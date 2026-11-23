/*
 *	Copyright (c) 2026, Signaloid.
 *
 *	Permission is hereby granted, free of charge, to any person obtaining a copy
 *	of this software and associated documentation files (the "Software"), to deal
 *	in the Software without restriction, including without limitation the rights
 *	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *	copies of the Software, and to permit persons to whom the Software is
 *	furnished to do so, subject to the following conditions:
 *
 *	The above copyright notice and this permission notice shall be included in all
 *	copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *	SOFTWARE.
 */

/*
 *	UxHw execution-mode kernel.
 *
 *	This file holds the distributional half of the xVA engine. Where the
 *	Monte Carlo kernel in `xva-monte-carlo.c` draws an ensemble of scalar
 *	paths and averages across it, this kernel simulates one path whose
 *	values carry a complete probability distribution, and integrates that
 *	distribution out in a single pass.
 *
 *	The two kernels share the GBM step constants in `xva-gbm.c`, the
 *	mark-to-market in `xva-exposure.c` and the mode-agnostic xVA metric
 *	computations in `xva-cva.c` and friends.
 */

#include "xva-config.h"
#include "xva-uxhw.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

int
xvaComputeExposureProfileUxHw(
	XvaExposureProfile *        profile,
	const XvaMarketData *       market,
	double                      strike)
{
	XvaGbmPathParameters    parameters;
	double                  riskFreeRate;
	double                  dividendYield;
	double                  maturity;
	double                  dt;
	double                  spotPrice;

	if (profile == NULL || market == NULL)
	{
		return EXIT_FAILURE;
	}

	xvaGetGbmPathParameters(market, &parameters);

	riskFreeRate    = market->riskFreeRate;
	dividendYield   = market->dividendYield;
	maturity        = market->maturity;
	dt              = parameters.dt;

	profile->numTimeSteps = parameters.numSteps + 1;  /*	Include t=0 */

	/*
	 *	A single distributional path represents the whole ensemble. It is
	 *	advanced one step at a time and each step is folded into the
	 *	exposure profile straight away, so only the current spot is held
	 *	rather than the whole path.
	 */
	spotPrice = market->spotPrice;

	/*
	 *	Compute exposure at each time step
	 */
	for (size_t ii = 0; ii < profile->numTimeSteps; ii++)
	{
		double  mtmValue;
		double  currentTime;
		double  timeToMaturity;
		double  probPositive;
		double  positiveExposureDistribution;
		double  negatedMtm;
		double  probNegative;

		/*
		 *	Time to maturity at this step
		 */
		currentTime     = ii * dt;
		timeToMaturity  = maturity - currentTime;

		if (timeToMaturity < 0.0)
		{
			timeToMaturity = 0.0;
		}

		/*
		 *	Advance the path from the previous step using exact GBM
		 *	simulation. Step 0 is the deterministic spot.
		 */
		if (ii > 0)
		{
			double  zz = UxHwDoubleGaussDist(0.0, 1.0);

			/*
			 *	Break dependency tracking between time steps. Each
			 *	step depends only on the previous step's value, so
			 *	take an independent copy of the previous price before
			 *	propagating it forward.
			 */
			spotPrice = UxHwDoubleGetIndependentCopy(spotPrice) *
			            exp(parameters.driftDt + parameters.diffusionSqrtDt * zz);
		}

		/*
		 *	Compute mark-to-market value.
		 *	Using forward contract as example (easily extensible)
		 */
		mtmValue    = xvaComputeForwardMtm(spotPrice, strike, riskFreeRate, dividendYield, timeToMaturity);

		/*
		 *	Exposure metrics are expectations: EE = E[V],
		 *	EPE = E[max(V, 0)], ENE = E[min(V, 0)]. They are scalars
		 *	per time step.
		 *
		 *	Here `mtmValue` is a distribution over scenarios. A scalar
		 *	comparison would collapse it, so we build max(V, 0) and
		 *	min(V, 0) distributionally as a truncated-plus-Dirac
		 *	mixture (limit the distribution to the relevant half-line
		 *	and mix with a Dirac at 0 weighted by the probability of
		 *	being on that side), then extract the first moment so the
		 *	stored exposure is the scalar expectation. This integrates
		 *	out the path scatter in a single pass, exactly as the Monte
		 *	Carlo kernel does by averaging across its ensemble.
		 *
		 *	EE = E[V].
		 */
		profile->expectedExposure[ii] = UxHwDoubleNthMoment(mtmValue, 1);

		/*
		 *	Positive-exposure distribution max(V, 0): limit V to its
		 *	positive half-line and mix with a Dirac at 0 weighted by
		 *	P(V > 0). EPE(t) is its mean; PFE(t) is its high (97.5%)
		 *	quantile. Both are read off the same distribution, so the
		 *	tail exposure is preserved here rather than being lost when
		 *	we collapse EPE to a scalar -- matching the cross-sectional
		 *	quantile a Monte Carlo ensemble would compute.
		 */
		probPositive                    = UxHwDoubleProbabilityGT(mtmValue, 0.0);
		positiveExposureDistribution    = UxHwDoubleMixture(
			UxHwDoubleLimitDistributionSupport(mtmValue, 0.0, DBL_MAX),
			0.0,
			probPositive
		);

		profile->expectedPositiveExposure[ii]   = UxHwDoubleNthMoment(positiveExposureDistribution, 1);
		profile->potentialFutureExposure[ii]    = UxHwDoubleQuantile(positiveExposureDistribution, kXvaPfeQuantile);

		/*
		 *	ENE = E[min(V, 0)] = -E[max(-V, 0)]. We reuse the positive
		 *	[0, DBL_MAX] truncation on -V (rather than a [-DBL_MAX, 0]
		 *	truncation on V, which can yield a NaN-tainted mixture when
		 *	V has little mass below zero) and negate the mean.
		 */
		negatedMtm      = -mtmValue;
		probNegative    = UxHwDoubleProbabilityGT(negatedMtm, 0.0);

		/*
		 *	When there is no mass below zero (for example the
		 *	deterministic spot at step 0), the truncation of -V has no
		 *	support and returns NaN. In that case max(-V, 0) is the Dirac
		 *	at 0, so ENE is exactly 0.
		 */
		if (probNegative == 0.0)
		{
			profile->expectedNegativeExposure[ii] = 0.0;

			continue;
		}

		profile->expectedNegativeExposure[ii] = -UxHwDoubleNthMoment(
			UxHwDoubleMixture(
				UxHwDoubleLimitDistributionSupport(negatedMtm, 0.0, DBL_MAX),
				0.0,
				probNegative
			),
			1
		);
	}

	return 0;
}

double
xvaEuropeanPutUxHw(
	const XvaMarketData *   market,
	double                  strikePrice)
{
	double putPrice;

	putPrice = xvaBlackScholesPut(
		market->spotPrice,
		strikePrice,
		market->riskFreeRate,
		market->dividendYield,
		market->volatility,
		market->maturity
	);

	/*
	 *	Collapse to the expectation. With a deterministic volatility this
	 *	is the identity; with a distributional one it is the mean of the
	 *	resulting put distribution, matching how the exposure profile
	 *	collapses EE, EPE and ENE for every other output.
	 */
	return UxHwDoubleNthMoment(putPrice, 1);
}
