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
 *	Monte Carlo execution-mode kernel.
 *
 *	This file holds the sampling half of the xVA engine: the path
 *	simulation and the exposure profile as a Monte Carlo estimator
 *	evaluates them. Every quantity here is a scalar sample, and the
 *	expectations that define the exposure metrics are formed by averaging
 *	across the path ensemble.
 *
 *	The distributional counterpart lives in `xva-uxhw.c`; the two share
 *	the GBM step constants in `xva-gbm.c`, the mark-to-market in
 *	`xva-exposure.c` and the mode-agnostic xVA metric computations in
 *	`xva-cva.c` and friends.
 */

#include "xva-config.h"
#include "xva-monte-carlo.h"
#include <math.h>
#include <stdlib.h>

int
xvaComputeExposureProfileMonteCarlo(
	XvaExposureProfile *        profile,
	const XvaMarketData *       market,
	double                      strike,
	size_t                      numPaths)
{
	XvaGbmPathParameters    parameters;
	double                  riskFreeRate;
	double                  dividendYield;
	double                  maturity;
	double *                spotPrices;
	double *                positiveExposures;

	if (profile == NULL || market == NULL || numPaths == 0)
	{
		return EXIT_FAILURE;
	}

	xvaGetGbmPathParameters(market, &parameters);

	riskFreeRate    = market->riskFreeRate;
	dividendYield   = market->dividendYield;
	maturity        = market->maturity;

	profile->numTimeSteps = parameters.numSteps + 1;  /*	Include t=0 */

	/*
	 *	The paths are advanced one step at a time and each step is folded
	 *	into the exposure profile straight away, so only the current spot
	 *	of each path is kept, rather than the whole
	 *	`numPaths` x `numSteps + 1` path matrix. A second array of the same
	 *	size holds the positive exposures for the cross-sectional PFE
	 *	quantile; it is separate because the quantile sorts it in place.
	 *	Both arrays come from one allocation.
	 */
	spotPrices = (double *) malloc(2 * numPaths * sizeof(double));

	if (spotPrices == NULL)
	{
		return EXIT_FAILURE;
	}

	positiveExposures = spotPrices + numPaths;

	for (size_t jj = 0; jj < numPaths; jj++)
	{
		spotPrices[jj] = market->spotPrice;
	}

	/*
	 *	Compute exposure at each time step
	 */
	for (size_t ii = 0; ii < profile->numTimeSteps; ii++)
	{
		double  sumMtm      = 0.0;
		double  sumPositive = 0.0;
		double  sumNegative = 0.0;
		double  currentTime;
		double  timeToMaturity;
		double  pathCount;

		/*
		 *	Time to maturity at this step
		 */
		currentTime     = ii * parameters.dt;
		timeToMaturity  = maturity - currentTime;

		if (timeToMaturity < 0.0)
		{
			timeToMaturity = 0.0;
		}

		/*
		 *	Process each path.
		 *	THIS LOOP AGGREGATES UNCERTAINTY FROM ALL SIMULATED PATHS
		 */
		for (size_t jj = 0; jj < numPaths; jj++)
		{
			double  mtmValue;
			double  positiveExp;
			double  negativeExp;

			/*
			 *	Advance the path from the previous step using exact
			 *	GBM simulation. Each path is an independent scalar
			 *	realisation: UxHwDoubleGaussDist() returns a fresh
			 *	pseudo-random standard normal sample on every call
			 *	under the native compatibility layer. Step 0 is the
			 *	deterministic spot.
			 */
			if (ii > 0)
			{
				double zz = UxHwDoubleGaussDist(0.0, 1.0);

				spotPrices[jj] *= exp(parameters.driftDt + parameters.diffusionSqrtDt * zz);
			}

			/*
			 *	Compute mark-to-market value.
			 *	Using forward contract as example (easily extensible)
			 */
			mtmValue = xvaComputeForwardMtm(spotPrices[jj], strike, riskFreeRate, dividendYield, timeToMaturity);

			/*
			 *	`mtmValue` is a single scalar sample, so max(V, 0)
			 *	and min(V, 0) are plain comparisons and the
			 *	expectations are the sample means taken below.
			 */
			sumMtm      += mtmValue;
			positiveExp = mtmValue > 0.0 ? mtmValue : 0.0;
			negativeExp = mtmValue < 0.0 ? mtmValue : 0.0;

			sumPositive             += positiveExp;
			positiveExposures[jj]   = positiveExp;

			sumNegative += negativeExp;
		}

		/*
		 *	Compute expected values (sample means)
		 *
		 *	EE(t) = (1/N) * SUM MtM_i
		 *	EPE(t) = (1/N) * SUM max(MtM_i, 0)
		 *	ENE(t) = (1/N) * SUM min(MtM_i, 0)
		 *
		 *	Reference: Gregory (2020), Equations (10.1)-(10.3)
		 */
		pathCount = (double) numPaths;
		profile->expectedExposure[ii]           = sumMtm / pathCount;
		profile->expectedPositiveExposure[ii]   = sumPositive / pathCount;
		profile->expectedNegativeExposure[ii]   = sumNegative / pathCount;

		/*
		 *	Compute PFE as the cross-sectional quantile of the
		 *	positive exposures.
		 *
		 *	PFE(t, q) = Quantile_q{max(MtM, 0)}
		 *
		 *	Typically q = 0.975 (97.5th percentile)
		 *
		 *	Reference: Gregory (2020), Section 10.2
		 */
		profile->potentialFutureExposure[ii] =
			xvaQuantile(positiveExposures, numPaths, kXvaPfeQuantile);
	}

	free(spotPrices);

	return 0;
}

double
xvaEuropeanPutMonteCarlo(
	const XvaMarketData *   market,
	double                  strikePrice)
{
	return xvaBlackScholesPut(
		market->spotPrice,
		strikePrice,
		market->riskFreeRate,
		market->dividendYield,
		market->volatility,
		market->maturity
	);
}
