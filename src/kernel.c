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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kernel.h"
#include "xva-config.h"
#include "xva-monte-carlo.h"
#include "xva-uxhw.h"

/*
 *	The model data (market, credit, funding, capital) and the xVA metrics
 *	themselves are mode-agnostic: they are ordinary arithmetic over the
 *	exposure profile and are identical whether that profile came from a
 *	Monte Carlo ensemble or from a single distributional path. Only the
 *	path simulation and the exposure profile differ between modes, and
 *	those live in `xva-monte-carlo.c` and `xva-uxhw.c`.
 *
 *	The two shared stages below therefore bracket the mode-specific
 *	kernels: `buildModelData()` runs before them and
 *	`computeXvaMetrics()` after.
 */
typedef struct
{
	XvaMarketData   market;
	XvaCreditData   cptyCredit;
	XvaCreditData   ownCredit;
	XvaFundingData  funding;
	XvaCapitalData  capital;
} XvaModelData;

/**
 *	@brief	Build the market, credit, funding, and capital data from the
 *		model parameters and input variables.
 *
 *	@param	parameters	: The model parameters and input variables.
 *	@param	model		: Output: the assembled model data.
 */
static void
buildModelData(
	const KernelParameters *    parameters,
	XvaModelData *              model)
{
	/*
	 *	Build market data from the model parameters and input variables.
	 *	Volatility comes from `parameters->inputVariables`. It is a
	 *	scalar unless the UxHw kernel is run with a CSV of volatility
	 *	samples (`-i`), in which case it is a distribution.
	 */
	model->market.spotPrice     = parameters->spotPrice;
	model->market.riskFreeRate  = parameters->riskFreeRate;
	model->market.dividendYield = parameters->dividendYield;
	model->market.volatility    = parameters->inputVariables[kXvaInputIndexVolatility];
	model->market.maturity      = parameters->maturity;

	/*
	 *	Build counterparty credit data
	 */
	model->cptyCredit.creditSpread  = parameters->creditSpread;
	model->cptyCredit.recoveryRate  = kXvaDemoDefaultRecoveryRate;
	model->cptyCredit.hazardRate    = xvaHazardRateFromSpread(
		model->cptyCredit.creditSpread,
		model->cptyCredit.recoveryRate
	);

	/*
	 *	Build own credit data
	 */
	model->ownCredit.creditSpread   = parameters->ownCreditSpread;
	model->ownCredit.recoveryRate   = kXvaDemoDefaultRecoveryRate;
	model->ownCredit.hazardRate     = xvaHazardRateFromSpread(
		model->ownCredit.creditSpread,
		model->ownCredit.recoveryRate
	);

	/*
	 *	Build funding data. Lending at the risk-free rate gives a zero
	 *	lending spread, so the funding benefit (FBA) term of FVA vanishes
	 *	and FVA reduces to the funding cost (FCA).
	 */
	model->funding.fundingSpread    = kXvaDemoDefaultFundingSpread;
	model->funding.borrowingRate    = parameters->riskFreeRate + model->funding.fundingSpread;
	model->funding.lendingRate      = parameters->riskFreeRate;

	/*
	 *	Build capital data
	 */
	model->capital.capitalRatio = (double) kXvaCapitalRatioPercent / 100.0;
	model->capital.hurdleRate   = (double) kXvaHurdleRatePercent / 100.0;
	model->capital.rwaFactor    = kXvaDefaultRiskWeight;

	return;
}

/**
 *	@brief	Compute the selected xVA metric(s) from an exposure profile. Each
 *		output is a single converged production figure, because the
 *		expectation across scenarios was already formed when the
 *		exposure profile was built.
 *
 *	@param	outputSelect	: Which output to compute, or `kXvaOutputIndexMax` for all.
 *	@param	model		: The market, credit, funding, and capital data.
 *	@param	profile		: The exposure profile.
 *	@param	outputVariables	: The output variables.
 */
static void
computeXvaMetrics(
	size_t                      outputSelect,
	const XvaModelData *        model,
	const XvaExposureProfile *  profile,
	double *                    outputVariables)
{
	bool calculateAllOutputs = (outputSelect == kXvaOutputIndexMax);

	if (calculateAllOutputs || (outputSelect == kXvaOutputIndexCva))
	{
		outputVariables[kXvaOutputIndexCva] =
			xvaComputeCva(profile, &model->cptyCredit, &model->market);
	}

	if (calculateAllOutputs || (outputSelect == kXvaOutputIndexDva))
	{
		outputVariables[kXvaOutputIndexDva] =
			xvaComputeDva(profile, &model->ownCredit, &model->market);
	}

	if (calculateAllOutputs || (outputSelect == kXvaOutputIndexBilateralCva))
	{
		outputVariables[kXvaOutputIndexBilateralCva] =
			xvaComputeBilateralCva(profile, &model->cptyCredit, &model->ownCredit, &model->market);
	}

	if (calculateAllOutputs || (outputSelect == kXvaOutputIndexFva))
	{
		outputVariables[kXvaOutputIndexFva] =
			xvaComputeFva(profile, &model->funding, &model->cptyCredit, &model->ownCredit, &model->market);
	}

	if (calculateAllOutputs || (outputSelect == kXvaOutputIndexMva))
	{
		outputVariables[kXvaOutputIndexMva] =
			xvaComputeMva(profile, &model->funding, &model->capital, &model->cptyCredit, &model->market);
	}

	if (calculateAllOutputs || (outputSelect == kXvaOutputIndexKva))
	{
		outputVariables[kXvaOutputIndexKva] =
			xvaComputeKva(profile, &model->capital, &model->market);
	}

	if (calculateAllOutputs || (outputSelect == kXvaOutputIndexTotalXva))
	{
		double  cva = xvaComputeCva(profile, &model->cptyCredit, &model->market);
		double  dva = xvaComputeDva(profile, &model->ownCredit, &model->market);
		double  fva = xvaComputeFva(profile, &model->funding, &model->cptyCredit, &model->ownCredit, &model->market);
		double  mva = xvaComputeMva(profile, &model->funding, &model->capital, &model->cptyCredit, &model->market);
		double  kva = xvaComputeKva(profile, &model->capital, &model->market);

		outputVariables[kXvaOutputIndexTotalXva] =
			cva - dva + fva + mva + kva;
	}

	return;
}

/**
 *	@brief	Whether an output index is selected, either on its own or as part
 *		of a request for all outputs.
 *
 *	@param	outputSelect	: The selected output, or `kXvaOutputIndexMax` for all.
 *	@param	outputIndex	: The 0-indexed output.
 *	@return			: `true` if the output should be computed.
 */
static bool
isOutputSelected(
	size_t  outputSelect,
	size_t  outputIndex)
{
	return (outputSelect == kXvaOutputIndexMax) ||
	       (outputSelect == outputIndex);
}

void
calculateOutputMonteCarlo(
	const KernelParameters *    parameters,
	size_t                      outputSelect,
	size_t                      numberOfMonteCarloIterations,
	double *                    outputVariables)
{
	XvaModelData        model;
	XvaExposureProfile  profile;

	buildModelData(parameters, &model);

	/*
	 *	Simulate `numberOfMonteCarloIterations` independent Brownian paths.
	 *	The exposure profile is their cross-sectional average, so every xVA
	 *	metric below is computed once from the expected exposure profile.
	 *	That ordering matters for the nonlinear metrics: KVA's Effective EPE
	 *	is a running maximum over the *expected* exposure profile, so the
	 *	expectation across paths must be formed before the running max.
	 *	Averaging per-path xVAs would instead compute E[max(...)], a
	 *	different, biased value.
	 */
	if (xvaComputeExposureProfileMonteCarlo(&profile, &model.market, parameters->strikePrice, numberOfMonteCarloIterations) != 0)
	{
		fprintf(stderr, "Error: Exposure profile computation failed.\n");

		return;
	}

	computeXvaMetrics(outputSelect, &model, &profile, outputVariables);

	/*
	 *	The European put never passes through the exposure profile, so it
	 *	is computed per mode rather than in the shared helper above.
	 */
	if (isOutputSelected(outputSelect, kXvaOutputIndexEuropeanPut))
	{
		outputVariables[kXvaOutputIndexEuropeanPut] =
			xvaEuropeanPutMonteCarlo(&model.market, parameters->strikePrice);
	}

	return;
}

void
calculateOutputUxHw(
	const KernelParameters *    parameters,
	size_t                      outputSelect,
	double *                    outputVariables)
{
	XvaModelData        model;
	XvaExposureProfile  profile;

	buildModelData(parameters, &model);

	/*
	 *	One distributional path represents the whole ensemble: the exposure
	 *	profile integrates the scenario distribution out in a single pass,
	 *	giving the same expected exposure profile that the Monte Carlo
	 *	kernel reaches by averaging, and every xVA metric below is then
	 *	computed once from it.
	 */
	if (xvaComputeExposureProfileUxHw(&profile, &model.market, parameters->strikePrice) != 0)
	{
		fprintf(stderr, "Error: Exposure profile computation failed.\n");

		return;
	}

	computeXvaMetrics(outputSelect, &model, &profile, outputVariables);

	/*
	 *	The European put never passes through the exposure profile, so
	 *	nothing has collapsed it: `xvaEuropeanPutUxHw()` takes its first
	 *	moment so that a distributional volatility input does not leak a
	 *	distribution into an output declared scalar.
	 */
	if (isOutputSelected(outputSelect, kXvaOutputIndexEuropeanPut))
	{
		outputVariables[kXvaOutputIndexEuropeanPut] =
			xvaEuropeanPutUxHw(&model.market, parameters->strikePrice);
	}

	return;
}
