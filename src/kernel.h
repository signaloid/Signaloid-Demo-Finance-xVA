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

#pragma once

#include <stdlib.h>
#include "xva-config.h"

/**
 *	@brief	Market, credit and trade parameters and input variables consumed
 *		by the kernels. The command-line parser populates it directly, and
 *		the caller overwrites `inputVariables` when reading inputs from a
 *		CSV file, so the kernels do not depend on the command-line
 *		handling.
 */
typedef struct
{
	double  spotPrice;
	double  riskFreeRate;
	double  dividendYield;
	double  maturity;
	double  creditSpread;
	double  ownCreditSpread;
	double  strikePrice;

	/*
	 *	Input variables, indexed by `XvaInputIndex`. The volatility
	 *	`inputVariables[kXvaInputIndexVolatility]` is a calibrated market
	 *	input, held at a fixed deterministic value (`-V`, or the demo
	 *	default) so that the Monte Carlo and UxHw modes compute the same
	 *	production xVA. It is a distribution only when read from an input
	 *	CSV file (`-i`). The stochasticity that drives exposure is the
	 *	Brownian path ensemble, not the volatility itself.
	 */
	double inputVariables[kXvaInputIndexMax];
} KernelParameters;

/**
 *	@brief	Monte Carlo calculation kernel. Simulates
 *		`numberOfMonteCarloIterations` independent
 *		Brownian paths, forms the exposure profile as their
 *		cross-sectional average, and computes the selected xVA metric(s)
 *		from it using scalar sample statistics.
 *
 *	@param	parameters			: The model parameters.
 *	@param	outputSelect			: Which output to compute, or `kXvaOutputIndexMax` for all.
 *	@param	numberOfMonteCarloIterations	: Number of Brownian paths to simulate.
 *	@param	outputVariables			: The output variables.
 */
void
calculateOutputMonteCarlo(
	const KernelParameters *    parameters,
	size_t                      outputSelect,
	size_t                      numberOfMonteCarloIterations,
	double *                    outputVariables);

/**
 *	@brief	UxHw calculation kernel. Simulates a single distributional
 *		Brownian path, forms the exposure profile by integrating that
 *		distribution out with UxHw distributional arithmetic, and
 *		computes the selected xVA metric(s) from it.
 *
 *	@param	parameters	: The model parameters.
 *	@param	outputSelect	: Which output to compute, or `kXvaOutputIndexMax` for all.
 *	@param	outputVariables	: The output variables.
 */
void
calculateOutputUxHw(
	const KernelParameters *    parameters,
	size_t                      outputSelect,
	double *                    outputVariables);
