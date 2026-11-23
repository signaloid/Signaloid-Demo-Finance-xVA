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
 *	@brief	Simulate a Monte Carlo ensemble of Geometric Brownian Motion
 *		paths using the exact (closed-form) solution and compute their
 *		exposure profile in the same pass. Each of the `numPaths` paths
 *		is an independent scalar realisation: the
 *		`UxHwDoubleGaussDist()` call in the stepping loop returns a
 *		fresh pseudo-random sample per call under the native
 *		compatibility layer. EE, EPE and ENE are cross-sectional sample
 *		means over the paths at each time step, and PFE is the
 *		cross-sectional sample quantile of the positive exposures. Only
 *		the current step's spot of each path is held, not the whole
 *		path matrix.
 *
 *	@param	profile		: Output: the exposure profile.
 *	@param	market		: Market data (spot, rates, volatility, maturity).
 *	@param	strike		: Strike price of the derivative.
 *	@param	numPaths	: Number of Monte Carlo paths to simulate.
 *	@return	int		: Returns 0 on success, `EXIT_FAILURE` on error.
 */
int
xvaComputeExposureProfileMonteCarlo(
	XvaExposureProfile *        profile,
	const XvaMarketData *       market,
	double                      strike,
	size_t                      numPaths);

/**
 *	@brief	Price the European put under Monte Carlo execution. The
 *		Black-Scholes formula is a closed form over the market scalars, so
 *		the result is already a plain scalar and needs no collapse.
 *
 *	@param	market		: Market data (spot, rates, volatility, maturity).
 *	@param	strikePrice	: Strike price of the option.
 *	@return	double		: Returns the European put price.
 */
double
xvaEuropeanPutMonteCarlo(
	const XvaMarketData *   market,
	double                  strikePrice);
