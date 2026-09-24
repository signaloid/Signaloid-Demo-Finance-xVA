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
 *	@brief	Simulate a single distributional Geometric Brownian Motion path
 *		using the exact (closed-form) solution. On the Signaloid Cloud
 *		Compute Engine the `UxHwDoubleGaussDist()` call in the stepping
 *		loop yields the full standard normal distribution, so the one
 *		path carries the whole scenario ensemble.
 *
 *	@param	paths	: Output: a one-path ensemble holding the distributional path.
 *	@param	market	: Market data (spot, rates, volatility, maturity).
 *	@return	int	: Returns 0 on success, `EXIT_FAILURE` on error.
 */
int
xvaSimulateGbmPathsUxHw(
	XvaSimulationPaths *    paths,
	const XvaMarketData *   market);

/**
 *	@brief	Compute the exposure profile of a single distributional path.
 *		EE, EPE and ENE are read off the mark-to-market distribution as
 *		first moments rather than being averaged across paths, and PFE
 *		is the high quantile of the positive-exposure distribution.
 *
 *	@param	profile	: Output: the exposure profile.
 *	@param	paths	: The one-path distributional ensemble.
 *	@param	market	: Market data (rates, dividend yield, maturity).
 *	@param	strike	: Strike price of the derivative.
 *	@return	int	: Returns 0 on success, `EXIT_FAILURE` on error.
 */
int
xvaComputeExposureProfileUxHw(
	XvaExposureProfile *        profile,
	const XvaSimulationPaths *  paths,
	const XvaMarketData *       market,
	double                      strike);

/**
 *	@brief	Price the European put under UxHw execution. Unlike the other
 *		xVA outputs, the put is a closed form over the market data and
 *		never passes through the exposure profile, so nothing has collapsed
 *		it. When the volatility input is distributional (for example loaded
 *		from a CSV) the Black-Scholes result is a distribution too, so take
 *		its first moment here to honour the scalar output contract declared
 *		in `main.c` and `signaloid.yaml`.
 *
 *	@param	market		: Market data (spot, rates, volatility, maturity).
 *	@param	strikePrice	: Strike price of the option.
 *	@return	double		: Returns the European put price as a scalar.
 */
double
xvaEuropeanPutUxHw(
	const XvaMarketData *   market,
	double                  strikePrice);
