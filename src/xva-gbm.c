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
 *
 *	xva-gbm.c
 *
 *	Geometric Brownian Motion Simulation for xVA Monte Carlo Framework
 *
 *	This module implements path simulation for assets following Geometric
 *	Brownian Motion (GBM), the foundation of exposure simulation in xVA.
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	GBM is the core stochastic process for modeling asset prices.
 *	The uncertainty enters through the diffusion term:
 *
 *	dS(t) = mu * S(t) * dt + sigma * S(t) * dW(t)
 *
 *	Where:
 *	- mu = drift (risk-free rate - dividend yield under Q-measure)
 *	- sigma = volatility (fixed calibrated input)
 *	- dW(t) = Wiener process increment ~ N(0, dt), the source of the
 *	  uncertainty, drawn with UxHwDoubleGaussDist() in the kernels
 *
 *	Under the risk-neutral measure Q:
 *	mu = r - q (risk-free rate minus dividend yield)
 *
 *	The volatility sigma scales the increments and so controls the width
 *	of the distribution of future asset prices, and hence of exposures.
 *
 *	ACADEMIC REFERENCES:
 *
 *	[1] Black, F. & Scholes, M. (1973). "The Pricing of Options and
 *	    Corporate Liabilities." Journal of Political Economy, 81(3).
 *	    - Foundation of risk-neutral pricing under GBM
 *
 *	[2] Glasserman, P. (2004). "Monte Carlo Methods in Financial Engineering."
 *	    Springer, Chapter 3: Generating Sample Paths.
 *	    - Euler-Maruyama discretization, exact simulation
 *
 *	[3] Gregory, J. (2020). "The xVA Challenge." Chapter 7: Exposure
 *	    and Credit Risk Models.
 *	    - Application of path simulation to xVA
 *
 *	[4] Kloeden, P.E. & Platen, E. (1992). "Numerical Solution of
 *	    Stochastic Differential Equations." Springer.
 *	    - Rigorous treatment of SDE discretization
 *
 *	PUBLIC IMPLEMENTATIONS REFERENCED:
 *
 *	- QuantLib: Monte Carlo simulation engines
 *	  https://github.com/lballabio/QuantLib
 *	- Open Source Risk Engine: Exposure simulation
 *	  https://github.com/opensourcerisk/Engine
 */

#include "xva-config.h"
#include <math.h>

/*
 *	xvaGetGbmPathParameters
 *
 *	Derive the per-step constants shared by the Monte Carlo and UxHw path
 *	simulators.
 *
 *	ALGORITHM:
 *	GBM has a closed-form solution:
 *
 *	S(t) = S(0) * exp((mu - sigma^2/2) * t + sigma * W(t))
 *
 *	For discrete time steps, this becomes:
 *
 *	S(t + dt) = S(t) * exp((mu - sigma^2/2) * dt + sigma * sqrt(dt) * Z)
 *
 *	where Z ~ N(0,1)
 *
 *	This is the EXACT simulation method (no discretization error)
 *	for GBM because we use the analytical solution.
 *
 *	Reference: Glasserman (2004), Section 3.2
 *
 *	This routine only computes `driftDt`, `diffusionSqrtDt` and `dt`. The
 *	paths are not stored: each mode-specific kernel (`xva-monte-carlo.c`
 *	and `xva-uxhw.c`) advances its paths one step at a time and folds each
 *	step straight into the exposure profile, so only the current spot of
 *	each path is ever held in memory.
 *
 *	Parameters:
 *	  market     - Market data (rates, volatility, maturity)
 *	  parameters - Output: per-step drift and diffusion constants
 */
void
xvaGetGbmPathParameters(
	const XvaMarketData *   market,
	XvaGbmPathParameters *  parameters)
{
	double  rr          = market->riskFreeRate;
	double  qq          = market->dividendYield;
	double  sigma       = market->volatility;
	double  maturity    = market->maturity;

	size_t  numSteps    = (size_t) kXvaNumTimeSteps;
	double  dt          = maturity / (double) numSteps;

	/*
	 *	Pre-compute constants for the exact GBM solution
	 */
	parameters->driftDt         = (rr - qq - 0.5 * sigma * sigma) * dt;
	parameters->diffusionSqrtDt = sigma * sqrt(dt);
	parameters->numSteps        = numSteps;
	parameters->dt              = dt;

	return;
}
