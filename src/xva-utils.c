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
 *	xva-utils.c
 *
 *	Utility Functions for xVA Monte Carlo Framework
 *
 *	This module provides common utility functions used across the
 *	xVA calculation modules, including mathematical helpers and
 *	data manipulation routines.
 *
 *	ACADEMIC REFERENCES:
 *
 *	[1] Glasserman, P. (2004). "Monte Carlo Methods in Financial
 *	    Engineering." Springer, Appendix A.
 *	    - Numerical methods and utilities
 *
 *	[2] Abramowitz, M. & Stegun, I.A. (1964). "Handbook of Mathematical
 *	    Functions." National Bureau of Standards.
 *	    - Standard mathematical approximations
 */

#include "xva-config.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdio.h>

/*
 *	xvaDiscountFactor
 *
 *	Compute risk-free discount factor.
 *
 *	FORMULA:
 *	DF(t) = exp(-r * t)
 */
double
xvaDiscountFactor(double rate, double time)
{
	if (time < 0.0)
	{
		return 1.0;
	}

	return exp(-rate * time);
}

/*
 *	xvaSurvivalProbability
 *
 *	Compute survival probability to time t given constant hazard rate.
 *
 *	FORMULA:
 *	Q(t) = P(tau > t) = exp(-lambda * t)
 *
 *	Reference: Gregory (2020), Chapter 8
 */
double
xvaSurvivalProbability(double hazardRate, double time)
{
	if (hazardRate < 0.0 || time < 0.0)
	{
		return 1.0;
	}

	return exp(-hazardRate * time);
}

/*
 *	xvaDefaultProbability
 *
 *	Compute probability of default between times t1 and t2.
 *
 *	FORMULA:
 *	P(tau in [t1, t2]) = Q(t1) - Q(t2) = exp(-lambda*t1) - exp(-lambda*t2)
 *
 *	Reference: Gregory (2020), Chapter 8
 */
double
xvaDefaultProbability(double hazardRate, double t1, double t2)
{
	if (t2 <= t1)
	{
		return 0.0;
	}

	return xvaSurvivalProbability(hazardRate, t1) -
	       xvaSurvivalProbability(hazardRate, t2);
}

/*
 *	xvaStandardNormalCdf
 *
 *	Compute the cumulative distribution function of the standard normal.
 *
 *	Evaluates Phi(x) = (1 + erf(x / sqrt(2))) / 2 with the Abramowitz &
 *	Stegun rational approximation to erf (formula 7.1.26), whose absolute
 *	error is at most 1.5e-7.
 *
 *	Reference: Abramowitz & Stegun (1964), Section 7.1
 *
 *	Parameters:
 *	  x - Point at which to evaluate
 *
 *	Returns: P(Z <= x) where Z ~ N(0,1)
 */
double
xvaStandardNormalCdf(double x)
{
	static const double coeffA1 = 0.254829592;
	static const double coeffA2 = -0.284496736;
	static const double coeffA3 = 1.421413741;
	static const double coeffA4 = -1.453152027;
	static const double coeffA5 = 1.061405429;
	static const double coeffP  = 0.3275911;

	int     sign    = (x < 0) ? -1 : 1;
	double  absX    = fabs(x) / sqrt(2.0);
	double  tValue  = 1.0 / (1.0 + coeffP * absX);
	double  yValue  = 1.0 - (((((coeffA5 * tValue + coeffA4) * tValue) + coeffA3) * tValue + coeffA2) * tValue + coeffA1) * tValue *
	                  exp(-absX * absX);

	return 0.5 * (1.0 + sign * yValue);
}

/*
 *	xvaBlackScholesPut
 *
 *	Compute Black-Scholes-Merton price for a European put option.
 *
 *	P = K * exp(-r*T) * N(-d2) - S * exp(-q*T) * N(-d1)
 *
 *	Reference: Black & Scholes (1973), Equation 27, extended with
 *	continuous dividend yield per Merton (1973).
 *
 *	The original Black-Scholes formula (q = 0):
 *	P = K * exp(-r*T) * N(-d2) - S * N(-d1)
 *
 *	Parameters:
 *	  spotPrice, strikePrice, riskFreeRate, dividendYield, volatility, timeToMaturity
 *	    - Same as xvaBlackScholesCall
 *
 *	Returns: Put option price
 */
double
xvaBlackScholesPut(
	double  spotPrice,
	double  strikePrice,
	double  riskFreeRate,
	double  dividendYield,
	double  volatility,
	double  timeToMaturity)
{
	if (timeToMaturity <= 0.0 || volatility <= 0.0 || spotPrice <= 0.0 || strikePrice <= 0.0)
	{
		return 0.0;
	}

	double  sqrtT   = sqrt(timeToMaturity);
	double  d1      = (log(spotPrice / strikePrice) + (riskFreeRate - dividendYield + 0.5 * volatility * volatility) * timeToMaturity) /
	                  (volatility * sqrtT);
	double  d2          = d1 - volatility * sqrtT;
	double  putPrice    = strikePrice * exp(-riskFreeRate * timeToMaturity) * xvaStandardNormalCdf(-d2) -
	                      spotPrice * exp(-dividendYield * timeToMaturity) * xvaStandardNormalCdf(-d1);

	return putPrice;
}

/*
 *	xvaHazardRateFromSpread
 *
 *	Compute the hazard rate from a credit spread.
 *
 *	lambda approximately equals spread / (1 - R)
 *
 *	where R is the recovery rate.
 *
 *	Reference: Gregory (2020), Chapter 8
 *
 *	Parameters:
 *	  spread       - Credit spread (decimal, e.g., 0.01 = 100 bps)
 *	  recoveryRate - Recovery rate (0 to 1)
 *
 *	Returns: Hazard rate (instantaneous default intensity)
 */
double
xvaHazardRateFromSpread(double spread, double recoveryRate)
{
	if (recoveryRate >= 1.0)
	{
		return 0.0;  /*	No loss given default */
	}

	return spread / (1.0 - recoveryRate);
}
