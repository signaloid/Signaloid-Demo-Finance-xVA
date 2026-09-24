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
 *	xva-exposure.c
 *
 *	Exposure Calculation Module for xVA Monte Carlo Framework
 *
 *	This module computes exposure metrics from simulated paths, including:
 *	- Expected Exposure (EE)
 *	- Expected Positive Exposure (EPE)
 *	- Expected Negative Exposure (ENE)
 *	- Potential Future Exposure (PFE)
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	Exposure metrics capture the uncertainty in the value that could
 *	be lost upon counterparty default. The metrics are computed from
 *	the distribution of mark-to-market (MtM) values:
 *
 *	- EPE(t) = E[max(V(t), 0)] - Expected positive exposure at time t
 *	  Captures the expected credit exposure (what we could lose if
 *	  counterparty defaults)
 *
 *	- ENE(t) = E[min(V(t), 0)] - Expected negative exposure at time t
 *	  Captures what counterparty could lose (relevant for DVA)
 *
 *	- PFE(t, q) = Quantile_q(max(V(t), 0)) - Potential future exposure
 *	  A tail risk measure (e.g., 97.5th percentile)
 *
 *	The uncertainty in these metrics flows directly from the uncertainty
 *	in the underlying asset paths (GBM simulation).
 *
 *	ACADEMIC REFERENCES:
 *
 *	[1] Gregory, J. (2020). "The xVA Challenge." Chapter 10: Exposure.
 *	    - Comprehensive treatment of exposure metrics
 *	    - EPE, ENE, PFE definitions and computation
 *
 *	[2] Pykhtin, M. & Zhu, S. (2007). "A Guide to Modelling Counterparty
 *	    Credit Risk." GARP Risk Review, July/August.
 *	    - Industry standard exposure modeling practices
 *
 *	[3] Canabarro, E. & Duffie, D. (2003). "Measuring and Marking
 *	    Counterparty Risk." In Asset/Liability Management for
 *	    Financial Institutions. Institutional Investor Books.
 *	    - Foundational paper on counterparty exposure
 *
 *	[4] BCBS (2015). "SA-CCR: Standardised Approach for Measuring
 *	    Counterparty Credit Risk Exposures." Basel Committee.
 *	    - Regulatory exposure calculations
 *
 *	PUBLIC IMPLEMENTATIONS REFERENCED:
 *
 *	- Open Source Risk Engine (ORE):
 *	  https://github.com/opensourcerisk/Engine
 *	  File: OREAnalytics/orea/aggregation/exposurecalculator.cpp
 *
 *	- xVA R Package: https://github.com/xvacode/xVA
 */

#include "xva-config.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/*
 *	xvaCompareDoubles
 *
 *	Comparison function for qsort, used in quantile calculation.
 */
int
xvaCompareDoubles(const void * aa, const void *  bb)
{
	double  valueA  = *(const double *) aa;
	double  valueB  = *(const double *) bb;

	if (valueA < valueB)
	{
		return -1;
	}

	if (valueA > valueB)
	{
		return 1;
	}

	return 0;
}

/*
 *	xvaQuantile
 *
 *	Compute the q-th quantile of a data array.
 *
 *	ALGORITHM:
 *	Uses linear interpolation between order statistics.
 *
 *	For quantile q in [0,1] with n data points:
 *	index = q * (n - 1)
 *	lower = floor(index)
 *	upper = ceil(index)
 *	quantile = data[lower] + (index - lower) * (data[upper] - data[lower])
 *
 *	Reference: Glasserman (2004), Appendix
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	The sample quantile is itself a random variable with uncertainty
 *	that depends on the sample size and the underlying distribution.
 *	For large samples, the standard error of the q-th quantile is
 *	approximately:
 *
 *	SE(q-quantile) approx sqrt(q*(1-q)/(n*f(F^{-1}(q))^2))
 *
 *	where f is the density and F is the CDF.
 *
 *	Parameters:
 *	  data - Array of values (will be sorted)
 *	  n    - Number of elements
 *	  q    - Quantile level in [0, 1]
 *
 *	Returns: Estimated quantile value
 */
double
xvaQuantile(double * data, size_t n, double q)
{
	double  indexValue;
	size_t  lowerIndex;
	size_t  upperIndex;
	double  fraction;

	if (data == NULL || n == 0)
	{
		return 0.0;
	}

	if (q <= 0.0)
	{
		return data[0];
	}

	if (q >= 1.0)
	{
		return data[n - 1];
	}

	/*
	 *	Sort the data
	 */
	qsort(data, n, sizeof(double), xvaCompareDoubles);

	/*
	 *	Compute quantile with linear interpolation
	 */
	indexValue  = q * (double) (n - 1);
	lowerIndex  = (size_t) floor(indexValue);
	upperIndex  = (size_t) ceil(indexValue);

	if (lowerIndex == upperIndex)
	{
		return data[lowerIndex];
	}

	fraction = indexValue - (double) lowerIndex;

	return data[lowerIndex] + fraction * (data[upperIndex] - data[lowerIndex]);
}

/*
 *	xvaComputeForwardMtm
 *
 *	Compute the mark-to-market value of a forward contract.
 *
 *	FORMULA:
 *	V(t) = S(t) * exp(-q*(T-t)) - K * exp(-r*(T-t))
 *
 *	where K is the forward price set at inception.
 *
 *	Reference: Hull (2017), Chapter 5
 *
 *	Parameters:
 *	  spot           - Current spot price S(t)
 *	  strike         - Forward strike price K
 *	  riskFreeRate   - Risk-free rate
 *	  dividendYield  - Dividend yield
 *	  timeToMaturity - Time to maturity (T - t)
 *
 *	Returns: Mark-to-market value
 */
double
xvaComputeForwardMtm(
	double  spot,
	double  strike,
	double  riskFreeRate,
	double  dividendYield,
	double  timeToMaturity)
{
	if (timeToMaturity <= 0.0)
	{
		/*
		 *	At maturity, payoff is simply S - K
		 */
		return spot - strike;
	}

	return spot * exp(-dividendYield * timeToMaturity) - strike * exp(-riskFreeRate * timeToMaturity);
}

/*
 *	xvaAllocateExposureProfile
 *
 *	Allocate the four exposure arrays (EE, EPE, ENE, PFE) of an exposure
 *	profile and zero them.
 *
 *	The exposure metrics themselves are:
 *
 *	1. EE(t)  = mean(MtM)              - expected exposure
 *	2. EPE(t) = mean(max(MtM, 0))      - positive exposure
 *	3. ENE(t) = mean(min(MtM, 0))      - negative exposure
 *	4. PFE(t) = quantile(max(MtM, 0))  - tail exposure
 *
 *	Reference: Gregory (2020), Chapter 10, Equations (10.1)-(10.4)
 *
 *	How those means and quantiles are formed differs fundamentally between
 *	execution modes -- Monte Carlo averages across a path ensemble, while
 *	UxHw integrates a single distributional path -- so the profile
 *	computation lives in the mode-specific kernels (`xva-monte-carlo.c`
 *	and `xva-uxhw.c`) and only the allocation is shared here.
 *
 *	Parameters:
 *	  profile  - Output: XvaExposureProfile structure to populate
 *	  numSteps - Number of time steps in the profile
 *
 *	Returns: 0 on success, EXIT_FAILURE on error
 */
int
xvaAllocateExposureProfile(
	XvaExposureProfile *    profile,
	size_t                  numSteps)
{
	if (profile == NULL)
	{
		return EXIT_FAILURE;
	}

	profile->numTimeSteps               = numSteps;
	profile->expectedExposure           = (double *) calloc(numSteps, sizeof(double));
	profile->expectedPositiveExposure   = (double *) calloc(numSteps, sizeof(double));
	profile->expectedNegativeExposure   = (double *) calloc(numSteps, sizeof(double));
	profile->potentialFutureExposure    = (double *) calloc(numSteps, sizeof(double));

	if (profile->expectedExposure == NULL ||
	    profile->expectedPositiveExposure == NULL ||
	    profile->expectedNegativeExposure == NULL ||
	    profile->potentialFutureExposure == NULL)
	{
		xvaFreeExposureProfile(profile);

		return EXIT_FAILURE;
	}

	return 0;
}

/*
 *	xvaFreeExposureProfile
 *
 *	Free memory allocated for an exposure profile.
 *
 *	Parameters:
 *	  profile - XvaExposureProfile structure to free
 */
void
xvaFreeExposureProfile(XvaExposureProfile * profile)
{
	if (profile == NULL)
	{
		return;
	}

	if (profile->expectedExposure != NULL)
	{
		free(profile->expectedExposure);
		profile->expectedExposure = NULL;
	}

	if (profile->expectedPositiveExposure != NULL)
	{
		free(profile->expectedPositiveExposure);
		profile->expectedPositiveExposure = NULL;
	}

	if (profile->expectedNegativeExposure != NULL)
	{
		free(profile->expectedNegativeExposure);
		profile->expectedNegativeExposure = NULL;
	}

	if (profile->potentialFutureExposure != NULL)
	{
		free(profile->potentialFutureExposure);
		profile->potentialFutureExposure = NULL;
	}

	profile->numTimeSteps = 0;

	return;
}
