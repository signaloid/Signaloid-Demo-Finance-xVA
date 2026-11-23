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
 *	xva-cva.c
 *
 *	Credit Valuation Adjustment (CVA) Calculation Module
 *
 *	This module implements the calculation of Credit Valuation Adjustment,
 *	which represents the market value of counterparty credit risk.
 *
 *	CVA DEFINITION:
 *	================================
 *	CVA is the price of counterparty credit risk - the expected loss
 *	from counterparty default weighted by the probability of default.
 *
 *	CVA = (1 - R) * integral_0^T EPE(t) * lambda(t) * exp(-integral_0^t lambda(s)ds) * DF(t) dt
 *
 *	where:
 *	- R = recovery rate
 *	- EPE(t) = expected positive exposure E[max(V(t), 0)] at time t
 *	- lambda(t) = hazard rate (default intensity) at time t
 *	- DF(t) = discount factor to time t
 *
 *	ACADEMIC REFERENCES:
 *
 *	[1] Gregory, J. (2020). "The xVA Challenge." Chapter 14: CVA.
 *	[2] Pykhtin, M. & Zhu, S. (2007). "A Guide to Modelling Counterparty Credit Risk."
 *	[3] Brigo, D. & Capponi, A. (2010). "Bilateral Counterparty Risk."
 *	[4] Hull, J. & White, A. (2012). "CVA and Wrong-Way Risk."
 */

#include "xva-config.h"
#include <math.h>
#include <stdlib.h>

/*
 *	xvaComputeCva
 *
 *	Compute Credit Valuation Adjustment using the unilateral formula.
 *
 *	FORMULA (Discrete approximation):
 *
 *	CVA = (1 - R) * sum_i EPE(t_i) * DF(t_i) * PD(t_{i-1}, t_i)
 *
 *	Reference: Gregory (2020), Equation (14.3)
 */
double
xvaComputeCva(
	const XvaExposureProfile *  profile,
	const XvaCreditData *       credit,
	const XvaMarketData *       market)
{
	if (profile == NULL || credit == NULL || market == NULL)
	{
		return 0.0;
	}

	double  rr          = market->riskFreeRate;
	double  lambda      = credit->hazardRate;
	double  recovery    = credit->recoveryRate;
	double  maturity    = market->maturity;
	double  dt          = maturity / (double) (profile->numTimeSteps - 1);

	/*
	 *	Loss given default
	 */
	double lgd = 1.0 - recovery;

	/*
	 *	Accumulate CVA using discrete summation
	 */
	double cva = 0.0;

	for (size_t ii = 1; ii < profile->numTimeSteps; ii++)
	{
		double  tPrev   = (double) (ii - 1) * dt;
		double  tCurr   = (double) ii * dt;

		/*
		 *	Expected positive exposure at current time
		 */
		double epe = profile->expectedPositiveExposure[ii];

		/*
		 *	Discount factor
		 */
		double df = xvaDiscountFactor(rr, tCurr);

		/*
		 *	Default probability in current interval
		 */
		double pd = xvaDefaultProbability(lambda, tPrev, tCurr);

		/*
		 *	Contribution to CVA
		 */
		cva += epe * df * pd;
	}

	return lgd * cva;
}
