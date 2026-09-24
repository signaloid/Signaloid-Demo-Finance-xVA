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
 *	xva-mva.c
 *
 *	Margin Valuation Adjustment (MVA) Calculation Module
 *
 *	This module implements the calculation of Margin Valuation Adjustment,
 *	which captures the cost of posting initial margin.
 *
 *	MVA DEFINITION:
 *	================================
 *	MVA is the cost of funding initial margin (IM) that must be posted
 *	under bilateral margin agreements or for clearing.
 *
 *	MVA = integrate_0^T IM(t) * s_f(t) * Q(t) * DF(t) dt
 *
 *	where:
 *	- IM(t) = expected initial margin at time t
 *	- s_f = funding spread
 *	- Q(t) = survival probability
 *	- DF(t) = discount factor
 *
 *	Initial margin is a regulatory requirement (UMR - Uncleared Margin
 *	Rules) that requires parties to post collateral covering potential
 *	future exposure over a margin period of risk (typically 10 days).
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	MVA uncertainty comes from:
 *
 *	1. INITIAL MARGIN UNCERTAINTY:
 *	   - IM is typically calculated as VaR or ES over the MPOR
 *	   - Higher volatility => higher IM => higher MVA
 *	   - IM(t) = VaR_{99%}(DeltaV | t, MPOR)
 *
 *	2. FUNDING COST UNCERTAINTY:
 *	   - s_f varies with bank's credit quality
 *	   - Funding the IM position requires borrowing
 *
 *	3. IM MODEL UNCERTAINTY:
 *	   - SIMM vs. schedule-based IM
 *	   - Different IM models give different amounts
 *
 *	ACADEMIC REFERENCES:
 *
 *	[1] Gregory, J. (2020). "The xVA Challenge." Chapter 17: MVA.
 *	    - Comprehensive MVA framework
 *
 *	[2] ISDA (2016). "ISDA SIMM: Standard Initial Margin Model."
 *	    - Industry standard IM calculation
 *
 *	[3] Andersen, L., Pykhtin, M. & Sokol, A. (2017). "Does Initial
 *	    Margin Eliminate Counterparty Risk?" Risk Magazine.
 *	    - Analysis of IM effectiveness
 *
 *	[4] BCBS-IOSCO (2015). "Margin Requirements for Non-Centrally
 *	    Cleared Derivatives." Basel Committee.
 *	    - Regulatory framework for IM
 *
 *	[5] Green, A. & Kenyon, C. (2015). "MVA: Initial Margin Valuation
 *	    Adjustment by Replication and Regression." Risk Magazine.
 *	    - MVA calculation methodology
 *
 *	PUBLIC IMPLEMENTATIONS REFERENCED:
 *
 *	- Open Source Risk Engine (ORE):
 *	  https://github.com/opensourcerisk/Engine
 *	  OREAnalytics/orea/aggregation/dimcalculator.cpp
 */

#include "xva-config.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
 *	xvaComputeMva
 *
 *	Compute Margin Valuation Adjustment.
 *
 *	MVA is the cost of funding the initial margin requirement over
 *	the life of the trade.
 *
 *	FORMULA:
 *
 *	MVA = integrate_0^T E[IM(t)] * s_f * Q(t) * DF(t) dt
 *
 *	Reference: Gregory (2020), Equation (17.1)
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	MVA combines:
 *
 *	1. IM uncertainty (from VaR calculation)
 *	   - Higher volatility => higher IM => higher MVA
 *
 *	2. Funding cost uncertainty
 *	   - Higher funding spread => higher MVA
 *
 *	3. Credit uncertainty (survival probability)
 *	   - IM only needs funding while parties survive
 *
 *	Parameters:
 *	  profile - Exposure profile
 *	  funding - Funding data (funding spread)
 *	  capital - Capital data (for regulatory parameters)
 *	  credit  - Credit data (for survival probability)
 *	  market  - Market data
 *
 *	Returns: Margin Valuation Adjustment
 */
double
xvaComputeMva(
	const XvaExposureProfile *  profile,
	const XvaFundingData *      funding,
	const XvaCapitalData *      capital,
	const XvaCreditData *       credit,
	const XvaMarketData *       market)
{
	double  riskFreeRate;
	double  maturity;
	double  fundingSpread;
	size_t  numSteps;
	double  dt;
	double  lambda;
	double  mva;

	if (profile == NULL || funding == NULL || market == NULL)
	{
		return 0.0;
	}

	/*
	 *	Capital parameter reserved for future SIMM-based IM calculation
	 */
	(void) capital;

	riskFreeRate    = market->riskFreeRate;
	maturity        = market->maturity;
	fundingSpread   = funding->fundingSpread;
	numSteps        = profile->numTimeSteps;
	dt              = maturity / (double) (numSteps - 1);

	/*
	 *	Hazard rate for survival weighting.
	 *	Uses actual credit data when available; falls back to funding spread.
	 */
	lambda = (credit != NULL) ? credit->hazardRate : fundingSpread;

	/*
	 *	Compute IM profile
	 *	In practice, this would use PFE or SIMM-based calculation
	 *	Here we approximate using PFE as proxy for IM
	 */
	mva = 0.0;

	for (size_t ii = 1; ii < numSteps; ii++)
	{
		double tt = (double) ii * dt;

		/*
		 *	Approximate IM using PFE (a common simplification)
		 *	IM approx PFE scaled for MPOR
		 *
		 *	Reference: Green & Kenyon (2015)
		 */
		double pfe = profile->potentialFutureExposure[ii];

		/*
		 *	Scale PFE to MPOR (10 days vs. full tenor)
		 *	IM ~ PFE * sqrt(MPOR/252)
		 */
		double  mporDays    = (double) kXvaImPeriodDays;
		double  mporScale   = sqrt(mporDays / (double) kXvaTradingDaysPerYear);
		double  im          = pfe * mporScale;

		/*
		 *	Survival probability
		 */
		double qq = xvaSurvivalProbability(lambda, tt);

		/*
		 *	Discount factor
		 */
		double df = xvaDiscountFactor(riskFreeRate, tt);

		/*
		 *	MVA contribution: IM * funding spread * survival * discount
		 */
		mva += im * fundingSpread * qq * df * dt;
	}

	return mva;
}
