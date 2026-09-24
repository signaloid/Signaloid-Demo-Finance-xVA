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
 *	xva-fva.c
 *
 *	Funding Valuation Adjustment (FVA) Calculation Module
 *
 *	This module implements the calculation of Funding Valuation Adjustment,
 *	which captures the cost/benefit of funding derivative positions.
 *
 *	FVA DEFINITION:
 *	================================
 *	FVA represents the market value of funding costs (or benefits)
 *	associated with uncollateralized derivative positions.
 *
 *	FVA = FCA - FBA
 *
 *	where:
 *	- FCA (Funding Cost Adjustment) = Cost of funding positive exposure
 *	- FBA (Funding Benefit Adjustment) = Benefit from negative exposure
 *
 *	FORMULAS:
 *
 *	FCA = integral_0^T EPE(t) * Q_c(t) * s_f(t) * DF(t) dt
 *	FBA = integral_0^T ENE(t) * Q_b(t) * s_b(t) * DF(t) dt
 *
 *	where:
 *	- s_f = funding spread (cost of borrowing)
 *	- s_b = lending spread (benefit of lending)
 *	- Q_c, Q_b = survival probabilities
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	FVA uncertainty comes from:
 *
 *	1. EXPOSURE UNCERTAINTY (EPE, ENE):
 *	   - Flows from the simulated exposure profile
 *	   - Higher volatility => higher exposures => higher FVA magnitude
 *
 *	2. FUNDING SPREAD UNCERTAINTY:
 *	   - Bank's funding costs vary with credit quality and market conditions
 *	   - Spread = borrowing rate - risk-free rate
 *
 *	3. SURVIVAL PROBABILITY UNCERTAINTY:
 *	   - FVA is weighted by survival (no funding needed post-default)
 *
 *	ACADEMIC REFERENCES:
 *
 *	[1] Gregory, J. (2020). "The xVA Challenge." Chapter 16: FVA.
 *	    - Comprehensive FVA framework
 *
 *	[2] Burgard, C. & Kjaer, M. (2013). "Funding Costs, Funding
 *	    Strategies." Risk Magazine.
 *	    - FVA derivation and debate
 *
 *	[3] Hull, J. & White, A. (2012). "The FVA Debate." Risk Magazine.
 *	    - Critical perspective on FVA
 *
 *	[4] Albanese, C., Andersen, L. & Iabichino, S. (2015). "FVA Accounting,
 *	    Risk Management, and Collateral Trading." Risk.
 *	    - Practical FVA implementation
 *
 *	[5] Crepey, S. (2015). "Bilateral Counterparty Risk under Funding
 *	    Constraints." Part I & II. Mathematical Finance.
 *	    - Rigorous mathematical treatment
 *
 *	PUBLIC IMPLEMENTATIONS REFERENCED:
 *
 *	- Open Source Risk Engine (ORE):
 *	  https://github.com/opensourcerisk/Engine
 *	  OREAnalytics/orea/aggregation/fvacalculator.cpp
 */

#include "xva-config.h"
#include <math.h>
#include <stdlib.h>

/*
 *	xvaComputeFca
 *
 *	Compute Funding Cost Adjustment.
 *
 *	FCA is the cost of funding POSITIVE exposures. When we have a
 *	positive MtM (counterparty owes us), we may need to borrow funds
 *	to finance our position until it pays off.
 *
 *	FORMULA:
 *
 *	FCA = integral_0^T EPE(t) * Q_c(t) * s_f * DF(t) dt
 *
 *	where s_f is the funding spread above the risk-free rate.
 *
 *	Reference: Gregory (2020), Equation (16.1)
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	FCA captures the COST OF BORROWING to fund positive exposures:
 *
 *	- EPE(t) determines HOW MUCH funding is needed
 *	- s_f determines the COST of that funding
 *	- Q_c(t) weights by counterparty survival (if they default,
 *	  we lose the exposure so funding becomes irrelevant)
 *
 *	Parameters:
 *	  profile           - Exposure profile with EPE
 *	  funding           - Funding data (spreads)
 *	  counterpartyCredit - Counterparty credit data (for survival probability)
 *	  market            - Market data
 *
 *	Returns: Funding Cost Adjustment (positive = cost)
 */
double
xvaComputeFca(
	const XvaExposureProfile *  profile,
	const XvaFundingData *      funding,
	const XvaCreditData *       counterpartyCredit,
	const XvaMarketData *       market)
{
	double  rr;
	double  maturity;
	double  dt;
	double  fundingSpread;
	double  lambdaC;
	double  fca;

	if (profile == NULL || funding == NULL || market == NULL)
	{
		return 0.0;
	}

	if (profile->expectedPositiveExposure == NULL)
	{
		return 0.0;
	}

	rr              = market->riskFreeRate;
	maturity        = market->maturity;
	dt              = maturity / (double) (profile->numTimeSteps - 1);
	fundingSpread   = funding->fundingSpread;

	/*
	 *	Counterparty hazard rate for survival weighting.
	 *	Uses actual credit data when available; falls back to funding spread.
	 */
	lambdaC = (counterpartyCredit != NULL) ? counterpartyCredit->hazardRate : fundingSpread;

	fca = 0.0;

	/*
	 *	Accumulate FCA using discrete summation
	 *
	 *	FCA = sum EPE(t_i) * Q_c(t_i) * s_f * DF(t_i) * dt
	 *
	 *	The product EPE * s_f represents the FUNDING COST per unit time
	 */
	for (size_t ii = 1; ii < profile->numTimeSteps; ii++)
	{
		double  tt;
		double  epe;
		double  qc;
		double  df;

		tt = (double) ii * dt;

		/*
		 *	Expected positive exposure
		 *	(From the exposure profile - captures MARKET UNCERTAINTY)
		 */
		epe = profile->expectedPositiveExposure[ii];

		/*
		 *	Survival probability (counterparty)
		 */
		qc = xvaSurvivalProbability(lambdaC, tt);

		/*
		 *	Discount factor
		 */
		df = xvaDiscountFactor(rr, tt);

		/*
		 *	FCA contribution: funding cost weighted by survival
		 */
		fca += epe * qc * fundingSpread * df * dt;
	}

	return fca;
}

/*
 *	xvaComputeFba
 *
 *	Compute Funding Benefit Adjustment.
 *
 *	FBA is the benefit from NEGATIVE exposures. When we have a
 *	negative MtM (we owe counterparty), we receive cash that can
 *	be invested or used to reduce borrowing.
 *
 *	FORMULA:
 *
 *	FBA = integral_0^T |ENE(t)| * Q_b(t) * s_b * DF(t) dt
 *
 *	where s_b is the lending/benefit spread.
 *
 *	Note: FBA is often assumed to be zero (s_b = 0) in practice
 *	because received cash may not be investable at a meaningful spread.
 *
 *	Reference: Gregory (2020), Equation (16.2)
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	FBA captures the BENEFIT of having negative exposure:
 *
 *	- |ENE(t)| determines the "free" funding received
 *	- s_b is typically < s_f (asymmetric funding)
 *	- Q_b(t) weights by own survival
 *
 *	Parameters:
 *	  profile   - Exposure profile with ENE
 *	  funding   - Funding data
 *	  ownCredit - Own bank credit data (for survival probability)
 *	  market    - Market data
 *
 *	Returns: Funding Benefit Adjustment (positive = benefit)
 */
double
xvaComputeFba(
	const XvaExposureProfile *  profile,
	const XvaFundingData *      funding,
	const XvaCreditData *       ownCredit,
	const XvaMarketData *       market)
{
	double  rr;
	double  maturity;
	double  dt;
	double  lendingSpread;
	double  lambdaB;
	double  fba;

	if (profile == NULL || funding == NULL || market == NULL)
	{
		return 0.0;
	}

	if (profile->expectedNegativeExposure == NULL)
	{
		return 0.0;
	}

	rr          = market->riskFreeRate;
	maturity    = market->maturity;
	dt          = maturity / (double) (profile->numTimeSteps - 1);

	/*
	 *	Lending spread (benefit rate)
	 *	Often set to zero in conservative calculations
	 */
	lendingSpread = funding->lendingRate - rr;

	if (lendingSpread < 0.0)
	{
		lendingSpread = 0.0;  /*	Can't have negative funding benefit */
	}

	/*
	 *	Own hazard rate for survival weighting.
	 *	Uses actual credit data when available; falls back to funding spread.
	 */
	lambdaB = (ownCredit != NULL) ? ownCredit->hazardRate : funding->fundingSpread;

	fba = 0.0;

	/*
	 *	Accumulate FBA using discrete summation
	 */
	for (size_t ii = 1; ii < profile->numTimeSteps; ii++)
	{
		double  tt;
		double  ene;
		double  qb;
		double  df;

		tt = (double) ii * dt;

		/*
		 *	Expected negative exposure (take absolute value)
		 */
		ene = -profile->expectedNegativeExposure[ii];

		if (ene < 0.0)
		{
			ene = 0.0;
		}

		/*
		 *	Own survival probability
		 */
		qb = xvaSurvivalProbability(lambdaB, tt);

		/*
		 *	Discount factor
		 */
		df = xvaDiscountFactor(rr, tt);

		/*
		 *	FBA contribution: funding benefit weighted by survival
		 */
		fba += ene * qb * lendingSpread * df * dt;
	}

	return fba;
}

/*
 *	xvaComputeFva
 *
 *	Compute total Funding Valuation Adjustment.
 *
 *	FORMULA:
 *	FVA = FCA - FBA
 *
 *	In most practical implementations, FBA is approximately 0 due to asymmetric
 *	funding (banks can borrow at spreads but not invest at spreads).
 *
 *	Reference: Gregory (2020), Equation (16.3)
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	FVA combines exposure and funding uncertainties:
 *
 *	1. Exposure uncertainty (from the exposure profile)
 *	2. Funding spread uncertainty (from bank's credit)
 *	3. Asymmetry: FCA typically dominates (borrowing > lending spread)
 *
 *	FVA can be SIGNIFICANT for large uncollateralized portfolios.
 *
 *	Parameters:
 *	  profile            - Exposure profile
 *	  funding            - Funding data
 *	  counterpartyCredit - Counterparty credit data (for FCA survival)
 *	  ownCredit          - Own bank credit data (for FBA survival)
 *	  market             - Market data
 *
 *	Returns: Total FVA (FCA - FBA)
 */
double
xvaComputeFva(
	const XvaExposureProfile *  profile,
	const XvaFundingData *      funding,
	const XvaCreditData *       counterpartyCredit,
	const XvaCreditData *       ownCredit,
	const XvaMarketData *       market)
{
	double  fca;
	double  fba;

	fca = xvaComputeFca(profile, funding, counterpartyCredit, market);
	fba = xvaComputeFba(profile, funding, ownCredit, market);

	/*
	 *	FVA = Funding Cost - Funding Benefit
	 *
	 *	Positive FVA = Net cost to us
	 *	Negative FVA = Net benefit (rare in practice)
	 */
	return fca - fba;
}
