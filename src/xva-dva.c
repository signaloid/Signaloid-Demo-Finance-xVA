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
 *	xva-dva.c
 *
 *	Debit Valuation Adjustment (DVA) Calculation Module
 *
 *	This module implements the calculation of Debit Valuation Adjustment
 *	and Bilateral CVA, which account for both parties' credit risk.
 *
 *	DVA DEFINITION:
 *	================================
 *	DVA is the market value of the bank's own default risk from the
 *	counterparty's perspective. It represents the expected gain to
 *	the bank from its own default.
 *
 *	DVA = (1 - R_bank) * integral_0^T ENE(t) * lambda_bank(t) * Q_bank(t) * DF(t) dt
 *
 *	where:
 *	- R_bank = bank's recovery rate
 *	- ENE(t) = expected negative exposure at time t
 *	- lambda_bank(t) = bank's hazard rate
 *	- Q_bank(t) = bank's survival probability
 *
 *	DVA is controversial because it implies a bank benefits from its
 *	own deteriorating credit quality.
 *
 *	BILATERAL CVA (BCVA):
 *	================================
 *	BCVA accounts for both parties' credit risk:
 *
 *	BCVA = CVA - DVA
 *
 *	This represents the symmetric treatment of counterparty risk
 *	where both parties face potential losses.
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	DVA uncertainty mirrors CVA but focuses on NEGATIVE exposures:
 *
 *	1. ENE(t) - Market risk (what counterparty could lose)
 *	   - ENE is the mirror image of EPE
 *	   - For a two-way trade: ENE_us = -EPE_counterparty
 *
 *	2. lambda_bank - Bank's own default probability
 *	   - Derived from bank's credit spread
 *	   - Higher spread => higher DVA (counterintuitive "benefit")
 *
 *	ACADEMIC REFERENCES:
 *
 *	[1] Gregory, J. (2020). "The xVA Challenge." Chapter 14: DVA.
 *	    - Comprehensive DVA framework
 *
 *	[2] Burgard, C. & Kjaer, M. (2011). "Partial Differential Equation
 *	    Representations of Derivatives with Bilateral Counterparty Risk
 *	    and Funding Costs." Journal of Credit Risk, 7(3).
 *	    - PDE approach to bilateral valuation
 *
 *	[3] Brigo, D., Morini, M. & Pallavicini, A. (2013). "Counterparty
 *	    Credit Risk, Collateral and Funding." Wiley.
 *	    - Rigorous bilateral framework
 *
 *	[4] Hull, J. & White, A. (2014). "Valuing Derivatives: Funding
 *	    Value Adjustments and Fair Value." Financial Analysts Journal.
 *	    - Industry perspective on DVA
 *
 *	PUBLIC IMPLEMENTATIONS REFERENCED:
 *
 *	- Open Source Risk Engine (ORE):
 *	  https://github.com/opensourcerisk/Engine
 *
 *	- QuantLib CVA/DVA calculators
 */

#include "xva-config.h"
#include <math.h>
#include <stdlib.h>

/*
 *	xvaComputeDva
 *
 *	Compute Debit Valuation Adjustment using the unilateral formula.
 *
 *	FORMULA (Discrete approximation):
 *
 *	DVA = (1 - R_own) * Sigma_i |ENE(t_i)| * DF(t_i) * PD_own(t_{i-1}, t_i)
 *
 *	Note: ENE is typically negative, so we use |ENE| to get positive DVA.
 *
 *	Reference: Gregory (2020), Equation (14.5)
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	DVA captures the expected benefit from the bank's own default:
 *
 *	1. |ENE(t_i)| - The exposure that counterparty has to us
 *	   - Read from the exposure profile (Monte Carlo ensemble or UxHw path)
 *	   - ENE = E[min(V, 0)] where V is our mark-to-market
 *
 *	2. PD_own - Our probability of default
 *	   - Higher own credit spread => higher DVA
 *	   - This is the controversial aspect of DVA
 *
 *	Parameters:
 *	  profile    - Exposure profile with ENE at each time step
 *	  ownCredit  - Bank's own credit data
 *	  market     - Market data
 *
 *	Returns: Computed DVA value (positive = benefit to us)
 */
double
xvaComputeDva(
	const XvaExposureProfile *  profile,
	const XvaCreditData *       ownCredit,
	const XvaMarketData *       market)
{
	double  riskFreeRate;
	double  lambda;
	double  recovery;
	double  maturity;
	double  dt;
	double  lgd;
	double  dva;

	if (profile == NULL || ownCredit == NULL || market == NULL)
	{
		return 0.0;
	}

	if (profile->expectedNegativeExposure == NULL)
	{
		return 0.0;
	}

	riskFreeRate    = market->riskFreeRate;
	lambda          = ownCredit->hazardRate;
	recovery        = ownCredit->recoveryRate;
	maturity        = market->maturity;
	dt              = maturity / (double) (profile->numTimeSteps - 1);

	/*
	 *	Loss given default (from counterparty's perspective)
	 */
	lgd = 1.0 - recovery;

	/*
	 *	Accumulate DVA using discrete summation
	 *
	 *	DVA = LGD * Sigma |ENE(t_i)| * DF(t_i) * PD_own(t_{i-1}, t_i)
	 */
	dva = 0.0;

	for (size_t ii = 1; ii < profile->numTimeSteps; ii++)
	{
		double  tPrev   = (double) (ii - 1) * dt;
		double  tCurr   = (double) ii * dt;
		double  ene;
		double  df;
		double  pd;

		/*
		 *	Expected negative exposure at current time
		 *	ENE is negative, so we negate to get positive contribution
		 *
		 *	This represents WHAT THE COUNTERPARTY COULD LOSE if we default
		 */
		ene = -profile->expectedNegativeExposure[ii];  /*	Negate */

		if (ene < 0.0)
		{
			ene = 0.0;  /*	Safety check */
		}

		/*
		 *	Discount factor
		 */
		df = xvaDiscountFactor(riskFreeRate, tCurr);

		/*
		 *	Our probability of default in current interval
		 *	(Captures OUR OWN CREDIT UNCERTAINTY)
		 */
		pd = xvaDefaultProbability(lambda, tPrev, tCurr);

		/*
		 *	Contribution to DVA
		 */
		dva += ene * df * pd;
	}

	return lgd * dva;
}

/*
 *	xvaComputeBilateralCva
 *
 *	Compute Bilateral CVA (BCVA) accounting for both parties' credit risk.
 *
 *	FORMULA:
 *	BCVA = CVA - DVA
 *
 *	This is the net credit adjustment considering that:
 *	- We may lose if counterparty defaults (CVA, positive)
 *	- We may "gain" if we default (DVA, subtracted)
 *
 *	Reference: Gregory (2020), Equation (14.7)
 *
 *	ACCOUNTING NOTE:
 *	Under IFRS 13, both CVA and DVA should be recognized, though
 *	DVA creates paradoxical P&L (profit when credit deteriorates).
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	BCVA combines uncertainties from both credit perspectives:
 *
 *	- Counterparty credit risk (CVA component)
 *	- Own credit risk (DVA component)
 *	- These may be correlated (systemic risk)
 *
 *	Parameters:
 *	  profile            - Exposure profile
 *	  counterpartyCredit - Counterparty credit data
 *	  ownCredit          - Bank's own credit data
 *	  market             - Market data
 *
 *	Returns: Bilateral CVA (CVA - DVA)
 */
double
xvaComputeBilateralCva(
	const XvaExposureProfile *  profile,
	const XvaCreditData *       counterpartyCredit,
	const XvaCreditData *       ownCredit,
	const XvaMarketData *       market)
{
	double  cva;
	double  dva;

	/*
	 *	Import CVA computation
	 */
	extern double
	xvaComputeCva(
		const XvaExposureProfile *,
		const XvaCreditData *,
		const XvaMarketData *);

	if (profile == NULL || counterpartyCredit == NULL ||
	    ownCredit == NULL || market == NULL)
	{
		return 0.0;
	}

	/*
	 *	Compute unilateral CVA (counterparty risk to us)
	 */
	cva = xvaComputeCva(profile, counterpartyCredit, market);

	/*
	 *	Compute DVA (our risk to counterparty)
	 */
	dva = xvaComputeDva(profile, ownCredit, market);

	/*
	 *	Bilateral CVA = CVA - DVA
	 *
	 *	CVA is a cost (positive)
	 *	DVA is a benefit (subtracted)
	 */
	return cva - dva;
}
