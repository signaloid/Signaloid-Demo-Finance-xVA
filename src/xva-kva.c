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
 *	xva-kva.c
 *
 *	Capital Valuation Adjustment (KVA) Calculation Module
 *
 *	This module implements the calculation of Capital Valuation Adjustment,
 *	which captures the cost of holding regulatory capital against
 *	counterparty credit risk.
 *
 *	KVA DEFINITION:
 *	================================
 *	KVA is the cost of capital that must be held against counterparty
 *	credit risk exposure over the life of the trade.
 *
 *	KVA = integral_0^T RC(t) * h * Q(t) * DF(t) dt
 *
 *	where:
 *	- RC(t) = regulatory capital requirement at time t
 *	- h = hurdle rate (cost of equity)
 *	- Q(t) = survival probability
 *	- DF(t) = discount factor
 *
 *	SIMPLIFICATIONS:
 *	================================
 *	This implementation uses a simplified "unilateral" KVA model:
 *
 *	1. Single survival probability with a fixed hazard rate:
 *	   The full Green-Kenyon-Dennis (2014) formula uses bilateral survival:
 *	   e^(-integral(r + lambda_B + lambda_C)ds) where lambda_B is the bank's
 *	   own hazard rate and lambda_C is the counterparty's. This implementation
 *	   uses a single survival factor with the fixed hazard rate
 *	   kXvaUnilateralHazardRate. It does not take the counterparty's
 *	   calibrated hazard rate as an input.
 *
 *	2. CCR capital only:
 *	   Green-Kenyon-Dennis (2014) Equation 34 decomposes capital as:
 *	   K = K_MR + K_CCR + K_CVA (market risk + CCR + CVA risk capital).
 *	   This implementation computes only the CCR component.
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	================================
 *	KVA uncertainty comes from:
 *
 *	1. EXPOSURE UNCERTAINTY:
 *	   - EAD is derived from the EPE profile
 *	   - Higher volatility => higher EPE => higher capital => higher KVA
 *
 *	2. CREDIT RISK PARAMETERS:
 *	   - PD, LGD uncertainty affects capital calculation
 *
 *	3. HURDLE RATE UNCERTAINTY:
 *	   - Cost of equity varies by institution (typically 8-15%)
 *
 *	ACADEMIC REFERENCES:
 *
 *	[1] Gregory, J. (2020). "The xVA Challenge." Chapter 18: KVA.
 *	[2] BCBS (2017). "Basel III: Finalising Post-Crisis Reforms."
 *	[3] Albanese, C., Caenazzo, S. & Crepey, S. (2017). "Capital
 *	    Valuation Adjustment and Funding Valuation Adjustment."
 *	[4] Green, A., Kenyon, C. & Dennis, C.R. (2014). "KVA: Capital
 *	    Valuation Adjustment by Replication." Risk Magazine.
 *
 *	PUBLIC IMPLEMENTATIONS REFERENCED:
 *
 *	- Open Source Risk Engine (ORE): https://github.com/opensourcerisk/Engine
 */

#include "xva-config.h"
#include <math.h>
#include <stdlib.h>

/*
 *	xvaComputeKva
 *
 *	Compute Capital Valuation Adjustment.
 *
 *	FORMULA:
 *	KVA = integral_0^T RC(t) * h * Q(t) * DF(t) dt
 *
 *	Reference: Gregory (2020), Equation (18.1)
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	KVA captures the COST OF CAPITAL UNCERTAINTY:
 *	- RC(t) depends on the EPE profile
 *	- h (hurdle rate) varies by institution
 *	- Q(t) captures credit survival uncertainty
 */
double
xvaComputeKva(
	const XvaExposureProfile *  profile,
	const XvaCapitalData *      capital,
	const XvaMarketData *       market)
{
	double  riskFreeRate;
	double  maturity;
	size_t  numSteps;
	double  dt;
	double  hurdleRate;
	double  rwFactor;
	double  capitalRatio;
	/*
	 *	Unilateral KVA with the fixed hazard rate kXvaUnilateralHazardRate
	 *	standing in for the counterparty's. For bilateral KVA per
	 *	Green-Kenyon-Dennis (2014), the bank's own hazard rate lambda_B
	 *	would also enter the discount factor.
	 */
	double  lambda      = kXvaUnilateralHazardRate;
	double  runningMax  = 0.0;
	double  kva         = 0.0;

	if (profile == NULL || capital == NULL || market == NULL)
	{
		return 0.0;
	}

	riskFreeRate    = market->riskFreeRate;
	maturity        = market->maturity;
	numSteps        = profile->numTimeSteps;
	dt              = maturity / (double) (numSteps - 1);

	hurdleRate = capital->hurdleRate;

	if (hurdleRate <= 0.0)
	{
		hurdleRate = (double) kXvaHurdleRatePercent / 100.0;
	}

	/*
	 *	rwaFactor is the risk weight w from standardised approach.
	 *	Default to 100% (BB-rated counterparty) if not specified.
	 */
	rwFactor        = capital->rwaFactor > 0 ? capital->rwaFactor : kXvaDefaultRiskWeight;
	capitalRatio    = capital->capitalRatio > 0 ?
	                  capital->capitalRatio :
	                  (double) kXvaCapitalRatioPercent / 100.0;

	for (size_t ii = 1; ii < numSteps; ii++)
	{
		double  time    = (double) ii * dt;
		double  epe     = profile->expectedPositiveExposure[ii];
		double  eadAtTime;
		double  rwa;
		double  regulatoryCapital;
		double  survivalProb;
		double  discountFactor;

		if (epe > runningMax)
		{
			runningMax = epe;
		}

		/*
		 *	Basel standardised: RWA = w × 12.5 × EAD
		 *	Capital = capitalRatio × RWA
		 *	(Green-Kenyon-Dennis 2014, Equation 35)
		 */
		eadAtTime           = kXvaBaselAlpha * runningMax;
		rwa                 = rwFactor * kXvaBaselRwaMultiplier * eadAtTime;
		regulatoryCapital   = rwa * capitalRatio;
		survivalProb        = xvaSurvivalProbability(lambda, time);
		discountFactor      = xvaDiscountFactor(riskFreeRate, time);

		kva += regulatoryCapital * hurdleRate * survivalProb * discountFactor * dt;
	}

	return kva;
}
