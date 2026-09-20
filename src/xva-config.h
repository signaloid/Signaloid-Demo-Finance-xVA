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
 *	xva-config.h
 *
 *	Header file for the xVA simulation framework
 *
 *	This header defines the configuration parameters, data structures,
 *	and function declarations for computing various valuation adjustments
 *	(xVA) from a simulated exposure profile, whether that profile comes
 *	from a Monte Carlo path ensemble or from a single UxHw distributional
 *	path.
 *
 *	ACADEMIC REFERENCES:
 *
 *	[1] Glasserman, P. (2004). "Monte Carlo Methods in Financial Engineering."
 *	    Springer, New York. ISBN: 978-0387004518
 *	    - Foundation for Monte Carlo simulation in finance
 *	    - Chapters 1-4: Random number generation, variance reduction
 *
 *	[2] Gregory, J. (2020). "The xVA Challenge: Counterparty Risk, Funding,
 *	    Collateral, Capital and Initial Margin." 4th Edition, Wiley.
 *	    ISBN: 978-1119508977
 *	    - Comprehensive xVA framework and formulas
 *
 *	[3] Crepey, S. (2014). "XVA: About CVA, DVA, FVA and Other Market
 *	    Adjustments." Opinion and Debates, No. 5.
 *	    - Mathematical formulation of interdependent xVA terms
 *
 *	[4] Brigo, D., Morini, M., & Pallavicini, A. (2013). "Counterparty Credit
 *	    Risk, Collateral and Funding." Wiley.
 *	    - Advanced credit risk modeling
 *
 *	[5] Abbas-Turki, L.A., Crepey, S., & Saadeddine, B. (2023). "Pathwise
 *	    CVA Regressions with Oversimulated Defaults."
 *	    Mathematical Finance, 33(2), pp. 274-307.
 *	    - Modern neural network approaches to xVA
 *
 *	[6] Hull, J.C. (2017). "Options, Futures, and Other Derivatives."
 *	    10th Edition, Pearson. ISBN: 978-0134472089
 *	    - Standard derivative pricing and risk management
 *
 *	PUBLIC IMPLEMENTATIONS REFERENCED:
 *
 *	- Open Source Risk Engine (ORE): https://github.com/opensourcerisk/Engine
 *	- QuantLib: https://www.quantlib.org/
 *	- DROP (Lakshmi DRIP): https://github.com/lakshmiDRIP/DROP
 *	- xVA R Package: https://github.com/xvacode/xVA
 *	- NeuralXVA: https://github.com/BouazzaSE/NeuralXVA
 *	- Monte Carlo Risk Engine (PyTorch):
 *	  https://github.com/konstantineder/montecarlo-risk-engine
 */

#ifndef XVA_CONFIG_H
#define XVA_CONFIG_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <uxhw.h>

/*
 *	Configuration enumeration defining all integer constants used throughout
 *	the xVA Monte Carlo framework. Using enum avoids magic numbers
 *	and provides centralized configuration management.
 *
 *	Naming convention: kXva prefix for all configuration values
 */
typedef enum
{
	/*
	 *	Simulation time grid
	 *	Reference: Glasserman (2004), Chapter 1
	 */
	kXvaNumTimeSteps        = 252,      /*	Time steps (trading days/year) */
	kXvaTradingDaysPerYear  = 252,      /*	Business days per year */

	/*
	 *	xVA calculation parameters
	 *	Reference: Gregory (2020), Chapters 14-19
	 */
	kXvaCapitalRatioPercent = 8,        /*	Regulatory capital ratio */
	kXvaHurdleRatePercent   = 10,       /*	Cost of capital hurdle rate */
	kXvaImPeriodDays        = 10        /*	Initial margin period of risk */
} XvaConfig;

/*
 *	Index of each output the demo can compute. `-S` selects one of these,
 *	and `kXvaOutputIndexMax` selects all of them. The order must match the
 *	`outputVariables[...]` entries in `signaloid.yaml`.
 */
typedef enum
{
	kXvaOutputIndexCva = 0,
	kXvaOutputIndexDva,
	kXvaOutputIndexBilateralCva,
	kXvaOutputIndexFva,
	kXvaOutputIndexMva,
	kXvaOutputIndexKva,
	kXvaOutputIndexTotalXva,
	kXvaOutputIndexEuropeanPut,
	kXvaOutputIndexMax
} XvaOutputIndex;

/*
 *	Index of each uncertain input the demo accepts. The demo takes the
 *	volatility alone, either from `-f` or from the input CSV read with `-i`.
 */
typedef enum
{
	kXvaInputIndexVolatility = 0,
	kXvaInputIndexMax
} XvaInputIndex;

/*
 *	Default values for the market, credit and contract parameters the demo
 *	exposes on the command line. These are floating point and so cannot live
 *	in `XvaConfig` above.
 */
#define kXvaDemoDefaultSpotPrice        (100.0)
#define kXvaDemoDefaultRiskFreeRate     (0.05)
#define kXvaDemoDefaultDividendYield    (0.02)
#define kXvaDemoDefaultVolatility       (0.20)
#define kXvaDemoDefaultMaturity         (5.0)
#define kXvaDemoDefaultCreditSpread     (0.015)
#define kXvaDemoDefaultOwnCreditSpread  (0.01)
#define kXvaDemoDefaultStrikePrice      (100.0)
#define kXvaDemoDefaultRecoveryRate     (0.40)
#define kXvaDemoDefaultFundingSpread    (0.01)

/*
 *	Floating-point constants (cannot be in enum)
 *	Naming convention: kXva prefix
 */

/*
 *	PFE quantile level (97.5th percentile)
 *	Reference: Gregory (2020), Section 10.2
 */
#define kXvaPfeQuantile 0.975

/*
 *	Basel III regulatory constants
 *	Reference: BCBS (2017), Green-Kenyon-Dennis (2014)
 */
#define kXvaBaselAlpha          1.4             /*	EAD multiplier */
#define kXvaBaselRwaMultiplier  12.5            /*	RWA = w * 12.5 * EAD */
#define kXvaDefaultRiskWeight   1.0             /*	100% for BB-rated */

/*
 *	Unilateral KVA default hazard rate
 *	Reference: Gregory (2020), Chapter 18
 */
#define kXvaUnilateralHazardRate 0.02

/*
 *	Market data structure containing all necessary inputs for simulation
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	The volatility (sigma) scales the Brownian increments that drive the
 *	diffusion. In this demo it is a fixed calibrated input, and the
 *	uncertainty comes from the increments themselves.
 *	Reference: Glasserman (2004), Section 3.1
 */
typedef struct
{
	double  spotPrice;          /*	Current asset price S(0) */
	double  riskFreeRate;       /*	Continuous risk-free rate r */
	double  dividendYield;      /*	Continuous dividend yield q */
	double  volatility;         /*	Annualized volatility sigma */
	double  maturity;           /*	Time to maturity T in years */
} XvaMarketData;

/*
 *	Counterparty credit data for CVA/DVA calculations
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	Hazard rates (derived from credit spreads) introduce credit
 *	uncertainty into the xVA calculations.
 *	Reference: Gregory (2020), Chapter 8
 */
typedef struct
{
	double  creditSpread;       /*	Credit spread in decimal (e.g., 0.01 = 100bps) */
	double  recoveryRate;       /*	Expected recovery rate (0 to 1) */
	double  hazardRate;         /*	Instantaneous default intensity lambda */
} XvaCreditData;

/*
 *	Funding parameters for FVA calculations
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	Funding spreads vary with market conditions and bank credit quality.
 *	Reference: Gregory (2020), Chapter 16
 */
typedef struct
{
	double  fundingSpread;      /*	Spread over risk-free rate (decimal) */
	double  borrowingRate;      /*	Rate for borrowing funds */
	double  lendingRate;        /*	Rate for lending/investing funds */
} XvaFundingData;

/*
 *	Capital requirements for KVA calculations
 *
 *	Reference: Gregory (2020), Chapter 18
 *	Basel III/IV regulatory framework
 */
typedef struct
{
	double  capitalRatio;       /*	Required capital ratio */
	double  hurdleRate;         /*	Cost of capital hurdle rate */
	double  rwaFactor;          /*	Risk-weighted asset factor */
} XvaCapitalData;

/*
 *	Exposure profile containing simulated exposure metrics
 *
 *	UNCERTAINTY QUANTIFICATION:
 *	EPE and PFE are computed from the distribution of simulated exposures,
 *	directly capturing market risk uncertainty.
 *	Reference: Gregory (2020), Chapter 10
 */
typedef struct
{
	double *    expectedExposure;           /*	EE[t] for each time step */
	double *    expectedPositiveExposure;   /*	EPE[t] = E[max(V,0)] */
	double *    expectedNegativeExposure;   /*	ENE[t] = E[min(V,0)] */
	double *    potentialFutureExposure;    /*	PFE[t] at specified quantile */
	size_t      numTimeSteps;               /*	Number of time steps */
} XvaExposureProfile;

/*
 *	Simulation paths storage structure
 */
typedef struct
{
	double * *  paths;          /*	paths[pathIndex][timeStep] */
	size_t      numPaths;       /*	Number of simulation paths */
	size_t      numSteps;       /*	Number of time steps per path */
	double      dt;             /*	Time step size */
} XvaSimulationPaths;

/*
 *	Per-step constants of the exact GBM solution, shared by the Monte Carlo
 *	and UxHw path simulators.
 */
typedef struct
{
	double  driftDt;            /*	(r - q - sigma^2/2) * dt */
	double  diffusionSqrtDt;    /*	sigma * sqrt(dt) */
	size_t  numSteps;           /*	Number of steps, excluding t=0 */
} XvaGbmPathParameters;

/* ========================================================================
 *	FUNCTION DECLARATIONS
 * ======================================================================== */

/*
 *	Geometric Brownian Motion Simulation (xva-gbm.c)
 *	Reference: Glasserman (2004), Chapter 3
 *
 *	Stochasticity is supplied directly by the UxHw distributional
 *	primitives (e.g. UxHwDoubleGaussDist). In native Monte Carlo mode
 *	the compatibility layer draws a fresh pseudo-random sample on each
 *	call; on the Signaloid Cloud Compute Engine the same call yields a
 *	distributional value carried through the computation.
 */
int
xvaAllocateGbmPaths(
	XvaSimulationPaths *    paths,
	const XvaMarketData *   market,
	size_t                  numPaths,
	XvaGbmPathParameters *  parameters);

void
xvaFreeSimulationPaths(XvaSimulationPaths * paths);

/*
 *	Exposure Calculations (xva-exposure.c)
 *	Reference: Gregory (2020), Chapter 10
 */
int
xvaAllocateExposureProfile(
	XvaExposureProfile *    profile,
	size_t                  numSteps);

double
xvaComputeForwardMtm(
	double  spot,
	double  strike,
	double  riskFreeRate,
	double  dividendYield,
	double  timeToMaturity);

void
xvaFreeExposureProfile(XvaExposureProfile * profile);

/*
 *	Credit Valuation Adjustment (xva-cva.c)
 *	Reference: Gregory (2020), Chapter 14
 */
double
xvaComputeCva(
	const XvaExposureProfile *  profile,
	const XvaCreditData *       credit,
	const XvaMarketData *       market);

/*
 *	Debit Valuation Adjustment (xva-dva.c)
 *	Reference: Gregory (2020), Chapter 14
 */
double
xvaComputeDva(
	const XvaExposureProfile *  profile,
	const XvaCreditData *       ownCredit,
	const XvaMarketData *       market);

double
xvaComputeBilateralCva(
	const XvaExposureProfile *  profile,
	const XvaCreditData *       counterpartyCredit,
	const XvaCreditData *       ownCredit,
	const XvaMarketData *       market);

/*
 *	Funding Valuation Adjustment (xva-fva.c)
 *	Reference: Gregory (2020), Chapter 16
 */
double
xvaComputeFva(
	const XvaExposureProfile *  profile,
	const XvaFundingData *      funding,
	const XvaCreditData *       counterpartyCredit,
	const XvaCreditData *       ownCredit,
	const XvaMarketData *       market);

double
xvaComputeFca(
	const XvaExposureProfile *  profile,
	const XvaFundingData *      funding,
	const XvaCreditData *       counterpartyCredit,
	const XvaMarketData *       market);

double
xvaComputeFba(
	const XvaExposureProfile *  profile,
	const XvaFundingData *      funding,
	const XvaCreditData *       ownCredit,
	const XvaMarketData *       market);

/*
 *	Margin Valuation Adjustment (xva-mva.c)
 *	Reference: Gregory (2020), Chapter 17
 */
double
xvaComputeMva(
	const XvaExposureProfile *  profile,
	const XvaFundingData *      funding,
	const XvaCapitalData *      capital,
	const XvaCreditData *       credit,
	const XvaMarketData *       market);

/*
 *	Capital Valuation Adjustment (xva-kva.c)
 *	Reference: Gregory (2020), Chapter 18
 */
double
xvaComputeKva(
	const XvaExposureProfile *  profile,
	const XvaCapitalData *      capital,
	const XvaMarketData *       market);

/*
 *	Utility functions (xva-utils.c)
 */
double
xvaDiscountFactor(double rate, double time);

double
xvaSurvivalProbability(double hazardRate, double time);

double
xvaDefaultProbability(double hazardRate, double t1, double t2);

int
xvaCompareDoubles(const void * aa, const void *  bb);

double
xvaQuantile(double * data, size_t nn, double qq);

double
xvaStandardNormalCdf(double x);

double
xvaBlackScholesPut(
	double  spotPrice,
	double  strikePrice,
	double  riskFreeRate,
	double  dividendYield,
	double  volatility,
	double  timeToMaturity);

double
xvaHazardRateFromSpread(double spread, double recoveryRate);

#endif /* XVA_CONFIG_H */
