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

#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include <assert.h>
#include "utilities.h"
#include "common.h"

void
printUsage(void)
{
	fprintf(stderr, "xVA Monte Carlo Simulation with Signaloid UxHw Support\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Computes valuation adjustments (xVA) for counterparty credit risk.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Output indices (-S flag):\n");
	fprintf(stderr, "  0: CVA  (Credit Valuation Adjustment)\n");
	fprintf(stderr, "  1: DVA  (Debit Valuation Adjustment)\n");
	fprintf(stderr, "  2: BCVA (Bilateral CVA = CVA - DVA)\n");
	fprintf(stderr, "  3: FVA  (Funding Valuation Adjustment)\n");
	fprintf(stderr, "  4: MVA  (Margin Valuation Adjustment)\n");
	fprintf(stderr, "  5: KVA  (Capital Valuation Adjustment)\n");
	fprintf(stderr, "  6: Total xVA (CVA - DVA + FVA + MVA + KVA)\n");
	fprintf(stderr, "  7: European Put (Black-Scholes)\n");
	fprintf(stderr, "\n");
	printCommonUsage();
	fprintf(
		stderr,
		"	[-f, --volatility <annualized volatility: double (Default: %lf)>]\n"
		"	[-s, --spot-price <spot price: double (Default: %lf)>]\n"
		"	[-r, --risk-free-rate <risk-free rate: double (Default: %lf)>]\n"
		"	[-d, --dividend-yield <dividend yield: double (Default: %lf)>]\n"
		"	[-m, --maturity <time to maturity in years: double (Default: %lf)>]\n"
		"	[-c, --credit-spread <counterparty credit spread: double (Default: %lf)>]\n"
		"	[-e, --own-credit-spread <own credit spread: double (Default: %lf)>]\n"
		"	[-k, --strike-price <strike price: double (Default: %lf)>]\n",
		kXvaDemoDefaultVolatility,
		kXvaDemoDefaultSpotPrice,
		kXvaDemoDefaultRiskFreeRate,
		kXvaDemoDefaultDividendYield,
		kXvaDemoDefaultMaturity,
		kXvaDemoDefaultCreditSpread,
		kXvaDemoDefaultOwnCreditSpread,
		kXvaDemoDefaultStrikePrice
	);
	fprintf(stderr, "\n");

	return;
}

/**
 *	@brief	Set the default values for the command line arguments.
 *
 *	@param	arguments	: command line arguments pointer.
 *	@return			: `kCommonConstantReturnTypeSuccess` if successful, else `kCommonConstantReturnTypeError`.
 */
CommonConstantReturnType
setDefaultCommandLineArguments(CommandLineArguments * arguments)
{
	if (arguments == NULL)
	{
		fprintf(stderr, "Error: The provided pointer to arguments is NULL.\n");

		return kCommonConstantReturnTypeError;
	}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
	*arguments = (CommandLineArguments) {
		.common             = (CommonCommandLineArguments) { 0 },
		.kernelParameters   = (KernelParameters) {
			.spotPrice          = kXvaDemoDefaultSpotPrice,
			.riskFreeRate       = kXvaDemoDefaultRiskFreeRate,
			.dividendYield      = kXvaDemoDefaultDividendYield,
			.maturity           = kXvaDemoDefaultMaturity,
			.creditSpread       = kXvaDemoDefaultCreditSpread,
			.ownCreditSpread    = kXvaDemoDefaultOwnCreditSpread,
			.strikePrice        = kXvaDemoDefaultStrikePrice,
			.inputVariables     = {
				[kXvaInputIndexVolatility] = kXvaDemoDefaultVolatility,
			},
		},
		.volatility         = kXvaDemoDefaultVolatility,
	};
#pragma GCC diagnostic pop

	return kCommonConstantReturnTypeSuccess;
}

#ifdef NO_OS_AVAILABLE
CommonConstantReturnType
setNoOSCommandLineArguments(CommandLineArguments * arguments)
{
	if (arguments == NULL)
	{
		fprintf(stderr, "Error: The provided pointer to arguments is NULL.\n");

		return kCommonConstantReturnTypeError;
	}

	/*
	 *	Start from the defaults so that every demo-specific field is
	 *	initialized, then override the ones the no-OS build fixes. The
	 *	`common` sub-struct is zeroed by this call, so it is set explicitly
	 *	below.
	 */
	if (setDefaultCommandLineArguments(arguments) != kCommonConstantReturnTypeSuccess)
	{
		return kCommonConstantReturnTypeError;
	}

	arguments->common.numberOfMonteCarloIterations  = 1;
	arguments->common.outputSelect                  = kXvaOutputIndexMax;
	arguments->common.isTimingEnabled               = false;
	arguments->common.isMonteCarloMode              = false;
	arguments->common.isOutputJSONMode              = false;
	arguments->common.isWriteToFileEnabled          = false;

	return kCommonConstantReturnTypeSuccess;
}
#endif

CommonConstantReturnType
getCommandLineArguments(int argc, char *  argv[], CommandLineArguments *  arguments)
{
#ifdef NO_OS_AVAILABLE
	/*
	 *	The no-OS build has no command line to parse, so ignore `argc` and
	 *	`argv` and use the hard-coded configuration instead.
	 */
	(void) argc;
	(void) argv;

	puts("Using hard coded command line arguments");

	return setNoOSCommandLineArguments(arguments);

#else
	const char *    volatilityArg       = NULL;
	const char *    spotPriceArg        = NULL;
	const char *    riskFreeRateArg     = NULL;
	const char *    dividendYieldArg    = NULL;
	const char *    maturityArg         = NULL;
	const char *    creditSpreadArg     = NULL;
	const char *    ownCreditSpreadArg  = NULL;
	const char *    strikePriceArg      = NULL;

	if (arguments == NULL)
	{
		fprintf(stderr, "Error: The provided pointer to arguments is NULL.\n");

		return kCommonConstantReturnTypeError;
	}

	if (setDefaultCommandLineArguments(arguments) != kCommonConstantReturnTypeSuccess)
	{
		return kCommonConstantReturnTypeError;
	}

	DemoOption options[] = {
		{ .opt = "f", .optAlternative = "volatility",        .hasArg = true, .foundArg = &volatilityArg,      .foundOpt = NULL },
		{ .opt = "s", .optAlternative = "spot-price",        .hasArg = true, .foundArg = &spotPriceArg,       .foundOpt = NULL },
		{ .opt = "r", .optAlternative = "risk-free-rate",    .hasArg = true, .foundArg = &riskFreeRateArg,    .foundOpt = NULL },
		{ .opt = "d", .optAlternative = "dividend-yield",    .hasArg = true, .foundArg = &dividendYieldArg,   .foundOpt = NULL },
		{ .opt = "m", .optAlternative = "maturity",          .hasArg = true, .foundArg = &maturityArg,        .foundOpt = NULL },
		{ .opt = "c", .optAlternative = "credit-spread",     .hasArg = true, .foundArg = &creditSpreadArg,    .foundOpt = NULL },
		{ .opt = "e", .optAlternative = "own-credit-spread", .hasArg = true, .foundArg = &ownCreditSpreadArg, .foundOpt = NULL },
		{ .opt = "k", .optAlternative = "strike-price",      .hasArg = true, .foundArg = &strikePriceArg,     .foundOpt = NULL },
		{ 0 },
	};

	if (parseArgs(argc, argv, &arguments->common, options) != kCommonConstantReturnTypeSuccess)
	{
		fprintf(stderr, "Error: Parsing command line arguments failed.\n");
		printUsage();

		return kCommonConstantReturnTypeError;
	}

	if (arguments->common.isHelpEnabled)
	{
		printUsage();

		exit(EXIT_SUCCESS);
	}

	/*
	 *	If no output is selected, set `outputSelect` to `kXvaOutputIndexMax`.
	 *	This triggers the demo to compute all outputs.
	 */
	if (!arguments->common.isOutputSelected)
	{
		arguments->common.outputSelect = kXvaOutputIndexMax;
	}

	/*
	 *	When `outputSelect` is set to `kXvaOutputIndexMax`, we cannot be
	 *	in Monte Carlo mode.
	 */
	if (arguments->common.outputSelect == kXvaOutputIndexMax)
	{
		if (arguments->common.isMonteCarloMode)
		{
			fprintf(stderr, "Error: Please select a single output when in Monte Carlo mode.\n");

			return kCommonConstantReturnTypeError;
		}
	}
	else if (arguments->common.outputSelect > kXvaOutputIndexMax)
	{
		fprintf(stderr, "Error: Wrong output selection.\n");

		return kCommonConstantReturnTypeError;
	}

	/*
	 *	Monte Carlo mode does not support input from file.
	 */
	if ((arguments->common.isMonteCarloMode) && (arguments->common.isInputFromFileEnabled))
	{
		fprintf(stderr, "Error: Monte Carlo mode does not support input from file.\n");

		return kCommonConstantReturnTypeError;
	}

	if (arguments->common.isVerbose)
	{
		fprintf(stderr, "Warning: Verbose mode not supported. Continuing in non-verbose mode.\n");
	}

	if (volatilityArg != NULL)
	{
		double volatility;

		if (parseDoubleChecked(volatilityArg, &volatility) != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The volatility must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		if (volatility <= 0.0)
		{
			fprintf(stderr, "Error: The volatility must be positive.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->volatility = volatility;
		arguments->kernelParameters.inputVariables[kXvaInputIndexVolatility] = volatility;
	}

	if (spotPriceArg != NULL)
	{
		double spotPrice;

		if (parseDoubleChecked(spotPriceArg, &spotPrice) != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The spot price must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		if (spotPrice <= 0.0)
		{
			fprintf(stderr, "Error: The spot price must be positive.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->kernelParameters.spotPrice = spotPrice;
	}

	if (riskFreeRateArg != NULL)
	{
		double riskFreeRate;

		if (parseDoubleChecked(riskFreeRateArg, &riskFreeRate) != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The risk-free rate must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->kernelParameters.riskFreeRate = riskFreeRate;
	}

	if (dividendYieldArg != NULL)
	{
		double dividendYield;

		if (parseDoubleChecked(dividendYieldArg, &dividendYield) != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The dividend yield must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->kernelParameters.dividendYield = dividendYield;
	}

	if (maturityArg != NULL)
	{
		double maturity;

		if (parseDoubleChecked(maturityArg, &maturity) != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The maturity must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		if (maturity <= 0.0)
		{
			fprintf(stderr, "Error: The maturity must be positive.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->kernelParameters.maturity = maturity;
	}

	if (creditSpreadArg != NULL)
	{
		double creditSpread;

		if (parseDoubleChecked(creditSpreadArg, &creditSpread) != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The credit spread must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		if (creditSpread < 0.0)
		{
			fprintf(stderr, "Error: The credit spread must be non-negative.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->kernelParameters.creditSpread = creditSpread;
	}

	if (ownCreditSpreadArg != NULL)
	{
		double ownCreditSpread;

		if (parseDoubleChecked(ownCreditSpreadArg, &ownCreditSpread) != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The own credit spread must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		if (ownCreditSpread < 0.0)
		{
			fprintf(stderr, "Error: The own credit spread must be non-negative.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->kernelParameters.ownCreditSpread = ownCreditSpread;
	}

	if (strikePriceArg != NULL)
	{
		double strikePrice;

		if (parseDoubleChecked(strikePriceArg, &strikePrice) != kCommonConstantReturnTypeSuccess)
		{
			fprintf(stderr, "Error: The strike price must be a real number.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		if (strikePrice <= 0.0)
		{
			fprintf(stderr, "Error: The strike price must be positive.\n");
			printUsage();

			return kCommonConstantReturnTypeError;
		}

		arguments->kernelParameters.strikePrice = strikePrice;
	}

	return kCommonConstantReturnTypeSuccess;

#endif
}
