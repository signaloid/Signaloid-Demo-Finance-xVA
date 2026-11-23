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

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "utilities.h"
#include "kernel.h"

#ifdef NO_OS_AVAILABLE

void
returnZeroNoOS(void);
#endif

int
main(int argc, char *  argv[])
{
	CommandLineArguments    arguments               = (CommandLineArguments) { 0 };
	const char *            applicationDescription  =
		"Computes credit risk valuation adjustments (xVA) using "
		"Geometric Brownian Motion path simulation.";
	const char *    expectedInputHeaders[kXvaInputIndexMax] = {
		"volatility"
	};
	double          outputVariables[kXvaOutputIndexMax] = { 0 };
	const char *    outputVariableNames[kXvaOutputIndexMax] = {
		"CVA",
		"DVA",
		"BCVA",
		"FVA",
		"MVA",
		"KVA",
		"TotalXVA",
		"EuropeanPut"
	};
	const char *    outputVariableDescriptions[kXvaOutputIndexMax] = {
		"Credit Valuation Adjustment",
		"Debit Valuation Adjustment",
		"Bilateral CVA (CVA - DVA)",
		"Funding Valuation Adjustment",
		"Margin Valuation Adjustment",
		"Capital Valuation Adjustment",
		"Total xVA (CVA - DVA + FVA + MVA + KVA)",
		"European Put (Black-Scholes)"
	};
	/*
	 *	Every xVA output is a single converged scalar.
	 */
	kOutputVariableTypeIndex    outputVariableTypes[kXvaOutputIndexMax] = {
		kOutputVariableTypeScalar, kOutputVariableTypeScalar, kOutputVariableTypeScalar, kOutputVariableTypeScalar,
		kOutputVariableTypeScalar, kOutputVariableTypeScalar, kOutputVariableTypeScalar, kOutputVariableTypeScalar,
	};
	double *                    monteCarloOutputSamples = NULL;
	clock_t                     start                   = 0;
	clock_t                     end                     = 0;
	double                      cpuTimeUsedInSeconds    = 0.0;

	/*
	 *	Get command line arguments.
	 */
	if (getCommandLineArguments(argc, argv, &arguments) != kCommonConstantReturnTypeSuccess)
	{
		return EXIT_FAILURE;
	}

	/*
	 *	Read input distributions from CSV if input from file is enabled.
	 */
	if (arguments.common.isInputFromFileEnabled)
	{
		if (readInputDoubleDistributionsFromCSV(
				arguments.common.inputFilePath,
				expectedInputHeaders,
				arguments.kernelParameters.inputVariables,
				kXvaInputIndexMax
		))
		{
			fprintf(stderr, "Error: Could not read from input CSV file \"%s\".\n", arguments.common.inputFilePath);

			return EXIT_FAILURE;
		}
	}

	/*
	 *	Start timing if timing is enabled.
	 */
	if (arguments.common.isTimingEnabled)
	{
		start = clock();
	}

	/*
	 *	Dispatch to the mode-specific kernel.
	 */
	if (arguments.common.isMonteCarloMode)
	{
		calculateOutputMonteCarlo(
			&arguments.kernelParameters,
			arguments.common.outputSelect,
			arguments.common.numberOfMonteCarloIterations,
			outputVariables
		);
	}
	else
	{
		calculateOutputUxHw(
			&arguments.kernelParameters,
			arguments.common.outputSelect,
			outputVariables
		);
	}


	/*
	 *	Stop timing if timing is enabled.
	 */
	if (arguments.common.isTimingEnabled)
	{
		end = clock();
		cpuTimeUsedInSeconds = ((double) (end - start)) / CLOCKS_PER_SEC;
	}


	/*
	 *	`outputSelect` is `kXvaOutputIndexMax` only when computing
	 *	all outputs, which the argument parser rejects in Monte Carlo mode,
	 *	so the types array is never indexed out of range here.
	 */
	bool isSelectedOutputScalar = (arguments.common.outputSelect != kXvaOutputIndexMax) &&
	                              (outputVariableTypes[arguments.common.outputSelect] == kOutputVariableTypeScalar);

	CommonCommandLineArguments printArguments = arguments.common;

	/*
	 *	Persist the selected value to data.out and route every printer
	 *	through the scalar code path. Unlike the demos that carry a
	 *	per-iteration sample array, there is exactly one value to save
	 *	here, so this is gated on the selected output being scalar.
	 */
	if (arguments.common.isMonteCarloMode && isSelectedOutputScalar)
	{
		double selectedOutput = outputVariables[arguments.common.outputSelect];

		saveMonteCarloDoubleDataToDataDotOutFile(
			&selectedOutput,
			(uint64_t) (cpuTimeUsedInSeconds * 1000000),
			1
		);

		printArguments.isMonteCarloMode             = false;
		printArguments.numberOfMonteCarloIterations = 1;
	}

	/*
	 *	Print json outputs if in JSON output mode.
	 */
	if (arguments.common.isOutputJSONMode)
	{
		printJSONFormattedOutput(
			&printArguments,
			monteCarloOutputSamples,
			outputVariables,
			outputVariableNames,
			kXvaOutputIndexMax,
			applicationDescription,
			outputVariableTypes
		);
	}
	/*
	 *	Print human-consumable output if not in JSON output mode.
	 */
	else
	{
		printHumanConsumableOutput(
			&printArguments,
			kXvaOutputIndexMax,
			outputVariables,
			outputVariableNames,
			outputVariableDescriptions,
			monteCarloOutputSamples,
			outputVariableTypes
		);
	}

	/*
	 *	Print timing if timing is enabled.
	 */
	if (arguments.common.isTimingEnabled)
	{
		printf("\nCPU time used: %" SignaloidParticleModifier "lf seconds\n", cpuTimeUsedInSeconds);
	}

	/*
	 *	Save outputs to file if write to file is enabled.
	 */
	if (arguments.common.isWriteToFileEnabled)
	{
		if (writeOutputDoubleDistributionsToCSV(
				arguments.common.outputFilePath,
				outputVariables,
				outputVariableNames,
				kXvaOutputIndexMax
		))
		{
			fprintf(stderr, "Error: Could not write to output CSV file \"%s\".\n", arguments.common.outputFilePath);

			return EXIT_FAILURE;
		}
	}

#ifdef NO_OS_AVAILABLE
	returnZeroNoOS();
#else

	return EXIT_SUCCESS;

#endif
}
