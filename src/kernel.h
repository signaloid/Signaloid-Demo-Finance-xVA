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

#pragma once

#include "utilities.h"
#include "xva-config.h"

/**
 *	@brief	Monte Carlo calculation kernel. Simulates
 *		`arguments->common.numberOfMonteCarloIterations` independent
 *		Brownian paths, forms the exposure profile as their
 *		cross-sectional average, and computes the selected xVA metric(s)
 *		from it using scalar sample statistics.
 *
 *	@param	arguments	: Pointer to command line arguments struct.
 *	@param	inputVariables	: The input variables.
 *	@param	outputVariables	: The output variables.
 */
void
calculateOutputMonteCarlo(
	CommandLineArguments *  arguments,
	double *                inputVariables,
	double *                outputVariables);

/**
 *	@brief	UxHw calculation kernel. Simulates a single distributional
 *		Brownian path, forms the exposure profile by integrating that
 *		distribution out with UxHw distributional arithmetic, and
 *		computes the selected xVA metric(s) from it.
 *
 *	@param	arguments	: Pointer to command line arguments struct.
 *	@param	inputVariables	: The input variables.
 *	@param	outputVariables	: The output variables.
 */
void
calculateOutputUxHw(
	CommandLineArguments *  arguments,
	double *                inputVariables,
	double *                outputVariables);
