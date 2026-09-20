# Source code

## main.c
Entry point. Parses the command line, reads the volatility from the input CSV
when `-i` is given, dispatches to the Monte Carlo or UxHw kernel on
`isMonteCarloMode`, and prints or saves the selected outputs.

## kernel.c/h
The two calculation kernels, `calculateOutputMonteCarlo()` and
`calculateOutputUxHw()`. Both build the same model data and compute the same
xVA metrics from an exposure profile. They differ only in how that profile is
formed.

## xva-monte-carlo.c/h and xva-uxhw.c/h
The mode-specific halves of the engine, covering GBM path simulation and the
exposure profile. The Monte Carlo half averages a scalar path ensemble. The
UxHw half integrates one distributional path in a single pass.

## xva-config.h
Central header for the xVA framework. It holds the `kXva` constants, the `Xva`
data structures and the prototypes of the `xva` functions in the modules below.

## xva-gbm.c, xva-exposure.c, xva-utils.c
Shared helpers. Allocation of paths and exposure profiles, the forward
mark-to-market, discounting, survival probabilities, Black-Scholes, quantiles
and other numerical routines.

## xva-cva.c, xva-dva.c, xva-fva.c, xva-mva.c, xva-kva.c
The mode-agnostic valuation adjustments, each computed from an exposure
profile.

## utilities.c/h
Demo-specific command-line arguments (`-f`, `-s`, `-r`, `-d`, `-m`, `-c`,
`-e`, `-k`) and their defaults, layered on the common arguments handled by
`common.c`.

## common.c/h
Utility methods for parsing, setting, and reporting the usage of command-line
arguments common to all of our C/C++ demo applications, as well as other
methods commonly used across them, such as CSV and JSON I/O. These source
files are symlinks to the original files contained in the repository
[Signaloid-Demo-CommonUtilityRoutines](https://github.com/signaloid/Signaloid-Demo-CommonUtilityRoutines),
which is included as a submodule in `submodules/common`.

## uxhw.c/h
Methods that implement the probabilistic versions of the methods in the UxHw
API (e.g., `UxHwDoubleGaussDist`) using the GNU Scientific Library (GSL)
random number generators. This allows building our C/C++ demo applications
natively (i.e., on conventional architectures) and running native Monte Carlo
evaluations without modifying the source code. These source files are symlinks
to the original files contained in the repository
[Signaloid-Demo-UxHwCompatibilityForNativeExecution](https://github.com/signaloid/Signaloid-Demo-UxHwCompatibilityForNativeExecution),
which is included as a submodule in `submodules/compat`.

## config.mk
The source list shared by the Signaloid cloud build and the native build. New
source files must be added here.

# Building natively
Use the top-level `Makefile`, which reads `config.mk`, adds `uxhw.c` and links
GSL:

```bash
make local-build
```
