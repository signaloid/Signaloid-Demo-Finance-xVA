[<img src="https://assets.signaloid.io/add-to-signaloid-cloud-logo-dark-latest.png#gh-dark-mode-only" alt="[Add to signaloid.io]" height="30">](https://signaloid.io/repositories?connect=https://github.com/signaloid/Signaloid-Demo-Finance-xVA#gh-dark-mode-only)
[<img src="https://assets.signaloid.io/add-to-signaloid-cloud-logo-light-latest.png#gh-light-mode-only" alt="[Add to signaloid.io]" height="30">](https://signaloid.io/repositories?connect=https://github.com/signaloid/Signaloid-Demo-Finance-xVA#gh-light-mode-only)


# xVA Monte Carlo Simulation with Signaloid UxHw® Support

This application computes counterparty credit risk valuation adjustments (xVA)
using Monte Carlo simulation, with support for Signaloid's distributional
arithmetic (UxHw) as an alternative to traditional Monte Carlo sampling. The
framework computes CVA, DVA, FVA, MVA, and KVA for a forward contract on an
equity underlying, using Geometric Brownian Motion path simulation.  The
repository has been structured to fit within the
[Signaloid UxHw benchmarking template](https://github.com/signaloid),
allowing direct comparison between traditional Monte Carlo sampling and
Signaloid's single-pass uncertainty tracking on the Signaloid Cloud
Compute Engine (SCCE).


## Cloning the repository

```bash
git clone --recursive <repository-url>
```

To update all submodules:
```bash
git pull --recurse-submodules
git submodule update --remote --recursive
```

If you forgot to clone with `--recursive`:
```bash
git submodule update --init --recursive
```


## Quickstart

### Native build and run
```bash
make local-build              # Build the executable
make local-run                # Build and run all 8 outputs with timing

# Run individual outputs:
./demo-native-mc -M 10000 -S 0 -T    # CVA with 10000 MC iterations and timing
./demo-native-mc -M 10000 -S 7 -T    # European Put
./demo-native-mc -S 0                 # UxHw kernel (natively one random path)
./demo-native-mc -M 10000 -S 0 -j    # JSON output
```

## Overview

In traditional Monte Carlo mode (`-M N`), the application simulates N GBM asset
paths, averages them into one expected exposure profile, and computes the
selected xVA from that profile. `-M 1` is a one-path Monte Carlo run, not the
UxHw mode.

In UxHw single-shot mode (no `-M`, on the Signaloid Cloud Compute Engine), the
application simulates one path whose values carry full distributional
information. Signaloid's uncertainty-tracking hardware propagates that
distribution through the exposure computation in a single pass, and the
expectation is then extracted so that each output is a single scalar. Natively
the same mode runs against the compatibility layer and yields one random path.

## Modeling framework

The demo values a forward contract on a dividend-paying equity underlying
against a single counterparty, then charges the five valuation adjustments
against that position. Every formula below is the one the code evaluates, on a
uniform time grid $0 = t_0 < t_1 < \dots < t_N = T$ with step $\Delta t = T/N$.

### Underlying and exposure

The underlying follows Geometric Brownian Motion under the risk-neutral
measure, simulated with the exact log-Euler scheme in
[`xva-gbm.c`](src/xva-gbm.c):

$$
dS(t) = (r - q)\,S(t)\,dt + \sigma\,S(t)\,dW(t)
$$

$$
S(t_{i+1}) = S(t_i)\exp\left[\left(r - q - \tfrac{1}{2}\sigma^{2}\right)\Delta t + \sigma\sqrt{\Delta t}\,Z_i\right],
\qquad Z_i \sim N(0, 1)
$$

The mark-to-market value of the forward at time $t$, in
[`xva-exposure.c`](src/xva-exposure.c), is

$$
V(t) = S(t)e^{-q(T-t)} - Ke^{-r(T-t)}
$$

From $V$ the code builds the four exposure measures that every adjustment is
derived from:

$$
\mathrm{EE}(t) = \mathbb{E}[V(t)], \qquad
\mathrm{EPE}(t) = \mathbb{E}[\max(V(t), 0)]
$$

$$
\mathrm{ENE}(t) = \mathbb{E}[\min(V(t), 0)], \qquad
\mathrm{PFE}(t) = Q_{0.975}\!\left(\max(V(t), 0)\right)
$$

This is the single point where the two execution modes differ. Monte Carlo mode
forms each expectation as a cross-sectional average over the path ensemble,
while UxHw mode integrates one distributional path out in a single pass.

### Credit and discounting

Hazard rates are bootstrapped from credit spreads under the standard credit
triangle, and survival is exponential:

$$
\lambda = \frac{s}{1 - R}, \qquad
Q(t) = e^{-\lambda t}, \qquad
\mathrm{DF}(t) = e^{-rt}
$$

The default probability over a grid interval is $Q(t_{i-1}) - Q(t_i)$.

### Valuation adjustments

With loss given default $\mathrm{LGD} = 1 - R$, the counterparty carrying
subscript $C$ and the bank subscript $B$, the five adjustments are

$$
\mathrm{CVA} = \mathrm{LGD}_C \sum_{i=1}^{N} \mathrm{EPE}(t_i)\,\mathrm{DF}(t_i)\left[Q_C(t_{i-1}) - Q_C(t_i)\right]
$$

$$
\mathrm{DVA} = \mathrm{LGD}_B \sum_{i=1}^{N} \left|\mathrm{ENE}(t_i)\right|\mathrm{DF}(t_i)\left[Q_B(t_{i-1}) - Q_B(t_i)\right]
$$

$$
\mathrm{BCVA} = \mathrm{CVA} - \mathrm{DVA}
$$

FVA is the funding cost net of the funding benefit, each weighted by the
survival of the party that bears it:

$$
\mathrm{FVA} = \underbrace{\sum_{i} \mathrm{EPE}(t_i)\,Q_C(t_i)\,s_f\,\mathrm{DF}(t_i)\,\Delta t}_{\text{FCA}}
- \underbrace{\sum_{i} \left|\mathrm{ENE}(t_i)\right|Q_B(t_i)\,s_l\,\mathrm{DF}(t_i)\,\Delta t}_{\text{FBA}}
$$

The demo's funding data lends at the risk-free rate, so $s_l = 0$ and only the
FCA term contributes.

MVA charges the funding spread against posted initial margin. The demo
approximates the initial margin from the PFE profile scaled to the ten-day
margin period of risk, and weights it by counterparty survival:

$$
\mathrm{MVA} = \sum_{i} \mathrm{IM}(t_i)\,s_f\,Q_C(t_i)\,\mathrm{DF}(t_i)\,\Delta t,
\qquad
\mathrm{IM}(t_i) = \mathrm{PFE}(t_i)\sqrt{10/252}
$$

KVA charges the hurdle rate against regulatory capital. Capital is driven by
the Basel effective expected positive exposure, which is a *running maximum*
over the expected profile rather than a per-path quantity:

$$
\mathrm{EEPE}(t_i) = \max_{j \le i} \mathrm{EPE}(t_j), \qquad
\mathrm{EAD}(t_i) = \alpha \cdot \mathrm{EEPE}(t_i)
$$

$$
\mathrm{KVA} = \sum_{i} \underbrace{w \cdot 12.5 \cdot \mathrm{EAD}(t_i) \cdot c}_{\text{regulatory capital}} \cdot h \cdot Q(t_i)\,\mathrm{DF}(t_i)\,\Delta t
$$

with Basel multiplier $\alpha = 1.4$, risk weight $w = 100\%$, capital ratio
$c = 8\%$ and hurdle rate $h = 10\%$. The survival weight $Q(t)$ in KVA uses a
fixed hazard rate of $2\%$ (`kXvaUnilateralHazardRate`) rather than the
counterparty's calibrated rate. The running maximum is why the expectation has
to be formed before the adjustment, as the [Dual execution
modes](#dual-execution-modes) section explains.

The total is the sum of the five, and the eighth output is the Black-Scholes
European put written on the same underlying, which bypasses the exposure
profile entirely:

$$
\mathrm{xVA} = \mathrm{CVA} - \mathrm{DVA} + \mathrm{FVA} + \mathrm{MVA} + \mathrm{KVA}
$$

$$
P = Ke^{-rT}\Phi(-d_2) - S_0e^{-qT}\Phi(-d_1), \qquad
d_{1,2} = \frac{\ln(S_0/K) + (r - q \pm \tfrac{1}{2}\sigma^{2})T}{\sigma\sqrt{T}}
$$

## Output indices

| `-S` | Name | Description |
|------|------|-------------|
| 0 | CVA | Credit Valuation Adjustment |
| 1 | DVA | Debit Valuation Adjustment |
| 2 | BCVA | Bilateral CVA (CVA - DVA) |
| 3 | FVA | Funding Valuation Adjustment |
| 4 | MVA | Margin Valuation Adjustment |
| 5 | KVA | Capital Valuation Adjustment |
| 6 | Total xVA | Total xVA (CVA - DVA + FVA + MVA + KVA) |
| 7 | European Put | European Put (Black-Scholes) |

## Building and running

### Prerequisites (native build)
- GCC with C11 support
- GSL (GNU Scientific Library): `sudo apt install libgsl-dev`

### Native build and run
```bash
make local-build              # Build the executable
make local-run                # Build and run all 8 outputs with timing

# Run individual outputs:
./demo-native-mc -M 10000 -S 0 -T    # CVA with 10000 MC iterations and timing
./demo-native-mc -M 10000 -S 7 -T    # European Put
./demo-native-mc -S 0                 # UxHw kernel (natively one random path)
./demo-native-mc -M 10000 -S 0 -j    # JSON output
```

### Command-line options
```
Common options:
  -S <index>       Select output (0-7, see table above). 8, or no -S,
                   computes all outputs and is allowed in UxHw mode only
  -M <iterations>  Monte Carlo mode with N paths. N is required and has no
                   default. Without -M the UxHw kernel runs
  -i <path>        Read the volatility distribution from a CSV file
                   (UxHw mode only, see inputs/README.md)
  -o <path>        Write the outputs to a CSV file
  -T               Enable timing
  -b               Benchmarking output format
  -j               JSON output mode
  -v               Verbose (accepted but not supported by this demo)
  -h               Help

xVA-specific options:
  -f <volatility>        Annualized volatility (default: 0.20)
  -s <spot-price>        Spot price (default: 100.0)
  -r <risk-free-rate>    Risk-free rate (default: 0.05)
  -d <dividend-yield>    Dividend yield (default: 0.02)
  -m <maturity>          Time to maturity in years (default: 5.0)
  -c <credit-spread>     Counterparty credit spread (default: 0.015)
  -e <own-credit-spread> Own credit spread (default: 0.01)
  -k <strike-price>      Strike price (default: 100.0)
```

## Architecture

### Dual execution modes
The two modes are implemented as separate kernels, dispatched from `main.c` on
`isMonteCarloMode`. Both build the same model data and compute the same xVA
metrics from an exposure profile. They differ only in how that profile is
formed.

- **Monte Carlo mode** (`-M N`): `calculateOutputMonteCarlo()` in `kernel.c`, backed by `xva-monte-carlo.c`. Simulates N independent scalar GBM paths and forms the exposure profile as their cross-sectional average, with PFE as the cross-sectional sample quantile.
- **UxHw single-shot mode** (no `-M`): `calculateOutputUxHw()` in `kernel.c`, backed by `xva-uxhw.c`. Simulates one distributional GBM path on SCCE, where `UxHwDoubleGaussDist()` yields a full distribution and `UxHwDoubleGetIndependentCopy()` breaks inter-step correlations. The exposure profile integrates that distribution out in a single pass: EE, EPE and ENE are first moments and PFE is the high quantile of the positive-exposure distribution.

Forming the expectation before the xVA metrics matters for the nonlinear ones.
KVA's effective EPE is a running maximum over the *expected* exposure profile,
so averaging per-path xVAs would instead compute a different, biased
$\mathbb{E}[\max(\cdot, 0)]$.

### Module flow
```
Input (volatility)
  -> GBM Path Simulation      (xva-monte-carlo.c | xva-uxhw.c, allocation in xva-gbm.c)
    -> Exposure Profile       (xva-monte-carlo.c | xva-uxhw.c, allocation in xva-exposure.c)
      -> xVA Computations     (xva-cva.c, xva-dva.c, xva-fva.c, xva-mva.c, xva-kva.c)
        -> Output
```

### Uncertainty model
Volatility is a calibrated market input, held at a fixed deterministic value
so that the Monte Carlo and UxHw modes compute the same production xVA. The
stochasticity that drives exposure is the Brownian path ensemble, meaning the
per-step `UxHwDoubleGaussDist()` draws in the GBM simulation, rather than the
volatility itself.

A distributional volatility can still be supplied with `-i` in UxHw mode,
since the CSV samples become a distribution via
`UxHwDoubleDistFromSamples()`. The CSV needs a single column headed
`volatility` (see [`inputs/README.md`](inputs/README.md)) and Monte Carlo mode
rejects `-i`. Every output stays scalar in that case. The seven
profile-derived metrics are collapsed by the first-moment extraction in
`xva-uxhw.c`, and the European put, which never passes through the exposure
profile, is collapsed by `xvaEuropeanPutUxHw()`.


## References

The xVA framework implemented in this repository draws on the following
academic papers, regulatory documents, and textbooks. Entries link to a
canonical source where one is freely resolvable, and carry full bibliographic
detail otherwise. Nothing is redistributed here.

### Monte Carlo and numerical methods
Gaussian variates come from `UxHwDoubleGaussDist()`. On SCCE this is a
distributional primitive, and natively the compatibility submodule draws them
with GSL, so the repository carries no random number generator of its own.
- Glasserman, P. (2004). *Monte Carlo Methods in Financial Engineering.* Springer.
  [Source](https://doi.org/10.1007/978-0-387-21617-1)
  Used in: `xva-gbm.c` (GBM path simulation).

### Options pricing and credit risk
- Black, F. & Scholes, M. (1973). "The Pricing of Options and Corporate Liabilities." *Journal of Political Economy*, 81(3), 637-654.
  [Source](https://doi.org/10.1086/260062)
  Used in: `xva-utils.c` (Black-Scholes formula for European options).
- Merton, R.C. (1974). "On the Pricing of Corporate Debt: The Risk Structure of Interest Rates." *Journal of Finance*, 29(2), 449-470.
  [Source](https://doi.org/10.1111/j.1540-6261.1974.tb03058.x)
  Background only. The code uses the reduced-form hazard-rate model, not Merton's structural model.
- Hull, J.C. (2017). *Options, Futures, and Other Derivatives.* 10th Ed., Pearson.
  Used in: `xva-exposure.c` (forward contract valuation, Chapter 5).

### xVA and valuation adjustments
- Gregory, J. (2020). *The xVA Challenge: Counterparty Risk, Funding, Collateral, Capital and Initial Margin.* 4th Ed., Wiley.
  Primary reference for: `xva-cva.c`, `xva-dva.c`, `xva-fva.c`, `xva-mva.c`, `xva-kva.c`, `xva-exposure.c`.
- Albanese, C., Caenazzo, S. & Crepey, S. (2016). "Capital Valuation Adjustment and Funding Valuation Adjustment." arXiv:1603.03012.
  [Source](https://arxiv.org/abs/1603.03012)
  Used in: `xva-kva.c` (KVA and FVA framework).
- Green, A., Kenyon, C. & Dennis, C. (2014). "KVA: Capital Valuation Adjustment by Replication." *Risk Magazine*, December 2014.
  Used in: `xva-kva.c` (KVA computation).
- Burgard, C. & Kjaer, M. (2011). "Partial Differential Equation Representations of Derivatives with Bilateral Counterparty Risk and Funding Costs." *Journal of Credit Risk*, 7(3), 1-19.
  Used in: `xva-dva.c` (bilateral counterparty risk).
- Pykhtin, M. & Zhu, S. (2007). "A Guide to Modeling Counterparty Credit Risk." *GARP Risk Review*, July/August 2007.
  Used in: `xva-exposure.c`, `xva-cva.c` (exposure and CVA methodology).

### Regulatory documents
- BCBS (2011). "Basel III: A Global Regulatory Framework for More Resilient Banks and Banking Systems." Bank for International Settlements (BCBS 189).
  [Source](https://www.bis.org/publ/bcbs189.htm)
  Used in: `xva-kva.c` (regulatory capital requirements).
- BCBS-IOSCO (2015). "Margin Requirements for Non-Centrally Cleared Derivatives." Bank for International Settlements.
  [Source](https://www.bis.org/bcbs/publ/d317.htm)
  Used in: `xva-mva.c` (initial margin requirements).
- ISDA (2016). "ISDA SIMM: Standard Initial Margin Model." International Swaps and Derivatives Association.
  Used in: `xva-mva.c` (industry-standard initial margin, cited as background for the PFE-based approximation).

### Open source implementations
- Acadia Inc. (2024). "Open Source Risk Engine Design." [github.com/OpenSourceRisk/Engine](https://github.com/OpenSourceRisk/Engine).
  Reference implementation for xVA computation pipelines.

### Signaloid documentation
- [Signaloid Developer Documentation](https://docs.signaloid.io)
- [UxHw API Reference](https://docs.signaloid.io/docs/uxhw-api/)

