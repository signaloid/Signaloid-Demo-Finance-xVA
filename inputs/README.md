# Input Files

## input.csv
Example input for the `-i` option. The demo has a single uncertain input, the
annualised volatility, so the CSV has one column headed `volatility`. The
header must match `expectedInputHeaders` in `src/main.c`. Each row is one
sample, and the samples become a distribution via
`UxHwDoubleDistFromSamples()`.

The shipped file holds 10,000 samples drawn from a normal distribution with
mean 0.20 and standard deviation 0.02.

`-i` is accepted only in UxHw single-shot mode (no `-M`). Monte Carlo mode
rejects it. On the Signaloid Cloud Compute Engine the volatility is then a
distribution and every output is collapsed to its mean. Natively, the
compatibility layer picks one sample at random.

```bash
./demo-native-mc -S 7 -i inputs/input.csv
```
