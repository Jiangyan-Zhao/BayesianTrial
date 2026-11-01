
<!-- README.md is generated from README.Rmd. Please edit README.Rmd only. -->

# BayesianTrial

[![R-CMD-check](https://github.com/Jiangyan-Zhao/BayesianTrial/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jiangyan-Zhao/BayesianTrial/actions/workflows/R-CMD-check.yaml)

**Bayesian decision-making for clinical trials in R.**

BayesianTrial provides posterior and predictive probability calculations
and simulation-based operating characteristics for single-arm phase II trials
with binary efficacy endpoints. The current implementation supports two-stage
designs with posterior-probability-based stopping rules.

## Installation

Install the package from GitHub using [pak](https://pak.r-lib.org/):

```r
install.packages("pak")
pak::pak("Jiangyan-Zhao/BayesianTrial")
```

## Current functionality

- `posterior_predictive_singlearm_binary()` calculates posterior summaries
  and predictive probabilities of meeting final posterior success criteria.
- `oc_singlearm_binary()` estimates overall trial success, early efficacy
  stopping, and early futility stopping through Monte Carlo simulation.

## Quick start

Consider a hypothetical trial with a historical response rate of 20%,
an interim analysis after 20 patients, and a maximum sample size of 40.
Nine patients respond at interim. We use the default Beta(0.2, 0.8) prior.

### Posterior and predictive probabilities

```r
library(BayesianTrial)

prob <- posterior_predictive_singlearm_binary(
  x_obs = 9, n_current = 20, N_total = 40,
  type = "both"
)

prob$p_null
prob$ppos_null
```

`p_null` is the current posterior probability that the response rate
exceeds 0.20. `ppos_null` is the predictive probability of meeting the
final success criterion if the planned sample size is completed.
The default `type = "posterior"` returns posterior summaries; use
`type = "predictive"` for predictive probabilities or `type = "both"`
for both.

### Operating characteristics

```r
oc <- oc_singlearm_binary(
  theta_true = 0.20, n_current = 20, N_total = 40
)

oc
```

The overall success probability (`power_null_post`) estimates the type I
error rate when `theta_true = 0.20` and statistical power when
`theta_true = 0.40`. Early stopping probabilities are also returned.

For more details, see the
[single-arm binary endpoint vignette](vignettes/singlearm-binary.Rmd).

## Planned extensions

Future versions may extend to other trial designs and endpoint types.

## License

GPL (\>= 3).
