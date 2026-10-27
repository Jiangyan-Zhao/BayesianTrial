# BayesianTrial: Bayesian Decision-Making for Clinical Trials

BayesianTrial currently supports Bayesian decision calculations for
two-stage single-arm trials with binary endpoints:

- `posterior_predictive_singlearm_binary()` calculates posterior summaries
  and predictive probabilities of meeting final success criteria.
- `oc_singlearm_binary()` estimates overall trial success and interim
  stopping probabilities by simulating complete trials.

## Usage guide

The [single-arm phase II trial guide](vignettes/singlearm-binary.Rmd)
explains how to specify binary response data and priors, calculate and
interpret posterior and predictive probabilities, apply posterior stopping
rules, and evaluate type I error, power, and early stopping probabilities.
It includes a hypothetical response-rate example and explains how to adapt
the inputs to another study.

Install the package with its vignette:

```r
install.packages(c("remotes", "knitr", "rmarkdown"))
remotes::install_github(
  "Jiangyan-Zhao/BayesianTrial",
  build_vignettes = TRUE
)
vignette("singlearm-binary", package = "BayesianTrial")
```
