# BayesianTrial: Bayesian Decision-Making for Clinical Trials

BayesianTrial currently supports Bayesian decision calculations for
two-stage single-arm trials with binary endpoints:

- `posterior_predictive_singlearm_binary()` calculates posterior summaries
  and predictive probabilities of meeting final success criteria.
- `oc_singlearm_binary()` estimates overall trial success and interim
  stopping probabilities by simulating complete trials.

## Reproduce the manuscript example

The [paper reproduction vignette](vignettes/paper-reproduction.Rmd)
reproduces the interim analysis and Table 1 in Zhao and Mu,
*Introduction to Bayesian Decision Method and Its Application in Medical
Research* (manuscript). It uses the hypothetical ORR example with 9 responses
among 20 patients, a maximum of 40 patients, and a Beta(0.2, 0.8) prior.

Install the optional `rpact` package before building the vignette to also
run the manuscript's frequentist group sequential comparison:

```r
install.packages(c("remotes", "knitr", "rmarkdown", "rpact"))
remotes::install_github(
  "Jiangyan-Zhao/BayesianTrial",
  build_vignettes = TRUE
)
vignette("paper-reproduction", package = "BayesianTrial")
```

The Bayesian calculations use 10,000 Monte Carlo replications and seed 666,
matching the reference script. The vignette reports the settings, decision
rules, generated results, and R session information.

