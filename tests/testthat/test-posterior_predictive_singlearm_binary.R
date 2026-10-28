
# 1. Validate Beta-Binomial conjugate posterior calculations
test_that("posterior calculations are correct", {

  result <- posterior_predictive_singlearm_binary(
    x_obs = 9,
    n_current = 20,
    N_total = 40,
    S = 100,
    seed = 666
  )

  # Posterior parameters after 9 successes and 11 failures
  expect_equal(result$posterior$alpha, 9.2)
  expect_equal(result$posterior$beta, 11.8)

  # Analytical posterior mean
  expect_equal(result$posterior_mean, 9.2 / 21)

  # Exact posterior probability based on the Beta distribution
  p_exact <- pbeta(
    0.20, 9.2, 11.8,
    lower.tail = FALSE
  )

  expect_equal(result$p_null, p_exact)

  expect_equal(
    result$p_null,
    0.992266339277512,
    tolerance = 1e-10
  )
})


# 2. Validate predictive probabilities against exact Beta-Binomial calculations
test_that("predictive probabilities match exact calculations", {

  S <- 10000

  result <- posterior_predictive_singlearm_binary(
    x_obs = 9,
    n_current = 20,
    N_total = 40,
    S = S,
    seed = 666
  )

  # Current posterior distribution: Beta(9.2, 11.8)
  a <- 9.2
  b <- 11.8
  m <- 20
  y <- 0:m

  # Exact Beta-Binomial probability mass function
  prob_y <- exp(
    lchoose(m, y) +
      lbeta(a + y, b + m - y) -
      lbeta(a, b)
  )

  # Final posterior probabilities for each possible future outcome
  p_final_null <- pbeta(
    0.20, a + y, b + m - y,
    lower.tail = FALSE
  )

  p_final_clin <- pbeta(
    0.40, a + y, b + m - y,
    lower.tail = FALSE
  )

  # Exact posterior predictive probabilities of success
  exact_null <- sum(
    prob_y * (p_final_null >= 0.95)
  )

  exact_clin <- sum(
    prob_y * (p_final_clin >= 0.90)
  )

  # Verify that the predictive probability mass sums to one
  expect_equal(
    sum(prob_y),
    1,
    tolerance = 1e-10
  )

  # Verify the exact predictive probabilities
  expect_equal(
    exact_null,
    0.964021691178353,
    tolerance = 1e-10
  )

  expect_equal(
    exact_clin,
    0.188147586887967,
    tolerance = 1e-10
  )

  # Calculate Monte Carlo standard errors
  se_null <- sqrt(
    exact_null * (1 - exact_null) / S
  )

  se_clin <- sqrt(
    exact_clin * (1 - exact_clin) / S
  )

  # Compare Monte Carlo estimates with exact probabilities
  # using five Monte Carlo standard errors as the tolerance
  expect_lte(
    abs(result$ppos_null - exact_null),
    5 * se_null
  )

  expect_lte(
    abs(result$ppos_clin - exact_clin),
    5 * se_clin
  )
})


# 3. Verify deterministic predictive outcomes at the final analysis
test_that("predictive outcome is deterministic at final analysis", {

  result <- posterior_predictive_singlearm_binary(
    x_obs = 13,
    n_current = 40,
    N_total = 40,
    S = 1,
    seed = 666
  )

  # No future patients remain at the final analysis
  expect_equal(
    result$ppos_null,
    as.numeric(result$p_null >= 0.95)
  )
})


# 4. Verify input validation
test_that("invalid inputs are rejected", {

  # Observed responses cannot exceed the sample size
  expect_error(
    posterior_predictive_singlearm_binary(
      x_obs = 21,
      n_current = 20,
      N_total = 40
    )
  )

  # Observed responses cannot be negative
  expect_error(
    posterior_predictive_singlearm_binary(
      x_obs = -1,
      n_current = 20,
      N_total = 40
    )
  )

  # Beta prior shape parameters must be positive
  expect_error(
    posterior_predictive_singlearm_binary(
      x_obs = 9,
      n_current = 20,
      N_total = 40,
      alpha0 = 0
    )
  )

  # The number of Monte Carlo replications must be positive
  expect_error(
    posterior_predictive_singlearm_binary(
      x_obs = 9,
      n_current = 20,
      N_total = 40,
      S = 0
    )
  )
})
