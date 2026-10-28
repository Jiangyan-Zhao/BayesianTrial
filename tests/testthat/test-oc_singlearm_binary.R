
# 1. Validate operating characteristics against exact binomial calculations
test_that("operating characteristics match exact probabilities", {

  # Reference values obtained by exact binomial enumeration
  # Design: Beta(0.2, 0.8) prior, interim n = 20, maximum N = 40
  # Thresholds: interim efficacy = 0.98, futility = 0.50,
  # and final efficacy = 0.95

  reference <- data.frame(
    theta = c(0.20, 0.40),
    success = c(0.0439384767, 0.8569076230),
    futility = c(0.6296482639, 0.0509519532),
    efficacy = c(0.0099817863, 0.4044012747)
  )

  ntrial <- 10000

  for (i in seq_len(nrow(reference))) {

    result <- oc_singlearm_binary(
      theta_true = reference$theta[i],
      n_current = 20,
      N_total = 40,
      alpha0 = 0.2,
      beta0 = 0.8,
      theta0_null = 0.20,
      gamma_early = 0.98,
      gamma_futility = 0.50,
      gamma_final = 0.95,
      ntrial = ntrial,
      seed = 666
    )

    # Calculate Monte Carlo standard errors
    se_success <- sqrt(
      reference$success[i] *
        (1 - reference$success[i]) / ntrial
    )

    se_futility <- sqrt(
      reference$futility[i] *
        (1 - reference$futility[i]) / ntrial
    )

    se_efficacy <- sqrt(
      reference$efficacy[i] *
        (1 - reference$efficacy[i]) / ntrial
    )

    # Compare overall success probability with the exact value
    expect_lte(
      abs(
        result$power_null_post -
          reference$success[i]
      ),
      5 * se_success
    )

    # Compare futility stopping probability with the exact value
    expect_lte(
      abs(
        result$early_stop_rate_futility -
          reference$futility[i]
      ),
      5 * se_futility
    )

    # Compare efficacy stopping probability with the exact value
    expect_lte(
      abs(
        result$early_stop_rate_efficacy -
          reference$efficacy[i]
      ),
      5 * se_efficacy
    )

    # Early efficacy and futility stopping are mutually exclusive
    expect_true(
      result$early_stop_rate_futility +
        result$early_stop_rate_efficacy <= 1
    )

    # Overall success includes all early efficacy successes
    expect_true(
      result$power_null_post >=
        result$early_stop_rate_efficacy
    )
  }
})


# 2. Verify reproducibility with a fixed random seed
test_that("OC simulation is reproducible", {

  result1 <- oc_singlearm_binary(
    theta_true = 0.20,
    n_current = 20,
    N_total = 40,
    ntrial = 1000,
    seed = 123
  )

  result2 <- oc_singlearm_binary(
    theta_true = 0.20,
    n_current = 20,
    N_total = 40,
    ntrial = 1000,
    seed = 123
  )

  # Identical seeds and inputs must produce identical outputs
  expect_identical(result1, result2)
})


# 3. Verify input validation for trial design parameters
test_that("invalid design parameters are rejected", {

  # The true response probability must lie between 0 and 1
  expect_error(
    oc_singlearm_binary(
      theta_true = -0.1,
      n_current = 20,
      N_total = 40
    )
  )

  # The interim sample size cannot exceed the maximum sample size
  expect_error(
    oc_singlearm_binary(
      theta_true = 0.20,
      n_current = 40,
      N_total = 20
    )
  )

  # The number of simulated trials must be positive
  expect_error(
    oc_singlearm_binary(
      theta_true = 0.20,
      n_current = 20,
      N_total = 40,
      ntrial = 0
    )
  )
})


# Verify strict success thresholds and inclusive futility threshold
test_that("interim and final threshold equality is handled correctly", {

  # All response counts are zero when the true response probability is zero
  p_interim <- pbeta(
    0.20, 0.2, 20.8,
    lower.tail = FALSE
  )

  p_final <- pbeta(
    0.20, 0.2, 40.8,
    lower.tail = FALSE
  )

  # Equality with early and final success thresholds is not success
  result <- oc_singlearm_binary(
    theta_true = 0,
    n_current = 20,
    N_total = 40,
    gamma_early = p_interim,
    gamma_futility = 0,
    gamma_final = p_final,
    ntrial = 100,
    seed = 123
  )

  expect_equal(result$early_stop_rate_efficacy, 0)
  expect_equal(result$early_stop_rate_futility, 0)
  expect_equal(result$power_null_post, 0)

  # Equality with the futility threshold causes early futility stopping
  futility_result <- oc_singlearm_binary(
    theta_true = 0,
    n_current = 20,
    N_total = 40,
    gamma_early = 1,
    gamma_futility = p_interim,
    gamma_final = p_final,
    ntrial = 100,
    seed = 123
  )

  expect_equal(futility_result$early_stop_rate_futility, 1)

  # Futility and efficacy thresholds must define disjoint decisions
  expect_error(
    oc_singlearm_binary(
      theta_true = 0.20,
      n_current = 20,
      N_total = 40,
      gamma_early = 0.50,
      gamma_futility = 0.50,
      ntrial = 10
    )
  )
})
