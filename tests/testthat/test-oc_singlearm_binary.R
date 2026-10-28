
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
      gamma_s_1 = 0.98,
      gamma_c = 0.50,
      gamma_f_1 = 0.95,
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
