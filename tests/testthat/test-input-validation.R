# Verify informative and consistent errors for invalid numeric inputs

test_that("OC validates scalar probabilities", {
  base <- list(theta_true = 0.20, n_current = 20, N_total = 40, ntrial = 1)

  invalid_values <- list(NA_real_, NaN, Inf, -Inf, numeric(0),
                         c(0.20, 0.40), "0.20", TRUE, matrix(0.20))

  for (value in invalid_values) {
    args <- base
    args$theta_true <- value
    expect_error(
      do.call(oc_singlearm_binary, args),
      "`theta_true` must be a single finite number"
    )
  }

  for (name in c("theta0_null", "gamma_early",
                 "gamma_futility", "gamma_final")) {
    args <- base
    args[[name]] <- NaN
    expect_error(do.call(oc_singlearm_binary, args), name)
  }

  for (value in c(-0.01, 1.01)) {
    args <- base
    args$theta_true <- value
    expect_error(do.call(oc_singlearm_binary, args), "`theta_true`")
  }

  base$gamma_futility <- 0.98
  expect_error(
    do.call(oc_singlearm_binary, base),
    "`gamma_futility` must be strictly less than `gamma_early`"
  )
})

test_that("OC validates integer counts, positive priors, and seeds", {
  base <- list(theta_true = 0.20, n_current = 20, N_total = 40, ntrial = 1)

  cases <- list(
    list(name = "n_current", value = 10.5),
    list(name = "n_current", value = NA_real_),
    list(name = "N_total", value = Inf),
    list(name = "N_total", value = c(40, 50)),
    list(name = "ntrial", value = 0),
    list(name = "ntrial", value = 1.5),
    list(name = "ntrial", value = .Machine$integer.max + 1),
    list(name = "alpha0", value = 0),
    list(name = "alpha0", value = NaN),
    list(name = "beta0", value = -1),
    list(name = "beta0", value = c(0.8, 1)),
    list(name = "seed", value = NA_real_),
    list(name = "seed", value = 1.5),
    list(name = "seed", value = c(1, 2)),
    list(name = "seed", value = "123")
  )

  for (case in cases) {
    args <- base
    args[[case$name]] <- case$value
    expect_error(
      do.call(oc_singlearm_binary, args),
      paste0("`", case$name, "`")
    )
  }

  args <- base
  args$N_total <- 10
  expect_error(
    do.call(oc_singlearm_binary, args),
    "`N_total` must be greater than or equal to `n_current`"
  )

  # A NULL seed is a supported request to continue the current RNG state
  args <- base
  args$seed <- NULL
  expect_type(do.call(oc_singlearm_binary, args), "list")
})

test_that("posterior-predictive function validates its inputs", {
  base <- list(x_obs = 9, n_current = 20, N_total = 40, S = 1)

  cases <- list(
    list(name = "x_obs", value = NA_real_),
    list(name = "x_obs", value = NaN),
    list(name = "x_obs", value = 2.5),
    list(name = "x_obs", value = c(9, 10)),
    list(name = "x_obs", value = -1),
    list(name = "n_current", value = Inf),
    list(name = "n_current", value = -1),
    list(name = "N_total", value = NA_real_),
    list(name = "alpha0", value = -1),
    list(name = "beta0", value = Inf),
    list(name = "theta0_null", value = NA_real_),
    list(name = "theta0_clin", value = c(0.40, 0.45)),
    list(name = "gamma_final", value = 1.5),
    list(name = "gamma_clin_final", value = NaN),
    list(name = "S", value = 0),
    list(name = "S", value = 1.1),
    list(name = "seed", value = NaN),
    list(name = "seed", value = 1.1),
    list(name = "seed", value = c(1, 2))
  )

  for (case in cases) {
    args <- base
    args[[case$name]] <- case$value
    expect_error(
      do.call(posterior_predictive_singlearm_binary, args),
      paste0("`", case$name, "`")
    )
  }

  args <- base
  args$x_obs <- 21
  expect_error(
    do.call(posterior_predictive_singlearm_binary, args),
    "`x_obs` must not exceed `n_current`"
  )

  args <- base
  args$N_total <- 19
  expect_error(
    do.call(posterior_predictive_singlearm_binary, args),
    "`N_total` must be greater than or equal to `n_current`"
  )

  args <- base
  args$seed <- NULL
  expect_type(do.call(posterior_predictive_singlearm_binary, args), "list")
})

test_that("input errors do not change the random number generator state", {
  set.seed(42)
  initial_rng <- .Random.seed

  expect_error(
    oc_singlearm_binary(
      theta_true = 0.20,
      n_current = 20,
      N_total = 40,
      ntrial = 1,
      seed = NA_real_
    ),
    "`seed`"
  )
  expect_identical(.Random.seed, initial_rng)

  expect_error(
    posterior_predictive_singlearm_binary(
      x_obs = 9,
      n_current = 20,
      N_total = 40,
      S = 1,
      seed = Inf
    ),
    "`seed`"
  )
  expect_identical(.Random.seed, initial_rng)
})
