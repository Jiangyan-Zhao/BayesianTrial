# Internal input validation helpers for single-arm binary trials

.is_finite_numeric_scalar <- function(x) {
  is.numeric(x) &&
    !is.complex(x) &&
    length(x) == 1L &&
    is.null(dim(x)) &&
    is.finite(x)
}

.validate_probability <- function(x, name) {
  if (!.is_finite_numeric_scalar(x) || x < 0 || x > 1) {
    stop(
      sprintf("`%s` must be a single finite number between 0 and 1.", name),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.validate_positive_number <- function(x, name) {
  if (!.is_finite_numeric_scalar(x) || x <= 0) {
    stop(
      sprintf("`%s` must be a single finite positive number.", name),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.validate_count <- function(x, name, min_value = 0L) {
  if (!.is_finite_numeric_scalar(x) ||
      x < min_value ||
      x > .Machine$integer.max ||
      x != floor(x)) {
    requirement <- if (min_value == 0L) "non-negative" else "positive"
    stop(
      sprintf(
        "`%s` must be a single finite %s integer no greater than %d.",
        name, requirement, .Machine$integer.max
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(NULL))
  }
  if (!.is_finite_numeric_scalar(seed) ||
      seed != floor(seed) ||
      abs(seed) > .Machine$integer.max) {
    stop(
      "`seed` must be NULL or a single finite integer within the supported range.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.validate_sample_sizes <- function(n_current, N_total) {
  .validate_count(n_current, "n_current")
  .validate_count(N_total, "N_total")
  if (N_total < n_current) {
    stop(
      "`N_total` must be greater than or equal to `n_current`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.validate_prior <- function(alpha0, beta0) {
  .validate_positive_number(alpha0, "alpha0")
  .validate_positive_number(beta0, "beta0")
  if (!is.finite(alpha0 + beta0)) {
    stop(
      "The sum of `alpha0` and `beta0` must be finite.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
