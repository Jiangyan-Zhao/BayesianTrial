#' Operating Characteristics of a Two-Stage Single-Arm Phase II Trial with a Binary Endpoint
#'
#' Evaluate a prespecified two-stage Bayesian decision rule for a single-arm
#' phase II trial with a binary efficacy endpoint. Monte Carlo simulation
#' estimates the probability of a positive trial conclusion and the
#' probabilities of early termination for efficacy or futility at an assumed
#' true response rate.
#'
#' @param theta_true Assumed true treatment response rate for generating trial outcomes; a
#'   probability between 0 and 1 (inclusive). Specify the null-scenario rate
#'   to assess the type I error rate or a clinically promising alternative
#'   rate to assess power.
#' @param n_current Planned number of patients with evaluable outcomes at the single interim
#'   analysis; a non-negative integer not exceeding \code{N_total}. A
#'   conventional two-stage design requires \code{0 < n_current < N_total}.
#' @param N_total Planned maximum number of evaluable patients; a non-negative integer at
#'   least as large as \code{n_current}.
#' @param alpha0,beta0 Positive, finite shape parameters of the Beta prior for the response
#'   rate. Their sum must also be finite.
#' @param theta0_null Prespecified unacceptable or historical reference response rate used as
#'   the primary efficacy benchmark; a probability between 0 and 1
#'   (inclusive).
#' @param gamma_early Interim posterior probability cutoff for early declaration of efficacy; a
#'   probability between 0 and 1 (inclusive). An early positive result is
#'   declared if the posterior probability of exceeding \code{theta0_null} is
#'   strictly greater than this cutoff.
#' @param gamma_futility Interim posterior probability cutoff for stopping due to futility; a
#'   probability between 0 and 1 (inclusive) and strictly less than
#'   \code{gamma_early}. Futility stopping occurs if the posterior probability
#'   of exceeding \code{theta0_null} is less than or equal to this cutoff.
#'   This is not a predictive probability cutoff.
#' @param gamma_final Posterior probability cutoff for declaring a positive result at the final
#'   analysis if the trial continues; a probability between 0 and 1
#'   (inclusive).
#' @param ntrial Positive integer number of independently simulated trials.
#' @param seed Integer random seed, or \code{NULL} to use the current random number
#'   generator state. A non-\code{NULL} value calls \code{set.seed()} and
#'   changes the session's random number generator state.
#'
#' @details
#' For each simulated trial, binary efficacy outcomes are generated at the
#' assumed true response rate \code{theta_true}. After \code{n_current}
#' patients have evaluable outcomes, the Beta posterior is updated and the
#' posterior probability that the unknown response rate exceeds
#' \code{theta0_null} is calculated.
#'
#' At the interim analysis, a positive result is declared and the trial stops
#' early for efficacy if this posterior probability is strictly greater than
#' \code{gamma_early}. The trial stops early for futility if the probability
#' is less than or equal to \code{gamma_futility}. Otherwise, the trial
#' continues to \code{N_total} evaluable patients. At the final analysis, a
#' positive result is declared only if the updated posterior probability is
#' strictly greater than \code{gamma_final}. Early efficacy declarations are
#' included in the overall probability of a positive trial result.
#'
#' The overall positive-result probability is estimated as
#' \code{power_null_post}. When \code{theta_true} equals the prespecified
#' null-scenario response rate, this quantity estimates the type I error
#' probability at that design point; when \code{theta_true} is a specified
#' promising alternative rate, it estimates statistical power. The other
#' outputs quantify early termination for futility and efficacy. The function
#' estimates operating characteristics for the supplied decision rule; it
#' neither calibrates the thresholds nor guarantees type I error control or
#' target power.
#'
#' This function implements one interim analysis with posterior-probability
#' stopping rules, not interim monitoring based on predictive probability or
#' the optimized stopping boundaries of the BOP2 design. It accepts boundary
#' sample sizes for computational purposes, but a genuine two-stage design has
#' \code{0 < n_current < N_total}. For reproducibility, responses for the
#' remaining patients are also simulated after an interim stop; these outcomes
#' are ignored in the trial decision.
#'
#' @return A named list containing:
#' \describe{
#'   \item{power_null_post}{Estimated overall probability of declaring a positive trial result,
#'     including early efficacy stopping. Interpreted as type I error
#'     probability under a null scenario or as power under a specified
#'     alternative.}
#'   \item{early_stop_rate_futility}{Estimated probability of stopping the trial at the interim analysis
#'     for futility.}
#'   \item{early_stop_rate_efficacy}{Estimated probability of stopping the trial at the interim analysis
#'     for efficacy.}
#' }
#'
#' @references
#' Thall, P. F., and Simon, R. M. (1994). Practical Bayesian guidelines for
#' phase IIB clinical trials. \emph{Biometrics}, 50, 337--349.
#' \doi{10.2307/2533377}.
#'
#' Lee, J. J., and Liu, D. D. (2008). A predictive probability design for
#' phase II cancer clinical trials. \emph{Clinical Trials}, 5, 93--106.
#' \doi{10.1177/1740774508089279}.
#'
#' Zhao, J., and Mu, R. (2026).
#' Introduction to Bayesian Decision Theory and Its Application
#' in Medical Research. Manuscript in preparation (in Chinese).
#'
#' @seealso \code{\link{posterior_predictive_singlearm_binary}}
#' @examples
#' null_oc <- oc_singlearm_binary(
#'   theta_true = 0.20, n_current = 20, N_total = 40,
#'   ntrial = 1000
#' )
#' alternative_oc <- oc_singlearm_binary(
#'   theta_true = 0.40, n_current = 20, N_total = 40,
#'   ntrial = 1000
#' )
#' null_oc$power_null_post  # Type I error probability at the null rate
#' alternative_oc$power_null_post  # Power at the alternative rate
#' @export

oc_singlearm_binary <- function(
    theta_true,           # True response rate (0.20 under the null; 0.40 under the alternative)
    n_current,            # Planned interim sample size
    N_total,              # Maximum planned sample size
    alpha0 = 0.2,         # First shape parameter of the Beta prior
    beta0 = 0.8,          # Second shape parameter of the Beta prior
    theta0_null = 0.20,   # Historical reference response rate
    gamma_early = 0.98,     # Interim posterior probability threshold for efficacy stopping
    gamma_futility = 0.50,       # Interim posterior probability threshold for futility stopping
    gamma_final = 0.95,     # Final posterior probability threshold for success
    ntrial = 10000,       # Number of Monte Carlo trial simulations
    seed = 666            # Random seed; NULL uses the current RNG state
) {

  # Validate all inputs before changing the RNG state
  .validate_probability(theta_true, "theta_true")
  .validate_sample_sizes(n_current, N_total)
  .validate_prior(alpha0, beta0)
  .validate_probability(theta0_null, "theta0_null")
  .validate_probability(gamma_early, "gamma_early")
  .validate_probability(gamma_futility, "gamma_futility")
  .validate_probability(gamma_final, "gamma_final")
  if (gamma_futility >= gamma_early) {
    stop(
      "`gamma_futility` must be strictly less than `gamma_early`.",
      call. = FALSE
    )
  }
  .validate_count(ntrial, "ntrial", min_value = 1L)
  .validate_seed(seed)

  if (!is.null(seed)) set.seed(seed)

  # Initialize the remaining sample size and outcome counters
  m_remain <- N_total - n_current  # Number of patients remaining after the interim analysis
  count_null_post <- 0             # Number of trials declared successful
  count_early_stop_f <- 0          # Number of trials stopped early for futility
  count_early_stop_e <- 0          # Number of trials stopped early for efficacy

  # Simulate trials using Monte Carlo sampling
  for (itrial in seq_len(ntrial)) {

    # Generate interim response data and update the posterior
    x_obs <- rbinom(1, n_current, theta_true)
    post <- list(
      alpha = alpha0 + x_obs,               # First shape parameter of the posterior
      beta  = beta0 + (n_current - x_obs)  # Second shape parameter of the posterior
    )

    # Compute the posterior probability of exceeding the historical rate
    p_null <- pbeta(theta0_null, post$alpha, post$beta, lower.tail = FALSE)

    # Generate future response data and update the final posterior
    y_fut <- rbinom(1, size = m_remain, prob = theta_true)
    alpha_f <- post$alpha + y_fut
    beta_f <- post$beta + (m_remain - y_fut)

    # Apply the prespecified interim and final decision thresholds
    if (p_null > gamma_early) {

      # Declare early success when the interim efficacy criterion is met
      count_null_post <- count_null_post + 1
      count_early_stop_e <- count_early_stop_e + 1

    } else {

      if (p_null > gamma_futility) {

        # Continue to the final analysis when neither stopping criterion is met
        p_final <- pbeta(theta0_null, alpha_f, beta_f, lower.tail = FALSE)

        # Declare success when the final efficacy criterion is met
        if (p_final > gamma_final) {
          count_null_post <- count_null_post + 1
        }

      } else {

        # Stop early for futility when the continuation criterion is not met
        count_early_stop_f <- count_early_stop_f + 1
      }
    }
  }

  # Return the estimated operating characteristics
  list(
    power_null_post = count_null_post / ntrial,             # Overall success probability (type I error or power)
    early_stop_rate_futility = count_early_stop_f / ntrial, # Probability of early futility stopping
    early_stop_rate_efficacy = count_early_stop_e / ntrial  # Probability of early efficacy stopping
  )
}
