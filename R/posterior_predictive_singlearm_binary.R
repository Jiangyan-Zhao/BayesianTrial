#' Posterior and Predictive Probabilities for a Single-Arm Phase II Trial
#'
#' Calculates posterior summaries of the treatment response rate and
#' predictive probabilities of meeting prespecified final-analysis efficacy
#' criteria in single-arm phase II trials with binary endpoints.
#'
#' @param x_obs Number of observed responses among \code{n_current} patients;
#'   a non-negative integer no greater than \code{n_current}.
#' @param n_current Number of patients with evaluable outcomes at the
#'   current analysis; a non-negative integer no greater than \code{N_total}.
#' @param N_total Planned maximum number of patients with evaluable outcomes;
#'   a non-negative integer at least as large as \code{n_current}.
#' @param alpha0,beta0 Positive shape parameters of the Beta prior for the
#'   response rate.
#' @param theta0_null Prespecified historical reference or unacceptable
#'   response rate; a probability between 0 and 1.
#' @param theta0_clin Clinically meaningful target response rate for a
#'   supplementary predictive assessment; a probability between 0 and 1.
#' @param gamma_final Final posterior probability threshold for trial success
#'   relative to \code{theta0_null}; a probability between 0 and 1.
#' @param gamma_clin_final Final posterior probability threshold for the
#'   supplementary criterion relative to \code{theta0_clin}; a probability
#'   between 0 and 1.
#' @param S Positive integer number of Monte Carlo replications used to
#'   estimate predictive probabilities.
#' @param type Results to return: \code{"posterior"} (default),
#'   \code{"predictive"}, or \code{"both"}.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A named list containing the components specified by \code{type}.
#' The default returns posterior summaries; \code{"predictive"} returns
#' predictive probabilities; \code{"both"} returns all components:
#' \describe{
#'   \item{posterior}{Shape parameters \code{alpha} and \code{beta} of the
#'     current Beta posterior.
#'   }
#'   \item{posterior_mean}{Posterior mean of the response rate.
#'   }
#'   \item{posterior_sd}{Posterior standard deviation of the response rate.
#'   }
#'   \item{p_null}{Posterior probability that the response rate exceeds
#'     \code{theta0_null}.
#'   }
#'   \item{ppos_null}{Predictive probability of meeting the final efficacy
#'     criterion defined by \code{theta0_null} and \code{gamma_final}.
#'   }
#'   \item{ppos_clin}{Predictive probability of meeting the supplementary
#'     criterion defined by \code{theta0_clin} and \code{gamma_clin_final}.
#'   }
#' }
#'
#' Predictive probabilities are estimated by Monte Carlo simulation,
#' conditional on the observed data and completion of the planned
#' sample size \code{N_total}.
#'
#' @references
#' Zhao, J., and Mu, R. (2026).
#' Introduction to Bayesian Decision Theory and Its Application
#' in Medical Research. Manuscript in preparation (in Chinese).
#'
#' Thall, P. F., and Simon, R. M. (1994). Practical Bayesian guidelines for
#' phase IIB clinical trials. \emph{Biometrics}, 50, 337--349.
#' \doi{10.2307/2533377}.
#'
#' Lee, J. J., and Liu, D. D. (2008). A predictive probability design for
#' phase II cancer clinical trials. \emph{Clinical Trials}, 5, 93--106.
#' \doi{10.1177/1740774508089279}.
#'
#' @seealso \code{\link{oc_singlearm_binary}}
#' @examples
#' posterior <- posterior_predictive_singlearm_binary(
#'   x_obs = 9, n_current = 20, N_total = 40
#' )
#' posterior$p_null
#'
#' predictive <- posterior_predictive_singlearm_binary(
#'   x_obs = 9, n_current = 20, N_total = 40,
#'   type = "predictive", S = 1000
#' )
#' predictive$ppos_null
#' @importFrom stats pbeta rbeta rbinom
#' @export

posterior_predictive_singlearm_binary <- function(
    x_obs,               # Number of observed responses
    n_current,           # Number of evaluated patients at the current analysis
    N_total,             # Maximum planned sample size
    alpha0 = 0.2,        # First shape parameter of the Beta prior
    beta0 = 0.8,         # Second shape parameter of the Beta prior
    theta0_null = 0.20,  # Historical reference response rate
    theta0_clin = 0.40,  # Prespecified clinically relevant response-rate threshold
    gamma_final = 0.95,   # Final posterior probability threshold for exceeding the historical rate
    gamma_clin_final = 0.90,   # Final posterior probability threshold for exceeding the clinical target rate
    S = 10000,           # Number of Monte Carlo replications
    type = c("posterior", "predictive", "both"),  # Requested output
    seed = 666          # Random seed
) {
  # Select the requested output before validating the inputs
  type <- match.arg(type)

  # Validate all inputs before changing the RNG state
  .validate_sample_sizes(n_current, N_total)
  .validate_count(x_obs, "x_obs")
  if (x_obs > n_current) {
    stop("`x_obs` must not exceed `n_current`.", call. = FALSE)
  }
  .validate_prior(alpha0, beta0)
  .validate_probability(theta0_null, "theta0_null")
  .validate_probability(theta0_clin, "theta0_clin")
  .validate_probability(gamma_final, "gamma_final")
  .validate_probability(gamma_clin_final, "gamma_clin_final")
  .validate_count(S, "S", min_value = 1L)
  .validate_seed(seed)

  # Update the posterior distribution
  m_remain <- N_total - n_current  # Number of remaining patients
  post <- list(
    alpha = alpha0 + x_obs,                 # First shape parameter of the posterior
    beta  = beta0 + (n_current - x_obs)    # Second shape parameter of the posterior
  )

  # Return exact posterior summaries without running predictive simulations
  if (type != "predictive") {
    post_mean <- post$alpha / (post$alpha + post$beta)
    post_sd <- sqrt(post$alpha * post$beta /
                      ((post$alpha + post$beta)^2 * (post$alpha + post$beta + 1)))
    p_null <- pbeta(theta0_null, post$alpha, post$beta, lower.tail = FALSE)

    posterior_result <- list(
      posterior = post,
      posterior_mean = post_mean,
      posterior_sd = post_sd,
      p_null = p_null
    )

    if (type == "posterior") {
      return(posterior_result)
    }
  }

  # Set the random seed only when predictive simulation is required
  if (!is.null(seed)) set.seed(seed)

  # Estimate predictive probabilities for both final success criteria
  success_null <- numeric(S)
  success_clin <- numeric(S)

  for (s in seq_len(S)) {

    # Draw a response rate from the current posterior
    pi_s <- rbeta(1, post$alpha, post$beta)

    # Simulate the number of responses among the remaining patients
    y_fut <- rbinom(1, size = m_remain, prob = pi_s)

    # Update the final posterior using current and simulated future responses
    alpha_f <- post$alpha + y_fut
    beta_f <- post$beta + (m_remain - y_fut)

    # Compute the final posterior probabilities
    p_final_null <- pbeta(theta0_null, alpha_f, beta_f, lower.tail = FALSE)
    p_final_clin <- pbeta(theta0_clin, alpha_f, beta_f, lower.tail = FALSE)

    # Evaluate the prespecified final success criteria
    success_null[s] <- as.numeric(p_final_null > gamma_final)
    success_clin[s] <- as.numeric(p_final_clin > gamma_clin_final)
  }

  # Estimate predictive probabilities as proportions of successful replications
  ppos_null <- mean(success_null)
  ppos_clin <- mean(success_clin)

  # Return only the requested components, retaining the original names
  predictive_result <- list(
    ppos_null = ppos_null,
    ppos_clin = ppos_clin
  )

  if (type == "predictive") {
    return(predictive_result)
  }

  c(posterior_result, predictive_result)
}
