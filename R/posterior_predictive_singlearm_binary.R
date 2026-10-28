#' Posterior and Predictive Probabilities for a Single-Arm Binary Trial
#'
#' Calculate the Beta-binomial posterior distribution and the posterior
#' predictive probabilities of meeting two prespecified final-analysis
#' success criteria after enrolling the remaining patients.
#'
#' @param x_obs Non-negative integer number of observed responses, not
#'   exceeding \code{n_current}.
#' @param n_current Non-negative integer number of evaluated patients at
#'   the current analysis, not exceeding \code{N_total}.
#' @param N_total Non-negative integer maximum planned sample size,
#'   including current patients.
#' @param alpha0,beta0 Finite positive shape parameters of the Beta prior;
#'   their sum must also be finite.
#' @param theta0_null Finite historical reference response rate, between
#'   0 and 1.
#' @param theta0_clin Finite clinically relevant response-rate threshold,
#'   between 0 and 1.
#' @param gamma_final Finite final posterior probability threshold for success
#'   relative to \code{theta0_null}, between 0 and 1.
#' @param gamma_clin_final Finite final posterior probability threshold for
#'   the additional clinical-target criterion relative to
#'   \code{theta0_clin}, between 0 and 1.
#' @param S Positive integer number of Monte Carlo replications.
#' @param seed Integer random seed, or \code{NULL} to use the current random
#'   number generator state. Providing a seed resets the session RNG state.
#'
#' @details
#' Given \code{x_obs} responses among \code{n_current} patients, the updated
#' Beta posterior has shape parameters \code{alpha0 + x_obs} and
#' \code{beta0 + n_current - x_obs}.
#'
#' For each Monte Carlo replication, a response rate is drawn from the
#' current Beta posterior, future binary responses are simulated, and
#' the posterior is updated using the combined data. Both predictive
#' probabilities use the same simulated future data. Final success requires
#' the corresponding posterior probability to be strictly greater than
#' its threshold; equality does not constitute success. These are Monte Carlo
#' estimates, not exact Beta-binomial sums.
#'
#' The clinical-target predictive probability \code{ppos_clin} is an
#' additional summary. It is not part of the trial success rule implemented
#' by \code{\link{oc_singlearm_binary}}.
#'
#' @return A named list containing:
#' \describe{
#'   \item{posterior}{A list with the posterior Beta shape parameters
#'     \code{alpha} and \code{beta}.}
#'   \item{posterior_mean}{Posterior mean of the response rate.}
#'   \item{posterior_sd}{Posterior standard deviation of the response rate.}
#'   \item{p_null}{Posterior probability that the response rate exceeds
#'     \code{theta0_null}.}
#'   \item{ppos_null}{Predictive probability that the final posterior
#'     probability of exceeding \code{theta0_null} is strictly greater than
#'     \code{gamma_final}.}
#'   \item{ppos_clin}{Predictive probability that the final posterior
#'     probability of exceeding \code{theta0_clin} is strictly greater than
#'     \code{gamma_clin_final}.}
#' }
#' @seealso \code{\link{oc_singlearm_binary}}
#' @examples
#' result <- posterior_predictive_singlearm_binary(
#'   x_obs = 9, n_current = 20, N_total = 40, S = 1000
#' )
#' result$p_null
#' result$ppos_null
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
    seed = 666           # Random seed
) {
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

  if (!is.null(seed)) set.seed(seed)
  
  # Update the posterior distribution
  m_remain <- N_total - n_current  # Number of remaining patients
  post <- list(
    alpha = alpha0 + x_obs,                 # First shape parameter of the posterior
    beta  = beta0 + (n_current - x_obs)    # Second shape parameter of the posterior
  )
  
  # Compute the posterior mean and standard deviation
  post_mean <- post$alpha / (post$alpha + post$beta)
  post_sd <- sqrt(post$alpha * post$beta /
                    ((post$alpha + post$beta)^2 * (post$alpha + post$beta + 1)))
  
  # Compute the current posterior probability
  p_null <- pbeta(theta0_null, post$alpha, post$beta, lower.tail = FALSE)
  
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
  
  # Return posterior summaries and decision probabilities
  list(
    posterior = post,           # Shape parameters of the Beta posterior
    posterior_mean = post_mean, # Posterior mean
    posterior_sd = post_sd,     # Posterior standard deviation
    p_null = p_null,            # Posterior probability of exceeding the historical rate
    ppos_null = ppos_null,      # Predictive probability of meeting the historical-rate success criterion
    ppos_clin = ppos_clin       # Predictive probability of meeting the clinical-target success criterion
  )
}
