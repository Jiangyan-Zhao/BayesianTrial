#' Posterior and Predictive Probabilities for a Single-Arm Binary Trial
#'
#' Calculate the Beta-binomial posterior distribution and the posterior
#' predictive probabilities of meeting two prespecified final-analysis
#' success criteria after enrolling the remaining patients.
#'
#' @param x_obs Number of observed responses at the current analysis.
#' @param n_current Number of evaluated patients at the current analysis.
#' @param N_total Maximum planned sample size, including current patients.
#' @param alpha0,beta0 Positive shape parameters of the Beta prior.
#' @param theta0_null Historical reference response rate, between 0 and 1.
#' @param theta0_clin Clinically relevant response-rate threshold, between 0 and 1.
#' @param gamma_null Final posterior probability threshold for success
#'   relative to \code{theta0_null}, between 0 and 1.
#' @param gamma_clin Final posterior probability threshold for success
#'   relative to \code{theta0_clin}, between 0 and 1.
#' @param S Positive integer number of Monte Carlo replications.
#' @param seed Integer random seed, or \code{NULL} to use the current random
#'   number generator state. Providing a seed resets the session RNG state.
#'
#' @details
#' Let \eqn{a = alpha0 + x_obs} and \eqn{b = beta0 + n_current - x_obs}.
#' For every Monte Carlo replication, draw a response rate from Beta(a, b),
#' generate the remaining responses, update the posterior, and evaluate the
#' two final success criteria. Both predictive probabilities use the same
#' simulated future data. The final success comparisons use \code{>=},
#' matching the supplied source script. The predictive probabilities are
#' Monte Carlo estimates, not exact Beta-binomial sums.
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
#'     probability of exceeding \code{theta0_null} is at least
#'     \code{gamma_null}.}
#'   \item{ppos_clin}{Predictive probability that the final posterior
#'     probability of exceeding \code{theta0_clin} is at least
#'     \code{gamma_clin}.}
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
    gamma_null = 0.95,   # Final posterior probability threshold for exceeding the historical rate
    gamma_clin = 0.90,   # Final posterior probability threshold for exceeding the clinical target rate
    S = 10000,           # Number of Monte Carlo replications
    seed = 666           # Random seed
) {
  # Validate inputs
  stopifnot(
    length(x_obs) == 1L, x_obs >= 0, x_obs <= n_current,
    x_obs == as.integer(x_obs),
    n_current >= 0, N_total >= n_current,
    n_current == as.integer(n_current), N_total == as.integer(N_total),
    alpha0 > 0, beta0 > 0, S >= 1, S == as.integer(S),
    theta0_null >= 0, theta0_null <= 1,
    theta0_clin >= 0, theta0_clin <= 1,
    gamma_null >= 0, gamma_null <= 1,
    gamma_clin >= 0, gamma_clin <= 1
  )
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
  p_null <- 1 - pbeta(theta0_null, post$alpha, post$beta)
  
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
    p_final_null <- 1 - pbeta(theta0_null, alpha_f, beta_f)
    p_final_clin <- 1 - pbeta(theta0_clin, alpha_f, beta_f)
    
    # Evaluate the prespecified final success criteria
    success_null[s] <- as.numeric(p_final_null >= gamma_null)
    success_clin[s] <- as.numeric(p_final_clin >= gamma_clin)
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
