#' Posterior and Predictive Probabilities for a Single-Arm Phase II Trial with a Binary Endpoint
#'
#' Summarize posterior evidence of treatment activity from observed binary
#' efficacy outcomes in a single-arm phase II trial, and estimate the
#' predictive probabilities of meeting prespecified posterior efficacy
#' criteria at the planned final analysis.
#'
#' @param x_obs Number of patients with an observed response among the \code{n_current}
#'   patients with evaluable outcomes; an integer from 0 to \code{n_current}.
#' @param n_current Number of patients with evaluable binary efficacy outcomes at the current
#'   analysis; a non-negative integer not exceeding \code{N_total}.
#' @param N_total Planned maximum number of evaluable patients, including those assessed at
#'   the current analysis; a non-negative integer at least as large as
#'   \code{n_current}.
#' @param alpha0,beta0 Positive, finite shape parameters of the Beta prior for the response
#'   rate. Their sum must also be finite.
#' @param theta0_null Prespecified unacceptable or historical reference response rate used as
#'   the primary efficacy benchmark; a probability between 0 and 1
#'   (inclusive).
#' @param theta0_clin Prespecified clinically meaningful target response rate for a
#'   supplementary predictive assessment; a probability between 0 and 1
#'   (inclusive). This target is not an additional requirement for the primary
#'   trial success decision.
#' @param gamma_final Posterior probability cutoff for declaring a positive trial result at the
#'   final analysis relative to \code{theta0_null}; a probability between 0
#'   and 1 (inclusive).
#' @param gamma_clin_final Posterior probability cutoff for the supplementary final-analysis
#'   criterion relative to \code{theta0_clin}; a probability between 0 and 1
#'   (inclusive).
#' @param S Positive integer number of Monte Carlo replications used to estimate the
#'   predictive probabilities.
#' @param seed Integer random seed, or \code{NULL} to use the current random number
#'   generator state. A non-\code{NULL} value calls \code{set.seed()} and
#'   changes the session's random number generator state.
#'
#' @details
#' Let \eqn{\theta} denote the unknown treatment response rate, with a
#' \eqn{\mathrm{Beta}(\alpha_0,\beta_0)} prior. Given \eqn{x} responses among
#' \eqn{n} patients with evaluable outcomes, its posterior distribution is
#' \eqn{\mathrm{Beta}(\alpha_0+x,\beta_0+n-x)}. The quantity \code{p_null} is
#' the current posterior probability that \eqn{\theta} exceeds the reference
#' rate \code{theta0_null}.
#'
#' To assess the prospect of a positive result after completing the planned
#' enrollment, the function simulates responses for the \code{N_total -
#' n_current} remaining patients from the Beta-binomial posterior predictive
#' distribution. In each replication, the final posterior is updated using the
#' observed and simulated responses, and the prespecified final-analysis
#' criterion is evaluated. The proportion of replications meeting the
#' criterion is reported as \code{ppos_null}. Both predictive assessments use
#' the same simulated future outcomes.
#'
#' The predictive probability is conditional on the current data and on
#' reaching \code{N_total} evaluable patients. It is not the frequentist
#' statistical power of the trial design and does not incorporate interim
#' stopping decisions. Final-analysis success requires the posterior
#' probability to be strictly greater than the corresponding cutoff; equality
#' does not satisfy the criterion.
#'
#' The additional quantity \code{ppos_clin} evaluates the probability of
#' meeting a separate posterior criterion for the clinically meaningful target
#' rate \code{theta0_clin}. It is a supplementary measure and is not included
#' in the trial success rule evaluated by \code{\link{oc_singlearm_binary}}.
#' Predictive probabilities are estimated by Monte Carlo simulation rather
#' than exact Beta-binomial summation.
#'
#' @return A named list containing:
#' \describe{
#'   \item{posterior}{A list containing the shape parameters \code{alpha} and \code{beta}
#'     of the current Beta posterior.}
#'   \item{posterior_mean}{Posterior mean of the treatment response rate at the current
#'     analysis.}
#'   \item{posterior_sd}{Posterior standard deviation of the treatment response rate at the
#'     current analysis.}
#'   \item{p_null}{Current posterior probability that the response rate exceeds
#'     \code{theta0_null}.}
#'   \item{ppos_null}{Posterior predictive probability of satisfying the final efficacy
#'     criterion based on \code{theta0_null} and \code{gamma_final},
#'     conditional on the observed data and completing the planned sample
#'     size.}
#'   \item{ppos_clin}{Supplementary posterior predictive probability of satisfying the
#'     final criterion based on \code{theta0_clin} and
#'     \code{gamma_clin_final}.}
#' }
#'
#' @references
#' Lee, J. J., and Liu, D. D. (2008). A predictive probability design for
#' phase II cancer clinical trials. \emph{Clinical Trials}, 5, 93--106.
#' \doi{10.1177/1740774508089279}.
#'
#' Zhou, H., Lee, J. J., and Yuan, Y. (2017). BOP2: Bayesian optimal design
#' for phase II clinical trials with simple and complex endpoints.
#' \emph{Statistics in Medicine}, 36, 3302--3314. \doi{10.1002/sim.7338}.
#'
#' @seealso \code{\link{oc_singlearm_binary}}
#' @examples
#' interim <- posterior_predictive_singlearm_binary(
#'   x_obs = 9, n_current = 20, N_total = 40,
#'   theta0_null = 0.20, gamma_final = 0.95, S = 1000
#' )
#' interim$p_null
#' interim$ppos_null
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
