#' Operating Characteristics of a Two-Stage Single-Arm Binary Trial
#'
#' Simulate the complete two-stage trial under a specified true response rate
#' and estimate the overall success probability and interim stopping rates.
#'
#' @param theta_true True response probability used to generate trial data.
#' @param n_current Planned interim sample size.
#' @param N_total Maximum planned sample size.
#' @param alpha0,beta0 Positive shape parameters of the Beta prior.
#' @param theta0_null Historical reference response rate, between 0 and 1.
#' @param gamma_early Interim posterior probability threshold for early
#'   efficacy stopping.
#' @param gamma_futility Interim posterior probability threshold for futility
#'   stopping. This applies to the \emph{current posterior probability},
#'   not to a predictive probability.
#' @param gamma_final Final posterior probability threshold for success.
#' @param ntrial Positive integer number of simulated trials.
#' @param seed Integer random seed, or \code{NULL} to use the current random
#'   number generator state. Providing a seed resets the session RNG state.
#'
#' @details
#' For each simulated trial, the interim posterior probability that the
#' true response rate exceeds \code{theta0_null} is evaluated. The trial
#' declares early success when this probability is strictly greater than
#' \code{gamma_early}; stops for futility when it is less than or equal to
#' \code{gamma_futility}; otherwise the final probability is compared strictly
#' against \code{gamma_final}. The routine draws future data even for trials
#' that stop early, preserving the random number sequence of the source
#' script; future data do not affect decisions in those trials.
#'
#' The reported \code{power_null_post} is the overall probability of
#' declaring trial success. Under a null response rate, it estimates the
#' type I error; under an alternative, it estimates statistical power.
#'
#' @return A named list containing:
#' \describe{
#'   \item{power_null_post}{Estimated overall trial success probability.}
#'   \item{early_stop_rate_futility}{Estimated probability of interim
#'     futility stopping.}
#'   \item{early_stop_rate_efficacy}{Estimated probability of interim
#'     efficacy stopping.}
#' }
#' @seealso \code{\link{posterior_predictive_singlearm_binary}}
#' @examples
#' null_oc <- oc_singlearm_binary(theta_true = 0.20, n_current = 20,
#'                       N_total = 40, ntrial = 1000)
#' alternative_oc <- oc_singlearm_binary(theta_true = 0.40, n_current = 20,
#'                              N_total = 40, ntrial = 1000)
#' null_oc
#' alternative_oc
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
  
  # Validate inputs
  stopifnot(
    theta_true >= 0, theta_true <= 1,
    n_current >= 0, N_total >= n_current,
    n_current == as.integer(n_current), N_total == as.integer(N_total),
    alpha0 > 0, beta0 > 0, theta0_null >= 0, theta0_null <= 1,
    gamma_early >= 0, gamma_early <= 1,
    gamma_futility >= 0, gamma_futility <= 1,
    gamma_final >= 0, gamma_final <= 1,
    gamma_futility < gamma_early,
    ntrial >= 1, ntrial == as.integer(ntrial)
  )
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
