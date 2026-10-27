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
#' @seealso \code{\link{calc_power}}
#' @examples
#' result <- calc_posterior_predictive(
#'   x_obs = 9, n_current = 20, N_total = 40, S = 1000
#' )
#' result$p_null
#' result$ppos_null
#' @importFrom stats pbeta rbeta rbinom
#' @export

calc_posterior_predictive <- function(
    x_obs,               # 当前观察到的缓解例数
    n_current,           # 当前样本量（期中分析时）
    N_total,             # 计划总样本量
    alpha0 = 0.2,        # Beta 先验分布的第一个形状参数
    beta0 = 0.8,         # Beta 先验分布的第二个形状参数
    theta0_null = 0.20,  # 历史对照缓解率
    theta0_clin = 0.40,  # 预设的临床目标缓解率
    gamma_null = 0.95,   # 最终判定优于历史对照的后验概率阈值
    gamma_clin = 0.90,   # 最终判定达到临床目标的后验概率阈值
    S = 10000,           # 蒙特卡洛模拟次数
    seed = 666           # 随机数种子
) {
  # 检查输入参数的有效性
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
  
  # 更新后验分布
  m_remain <- N_total - n_current  # 剩余样本量
  post <- list(
    alpha = alpha0 + x_obs,                 # 后验分布的第一个形状参数
    beta  = beta0 + (n_current - x_obs)    # 后验分布的第二个形状参数
  )
  
  # 计算后验均值和标准差
  post_mean <- post$alpha / (post$alpha + post$beta)
  post_sd <- sqrt(post$alpha * post$beta /
                    ((post$alpha + post$beta)^2 * (post$alpha + post$beta + 1)))
  
  # 计算当前后验概率
  p_null <- 1 - pbeta(theta0_null, post$alpha, post$beta)
  
  # 通过蒙特卡洛模拟计算两种最终成功标准的预测成功概率
  success_null <- numeric(S)
  success_clin <- numeric(S)
  
  for (s in seq_len(S)) {
    
    # 从当前后验分布中抽取真实缓解率
    pi_s <- rbeta(1, post$alpha, post$beta)
    
    # 模拟剩余患者的缓解例数
    y_fut <- rbinom(1, size = m_remain, prob = pi_s)
    
    # 根据当前及未来数据更新最终后验分布
    alpha_f <- post$alpha + y_fut
    beta_f <- post$beta + (m_remain - y_fut)
    
    # 计算最终后验概率
    p_final_null <- 1 - pbeta(theta0_null, alpha_f, beta_f)
    p_final_clin <- 1 - pbeta(theta0_clin, alpha_f, beta_f)
    
    # 判断是否达到预设的最终成功标准
    success_null[s] <- as.numeric(p_final_null >= gamma_null)
    success_clin[s] <- as.numeric(p_final_clin >= gamma_clin)
  }
  
  # 计算达到最终成功标准的模拟比例
  ppos_null <- mean(success_null)
  ppos_clin <- mean(success_clin)
  
  # 返回后验分布及相关决策概率
  list(
    posterior = post,           # Beta 后验分布的形状参数
    posterior_mean = post_mean, # 后验均值
    posterior_sd = post_sd,     # 后验标准差
    p_null = p_null,            # 缓解率超过历史对照水平的后验概率
    ppos_null = ppos_null,      # 最终判定优于历史对照的预测成功概率
    ppos_clin = ppos_clin       # 最终判定达到临床目标的预测成功概率
  )
}
