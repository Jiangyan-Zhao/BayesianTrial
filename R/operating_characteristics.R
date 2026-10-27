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
#' @param gamma_s_1 Interim posterior probability threshold for early
#'   efficacy stopping.
#' @param gamma_c Interim posterior probability threshold for futility
#'   stopping. This applies to the \emph{current posterior probability},
#'   not to a predictive probability.
#' @param gamma_f_1 Final posterior probability threshold for success.
#' @param ntrial Positive integer number of simulated trials.
#' @param seed Integer random seed, or \code{NULL} to use the current random
#'   number generator state. Providing a seed resets the session RNG state.
#'
#' @details
#' For each simulated trial, the interim posterior probability that the
#' true response rate exceeds \code{theta0_null} is evaluated. The trial
#' declares early success when this probability is strictly greater than
#' \code{gamma_s_1}; stops for futility when it is less than or equal to
#' \code{gamma_c}; otherwise the final probability is compared strictly
#' against \code{gamma_f_1}. The routine draws future data even for trials
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
#' @seealso \code{\link{calc_posterior_predictive}}
#' @examples
#' null_oc <- calc_power(theta_true = 0.20, n_current = 20,
#'                       N_total = 40, ntrial = 1000)
#' alternative_oc <- calc_power(theta_true = 0.40, n_current = 20,
#'                              N_total = 40, ntrial = 1000)
#' null_oc
#' alternative_oc
#' @export

calc_power <- function(
    theta_true,           # 真实缓解率（0.20 对应原假设，0.40 对应备择假设）
    n_current,            # 期中分析时的样本量
    N_total,              # 计划总样本量
    alpha0 = 0.2,         # Beta 先验分布的第一个形状参数
    beta0 = 0.8,          # Beta 先验分布的第二个形状参数
    theta0_null = 0.20,   # 历史对照缓解率
    gamma_s_1 = 0.98,     # 期中有效性提前停止的后验概率阈值
    gamma_c = 0.50,       # 期中无效性提前停止的后验概率阈值
    gamma_f_1 = 0.95,     # 最终有效性判定的后验概率阈值
    ntrial = 10000,       # 蒙特卡洛模拟的试验次数
    seed = 666            # 随机数种子，NULL 表示不重新设置
) {
  
  # 检查输入参数的有效性
  stopifnot(
    theta_true >= 0, theta_true <= 1,
    n_current >= 0, N_total >= n_current,
    n_current == as.integer(n_current), N_total == as.integer(N_total),
    alpha0 > 0, beta0 > 0, theta0_null >= 0, theta0_null <= 1,
    gamma_s_1 >= 0, gamma_s_1 <= 1,
    gamma_c >= 0, gamma_c <= 1,
    gamma_f_1 >= 0, gamma_f_1 <= 1,
    ntrial >= 1, ntrial == as.integer(ntrial)
  )
  if (!is.null(seed)) set.seed(seed)
  
  # 初始化样本量及计数器
  m_remain <- N_total - n_current  # 期中分析后的剩余样本量
  count_null_post <- 0             # 宣布试验成功的次数
  count_early_stop_f <- 0          # 因无效性提前停止的次数
  count_early_stop_e <- 0          # 因有效性提前停止的次数
  
  # 蒙特卡洛模拟
  for (itrial in seq_len(ntrial)) {
    
    # 生成期中分析数据并更新后验分布
    x_obs <- rbinom(1, n_current, theta_true)
    post <- list(
      alpha = alpha0 + x_obs,               # 后验分布的第一个形状参数
      beta  = beta0 + (n_current - x_obs)  # 后验分布的第二个形状参数
    )
    
    # 计算当前缓解率超过历史对照水平的后验概率
    p_null <- 1 - pbeta(theta0_null, post$alpha, post$beta)
    
    # 生成未来数据并更新最终后验分布
    y_fut <- rbinom(1, size = m_remain, prob = theta_true)
    alpha_f <- post$alpha + y_fut
    beta_f <- post$beta + (m_remain - y_fut)
    
    # 根据预设阈值执行期中及最终决策
    if (p_null > gamma_s_1) {
      
      # 期中达到有效性标准，提前宣布成功
      count_null_post <- count_null_post + 1
      count_early_stop_e <- count_early_stop_e + 1
      
    } else {
      
      if (p_null > gamma_c) {
        
        # 期中未达到提前停止标准，继续至最终分析
        p_final <- 1 - pbeta(theta0_null, alpha_f, beta_f)
        
        # 最终达到有效性标准，宣布试验成功
        if (p_final > gamma_f_1) {
          count_null_post <- count_null_post + 1
        }
        
      } else {
        
        # 期中未达到继续试验标准，因无效性提前停止
        count_early_stop_f <- count_early_stop_f + 1
      }
    }
  }
  
  # 返回各项操作特征
  list(
    power_null_post = count_null_post / ntrial,             # 试验成功率（一类错误率或检验功效）
    early_stop_rate_futility = count_early_stop_f / ntrial, # 无效性提前停止概率
    early_stop_rate_efficacy = count_early_stop_e / ntrial  # 有效性提前停止概率
  )
}
