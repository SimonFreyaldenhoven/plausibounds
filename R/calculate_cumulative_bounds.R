#' Calculate Cumulative Bounds (Average Treatment Effect)
#'
#' This function calculates bounds for the average treatment effect (ATE) from
#' a vector of estimates. It can calculate both pointwise and simultaneous (sup-t) bounds.
#' Supports pre-treatment periods for event study designs.
#'
#' @param estimates A vector of point estimates. If preperiods > 0, the first preperiods
#'   elements are pre-treatment estimates, followed by post-treatment estimates.
#' @param var The variance-covariance matrix of the estimates
#' @param alpha Significance level (default: 0.05)
#' @param preperiods Number of pre-treatment periods (default: 0). Period 0 is assumed
#'   to be normalized and not included in estimates.
#' @param include_pointwise Whether to include pointwise bounds (default: TRUE)
#' @param include_supt Whether to include sup-t bounds (default: TRUE)
#'
#' @return A list containing:
#'   \item{cumulative_bounds}{A data frame with columns for horizon (event time), coefficients, and bounds}
#'   \item{ate}{Average treatment effect with standard error}
#'   \item{metadata}{A list with metadata about the calculation}
#'
#' @examples
#' # Example with constant estimates and IID errors
#' data(estimates_constant_iid)
#' data(var_constant_iid)
#' cumul_bounds <- calculate_cumulative_bounds(estimates_constant_iid, var_constant_iid)
#'
#' # Example with wiggly estimates and strong correlation
#' data(estimates_wiggly_strong_corr)
#' data(var_wiggly_strong_corr)
#' cumul_bounds_complex <- calculate_cumulative_bounds(
#'   estimates_wiggly_strong_corr,
#'   var_wiggly_strong_corr
#' )
#'
#' @export
calculate_cumulative_bounds <- function(estimates, var, alpha = 0.05,
                                       preperiods = 0,
                                       include_pointwise = TRUE, include_supt = TRUE) {
  # Check inputs
  if (!is.numeric(estimates) || !is.vector(estimates)) {
    stop("estimates must be a numeric vector")
  }
  if (any(is.na(estimates))) {
    stop("estimates cannot contain NA values")
  }
  if (any(!is.finite(estimates))) {
    stop("estimates cannot contain infinite values")
  }
  if (!is.matrix(var) || nrow(var) != length(estimates) || ncol(var) != length(estimates)) {
    stop("var must be a square matrix with dimensions matching the length of estimates")
  }
  if (any(is.na(var))) {
    stop("var cannot contain NA values")
  }
  if (any(!is.finite(var))) {
    stop("var cannot contain infinite values")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number between 0 and 1")
  }
  if (!is.numeric(preperiods) || length(preperiods) != 1 || preperiods < 0 || preperiods != floor(preperiods)) {
    stop("preperiods must be a non-negative integer")
  }
  if (preperiods >= length(estimates)) {
    stop("preperiods must be less than the length of estimates")
  }

  # Store full data
  estimatesAll <- estimates
  varAll <- var
  n_all <- length(estimatesAll)

  # Extract post-period data
  if (preperiods > 0) {
    estimates_pre <- estimatesAll[1:preperiods]
    var_pre <- varAll[1:preperiods, 1:preperiods, drop = FALSE]
    estimates_post <- estimatesAll[(preperiods + 1):n_all]
    var_post <- varAll[(preperiods + 1):n_all, (preperiods + 1):n_all, drop = FALSE]
  } else {
    estimates_post <- estimates
    var_post <- var
  }

  p <- length(estimates_post)  # p is the number of post-periods

  # Calculate Average Treatment Effect (ATE) with standard error
  ate_hat <- mean(estimates_post)
  se_ate <- sqrt(sum(var_post)) / p
  ate <- c(estimate = ate_hat, se = se_ate)

  # Calculate bounds using simple normal CI
  z_crit <- stats::qnorm(1 - alpha / 2)
  lb <- ate_hat - z_crit * se_ate
  ub <- ate_hat + z_crit * se_ate

  # Create data frame with event time
  if (preperiods > 0) {
    # Post-period bounds
    bounds_df_post <- data.frame(
      horizon = 1:p,
      coef = estimates_post,
      lower = lb,
      upper = ub
    )
    # Pre-period (no ATE bounds for pre-periods)
    bounds_df_pre <- data.frame(
      horizon = -preperiods:-1,
      coef = estimates_pre,
      lower = NA_real_,
      upper = NA_real_
    )
    bounds_df <- rbind(bounds_df_pre, bounds_df_post)
  } else {
    bounds_df <- data.frame(
      horizon = 1:p,
      coef = estimates_post,
      lower = lb,
      upper = ub
    )
  }

  # Calculate width
  width <- ub - lb

  # Return list with bounds data frame and metadata
  result <- list(
    cumulative_bounds = bounds_df,
    ate = ate
  )

  # Calculate pointwise bounds if requested (for ALL periods)
  if (include_pointwise) {
    if (preperiods > 0) {
      pw_pre <- calculate_pointwise_bounds(estimates_pre, var_pre, alpha)
      pw_post <- calculate_pointwise_bounds(estimates_post, var_post, alpha)
      result$pointwise_bounds <- list(
        lower = c(pw_pre$lower, pw_post$lower),
        upper = c(pw_pre$upper, pw_post$upper)
      )
    } else {
      result$pointwise_bounds <- calculate_pointwise_bounds(estimates_post, var_post, alpha)
    }
  }

  # Calculate sup-t bounds if requested (for ALL periods)
  if (include_supt) {
    if (preperiods > 0) {
      # Compute sup-t critical value using ALL periods
      stdv_all <- sqrt(diag(varAll))
      d_nonzero <- stdv_all > .Machine$double.eps
      cV_all <- diag(1 / stdv_all[d_nonzero]) %*% varAll[d_nonzero, d_nonzero] %*% diag(1 / stdv_all[d_nonzero])

      set.seed(42)
      kk <- 10000
      eig_c <- eigen(cV_all, symmetric = TRUE)
      Corrmat_sqrt <- eig_c$vectors %*% diag(sqrt(pmax(eig_c$values, 0))) %*% t(eig_c$vectors)
      t_stat <- abs(matrix(stats::rnorm(kk * sum(d_nonzero)), nrow = kk) %*% Corrmat_sqrt)
      supt_critval <- as.numeric(stats::quantile(apply(t_stat, 1, max), 1 - alpha))

      supt_lower_all <- estimatesAll - supt_critval * sqrt(diag(varAll))
      supt_upper_all <- estimatesAll + supt_critval * sqrt(diag(varAll))
      result$supt_bounds <- list(
        lower = supt_lower_all,
        upper = supt_upper_all
      )
    } else {
      result$supt_bounds <- calculate_supt_bounds(estimates_post, var_post, alpha)
    }
  }

  result$metadata <- list(
    alpha = alpha,
    preperiods = preperiods,
    width = width,
    lb = lb,
    ub = ub
  )

  class(result) <- c("cumulative_bounds", "plausible_bounds_result")
  return(result)
}