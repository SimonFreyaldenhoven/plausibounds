# Internal utility functions

# Calculate sup-t bands
bands_plugin <- function(delta, var, p, nsim = 2000, level = 0.95) {
  # Simulate from N(0, var), calculate sup-t bands, return LB and UB
  rd <- MASS::mvrnorm(nsim, mu = rep(0, p), Sigma = var)
  std_devs <- sqrt(diag(var))
  rd_standardized <- rd / matrix(std_devs, nrow = nsim, ncol = length(std_devs), byrow = TRUE)
  sup_t <- as.numeric(stats::quantile(apply(abs(rd_standardized), 1, max), level))
  
  list(
    LB = delta - sup_t * sqrt(diag(var)),
    UB = delta + sup_t * sqrt(diag(var)),
    sup_t = sup_t
  )
}

#' Calculate pointwise bounds
#'
#' @param estimates A vector of point estimates
#' @param var The variance-covariance matrix of the estimates
#' @param alpha Significance level
#'
#' @return A list containing lower bounds, upper bounds, and critical value
#' @keywords internal
calculate_pointwise_bounds <- function(estimates, var, alpha) {
  pointwise_critical <- stats::qnorm(1 - alpha/2)
  pointwise_lower <- estimates - pointwise_critical * sqrt(diag(var))
  pointwise_upper <- estimates + pointwise_critical * sqrt(diag(var))
  
  list(
    lower = pointwise_lower,
    upper = pointwise_upper,
    critval = pointwise_critical
  )
}

#' Calculate sup-t bounds
#'
#' @param estimates A vector of point estimates
#' @param var The variance-covariance matrix of the estimates
#' @param alpha Significance level
#'
#' @return A list containing lower bounds, upper bounds, and critical value
#' @keywords internal
calculate_supt_bounds <- function(estimates, var, alpha) {
  p <- length(estimates)
  result <- bands_plugin(estimates, var, p, nsim = 2000, level = 1 - alpha)
  
  list(
    lower = result$LB,
    upper = result$UB,
    critval = result$sup_t
  )
}

# Declare global variables used in dplyr/ggplot2 operations
utils::globalVariables(c(
  "horizon", "lower", "upper", "coef", "surrogate", 
  "restricted_lower", "restricted_upper", "Var1"
))
