# Internal helper function for Wald bounds
Wald_bounds <- function(dhat, Vhat, alpha, df = 1) {
  h <- length(dhat)
  critval <- stats::qchisq(1 - alpha, df)
  
  lambda1 <- sqrt(sum(Vhat) / (4 * critval))
  lambda2 <- -sqrt(sum(Vhat) / (4 * critval))
  
  UB <- t(dhat) + (1 / (2 * lambda1)) * (t(rep(1, h)) %*% Vhat)
  LB <- t(dhat) + (1 / (2 * lambda2)) * (t(rep(1, h)) %*% Vhat)
  
  return(list(LB = t(LB), UB = t(UB)))
}

#' Calculate Cumulative Bounds
#'
#' This function calculates cumulative bounds for a vector of estimates.
#' It can calculate both pointwise and simultaneous (sup-t) bounds.
#'
#' @param estimates A vector of point estimates
#' @param var The variance-covariance matrix of the estimates
#' @param alpha Significance level (default: 0.05)
#' @param include_pointwise Whether to include pointwise bounds (default: TRUE)
#' @param include_supt Whether to include sup-t bounds (default: TRUE)
#'
#' @return A list containing:
#'   \item{bounds}{A data frame with columns for horizon, coefficients, and bounds}
#'   \item{type}{The type of bounds ("cumulative")}
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
  
  p <- length(estimates)

  # Calculate Wald bounds for cumulative effect
  wald_bounds <- Wald_bounds(estimates, var, alpha)
  lb <- sum(wald_bounds$LB)
  ub <- sum(wald_bounds$UB)
  
  # Create data frame with horizon and coefficients
  bounds_df <- data.frame(
    horizon = 1:p,
    coef = estimates,
    lower = lb / p, 
    upper = ub / p
  )

  # Calculate width
  width <- mean(bounds_df$upper - bounds_df$lower)
  
  # Return list with bounds data frame and metadata
  result <- list(
    cumulative_bounds = bounds_df
  )

  # Calculate pointwise bounds if requested
  if (include_pointwise) {
    result$pointwise_bounds <- calculate_pointwise_bounds(estimates, var, alpha)
  }
  
  # Calculate sup-t bounds if requested
  if (include_supt) {
    result$supt_bounds <- calculate_supt_bounds(estimates, var, alpha)
  }
  
  
  result$metadata <- list(
      alpha = alpha,
      width = width,
      individual_upper = wald_bounds$UB,
      individual_lower = wald_bounds$LB
    )
  
  class(result) <- c("cumulative_bounds", "plausible_bounds_result")
  return(result)
}