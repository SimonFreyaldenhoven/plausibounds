# Internal helper function for Wald bounds
Wald_bounds <- function(dhat, Vhat, alpha, df = 1) {
  h <- length(dhat)
  critval <- qchisq(1 - alpha, df)
  
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
#' @export
calculate_cumulative_bounds <- function(estimates, var, alpha = 0.05, 
                                       include_pointwise = TRUE, include_supt = TRUE) {
  # Check inputs
  if (!is.numeric(estimates) || !is.vector(estimates)) {
    stop("estimates must be a numeric vector")
  }
  if (!is.matrix(var) || nrow(var) != length(estimates) || ncol(var) != length(estimates)) {
    stop("var must be a square matrix with dimensions matching the length of estimates")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
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
    lower = lb, 
    upper = ub
  )
  
  # Calculate pointwise bounds if requested
  if (include_pointwise) {
    pointwise_critical <- qnorm(1 - alpha/2)
    bounds_df$pointwise_lower <- estimates - pointwise_critical * sqrt(diag(var))
    bounds_df$pointwise_upper <- estimates + pointwise_critical * sqrt(diag(var))
  }
  
  # Calculate sup-t bounds if requested
  if (include_supt) {
    supt_bands <- bands_plugin(estimates, var, p, nsim = 2000, level = 1 - alpha)
    bounds_df$supt_lower <- supt_bands$LB
    bounds_df$supt_upper <- supt_bands$UB
    supt_critval <- supt_bands$sup_t
  } else {
    supt_critval <- NA
  }
  
  
  # Calculate width
  width <- mean(bounds_df$upper - bounds_df$lower)
  
  # Return list with bounds data frame and metadata
  result <- list(
    bounds = bounds_df,
    type = "cumulative",
    metadata = list(
      alpha = alpha,
      width = width,
      individual_upper = wald_bounds$UB,
      individual_lower = wald_bounds$LB
    )
  )
  
  class(result) <- c("cumulative_bounds", "plausible_bounds_result")
  return(result)
}