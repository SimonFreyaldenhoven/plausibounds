#' Calculate Plausible Bounds
#'
#' This function calculates both cumulative and restricted bounds for a vector of estimates.
#' It calls both calculate_cumulative_bounds and calculate_restricted_bounds functions.
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
#' @param parallel Whether to use parallel processing for restricted bounds calculation (default: FALSE)
#'
#' @return A list containing:
#'   \item{cumulative_bounds}{Results from calculate_cumulative_bounds}
#'   \item{restricted_bounds}{Results from calculate_restricted_bounds}
#'   \item{ate}{Average treatment effect with standard error}
#'   \item{Wpre}{Wald test for pre-trends (if preperiods > 0)}
#'   \item{Wpost}{Wald test for no treatment effect}
#'
#' @examples
#' # Example with constant estimates and IID errors (simple case)
#' data(estimates_constant_iid)
#' data(var_constant_iid)
#' result1 <- plausible_bounds(estimates_constant_iid[1:4], var_constant_iid[1:4, 1:4])
#' print(result1)
#'
#'
#' @export
plausible_bounds <- function(estimates, var, alpha = 0.05,
                            preperiods = 0,
                            include_pointwise = TRUE, include_supt = TRUE,
                            parallel = FALSE) {
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
  if (!is.numeric(preperiods) || length(preperiods) != 1 || preperiods < 0 || preperiods != floor(preperiods)) {
    stop("preperiods must be a non-negative integer")
  }
  if (preperiods >= length(estimates)) {
    stop("preperiods must be less than the length of estimates")
  }

  # Calculate cumulative bounds
  cumul_bd <- calculate_cumulative_bounds(estimates, var, alpha,
                                         preperiods = preperiods)

  # Calculate restricted bounds
  restr_bd <- calculate_restricted_bounds(estimates, var, alpha,
                                         preperiods = preperiods,
                                         parallel = parallel)

  # Return combined results
  result <- list(
    cumulative_bounds = cumul_bd$cumulative_bounds,
    restricted_bounds = restr_bd$restricted_bounds,
    ate = cumul_bd$ate,
    Wpost = restr_bd$Wpost,
    cumulative_metadata = cumul_bd$metadata,
    restricted_metadata = restr_bd$metadata
  )

  # Add Wpre if preperiods > 0
  if (!is.null(restr_bd$Wpre)) {
    result$Wpre <- restr_bd$Wpre
  }

  # Calculate pointwise bounds if requested (for ALL periods)
  if (include_pointwise) {
    if (preperiods > 0) {
      # Split estimates and var for pre and post periods
      estimates_pre <- estimates[1:preperiods]
      estimates_post <- estimates[(preperiods + 1):length(estimates)]
      var_pre <- var[1:preperiods, 1:preperiods, drop = FALSE]
      var_post <- var[(preperiods + 1):nrow(var), (preperiods + 1):ncol(var), drop = FALSE]

      pw_pre <- calculate_pointwise_bounds(estimates_pre, var_pre, alpha)
      pw_post <- calculate_pointwise_bounds(estimates_post, var_post, alpha)
      result$pointwise_bounds <- list(
        lower = c(pw_pre$lower, pw_post$lower),
        upper = c(pw_pre$upper, pw_post$upper)
      )
    } else {
      result$pointwise_bounds <- calculate_pointwise_bounds(estimates, var, alpha)
    }
  }

  # Calculate sup-t bounds if requested (for ALL periods, using full variance)
  if (include_supt) {
    if (preperiods > 0) {
      # Use the sup-t critical value computed from ALL periods (from restricted metadata)
      supt_critval <- restr_bd$metadata$supt_critval
      supt_lower_all <- estimates - supt_critval * sqrt(diag(var))
      supt_upper_all <- estimates + supt_critval * sqrt(diag(var))
      result$supt_bounds <- list(
        lower = supt_lower_all,
        upper = supt_upper_all
      )
    } else {
      result$supt_bounds <- calculate_supt_bounds(estimates, var, alpha)
    }
  }

  class(result) <- "plausible_bounds"
  return(result)
}

#' Summary method for plausible_bounds objects
#'
#' @param object A plausible_bounds object
#' @param ... Additional arguments passed to summary
#'
#' @export
summary.plausible_bounds <- function(object, ...) {

  cumul_df <- object$cumulative_bounds %>%
    dplyr::rename(cumul_lower = lower, cumul_upper = upper) 

  restr_df <- object$restricted_bounds %>%
    dplyr::rename(restr_lower = lower, restr_upper = upper) %>%
    dplyr::select(-coef)

  result <- dplyr::left_join(cumul_df, restr_df, by = "horizon")
  
  class(result) <- c("summary.plausible_bounds", class(result))
  return(result)
}

#' Print method for summary.plausible_bounds objects
#'
#' @param x A summary.plausible_bounds object
#' @param ... Additional arguments passed to print
#'
#' @export
print.summary.plausible_bounds <- function(x, ...) {
  cat("Summary of Plausible Bounds Results\n")
  cat("----------------------------------\n")
  
  cat("\nBounds data.frame:\n")
  class(x) <- setdiff(class(x), "summary.plausible_bounds")
  print(x)

  invisible(x)

}