#' Calculate Plausible Bounds
#'
#' This function calculates both cumulative and restricted bounds for a vector of estimates.
#' It calls both calculate_cumulative_bounds and calculate_restricted_bounds functions.
#'
#' @param estimates A vector of point estimates
#' @param var The variance-covariance matrix of the estimates
#' @param alpha Significance level (default: 0.05)
#' @param include_pointwise Whether to include pointwise bounds (default: TRUE)
#' @param include_supt Whether to include sup-t bounds (default: TRUE)
#' @param parallel Whether to use parallel processing for restricted bounds calculation (default: FALSE)
#'
#' @return A list containing:
#'   \item{cumulative_bounds}{Results from calculate_cumulative_bounds}
#'   \item{restricted_bounds}{Results from calculate_restricted_bounds}
#'
#' @examples
#' # Example with constant estimates and IID errors (simple case)
#' data(estimates_constant_iid)
#' data(var_constant_iid)
#' result1 <- plausible_bounds(estimates_constant_iid, var_constant_iid)
#' print(result1)
#'
#' # Example with wiggly estimates and strong correlation (complex case)
#' data(estimates_wiggly_strong_corr)
#' data(var_wiggly_strong_corr)
#' result2 <- plausible_bounds(estimates_wiggly_strong_corr, var_wiggly_strong_corr)
#' print(result2)
#'
#' @export
plausible_bounds <- function(estimates, var, alpha = 0.05,
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
  
  # Calculate cumulative bounds
  cumul_bd <- calculate_cumulative_bounds(estimates, var, alpha, 
                                         include_pointwise, include_supt)
  
  # Calculate restricted bounds
  restr_bd <- calculate_restricted_bounds(estimates, var, alpha,
                                         include_pointwise = FALSE, include_supt = FALSE,
                                         parallel = parallel)
  
  # Return combined results
  result <- list(
    cumulative_bounds = cumul_bd$cumulative_bounds,
    restricted_bounds = restr_bd$restricted_bounds,
    cumulative_metadata = cumul_bd$metadata,
    restricted_metadata = restr_bd$metadata
  )

  if(include_pointwise) {
    result$pointwise_bounds <- cumul_bd$pointwise_bounds
  }

  if(include_supt) {
    result$supt_bounds <- cumul_bd$supt_bounds
  }
  
  class(result) <- "plausible_bounds"
  return(result)
}

#' Print method for plausible_bounds objects
#'
#' @param x A plausible_bounds object
#' @param ... Additional arguments passed to print
#'
#' @export
print.plausible_bounds <- function(x, ...) {
  cat("Plausible Bounds Analysis\n")
  cat("------------------------\n")
  
  cat("\nCumulative Bounds:\n")
  cat("  Width:", x$cumulative_metadata$width, "\n")
  cat("  Alpha:", x$cumulative_metadata$alpha, "\n")
  
  cat("\nRestricted Bounds:\n")
  cat("  Width:", x$restricted_metadata$width, "\n")
  cat("  Alpha:", x$restricted_metadata$alpha, "\n")
  cat("  Surrogate class:", x$restricted_metadata$surrogate_class, "\n")
  cat("  Degrees of freedom:", x$restricted_metadata$df, "\n")
  
  cat("\nUse create_plot() to visualize the results.\n")
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