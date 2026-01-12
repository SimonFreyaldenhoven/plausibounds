#' Calculate Plausible Bounds
#'
#' This function calculates the plausible bounds for a vector of estimates along with
#' average treatment effect, Wald tests, and optional pointwise/sup-t bounds.
#' Supports pre-treatment periods for event study designs.
#'
#' @param estimates A numeric vector or single-row/single-column matrix of point estimates.
#'   If preperiods > 0, the first preperiods elements are pre-treatment estimates,
#'   followed by post-treatment estimates.
#' @param var The variance-covariance matrix of the estimates
#' @param alpha Significance level (default: 0.05)
#' @param preperiods Number of pre-treatment periods (default: 0). Period 0 is assumed
#'   to be normalized and not included in estimates.
#' @param include_pointwise Whether to include pointwise bounds (default: TRUE)
#' @param include_supt Whether to include sup-t bounds (default: TRUE)
#' @param parallel Whether to use parallel processing for restricted bounds calculation (default: FALSE)
#' @param n_cores Number of cores to use for parallel processing (default: NULL, which uses
#'   detectCores() - 1). Only used when parallel = TRUE.
#'
#' @return A list containing:
#'   \item{alpha}{Significance level}
#'   \item{preperiods}{Number of pre-treatment periods}
#'   \item{wald_test}{List with post (and pre if preperiods > 0) Wald test results}
#'   \item{restricted_bounds}{Data frame with horizon, coef, surrogate, lower, upper}
#'   \item{restricted_bounds_metadata}{List with supt_critval, supt_b,
#'     degrees_of_freedom, K, lambda1, lambda2, surrogate_class, best_fit_model}
#'   \item{avg_treatment_effect}{List with estimate, se, lower, upper}
#'   \item{pointwise_bounds}{List with lower and upper vectors (if include_pointwise = TRUE)}
#'   \item{supt_bounds}{List with lower and upper vectors (if include_supt = TRUE)}
#'
#' @examples
#' # Example with constant estimates and IID errors (simple case)
#' data(estimates_constant)
#' data(var_iid)
#' result1 <- plausible_bounds(estimates_constant[1:4], var_iid[1:4, 1:4])
#' print(result1)
#'
#'
#' @export
plausible_bounds <- function(estimates, var, alpha = 0.05,
                            preperiods = 0,
                            include_pointwise = TRUE, include_supt = TRUE,
                            parallel = FALSE, n_cores = NULL) {
  # Coerce estimates to plain numeric vector
  if (is.matrix(estimates)) {
    if (ncol(estimates) == 1) {
      estimates <- as.vector(estimates)
    } else if (nrow(estimates) == 1) {
      estimates <- as.vector(t(estimates))
    } else {
      stop("If estimates is a matrix, it must have exactly one row or one column")
    }
  } else if (is.numeric(estimates) && length(estimates) > 0) {
    # Accept any numeric object with length (e.g., vectors, numeric with attributes)
    estimates <- as.vector(estimates)
  } else {
    stop("estimates must be a numeric vector or single-row/single-column matrix")
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

  # Calculate cumulative bounds (for ATE)
  cumul_bd <- calculate_cumulative_bounds(estimates, var, alpha,
                                         preperiods = preperiods)

  # Calculate restricted bounds
  restr_bd <- calculate_restricted_bounds(estimates, var, alpha,
                                         preperiods = preperiods,
                                         parallel = parallel,
                                         n_cores = n_cores)

  # Build wald_test
  wald_test <- list(
    post = list(
      statistic = unname(restr_bd$Wpost["statistic"]),
      p_value = unname(restr_bd$Wpost["pvalue"])
    )
  )
  if (!is.null(restr_bd$Wpre)) {
    wald_test$pre <- list(
      statistic = unname(restr_bd$Wpre["statistic"]),
      p_value = unname(restr_bd$Wpre["pvalue"])
    )
  }

  # Build restricted_bounds_metadata
  restricted_bounds_metadata <- list(
    supt_critval = restr_bd$metadata$supt_critval,
    supt_b = restr_bd$metadata$suptb,
    degrees_of_freedom = restr_bd$metadata$df,
    K = restr_bd$metadata$K,
    lambda1 = restr_bd$metadata$lambda1,
    lambda2 = restr_bd$metadata$lambda2,
    surrogate_class = restr_bd$metadata$surrogate_class,
    best_fit_model = restr_bd$metadata$best_fit_model
  )

  # Build avg_treatment_effect
  avg_treatment_effect <- list(
    estimate = unname(cumul_bd$ate["estimate"]),
    se = unname(cumul_bd$ate["se"]),
    lower = cumul_bd$metadata$lb,
    upper = cumul_bd$metadata$ub
  )

  # Build result list

  result <- list(
    alpha = alpha,
    preperiods = preperiods,
    wald_test = wald_test,
    restricted_bounds = restr_bd$restricted_bounds,
    restricted_bounds_metadata = restricted_bounds_metadata,
    avg_treatment_effect = avg_treatment_effect
  )

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
  # Build summary from restricted_bounds (which has coef, surrogate, lower, upper)
  result <- object$restricted_bounds %>%
    dplyr::rename(restr_lower = lower, restr_upper = upper)

  class(result) <- c("summary.plausible_bounds", class(result))
  return(result)
}

#' Print method for plausible_bounds objects
#'
#' @param x A plausible_bounds object
#' @param ... Additional arguments passed to print
#'
#' @export
print.plausible_bounds <- function(x, ...) {
  cat("Plausible Bounds Results\n")
  cat("========================\n\n")

  # Wald tests
  cat("Wald Tests:\n")
  cat(sprintf("  Post-treatment (H0: no effect): stat = %.3f, p = %.4f\n",
              x$wald_test$post$statistic, x$wald_test$post$p_value))
  if (!is.null(x$wald_test$pre)) {
    cat(sprintf("  Pre-treatment (H0: no pre-trends): stat = %.3f, p = %.4f\n",
                x$wald_test$pre$statistic, x$wald_test$pre$p_value))
  }

  # ATE
  cat("\nAverage Treatment Effect:\n")
  cat(sprintf("  Estimate: %.4f (SE: %.4f)\n",
              x$avg_treatment_effect$estimate, x$avg_treatment_effect$se))
  cat(sprintf("  %.0f%% CI: [%.4f, %.4f]\n",
              (1 - x$alpha) * 100,
              x$avg_treatment_effect$lower, x$avg_treatment_effect$upper))

  # Restricted bounds info
  cat("\nRestricted Bounds:\n")
  cat(sprintf("  Surrogate class: %s\n", x$restricted_bounds_metadata$surrogate_class))
  cat(sprintf("  Degrees of freedom: %.2f\n", x$restricted_bounds_metadata$degrees_of_freedom))

  cat("\nBounds by horizon:\n")
  print(x$restricted_bounds, row.names = FALSE)

  invisible(x)
}

#' Print method for summary.plausible_bounds objects
#'
#' @param x A summary.plausible_bounds object
#' @param ... Additional arguments passed to print
#'
#' @export
print.summary.plausible_bounds <- function(x, ...) {
  cat("Summary of Plausible Bounds Results\n")
  cat("-----------------------------------\n\n")

  class(x) <- setdiff(class(x), "summary.plausible_bounds")
  print(x, row.names = FALSE)

  invisible(x)
}