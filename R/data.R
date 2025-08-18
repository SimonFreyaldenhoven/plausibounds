#' Constant Estimates with IID Errors
#'
#' A dataset containing constant estimates with independent and identically distributed errors.
#' This represents a simple case where the true effect is constant across all horizons,
#' and the errors are independent with variance increasing slightly with horizon.
#'
#' @format A numeric vector with 12 elements
#' @source Generated from simulation with constant mean and IID errors
#' @examples
#' data(estimates_constant_iid)
#' data(var_constant_iid)
#' result <- plausible_bounds(estimates_constant_iid, var_constant_iid)
#' create_plot(result)
"estimates_constant_iid"

#' Variance Matrix for Constant Estimates with IID Errors
#'
#' A variance-covariance matrix for the constant estimates with independent and identically 
#' distributed errors. This is a diagonal matrix where the variance increases slightly with horizon.
#'
#' @format A 12 x 12 diagonal matrix
#' @source Generated from simulation with constant mean and IID errors
#' @examples
#' data(estimates_constant_iid)
#' data(var_constant_iid)
#' result <- plausible_bounds(estimates_constant_iid, var_constant_iid)
#' create_plot(result)
"var_constant_iid"

#' Wiggly Estimates with Strong Correlation in Errors
#'
#' A dataset containing wiggly estimates with strong correlation in errors.
#' This represents a complex case where the true effect varies across horizons,
#' and the errors have strong correlation structure.
#'
#' @format A numeric vector with 12 elements
#' @source Generated from simulation with wiggly mean and strongly correlated errors
#' @examples
#' data(estimates_wiggly_strong_corr)
#' data(var_wiggly_strong_corr)
#' result <- plausible_bounds(estimates_wiggly_strong_corr, var_wiggly_strong_corr)
#' create_plot(result)
"estimates_wiggly_strong_corr"

#' Variance Matrix for Wiggly Estimates with Strong Correlation in Errors
#'
#' A variance-covariance matrix for the wiggly estimates with strong correlation in errors.
#' This is a full matrix with non-zero off-diagonal elements, indicating correlation between
#' the errors at different horizons.
#'
#' @format A 12 x 12 matrix
#' @source Generated from simulation with wiggly mean and strongly correlated errors
#' @examples
#' data(estimates_wiggly_strong_corr)
#' data(var_wiggly_strong_corr)
#' result <- plausible_bounds(estimates_wiggly_strong_corr, var_wiggly_strong_corr)
#' create_plot(result)
"var_wiggly_strong_corr"