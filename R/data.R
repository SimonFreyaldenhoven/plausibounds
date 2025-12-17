#' Constant Estimates with IID Errors
#'
#' A dataset containing estimates from a constant design with independent and
#' identically distributed errors. The true effect is approximately constant,
#' representing a simple treatment effect scenario.
#'
#' @format A numeric vector with 12 elements
#' @source Generated from simulation with constant design and IID errors
#' @examples
#' data(estimates_constant)
#' data(var_iid)
#' result <- calculate_cumulative_bounds(estimates_constant, var_iid)
#' create_plot(result)
"estimates_constant"

#' Variance Matrix for Constant Estimates with IID Errors
#'
#' A variance-covariance matrix for the constant estimates with independent
#' and identically distributed errors. This is a diagonal matrix representing
#' no correlation across time periods.
#'
#' @format A 12 x 12 diagonal matrix
#' @source Generated from simulation with constant design and IID errors
#' @examples
#' data(estimates_constant)
#' data(var_iid)
#' result <- calculate_cumulative_bounds(estimates_constant, var_iid)
#' create_plot(result)
"var_iid"

#' Wiggly Estimates with Correlated Errors
#'
#' A dataset containing estimates from a wiggly (oscillating) design with
#' correlated errors. The true effect follows an oscillating path,
#' representing a treatment effect with complex dynamics. The errors have
#' strong correlation structure.
#'
#' @format A numeric vector with 12 elements
#' @source Generated from simulation with wiggly design and correlated errors
#' @examples
#' data(estimates_wiggly)
#' data(var_corr)
#' result <- calculate_cumulative_bounds(estimates_wiggly, var_corr)
#' create_plot(result)
"estimates_wiggly"

#' Variance Matrix for Wiggly Estimates with Correlated Errors
#'
#' A variance-covariance matrix for the wiggly estimates with correlated errors.
#' This is a full matrix with non-zero off-diagonal elements representing
#' correlation across time periods.
#'
#' @format A 12 x 12 matrix
#' @source Generated from simulation with wiggly design and correlated errors
#' @examples
#' data(estimates_wiggly)
#' data(var_corr)
#' result <- calculate_cumulative_bounds(estimates_wiggly, var_corr)
#' create_plot(result)
"var_corr"