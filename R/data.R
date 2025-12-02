#' Smooth Eventually Flat Estimates with IID Errors
#'
#' A dataset containing estimates from a quadratic (hump-shaped, eventually flat) design
#' with independent and identically distributed errors. The true effect follows a quadratic
#' path that peaks and then becomes flat, representing a treatment effect that stabilizes
#' after an initial dynamic period. The errors are independent with variance increasing
#' slightly with horizon.
#'
#' @format A numeric vector with 12 elements
#' @source Generated from simulation with quadratic design and IID errors (rho = 0)
#' @examples
#' data(estimates_smooth_iid)
#' data(var_smooth_iid)
#' result <- calculate_cumulative_bounds(estimates_smooth_iid, var_smooth_iid)
#' create_plot(result)
"estimates_smooth_iid"

#' Variance Matrix for Smooth Eventually Flat Estimates with IID Errors
#'
#' A variance-covariance matrix for the smooth eventually flat estimates with independent
#' and identically distributed errors. This is a diagonal matrix where the variance
#' increases slightly with horizon, representing heteroskedasticity but no correlation
#' across time periods.
#'
#' @format A 12 x 12 diagonal matrix
#' @source Generated from simulation with quadratic design and IID errors (rho = 0)
#' @examples
#' data(estimates_smooth_iid)
#' data(var_smooth_iid)
#' result <- calculate_cumulative_bounds(estimates_smooth_iid, var_smooth_iid)
#' create_plot(result)
"var_smooth_iid"

#' Wiggly Estimates with Moderate Correlation in Errors
#'
#' A dataset containing estimates from a wiggly (oscillating) design with moderate
#' correlation in errors. The true effect follows an oscillating path with noise that
#' eventually becomes flat, representing a treatment effect with complex dynamics before
#' stabilization. The errors have moderate correlation structure (rho = 0.8), with
#' correlation decreasing as the time between periods increases.
#'
#' @format A numeric vector with 12 elements
#' @source Generated from simulation with wiggly design and moderate correlation (rho = 0.8)
#' @examples
#' data(estimates_wiggly_corr)
#' data(var_wiggly_corr)
#' result <- calculate_cumulative_bounds(estimates_wiggly_corr, var_wiggly_corr)
#' create_plot(result)
"estimates_wiggly_corr"

#' Variance Matrix for Wiggly Estimates with Moderate Correlation in Errors
#'
#' A variance-covariance matrix for the wiggly estimates with moderate correlation in errors.
#' This is a full matrix with non-zero off-diagonal elements following a Toeplitz structure
#' where correlation decays as rho^|i-j| (with rho = 0.8). This represents a realistic
#' scenario where errors across time periods are correlated, with nearby periods more
#' strongly correlated than distant ones.
#'
#' @format A 12 x 12 matrix
#' @source Generated from simulation with wiggly design and moderate correlation (rho = 0.8)
#' @examples
#' data(estimates_wiggly_corr)
#' data(var_wiggly_corr)
#' result <- calculate_cumulative_bounds(estimates_wiggly_corr, var_wiggly_corr)
#' create_plot(result)
"var_wiggly_corr"