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
#' result <- plausible_bounds(estimates_constant[1:7], var_iid[1:7, 1:7])
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
#' result <- plausible_bounds(estimates_constant[1:7], var_iid[1:7, 1:7])
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
#' result <- plausible_bounds(estimates_wiggly[1:7], var_corr[1:7, 1:7])
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
#' result <- plausible_bounds(estimates_wiggly[1:7], var_corr[1:7, 1:7])
#' create_plot(result)
"var_corr"

#' Pretrends Estimates with Significant Pre-Treatment Trends
#'
#' A dataset containing estimates with significant pretrends in the first 6
#' pre-treatment periods followed by 12 post-treatment periods. This dataset
#' exhibits a linear trend before treatment that violates the parallel trends
#' assumption, making it useful for testing pretrend detection methods.
#' Statistical tests would reject the null hypothesis of no pretrends.
#'
#' @format A numeric vector with 18 elements (6 pre-treatment, 12 post-treatment)
#' @source Generated from simulation with linear pretrend design and moderate correlation (rho = 0.4)
#' @examples
#' data(estimates_pretrends)
#' data(var_pretrends)
#' # Use first 7 post-treatment estimates (indices 1:13 = 6 pre + 7 post)
#' result <- plausible_bounds(estimates_pretrends[1:13], var_pretrends[1:13, 1:13], preperiods = 6)
#' create_plot(result)
"estimates_pretrends"

#' Variance Matrix for Pretrends Estimates
#'
#' A variance-covariance matrix for the pretrends estimates with moderate
#' correlation structure (rho = 0.4). This matrix has non-zero off-diagonal
#' elements representing moderate correlation across time periods.
#'
#' @format A 18 x 18 matrix
#' @source Generated from simulation with linear pretrend design and moderate correlation (rho = 0.4)
#' @examples
#' data(estimates_pretrends)
#' data(var_pretrends)
#' # Use first 7 post-treatment estimates (indices 1:13 = 6 pre + 7 post)
#' result <- plausible_bounds(estimates_pretrends[1:13], var_pretrends[1:13, 1:13], preperiods = 6)
#' create_plot(result)
"var_pretrends"