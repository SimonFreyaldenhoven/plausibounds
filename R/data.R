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
#' result <- plausible_bounds(estimates_constant[1:4], var_iid[1:4, 1:4])
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
#' result <- plausible_bounds(estimates_constant[1:4], var_iid[1:4, 1:4])
#' create_plot(result)
"var_iid"

#' Smooth Estimates from Empirical Data
#'
#' A dataset containing smooth treatment effect estimates from empirical data.
#' This dataset represents real-world treatment effects with gradual changes
#' over time. The first 8 observations comprise the preperiods.
#'
#' @format A numeric vector with estimates over time (starting after preperiod 8)
#' @source Point estimates in Figure 1 of Freyaldenhoven and Hansen (2026)
#' @examples
#' data(estimates_smooth)
#' data(var_smooth)
#' result <- plausible_bounds(estimates_smooth[9:13], var_smooth[9:13, 9:13])
#' create_plot(result)
"estimates_smooth"

#' Variance Matrix for Smooth Estimates
#'
#' A variance-covariance matrix for the smooth estimates. This matrix captures
#' the correlation structure of the estimation errors from the empirical data.
#' The first 8 rows and 8 columns comprise the variance matrix for the 8 preperiods.
#'
#' @format A square matrix matching the length of estimates_smooth
#' @source Variance Matrix in Figure 1 of Freyaldenhoven and Hansen (2026)
#' @examples
#' data(estimates_smooth)
#' data(var_smooth)
#' result <- plausible_bounds(estimates_smooth[9:13], var_smooth[9:13, 9:13])
#' create_plot(result)
"var_smooth"

#' Sinusoidal Estimates with Moderate Correlation
#'
#' A dataset containing estimates with a curved sinusoidal pattern in the first
#' 6 periods that then converges to zero for the remaining 30 periods. The
#' initial effect follows the formula -0.35 - 0.35*sin(3/2*(1:6)*pi/6),
#' creating a smooth curved trajectory before the treatment effect dissipates
#' completely. Generated with moderate correlation (rho = 0.5).
#'
#' @format A numeric vector with 36 elements (6 sinusoidal periods + 30 zero periods)
#' @source Generated from simulation with sinusoidal design and moderate correlation (rho = 0.5)
#' @examples
#' data(estimates_bighump)
#' data(var_bighump)
#' result <- plausible_bounds(estimates_bighump[1:4], var_bighump[1:4, 1:4])
#' create_plot(result)
"estimates_bighump"

#' Variance Matrix for Sinusoidal Estimates
#'
#' A variance-covariance matrix for the sinusoidal estimates with moderate
#' correlation structure (rho = 0.5). This matrix has non-zero off-diagonal
#' elements representing moderate correlation across time periods.
#'
#' @format A 36 x 36 matrix
#' @source Generated from simulation with sinusoidal design and moderate correlation (rho = 0.5)
#' @examples
#' data(estimates_bighump)
#' data(var_bighump)
#' result <- plausible_bounds(estimates_bighump[1:4], var_bighump[1:4, 1:4])
#' create_plot(result)
"var_bighump"
