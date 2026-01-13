#' Constant Estimates with IID Errors
#'
#' A dataset containing estimates from a constant design, independent across horizons.
#' The true effect is approximately constant, representing a simple treatment effect scenario.
#'
#' @format A numeric vector with 12 elements
#' @source Generated from simulation with constant design and IID errors
#' @examples
#' data(estimates_constant)
#' data(var_constant)
#' result <- plausible_bounds(estimates_constant[1:4], var_constant[1:4, 1:4])
#' create_plot(result)
"estimates_constant"

#' Variance Matrix for Constant Estimates with IID Errors
#'
#' A variance matrix for the constant, independent, and identically distributed estimates.
#' This is a diagonal matrix representing no correlation across time periods.
#'
#' @format A 12 x 12 matrix
#' @source Generated from simulation with constant design and IID errors
#' @examples
#' data(estimates_constant)
#' data(var_constant)
#' result <- plausible_bounds(estimates_constant[1:4], var_constant[1:4, 1:4])
#' create_plot(result)
"var_constant"

#' Smooth Estimates from Freyaldenhoven and Hansen (2026)
#'
#' A dataset containing smooth treatment effect estimates that dip down and then level off.
#' The first 8 observations comprise the preperiods, the next 36 are post-period.
#'
#' @format A numeric vector with 44 elements (8 preperiods, 36 postperiods)
#' @source Point estimates from in Figure 1 of Freyaldenhoven and Hansen (2026)
#' @examples
#' data(estimates_smooth)
#' data(var_smooth)
#' result <- plausible_bounds(estimates_smooth[9:13], var_smooth[9:13, 9:13])
#' create_plot(result)
"estimates_smooth"

#' Variance Matrix for Smooth Estimates
#'
#' A variance matrix for the smooth estimates. This matrix captures
#' the correlation structure of the estimates.
#' The first 8 rows and 8 columns comprise the variance matrix for the 8 preperiods.
#'
#' @format A 44 x 44 matrix 
#' @source Variance Matrix of Figure 1 of Freyaldenhoven and Hansen (2026)
#' @examples
#' data(estimates_smooth)
#' data(var_smooth)
#' result <- plausible_bounds(estimates_smooth[9:13], var_smooth[9:13, 9:13])
#' create_plot(result)
"var_smooth"

#' Sinusoidal Estimates with Moderate Correlation
#'
#' A dataset containing estimates with a curved sinusoidal pattern in the first
#' 6 periods that then converges to zero for the remaining 30 periods. The effect
#' is a smooth curved trajectory with a large dip before the treatment effect quickly returns to 0.
#' Generated with moderate correlation (rho = 0.5).
#'
#' @format A numeric vector with 36 elements
#' @source Generated from simulation with sinusoidal design and moderate correlation (rho = 0.5)
#' @examples
#' data(estimates_bighump)
#' data(var_bighump)
#' result <- plausible_bounds(estimates_bighump[1:4], var_bighump[1:4, 1:4])
#' create_plot(result)
"estimates_bighump"

#' Variance Matrix for Sinusoidal Estimates
#'
#' A variance matrix for the sinusoidal estimates with moderate
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
