#' plausibounds: Plausible Bounds for Treatment Path Estimates
#'
#' Calculate and visualize plausible bounds as proposed in Freyaldenhoven (2025)
#' for dynamic treatment path estimates. Given a treatment path, variance matrix, and significance level, plot various uncertainty measures.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{plausible_bounds}}: Calculate restricted bounds
#'   \item \code{\link{create_plot}}: Create plots of the bounds
#' }
#'
#' @section Example Datasets:
#' The package includes example datasets to demonstrate the functionality:
#' \itemize{
#'   \item \code{\link{estimates_constant}} and \code{\link{var_iid}}:
#'         A simple case with constant estimates and independent errors
#'   \item \code{\link{estimates_wiggly}} and \code{\link{var_corr}}:
#'         A case with wiggly estimates and strongly correlated errors
#'   \item \code{\link{estimates_pretrends}} and \code{\link{var_pretrends}}:
#'         A case with 6 preperiods, 12 post
#' }
#'
#' 
#' @name plausibounds-package
#' @aliases plausibounds
NULL