#' plausibounds: Plausible Bounds for Treatment Path Estimates
#'
#' Implements the additions to dynamic effect plots suggested in Freyaldenhoven and Hansen (2026). 
#' Data-driven smoothing delivers a smooth estimated path with potentially improved point estimation properties and confidence regions covering a surrogate that can be substantially tighter than conventional pointwise or uniform bands.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{plausible_bounds}}: Calculate plausible bounds
#'   \item \code{\link{create_plot}}: Create plots of the bounds
#' }
#'
#' @section Example Datasets:
#' The package includes example datasets to demonstrate the functionality:
#' \itemize{
#'   \item \code{\link{estimates_constant}} and \code{\link{var_iid}}:
#'         A simple case with constant estimates and no correlation
#'   \item \code{\link{estimates_bighump}} and \code{\link{var_bighump}}:
#'         A case with sinusoidal estimates and moderate correlation
#'   \item \code{\link{estimates_smooth}} and \code{\link{var_smooth}}:
#'         A smooth case with effects that slowly level off and no correlation, from Figure 1 of Freyaldenhoven and Hansen (2026)
#' }
#'
#' 
#' @name plausibounds-package
#' @aliases plausibounds
NULL