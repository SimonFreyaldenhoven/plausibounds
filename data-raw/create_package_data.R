## Code to prepare built-in datasets for the plausibounds package
## This script generates example data for demonstrating package functionality

# Helper functions for data generation ------------------------------------

#' Create Dynamic Path (Treatment Effect Path)
create_dynamic_path <- function(design_name, p = 12) {

  if (design_name == "constant") {
    # Smooth, flat (constant effect)
    delta <- rep(-0.4, p)

  } else if (design_name == "bighump") {
    # Big hump design: sinusoidal pattern in first 6 periods, then flat at zero
    delta <- -0.35 - 0.35 * sin(3/2 * (1:6) * pi / 6)
    # Pad with zeros for remaining periods
    if (p > 6) {
      delta <- c(delta, rep(0, p - 6))
    }

  } else {
    stop("Invalid design name. Choose from: 'constant', 'bighump'")
  }

  return(delta)
}

#' Create Variance-Covariance Matrix
create_var <- function(rho, se = 0.014, p = 12) {

  # Create Toeplitz correlation matrix with rho^|i-j| structure
  indices <- 0:(p-1)
  corr_mat <- outer(indices, indices, function(i, j) rho^abs(i - j))

  # Create heteroskedastic standard deviation vector
  sd_vector <- (100 + (1:p)) / 100

  # Construct variance-covariance matrix
  V <- se * diag(sd_vector) %*% corr_mat %*% diag(sd_vector)

  return(V)
}

#' Generate Simulated Data from DGP
simulate_dgp <- function(design_name, rho = 0.0, se = 0.014, p = 12) {

  # Generate true treatment path
  delta <- create_dynamic_path(design_name, p)

  # Generate variance-covariance matrix
  Vhat <- create_var(rho, se, p)

  # Generate noisy estimates
  dhat <- delta + MASS::mvrnorm(n = 1, mu = rep(0, p), Sigma = Vhat)

  return(list(
    delta = delta,
    dhat = as.vector(dhat),
    Vhat = Vhat,
    design = design_name,
    rho = rho,
    p = p
  ))
}

# Generate package datasets ------------------------------------------------

# Set seed for reproducibility
set.seed(916)

# Generate constant estimates with IID errors (rho = 0)
sim_constant <- simulate_dgp("constant", rho = 0, se = 0.014, p = 12)
estimates_constant <- sim_constant$dhat
var_iid <- sim_constant$Vhat

# Generate bighump estimates with moderate correlation (rho = 0.5)
# Big hump in first 6 periods, then flat at zero for remaining 30 periods = 36 total
sim_bighump <- simulate_dgp("bighump", rho = 0.5, se = 0.014, p = 36)
estimates_bighump <- sim_bighump$dhat
var_bighump <- sim_bighump$Vhat

# Save datasets to data/ directory -----------------------------------------
library(R.matlab)
smooth <- readMat("analysis/preperiod_data.mat")

var_smooth <- ex$Vall
estimates_smooth <- as.numeric(ex$dall)

usethis::use_data(estimates_smooth, overwrite = TRUE)
usethis::use_data(var_smooth, overwrite = TRUE)
usethis::use_data(estimates_constant, overwrite = TRUE)
usethis::use_data(var_iid, overwrite = TRUE)
usethis::use_data(estimates_bighump, overwrite = TRUE)
usethis::use_data(var_bighump, overwrite = TRUE)

