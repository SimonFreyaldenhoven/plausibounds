## Code to prepare built-in datasets for the plausibounds package
## This script generates example data for demonstrating package functionality

# Helper functions for data generation ------------------------------------

#' Create Dynamic Path (Treatment Effect Path)
create_dynamic_path <- function(design_name, p = 12) {

  if (design_name == "quadratic") {
    # Hump-shaped, eventually flat
    peak_period <- ceiling(p * 16/36)
    delta <- -0.3 + (peak_period - (1:p))^2 / 150
    if (peak_period <= p) {
      delta[peak_period:p] <- -0.3
    }

  } else if (design_name == "wiggly") {
    # Wiggly (oscillating) path with noise, eventually flat
    set.seed(12082023)
    flat_period <- ceiling(p * 20/36)
    delta <- -0.4 * sin((0:(p-1)) * pi / (p-1))
    if (flat_period <= p) {
      delta[flat_period:p] <- -0.4
    }
    noise <- 0.1 * rnorm(p)
    delta <- delta + noise - mean(noise)

  } else if (design_name == "constant") {
    # Smooth, flat (constant effect)
    delta <- rep(-0.4, p)

  } else {
    stop("Invalid design name. Choose from: 'quadratic', 'wiggly', 'constant'")
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

# Generate smooth, eventually flat estimates with IID errors (rho = 0)
sim_eventually_flat <- simulate_dgp("quadratic", rho = 0, se = 0.014, p = 12)
estimates_smooth_iid <- sim_eventually_flat$dhat
var_smooth_iid <- sim_eventually_flat$Vhat

# Generate wiggly estimates with moderate correlation (rho = 0.8)
sim_wiggly_corr <- simulate_dgp("wiggly", rho = 0.8, se = 0.014, p = 12)
estimates_wiggly_corr <- sim_wiggly_corr$dhat
var_wiggly_corr <- sim_wiggly_corr$Vhat

# Save datasets to data/ directory -----------------------------------------

usethis::use_data(estimates_smooth_iid, overwrite = TRUE)
usethis::use_data(var_smooth_iid, overwrite = TRUE)
usethis::use_data(estimates_wiggly_corr, overwrite = TRUE)
usethis::use_data(var_wiggly_corr, overwrite = TRUE)

