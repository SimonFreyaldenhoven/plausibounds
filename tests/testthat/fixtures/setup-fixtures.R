# Setup script to generate test fixtures for create_plot tests
# This script generates fixtures that will be used by test-create_plot.R
# It will overwrite existing fixtures when run

set.seed(20)

# Create fixtures directory if it doesn't exist
fixtures_dir <- "tests/testthat/fixtures"
if (!dir.exists(fixtures_dir)) {
  dir.create(fixtures_dir, recursive = TRUE)
}

# Load required data
data(estimates_constant, envir = environment())
data(var_iid, envir = environment())
data(estimates_wiggly, envir = environment())
data(var_corr, envir = environment())

# Generate and save fixtures

# Simple cumulative bounds
simple_cumulative <- calculate_cumulative_bounds(
  estimates_constant, var_iid
)
saveRDS(simple_cumulative, file.path(fixtures_dir, "simple_cumulative.rds"))

# Simple plausible bounds
simple_plausible <- plausible_bounds(
  estimates_constant, var_iid
)
saveRDS(simple_plausible, file.path(fixtures_dir, "simple_plausible.rds"))

# Complex cumulative bounds
complex_cumulative <- calculate_cumulative_bounds(
  estimates_wiggly, var_corr
)
saveRDS(complex_cumulative, file.path(fixtures_dir, "complex_cumulative.rds"))

# Complex restricted bounds (expensive computation)
complex_restricted <- calculate_restricted_bounds(
  estimates_wiggly, var_corr
)
saveRDS(complex_restricted, file.path(fixtures_dir, "complex_restricted.rds"))

# Complex plausible bounds
complex_plausible <- plausible_bounds(
  estimates_wiggly, var_corr
)
saveRDS(complex_plausible, file.path(fixtures_dir, "complex_plausible.rds"))


# Preperiods -------

create_var <- function(rho, se, p) {
  corr_mat  <- stats::toeplitz(rho^(0:(p - 1)))
  sd_vector <- (100 + seq_len(p)) / 100
  se * diag(sd_vector) %*% corr_mat %*% diag(sd_vector)
}

block_diag <- function(A, B) {
  out <- matrix(0, nrow(A) + nrow(B), ncol(A) + ncol(B))
  out[1:nrow(A), 1:ncol(A)] <- A
  out[(nrow(A)+1):(nrow(A)+nrow(B)), (ncol(A)+1):(ncol(A)+ncol(B))] <- B
  out
}

combine_pre_post <- function(pre, Vpre, post, Vpost) {
  stopifnot(length(pre)  == nrow(Vpre),  nrow(Vpre)  == ncol(Vpre))
  stopifnot(length(post) == nrow(Vpost), nrow(Vpost) == ncol(Vpost))
  list(d = c(pre, post), V = block_diag(Vpre, Vpost))
}

draw_mvn <- function(mu, Sigma) {
  U <- chol(Sigma)
  as.numeric(mu + t(U) %*% rnorm(length(mu)))
}

post <- as.numeric(estimates_constant)
p <- length(post)

se_values   <- c(0.014, 0.014 * 0.5^7, 0.014)
corr_values <- c(0.0, 0.0, 0.95)

npre <- 8
se_pre <- se_values[1]

Vpost <- create_var(corr_values[1], se_values[1], p)

Vpre_base  <- create_var(rho = 0, se = se_pre, p = npre)
sigma2_pre <- mean(diag(Vpre_base))

Vpre_iid  <- sigma2_pre * diag(npre)
Vpre_corr <- create_var(rho = corr_values[3], se = se_pre, p = npre)
Vpre_reject <- sigma2_pre * stats::toeplitz(0.95^(0:(npre - 1)))

pre_mean0 <- rnorm(npre, mean = 0, sd = sqrt(sigma2_pre))

pre_fixed <- c(
  0.395347730453505,
  0.363670223170499,
  0.328916662914659,
  0.300151381192354,
  0.288201567338562,
  0.276213177240659,
  0.261694387925979,
  0.241136813588603
)

pre_corr <- rep(0.3, npre) + draw_mvn(rep(0, npre), Vpre_corr)

pre_reject <- rnorm(npre, mean = 0.08, sd = sqrt(sigma2_pre))

combo_reject <- combine_pre_post(pre_reject, Vpre_reject, post, Vpost)

saveRDS(pre_mean0,  file.path(fixtures_dir, "preperiods_mean0.rds"))
saveRDS(pre_fixed,  file.path(fixtures_dir, "preperiods_fixed.rds"))
saveRDS(pre_corr,   file.path(fixtures_dir, "preperiods_corr.rds"))
saveRDS(pre_reject, file.path(fixtures_dir, "preperiods_reject.rds"))

saveRDS(Vpre_iid,    file.path(fixtures_dir, "Vpre_iid.rds"))
saveRDS(Vpre_corr,   file.path(fixtures_dir, "Vpre_corr.rds"))
saveRDS(Vpre_reject, file.path(fixtures_dir, "Vpre_reject.rds"))



