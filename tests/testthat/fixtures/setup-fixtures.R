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
data(var_constant, envir = environment())
data(estimates_bighump, envir = environment())
data(var_bighump, envir = environment())

# Generate and save fixtures

# Simple plausible bounds
simple_plausible <- plausible_bounds(
  estimates_constant, var_constant
)
saveRDS(simple_plausible, file.path(fixtures_dir, "simple_plausible.rds"))

# Complex plausible bounds
complex_plausible <- plausible_bounds(
  estimates_bighump, var_bighump
)
saveRDS(complex_plausible, file.path(fixtures_dir, "complex_plausible.rds"))


# Preperiods -------

create_var <- function(rho, se, p) {
  corr_mat  <- stats::toeplitz(rho^(0:(p - 1)))
  sd_vector <- (100 + seq_len(p)) / 100
  se * diag(sd_vector) %*% corr_mat %*% diag(sd_vector)
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
Vpre_reject <- sigma2_pre * stats::toeplitz(0.95^(0:(npre - 1)))

pre_mean0 <- rnorm(npre, mean = 0, sd = sqrt(sigma2_pre))

pre_reject <- rnorm(npre, mean = 0.08, sd = sqrt(sigma2_pre))

saveRDS(pre_mean0,  file.path(fixtures_dir, "preperiods_mean0.rds"))
saveRDS(pre_reject, file.path(fixtures_dir, "preperiods_reject.rds"))

saveRDS(Vpre_iid,    file.path(fixtures_dir, "Vpre_iid.rds"))
saveRDS(Vpre_reject, file.path(fixtures_dir, "Vpre_reject.rds"))



