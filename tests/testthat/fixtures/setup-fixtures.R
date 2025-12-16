# Setup script to generate test fixtures for create_plot tests
# This script generates fixtures that will be used by test-create_plot.R
# It will overwrite existing fixtures when run

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