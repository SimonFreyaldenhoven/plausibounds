
# Test basic functionality
test_that("plausible_bounds functions work with simple test data", {
  # Create test data
  set.seed(123)
  p <- 4
  estimates <- stats::rnorm(p)
  var <- diag(p) * 0.1

  
  # Test plausible_bounds
  pb <- plausible_bounds(estimates, var)
  expect_s3_class(pb, "plausible_bounds")
  
  # Test create_plot (just check that it runs without error)
  expect_error(create_plot(pb), NA)
  expect_error(create_plot(pb, show_supt = TRUE, show_pointwise = FALSE), NA)
  expect_error(create_plot(pb, show_supt = FALSE, show_pointwise = TRUE), NA)
})

# Test comparison with original function using pre-computed fixtures
test_that("new functions match original full_eventplot_l2tf results", {
  # Skip this test if the fixture files don't exist
  skip_if_not(file.exists(test_path("fixtures", "dhat.csv")) &&
              file.exists(test_path("fixtures", "vhat.csv")) &&
              file.exists(test_path("fixtures", "full_eventplot_l2tf_results.rds")))
  
  skip_on_cran()
  
  # Load test data
  delta <- as.matrix(read.csv(test_path("fixtures", "dhat.csv"), header = FALSE))
  vhat <- as.matrix(read.csv(test_path("fixtures", "vhat.csv"), header = FALSE))
  
  # Load pre-computed original results from fixture
  original_results <- readRDS(test_path("fixtures", "full_eventplot_l2tf_results.rds"))
  
  # Run new functions
  pb <- plausible_bounds(as.vector(delta), vhat)
  
  # Compare results
  # Average treatment effect bounds
  expect_equal(
    pb$avg_treatment_effect$bounds$lb,
    original_results$lb/length(delta),
    tolerance = 1e-2
  )
  expect_equal(
    pb$avg_treatment_effect$bounds$ub,
    original_results$ub/length(delta),
    tolerance = 1e-2
  )
  
  # Restricted bounds
  expect_equal(
    pb$restricted_bounds$surrogate,
    as.vector(original_results$surrogate),
    tolerance = 1e-2
  )
  expect_equal(
    pb$restricted_bounds$lower,
    as.vector(original_results$restricted_LB),
    tolerance = 1e-2
  )
  expect_equal(
    pb$restricted_bounds$upper,
    as.vector(original_results$restricted_UB),
    tolerance = 1e-2
  )
})

# Test parallel vs non-parallel results
test_that("parallel and non-parallel results are identical", {
  skip_on_cran()
  l <- 6

  # Test with wiggly estimates and strong correlation
  data(estimates_wiggly)
  data(var_corr)

  # Run with parallel = TRUE
  set.seed(42)  # Set seed for reproducibility
  result_parallel <- plausible_bounds(
    estimates_wiggly[1:l],
    var_corr[1:l,1:l],
    parallel = TRUE
  )

  # Run with parallel = FALSE
  set.seed(42)  # Reset seed to ensure same random numbers
  result_sequential <- plausible_bounds(
    estimates_wiggly[1:l],
    var_corr[1:l, 1:l],
    parallel = FALSE
  )
  
  # Compare results - restricted bounds
  expect_equal(
    result_parallel$restricted_bounds,
    result_sequential$restricted_bounds,
    tolerance = 1e-10
  )
  
  # Compare metadata
  expect_equal(
    result_parallel$restricted_bounds_metadata,
    result_sequential$restricted_bounds_metadata,
    tolerance = 1e-10
  )

  expect_equal(
    result_parallel$restricted_bounds_metadata$width,
    result_sequential$restricted_bounds_metadata$width,
    tolerance = 1e-10
  )
  
  # Test with constant estimates and IID errors
  data(estimates_constant)
  data(var_iid)
  
  # Run with parallel = TRUE
  set.seed(42)  # Set seed for reproducibility
  result_parallel_const <- plausible_bounds(
    estimates_constant[1:l],
    var_iid[1:l,1:l],
    parallel = TRUE
  )
  
  # Run with parallel = FALSE
  set.seed(42)  # Reset seed to ensure same random numbers
  result_sequential_const <- plausible_bounds(
    estimates_constant[1:l],
    var_iid[1:l,1:l],
    parallel = FALSE
  )
  
  # Compare results - restricted bounds
  expect_equal(
    result_parallel_const$restricted_bounds,
    result_sequential_const$restricted_bounds,
    tolerance = 1e-10
  )
  
  # Compare metadata
  expect_equal(
    result_parallel_const$restricted_bounds_metadata,
    result_sequential_const$restricted_bounds_metadata,
    tolerance = 1e-10
  )

  expect_equal(
    result_parallel_const$restricted_bounds_metadata$width,
    result_sequential_const$restricted_bounds_metadata$width,
    tolerance = 1e-10
  )
})

# Test n_cores parameter validation and edge cases
test_that("plausible_bounds validates n_cores parameter", {
  skip_if_not_installed("parallel")
  skip_on_cran()

  set.seed(123)
  p <- 6
  estimates <- stats::rnorm(p)
  var <- diag(p) * 0.1

  # Invalid: non-numeric
  expect_error(
    plausible_bounds(estimates, var, parallel = TRUE, n_cores = "two"),
    "n_cores must be a positive integer"
  )

  # Invalid: negative
  expect_error(
    plausible_bounds(estimates, var, parallel = TRUE, n_cores = -1),
    "n_cores must be a positive integer"
  )

  # Invalid: zero
  expect_error(
    plausible_bounds(estimates, var, parallel = TRUE, n_cores = 0),
    "n_cores must be a positive integer"
  )

  # Invalid: non-integer
  expect_error(
    plausible_bounds(estimates, var, parallel = TRUE, n_cores = 3.5),
    "n_cores must be a positive integer"
  )

  # Invalid: vector
  expect_error(
    plausible_bounds(estimates, var, parallel = TRUE, n_cores = c(2, 4)),
    "n_cores must be a positive integer"
  )
})

test_that("plausible_bounds works with specific n_cores values", {
  skip_if_not_installed("parallel")
  skip_on_cran()

  set.seed(456)
  p <- 6
  estimates <- stats::rnorm(p)
  var <- diag(p) * 0.1

  # Test with n_cores = 1
  set.seed(42)
  result_1core <- plausible_bounds(estimates, var, parallel = TRUE, n_cores = 1)

  # Test with n_cores = 2
  set.seed(42)
  result_2cores <- plausible_bounds(estimates, var, parallel = TRUE, n_cores = 2)

  # Test with parallel = FALSE for comparison
  set.seed(42)
  result_seq <- plausible_bounds(estimates, var, parallel = FALSE)

  # All should produce valid results
  expect_s3_class(result_1core, "plausible_bounds")
  expect_s3_class(result_2cores, "plausible_bounds")

  # Results should be very similar across different n_cores
  expect_equal(result_1core$restricted_bounds$horizon,
               result_2cores$restricted_bounds$horizon)
  expect_equal(result_1core$restricted_bounds$coef,
               result_2cores$restricted_bounds$coef)
  expect_equal(result_1core$restricted_bounds$surrogate,
               result_2cores$restricted_bounds$surrogate,
               tolerance = 0.01)
})
