
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
  expect_error(create_plot(pb, show_cumulative = TRUE, show_restricted = FALSE), NA)
  expect_error(create_plot(pb, show_cumulative = FALSE, show_restricted = TRUE), NA)
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
  # Cumulative bounds
  expect_equal(
    pb$cumulative_bounds$lower[1],
    original_results$lb/length(delta),
    tolerance = 1e-2
  )
  expect_equal(
    pb$cumulative_bounds$upper[1],
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
  l <- 8

  # Test with wiggly estimates and strong correlation
  data(estimates_wiggly_strong_corr)
  data(var_wiggly_strong_corr)
  
  # Run with parallel = TRUE
  set.seed(42)  # Set seed for reproducibility
  result_parallel <- plausible_bounds(
    estimates_wiggly_strong_corr[1:l],
    var_wiggly_strong_corr[1:l,1:l],
    parallel = TRUE
  )
  
  # Run with parallel = FALSE
  set.seed(42)  # Reset seed to ensure same random numbers
  result_sequential <- plausible_bounds(
    estimates_wiggly_strong_corr[1:l],
    var_wiggly_strong_corr[1:l, 1:l],
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
    result_parallel$restricted_metadata,
    result_sequential$restricted_metadata,
    tolerance = 1e-10
  )
  
  expect_equal(
    result_parallel$restricted_metadata$width,
    result_sequential$restricted_metadata$width,
    tolerance = 1e-10
  )
  
  # Test with constant estimates and IID errors
  data(estimates_constant_iid)
  data(var_constant_iid)
  
  # Run with parallel = TRUE
  set.seed(42)  # Set seed for reproducibility
  result_parallel_const <- plausible_bounds(
    estimates_constant_iid[1:l],
    var_constant_iid[1:l,1:l],
    parallel = TRUE
  )
  
  # Run with parallel = FALSE
  set.seed(42)  # Reset seed to ensure same random numbers
  result_sequential_const <- plausible_bounds(
    estimates_constant_iid[1:l],
    var_constant_iid[1:l,1:l],
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
    result_parallel_const$restricted_metadata,
    result_sequential_const$restricted_metadata,
    tolerance = 1e-10
  )
  
  expect_equal(
    result_parallel_const$restricted_metadata$width,
    result_sequential_const$restricted_metadata$width,
    tolerance = 1e-10
  )
})
