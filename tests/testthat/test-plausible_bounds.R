library(testthat)

# Test data
test_that("plausible_bounds functions work with simple test data", {
  # Create test data
  set.seed(123)
  p <- 10
  estimates <- rnorm(p)
  var <- diag(p) * 0.1
  
  # Test calculate_cumulative_bounds
  cumul_bounds <- calculate_cumulative_bounds(estimates, var)
  expect_s3_class(cumul_bounds, "cumulative_bounds")
  expect_s3_class(cumul_bounds, "plausible_bounds_result")
  expect_equal(nrow(cumul_bounds$cumulative_bounds), p)
  expect_true("horizon" %in% names(cumul_bounds$cumulative_bounds))
  expect_true("coef" %in% names(cumul_bounds$cumulative_bounds))
  expect_true("lower" %in% names(cumul_bounds$cumulative_bounds))
  expect_true("upper" %in% names(cumul_bounds$cumulative_bounds))
  expect_true(!is.null(cumul_bounds$pointwise_bounds))
  expect_true(!is.null(cumul_bounds$supt_bounds))
  expect_true(!is.null(cumul_bounds$metadata))

  # Test calculate_restricted_bounds
  restr_bounds <- calculate_restricted_bounds(estimates, var)
  expect_s3_class(restr_bounds, "restricted_bounds")
  expect_s3_class(restr_bounds, "plausible_bounds_result")
  expect_equal(nrow(restr_bounds$restricted_bounds), p)
  expect_true("horizon" %in% names(restr_bounds$restricted_bounds))
  expect_true("coef" %in% names(restr_bounds$restricted_bounds))
  expect_true("surrogate" %in% names(restr_bounds$restricted_bounds))
  expect_true("lower" %in% names(restr_bounds$restricted_bounds))
  expect_true("upper" %in% names(restr_bounds$restricted_bounds))
  expect_true(!is.null(restr_bounds$metadata))
  
  # Test plausible_bounds
  pb <- plausible_bounds(estimates, var)
  expect_s3_class(pb, "plausible_bounds")
  expect_s3_class(pb$cumulative_bounds, "cumulative_bounds")
  expect_s3_class(pb$restricted_bounds, "restricted_bounds")
  
  # Test create_plot (just check that it runs without error)
  expect_error(create_plot(pb), NA)
  expect_error(create_plot(pb, show_cumulative = TRUE, show_restricted = FALSE), NA)
  expect_error(create_plot(pb, show_cumulative = FALSE, show_restricted = TRUE), NA)
  expect_error(create_plot(cumul_bounds), NA)
  expect_error(create_plot(restr_bounds), NA)
  expect_error(create_plot(cumul_bounds, restr_bounds), NA)
})

# Test comparison with original function
test_that("new functions match original full_eventplot_l2tf results", {
  # Skip this test if the test data files don't exist
  skip_if_not(file.exists("../delta.csv") && file.exists("../vhat.csv"))
  
  # Load test data
  delta <- as.matrix(read.csv("tests/delta.csv", header = FALSE))
  vhat <- as.matrix(read.csv("tests/vhat.csv", header = FALSE))
  
  # Run original function (if available)
  if (exists("full_eventplot_l2tf")) {
    original_results <- full_eventplot_l2tf(delta, vhat, nsim = 1000)
    
    # Run new functions
    pb <- plausible_bounds(as.vector(delta), vhat)
    
    # Compare results
    # Cumulative bounds
    expect_equal(
      pb$cumulative_bounds$lower[1],
      original_results$lb,
      tolerance = 1e-6
    )
    expect_equal(
      pb$cumulative_bounds$upper[1],
      original_results$ub,
      tolerance = 1e-6
    )
    
    # Restricted bounds
    expect_equal(
      pb$restricted_bounds$surrogate,
      as.vector(original_results$surrogate),
      tolerance = 1e-6
    )
    expect_equal(
      pb$restricted_bounds$lower,
      as.vector(original_results$restricted_LB),
      tolerance = 1e-6
    )
    expect_equal(
      pb$restricted_bounds$upper,
      as.vector(original_results$restricted_UB),
      tolerance = 1e-6
    )
  }
})

# Test parallel vs non-parallel results
test_that("parallel and non-parallel results are identical", {
  # Test with wiggly estimates and strong correlation
  data(estimates_wiggly_strong_corr)
  data(var_wiggly_strong_corr)
  
  # Run with parallel = TRUE
  set.seed(42)  # Set seed for reproducibility
  result_parallel <- plausible_bounds(
    estimates_wiggly_strong_corr,
    var_wiggly_strong_corr,
    parallel = TRUE
  )
  
  # Run with parallel = FALSE
  set.seed(42)  # Reset seed to ensure same random numbers
  result_sequential <- plausible_bounds(
    estimates_wiggly_strong_corr,
    var_wiggly_strong_corr,
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
    estimates_constant_iid,
    var_constant_iid,
    parallel = TRUE
  )
  
  # Run with parallel = FALSE
  set.seed(42)  # Reset seed to ensure same random numbers
  result_sequential_const <- plausible_bounds(
    estimates_constant_iid,
    var_constant_iid,
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