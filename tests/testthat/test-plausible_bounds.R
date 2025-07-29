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
  expect_equal(nrow(cumul_bounds$bounds), p)
  expect_true("horizon" %in% names(cumul_bounds$bounds))
  expect_true("coef" %in% names(cumul_bounds$bounds))
  expect_true("lower" %in% names(cumul_bounds$bounds))
  expect_true("upper" %in% names(cumul_bounds$bounds))
  expect_true("pointwise_lower" %in% names(cumul_bounds$bounds))
  expect_true("pointwise_upper" %in% names(cumul_bounds$bounds))
  expect_true("supt_lower" %in% names(cumul_bounds$bounds))
  expect_true("supt_upper" %in% names(cumul_bounds$bounds))

  # Test calculate_restricted_bounds
  restr_bounds <- calculate_restricted_bounds(estimates, var)
  expect_s3_class(restr_bounds, "restricted_bounds")
  expect_s3_class(restr_bounds, "plausible_bounds_result")
  expect_equal(nrow(restr_bounds$bounds), p)
  expect_true("horizon" %in% names(restr_bounds$bounds))
  expect_true("coef" %in% names(restr_bounds$bounds))
  expect_true("surrogate" %in% names(restr_bounds$bounds))
  expect_true("lower" %in% names(restr_bounds$bounds))
  expect_true("upper" %in% names(restr_bounds$bounds))
  expect_true("pointwise_lower" %in% names(restr_bounds$bounds))
  expect_true("pointwise_upper" %in% names(restr_bounds$bounds))
  expect_true("supt_lower" %in% names(restr_bounds$bounds))
  expect_true("supt_upper" %in% names(restr_bounds$bounds))
  
  # Test plausible_bounds
  pb <- plausible_bounds(estimates, var)
  expect_s3_class(pb, "plausible_bounds")
  expect_s3_class(pb$cumulative_bounds, "cumulative_bounds")
  expect_s3_class(pb$restricted_bounds, "restricted_bounds")
  
  # Test plot_bounds (just check that it runs without error)
  expect_error(plot_bounds(pb), NA)
  expect_error(plot_bounds(pb, which = "cumulative"), NA)
  expect_error(plot_bounds(pb, which = "restricted"), NA)
  expect_error(plot_bounds(cumul_bounds), NA)
  expect_error(plot_bounds(restr_bounds), NA)
  expect_error(plot_bounds(cumul_bounds, restr_bounds), NA)
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
      pb$cumulative_bounds$bounds$lower[1],
      original_results$lb,
      tolerance = 1e-6
    )
    expect_equal(
      pb$cumulative_bounds$bounds$upper[1],
      original_results$ub,
      tolerance = 1e-6
    )
    
    # Restricted bounds
    expect_equal(
      pb$restricted_bounds$bounds$surrogate,
      as.vector(original_results$surrogate),
      tolerance = 1e-6
    )
    expect_equal(
      pb$restricted_bounds$bounds$lower,
      as.vector(original_results$restricted_LB),
      tolerance = 1e-6
    )
    expect_equal(
      pb$restricted_bounds$bounds$upper,
      as.vector(original_results$restricted_UB),
      tolerance = 1e-6
    )
  }
})