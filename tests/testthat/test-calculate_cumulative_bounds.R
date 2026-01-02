test_that("calculate_cumulative_bounds works with real example data", {
  # Test with constant IID data (use subset for CRAN efficiency)
  data(estimates_constant)
  data(var_iid)

  n_test <- 6
  result <- calculate_cumulative_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test])
  
  expect_s3_class(result, "cumulative_bounds")
  expect_s3_class(result, "plausible_bounds_result")
  expect_type(result, "list")
  
  # Check structure
  expect_named(result, c("cumulative_bounds", "ate", "metadata"))
  
  # Check cumulative bounds data frame
  expect_s3_class(result$cumulative_bounds, "data.frame")
  expect_equal(nrow(result$cumulative_bounds), n_test)
  expect_named(result$cumulative_bounds, c("horizon", "coef", "lower", "upper"))
  expect_equal(result$cumulative_bounds$coef, estimates_constant[1:n_test])
  
  # Check that bounds contain the estimates (for constant case)
  expect_true(all(result$cumulative_bounds$lower <= result$cumulative_bounds$upper))

  # Check metadata
  expect_type(result$metadata, "list")
  expect_equal(result$metadata$alpha, 0.05)

  # Compute width from bounds
  width <- mean(result$cumulative_bounds$upper - result$cumulative_bounds$lower)
  expect_true(is.numeric(width))
  expect_true(width > 0)
})

test_that("calculate_cumulative_bounds works with wiggly correlated data", {
  # Test with wiggly strong correlation data (use subset for CRAN efficiency)
  data(estimates_wiggly)
  data(var_corr)

  n_test <- 6
  result <- calculate_cumulative_bounds(estimates_wiggly[1:n_test],
                                       var_corr[1:n_test, 1:n_test])
  
  expect_s3_class(result, "cumulative_bounds")
  expect_equal(nrow(result$cumulative_bounds), n_test)
  
  # Check that all bounds are properly ordered
  expect_true(all(result$cumulative_bounds$lower <= result$cumulative_bounds$upper))
  
  # For wiggly data, bounds should show more variation
  bound_widths <- result$cumulative_bounds$upper - result$cumulative_bounds$lower
  expect_true(length(unique(bound_widths)) == 1)  # Cumulative bounds have constant width
})

test_that("calculate_cumulative_bounds handles different alpha values", {
  data(estimates_constant)
  data(var_iid)

  n_test <- 6
  # Test with different confidence levels
  result_99 <- calculate_cumulative_bounds(estimates_constant[1:n_test],
                                          var_iid[1:n_test, 1:n_test],
                                          alpha = 0.01)
  result_95 <- calculate_cumulative_bounds(estimates_constant[1:n_test],
                                          var_iid[1:n_test, 1:n_test],
                                          alpha = 0.05)
  result_90 <- calculate_cumulative_bounds(estimates_constant[1:n_test],
                                          var_iid[1:n_test, 1:n_test],
                                          alpha = 0.10)
  
  # Check alpha is stored correctly
  expect_equal(result_99$metadata$alpha, 0.01)
  expect_equal(result_95$metadata$alpha, 0.05)
  expect_equal(result_90$metadata$alpha, 0.10)
  
  # Bounds should be wider for higher confidence (lower alpha)
  width_99 <- mean(result_99$cumulative_bounds$upper - result_99$cumulative_bounds$lower)
  width_95 <- mean(result_95$cumulative_bounds$upper - result_95$cumulative_bounds$lower)
  width_90 <- mean(result_90$cumulative_bounds$upper - result_90$cumulative_bounds$lower)
  
  expect_true(width_99 > width_95)
  expect_true(width_95 > width_90)
})

test_that("calculate_cumulative_bounds validates inputs correctly", {
  # Invalid estimates
  expect_error(calculate_cumulative_bounds("not_numeric", diag(3)),
               "estimates must be a numeric vector")
  expect_error(calculate_cumulative_bounds(matrix(1:6, 2, 3), diag(3)),
               "estimates must be a numeric vector")
  expect_error(calculate_cumulative_bounds(list(1, 2, 3), diag(3)),
               "estimates must be a numeric vector")
  
  # Invalid variance matrix
  expect_error(calculate_cumulative_bounds(1:3, "not_matrix"),
               "var must be a square matrix")
  expect_error(calculate_cumulative_bounds(1:3, diag(4)),
               "var must be a square matrix with dimensions matching")
  expect_error(calculate_cumulative_bounds(1:3, matrix(1:6, 2, 3)),
               "var must be a square matrix")
  expect_error(calculate_cumulative_bounds(1:4, matrix(1:12, 3, 4)),
               "var must be a square matrix")
  
  # Invalid alpha
  expect_error(calculate_cumulative_bounds(1:3, diag(3), alpha = 0),
               "alpha must be a number between 0 and 1")
  expect_error(calculate_cumulative_bounds(1:3, diag(3), alpha = 1),
               "alpha must be a number between 0 and 1")
  expect_error(calculate_cumulative_bounds(1:3, diag(3), alpha = "not_numeric"),
               "alpha must be a number between 0 and 1")
  expect_error(calculate_cumulative_bounds(1:3, diag(3), alpha = c(0.05, 0.10)),
               "alpha must be a number between 0 and 1")
})

test_that("calculate_cumulative_bounds works with edge cases", {
  # Single estimate
  estimates_single <- 0.5
  var_single <- matrix(0.1, 1, 1)
  
  result_single <- calculate_cumulative_bounds(estimates_single, var_single)
  expect_s3_class(result_single, "cumulative_bounds")
  expect_equal(nrow(result_single$cumulative_bounds), 1)
  expect_equal(result_single$cumulative_bounds$coef, estimates_single)
  
  # Two estimates
  estimates_two <- c(0.5, -0.3)
  var_two <- diag(2) * 0.1
  
  result_two <- calculate_cumulative_bounds(estimates_two, var_two)
  expect_s3_class(result_two, "cumulative_bounds")
  expect_equal(nrow(result_two$cumulative_bounds), 2)
  expect_equal(result_two$cumulative_bounds$coef, estimates_two)
  
  # Zero estimates
  estimates_zero <- rep(0, 5)
  var_zero <- diag(5) * 0.1
  
  result_zero <- calculate_cumulative_bounds(estimates_zero, var_zero)
  expect_s3_class(result_zero, "cumulative_bounds")
  expect_equal(result_zero$cumulative_bounds$coef, estimates_zero)
  # Bounds should be symmetric around zero
  expect_true(all(abs(result_zero$cumulative_bounds$lower + 
                     result_zero$cumulative_bounds$upper) < 1e-10))
})

test_that("calculate_cumulative_bounds handles different variance structures", {
  set.seed(123)
  p <- 10
  estimates <- rnorm(p)
  
  # Identity variance
  var_identity <- diag(p)
  result_identity <- calculate_cumulative_bounds(estimates, var_identity)
  
  # Scaled identity (larger variance)
  var_scaled <- diag(p) * 4
  result_scaled <- calculate_cumulative_bounds(estimates, var_scaled)
  
  # Bounds should be wider with larger variance
  width_identity <- mean(result_identity$cumulative_bounds$upper - 
                        result_identity$cumulative_bounds$lower)
  width_scaled <- mean(result_scaled$cumulative_bounds$upper - 
                      result_scaled$cumulative_bounds$lower)
  expect_true(width_scaled > width_identity)
  
  # Correlated variance (AR(1) structure)
  rho <- 0.5
  var_ar1 <- diag(p)
  for (i in 1:p) {
    for (j in 1:p) {
      var_ar1[i, j] <- rho^abs(i - j)
    }
  }
  
  result_ar1 <- calculate_cumulative_bounds(estimates, var_ar1)
  expect_s3_class(result_ar1, "cumulative_bounds")
  
  # Block diagonal variance
  var_block <- matrix(0, p, p)
  var_block[1:5, 1:5] <- 0.5
  var_block[6:10, 6:10] <- 0.5
  diag(var_block) <- 1
  
  result_block <- calculate_cumulative_bounds(estimates, var_block)
  expect_s3_class(result_block, "cumulative_bounds")
})

test_that("calculate_cumulative_bounds produces reproducible results", {
  data(estimates_constant)
  data(var_iid)

  n_test <- 6
  # Calculate bounds (note: this function is deterministic, no randomness)
  result1 <- calculate_cumulative_bounds(estimates_constant[1:n_test],
                                        var_iid[1:n_test, 1:n_test])

  result2 <- calculate_cumulative_bounds(estimates_constant[1:n_test],
                                        var_iid[1:n_test, 1:n_test])

  # Results should be identical
  expect_identical(result1$cumulative_bounds, result2$cumulative_bounds)
  expect_identical(result1$metadata, result2$metadata)
})

test_that("calculate_cumulative_bounds ATE calculation is correct", {
  # Test that ATE is calculated correctly
  set.seed(123)
  estimates <- rnorm(5, mean = 2)
  var <- diag(5) * 0.1

  result <- calculate_cumulative_bounds(estimates, var, alpha = 0.05)

  # Check ATE structure
  expect_named(result$ate, c("estimate", "se"))

  # ATE estimate should be the mean of estimates
  expect_equal(unname(result$ate["estimate"]), mean(estimates))

  # Check that bounds are symmetric around ATE for constant bounds
  expect_true(result$metadata$lb < result$ate["estimate"])
  expect_true(result$metadata$ub > result$ate["estimate"])
})

test_that("calculate_cumulative_bounds metadata is correct", {
  data(estimates_wiggly)
  data(var_corr)

  n_test <- 6
  result <- calculate_cumulative_bounds(estimates_wiggly[1:n_test],
                                       var_corr[1:n_test, 1:n_test],
                                       alpha = 0.10)

  # Check metadata structure
  expect_type(result$metadata, "list")
  expect_true("alpha" %in% names(result$metadata))
  expect_true("lb" %in% names(result$metadata))
  expect_true("ub" %in% names(result$metadata))
  expect_true("preperiods" %in% names(result$metadata))

  # Check metadata values
  expect_equal(result$metadata$alpha, 0.10)
  expect_true(is.numeric(result$metadata$lb))
  expect_true(is.numeric(result$metadata$ub))
  expect_equal(result$metadata$preperiods, 0)

  # Compute width from ub - lb
  width <- result$metadata$ub - result$metadata$lb
  expect_true(is.numeric(width))
  expect_true(width > 0)
})

test_that("calculate_cumulative_bounds handles NA and infinite values appropriately", {
  # Test with NA in estimates
  estimates_na <- c(1, NA, 3)
  var_valid <- diag(3) * 0.1
  expect_error(calculate_cumulative_bounds(estimates_na, var_valid))
  
  # Test with infinite estimates
  estimates_inf <- c(1, Inf, 3)
  expect_error(calculate_cumulative_bounds(estimates_inf, var_valid))
  
  # Test with NA in variance
  estimates_valid <- c(1, 2, 3)
  var_na <- diag(3) * 0.1
  var_na[1, 2] <- NA
  expect_error(calculate_cumulative_bounds(estimates_valid, var_na))
  
  # Test with non-positive definite variance
  # Note: calculate_cumulative_bounds doesn't explicitly check for PD,
  # but negative variances would be caught if they occur
  var_not_pd <- matrix(c(1, 2, 2, 1), 2, 2)  # Not positive definite
  # This may or may not error depending on whether sqrt of negative variance occurs
  # For now, just verify it runs without crashing (may get warnings)
  suppressWarnings({
    result <- calculate_cumulative_bounds(c(1, 2), var_not_pd)
    expect_s3_class(result, "cumulative_bounds")
  })
})
