
test_that("calculate_restricted_bounds works with real example data", {
  # Use smaller subset for computational efficiency
  data(estimates_constant_iid)
  data(var_constant_iid)
  
  # Use first 10 observations for faster testing
  n_test <- min(10, length(estimates_constant_iid))
  estimates_subset <- estimates_constant_iid[1:n_test]
  var_subset <- var_constant_iid[1:n_test, 1:n_test]
  
  result <- calculate_restricted_bounds(estimates_subset, var_subset)
  
  expect_s3_class(result, "restricted_bounds")
  expect_s3_class(result, "plausible_bounds_result")
  expect_type(result, "list")
    
  # Check restricted bounds data frame
  expect_s3_class(result$restricted_bounds, "data.frame")
  expect_equal(nrow(result$restricted_bounds), n_test)
  expect_named(result$restricted_bounds, c("horizon", "coef", "surrogate", "lower", "upper"))
  expect_equal(result$restricted_bounds$coef, estimates_subset)
  
  # Check that bounds are properly ordered
  expect_true(all(result$restricted_bounds$lower <= result$restricted_bounds$upper))
  
  # Check metadata
  expect_type(result$metadata, "list")
  expect_equal(result$metadata$alpha, 0.05)
  expect_true(is.numeric(result$metadata$width))
  expect_true(result$metadata$width > 0)
})

test_that("calculate_restricted_bounds works with wiggly correlated data", {
  # Test with smaller subset of wiggly strong correlation data
  data(estimates_wiggly_strong_corr)
  data(var_wiggly_strong_corr)
  
  # Use first 8 observations for faster testing
  n_test <- min(8, length(estimates_wiggly_strong_corr))
  estimates_subset <- estimates_wiggly_strong_corr[1:n_test]
  var_subset <- var_wiggly_strong_corr[1:n_test, 1:n_test]
  
  result <- calculate_restricted_bounds(estimates_subset, var_subset)
  
  expect_s3_class(result, "restricted_bounds")
  expect_equal(nrow(result$restricted_bounds), n_test)
  
  # Check that surrogate is smooth (restricted)
  surrogate_diff <- diff(result$restricted_bounds$surrogate)
  # Surrogate should be relatively smooth compared to original estimates
  estimates_diff <- diff(estimates_subset)
  expect_true(var(surrogate_diff) <= var(estimates_diff))
  
  # Bounds should contain the surrogate
  expect_true(all(result$restricted_bounds$lower <= result$restricted_bounds$surrogate))
  expect_true(all(result$restricted_bounds$surrogate <= result$restricted_bounds$upper))
})

test_that("calculate_restricted_bounds handles different alpha values", {
  # Use small sample for speed
  set.seed(123)
  n <- 5
  estimates <- rnorm(n)
  var <- diag(n) * 0.1
  
  # Test with different confidence levels
  result_99 <- calculate_restricted_bounds(estimates, var, alpha = 0.01)
  result_95 <- calculate_restricted_bounds(estimates, var, alpha = 0.05)
  result_90 <- calculate_restricted_bounds(estimates, var, alpha = 0.10)
  
  # Check alpha is stored correctly
  expect_equal(result_99$metadata$alpha, 0.01)
  expect_equal(result_95$metadata$alpha, 0.05)
  expect_equal(result_90$metadata$alpha, 0.10)
  
  # Bounds should be wider for higher confidence (lower alpha)
  width_99 <- mean(result_99$restricted_bounds$upper - result_99$restricted_bounds$lower)
  width_95 <- mean(result_95$restricted_bounds$upper - result_95$restricted_bounds$lower)
  width_90 <- mean(result_90$restricted_bounds$upper - result_90$restricted_bounds$lower)
  
  expect_true(width_99 > width_95)
  expect_true(width_95 > width_90)
})

test_that("calculate_restricted_bounds handles parallel parameter", {
  skip_if_not_installed("parallel")
  
  # Use small sample for speed
  set.seed(456)
  n <- 10
  estimates <- rnorm(n)
  var <- diag(n) * 0.1
  
  # Test with parallel = FALSE
  set.seed(100)
  result_seq <- calculate_restricted_bounds(estimates, var, parallel = FALSE)
  
  # Test with parallel = TRUE
  set.seed(100)
  result_par <- calculate_restricted_bounds(estimates, var, parallel = TRUE)
  
  # Results should be very similar
  expect_equal(result_seq$restricted_bounds$horizon, 
               result_par$restricted_bounds$horizon)
  expect_equal(result_seq$restricted_bounds$coef, 
               result_par$restricted_bounds$coef)
  
  # Surrogate and bounds might have small numerical differences
  expect_equal(result_seq$restricted_bounds$surrogate,
               result_par$restricted_bounds$surrogate,
               tolerance = 0.01)
})

test_that("calculate_restricted_bounds validates inputs correctly", {
  # Single estimate should error (need at least 2 for restrictions)
  expect_error(calculate_restricted_bounds(0.5, matrix(.1)),
               "calculate_restricted_bounds requires at least 2 estimates")
  
  # Invalid estimates
  expect_error(calculate_restricted_bounds("not_numeric", diag(3)),
               "estimates must be a numeric vector")
  expect_error(calculate_restricted_bounds(matrix(1:6, 2, 3), diag(3)),
               "estimates must be a numeric vector")
  
  # Invalid variance matrix
  expect_error(calculate_restricted_bounds(1:3, "not_matrix"),
               "var must be a square matrix")
  expect_error(calculate_restricted_bounds(1:3, diag(4)),
               "var must be a square matrix with dimensions matching")
  expect_error(calculate_restricted_bounds(1:3, matrix(1:6, 2, 3)),
               "var must be a square matrix")
  
  # Invalid alpha
  expect_error(calculate_restricted_bounds(1:3, diag(3), alpha = 0),
               "alpha must be a number between 0 and 1")
  expect_error(calculate_restricted_bounds(1:3, diag(3), alpha = 1),
               "alpha must be a number between 0 and 1")
  expect_error(calculate_restricted_bounds(1:3, diag(3), alpha = "not_numeric"),
               "alpha must be a number between 0 and 1")
})

test_that("calculate_restricted_bounds works with edge cases", {

  # Two estimates
  estimates_two <- c(0.5, -0.3)
  var_two <- diag(2) * 0.1
  
  result_two <- calculate_restricted_bounds(estimates_two, var_two)
  expect_s3_class(result_two, "restricted_bounds")
  expect_equal(nrow(result_two$restricted_bounds), 2)
  expect_equal(result_two$restricted_bounds$coef, estimates_two)
  
  # Three estimates (minimum for meaningful restriction)
  estimates_three <- c(0.5, -0.3, 0.2)
  var_three <- diag(3) * 0.1
  
  result_three <- calculate_restricted_bounds(estimates_three, var_three)
  expect_s3_class(result_three, "restricted_bounds")
  expect_equal(nrow(result_three$restricted_bounds), 3)
})

test_that("calculate_restricted_bounds handles different variance structures", {
  set.seed(789)
  n <- 5  # Small sample for speed
  estimates <- rnorm(n)
  
  # Identity variance
  var_identity <- diag(n)
  result_identity <- calculate_restricted_bounds(estimates, var_identity)
  
  # Scaled identity (larger variance)
  var_scaled <- diag(n) * 4
  result_scaled <- calculate_restricted_bounds(estimates, var_scaled)
  
  # Bounds should be wider with larger variance
  width_identity <- mean(result_identity$restricted_bounds$upper - 
                        result_identity$restricted_bounds$lower)
  width_scaled <- mean(result_scaled$restricted_bounds$upper - 
                      result_scaled$restricted_bounds$lower)
  expect_true(width_scaled > width_identity)
  
  # Correlated variance
  var_corr <- diag(n) * 0.5
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      var_corr[i, j] <- var_corr[j, i] <- 0.2 * exp(-abs(i - j) / 2)
    }
  }
  
  result_corr <- calculate_restricted_bounds(estimates, var_corr)
  expect_s3_class(result_corr, "restricted_bounds")
})

test_that("calculate_restricted_bounds produces smooth surrogates", {
  # Generate wiggly estimates
  set.seed(111)
  n <- 8
  t <- seq(0, 1, length.out = n)
  estimates <- sin(4 * pi * t) + rnorm(n, 0, 0.2)
  var <- diag(n) * 0.1
  
  result <- calculate_restricted_bounds(estimates, var)
  
  # Surrogate should be smoother than original estimates
  surrogate_roughness <- sum(diff(result$restricted_bounds$surrogate)^2)
  estimates_roughness <- sum(diff(estimates)^2)
  
  expect_true(surrogate_roughness < estimates_roughness)
  
  # Check that surrogate is within reasonable range of estimates
  expect_true(all(abs(result$restricted_bounds$surrogate) <= 
                 max(abs(estimates)) + 2 * sqrt(max(diag(var)))))
})

test_that("calculate_restricted_bounds metadata is correct", {
  set.seed(222)
  n <- 6
  estimates <- rnorm(n)
  var <- diag(n) * 0.15
  
  result <- calculate_restricted_bounds(estimates, var, alpha = 0.10)
  
  # Check metadata structure
  expect_type(result$metadata, "list")
  expect_true("alpha" %in% names(result$metadata))
  expect_true("width" %in% names(result$metadata))
  expect_true("suptb" %in% names(result$metadata))
  expect_true("df" %in% names(result$metadata))
  
  # Check metadata values
  expect_equal(result$metadata$alpha, 0.10)
  expect_true(is.numeric(result$metadata$width))
  expect_true(result$metadata$width > 0)
  expect_true(is.numeric(result$metadata$suptb))
  expect_true(result$metadata$suptb > 0)
  expect_true(is.numeric(result$metadata$df))
  
  # Width should match the actual bounds
  actual_width <- mean(result$restricted_bounds$upper - result$restricted_bounds$lower)
  expect_equal(result$metadata$width, actual_width, tolerance = 1e-10)
})

test_that("calculate_restricted_bounds handles pointwise and supt parameters", {
  # Use small sample
  set.seed(333)
  n <- 5
  estimates <- rnorm(n)
  var <- diag(n) * 0.1
  
  # Test with both FALSE (default for restricted bounds)
  result_none <- calculate_restricted_bounds(estimates, var,
                                            include_pointwise = FALSE,
                                            include_supt = FALSE)
  expect_null(result_none$pointwise_bounds)
  expect_null(result_none$supt_bounds)
  
  # Test with pointwise only
  result_pw <- calculate_restricted_bounds(estimates, var,
                                          include_pointwise = TRUE,
                                          include_supt = FALSE)
  expect_type(result_pw$pointwise_bounds, "list")
  expect_named(result_pw$pointwise_bounds, c("lower", "upper", "critval"))
  expect_null(result_pw$supt_bounds)
  
  # Test with supt only
  result_supt <- calculate_restricted_bounds(estimates, var,
                                            include_pointwise = FALSE,
                                            include_supt = TRUE)
  expect_null(result_supt$pointwise_bounds)
  expect_type(result_supt$supt_bounds, "list")
  expect_named(result_supt$supt_bounds, c("lower", "upper", "critval"))
  
  # Test with both TRUE
  result_both <- calculate_restricted_bounds(estimates, var,
                                            include_pointwise = TRUE,
                                            include_supt = TRUE)
  expect_type(result_both$pointwise_bounds, "list")
  expect_type(result_both$supt_bounds, "list")
})

test_that("calculate_restricted_bounds optimization works correctly", {
  # Test that optimization finds a reasonable solution
  set.seed(444)
  n <- 6
  # Create estimates with clear trend
  estimates <- c(1, 1.5, 2, 2.3, 2.5, 2.6)
  var <- diag(n) * 0.1
  
  result <- calculate_restricted_bounds(estimates, var)
  
  # Surrogate should capture the trend
  expect_true(all(diff(result$restricted_bounds$surrogate) >= -0.5))
  
  # Check that metadata contains expected fields
  expect_true(!is.null(result$metadata$surrogate_class))
  expect_true(is.character(result$metadata$surrogate_class))
})

test_that("calculate_restricted_bounds handles constant estimates", {
  # All estimates the same
  n <- 5
  estimates_const <- rep(2, n)
  var <- diag(n) * 0.1
  
  result <- calculate_restricted_bounds(estimates_const, var)
  
  # Surrogate should be approximately constant
  expect_true(all(abs(diff(result$restricted_bounds$surrogate)) < 0.1))
  
  # Bounds should be relatively uniform
  widths <- result$restricted_bounds$upper - result$restricted_bounds$lower
  expect_true(sd(widths) / mean(widths) < 0.1)
})

test_that("calculate_restricted_bounds handles zero variance edge case", {
  # Very small variance
  n <- 4
  estimates <- rnorm(n)
  var_small <- diag(n) * 1e-7
  
  # Should still work but with very narrow bounds
  result <- calculate_restricted_bounds(estimates, var_small)
  expect_s3_class(result, "restricted_bounds")
  
  widths <- result$restricted_bounds$upper - result$restricted_bounds$lower
  expect_true(all(widths < 0.01))
})

test_that("calculate_restricted_bounds produces reproducible results", {
  set.seed(555)
  n <- 5
  estimates <- rnorm(n)
  var <- diag(n) * 0.1
  
  set.seed(999)
  result1 <- calculate_restricted_bounds(estimates, var)
  
  set.seed(999)
  result2 <- calculate_restricted_bounds(estimates, var)
  
  expect_identical(result1$restricted_bounds, result2$restricted_bounds)
  expect_identical(result1$metadata, result2$metadata)
})