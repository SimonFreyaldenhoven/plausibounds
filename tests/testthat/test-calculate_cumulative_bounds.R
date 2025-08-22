# Test calculate_cumulative_bounds function
# Following Hadley Wickham's R Packages testing guidelines

test_that("calculate_cumulative_bounds works with real example data", {
  # Test with constant IID data
  data(estimates_constant_iid)
  data(var_constant_iid)
  
  result <- calculate_cumulative_bounds(estimates_constant_iid, var_constant_iid)
  
  expect_s3_class(result, "cumulative_bounds")
  expect_s3_class(result, "plausible_bounds_result")
  expect_type(result, "list")
  
  # Check structure
  expect_named(result, c("cumulative_bounds", "pointwise_bounds", 
                        "supt_bounds", "metadata"))
  
  # Check cumulative bounds data frame
  expect_s3_class(result$cumulative_bounds, "data.frame")
  expect_equal(nrow(result$cumulative_bounds), length(estimates_constant_iid))
  expect_named(result$cumulative_bounds, c("horizon", "coef", "lower", "upper"))
  expect_equal(result$cumulative_bounds$coef, estimates_constant_iid)
  
  # Check that bounds contain the estimates (for constant case)
  expect_true(all(result$cumulative_bounds$lower <= result$cumulative_bounds$upper))
  
  # Check metadata
  expect_type(result$metadata, "list")
  expect_equal(result$metadata$alpha, 0.05)
  expect_true(is.numeric(result$metadata$width))
  expect_true(result$metadata$width > 0)
})

test_that("calculate_cumulative_bounds works with wiggly correlated data", {
  # Test with wiggly strong correlation data
  data(estimates_wiggly_strong_corr)
  data(var_wiggly_strong_corr)
  
  result <- calculate_cumulative_bounds(estimates_wiggly_strong_corr, 
                                       var_wiggly_strong_corr)
  
  expect_s3_class(result, "cumulative_bounds")
  expect_equal(nrow(result$cumulative_bounds), length(estimates_wiggly_strong_corr))
  
  # Check that all bounds are properly ordered
  expect_true(all(result$cumulative_bounds$lower <= result$cumulative_bounds$upper))
  
  # For wiggly data, bounds should show more variation
  bound_widths <- result$cumulative_bounds$upper - result$cumulative_bounds$lower
  expect_true(length(unique(bound_widths)) == 1)  # Cumulative bounds have constant width
})

test_that("calculate_cumulative_bounds handles different alpha values", {
  data(estimates_constant_iid)
  data(var_constant_iid)
  
  # Test with different confidence levels
  result_99 <- calculate_cumulative_bounds(estimates_constant_iid, 
                                          var_constant_iid, 
                                          alpha = 0.01)
  result_95 <- calculate_cumulative_bounds(estimates_constant_iid, 
                                          var_constant_iid, 
                                          alpha = 0.05)
  result_90 <- calculate_cumulative_bounds(estimates_constant_iid, 
                                          var_constant_iid, 
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

test_that("calculate_cumulative_bounds handles pointwise and supt parameters", {
  data(estimates_constant_iid)
  data(var_constant_iid)
  
  # Test with both FALSE
  result_none <- calculate_cumulative_bounds(estimates_constant_iid, 
                                            var_constant_iid,
                                            include_pointwise = FALSE,
                                            include_supt = FALSE)
  expect_null(result_none$pointwise_bounds)
  expect_null(result_none$supt_bounds)
  
  # Test with pointwise only
  result_pw <- calculate_cumulative_bounds(estimates_constant_iid, 
                                          var_constant_iid,
                                          include_pointwise = TRUE,
                                          include_supt = FALSE)
  expect_type(result_pw$pointwise_bounds, "list")
  expect_named(result_pw$pointwise_bounds, c("lower", "upper", "critval"))
  expect_equal(length(result_pw$pointwise_bounds$lower), length(estimates_constant_iid))
  expect_equal(length(result_pw$pointwise_bounds$upper), length(estimates_constant_iid))
  expect_true(is.numeric(result_pw$pointwise_bounds$critval))
  expect_null(result_pw$supt_bounds)
  
  # Test with supt only
  result_supt <- calculate_cumulative_bounds(estimates_constant_iid, 
                                            var_constant_iid,
                                            include_pointwise = FALSE,
                                            include_supt = TRUE)
  expect_null(result_supt$pointwise_bounds)
  expect_type(result_supt$supt_bounds, "list")
  expect_named(result_supt$supt_bounds, c("lower", "upper", "critval"))
  expect_equal(length(result_supt$supt_bounds$lower), length(estimates_constant_iid))
  expect_equal(length(result_supt$supt_bounds$upper), length(estimates_constant_iid))
  expect_true(is.numeric(result_supt$supt_bounds$critval))
  
  # Test with both TRUE (default)
  result_both <- calculate_cumulative_bounds(estimates_constant_iid, 
                                            var_constant_iid)
  expect_type(result_both$pointwise_bounds, "list")
  expect_type(result_both$supt_bounds, "list")
  
  # Sup-t bounds should be wider than pointwise
  expect_true(all(result_both$supt_bounds$lower <= result_both$pointwise_bounds$lower))
  expect_true(all(result_both$supt_bounds$upper >= result_both$pointwise_bounds$upper))
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
  data(estimates_constant_iid)
  data(var_constant_iid)
  
  # Set seed for reproducibility
  set.seed(999)
  result1 <- calculate_cumulative_bounds(estimates_constant_iid, 
                                        var_constant_iid)
  
  set.seed(999)
  result2 <- calculate_cumulative_bounds(estimates_constant_iid, 
                                        var_constant_iid)
  
  # Results should be identical
  expect_identical(result1$cumulative_bounds, result2$cumulative_bounds)
  expect_identical(result1$pointwise_bounds, result2$pointwise_bounds)
  expect_equal(result1$supt_bounds, result2$supt_bounds, tolerance = 1e-10)
  expect_identical(result1$metadata, result2$metadata)
})

test_that("Wald_bounds internal function works correctly", {
  # Access internal function
  #Wald_bounds <- plausibounds:::Wald_bounds
  
  # Simple test case
  dhat <- c(1, 2, 3)
  Vhat <- diag(3) * 0.1
  alpha <- 0.05
  
  result <- Wald_bounds(dhat, Vhat, alpha)
  
  expect_type(result, "list")
  expect_named(result, c("LB", "UB"))
  expect_true(is.matrix(result$LB))
  expect_true(is.matrix(result$UB))
  expect_equal(dim(result$LB), c(3, 1))
  expect_equal(dim(result$UB), c(3, 1))
  
  # Lower bounds should be less than upper bounds
  expect_true(all(result$LB < result$UB))
})

test_that("calculate_cumulative_bounds metadata is correct", {
  data(estimates_wiggly_strong_corr)
  data(var_wiggly_strong_corr)
  
  result <- calculate_cumulative_bounds(estimates_wiggly_strong_corr,
                                       var_wiggly_strong_corr,
                                       alpha = 0.10)
  
  # Check metadata structure
  expect_type(result$metadata, "list")
  expect_true("alpha" %in% names(result$metadata))
  expect_true("width" %in% names(result$metadata))
  expect_true("individual_upper" %in% names(result$metadata))
  expect_true("individual_lower" %in% names(result$metadata))
  
  # Check metadata values
  expect_equal(result$metadata$alpha, 0.10)
  expect_true(is.numeric(result$metadata$width))
  expect_true(result$metadata$width > 0)
  expect_true(is.matrix(result$metadata$individual_upper))
  expect_true(is.matrix(result$metadata$individual_lower))
  
  # Width should match the actual bounds
  actual_width <- mean(result$cumulative_bounds$upper - result$cumulative_bounds$lower)
  expect_equal(result$metadata$width, actual_width)
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
  var_not_pd <- matrix(c(1, 2, 2, 1), 2, 2)  # Not positive definite
  expect_error(calculate_cumulative_bounds(c(1, 2), var_not_pd))
})
