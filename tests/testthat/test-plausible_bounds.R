# Test basic functionality
test_that("plausible_bounds functions work with simple test data", {
  # Create test data
  set.seed(123)
  p <- 7
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

# Test estimates parameter validation and coercion
test_that("plausible_bounds accepts and coerces different estimate formats", {
  set.seed(123)
  p <- 5
  est_vec <- stats::rnorm(p)
  var_mat <- diag(p) * 0.1

  # Standard numeric vector (baseline)
  result_vec <- plausible_bounds(est_vec, var_mat)
  expect_s3_class(result_vec, "plausible_bounds")
  expect_equal(length(result_vec$restricted_bounds$unrestr_est), p)

  # Single-column matrix (should be coerced)
  est_col <- matrix(est_vec, ncol = 1)
  result_col <- plausible_bounds(est_col, var_mat)
  expect_s3_class(result_col, "plausible_bounds")
  expect_equal(result_col$restricted_bounds$unrestr_est, result_vec$restricted_bounds$unrestr_est)

  # Single-row matrix (should be transposed and coerced)
  est_row <- matrix(est_vec, nrow = 1)
  result_row <- plausible_bounds(est_row, var_mat)
  expect_s3_class(result_row, "plausible_bounds")
  expect_equal(result_row$restricted_bounds$unrestr_est, result_vec$restricted_bounds$unrestr_est)

  # Named numeric vector (should work)
  est_named <- est_vec
  names(est_named) <- paste0("t", 1:p)
  result_named <- plausible_bounds(est_named, var_mat)
  expect_s3_class(result_named, "plausible_bounds")
  expect_equal(result_named$restricted_bounds$unrestr_est, result_vec$restricted_bounds$unrestr_est)
})

test_that("plausible_bounds rejects invalid estimate formats", {
  set.seed(123)
  p <- 5
  var_mat <- diag(p) * 0.1

  # Multi-row, multi-column matrix
  est_matrix <- matrix(stats::rnorm(p * 3), nrow = 3, ncol = p)
  expect_error(
    plausible_bounds(est_matrix, var_mat),
    "If estimates is a matrix, it must have exactly one row or one column"
  )

  # Character vector
  est_char <- as.character(1:p)
  expect_error(
    plausible_bounds(est_char, var_mat),
    "estimates must be a numeric vector or single-row/single-column matrix"
  )

  # NULL
  expect_error(
    plausible_bounds(NULL, var_mat),
    "estimates must be a numeric vector or single-row/single-column matrix"
  )

  # Empty numeric vector
  expect_error(
    plausible_bounds(numeric(0), var_mat),
    "estimates must be a numeric vector or single-row/single-column matrix"
  )

  # List
  expect_error(
    plausible_bounds(list(1, 2, 3), var_mat),
    "estimates must be a numeric vector or single-row/single-column matrix"
  )

  # Data frame
  expect_error(
    plausible_bounds(data.frame(x = 1:p), var_mat),
    "estimates must be a numeric vector or single-row/single-column matrix"
  )
})

# Test parallel vs non-parallel results
test_that("parallel and non-parallel results are identical", {
  skip_on_cran()
  skip_if_not_installed("doParallel")
  skip_if_not_installed("foreach")

  l <- 8

  # Test with bighump estimates and corresponding variance
  data(estimates_bighump)
  data(var_bighump)

  # Run with parallel = TRUE
  set.seed(42)  # Set seed for reproducibility
  result_parallel <- plausible_bounds(
    estimates_bighump[1:l],
    var_bighump[1:l,1:l],
    parallel = TRUE
  )

  # Run with parallel = FALSE
  set.seed(42)  # Reset seed to ensure same random numbers
  result_sequential <- plausible_bounds(
    estimates_bighump[1:l],
    var_bighump[1:l, 1:l],
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
})

# Test n_cores parameter validation and edge cases
test_that("plausible_bounds validates n_cores parameter", {
  skip_if_not_installed("parallel")

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
  skip_on_cran()
  skip_if_not_installed("parallel")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("foreach")

  set.seed(456)
  p <- 8
  estimates <- stats::rnorm(p)
  var <- diag(p) * 0.1

  # Test with n_cores = 1
  set.seed(42)
  result_1core <- plausible_bounds(estimates, var, parallel = TRUE, n_cores = 1)

  # Test with n_cores = 2
  set.seed(42)
  result_2cores <- plausible_bounds(estimates, var, parallel = TRUE, n_cores = 2)


  # All should produce valid results
  expect_s3_class(result_1core, "plausible_bounds")
  expect_s3_class(result_2cores, "plausible_bounds")

  # Results should be very similar across different n_cores
  expect_equal(result_1core$restricted_bounds$horizon,
               result_2cores$restricted_bounds$horizon)
  expect_equal(result_1core$restricted_bounds$unrestr_est,
               result_2cores$restricted_bounds$unrestr_est)
  expect_equal(result_1core$restricted_bounds$restr_est,
               result_2cores$restricted_bounds$restr_est,
               tolerance = 0.01)
})
