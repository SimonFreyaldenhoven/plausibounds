# Tests for preperiod functionality
# These tests validate that the package correctly handles pre-treatment periods
# in event study designs, including horizon numbering, Wald tests, and plotting.

# Helper function to create block diagonal matrices
block_diag_matrix <- function(A, B) {
  out <- matrix(0, nrow(A) + nrow(B), ncol(A) + ncol(B))
  out[1:nrow(A), 1:ncol(A)] <- A
  out[(nrow(A)+1):(nrow(A)+nrow(B)), (ncol(A)+1):(ncol(A)+ncol(B))] <- B
  out
}

# Section 1: Input Validation Tests ----

test_that("preperiods parameter accepts valid values", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  n_test <- 7
  # Valid: zero preperiods
  expect_error(plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = 0), NA)

  # Valid: single preperiod
  pre_single <- rnorm(1)
  Vpre_single <- matrix(0.1, 1, 1)
  estimates <- c(pre_single, estimates_constant[1:n_test])
  var <- block_diag_matrix(Vpre_single, var_iid[1:n_test, 1:n_test])
  expect_error(plausible_bounds(estimates, var, preperiods = 1), NA)

  # Valid: multiple preperiods (8 pre + 7 post)
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  estimates_full <- c(pre_mean0, estimates_constant[1:7])
  var_full <- block_diag_matrix(Vpre_iid, var_iid[1:7, 1:7])
  expect_error(plausible_bounds(estimates_full, var_full, preperiods = 8), NA)
})

test_that("preperiods parameter rejects negative values", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  n_test <- 6
  expect_error(
    plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = -1),
    "preperiods must be a non-negative integer"
  )

  expect_error(
    plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = -5),
    "preperiods must be a non-negative integer"
  )
})

test_that("preperiods parameter rejects non-integer values", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  n_test <- 6
  expect_error(
    plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = 2.5),
    "preperiods must be a non-negative integer"
  )

  expect_error(
    plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = 3.7),
    "preperiods must be a non-negative integer"
  )
})

test_that("preperiods parameter rejects values >= length(estimates)", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  n_test <- 6
  expect_error(
    plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = n_test),
    "preperiods must be less than the length of estimates"
  )

  expect_error(
    plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = n_test + 1),
    "preperiods must be less than the length of estimates"
  )
})

test_that("preperiods parameter rejects invalid types", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  n_test <- 6
  # Non-numeric string
  expect_error(
    plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = "two"),
    "preperiods must be a non-negative integer"
  )

  # NULL
  expect_error(
    plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = NULL),
    "preperiods must be a non-negative integer"
  )

  # NA
  expect_error(
    plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = NA),
    "preperiods must be a non-negative integer"
  )
})

test_that("preperiods parameter rejects vector inputs", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  n_test <- 6
  expect_error(
    plausible_bounds(estimates_constant[1:n_test], var_iid[1:n_test, 1:n_test], preperiods = c(1, 2)),
    "preperiods must be a non-negative integer"
  )
})

# Section 2: Core Functionality Tests ----

test_that("core functionality with 8 preperiods (mean0 scenario)", {
  # Setup data once for all assertions
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  npre <- 8
  npost <- 7
  estimates <- c(pre_mean0[1:npre], estimates_constant[1:npost])
  var <- block_diag_matrix(Vpre_iid[1:npre, 1:npre], var_iid[1:npost, 1:npost])

  # Call functions once (used across multiple assertions below)
  result <- plausible_bounds(estimates, var, preperiods = npre)
  restr <- calculate_restricted_bounds(estimates, var, preperiods = npre)
  cumul <- calculate_cumulative_bounds(estimates, var, preperiods = npre)

  # Test 1: Function accepts multiple preperiods without error
  expect_s3_class(result, "plausible_bounds")

  # Test 2: Horizon numbering in restricted bounds
  expect_equal(min(result$restricted_bounds$horizon), -npre)
  expect_equal(max(result$restricted_bounds$horizon), npost)
  expect_true(all(result$restricted_bounds$horizon[1:npre] < 0))
  expect_true(all(result$restricted_bounds$horizon[(npre+1):(npre+npost)] > 0))
  expect_false(any(result$restricted_bounds$horizon == 0))  # No period 0
  expected_horizons <- c(-npre:-1, 1:npost)
  expect_equal(result$restricted_bounds$horizon, expected_horizons)

  # Test 3: Wpre is computed when preperiods > 0
  expect_true(!is.null(result$wald_test$pre))
  expect_true(is.numeric(result$wald_test$pre$statistic))
  expect_true(is.numeric(result$wald_test$pre$p_value))
  expect_true(result$wald_test$pre$p_value >= 0 && result$wald_test$pre$p_value <= 1)

  # Test 4: Restricted bounds are NA for preperiods but not for post-periods
  pre_rows <- result$restricted_bounds$horizon < 0
  expect_true(all(is.na(result$restricted_bounds$surrogate[pre_rows])))
  expect_true(all(is.na(result$restricted_bounds$lower[pre_rows])))
  expect_true(all(is.na(result$restricted_bounds$upper[pre_rows])))
  post_rows <- result$restricted_bounds$horizon > 0
  expect_false(any(is.na(result$restricted_bounds$surrogate[post_rows])))
  expect_false(any(is.na(result$restricted_bounds$lower[post_rows])))
  expect_false(any(is.na(result$restricted_bounds$upper[post_rows])))
  expect_false(any(is.na(result$restricted_bounds$coef)))

  # Test 5: Wpre passes for mean0 scenario (no pre-trends)
  expect_true(result$wald_test$pre$p_value >= 0.05)

  # Test 6: Horizon numbering consistent across all functions
  expect_equal(restr$restricted_bounds$horizon, expected_horizons)
  expect_equal(cumul$cumulative_bounds$horizon, expected_horizons)

  # Test 7: Cumulative bounds are NA for preperiods but not for post-periods
  pre_rows_cumul <- cumul$cumulative_bounds$horizon < 0
  expect_true(all(is.na(cumul$cumulative_bounds$lower[pre_rows_cumul])))
  expect_true(all(is.na(cumul$cumulative_bounds$upper[pre_rows_cumul])))
  post_rows_cumul <- cumul$cumulative_bounds$horizon > 0
  expect_false(any(is.na(cumul$cumulative_bounds$lower[post_rows_cumul])))
  expect_false(any(is.na(cumul$cumulative_bounds$upper[post_rows_cumul])))
  expect_false(any(is.na(cumul$cumulative_bounds$coef)))

  # Test 8: Pre and post periods treated as independent blocks
  var_block <- matrix(0, npre + npost, npre + npost)
  var_block[1:npre, 1:npre] <- Vpre_iid[1:npre, 1:npre]
  var_block[(npre + 1):(npre + npost), (npre + 1):(npre + npost)] <- var_iid[1:npost, 1:npost]
  estimates_block <- c(pre_mean0[1:npre], estimates_constant[1:npost])
  result_block <- plausible_bounds(estimates_block, var_block, preperiods = npre)
  expect_true(all(var_block[1:npre, (npre + 1):(npre + npost)] == 0))
  expect_true(all(var_block[(npre + 1):(npre + npost), 1:npre] == 0))
  post_rows_block <- result_block$restricted_bounds$horizon > 0
  expect_false(any(is.na(result_block$restricted_bounds$surrogate[post_rows_block])))

  # Test 9: Metadata correctly populated
  expect_equal(result$preperiods, 8)
  expect_equal(restr$metadata$preperiods, 8)
  expect_equal(cumul$metadata$preperiods, 8)

  # Test 10: Wpre calculation is mathematically correct
  expected_stat <- as.numeric(t(pre_mean0[1:npre]) %*% solve(Vpre_iid[1:npre, 1:npre]) %*% pre_mean0[1:npre])
  expected_pval <- 1 - pchisq(expected_stat, df = npre)
  expect_equal(unname(restr$Wpre["statistic"]), expected_stat, tolerance = 1e-6)
  expect_equal(unname(restr$Wpre["pvalue"]), expected_pval, tolerance = 1e-6)
})

# Section 3: Horizon Numbering Tests ----

# Section 3: NA Bounds Tests ----

# Section 4: Wald Test Coverage ----

test_that("zero preperiods functionality", {
  # Setup data once for all assertions
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  result <- plausible_bounds(estimates_constant[1:7], var_iid[1:7, 1:7], preperiods = 0)

  # Test 1: Wpre is NULL when preperiods = 0
  expect_null(result$wald_test$pre)
  expect_true(!is.null(result$wald_test$post))

  # Test 2: Zero preperiods metadata and horizon numbering
  expect_equal(result$preperiods, 0)
  expect_true(all(result$restricted_bounds$horizon > 0))
})

test_that("Wpre fails for reject scenario with pre-trends", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Test with reject scenario (designed to fail) - use 8 pre + 7 post
  pre_reject <- readRDS(test_path("fixtures", "preperiods_reject.rds"))
  Vpre_reject <- readRDS(test_path("fixtures", "Vpre_reject.rds"))
  estimates <- c(pre_reject, estimates_constant[1:7])
  var <- block_diag_matrix(Vpre_reject, var_iid[1:7, 1:7])

  result <- plausible_bounds(estimates, var, preperiods = 8)

  expect_true(!is.null(result$wald_test$pre))
  expect_true(result$wald_test$pre$p_value < 0.05)  # Should reject H0
})

test_that("Wpost is computed on post-periods only", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # With preperiods - use 8 pre + 7 post
  estimates <- c(pre_mean0, estimates_constant[1:7])
  var <- block_diag_matrix(Vpre_iid, var_iid[1:7, 1:7])
  result_with_pre <- plausible_bounds(estimates, var, preperiods = 8)

  # Without preperiods (same post data) - use 6 post
  result_no_pre <- plausible_bounds(estimates_constant[1:7], var_iid[1:7, 1:7], preperiods = 0)

  # Wpost should be the same (only depends on post-periods)
  expect_equal(
    result_with_pre$wald_test$post$statistic,
    result_no_pre$wald_test$post$statistic,
    tolerance = 1e-6
  )
  expect_equal(
    result_with_pre$wald_test$post$p_value,
    result_no_pre$wald_test$post$p_value,
    tolerance = 1e-6
  )
})

# Section 5: Block Structure Tests ----

test_that("ATE is computed from post-periods only", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Use 8 pre + 7 post
  estimates <- c(pre_mean0, estimates_constant[1:7])
  var <- block_diag_matrix(Vpre_iid, var_iid[1:7, 1:7])

  result_with_pre <- plausible_bounds(estimates, var, preperiods = 8)
  result_no_pre <- plausible_bounds(estimates_constant[1:7], var_iid[1:7, 1:7], preperiods = 0)

  # ATE should be the same (only depends on post-periods)
  expect_equal(
    result_with_pre$avg_treatment_effect$estimate,
    result_no_pre$avg_treatment_effect$estimate,
    tolerance = 1e-10
  )
  expect_equal(
    result_with_pre$avg_treatment_effect$se,
    result_no_pre$avg_treatment_effect$se,
    tolerance = 1e-10
  )
})

# Section 6: Pointwise and Sup-t Bounds ----

test_that("supt bounds use critical value from all periods", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Use 8 pre + 7 post
  estimates <- c(pre_mean0, estimates_constant[1:7])
  var <- block_diag_matrix(Vpre_iid, var_iid[1:7, 1:7])

  result <- plausible_bounds(estimates, var, preperiods = 8, include_supt = TRUE)

  expect_true(!is.null(result$supt_bounds))
  expect_equal(length(result$supt_bounds$lower), length(estimates))
  expect_equal(length(result$supt_bounds$upper), length(estimates))

  # No NAs in supt bounds
  expect_false(any(is.na(result$supt_bounds$lower)))
  expect_false(any(is.na(result$supt_bounds$upper)))

  # Check that supt critical value is stored in metadata
  expect_true(!is.null(result$restricted_bounds_metadata$supt_critval))
  expect_true(is.numeric(result$restricted_bounds_metadata$supt_critval))
})

# Section 7: Plotting Integration ----

test_that("plotting handles preperiods correctly", {
  # Setup data once for all assertions
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant[1:7])
  var <- block_diag_matrix(Vpre_iid, var_iid[1:7, 1:7])
  result <- plausible_bounds(estimates, var, preperiods = 8)

  # Call plot once (used across multiple assertions below)
  plot <- create_plot(result)

  # Test 1: create_plot returns ggplot object
  expect_s3_class(plot, "ggplot")

  # Test 2: X-axis includes negative values (pre-periods)
  built <- ggplot2::ggplot_build(plot)
  x_range <- range(built$layout$panel_params[[1]]$x.range)
  expect_true(x_range[1] < 0)
  expect_true(x_range[2] > 0)

  # Test 3: Vertical line at event time 0
  has_vline <- any(sapply(plot$layers, function(l) {
    inherits(l$geom, "GeomVline")
  }))
  expect_true(has_vline)
})

# Section 8: Integration Tests ----

# Section 9: Edge Cases and Boundary Conditions ----

test_that("single preperiod works correctly", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Single preperiod + 6 post
  pre_single <- rnorm(1, mean = 0, sd = 0.3)
  Vpre_single <- matrix(0.1, 1, 1)
  estimates <- c(pre_single, estimates_constant[1:7])
  var <- block_diag_matrix(Vpre_single, var_iid[1:7, 1:7])

  result <- plausible_bounds(estimates, var, preperiods = 1)
  expect_equal(result$preperiods, 1)
  expect_true(!is.null(result$wald_test$pre))
  expect_equal(min(result$restricted_bounds$horizon), -1)
})
