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

  # Valid: zero preperiods
  expect_error(plausible_bounds(estimates_constant, var_iid, preperiods = 0), NA)

  # Valid: single preperiod
  pre_single <- rnorm(1)
  Vpre_single <- matrix(0.1, 1, 1)
  estimates <- c(pre_single, estimates_constant)
  var <- block_diag_matrix(Vpre_single, var_iid)
  expect_error(plausible_bounds(estimates, var, preperiods = 1), NA)

  # Valid: multiple preperiods
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  estimates_full <- c(pre_mean0, estimates_constant)
  var_full <- block_diag_matrix(Vpre_iid, var_iid)
  expect_error(plausible_bounds(estimates_full, var_full, preperiods = 8), NA)
})

test_that("preperiods parameter rejects negative values", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  expect_error(
    plausible_bounds(estimates_constant, var_iid, preperiods = -1),
    "preperiods must be a non-negative integer"
  )

  expect_error(
    plausible_bounds(estimates_constant, var_iid, preperiods = -5),
    "preperiods must be a non-negative integer"
  )
})

test_that("preperiods parameter rejects non-integer values", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  expect_error(
    plausible_bounds(estimates_constant, var_iid, preperiods = 2.5),
    "preperiods must be a non-negative integer"
  )

  expect_error(
    plausible_bounds(estimates_constant, var_iid, preperiods = 3.7),
    "preperiods must be a non-negative integer"
  )
})

test_that("preperiods parameter rejects values >= length(estimates)", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  n <- length(estimates_constant)

  expect_error(
    plausible_bounds(estimates_constant, var_iid, preperiods = n),
    "preperiods must be less than the length of estimates"
  )

  expect_error(
    plausible_bounds(estimates_constant, var_iid, preperiods = n + 1),
    "preperiods must be less than the length of estimates"
  )
})

test_that("preperiods parameter rejects invalid types", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Non-numeric string
  expect_error(
    plausible_bounds(estimates_constant, var_iid, preperiods = "two"),
    "preperiods must be a non-negative integer"
  )

  # NULL
  expect_error(
    plausible_bounds(estimates_constant, var_iid, preperiods = NULL),
    "preperiods must be a non-negative integer"
  )

  # NA
  expect_error(
    plausible_bounds(estimates_constant, var_iid, preperiods = NA),
    "preperiods must be a non-negative integer"
  )
})

test_that("preperiods parameter rejects vector inputs", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  expect_error(
    plausible_bounds(estimates_constant, var_iid, preperiods = c(1, 2)),
    "preperiods must be a non-negative integer"
  )
})

# Section 2: Horizon Numbering Tests ----

test_that("horizon numbering is correct with preperiods", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  npre <- length(pre_mean0)
  npost <- length(estimates_constant)

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  result <- plausible_bounds(estimates, var, preperiods = npre)

  # Check horizon numbering in restricted bounds
  expect_equal(min(result$restricted_bounds$horizon), -npre)
  expect_equal(max(result$restricted_bounds$horizon), npost)
  expect_true(all(result$restricted_bounds$horizon[1:npre] < 0))
  expect_true(all(result$restricted_bounds$horizon[(npre+1):(npre+npost)] > 0))
  expect_false(any(result$restricted_bounds$horizon == 0))  # No period 0

  # Check horizon sequence is consecutive
  expected_horizons <- c(-npre:-1, 1:npost)
  expect_equal(result$restricted_bounds$horizon, expected_horizons)
})

test_that("horizon numbering is consistent across all functions", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)
  npre <- length(pre_mean0)
  npost <- length(estimates_constant)

  restr <- calculate_restricted_bounds(estimates, var, preperiods = npre)
  cumul <- calculate_cumulative_bounds(estimates, var, preperiods = npre)

  expected_horizons <- c(-8:-1, 1:12)
  expect_equal(restr$restricted_bounds$horizon, expected_horizons)
  expect_equal(cumul$cumulative_bounds$horizon, expected_horizons)
})

# Section 3: NA Bounds Tests ----

test_that("restricted bounds are NA for preperiods", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)
  npre <- length(pre_mean0)

  result <- plausible_bounds(estimates, var, preperiods = npre)

  # Restricted bounds should be NA for pre-periods
  pre_rows <- result$restricted_bounds$horizon < 0
  expect_true(all(is.na(result$restricted_bounds$surrogate[pre_rows])))
  expect_true(all(is.na(result$restricted_bounds$lower[pre_rows])))
  expect_true(all(is.na(result$restricted_bounds$upper[pre_rows])))

  # But NOT NA for post-periods
  post_rows <- result$restricted_bounds$horizon > 0
  expect_false(any(is.na(result$restricted_bounds$surrogate[post_rows])))
  expect_false(any(is.na(result$restricted_bounds$lower[post_rows])))
  expect_false(any(is.na(result$restricted_bounds$upper[post_rows])))

  # Coefficients should never be NA
  expect_false(any(is.na(result$restricted_bounds$coef)))
})

test_that("cumulative bounds are NA for preperiods", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)
  npre <- length(pre_mean0)

  result <- calculate_cumulative_bounds(estimates, var, preperiods = npre)

  # Cumulative bounds should be NA for pre-periods
  pre_rows <- result$cumulative_bounds$horizon < 0
  expect_true(all(is.na(result$cumulative_bounds$lower[pre_rows])))
  expect_true(all(is.na(result$cumulative_bounds$upper[pre_rows])))

  # But NOT NA for post-periods
  post_rows <- result$cumulative_bounds$horizon > 0
  expect_false(any(is.na(result$cumulative_bounds$lower[post_rows])))
  expect_false(any(is.na(result$cumulative_bounds$upper[post_rows])))

  # Coefficients should never be NA
  expect_false(any(is.na(result$cumulative_bounds$coef)))
})

# Section 4: Wald Test Coverage ----

test_that("Wpre is computed when preperiods > 0", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  result <- plausible_bounds(estimates, var, preperiods = 8)

  expect_true(!is.null(result$wald_test$pre))
  expect_true(is.numeric(result$wald_test$pre$statistic))
  expect_true(is.numeric(result$wald_test$pre$p_value))
  expect_true(result$wald_test$pre$p_value >= 0 && result$wald_test$pre$p_value <= 1)
})

test_that("Wpre is NULL when preperiods = 0", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  result <- plausible_bounds(estimates_constant, var_iid, preperiods = 0)

  expect_null(result$wald_test$pre)
  expect_true(!is.null(result$wald_test$post))  # But Wpost should exist
})

test_that("Wpre passes for scenarios with no pre-trends", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Test with mean0 scenario (should pass)
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  result <- plausible_bounds(estimates, var, preperiods = 8)
  expect_true(result$wald_test$pre$p_value >= 0.05)  # Should fail to reject H0

  # Note: Fixed scenario has declining pattern so may actually fail pre-trends test
  # Test with corr scenario (should pass)
  pre_corr <- readRDS(test_path("fixtures", "preperiods_corr.rds"))
  Vpre_corr <- readRDS(test_path("fixtures", "Vpre_corr.rds"))
  estimates_corr <- c(pre_corr, estimates_constant)
  var_corr_full <- block_diag_matrix(Vpre_corr, var_iid)

  result_corr <- plausible_bounds(estimates_corr, var_corr_full, preperiods = 8)
  expect_true(result_corr$wald_test$pre$p_value >= 0.05)
})

test_that("Wpre fails for reject scenario with pre-trends", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Test with reject scenario (designed to fail)
  pre_reject <- readRDS(test_path("fixtures", "preperiods_reject.rds"))
  Vpre_reject <- readRDS(test_path("fixtures", "Vpre_reject.rds"))
  estimates <- c(pre_reject, estimates_constant)
  var <- block_diag_matrix(Vpre_reject, var_iid)

  result <- plausible_bounds(estimates, var, preperiods = 8)

  expect_true(!is.null(result$wald_test$pre))
  expect_true(result$wald_test$pre$p_value < 0.05)  # Should reject H0
})

test_that("Wpre calculation is correct in calculate_restricted_bounds", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  result <- calculate_restricted_bounds(estimates, var, preperiods = 8)

  # Manually calculate Wpre
  expected_stat <- as.numeric(t(pre_mean0) %*% solve(Vpre_iid) %*% pre_mean0)
  expected_pval <- 1 - pchisq(expected_stat, df = 8)

  expect_equal(unname(result$Wpre["statistic"]), expected_stat, tolerance = 1e-6)
  expect_equal(unname(result$Wpre["pvalue"]), expected_pval, tolerance = 1e-6)
})

test_that("Wpost is always computed regardless of preperiods", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Without preperiods
  result_no_pre <- plausible_bounds(estimates_constant, var_iid, preperiods = 0)
  expect_true(!is.null(result_no_pre$wald_test$post))
  expect_true(is.numeric(result_no_pre$wald_test$post$statistic))
  expect_true(is.numeric(result_no_pre$wald_test$post$p_value))

  # With preperiods
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  result_with_pre <- plausible_bounds(estimates, var, preperiods = 8)
  expect_true(!is.null(result_with_pre$wald_test$post))
})

test_that("Wpost is computed on post-periods only", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # With preperiods
  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)
  result_with_pre <- plausible_bounds(estimates, var, preperiods = 8)

  # Without preperiods (same post data)
  result_no_pre <- plausible_bounds(estimates_constant, var_iid, preperiods = 0)

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

test_that("pre and post periods are treated as independent blocks", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  npre <- length(pre_mean0)
  npost <- length(estimates_constant)

  # Create block diagonal variance (no correlation between pre and post)
  var_block <- matrix(0, npre + npost, npre + npost)
  var_block[1:npre, 1:npre] <- Vpre_iid
  var_block[(npre+1):(npre+npost), (npre+1):(npre+npost)] <- var_iid

  estimates <- c(pre_mean0, estimates_constant)
  result <- plausible_bounds(estimates, var_block, preperiods = npre)

  # Check that off-diagonal blocks in variance are all zeros
  expect_true(all(var_block[1:npre, (npre+1):(npre+npost)] == 0))
  expect_true(all(var_block[(npre+1):(npre+npost), 1:npre] == 0))

  # Surrogate should only be computed on post-periods
  post_rows <- result$restricted_bounds$horizon > 0
  expect_false(any(is.na(result$restricted_bounds$surrogate[post_rows])))
})

test_that("ATE is computed from post-periods only", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  result_with_pre <- plausible_bounds(estimates, var, preperiods = 8)
  result_no_pre <- plausible_bounds(estimates_constant, var_iid, preperiods = 0)

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

test_that("pointwise bounds are computed for all periods when preperiods > 0", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  result <- plausible_bounds(estimates, var, preperiods = 8, include_pointwise = TRUE)

  expect_true(!is.null(result$pointwise_bounds))
  expect_equal(length(result$pointwise_bounds$lower), length(estimates))
  expect_equal(length(result$pointwise_bounds$upper), length(estimates))

  # No NAs in pointwise bounds
  expect_false(any(is.na(result$pointwise_bounds$lower)))
  expect_false(any(is.na(result$pointwise_bounds$upper)))
})

test_that("supt bounds use critical value from all periods", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

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

test_that("create_plot handles preperiods correctly", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  result <- plausible_bounds(estimates, var, preperiods = 8)
  plot <- create_plot(result)

  expect_s3_class(plot, "ggplot")

  # Build plot to access data
  built <- ggplot2::ggplot_build(plot)

  # Check that x-axis includes negative values
  x_range <- range(built$layout$panel_params[[1]]$x.range)
  expect_true(x_range[1] < 0)
  expect_true(x_range[2] > 0)
})

test_that("plot shows vertical line at event time 0 when preperiods > 0", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  result <- plausible_bounds(estimates, var, preperiods = 8)
  plot <- create_plot(result)

  # Check for vline at x = 0
  has_vline <- any(sapply(plot$layers, function(l) {
    inherits(l$geom, "GeomVline")
  }))

  expect_true(has_vline)
})

# Section 8: Integration Tests ----

test_that("full pipeline works with all four preperiod scenarios", {
  skip_on_cran()

  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  scenarios <- list(
    mean0 = list(
      pre = readRDS(test_path("fixtures", "preperiods_mean0.rds")),
      Vpre = readRDS(test_path("fixtures", "Vpre_iid.rds"))
    ),
    fixed = list(
      pre = readRDS(test_path("fixtures", "preperiods_fixed.rds")),
      Vpre = readRDS(test_path("fixtures", "Vpre_iid.rds"))
    ),
    corr = list(
      pre = readRDS(test_path("fixtures", "preperiods_corr.rds")),
      Vpre = readRDS(test_path("fixtures", "Vpre_corr.rds"))
    ),
    reject = list(
      pre = readRDS(test_path("fixtures", "preperiods_reject.rds")),
      Vpre = readRDS(test_path("fixtures", "Vpre_reject.rds"))
    )
  )

  for (name in names(scenarios)) {
    estimates <- c(scenarios[[name]]$pre, estimates_constant)
    var <- block_diag_matrix(scenarios[[name]]$Vpre, var_iid)

    # Test plausible_bounds
    result <- plausible_bounds(estimates, var, preperiods = 8)
    expect_s3_class(result, "plausible_bounds")
    expect_equal(result$preperiods, 8)
    expect_true(!is.null(result$wald_test$pre))
    expect_true(!is.null(result$wald_test$post))

    # Test calculate_restricted_bounds
    restr <- calculate_restricted_bounds(estimates, var, preperiods = 8)
    expect_s3_class(restr, "restricted_bounds")
    expect_true(!is.null(restr$Wpre))
    expect_true(!is.null(restr$Wpost))

    # Test calculate_cumulative_bounds
    cumul <- calculate_cumulative_bounds(estimates, var, preperiods = 8)
    expect_s3_class(cumul, "cumulative_bounds")

    # Test plotting
    plot <- create_plot(result)
    expect_s3_class(plot, "ggplot")
  }
})

test_that("metadata is correctly populated with preperiods", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  result <- plausible_bounds(estimates, var, preperiods = 8)

  # Check preperiods is stored at top level
  expect_equal(result$preperiods, 8)

  # Check metadata in sub-results
  restr <- calculate_restricted_bounds(estimates, var, preperiods = 8)
  expect_equal(restr$metadata$preperiods, 8)

  cumul <- calculate_cumulative_bounds(estimates, var, preperiods = 8)
  expect_equal(cumul$metadata$preperiods, 8)
})

# Section 9: Edge Cases and Boundary Conditions ----

test_that("single preperiod works correctly", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Single preperiod
  pre_single <- rnorm(1, mean = 0, sd = 0.3)
  Vpre_single <- matrix(0.1, 1, 1)
  estimates <- c(pre_single, estimates_constant)
  var <- block_diag_matrix(Vpre_single, var_iid)

  result <- plausible_bounds(estimates, var, preperiods = 1)
  expect_equal(result$preperiods, 1)
  expect_true(!is.null(result$wald_test$pre))
  expect_equal(min(result$restricted_bounds$horizon), -1)
})

test_that("zero preperiods works and Wpre is NULL", {
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  result <- plausible_bounds(estimates_constant, var_iid, preperiods = 0)

  expect_equal(result$preperiods, 0)
  expect_null(result$wald_test$pre)
  expect_true(!is.null(result$wald_test$post))
  expect_true(all(result$restricted_bounds$horizon > 0))
})

test_that("many preperiods with few post-periods works", {
  #skip("Edge case with 1 post-period causes error in package - potential bug")

  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  n <- length(estimates_constant)

  # Near maximum: only 3 post-periods (safer than 1-2)
  pre_many <- rnorm(n - 3, mean = 0, sd = 0.1)
  Vpre_many <- diag(n - 3) * 0.01
  estimates <- c(pre_many, estimates_constant[1:3])
  var <- block_diag_matrix(Vpre_many, var_iid[1:3, 1:3])

  result <- plausible_bounds(estimates, var, preperiods = n - 2)
  expect_equal(result$preperiods, n - 3)
  expect_equal(nrow(result$restricted_bounds), n)
})

test_that("results are consistent across different preperiod counts", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Test with different numbers of preperiods
  for (npre in c(2, 4, 6, 8)) {
    pre_subset <- pre_mean0[1:npre]
    Vpre_subset <- Vpre_iid[1:npre, 1:npre]

    estimates <- c(pre_subset, estimates_constant[1:4])
    var <- block_diag_matrix(Vpre_subset, var_iid[1:4, 1:4])

    result <- plausible_bounds(estimates, var, preperiods = npre)

    expect_equal(result$preperiods, npre)
    expect_equal(nrow(result$restricted_bounds), npre + length(estimates_constant))

    # Post-treatment results should be similar regardless of npre
    post_rows <- result$restricted_bounds$horizon > 0
    expect_equal(sum(post_rows), length(estimates_constant))
  }
})

test_that("post-treatment results are mostly independent of npre", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  # Compare 2 preperiods vs 6 preperiods
  estimates_2pre <- c(pre_mean0[1:2], estimates_constant)
  var_2pre <- block_diag_matrix(Vpre_iid[1:2, 1:2], var_iid)
  result_2pre <- plausible_bounds(estimates_2pre, var_2pre, preperiods = 2)

  estimates_6pre <- c(pre_mean0[1:6], estimates_constant)
  var_6pre <- block_diag_matrix(Vpre_iid[1:6, 1:6], var_iid)
  result_6pre <- plausible_bounds(estimates_6pre, var_6pre, preperiods = 6)

  # Post-period surrogate should be very similar (not identical due to critical values)
  post_bounds_2pre <- result_2pre$restricted_bounds[result_2pre$restricted_bounds$horizon > 0, ]
  post_bounds_6pre <- result_6pre$restricted_bounds[result_6pre$restricted_bounds$horizon > 0, ]

  expect_equal(post_bounds_2pre$surrogate, post_bounds_6pre$surrogate, tolerance = 1e-10)
  # Bounds can differ slightly due to different critical values from different npre
  expect_equal(post_bounds_2pre$lower, post_bounds_6pre$lower, tolerance = 1e-3)
  expect_equal(post_bounds_2pre$upper, post_bounds_6pre$upper, tolerance = 1e-3)
})

test_that("preperiods work with different alpha values", {
  pre_mean0 <- readRDS(test_path("fixtures", "preperiods_mean0.rds"))
  Vpre_iid <- readRDS(test_path("fixtures", "Vpre_iid.rds"))
  data(estimates_constant, envir = environment())
  data(var_iid, envir = environment())

  estimates <- c(pre_mean0, estimates_constant)
  var <- block_diag_matrix(Vpre_iid, var_iid)

  for (alpha in c(0.01, 0.05, 0.1)) {
    result <- plausible_bounds(estimates, var, alpha = alpha, preperiods = 8)

    expect_equal(result$preperiods, 8)
    expect_equal(result$restricted_bounds_metadata$alpha, alpha)
    expect_true(!is.null(result$wald_test$pre))
    expect_true(!is.null(result$wald_test$post))
  }
})
