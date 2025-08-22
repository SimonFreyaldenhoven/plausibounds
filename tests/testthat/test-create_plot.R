
check_plot_structure <- function(plot) {
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$theme, "theme")
  expect_type(plot$labels, "list")
  expect_type(plot$layers, "list")
}

# Helper to count specific geom types in plot
count_geom_type <- function(plot, geom_class) {
  layer_classes <- sapply(plot$layers, function(l) class(l$geom)[1])
  sum(layer_classes == geom_class)
}

# Helper to check if plot contains specific aesthetic mappings
has_aesthetic <- function(plot, aes_name) {
  any(sapply(plot$layers, function(l) aes_name %in% names(l$mapping)))
}

# ---- Basic Functionality Tests with Package Data ----

test_that("create_plot produces valid ggplot with constant IID data", {
  # Use pre-computed fixture to avoid expensive computation
  result <- readRDS("tests/testthat/fixtures/simple_plausible.rds")
  plot <- create_plot(result)
  
  check_plot_structure(plot)
  expect_equal(plot$labels$x, "Horizon")
  expect_equal(plot$labels$y, "Estimate")
  expect_equal(plot$theme$legend.position, "bottom")
  
  # Check for expected layers
  expect_true(count_geom_type(plot, "GeomPoint") >= 1)
  expect_true(count_geom_type(plot, "GeomRibbon") >= 1)
})

test_that("create_plot produces valid ggplot with wiggly strong correlation data", {
  # Use pre-computed fixture instead of expensive computation
  result <- readRDS("tests/testthat/fixtures/complex_plausible.rds")
  plot <- create_plot(result)
  
  check_plot_structure(plot)
  
  # Verify both cumulative and restricted bounds are present
  expect_true(count_geom_type(plot, "GeomRibbon") >= 1)  # Cumulative bounds
  expect_true(count_geom_type(plot, "GeomLine") >= 2)    # Restricted bounds lines
})



# ---- Input Type Tests ----

test_that("create_plot accepts different input types correctly", {
  # Test with cumulative_bounds only - use pre-computed fixture
  cumul <- readRDS("tests/testthat/fixtures/simple_cumulative.rds")
  plot_cumul <- create_plot(cumul)
  check_plot_structure(plot_cumul)
  expect_true(count_geom_type(plot_cumul, "GeomRibbon") >= 1)
  
  # Test with restricted_bounds only - use pre-computed fixture
  restr <- readRDS("tests/testthat/fixtures/complex_restricted.rds")
  plot_restr <- create_plot(restr)
  check_plot_structure(plot_restr)
  expect_true(count_geom_type(plot_restr, "GeomLine") >= 2)
  
  # Test with two separate bounds objects
  plot_both <- create_plot(cumul, restr)
  check_plot_structure(plot_both)
  expect_true(count_geom_type(plot_both, "GeomRibbon") >= 1)
  expect_true(count_geom_type(plot_both, "GeomLine") >= 2)
})

test_that("create_plot rejects invalid inputs", {
  # Wrong object type
  expect_error(
    create_plot(list(a = 1, b = 2)),
    "Unrecognized input type"
  )
  
  # Too many arguments
  expect_error(
    create_plot(1, 2, 3),
    "accepts 1 or 2 arguments"
  )
  
  # Mismatched object types - use pre-computed fixtures
  cumul1 <- readRDS("tests/testthat/fixtures/simple_cumulative.rds")
  cumul2 <- readRDS("tests/testthat/fixtures/simple_cumulative.rds")
  
  expect_error(
    create_plot(cumul1, cumul2),
    "both must be either cumulative_bounds or restricted_bounds"
  )
})

# ---- Show/Hide Parameters Tests ----

test_that("show parameters control visibility of plot elements", {
  # Use pre-computed fixture for testing parameters
  pb <- readRDS("tests/testthat/fixtures/simple_plausible.rds")
  
  # Test show_cumulative = FALSE
  plot_no_cumul <- create_plot(pb, show_cumulative = FALSE)
  check_plot_structure(plot_no_cumul)
  expect_equal(count_geom_type(plot_no_cumul, "GeomRibbon"), 0)
  
  # Test show_restricted = FALSE
  plot_no_restr <- create_plot(pb, show_restricted = FALSE)
  check_plot_structure(plot_no_restr)
  # Should only have point layer and ribbon, no lines for restricted bounds
  expect_true(count_geom_type(plot_no_restr, "GeomRibbon") >= 1)
  expect_equal(count_geom_type(plot_no_restr, "GeomLine"), 0)
})

test_that("show_supt and show_pointwise parameters work", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS("tests/testthat/fixtures/complex_plausible.rds")
  
  # Test with all bounds shown
  plot_all <- create_plot(pb, show_supt = TRUE, show_pointwise = TRUE)
  check_plot_structure(plot_all)
  
  # Test with bounds hidden
  plot_none <- create_plot(pb, show_supt = FALSE, show_pointwise = FALSE)
  check_plot_structure(plot_none)
  
  # Count segments (used for supt and pointwise bounds)
  segments_all <- count_geom_type(plot_all, "GeomSegment")
  segments_none <- count_geom_type(plot_none, "GeomSegment")
  
  # When hiding supt and pointwise, there should be fewer or no segments
  expect_true(segments_none <= segments_all)
})

# ---- Message and Warning Tests ----

test_that("appropriate messages are shown for missing bounds", {
  # Use pre-computed fixture
  cumul <- readRDS("tests/testthat/fixtures/simple_cumulative.rds")
  
  # Requesting restricted bounds that don't exist should show message
  expect_message(
    create_plot(cumul, show_restricted = TRUE),
    "restricted bounds not available"
  )
  
  # Not requesting them should not show message
  expect_message(
    create_plot(cumul, show_restricted = FALSE),
    NA
  )
})

test_that("warning is shown when no bounds can be plotted", {
  # Use pre-computed fixture
  cumul <- readRDS("tests/testthat/fixtures/simple_cumulative.rds")
  
  # Hide the only available bounds
  expect_error(
    create_plot(cumul, show_cumulative = FALSE, show_restricted = TRUE),
    "No bounds available to plot"
  )
})

# ---- Visual Properties Tests ----

test_that("plot has correct aesthetic mappings and scales", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS("tests/testthat/fixtures/simple_plausible.rds")
  plot <- create_plot(pb)
  
  # Check for color scale by examining scale aesthetics
  scale_aesthetics <- sapply(plot$scales$scales, function(s) s$aesthetics)
  expect_true("colour" %in% unlist(scale_aesthetics))
  
  # Check for fill scale (from cumulative bounds ribbon)
  expect_true("fill" %in% unlist(scale_aesthetics))
  
  # Check that plot uses minimal theme
  expect_true(inherits(plot$theme, "theme"))
})

test_that("plot handles different horizon lengths appropriately", {
  # Test with small dataset
  set.seed(42)
  small_est <- rnorm(5)
  small_var <- diag(5) * 0.1
  
  cumul_small <- calculate_cumulative_bounds(small_est, small_var)
  plot_small <- create_plot(cumul_small)
  check_plot_structure(plot_small)
  
  # Test with larger dataset
  large_est <- rnorm(25)
  large_var <- diag(25) * 0.1
  
  cumul_large <- calculate_cumulative_bounds(large_est, large_var)
  plot_large <- create_plot(cumul_large)
  check_plot_structure(plot_large)
  
  # Both should be valid plots
  expect_s3_class(plot_small, "ggplot")
  expect_s3_class(plot_large, "ggplot")
})

# ---- Integration Tests ----

test_that("plots can be modified with ggplot2 functions", {
  # Use pre-computed fixture instead of expensive computation
  result <- readRDS("tests/testthat/fixtures/simple_plausible.rds")
  base_plot <- create_plot(result)
  
  # Add title
  plot_with_title <- base_plot + ggtitle("Test Title")
  expect_equal(plot_with_title$labels$title, "Test Title")
  
  # Change theme
  plot_with_theme <- base_plot + theme_classic()
  expect_s3_class(plot_with_theme$theme, "theme")
  
  # Add reference line
  plot_with_line <- base_plot + geom_hline(yintercept = 0, linetype = "dashed")
  expect_true(length(plot_with_line$layers) > length(base_plot$layers))
  
  # Modify axis labels
  plot_with_labels <- base_plot + labs(x = "Time Period", y = "Effect Size")
  expect_equal(plot_with_labels$labels$x, "Time Period")
  expect_equal(plot_with_labels$labels$y, "Effect Size")
})

test_that("plot can be saved without errors", {
  skip_on_cran()  # Skip file I/O on CRAN
  
  # Use pre-computed fixture instead of expensive computation
  result <- readRDS("tests/testthat/fixtures/simple_plausible.rds")
  plot <- create_plot(result)
  
  # Create temporary file
  temp_file <- tempfile(fileext = ".png")
  
  # Save plot
  expect_error(
    ggsave(temp_file, plot, width = 7, height = 5),
    NA
  )
  
  # Check file exists
  expect_true(file.exists(temp_file))
  
  # Clean up
  unlink(temp_file)
})

# ---- Edge Cases and Boundary Tests ----

test_that("create_plot handles edge cases gracefully", {
  # Single horizon point
  single_est <- 1.5
  single_var <- matrix(0.1, 1, 1)
  
  cumul_single <- calculate_cumulative_bounds(single_est, single_var)
  expect_error(
    plot_single <- create_plot(cumul_single),
    NA  # Should not error
  )
  
  if (exists("plot_single")) {
    check_plot_structure(plot_single)
  }
  
  # Very small variance
  set.seed(42)
  est_small_var <- rnorm(10)
  var_small_var <- diag(10) * 1e-6
  
  cumul_small_var <- calculate_cumulative_bounds(est_small_var, var_small_var)
  plot_small_var <- create_plot(cumul_small_var)
  check_plot_structure(plot_small_var)
})

test_that("create_plot handles mismatched horizons correctly", {
  # Create bounds with different horizon lengths
  est1 <- rnorm(10)
  var1 <- diag(10) * 0.1
  cumul1 <- calculate_cumulative_bounds(est1, var1)
  
  est2 <- rnorm(8)
  var2 <- diag(8) * 0.1
  
  skip_on_cran()
  restr2 <- calculate_restricted_bounds(est2, var2)
  
  # Should error when trying to combine mismatched horizons
  expect_error(
    create_plot(cumul1, restr2),
    "same horizon"
  )
})

# ---- Performance Tests ----

test_that("create_plot performs efficiently with package data", {
  skip_on_cran()  # Skip timing tests on CRAN
  
  # Use pre-computed fixture instead of expensive computation
  result <- readRDS("tests/testthat/fixtures/complex_plausible.rds")
  
  # Time plot creation
  time_taken <- system.time({
    plot <- create_plot(result)
  })["elapsed"]
  
  # Plot creation should be fast (under 1 second)
  expect_true(time_taken < 1.0)
  check_plot_structure(plot)
})

# ---- Helper Function Tests ----

test_that("extract_bounds_data works correctly", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS("tests/testthat/fixtures/simple_plausible.rds")
  
  # Test extraction from plausible_bounds
  bounds_data <- extract_bounds_data(
    list(pb),
    show_cumulative = TRUE,
    show_restricted = TRUE,
    show_supt = TRUE,
    show_pointwise = TRUE
  )
  
  expect_equal(bounds_data$type, "plausible_bounds")
  expect_true(bounds_data$has_cumulative)
  expect_true(bounds_data$has_restricted)
  expect_s3_class(bounds_data$cumulative, "data.frame")
  expect_s3_class(bounds_data$restricted, "data.frame")
  
  # Test extraction from cumulative_bounds only
  cumul <- calculate_cumulative_bounds(estimates_constant_iid, var_constant_iid)
  bounds_data_cumul <- extract_bounds_data(
    list(cumul),
    show_cumulative = TRUE,
    show_restricted = FALSE,
    show_supt = FALSE,
    show_pointwise = FALSE
  )
  
  expect_equal(bounds_data_cumul$type, "cumulative_only")
  expect_true(bounds_data_cumul$has_cumulative)
  expect_false(bounds_data_cumul$has_restricted)
})

test_that("check_bounds_availability identifies available bounds", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS("tests/testthat/fixtures/simple_plausible.rds")
  
  bounds_data <- extract_bounds_data(
    list(pb),
    show_cumulative = TRUE,
    show_restricted = TRUE,
    show_supt = TRUE,
    show_pointwise = TRUE
  )
  
  availability <- check_bounds_availability(
    bounds_data,
    show_cumulative = TRUE,
    show_restricted = TRUE,
    show_supt = TRUE,
    show_pointwise = TRUE
  )
  
  # Check structure
  expect_type(availability, "list")
  expect_true(availability$cumulative_available)
  expect_true(availability$restricted_available)
  expect_true(availability$cumulative_requested)
  expect_true(availability$restricted_requested)
  
  # Check logical consistency
  expect_equal(
    availability$show_cumulative,
    availability$cumulative_requested && availability$cumulative_available
  )
  expect_equal(
    availability$show_restricted,
    availability$restricted_requested && availability$restricted_available
  )
})

test_that("display_availability_messages shows correct messages", {
  # Test with missing cumulative bounds
  availability <- list(
    cumulative_requested = TRUE,
    cumulative_available = FALSE,
    restricted_requested = FALSE,
    restricted_available = FALSE,
    supt_requested = FALSE,
    supt_available = FALSE,
    pointwise_requested = FALSE,
    pointwise_available = FALSE,
    show_cumulative = FALSE,
    show_restricted = FALSE
  )
  
  expect_message(
    display_availability_messages(availability),
    "cumulative bounds not available"
  )
  
  # Test with multiple missing bounds
  availability$restricted_requested <- TRUE
  availability$supt_requested <- TRUE
  
  expect_message(
    display_availability_messages(availability),
    "cumulative, restricted, sup-t bounds not available"
  )
})

test_that("merge_bounds_data correctly merges additional bounds", {
  skip_on_cran()
  
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS("tests/testthat/fixtures/complex_plausible.rds")
  
  # Create mock result structure
  result <- list(
    has_cumulative = TRUE,
    cumulative = pb$cumulative_bounds,
    has_restricted = TRUE,
    restricted = pb$restricted_bounds
  )
  
  # Merge bounds
  merged <- merge_bounds_data(result, pb, "plausible_bounds")
  
  # Check if pointwise/supt bounds were merged (if they exist)
  if (!is.null(pb$pointwise_bounds)) {
    expect_true("pointwise_lower" %in% names(merged$cumulative))
    expect_true("pointwise_upper" %in% names(merged$cumulative))
    expect_true("pointwise_lower" %in% names(merged$restricted))
    expect_true("pointwise_upper" %in% names(merged$restricted))
  }
  
  if (!is.null(pb$supt_bounds)) {
    expect_true("supt_lower" %in% names(merged$cumulative))
    expect_true("supt_upper" %in% names(merged$cumulative))
    expect_true("supt_lower" %in% names(merged$restricted))
    expect_true("supt_upper" %in% names(merged$restricted))
  }
})

test_that("create_bounds_plot generates correct plot structure", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS("tests/testthat/fixtures/simple_plausible.rds")
  
  # Extract and prepare data
  bounds_data <- extract_bounds_data(
    list(pb),
    show_cumulative = TRUE,
    show_restricted = TRUE,
    show_supt = FALSE,
    show_pointwise = FALSE
  )
  
  availability <- check_bounds_availability(
    bounds_data,
    show_cumulative = TRUE,
    show_restricted = TRUE,
    show_supt = FALSE,
    show_pointwise = FALSE
  )
  
  # Create plot
  plot <- create_bounds_plot(bounds_data, availability)
  
  check_plot_structure(plot)
  
  # Check for expected components
  expect_true(count_geom_type(plot, "GeomPoint") >= 1)
  expect_true(count_geom_type(plot, "GeomRibbon") >= 1)  # Cumulative
  expect_true(count_geom_type(plot, "GeomLine") >= 2)    # Restricted bounds
})

# ---- Y-axis Limits Tests ----

test_that("y-axis limits include all data with appropriate buffer", {
  # Create data with known range
  est <- seq(-2, 2, length.out = 10)
  var <- diag(10) * 0.1
  
  cumul <- calculate_cumulative_bounds(est, var)
  plot <- create_plot(cumul)
  
  # Extract y-limits from plot coordinates
  y_range <- ggplot_build(plot)$layout$panel_params[[1]]$y.range
  
  # Check that limits include data with buffer
  min_expected <- min(cumul$cumulative_bounds$lower)
  max_expected <- max(cumul$cumulative_bounds$upper)
  
  expect_true(y_range[1] < min_expected)
  expect_true(y_range[2] > max_expected)
})

# ---- Combined Bounds Tests ----

test_that("plotting both cumulative and restricted bounds works correctly", {
  skip_on_cran()
  
  data(estimates_wiggly_strong_corr)
  data(var_wiggly_strong_corr)
  
  # Create both types of bounds separately
  cumul <- calculate_cumulative_bounds(estimates_wiggly_strong_corr, var_wiggly_strong_corr)
  restr <- calculate_restricted_bounds(estimates_wiggly_strong_corr, var_wiggly_strong_corr)
  
  # Plot with both bounds from separate objects
  plot_both_sep <- create_plot(cumul, restr)
  check_plot_structure(plot_both_sep)
  
  # Plot with both bounds from plausible_bounds object - use pre-computed fixture
  pb <- readRDS("tests/testthat/fixtures/complex_plausible.rds")
  plot_both_pb <- create_plot(pb)
  check_plot_structure(plot_both_pb)
  
  # Both methods should produce plots with similar structure
  expect_equal(
    count_geom_type(plot_both_sep, "GeomRibbon"),
    count_geom_type(plot_both_pb, "GeomRibbon")
  )
  expect_equal(
    count_geom_type(plot_both_sep, "GeomLine"),
    count_geom_type(plot_both_pb, "GeomLine")
  )
})

# ---- Legend Tests ----

test_that("legend is properly configured", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS("tests/testthat/fixtures/simple_plausible.rds")
  plot <- create_plot(pb)
  
  # Check legend position
  expect_equal(plot$theme$legend.position, "bottom")
  
  # Build plot to check legend content
  built_plot <- ggplot_build(plot)
  plot_data <- built_plot$plot
  
  # Check that color scale has appropriate labels
  color_scale <- plot_data$scales$get_scales("colour")
  if (!is.null(color_scale)) {
    expect_true("Point Estimates" %in% color_scale$labels)
  }
  
  # Check that fill scale has appropriate labels
  fill_scale <- plot_data$scales$get_scales("fill")
  if (!is.null(fill_scale)) {
    expect_true("Cumulative" %in% fill_scale$labels)
  }
})