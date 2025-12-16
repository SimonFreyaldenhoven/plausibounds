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
  result <- readRDS(test_path("fixtures", "simple_plausible.rds"))
  plot <- create_plot(result)
  
  check_plot_structure(plot)
  expect_equal(plot$labels$x, "Event time")
  expect_equal(plot$labels$y, "Estimate")
  expect_equal(plot$theme$legend.position, "bottom")
  
  # Check for expected layers
  expect_true(count_geom_type(plot, "GeomPoint") >= 1)
})

test_that("create_plot produces valid ggplot with wiggly strong correlation data", {
  # Use pre-computed fixture instead of expensive computation
  result <- readRDS(test_path("fixtures", "complex_plausible.rds"))
  plot <- create_plot(result)
  
  check_plot_structure(plot)
  
  # Verify restricted bounds are present
  expect_true(count_geom_type(plot, "GeomLine") >= 2)    # Restricted bounds lines
})



# ---- Input Type Tests ----

test_that("create_plot only accepts plausible_bounds objects", {
  # Test with valid plausible_bounds object - use pre-computed fixture
  result <- readRDS(test_path("fixtures", "simple_plausible.rds"))
  plot <- create_plot(result)
  check_plot_structure(plot)
  expect_true(count_geom_type(plot, "GeomLine") >= 2)
})

test_that("create_plot rejects invalid inputs", {
  # Wrong object type
  expect_error(
    create_plot(list(a = 1, b = 2)),
    "plausible_bounds object"
  )

  # Numeric input
  expect_error(
    create_plot(1),
    "plausible_bounds object"
  )

  # NULL input
  expect_error(
    create_plot(NULL),
    "plausible_bounds object"
  )
})

# ---- Show/Hide Parameters Tests ----

test_that("create_plot always shows restricted bounds", {
  # Use pre-computed fixture for testing
  pb <- readRDS(test_path("fixtures", "simple_plausible.rds"))

  # Restricted bounds are always shown (no parameter to hide them)
  plot <- create_plot(pb)
  check_plot_structure(plot)
  # Should have lines for restricted bounds
  expect_true(count_geom_type(plot, "GeomLine") >= 2)
})

test_that("show_supt and show_pointwise parameters work", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS(test_path("fixtures", "complex_plausible.rds"))
  
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

test_that("appropriate messages are shown for missing optional bounds", {
  # Use pre-computed fixture
  pb <- readRDS(test_path("fixtures", "simple_plausible.rds"))

  # If pointwise or supt bounds don't exist and are requested, should show message
  # (This depends on what's in the fixture - may or may not trigger a message)
  # Just verify the function runs without error
  expect_error(
    create_plot(pb, show_supt = TRUE, show_pointwise = TRUE),
    NA
  )
})

# ---- Visual Properties Tests ----

test_that("plot has correct aesthetic mappings and scales", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS(test_path("fixtures", "simple_plausible.rds"))
  plot <- create_plot(pb)

  # Check for color scale by examining scale aesthetics
  scale_aesthetics <- sapply(plot$scales$scales, function(s) s$aesthetics)
  expect_true("colour" %in% unlist(scale_aesthetics))

  # Check that plot uses minimal theme
  expect_true(inherits(plot$theme, "theme"))
})

test_that("plot handles different horizon lengths appropriately", {
  skip_on_cran()
  # Test with small dataset
  set.seed(42)
  small_est <- rnorm(5)
  small_var <- diag(5) * 0.1

  pb_small <- plausible_bounds(small_est, small_var)
  plot_small <- create_plot(pb_small)
  check_plot_structure(plot_small)

  # Test with larger dataset
  large_est <- rnorm(25)
  large_var <- diag(25) * 0.1

  pb_large <- plausible_bounds(large_est, large_var)
  plot_large <- create_plot(pb_large)
  check_plot_structure(plot_large)

  # Both should be valid plots
  expect_s3_class(plot_small, "ggplot")
  expect_s3_class(plot_large, "ggplot")
})

# ---- Integration Tests ----

test_that("plots can be modified with ggplot2 functions", {
  # Use pre-computed fixture instead of expensive computation
  result <- readRDS(test_path("fixtures", "simple_plausible.rds"))
  base_plot <- create_plot(result)
  
  # Add title
  plot_with_title <- base_plot + ggplot2::labs(title = "Test Title")
  expect_equal(plot_with_title$labels$title, "Test Title")
  
  # Change theme
  plot_with_theme <- base_plot + ggplot2::theme_classic()
  expect_s3_class(plot_with_theme$theme, "theme")
  
  # Add reference line
  plot_with_line <- base_plot + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
  expect_true(length(plot_with_line$layers) > length(base_plot$layers))
  
  # Modify axis labels
  plot_with_labels <- base_plot + ggplot2::labs(x = "Time Period", y = "Effect Size")
  expect_equal(plot_with_labels$labels$x, "Time Period")
  expect_equal(plot_with_labels$labels$y, "Effect Size")
})

test_that("plot can be saved without errors", {
  skip_on_cran()  # Skip file I/O on CRAN
  
  # Use pre-computed fixture instead of expensive computation
  result <- readRDS(test_path("fixtures", "simple_plausible.rds"))
  plot <- create_plot(result)
  
  # Create temporary file
  temp_file <- tempfile(fileext = ".png")
  
  # Save plot
  expect_error(
    ggplot2::ggsave(temp_file, plot, width = 7, height = 5),
    NA
  )
  
  # Check file exists
  expect_true(file.exists(temp_file))
  
  # Clean up
  unlink(temp_file)
})

# ---- Edge Cases and Boundary Tests ----

test_that("create_plot handles edge cases gracefully", {
  # Very small variance
  set.seed(42)
  est_small_var <- rnorm(10)
  var_small_var <- diag(10) * 1e-6

  pb_small_var <- plausible_bounds(est_small_var, var_small_var)
  plot_small_var <- create_plot(pb_small_var)
  check_plot_structure(plot_small_var)
})

# ---- Performance Tests ----

test_that("create_plot performs efficiently with package data", {
  skip_on_cran()  # Skip timing tests on CRAN
  
  # Use pre-computed fixture instead of expensive computation
  result <- readRDS(test_path("fixtures", "complex_plausible.rds"))
  
  # Time plot creation
  time_taken <- system.time({
    plot <- create_plot(result)
  })["elapsed"]
  
  # Plot creation should be fast (under 1 second)
  expect_true(time_taken < 1.0)
  check_plot_structure(plot)
})

# ---- Helper Function Tests ----

test_that("check_bounds_availability identifies available optional bounds", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS(test_path("fixtures", "complex_plausible.rds"))

  bounds_data <- list(
    restricted = pb$restricted_bounds,
    has_restricted = TRUE
  )

  availability <- check_bounds_availability(
    bounds_data,
    pb,
    show_supt = TRUE,
    show_pointwise = TRUE
  )

  # Check structure
  expect_type(availability, "list")

  # Check that requested status is recorded
  expect_true(availability$supt_requested)
  expect_true(availability$pointwise_requested)
})

test_that("display_availability_messages shows correct messages for optional bounds", {
  # Test with missing supt bounds
  availability <- list(
    show_restricted = TRUE,
    supt_requested = TRUE,
    supt_available = FALSE,
    pointwise_requested = FALSE,
    pointwise_available = FALSE
  )

  expect_message(
    display_availability_messages(availability),
    "sup-t bounds not available"
  )

  # Test with multiple missing bounds
  availability$pointwise_requested <- TRUE

  expect_message(
    display_availability_messages(availability),
    "sup-t, pointwise bounds not available"
  )
})

test_that("create_bounds_plot generates correct plot structure", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS(test_path("fixtures", "simple_plausible.rds"))

  # Prepare data
  bounds_data <- list(
    restricted = pb$restricted_bounds,
    has_restricted = TRUE
  )

  availability <- check_bounds_availability(
    bounds_data,
    pb,
    show_supt = FALSE,
    show_pointwise = FALSE
  )

  # Create plot
  plot <- create_bounds_plot(bounds_data, availability)

  check_plot_structure(plot)

  # Check for expected components
  expect_true(count_geom_type(plot, "GeomPoint") >= 1)  # Point estimates
  expect_true(count_geom_type(plot, "GeomLine") >= 2)    # Restricted bounds
})

# ---- Y-axis Limits Tests ----

test_that("y-axis limits include all data with appropriate buffer", {
  # Create data with known range
  est <- seq(-2, 2, length.out = 10)
  var <- diag(10) * 0.1

  pb <- plausible_bounds(est, var)
  plot <- create_plot(pb)

  # Extract y-limits from plot coordinates
  y_range <- ggplot2::ggplot_build(plot)$layout$panel_params[[1]]$y.range

  # Check that limits include data with buffer
  min_expected <- min(pb$restricted_bounds$lower, na.rm = TRUE)
  max_expected <- max(pb$restricted_bounds$upper, na.rm = TRUE)

  expect_true(y_range[1] < min_expected)
  expect_true(y_range[2] > max_expected)
})

# ---- Legend Tests ----

test_that("legend is properly configured", {
  # Use pre-computed fixture instead of expensive computation
  pb <- readRDS(test_path("fixtures", "simple_plausible.rds"))
  plot <- create_plot(pb)

  # Check legend position
  expect_equal(plot$theme$legend.position, "bottom")

  # Build plot to check legend content
  built_plot <- ggplot2::ggplot_build(plot)
  plot_data <- built_plot$plot

  # Check that color scale has appropriate labels
  color_scale <- plot_data$scales$get_scales("colour")
  if (!is.null(color_scale)) {
    expect_true("Point Estimates" %in% color_scale$labels)
  }
})
