#' Create Plot for Plausible Bounds
#'
#' This function creates a plot of plausible bounds from a plausible_bounds object.
#' The plot displays plausible bounds as the main visualization, with optional
#' pointwise and sup-t bounds overlays.
#' Supports event study designs with pre-treatment periods.
#'
#' @param result A plausible_bounds object returned by the plausible_bounds() function
#' @param show_supt Whether to show sup-t bounds (default: TRUE)
#' @param show_pointwise Whether to show pointwise bounds (default: TRUE)
#' @param show_annotations Whether to show annotations with test statistics and ATE (default: TRUE)
#'
#' @return A ggplot2 object
#'
#' @examples
#'
#' # Example with bighump estimates and correlated errors
#' data(estimates_bighump)
#' data(var_bighump)
#' result_complex <- plausible_bounds(estimates_bighump[1:4], var_bighump[1:4, 1:4])
#' plot_complex <- create_plot(result_complex)
#'
#' @importFrom magrittr %>%
#' @export
create_plot <- function(result,
                       show_supt = TRUE,
                       show_pointwise = TRUE,
                       show_annotations = TRUE) {

  # Validate input
  if (!inherits(result, "plausible_bounds")) {
    stop("create_plot() requires a plausible_bounds object. Use plausible_bounds() to create one.")
  }

  # Check that restricted bounds are available
  if (is.null(result$restricted_bounds)) {
    stop("No restricted bounds available in the plausible_bounds object.")
  }

  # Prepare bounds data
  bounds_data <- list(
    restricted = result$restricted_bounds,
    has_restricted = TRUE
  )

  # Check availability of optional bounds
  chk <- check_bounds_availability(bounds_data, result, show_supt, show_pointwise)

  # Display informative messages about missing optional bounds
  bounds_data <- chk$bounds_data
  availability <- chk$availability
  
  display_availability_messages(availability)

  # Extract test statistics for annotations
  annotations <- NULL
  if (show_annotations) {
    annotations <- list()
    if (!is.null(result$avg_treatment_effect)) {
      annotations$ate <- result$avg_treatment_effect
    }
    if (!is.null(result$wald_test$pre)) {
      annotations$Wpre <- result$wald_test$pre
    }
    if (!is.null(result$wald_test$post)) {
      annotations$Wpost <- result$wald_test$post
    }
  }

  # Create the plot with available data
  create_bounds_plot(bounds_data, availability, annotations)
}

# Helper function to check what bounds are available and requested
check_bounds_availability <- function(bounds_data, result, show_supt, show_pointwise) {
  n <- nrow(bounds_data$restricted)
  
  pointwise_ok <- !is.null(result$pointwise_bounds) &&
    all(c("lower", "upper") %in% names(result$pointwise_bounds)) &&
    length(result$pointwise_bounds$lower) == n &&
    length(result$pointwise_bounds$upper) == n
  
  supt_ok <- !is.null(result$supt_bounds) &&
    all(c("lower", "upper") %in% names(result$supt_bounds)) &&
    length(result$supt_bounds$lower) == n &&
    length(result$supt_bounds$upper) == n
  
  if (pointwise_ok) {
    bounds_data$restricted$pointwise_lower <- result$pointwise_bounds$lower
    bounds_data$restricted$pointwise_upper <- result$pointwise_bounds$upper
  }
  
  if (supt_ok) {
    bounds_data$restricted$supt_lower <- result$supt_bounds$lower
    bounds_data$restricted$supt_upper <- result$supt_bounds$upper
  }
  
  availability <- list(
    show_restricted = TRUE,
    supt_requested = show_supt,
    supt_available = supt_ok,
    pointwise_requested = show_pointwise,
    pointwise_available = pointwise_ok,
    show_supt = show_supt && supt_ok,
    show_pointwise = show_pointwise && pointwise_ok
  )
  
  list(bounds_data = bounds_data, availability = availability)
}


# Helper function to display informative messages about missing optional bounds
display_availability_messages <- function(availability) {
  missing_types <- c()

  # Check for requested but unavailable optional bounds
  if (availability$supt_requested && !availability$supt_available) {
    missing_types <- c(missing_types, "sup-t")
  }

  if (availability$pointwise_requested && !availability$pointwise_available) {
    missing_types <- c(missing_types, "pointwise")
  }

  # Display message if any requested bounds are not available
  if (length(missing_types) > 0) {
    msg <- paste("Note:", paste(missing_types, collapse = ", "), "bounds not available in data.")

    param_names <- gsub("sup-t", "supt", missing_types)
    param_names <- paste0("show_", param_names, " = FALSE")

    msg <- paste(msg, "Set", paste(param_names, collapse = " and "), "to silence this message.")

    message(msg)
  }
}


# Helper function to create the actual plot
create_bounds_plot <- function(bounds_data, availability, annotations = NULL) {

  colors <- c(
      estimate = "black",
      surrogate = "cornflowerblue",
      surrogate_bounds = "cornflowerblue",
      pointwise = "darkgray",
      supt = "gray"
    )

  line_types <- c(
      estimate = "solid",
      surrogate = "solid",
      surrogate_bounds = "dashed",
      pointwise = "solid",
      supt = "solid"
  )

  # Use restricted bounds data
  df <- bounds_data$restricted

  # Check if we have pre-periods (negative horizon values)
  has_preperiods <- any(df$horizon < 0)

  # Separate post-period data for surrogate and bounds
  df_post <- df[df$horizon > 0, ]

  # Create base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    ggplot2::geom_point(ggplot2::aes(y = coef, color = "estimate"))

  # Add surrogate and bounds only for post-periods
  if (nrow(df_post) > 0 && !all(is.na(df_post$surrogate))) {
    p <- p +
      ggplot2::geom_line(data = df_post, ggplot2::aes(y = surrogate, color = "surrogate", linetype = "surrogate")) +
      ggplot2::geom_line(data = df_post, ggplot2::aes(y = lower, color = "surrogate_bounds", linetype = "surrogate_bounds")) +
      ggplot2::geom_line(data = df_post, ggplot2::aes(y = upper, color = "surrogate_bounds", linetype = "surrogate_bounds"))
  }

  # Add vertical line at event time 0 if we have pre-periods
  if (has_preperiods) {
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7)
  }
  
  # Add pointwise bounds if requested and available
  if (availability$show_pointwise && "pointwise_lower" %in% names(df) && "pointwise_upper" %in% names(df)) {
    for (i in 1:nrow(df)) {
      p <- p + ggplot2::geom_segment(
        x = df$horizon[i] - 0.2, xend = df$horizon[i] + 0.2,
        y = df$pointwise_lower[i], yend = df$pointwise_lower[i],
        color = colors["pointwise"], linewidth = 0.5
      ) +
      ggplot2::geom_segment(
        x = df$horizon[i] - 0.2, xend = df$horizon[i] + 0.2,
        y = df$pointwise_upper[i], yend = df$pointwise_upper[i],
        color = colors["pointwise"], linewidth = 0.5
      )
    }
  }
  
  # Add sup-t bounds if requested and available
  if (availability$show_supt && "supt_lower" %in% names(df) && "supt_upper" %in% names(df)) {
    for (i in 1:nrow(df)) {
      p <- p + ggplot2::geom_segment(
        x = df$horizon[i], xend = df$horizon[i],
        y = df$supt_lower[i], yend = df$supt_upper[i],
        color = colors["supt"], linewidth = 0.5
      )
    }
  }

  # Calculate y-axis limits to include all bounds
  y_min <- min(df$coef, df$lower, na.rm = TRUE)
  y_max <- max(df$coef, df$upper, na.rm = TRUE)
  
  # Include pointwise bounds in y-axis limits if available
  if ("pointwise_lower" %in% names(df) && "pointwise_upper" %in% names(df)) {
    y_min <- min(y_min, min(df$pointwise_lower, na.rm = TRUE))
    y_max <- max(y_max, max(df$pointwise_upper, na.rm = TRUE))
  }
  
  # Include sup-t bounds in y-axis limits if available
  if ("supt_lower" %in% names(df) && "supt_upper" %in% names(df)) {
    y_min <- min(y_min, min(df$supt_lower, na.rm = TRUE))
    y_max <- max(y_max, max(df$supt_upper, na.rm = TRUE))
  }
  
  # Add a small buffer to the y-axis limits
  y_buffer <- 0.05 * (y_max - y_min)
  y_min <- y_min - y_buffer
  y_max <- y_max + y_buffer
  
  # Add theme and labels
  p <- p +
    ggplot2::labs(
      x = "Event time",
      y = "Estimate"
    ) +
    # Set y-axis limits to include all bounds
    ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
    # Ensure x-axis has appropriate integer breaks (avoid clutter for large horizons)
    ggplot2::scale_x_continuous(breaks = function(x) {
      range <- ceiling(max(x)) - floor(min(x))
      if (range <= 10) {
        # For small ranges, show all integer breaks
        return(seq(floor(min(x)), ceiling(max(x)), by = 1))
      } else if (range <= 20) {
        # For medium ranges, show every other integer
        return(seq(floor(min(x)), ceiling(max(x)), by = 2))
      } else {
        # For large ranges, use approximately 8-10 breaks
        by_value <- ceiling(range / 8)
        return(seq(floor(min(x)), ceiling(max(x)), by = by_value))
      }
    }) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = colors,
      breaks = c("estimate", "surrogate"),
      labels = c("Point Estimates", "Restricted")
    ) +
    # Hide linetype from legend
    ggplot2::guides(linetype = "none")

  p <- p + ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom")

  # Add annotations if provided
  if (!is.null(annotations) && length(annotations) > 0) {
    annotation_parts <- c()

    # Add pretrends p-value if available
    if (!is.null(annotations$Wpre)) {
      annotation_parts <- c(annotation_parts,
        sprintf("Pretrends p-value: %.2f", annotations$Wpre$p_value))
    }

    # Add no effect p-value if available
    if (!is.null(annotations$Wpost)) {
      annotation_parts <- c(annotation_parts,
        sprintf("No effect p-value: %.2f", annotations$Wpost$p_value))
    }

    # Add ATE with CI if available
    if (!is.null(annotations$ate)) {
      annotation_parts <- c(annotation_parts,
        sprintf("Average effect [CI]: %.3g [%.3g, %.3g]",
                annotations$ate$estimate,
                annotations$ate$lower,
                annotations$ate$upper))
    }

    if (length(annotation_parts) > 0) {
      annotation_text <- paste(annotation_parts, collapse = " --- ")
      p <- p + ggplot2::labs(caption = annotation_text) +
        ggplot2::theme(
          plot.caption = ggplot2::element_text(hjust = 0.5, size = 8)
        )
    }
  }

  return(p)
}
