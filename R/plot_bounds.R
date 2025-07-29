#' Plot Plausible Bounds
#'
#' This function creates a plot of plausible bounds from the results of plausible_bounds,
#' calculate_cumulative_bounds, or calculate_restricted_bounds functions.
#'
#' @param ... One or two objects: either a plausible_bounds object or individual bounds objects
#' @param which Which bounds to plot: "both", "cumulative", or "restricted" (default: "both")
#' @param show_pointwise Whether to show pointwise bounds (default: TRUE)
#' @param show_supt Whether to show sup-t bounds (default: TRUE)
#' @param theme Optional ggplot2 theme to use
#' @param colors Optional vector of colors for different elements
#' @param line_types Optional vector of line types for different elements
#'
#' @return A ggplot2 object
#'
#' @export
plot_bounds <- function(..., which = c("both", "cumulative", "restricted"), 
                       show_pointwise = TRUE, show_supt = TRUE,
                       theme = NULL, colors = NULL, line_types = NULL) {
  which <- match.arg(which)
  
  # Collect all arguments
  args <- list(...)
  
  # Handle different input patterns
  if (length(args) == 1) {
    x <- args[[1]]
    bounds_data <- detect_and_extract(x, which)
  } else if (length(args) == 2) {
    # Do not assume particular order of args
    bounds_data <- combine_separate_results(args[[1]], args[[2]], which)
  } else {
    stop("plot_bounds() accepts 1 or 2 arguments")
  }
  
  # Create the plot
  create_bounds_plot(bounds_data, show_pointwise, show_supt, theme, colors, line_types)
}

# Helper function to detect and extract data from different object types
detect_and_extract <- function(x, which) {
  if (inherits(x, "plausible_bounds")) {
    # Extract from plausible_bounds object based on 'which'
    if (which == "both") {
      return(list(
        cumulative = x$cumulative_bounds$bounds,
        restricted = x$restricted_bounds$bounds,
        type = "both"
      ))
    } else if (which == "cumulative") {
      return(list(
        cumulative = x$cumulative_bounds$bounds,
        type = "cumulative"
      ))
    } else if (which == "restricted") {
      return(list(
        restricted = x$restricted_bounds$bounds,
        type = "restricted"
      ))
    }
  } else if (inherits(x, "cumulative_bounds")) {
    # Extract from cumulative_bounds object
    if (which == "restricted") {
      stop("Cannot show restricted bounds from a cumulative_bounds object")
    }
    return(list(
      cumulative = x$bounds,
      type = "cumulative"
    ))
  } else if (inherits(x, "restricted_bounds")) {
    # Extract from restricted_bounds object
    if (which == "cumulative") {
      stop("Cannot show cumulative bounds from a restricted_bounds object")
    }
    return(list(
      restricted = x$bounds,
      type = "restricted"
    ))
  } else {
    stop("Unrecognized input type")
  }
}

# Helper function to combine data from separate bounds objects
combine_separate_results <- function(x, y, which) {
  # Determine the types of x and y
  x_type <- if (inherits(x, "cumulative_bounds")) "cumulative" else if (inherits(x, "restricted_bounds")) "restricted" else NULL
  y_type <- if (inherits(y, "cumulative_bounds")) "cumulative" else if (inherits(y, "restricted_bounds")) "restricted" else NULL
  
  if (is.null(x_type) || is.null(y_type)) {
    stop("Both arguments must be either cumulative_bounds or restricted_bounds objects")
  }
  
  if (x_type == y_type) {
    stop("Cannot combine two objects of the same type")
  }
  
  # Extract data based on which
  if (which == "both") {
    result <- list(type = "both")
    if (x_type == "cumulative") {
      result$cumulative <- x$bounds
      result$restricted <- y$bounds
    } else {
      result$cumulative <- y$bounds
      result$restricted <- x$bounds
    }
    return(result)
  } else if (which == "cumulative") {
    if (x_type == "cumulative") {
      return(list(cumulative = x$bounds, type = "cumulative"))
    } else {
      return(list(cumulative = y$bounds, type = "cumulative"))
    }
  } else if (which == "restricted") {
    if (x_type == "restricted") {
      return(list(restricted = x$bounds, type = "restricted"))
    } else {
      return(list(restricted = y$bounds, type = "restricted"))
    }
  }
}

# Helper function to create the actual plot
create_bounds_plot <- function(bounds_data, show_pointwise, show_supt, 
                              theme = NULL, colors = NULL, line_types = NULL) {
  # Set default colors if not provided
  if (is.null(colors)) {
    colors <- c(
      estimate = "black",
      surrogate = "green",
      surrogate_bounds = "green",
      pointwise = "darkgray",
      supt = "gray",
      restricted = "green",
      cumulative = "red"
    )
  }
  
  # Set default line types if not provided
  if (is.null(line_types)) {
    line_types <- c(
      estimate = "solid",
      surrogate = "solid",
      surrogate_bounds = "dashed",
      pointwise = "solid",
      supt = "solid",
      restricted = "solid",
      cumulative = "solid"
    )
  }
  
  # Initialize the plot
  if (bounds_data$type == "cumulative") {
    # Plot cumulative bounds
    df <- bounds_data$cumulative %>%
      dplyr::mutate(
        lower = lower / max(horizon),
        upper = upper / max(horizon) 
      )
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon)) +
      ggplot2::geom_point(ggplot2::aes(y = coef, color = "estimate")) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = "cumulative"), alpha = 0.2)
    
    if (show_pointwise && "pointwise_lower" %in% names(df) && "pointwise_upper" %in% names(df)) {
      # Add horizontal ticks for pointwise intervals
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
    
    if (show_supt && "supt_lower" %in% names(df) && "supt_upper" %in% names(df)) {
      # Add vertical lines for sup-t bands
      for (i in 1:nrow(df)) {
        p <- p + ggplot2::geom_segment(
          x = df$horizon[i], xend = df$horizon[i],
          y = df$supt_lower[i], yend = df$supt_upper[i],
          color = colors["supt"], linewidth = 0.5
        )
      }
    }
    
  } else if (bounds_data$type == "restricted") {
    # Plot restricted bounds
    df <- bounds_data$restricted
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon)) +
      ggplot2::geom_point(ggplot2::aes(y = coef, color = "estimate")) +
      ggplot2::geom_line(ggplot2::aes(y = surrogate, color = "surrogate", linetype = "surrogate")) +
      ggplot2::geom_line(ggplot2::aes(y = lower, color = "surrogate_bounds", linetype = "surrogate_bounds")) +
      ggplot2::geom_line(ggplot2::aes(y = upper, color = "surrogate_bounds", linetype = "surrogate_bounds"))
    
    if (show_pointwise && "pointwise_lower" %in% names(df) && "pointwise_upper" %in% names(df)) {
      # Add horizontal ticks for pointwise intervals
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
    
    if (show_supt && "supt_lower" %in% names(df) && "supt_upper" %in% names(df)) {
      # Add vertical lines for sup-t bands
      for (i in 1:nrow(df)) {
        p <- p + ggplot2::geom_segment(
          x = df$horizon[i], xend = df$horizon[i],
          y = df$supt_lower[i], yend = df$supt_upper[i],
          color = colors["supt"], linewidth = 0.5
        )
      }
    }
    
  } else if (bounds_data$type == "both") {
    # Plot both cumulative and restricted bounds
    df_c <- bounds_data$cumulative %>%
     dplyr::mutate(
        lower = lower / max(horizon),
        upper = upper / max(horizon) 
      )
    df_r <- bounds_data$restricted
    
    # Ensure both data frames have the same horizon values
    if (!identical(df_c$horizon, df_r$horizon)) {
      stop("Cumulative and restricted bounds must have the same horizon values")
    }
    
    # Create a combined data frame
    df <- df_c
    df$surrogate <- df_r$surrogate
    df$restricted_lower <- df_r$lower
    df$restricted_upper <- df_r$upper
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = "cumulative"), alpha = 0.2) +
      ggplot2::geom_point(ggplot2::aes(y = coef, color = "estimate")) +
      ggplot2::geom_line(ggplot2::aes(y = surrogate, color = "surrogate", linetype = "surrogate")) +
      ggplot2::geom_line(ggplot2::aes(y = restricted_lower, color = "surrogate_bounds", linetype = "surrogate_bounds")) +
      ggplot2::geom_line(ggplot2::aes(y = restricted_upper, color = "surrogate_bounds", linetype = "surrogate_bounds"))
    
    if (show_pointwise && "pointwise_lower" %in% names(df) && "pointwise_upper" %in% names(df)) {
      # Add horizontal ticks for pointwise intervals
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
    
    if (show_supt && "supt_lower" %in% names(df) && "supt_upper" %in% names(df)) {
      # Add vertical lines for sup-t bands
      for (i in 1:nrow(df)) {
        p <- p + ggplot2::geom_segment(
          x = df$horizon[i], xend = df$horizon[i],
          y = df$supt_lower[i], yend = df$supt_upper[i],
          color = colors["supt"], linewidth = 0.5
        )
      }
    }
  }
  
  # Add theme and labels
  p <- p +
    ggplot2::labs(
      title = "Plausible Bounds",
      x = "Horizon",
      y = "Estimate"
    ) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = colors,
      breaks = c("estimate", "surrogate"),
      labels = c("Point Estimates", "Restricted")
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = colors,
      breaks = c("cumulative"),
      labels = c("Cumulative")
    ) +
    # Hide linetype from legend
    ggplot2::guides(linetype = "none")
  
  # Apply custom theme if provided
  if (!is.null(theme)) {
    p <- p + theme + ggplot2::theme(legend.position = "bottom")
  } else {
    p <- p + ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom")
  }
  
  return(p)
}
