#' Create Plot for Plausible Bounds
#'
#' This function creates a plot of plausible bounds from the results of plausible_bounds,
#' calculate_cumulative_bounds, or calculate_restricted_bounds functions.
#'
#' @param ... One or two objects: either a plausible_bounds object or individual bounds objects
#' @param show_cumulative Whether to show cumulative bounds (default: TRUE)
#' @param show_restricted Whether to show restricted bounds (default: TRUE)
#' @param show_supt Whether to show sup-t bounds (default: TRUE)
#' @param show_pointwise Whether to show pointwise bounds (default: TRUE)
#'
#' @return A ggplot2 object
#'
#' @examples
#' # Example with constant estimates and IID errors
#' data(estimates_constant_iid)
#' data(var_constant_iid)
#' result <- plausible_bounds(estimates_constant_iid, var_constant_iid)
#' plot <- create_plot(result)
#'
#' # Example with wiggly estimates and strong correlation
#' data(estimates_wiggly_strong_corr)
#' data(var_wiggly_strong_corr)
#' result_complex <- plausible_bounds(estimates_wiggly_strong_corr, var_wiggly_strong_corr)
#' plot_complex <- create_plot(result_complex)
#'
#' # Example showing only restricted bounds
#' plot_restricted <- create_plot(result_complex, show_cumulative = FALSE)
#'
#' @importFrom magrittr %>%
#' @export
create_plot <- function(..., 
                       show_cumulative = TRUE,
                       show_restricted = TRUE,
                       show_supt = TRUE,
                       show_pointwise = TRUE) {
  
  # Collect all arguments
  args <- list(...)
  
  # Handle different input patterns and extract bounds data
  bounds_data <- extract_bounds_data(args, show_cumulative, show_restricted, show_supt, show_pointwise)
  # Check what's available and what's requested
  availability <- check_bounds_availability(bounds_data, show_cumulative, show_restricted, show_supt, show_pointwise)
  # Display informative messages about missing bounds
  display_availability_messages(availability)
  
  # Create the plot with available data
  create_bounds_plot(bounds_data, availability)
}

# Helper function to merge pointwise and supt bounds with main bounds data frames
merge_bounds_data <- function(result, x, type) {
  # Process cumulative bounds if available
  if (result$has_cumulative && !is.null(result$cumulative)) {
    df <- result$cumulative
    
    # Check if pointwise bounds are available in the original object
    if (!is.null(x$pointwise_bounds)) {
      df$pointwise_lower <- x$pointwise_bounds$lower
      df$pointwise_upper <- x$pointwise_bounds$upper
    }
    
    # Check if supt bounds are available in the original object
    if (!is.null(x$supt_bounds)) {
      df$supt_lower <- x$supt_bounds$lower
      df$supt_upper <- x$supt_bounds$upper
    }
    
    result$cumulative <- df
  }
  
  # Process restricted bounds if available
  if (result$has_restricted && !is.null(result$restricted)) {
    df <- result$restricted
    
    # For plausible_bounds objects, we need to handle the structure differently
    if (type == "plausible_bounds") {
      # Check if pointwise bounds are available in the original object
      if (!is.null(x$pointwise_bounds)) {
        df$pointwise_lower <- x$pointwise_bounds$lower
        df$pointwise_upper <- x$pointwise_bounds$upper
      }
      
      # Check if supt bounds are available in the original object
      if (!is.null(x$supt_bounds)) {
        df$supt_lower <- x$supt_bounds$lower
        df$supt_upper <- x$supt_bounds$upper
      }
    } else {
      # For restricted_bounds objects
      # Check if pointwise bounds are available in the original object
      if (!is.null(x$pointwise_bounds)) {
        df$pointwise_lower <- x$pointwise_bounds$lower
        df$pointwise_upper <- x$pointwise_bounds$upper
      }
      
      # Check if supt bounds are available in the original object
      if (!is.null(x$supt_bounds)) {
        df$supt_lower <- x$supt_bounds$lower
        df$supt_upper <- x$supt_bounds$upper
      }
    }
    
    result$restricted <- df
  }
  
  return(result)
}

# Helper function to extract bounds data from input arguments
extract_bounds_data <- function(args, show_cumulative, show_restricted, show_supt, show_pointwise) {
  if (length(args) == 1) {
    x <- args[[1]]
    
    if (inherits(x, "plausible_bounds")) {
      # Extract from plausible_bounds object
      result <- list(type = "plausible_bounds")
      
      if (show_cumulative && !is.null(x$cumulative_bounds)) {
        result$cumulative <- x$cumulative_bounds  # Direct access, not $bounds
        result$has_cumulative <- TRUE
      } else {
        result$has_cumulative <- FALSE
      }
      
      if (show_restricted && !is.null(x$restricted_bounds)) {
        result$restricted <- x$restricted_bounds  # Direct access, not $bounds
        result$has_restricted <- TRUE
      } else {
        result$has_restricted <- FALSE
      }
      
      # Merge pointwise and supt bounds if requested
      if ((show_pointwise || show_supt) && 
          (!is.null(x$pointwise_bounds) || !is.null(x$supt_bounds))) {
        result <- merge_bounds_data(result, x, "plausible_bounds")
      }
      
      return(result)
      
    } else if (inherits(x, "cumulative_bounds")) {
      # Extract from cumulative_bounds object
      result <- list(
        type = "cumulative_only",
        has_cumulative = TRUE,
        has_restricted = FALSE
      )
      
      if (show_cumulative) {
        result$cumulative <- x$cumulative_bounds  # Direct access to cumulative_bounds
      }
      
      # Merge pointwise and supt bounds if requested
      if ((show_pointwise || show_supt) && 
          (!is.null(x$pointwise_bounds) || !is.null(x$supt_bounds))) {
        result <- merge_bounds_data(result, x, "cumulative_bounds")
      }
      
      return(result)
      
    } else if (inherits(x, "restricted_bounds")) {
      # Extract from restricted_bounds object
      result <- list(
        type = "restricted_only",
        has_cumulative = FALSE,
        has_restricted = TRUE
      )
      
      if (show_restricted) {
        result$restricted <- x$restricted_bounds  # Direct access to restricted_bounds
      }
      
      # Merge pointwise and supt bounds if requested
      if ((show_pointwise || show_supt) && 
          (!is.null(x$pointwise_bounds) || !is.null(x$supt_bounds))) {
        result <- merge_bounds_data(result, x, "restricted_bounds")
      }
      
      return(result)
      
    } else {
      stop("Unrecognized input type. Expected plausible_bounds, cumulative_bounds, or restricted_bounds object.")
    }
    
  } else if (length(args) == 2) {
    # Handle two separate bounds objects
    result <- list(type = "combined")
    result$has_cumulative <- FALSE
    result$has_restricted <- FALSE
    
    for (arg in args) {
      if (inherits(arg, "cumulative_bounds")) {
        if (show_cumulative) {
          result$cumulative <- arg$cumulative_bounds  # Direct access to cumulative_bounds
        }
        result$has_cumulative <- TRUE
        
        # Merge pointwise and supt bounds if requested
        if ((show_pointwise || show_supt) && 
            (!is.null(arg$pointwise_bounds) || !is.null(arg$supt_bounds))) {
          result <- merge_bounds_data(result, arg, "cumulative_bounds")
        }
      } else if (inherits(arg, "restricted_bounds")) {
        if (show_restricted) {
          result$restricted <- arg$restricted_bounds  # Direct access to restricted_bounds
        }
        result$has_restricted <- TRUE
        
        # Merge pointwise and supt bounds if requested
        if ((show_pointwise || show_supt) && 
            (!is.null(arg$pointwise_bounds) || !is.null(arg$supt_bounds))) {
          result <- merge_bounds_data(result, arg, "restricted_bounds")
        }
      } else {
        stop("When providing two arguments, both must be either cumulative_bounds or restricted_bounds objects.")
      }
    }
    
    if (!result$has_cumulative && !result$has_restricted) {
      stop("At least one valid bounds object must be provided.")
    }
    
    return(result)
    
  } else {
    stop("create_plot() accepts 1 or 2 arguments")
  }
}

# Helper function to check what bounds are available and requested
check_bounds_availability <- function(bounds_data, show_cumulative, show_restricted, show_supt, show_pointwise) {
  availability <- list(
    cumulative_requested = show_cumulative,
    cumulative_available = bounds_data$has_cumulative && !is.null(bounds_data$cumulative),
    restricted_requested = show_restricted,
    restricted_available = bounds_data$has_restricted && !is.null(bounds_data$restricted),
    supt_requested = show_supt,
    supt_available = FALSE,
    pointwise_requested = show_pointwise,
    pointwise_available = FALSE
  )
  
  # Check for sup-t availability in the merged data
  if (availability$cumulative_available && !is.null(bounds_data$cumulative)) {
    if ("supt_lower" %in% names(bounds_data$cumulative) && "supt_upper" %in% names(bounds_data$cumulative)) {
      availability$supt_available <- TRUE
    }
  }
  if (availability$restricted_available && !is.null(bounds_data$restricted)) {
    if ("supt_lower" %in% names(bounds_data$restricted) && "supt_upper" %in% names(bounds_data$restricted)) {
      availability$supt_available <- TRUE
    }
  }
  
  # Check for pointwise availability in the merged data
  if (availability$cumulative_available && !is.null(bounds_data$cumulative)) {
    if ("pointwise_lower" %in% names(bounds_data$cumulative) && "pointwise_upper" %in% names(bounds_data$cumulative)) {
      availability$pointwise_available <- TRUE
    }
  }
  if (availability$restricted_available && !is.null(bounds_data$restricted)) {
    if ("pointwise_lower" %in% names(bounds_data$restricted) && "pointwise_upper" %in% names(bounds_data$restricted)) {
      availability$pointwise_available <- TRUE
    }
  }
  
  # Determine what will actually be shown
  availability$show_cumulative <- availability$cumulative_requested && availability$cumulative_available
  availability$show_restricted <- availability$restricted_requested && availability$restricted_available
  availability$show_supt <- availability$supt_requested && availability$supt_available
  availability$show_pointwise <- availability$pointwise_requested && availability$pointwise_available
  
  return(availability)
}

# Helper function to display informative messages about missing bounds
display_availability_messages <- function(availability) {
  missing_types <- c()
  
  # Check for requested but unavailable bounds
  if (availability$cumulative_requested && !availability$cumulative_available) {
    missing_types <- c(missing_types, "cumulative")
  }
  
  if (availability$restricted_requested && !availability$restricted_available) {
    missing_types <- c(missing_types, "restricted")
  }
  
  if (availability$supt_requested && !availability$supt_available) {
    missing_types <- c(missing_types, "sup-t")
  }
  
  if (availability$pointwise_requested && !availability$pointwise_available) {
    missing_types <- c(missing_types, "pointwise")
  }
  
  # Display message if any requested bounds are not available
  if (length(missing_types) > 0) {
    message("Note: ", paste(missing_types, collapse = ", "), " bounds not available in data.")
  }
  
  # Check if nothing will be plotted
  if (!availability$show_cumulative && !availability$show_restricted) {
    warning("No bounds available to plot. Check your input data and parameters.")
  }
}


# Helper function to create the actual plot
create_bounds_plot <- function(bounds_data, availability) {
  
  colors <- c(
      estimate = "black",
      surrogate = "cornflowerblue",
      surrogate_bounds = "cornflowerblue",
      pointwise = "darkgray",
      supt = "gray",
      restricted = "cornflowerblue",
      cumulative = "orange"
    )
  
  line_types <- c(
      estimate = "solid",
      surrogate = "solid",
      surrogate_bounds = "dashed",
      pointwise = "solid",
      supt = "solid",
      restricted = "solid",
      cumulative = "solid"
  )
  

  p <- NULL
  df <- NULL
  
  # Determine which data to use for the base plot
  if (availability$show_cumulative && availability$show_restricted) {
    # Plot both cumulative and restricted bounds
    df_c <- bounds_data$cumulative
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
    
    # Copy pointwise and supt bounds from restricted to combined if they exist
    if ("pointwise_lower" %in% names(df_r) && "pointwise_upper" %in% names(df_r)) {
      df$pointwise_lower <- df_r$pointwise_lower
      df$pointwise_upper <- df_r$pointwise_upper
    }
    if ("supt_lower" %in% names(df_r) && "supt_upper" %in% names(df_r)) {
      df$supt_lower <- df_r$supt_lower
      df$supt_upper <- df_r$supt_upper
    }
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = "cumulative"), alpha = 0.2) +
      ggplot2::geom_point(ggplot2::aes(y = coef, color = "estimate")) +
      ggplot2::geom_line(ggplot2::aes(y = surrogate, color = "surrogate", linetype = "surrogate")) +
      ggplot2::geom_line(ggplot2::aes(y = restricted_lower, color = "surrogate_bounds", linetype = "surrogate_bounds")) +
      ggplot2::geom_line(ggplot2::aes(y = restricted_upper, color = "surrogate_bounds", linetype = "surrogate_bounds"))
    
  } else if (availability$show_cumulative) {
    # Plot only cumulative bounds
    df <- bounds_data$cumulative 
    p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon)) +
      ggplot2::geom_point(ggplot2::aes(y = coef, color = "estimate")) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = "cumulative"), alpha = 0.2)
    
  } else if (availability$show_restricted) {
    # Plot only restricted bounds
    df <- bounds_data$restricted
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = horizon)) +
      ggplot2::geom_point(ggplot2::aes(y = coef, color = "estimate")) +
      ggplot2::geom_line(ggplot2::aes(y = surrogate, color = "surrogate", linetype = "surrogate")) +
      ggplot2::geom_line(ggplot2::aes(y = lower, color = "surrogate_bounds", linetype = "surrogate_bounds")) +
      ggplot2::geom_line(ggplot2::aes(y = upper, color = "surrogate_bounds", linetype = "surrogate_bounds"))
    
  } else {
    stop("No bounds available to plot.")
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
      x = "Horizon",
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
    ggplot2::scale_fill_manual(
      name = NULL,
      values = colors,
      breaks = c("cumulative"),
      labels = c("Cumulative")
    ) +
    # Hide linetype from legend
    ggplot2::guides(linetype = "none")
    
  p <- p + ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom")
  
  
  return(p)
}
