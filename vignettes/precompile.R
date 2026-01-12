#!/usr/bin/env Rscript

# Precompile vignettes script
# This script pre-computes vignettes to avoid re-running code during package builds
# Based on recommendations from: https://ropensci.org/blog/2019/12/08/precompute-vignettes/

# Set working directory to vignettes folder if not already there
if (!grepl("vignettes$", getwd())) {
  if (file.exists("vignettes")) {
    setwd("vignettes")
  } else {
    stop("Please run this script from the package root or vignettes directory")
  }
}

# Function to precompile a single vignette
precompile_vignette <- function(input_file, output_file = NULL) {
  if (is.null(output_file)) {
    output_file <- sub("\\.orig$", "", input_file)
  }
  
  if (!file.exists(input_file)) {
    stop(paste("Input file not found:", input_file))
  }
  
  cat("Precompiling:", input_file, "->", output_file, "\n")
  
  # Set up knitr options for figure paths
  # Use named chunks for better figure naming
  knitr::opts_chunk$set(
    fig.path = "figs/",
    cache = FALSE
  )
  
  # Knit the vignette
  tryCatch({
    knitr::knit(input_file, output = output_file)
    cat("Successfully precompiled:", output_file, "\n")
  }, error = function(e) {
    cat("Error precompiling", input_file, ":\n")
    cat(conditionMessage(e), "\n")
    stop(e)
  })
}

# List of vignettes to precompile
vignettes_to_precompile <- list(
  list(input = "documentation.Rmd.orig", output = "documentation.Rmd")
  # Add more vignettes here as needed:
  # list(input = "another_vignette.Rmd.orig", output = "another_vignette.Rmd")
)

# Main execution
cat("========================================\n")
cat("Precompiling vignettes for plausibounds\n")
cat("========================================\n\n")

# Check if required packages are available
if (!requireNamespace("knitr", quietly = TRUE)) {
  stop("Package 'knitr' is required but not installed")
}

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Package 'rmarkdown' is required but not installed")
}

# Precompile each vignette
for (vignette in vignettes_to_precompile) {
  precompile_vignette(vignette$input, vignette$output)
  cat("\n")
}

# List any generated figure files
figure_files <- list.files(path = "figs", pattern = "\\.(png|pdf|jpg|jpeg|svg)$", 
                          recursive = FALSE)
if (length(figure_files) > 0) {
  cat("Generated figure files:\n")
  for (fig in figure_files) {
    cat("  -", fig, "\n")
  }
} else {
  cat("No figure files were generated.\n")
}

cat("\n========================================\n")
cat("Vignette precompilation complete!\n")
cat("========================================\n")
cat("\nREMINDER: Remember to re-run this script when:\n")
cat("  - The package code changes significantly\n")
cat("  - Before releasing a new version\n")
cat("  - If the vignette content needs updating\n")