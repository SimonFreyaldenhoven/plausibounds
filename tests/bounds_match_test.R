library(testthat)
library(readr)

# Source the R function
#source("R/full_eventplot_l2tf.R")

set.seed(42)


# Simple function to extract test case stems from directory
find_test_cases <- function(dir = "test_across_specs") {
  # Get all CSV files in directory
  files <- list.files(dir, pattern = "\\.csv$")

  stems <- gsub("^[^_]+_(.+)\\.csv$", "\\1", files)
  
  stems <- gsub("_suptbands$", "", stems)
  
  unique(stems)
}

# Initialize test log
init_test_log <- function(logfile = "test_results.log") {
  # Create header
  header <- paste("Test Log Started:", Sys.time(), "\n", 
                  paste(rep("=", 50), collapse = ""), "\n")
  writeLines(header, logfile)
  return(logfile)
}

# Modified test function with logging
test_r_matlab_diagnostics <- function(test_case, logfile = "test_results.log") {
  
  # Log test start
  cat(paste("Starting test case:", test_case, "at", Sys.time()), "\n", file = logfile, append = TRUE)
  
  tryCatch({
    # Read input files
    cat("Reading input files...\n")
    delta <- as.matrix(read_csv(paste0("test_across_specs/delta_", test_case, ".csv"), 
                                col_names = FALSE))
    
    
    dhat <- as.matrix(read_csv(paste0("test_across_specs/dhat_", test_case, ".csv"), 
                               col_names = FALSE))
    vhat_raw <- as.matrix(read_csv(paste0("test_across_specs/vhat_", test_case, ".csv"), 
                                   col_names = FALSE))
    
    # Convert vhat from CSV format to matrix
    p <- nrow(delta)
    vhat <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        vhat[i, j] <- vhat_raw[i, j]
      }
    }
    
    # Read Matlab diagnostics
    cat("Reading Matlab diagnostics...\n")
    matlab_diag <- read_csv(paste0("test_across_specs/diagnostics_", test_case, ".csv"))
    matlab_diag_values <- setNames(matlab_diag$Value, matlab_diag$Field)
    
    matlab_suptbands <- as.matrix(read_csv(paste0("test_across_specs/diagnostics_", 
                                                  test_case, "_suptbands.csv"), 
                                           col_names = FALSE))
    
    # Call R function
    cat("Running R function...\n")
    r_results <- full_eventplot_l2tf(dhat, vhat, nsim = 2000)
    
    # Extract diagnostics from R results
    # Calculate widths for comparison
    r_pw_width <- r_results$diagnostics$pw_width
    r_supt_width <- r_results$diagnostics$supt_width
    r_restricted_width <- mean(r_results$restricted_UB - r_results$restricted_LB)
    r_wald_width <- (r_results$ub - r_results$lb) / p
    
    # Get df from best_fit
    r_df <- if (!is.null(r_results$diagnostics$df)) r_results$diagnostics$df else length(r_results$best_fit$d)
    
    # Compare results
    cat("\nComparing diagnostics...\n")
    
    # Create a tolerance for floating point comparisons
    tol <- 1e-2  # Using a slightly larger tolerance due to simulation randomness
    
    # Track test results
    test_results <- list()
    
    # Test scalar diagnostics
    test_results$pw_width <- tryCatch({
      expect_equal(r_results$diagnostics$pw_width, unname(matlab_diag_values["pw_width"]), tolerance = tol)
      "PASSED"
    }, error = function(e) paste("FAILED:", e$message))
    
    test_results$supt_width <- tryCatch({
      expect_equal(r_results$diagnostics$supt_width, unname(matlab_diag_values["supt_width"]), tolerance = tol)
      "PASSED"
    }, error = function(e) paste("FAILED:", e$message))
    
    test_results$restricted_width <- tryCatch({
      expect_equal(r_restricted_width, unname(matlab_diag_values["restricted_width"]), tolerance = tol)
      "PASSED"
    }, error = function(e) paste("FAILED:", e$message))
    
    test_results$wald_width <- tryCatch({
      expect_equal(r_wald_width, unname(matlab_diag_values["wald_width"]), tolerance = tol)
      "PASSED"
    }, error = function(e) paste("FAILED:", e$message))
    
    test_results$df <- tryCatch({
      expect_equal(r_df, unname(matlab_diag_values["df"]), tolerance = tol)
      "PASSED"
    }, error = function(e) paste("FAILED:", e$message))
    
    # Log individual test results
    for (test_name in names(test_results)) {
      result_line <- paste("  ", test_name, ":", test_results[[test_name]])
      cat(result_line, "\n", file = logfile, append = TRUE)
    }
    
    # Overall result
    all_passed <- all(sapply(test_results, function(x) x == "PASSED"))
    overall_result <- if (all_passed) "ALL TESTS PASSED" else "SOME TESTS FAILED"
    
    cat(paste("Test case", test_case, "completed:", overall_result, "at", Sys.time()), "\n", 
        file = logfile, append = TRUE)
    cat(paste(rep("-", 30), collapse = ""), "\n", file = logfile, append = TRUE)
    
    
    
    cat("\nAll tests completed.\n")
    
    # Create a comparison table with Statistic, Matlab, and R columns
    summary_table <- data.frame(
      Statistic = c("pw_width", "supt_width", "restricted_width", "wald_width", "df"),
      Matlab = c(
        matlab_diag_values["pw_width"],
        matlab_diag_values["supt_width"], 
        matlab_diag_values["restricted_width"],
        matlab_diag_values["wald_width"],
        matlab_diag_values["df"]
      ),
      R = c(
        r_pw_width,
        r_supt_width,
        r_restricted_width,
        r_wald_width,
        r_df
      ),
      Test_Result = sapply(test_results, function(x) if(x == "PASSED") "✓" else "✗")
    )
    
    
    # Print the summary table to log
    cat("Summary Table:\n", file = logfile, append = TRUE)
    capture.output(print(summary_table), file = logfile, append = TRUE)
    cat("\n", file = logfile, append = TRUE)
    
    cat(paste(rep("-", 30), collapse = ""), "\n", file = logfile, append = TRUE)
    flush(file(logfile, open = "a"))  # Force write to disk
    

    # Return a summary of the comparison
    return(summary_table)
    
  }, error = function(e) {
    # Log error
    error_msg <- paste("ERROR in test case", test_case, ":", e$message, "at", Sys.time())
    cat(error_msg, "\n", file = logfile, append = TRUE)
    cat(paste(rep("-", 30), collapse = ""), "\n", file = logfile, append = TRUE)
    stop(e)
  })
}

init_test_log("test_results.log")


test_cases <- find_test_cases()
for(spec in test_cases){
  test_r_matlab_diagnostics(spec) 
}


