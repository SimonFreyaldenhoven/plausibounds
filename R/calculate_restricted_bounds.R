#' Calculate Plausible Bounds
#'
#' This function calculates the plausible (or "restricted") bounds for a vector of estimates using model selection.
#' Supports pre-treatment periods for event study designs.
#'
#'
#' @param estimates A vector of point estimates. If preperiods > 0, the first preperiods
#'   elements are pre-treatment estimates, followed by post-treatment estimates.
#' @param var The variance-covariance matrix of the estimates
#' @param alpha Significance level (default: 0.05)
#' @param preperiods Number of pre-treatment periods (default: 0). Period 0 is assumed
#'   to be normalized and not included in estimates.
#' @param parallel Whether to use parallel processing (default: FALSE)
#' @param n_cores Number of cores to use for parallel processing (default: NULL, which uses
#'   detectCores() - 1). Only used when parallel = TRUE.
#'
#' @return A list containing:
#'   \item{restricted_bounds}{A data frame with columns for horizon (event time),
#'     unrestricted estimates, restricted estimates, and plausible bounds}
#'   \item{Wpre}{Wald test for no pre-trends (statistic and p-value), if preperiods > 0}
#'   \item{Wpost}{Wald test for no treatment effect (statistic and p-value)}
#'   \item{metadata}{A list with metadata about the calculation}
#' @keywords internal

calculate_restricted_bounds <- function(estimates, var, alpha = 0.05,
                                       preperiods = 0,
                                       parallel = FALSE,
                                       n_cores = NULL) {
  if (!is.numeric(estimates) || !is.vector(estimates)) {
    stop("estimates must be a numeric vector")
  }
  if (length(estimates) - preperiods < 2) {
    stop("calculate_restricted_bounds requires at least 2 estimates")
  }
  if (any(is.na(estimates))) {
    stop("estimates cannot contain NA values")
  }
  if (any(!is.finite(estimates))) {
    stop("estimates cannot contain infinite values")
  }
  if (!is.matrix(var) || nrow(var) != length(estimates) || ncol(var) != length(estimates)) {
    stop("var must be a square matrix with dimensions matching the length of estimates")
  }
  if (any(is.na(var))) {
    stop("var cannot contain NA values")
  }
  if (any(!is.finite(var))) {
    stop("var cannot contain infinite values")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number between 0 and 1")
  }
  if (!is.numeric(preperiods) || length(preperiods) != 1 || preperiods < 0 || preperiods != floor(preperiods)) {
    stop("preperiods must be a non-negative integer")
  }
  if (preperiods >= length(estimates)) {
    stop("preperiods must be less than the length of estimates")
  }
  if (!is.null(n_cores)) {
    if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1 || n_cores != floor(n_cores)) {
      stop("n_cores must be a positive integer")
    }
  }

  # Store full data
  estimatesAll <- estimates
  varAll <- var
  n_all <- length(estimatesAll)


  # Extract post-period data (smoothing is done only on post-periods)
  if (preperiods > 0) {
    estimates_pre <- estimatesAll[1:preperiods]
    var_pre <- varAll[1:preperiods, 1:preperiods, drop = FALSE]
    estimates <- estimatesAll[(preperiods + 1):n_all]
    var <- varAll[(preperiods + 1):n_all, (preperiods + 1):n_all, drop = FALSE]
  }

  p <- length(estimates)  # p is now the number of POST-periods
  
  # Ensure variance matrices are positive definite

  # For post-period variance
  eigs <- eigen(var, symmetric = TRUE)
  if (min(eigs$values) < 0) {
    var <- eigs$vectors %*% (diag(eigs$values) + (abs(min(eigs$values)) + 1e-8) * diag(length(eigs$values))) %*% t(eigs$vectors)
  }

  # For full variance matrix (used for sup-t)
  eigs_all <- eigen(varAll, symmetric = TRUE)
  if (min(eigs_all$values) < 0) {
    varAll <- eigs_all$vectors %*% (diag(eigs_all$values) + (abs(min(eigs_all$values)) + 1e-8) * diag(length(eigs_all$values))) %*% t(eigs_all$vectors)
  }

  # Compute sup-t critical value using ALL periods (pre + post)
  stdv_all <- sqrt(diag(varAll))
  d_nonzero <- stdv_all > .Machine$double.eps
  cV_all <- diag(1 / stdv_all[d_nonzero]) %*% varAll[d_nonzero, d_nonzero] %*% diag(1 / stdv_all[d_nonzero])

  kk <- 10000

  Corrmat <- cV_all
  eig_c <- eigen(Corrmat, symmetric = TRUE)
  Corrmat_sqrt <- eig_c$vectors %*% diag(sqrt(pmax(eig_c$values, 0))) %*% t(eig_c$vectors)
  t_stat <- abs(matrix(stats::rnorm(kk * sum(d_nonzero)), nrow = kk) %*% Corrmat_sqrt)
  supt_critval <- as.numeric(stats::quantile(apply(t_stat, 1, max), 1 - alpha))

  # Random draws for POSI critical value (post-periods only)
  rd <- matrix(stats::rnorm(kk * p), nrow = kk)
  
  best_bic <- Inf
  best_fit <- NULL
  best_df <- NA
  best_J <- NULL
  best_obj <- NULL
  best_pval <- NULL
  restr_class <- "polynomial"
  
  Xtmp <- matrix(1, nrow = p, ncol = 1)
  res <- MDproj2(estimates, var, Xtmp)
  bic <- res$bic + log(p)
  best_fit <- res
  best_df <- ncol(Xtmp)
  best_obj <- res$bic
  best_J <- res$J
  best_pval <- res$model_fit_pval
  best_bic <- bic
  mbhsf <- apply(abs(rd %*% best_J), 1, max)

  # Limit polynomial degree based on number of observations
  max_degree <- min(3, p - 1)
  for (degree in 1:max_degree) {
    Xtmp <- cbind(Xtmp, (1:p)^degree)
    # Skip if we have more columns than rows (would be singular)
    if (ncol(Xtmp) > p) {
      break
    }
    res <- MDproj2(estimates, var, Xtmp)
    bic <- res$bic + log(p) * ncol(Xtmp)
    if (bic < best_bic) {
      best_fit <- res
      best_df <- ncol(Xtmp)
      best_obj <- res$bic
      best_J <- res$J
      best_pval <- res$model_fit_pval
      restr_class <- "polynomial"
      best_bic <- bic
    }
    mbhsf <- pmax(mbhsf, apply(abs(rd %*% res$J), 1, max))
  }

  bic_unrestricted <- log(p) * p
  if (bic_unrestricted < best_bic) {
    estimates_proj <- estimates
    var_proj <- var
    obj <- 0
    model_fit_pval <- 1
    stdVd <- sqrt(diag(var))
    cVdb <- diag(1 / stdVd) %*% var %*% diag(1 / stdVd)
    J <- tryCatch(chol(cVdb), error = function(e) diag(p))
    best_fit <- list(estimates_proj = estimates_proj, var_proj = var_proj, obj = obj, J = J, model_fit_pval = model_fit_pval)
    best_df <- p
    restr_class <- "unrestricted"
    best_bic <- bic_unrestricted
  }
  mbhsf <- pmax(mbhsf, apply(abs(rd %*% best_fit$J), 1, max))
  
  bestK <- NA
  bestlam1 <- NA
  bestlam2 <- NA
  
  if (p >= 6) {
    target_df <- p - 1
    lam_bounds <- find_lam_bounds(p, var, 4)
    loglam1_range <- lam_bounds$lam1_range
    loglam2_range <- lam_bounds$lam2_range
    
    # Set up parallel processing
    use_parallel <- FALSE
    pb_cluster <- NULL  # Initialize cluster variable
    if (parallel) {
      if (!requireNamespace("foreach", quietly = TRUE)) {
        warning("Package 'foreach' is required for parallel processing but not installed. Falling back to sequential processing.")
      } else if (!requireNamespace("doParallel", quietly = TRUE)) {
        warning("Package 'doParallel' is required for parallel processing but not installed. Falling back to sequential processing.")
      } else {
        use_parallel <- TRUE
        # Determine number of cores to use
        if (is.null(n_cores)) {
          # Default: leave one core free
          num_cores <- max(1, parallel::detectCores() - 1)
        } else {
          # Validate n_cores
          if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1 || n_cores != floor(n_cores)) {
            stop("n_cores must be a positive integer")
          }
          num_cores <- n_cores
        }

        `%dopar%` <- foreach::`%dopar%`

        tryCatch({
          # Clean up any existing parallel backends first
          if (foreach::getDoParRegistered()) {
            foreach::registerDoSEQ()  # Reset to sequential backend
          }

          pb_cluster <- parallel::makeCluster(num_cores)
          doParallel::registerDoParallel(pb_cluster)

          # Verify that the backend is registered
          if (!foreach::getDoParRegistered()) {
            stop("Failed to register parallel backend")
          }

          on.exit({
            if (!is.null(pb_cluster)) {
              tryCatch({
                parallel::stopCluster(pb_cluster)
                foreach::registerDoSEQ()  # Reset to sequential after cleanup
              }, error = function(e) NULL)
            }
          }, add = TRUE)
          message(paste("Parallelizing over", num_cores, "cores"))
        }, error = function(e) {
          warning(paste("Failed to set up parallel cluster:", e$message, "\nFalling back to sequential processing."))
          use_parallel <<- FALSE
        })
      }
    }
    
    # Initialize progress reporting
    if (!use_parallel) {
      # Use standard progress bar for sequential processing
      cli::cli_progress_bar(
        format = "Calculating restricted bounds: {cli::pb_bar} K = {K}/{p-1} [{cli::pb_percent}]",
        total = p-1,
        clear = FALSE
      )
    } else {
      message("Processing in parallel...")
    }
    
    # Initialize results variable
    results <- NULL

    # Process all K values (either in parallel or sequentially)
    if (use_parallel) {
      # Parallel processing using foreach
      tryCatch({
        results <- foreach::foreach(
          K = 1:(p-1),
          .packages = c("Matrix", "dplyr"),
          .export = c("setup_grid", "MDprojl2tf", "my_df", "diff_df", "process_K"),
          .errorhandling = "pass"
        ) %dopar% {
          # Add some progress indication even in parallel mode - this isn't working
          if (K %% 5 == 0 || K == 1 || K == (p-1)) {
            message(sprintf("Processing K = %d of %d\n", K, p-1))
          }
          process_K(K, estimates, var, loglam1_range, loglam2_range, target_df, p, rd)
        }

        # Check for errors in parallel execution
        error_indices <- which(sapply(results, inherits, "error"))
        if (length(error_indices) > 0) {
          warning(paste("Errors occurred during parallel processing for K values:",
                        paste(error_indices, collapse = ", ")))
          # Remove error results
          results <- results[!sapply(results, inherits, "error")]
          if (length(results) == 0) {
            stop("All parallel computations failed. Check the error messages above.")
          }
        }
      }, error = function(e) {
        warning(paste("Error in parallel execution:", e$message, "\nFalling back to sequential processing."))
        # Clean up cluster if it exists
        if (exists("pb_cluster") && !is.null(pb_cluster)) {
          tryCatch(parallel::stopCluster(pb_cluster), error = function(e) NULL)
        }
        use_parallel <<- FALSE
        results <<- NULL
      })
    }
    
    # Process sequentially if parallel is not available or failed
    if (!use_parallel || is.null(results)) {
      # Initialize progress bar if not already initialized (e.g., after fallback from parallel)
      # Check if we need to initialize by seeing if we're falling back from parallel
      if (is.null(results)) {
        cli::cli_progress_bar(
          format = "Calculating restricted bounds: {cli::pb_bar} K = {K}/{p-1} [{cli::pb_percent}]",
          total = p-1,
          clear = FALSE
        )
      }
      # Sequential processing with progress updates
      results <- list()
      for (K in 1:(p-1)) {
        cli::cli_progress_update(set = K)
        results[[K]] <- process_K(K, estimates, var, loglam1_range, loglam2_range, target_df, p, rd)
      }
      # Complete the progress bar
      cli::cli_progress_done()
    }

    # Combine results from all K values
    for (result in results) {
      if (result$best_bic < best_bic) {
        bestlam1 <- result$bestlam1
        bestlam2 <- result$bestlam2
        bestK <- result$K
        best_fit <- result$best_fit
        best_df <- result$best_df
        best_obj <- result$best_obj
        best_J <- result$best_J
        best_pval <- result$best_pval
        restr_class <- "M"
        best_bic <- result$best_bic
      }
      mbhsf <- pmax(mbhsf, result$mbhsf)
    }
  }
  
  suptb <- as.numeric(stats::quantile(mbhsf, 1 - alpha))

  # Restricted bounds (post-periods only)
  restricted_LB <- best_fit$estimates_proj - suptb * sqrt(diag(best_fit$var_proj))
  restricted_UB <- best_fit$estimates_proj + suptb * sqrt(diag(best_fit$var_proj))

  # Wald test for pre-trends (H0: no pre-treatment effects)
  Wpre <- NULL
  if (preperiods > 0) {
    Wpre_stat <- as.numeric(t(estimates_pre) %*% solve(var_pre) %*% estimates_pre)
    Wpre_pval <- 1 - stats::pchisq(Wpre_stat, df = preperiods)
    Wpre <- c(statistic = Wpre_stat, pvalue = Wpre_pval)
  }

  # Wald test for no treatment effect (H0: no post-treatment effects)
  Wpost_stat <- as.numeric(t(estimates) %*% solve(var) %*% estimates)
  Wpost_pval <- 1 - stats::pchisq(Wpost_stat, df = p)
  Wpost <- c(statistic = Wpost_stat, pvalue = Wpost_pval)

  # Create data frame with event time
  # Pre-periods: -preperiods to -1, Post-periods: 1 to p
  if (preperiods > 0) {
    # Post-period bounds
    bounds_df_post <- data.frame(
      horizon = 1:p,
      unrestr_est = estimates,
      restr_est = best_fit$estimates_proj,
      lower = restricted_LB,
      upper = restricted_UB
    )
    # Pre-period bounds (no surrogate for pre-periods)
    bounds_df_pre <- data.frame(
      horizon = -preperiods:-1,
      unrestr_est = estimates_pre,
      restr_est = NA_real_,
      lower = NA_real_,
      upper = NA_real_
    )
    bounds_df <- rbind(bounds_df_pre, bounds_df_post)
  } else {
    bounds_df <- data.frame(
      horizon = 1:p,
      unrestr_est = estimates,
      restr_est = best_fit$estimates_proj,
      lower = restricted_LB,
      upper = restricted_UB
    )
  }

  # Create result list with restricted bounds
  result <- list(
    restricted_bounds = bounds_df
  )

  # Add Wald tests
  if (!is.null(Wpre)) {
    result$Wpre <- Wpre
  }

  result$Wpost <- Wpost

  # Remove J from best_fit before storing in metadata
  best_fit_clean <- best_fit
  best_fit_clean$J <- NULL

  # Add metadata
  result$metadata <- list(
    alpha = alpha,
    preperiods = preperiods,
    supt_critval = supt_critval,
    suptb = suptb,
    df = best_df,
    K = bestK,
    lambda1 = bestlam1,
    lambda2 = bestlam2,
    restr_class = restr_class,
    individual_upper = restricted_UB,
    individual_lower = restricted_LB,
    best_fit_model = best_fit_clean
  )

  class(result) <- c("restricted_bounds", "plausible_bounds_result")
  return(result)
}


process_K <- function(K, estimates, var, loglam1_range, loglam2_range, target_df, p, rd) {
  # Initialize results for this K
  K_best_bic <- Inf
  K_bestlam1 <- NA
  K_bestlam2 <- NA
  K_best_fit <- NULL
  K_best_df <- NA
  K_best_obj <- NULL
  K_best_J <- NULL
  K_best_pval <- NULL
  K_mbhsf <- numeric(nrow(rd))
  
  loglambda_grid <- setup_grid(20, loglam1_range, loglam2_range, K, var, 4, target_df)
  
  for (jj in 1:nrow(loglambda_grid)) {
    res <- MDprojl2tf(estimates, var, exp(loglambda_grid[jj, 1]), exp(loglambda_grid[jj, 2]), K)
    bic <- res$bic + log(p) * res$df

    # Update best results for this K
    if (bic < K_best_bic) {
      K_bestlam1 <- exp(loglambda_grid[jj, 1])
      K_bestlam2 <- exp(loglambda_grid[jj, 2])
      K_best_fit <- res
      K_best_df <- res$df
      K_best_obj <- res$bic
      K_best_J <- res$J
      K_best_pval <- res$model_fit_pval
      K_best_bic <- bic
    }

    # Update mbhsf for this K
    K_mbhsf <- pmax(K_mbhsf, apply(abs(rd %*% res$J), 1, max))
  }
  
  # Return results for this K
  return(list(
    K = K,
    best_bic = K_best_bic,
    bestlam1 = K_bestlam1,
    bestlam2 = K_bestlam2,
    best_fit = K_best_fit,
    best_df = K_best_df,
    best_obj = K_best_obj,
    best_J = K_best_J,
    best_pval = K_best_pval,
    mbhsf = K_mbhsf
  ))
}


MDproj2 <- function(delta, V, X) {
  p <- length(delta)
  Vb <- solve(t(X) %*% solve(V) %*% X)
  beta <- Vb %*% (t(X) %*% solve(V, delta))
  estimates_proj <- X %*% beta
  var_proj <- X %*% Vb %*% t(X)
  stdVd <- sqrt(diag(var_proj))
  cVdb <- diag(1 / stdVd) %*% var_proj %*% diag(1 / stdVd)

  eig <- eigen(cVdb)
  eigVec <- eig$vectors
  D <- eig$values
  J <- eigVec %*% diag(sqrt(D * (Re(D) > 1e-12))) %*% t(eigVec)

  bic <- t(delta - estimates_proj) %*% solve(V, delta - estimates_proj)
  model_fit_pval <- 1 - stats::pchisq(bic, p - length(beta))

  return(list(estimates_proj = estimates_proj, var_proj = var_proj, bic = bic, J = J, model_fit_pval = model_fit_pval))
}

find_lam_bounds <- function(p, V, target_df = 4) {
  lim <- 10
  V <- V / mean(diag(V))
  
  f <- function(lam2) {
    diff_df(-lim, lam2, 1, V, target_df)
  }
  
  result <- stats::optim(par = 1, fn = f, method = "Brent", lower = -50, upper = 40)
                         
  lam2_upper <- result$par

  lam1_range <- matrix(0, nrow = p-1, ncol = 2)
  lam2_range <- matrix(0, nrow = p-1, ncol = 2)
  
  lam1_range[, 1] <- -lim
  lam1_range[, 2] <- lim
  lam2_range[, 1] <- -lim
  lam2_range[, 2] <- lam2_upper
  
  return(list(lam1_range = lam1_range, lam2_range = lam2_range))
}

setup_grid <- function(n_grid, loglam1_range, loglam2_range, K, V, lb, ub) {
  grid_ll2 <- seq(loglam2_range[K, 1], loglam2_range[K, 2], length.out = n_grid)
  grid_ll1 <- seq(loglam1_range[K, 1], loglam1_range[K, 2], length.out = n_grid)
  full_grid <- expand.grid(grid_ll1, grid_ll2) %>%
    dplyr::arrange(Var1)
  
  # Normalize V once outside the loop
  scaled_V <- V / mean(diag(V))
  
  grid_df <- numeric(n_grid^2)
  
  for (sim in 1:(n_grid^2)) {
    grid_df[sim] <- my_df(full_grid[sim, 1], full_grid[sim, 2], K, scaled_V)
  }
  
  legit <- (grid_df <= ub) & (grid_df >= lb)
  
  if (mean(legit) < 0.1) {
    warning('Few gridpoints in range, results may be suboptimal')
  }
  
  loglambda_grid <- full_grid[legit, ]
  return(loglambda_grid)
}

MDprojl2tf <- function(delta, V, lambda1, lambda2, K) {

  p <- length(delta)

  D1 <- matrix(0, nrow = p, ncol = p)
  D1[cbind(2:p, 1:(p-1))] <- 1
  D1 <- D1 - diag(p)
  D1 <- D1[-1, ]

  D2 <- D1[-1, -1]
  D3 <- D2[-1, -1] %*% D2 %*% D1

  diag_mean <- mean(diag(V))
  scaledV <- V / diag_mean
  iV <- solve(V)
  scalediV <- solve(scaledV)

  # Calculate vD1 and related matrices
  vD1 <- D1 %*% scaledV %*% t(D1)
  zeros_block <- matrix(0, nrow = K-1, ncol = K-1)
  partial_vD1 <- vD1[K:nrow(vD1), K:ncol(vD1)] %>% as.matrix()
  vD1_block <- partial_vD1 / mean(diag(partial_vD1))
  W1 <- as.matrix(Matrix::bdiag(zeros_block, vD1_block))

  # Calculate vD3 and W3
  vD3 <- D3 %*% scaledV %*% t(D3)
  W3 <- vD3 / mean(diag(vD3))

  # Solve system and calculate results
  den <- solve(scalediV + lambda1 * t(D1) %*% diag(diag(W1)) %*% D1 +
                 lambda2 * t(D3) %*% diag(diag(W3)) %*% D3)

  df <- sum(diag(den %*% scalediV))

  estimates_proj <- den %*% (scalediV %*% delta)

  var_proj <- den %*% scalediV %*% t(den)
  var_proj <- var_proj * diag_mean

  stdVd <- sqrt(diag(var_proj))
  cVdb <- diag(1/stdVd) %*% var_proj %*% diag(1/stdVd)

  eigen_result <- eigen(cVdb)
  eigVec <- eigen_result$vectors
  eigen_values <- eigen_result$values
  D_eigen <- diag(eigen_values)

  if (any(Im(eigen_values) > 0)) {
    warning("Complex eigenvalues detected in covariance matrix, results may be unreliable")
  }

  # Create J matrix more efficiently
  sqrt_D <- sqrt(pmax(D_eigen, 1e-12))
  sqrt_D[D_eigen <= 1e-12] <- 0
  J <- eigVec %*% sqrt_D %*% t(eigVec)

  residual <- delta - estimates_proj
  bic <- as.numeric(t(residual) %*% iV %*% residual)

  model_fit_pval <- 1 - stats::pchisq(bic, df = p - df)

  return(list(
    estimates_proj = as.vector(estimates_proj),
    var_proj = var_proj,
    bic = bic,
    J = J,
    model_fit_pval = model_fit_pval,
    df = df
  ))
}

diff_df <- function(loglam1, loglam2, K, V, targetdf) {
  p <- nrow(V)
  
  # Create D matrices more efficiently
  D1 <- matrix(0, nrow = p, ncol = p)
  D1[cbind(2:p, 1:(p-1))] <- 1
  D1 <- D1 - diag(p)
  D1 <- D1[-1, ]
  
  D2 <- D1[-1, -1]
  D3 <- D2[-1, -1] %*% D2 %*% D1
  
  # Calculate vD1 and related matrices
  vD1 <- D1 %*% V %*% t(D1)
  zeros_block <- matrix(0, nrow = K-1, ncol = K-1)
  vD1_subset <- vD1[K:nrow(vD1), K:ncol(vD1)]
  diag_mean <- mean(diag(vD1_subset))
  vD1_normalized <- vD1_subset / diag_mean
  
  W1 <- matrix(0, nrow = nrow(vD1), ncol = ncol(vD1))
  W1[K:nrow(W1), K:ncol(W1)] <- vD1_normalized
  
  # Calculate vD3 and W2
  vD3 <- D3 %*% V %*% t(D3)
  W2 <- vD3 / mean(diag(vD3))
  
  # Pre-compute exponentials
  exp_loglam1 <- exp(loglam1)
  exp_loglam2 <- exp(loglam2)
  
  # Solve system and calculate df
  iV <- solve(V)
  den <- solve(iV + exp_loglam1 * t(D1) %*% diag(diag(W1)) %*% D1 +
                 exp_loglam2 * t(D3) %*% diag(diag(W2)) %*% D3)
  
  df <- sum(diag(den %*% iV))
  
  # Calculate squared difference from target
  f <- (df - targetdf)^2
  
  return(f)
}

my_df <- function(loglam1, loglam2, K, V) {
  
  p <- nrow(V)
  
  # Create D matrices more efficiently
  D1 <- matrix(0, p, p)
  diag(D1[-1, ]) <- 1
  D1 <- D1 - diag(1, p)
  D1 <- D1[-1, ]
  
  D2 <- D1[-1, -1]
  D3 <- D2[-1, -1] %*% D2 %*% D1
  
  # Calculate vD1 and related matrices
  vD1 <- D1 %*% V %*% t(D1)
  zeros_block <- matrix(0, K-1, K-1)
  vD1_subset <- vD1[K:nrow(vD1), K:ncol(vD1)] %>% as.matrix()
  diag_mean <- mean(diag(vD1_subset))
  normalized_vD1 <- vD1_subset / diag_mean
  W1 <- Matrix::bdiag(zeros_block, normalized_vD1)
  W1 <- as.matrix(W1)
  
  # Calculate vD3 and W2
  vD3 <- D3 %*% V %*% t(D3)
  W2 <- vD3 / mean(diag(vD3))
  
  # Pre-compute exponentials
  exp_loglam1 <- exp(loglam1)
  exp_loglam2 <- exp(loglam2)
  
  # Solve system and calculate df
  iV <- solve(V)
  den <- solve(iV + exp_loglam1 * t(D1) %*% diag(diag(W1)) %*% D1 +
                 exp_loglam2 * t(D3) %*% diag(diag(W2)) %*% D3)
  
  df_result <- sum(diag(den %*% iV))
  
  return(df_result)
}
