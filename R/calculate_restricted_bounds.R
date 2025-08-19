#' Calculate Restricted Bounds
#'
#' This function calculates restricted bounds for a vector of estimates using model selection.
#' It can calculate both pointwise and simultaneous (sup-t) bounds.
#' 
#' 
#' @param estimates A vector of point estimates
#' @param var The variance-covariance matrix of the estimates
#' @param alpha Significance level (default: 0.05)
#' @param include_pointwise Whether to include pointwise bounds (default: TRUE)
#' @param include_supt Whether to include sup-t bounds (default: TRUE)
#' @param parallel Whether to use parallel processing (default: FALSE)
#'
#' @return A list containing:
#'   \item{bounds}{A data frame with columns for horizon, coefficients, surrogate values, and bounds}
#'   \item{type}{The type of bounds ("restricted")}
#'   \item{metadata}{A list with metadata about the calculation}
#'
#' @examples
#' # Example with constant estimates and IID errors
#' data(estimates_constant_iid)
#' data(var_constant_iid)
#' restr_bounds <- calculate_restricted_bounds(estimates_constant_iid, var_constant_iid)
#'
#' @export

calculate_restricted_bounds <- function(estimates, var, alpha = 0.05,
                                       include_pointwise = TRUE, include_supt = TRUE,
                                       parallel = FALSE) {
  if (!is.numeric(estimates) || !is.vector(estimates)) {
    stop("estimates must be a numeric vector")
  }
  if (!is.matrix(var) || nrow(var) != length(estimates) || ncol(var) != length(estimates)) {
    stop("var must be a square matrix with dimensions matching the length of estimates")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number between 0 and 1")
  }
  
  p <- length(estimates)
  
  eigs <- eigen(var, symmetric = TRUE)
  if (min(eigs$values) < 0) {
    var <- eigs$vectors %*% (diag(eigs$values) + (abs(min(eigs$values)) + 1e-8) * diag(length(eigs$values))) %*% t(eigs$vectors)
  }
  
  stdv <- sqrt(diag(var))
  cV <- diag(1 / stdv) %*% var %*% diag(1 / stdv)
  
  set.seed(42)
  kk <- 10000
  rd <- matrix(stats::rnorm(kk * p), nrow = kk)
  
  Corrmat <- cV
  eig_c <- eigen(Corrmat, symmetric = TRUE)
  Corrmat_sqrt <- eig_c$vectors %*% diag(sqrt(pmax(eig_c$values, 0))) %*% t(eig_c$vectors)
  t_stat <- abs(rd %*% Corrmat_sqrt)
  supt_critval <- as.numeric(stats::quantile(apply(t_stat, 1, max), 1 - alpha))
  
  best_bic <- Inf
  best_fit <- NULL
  best_df <- NA
  best_J <- NULL
  best_obj <- NULL
  best_pval <- NULL
  surrogate_class <- "polynomial"
  
  Xtmp <- matrix(1, nrow = p, ncol = 1)
  res <- MDproj2(estimates, var, Xtmp) 
  bic <- res$MD + log(p) 
  best_fit <- res
  best_df <- ncol(Xtmp)
  best_obj <- res$MD
  best_J <- res$J
  best_pval <- res$pval
  best_bic <- bic
  mbhsf <- apply(abs(rd %*% best_J), 1, max)
  
  for (degree in 1:3) {
    Xtmp <- cbind(Xtmp, (1:p)^degree)
    res <- MDproj2(estimates, var, Xtmp)
    bic <- res$MD + log(p) * ncol(Xtmp)
    if (bic < best_bic) {
      best_fit <- res
      best_df <- ncol(Xtmp)
      best_obj <- res$MD
      best_J <- res$J
      best_pval <- res$pval
      surrogate_class <- "polynomial"
      best_bic <- bic
    }
    mbhsf <- pmax(mbhsf, apply(abs(rd %*% res$J), 1, max))
  }
  
  bic_unrestricted <- log(p) * p
  if (bic_unrestricted < best_bic) {
    d <- estimates
    Vd <- var
    obj <- 0
    pval <- 1
    stdVd <- sqrt(diag(var))
    cVdb <- diag(1 / stdVd) %*% var %*% diag(1 / stdVd)
    J <- tryCatch(chol(cVdb), error = function(e) diag(p))
    best_fit <- list(d = d, Vd = Vd, obj = obj, J = J, pval = pval)
    best_df <- p
    surrogate_class <- "unrestricted"
    best_bic <- bic_unrestricted
  }
  mbhsf <- pmax(mbhsf, apply(abs(rd %*% best_fit$J), 1, max))
  
  bestK <- NA
  bestlam1 <- NA
  bestlam2 <- NA
  
  if (p >= 6) {
    target_df <- p - 1
    lam_bounds <- find_lam_bounds(p, var, target_df)
    loglam1_range <- lam_bounds$lam1_range
    loglam2_range <- lam_bounds$lam2_range
    
    # Set up parallel processing
    use_parallel <- FALSE
    if (parallel) {
      if (!requireNamespace("foreach", quietly = TRUE)) {
        warning("Package 'foreach' is required for parallel processing but not installed. Falling back to sequential processing.")
      } else if (!requireNamespace("doParallel", quietly = TRUE)) {
        warning("Package 'doParallel' is required for parallel processing but not installed. Falling back to sequential processing.")
      } else {
        use_parallel <- TRUE
        # Determine number of cores to use (leave one core free)
        num_cores <- max(1, parallel::detectCores() - 1)
        
        `%dopar%` <- foreach::`%dopar%` 
        # Limit to a reasonable number of cores to avoid system overload
        max_recommended_cores <- 16
        if (num_cores > max_recommended_cores) {
          message(paste("Limiting to", max_recommended_cores, "cores to avoid system overload"))
          num_cores <- max_recommended_cores
        }
        
        tryCatch({
          cl <- parallel::makeCluster(num_cores)
          doParallel::registerDoParallel(cl)
          on.exit(parallel::stopCluster(cl), add = TRUE)
          message(paste("Parallelizing over", num_cores, "cores"))
        }, error = function(e) {
          warning(paste("Failed to set up parallel cluster:", e$message, "\nFalling back to sequential processing."))
          use_parallel <- FALSE
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
      message("Processing in parallel. Progress updates will be limited.")
    }
    
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
        use_parallel <- FALSE
      })
    }
    
    # Process sequentially if parallel is not available or failed
    if (!use_parallel) {
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
        surrogate_class <- "M"
        best_bic <- result$best_bic
      }
      mbhsf <- pmax(mbhsf, result$mbhsf)
    }
  }
  
  suptb <- as.numeric(stats::quantile(mbhsf, 1 - alpha))
  
  restricted_LB <- best_fit$d - suptb * sqrt(diag(best_fit$Vd))
  restricted_UB <- best_fit$d + suptb * sqrt(diag(best_fit$Vd))
  
  # Create data frame with horizon and coefficients
  bounds_df <- data.frame(
    horizon = 1:p,
    coef = estimates,
    surrogate = best_fit$d,
    lower = restricted_LB,
    upper = restricted_UB
  )
  
  # Calculate width
  restricted_width <- mean(bounds_df$upper - bounds_df$lower)
  
  # Create result list with restricted bounds
  result <- list(
    restricted_bounds = bounds_df
  )
  
  # Calculate pointwise bounds if requested
  if (include_pointwise) {
    result$pointwise_bounds <- calculate_pointwise_bounds(estimates, var, alpha)
  }
  
  # Calculate sup-t bounds if requested
  if (include_supt) {
    result$supt_bounds <- calculate_supt_bounds(estimates, var, alpha)
  }
  
  # Add metadata
  result$metadata <- list(
    alpha = alpha,
    suptb = suptb,
    df = best_df,
    K = bestK,
    lambda1 = bestlam1,
    lambda2 = bestlam2,
    surrogate_class = surrogate_class,
    width = restricted_width,
    individual_upper = restricted_UB,
    individual_lower = restricted_LB,
    best_fit_model = best_fit
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
  
  # Setup grid for this K
  loglambda_grid <- setup_grid(20, loglam1_range, loglam2_range, K, var, 4, target_df)
  
  # Process each grid point
  for (jj in 1:nrow(loglambda_grid)) {
    res <- MDprojl2tf(estimates, var, exp(loglambda_grid[jj, 1]), exp(loglambda_grid[jj, 2]), K)
    bic <- res$MD + log(p) * res$df
    
    # Update best results for this K
    if (bic < K_best_bic) {
      K_bestlam1 <- exp(loglambda_grid[jj, 1])
      K_bestlam2 <- exp(loglambda_grid[jj, 2])
      K_best_fit <- res
      K_best_df <- res$df
      K_best_obj <- res$MD
      K_best_J <- res$J
      K_best_pval <- res$pval
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
  d <- X %*% beta                           
  Vd <- X %*% Vb %*% t(X)                   
  stdVd <- sqrt(diag(Vd))
  cVdb <- diag(1 / stdVd) %*% Vd %*% diag(1 / stdVd)
  
  eig <- eigen(cVdb)
  eigVec <- eig$vectors
  D <- eig$values
  J <- eigVec %*% diag(sqrt(D * (Re(D) > 1e-12))) %*% t(eigVec)
  
  MD <- t(delta - d) %*% solve(V, delta - d)
  pv <- 1 - stats::pchisq(MD, p - length(beta))
  
  return(list(d = d, Vd = Vd, MD = MD, J = J, pval = pv))
}

find_lam_bounds <- function(p, V, target_df) {
  lim <- 10
  V <- V / mean(diag(V))
  
  f <- function(lam2) {
    diff_df(-lim, lam2, 1, V, target_df)
  }
  
  result <- stats::optim(par = 1, fn = f, method = "Brent", lower = -50, upper = 50)
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
  
  # Create D matrices more efficiently
  D1 <- matrix(0, nrow = p, ncol = p)
  D1[cbind(2:p, 1:(p-1))] <- 1
  D1 <- D1 - diag(p)
  D1 <- D1[-1, ]
  
  D2 <- D1[-1, -1]
  D3 <- D2[-1, -1] %*% D2 %*% D1
  
  # Pre-compute values used multiple times
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
  
  d <- den %*% (scalediV %*% delta)
  
  Vd <- den %*% scalediV %*% t(den)
  Vd <- Vd * diag_mean
  
  stdVd <- sqrt(diag(Vd))
  cVdb <- diag(1/stdVd) %*% Vd %*% diag(1/stdVd)
  
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
  
  residual <- delta - d
  MD <- as.numeric(t(residual) %*% iV %*% residual)
  
  pv <- 1 - stats::pchisq(MD, df = p - df)
  
  return(list(
    d = as.vector(d),
    Vd = Vd,
    MD = MD,
    J = J,
    pval = pv,
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
