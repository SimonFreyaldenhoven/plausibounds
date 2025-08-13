#' Calculate Restricted Bounds
#'
#' This function calculates restricted bounds for a vector of estimates using model selection.
#' It can calculate both pointwise and simultaneous (sup-t) bounds.
#'
#' @param estimates A vector of point estimates
#' @param var The variance-covariance matrix of the estimates
#' @param alpha Significance level (default: 0.05)
#' @param include_pointwise Whether to include pointwise bounds (default: TRUE)
#' @param include_supt Whether to include sup-t bounds (default: TRUE)
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
#' # Example with wiggly estimates and strong correlation
#' data(estimates_wiggly_strong_corr)
#' data(var_wiggly_strong_corr)
#' restr_bounds_complex <- calculate_restricted_bounds(
#'   estimates_wiggly_strong_corr,
#'   var_wiggly_strong_corr
#' )
#'
#' @export
calculate_restricted_bounds <- function(estimates, var, alpha = 0.05, 
                                       include_pointwise = TRUE, include_supt = TRUE) {
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
  rd <- matrix(rnorm(kk * p), nrow = kk)
  
  Corrmat <- cV
  eig_c <- eigen(Corrmat, symmetric = TRUE)
  Corrmat_sqrt <- eig_c$vectors %*% diag(sqrt(pmax(eig_c$values, 0))) %*% t(eig_c$vectors)
  t_stat <- abs(rd %*% Corrmat_sqrt)
  supt_critval <- as.numeric(quantile(apply(t_stat, 1, max), 1 - alpha))
  
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
    
    for (K in 1:(p-1)) {
      loglambda_grid <- setup_grid(20, loglam1_range, loglam2_range, K, var, 4, target_df)
      for (jj in 1:nrow(loglambda_grid)) {
        res <- MDprojl2tf(estimates, var, exp(loglambda_grid[jj, 1]), exp(loglambda_grid[jj, 2]), K)
        bic <- res$MD + log(p) * res$df
        if (bic < best_bic) {
          bestlam1 <- exp(loglambda_grid[jj, 1])
          bestlam2 <- exp(loglambda_grid[jj, 2])
          bestK <- K
          best_fit <- res
          best_df <- res$df
          best_obj <- res$MD
          best_J <- res$J
          best_pval <- res$pval
          surrogate_class <- "M"
          best_bic <- bic
        }
        mbhsf <- pmax(mbhsf, apply(abs(rd %*% res$J), 1, max))
      }
    }
  }
  
  suptb <- as.numeric(quantile(mbhsf, 1 - alpha))
  
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
  pv <- 1 - pchisq(MD, p - length(beta))
  
  return(list(d = d, Vd = Vd, MD = MD, J = J, pval = pv))
}

find_lam_bounds <- function(p, V, target_df) {
  lim <- 10
  V <- V / mean(diag(V))
  
  f <- function(lam2) {
    diff_df(-lim, lam2, 1, V, target_df)
  }
  
  result <- optim(par = 1, fn = f, method = "Brent", lower = -50, upper = 50)
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
  
  grid_df <- numeric(n_grid^2)
  
  for (sim in 1:(n_grid^2)) {
    grid_df[sim] <- my_df(full_grid[sim, 1], full_grid[sim, 2], K, V / mean(diag(V)))
  }
  
  legit <- (grid_df <= ub) & (grid_df >= lb)
  
  if (mean(legit) < 0.1) {
    cat('Warning: few gridpoints in range\n')
  }
  
  loglambda_grid <- full_grid[legit, ]
  return(loglambda_grid)
}

MDprojl2tf <- function(delta, V, lambda1, lambda2, K) {
  library(Matrix)
  
  p <- length(delta)
  
  D1 <- matrix(0, nrow = p, ncol = p)
  D1[cbind(2:p, 1:(p-1))] <- 1  
  D1 <- D1 - diag(p)  
  D1 <- D1[-1, ]  
  
  D2 <- D1[-1, -1]  
  D3 <- D2[-1, -1]  
  D3 <- D3 %*% D2 %*% D1  
  
  iV <- solve(V)
  scaledV <- V / mean(diag(V))
  
  vD1 <- D1 %*% scaledV %*% t(D1)
  zeros_block <- matrix(0, nrow = K-1, ncol = K-1)
  partial_vD1 <- vD1[K:nrow(vD1), K:ncol(vD1)] %>% as.matrix()
  vD1_block <- partial_vD1 / mean(diag(partial_vD1))
  W1 <- as.matrix(bdiag(zeros_block, vD1_block))
  
  vD3 <- D3 %*% scaledV %*% t(D3)
  W3 <- vD3 / mean(diag(vD3))
  
  scalediV <- solve(scaledV)
  
  den <- solve(scalediV + lambda1 * t(D1) %*% diag(diag(W1)) %*% D1 + 
                 lambda2 * t(D3) %*% diag(diag(W3)) %*% D3)
  
  df <- sum(diag(den %*% scalediV))
  
  d <- den %*% (scalediV %*% delta)
  
  Vd <- den %*% scalediV %*% t(den)
  Vd <- Vd * mean(diag(V))
  
  stdVd <- sqrt(diag(Vd))
  cVdb <- diag(1/stdVd) %*% Vd %*% diag(1/stdVd)
  
  eigen_result <- eigen(cVdb)
  eigVec <- eigen_result$vectors
  D_eigen <- diag(eigen_result$values)
  
  if (any(Im(eigen_result$values) > 0)) {
    cat("hi - perhaps issue with covariance?\n")
  }
  
  sqrt_D <- sqrt(pmax(D_eigen, 1e-12))
  sqrt_D[D_eigen <= 1e-12] <- 0
  J <- eigVec %*% sqrt_D %*% t(eigVec)
  
  residual <- delta - d
  MD <- as.numeric(t(residual) %*% solve(V) %*% residual)
  
  pv <- 1 - pchisq(MD, df = p - df)
  
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
  
  D1 <- matrix(0, nrow = p, ncol = p)
  D1[cbind(2:p, 1:(p-1))] <- 1  
  D1 <- D1 - diag(p)  
  D1 <- D1[-1, ]  
  
  D2 <- D1[-1, -1]  
  D3 <- D2[-1, -1]  
  D3 <- D3 %*% D2 %*% D1  
  
  vD1 <- D1 %*% V %*% t(D1)
  
  zeros_block <- matrix(0, nrow = K-1, ncol = K-1)
  vD1_subset <- vD1[K:nrow(vD1), K:ncol(vD1)]
  vD1_normalized <- vD1_subset / mean(diag(vD1_subset))
  
  W1 <- matrix(0, nrow = nrow(vD1), ncol = ncol(vD1))
  W1[K:nrow(W1), K:ncol(W1)] <- vD1_normalized
  
  vD3 <- D3 %*% V %*% t(D3)
  W2 <- vD3 / mean(diag(vD3))
  
  iV <- solve(V)
  
  den <- solve(iV + exp(loglam1) * t(D1) %*% diag(diag(W1)) %*% D1 + 
                 exp(loglam2) * t(D3) %*% diag(diag(W2)) %*% D3)
  
  df <- sum(diag(den %*% iV))  
  
  f <- (df - targetdf)^2
  
  return(f)
}

my_df <- function(loglam1, loglam2, K, V) {
  library(Matrix)  
  
  p <- nrow(V)
  
  D1 <- matrix(0, p, p)
  diag(D1[-1, ]) <- 1  
  D1 <- D1 - diag(1, p)  
  D1 <- D1[-1, ]  
  
  D2 <- D1[-1, -1]  
  
  D3 <- D2[-1, -1]  
  
  D3 <- D3 %*% D2 %*% D1  
  
  vD1 <- D1 %*% V %*% t(D1)
  
  zeros_block <- matrix(0, K-1, K-1)
  vD1_subset <- vD1[K:nrow(vD1), K:ncol(vD1)] %>% as.matrix()
  normalized_vD1 <- vD1_subset / mean(diag(vD1_subset))
  W1 <- bdiag(zeros_block, normalized_vD1)
  W1 <- as.matrix(W1)  
  
  vD3 <- D3 %*% V %*% t(D3)
  W2 <- vD3 / mean(diag(vD3))
  
  iV <- solve(V)
  
  den <- solve(iV + exp(loglam1) * t(D1) %*% diag(diag(W1)) %*% D1 + 
                 exp(loglam2) * t(D3) %*% diag(diag(W2)) %*% D3)
  
  df_result <- sum(diag(den %*% iV))
  
  return(df_result)
}

