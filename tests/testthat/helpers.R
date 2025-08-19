
full_eventplot_l2tf <- function(delta, var, truth = NULL, nsim = 2000) {
  # Ensure positive-definite variance matrix (MATLAB: eig, A*B*A')
  eigs <- eigen(var, symmetric = TRUE)
  if (min(eigs$values) < 0) {
    var <- eigs$vectors %*% (diag(eigs$values) + (abs(min(eigs$values)) + 1e-8) * diag(length(eigs$values))) %*% t(eigs$vectors)
  }

  p <- length(delta)
  # --- SUP-t (simultaneous) bands ---
  supt_bands <- bands_plugin(delta, var, p, nsim = nsim, level = 0.95) 
  supt_LB <- supt_bands$LB
  supt_UB <- supt_bands$UB
  supt <- supt_bands$sup_t

  # --- Standardized covariance and correlation matrix (MATLAB: diag, sqrt, diag(1./stdV)*V*diag(1./stdV)) ---
  stdv <- sqrt(diag(var))
  cV <- diag(1 / stdv) %*% var %*% diag(1 / stdv)

  # --- Critical value simulation (MATLAB: randn) ---
  set.seed(42)
  kk <- 10000
  rd <- matrix(rnorm(kk * p), nrow = kk)
  
  # --- sup-t (MATLAB: max, prctile) ---
  Corrmat <- cV
  eig_c <- eigen(Corrmat, symmetric = TRUE)
  Corrmat_sqrt <- eig_c$vectors %*% diag(sqrt(pmax(eig_c$values, 0))) %*% t(eig_c$vectors)
  t_stat <- abs(rd %*% Corrmat_sqrt)
  supt_critval <- as.numeric(quantile(apply(t_stat, 1, max), 0.95)) # This may be the same as supt above

  # --- Model selection among polynomials (MATLAB: while, Xtmp, MDproj2) ---
  best_bic <- Inf
  best_fit <- NULL
  best_df <- NA
  best_J <- NULL
  best_obj <- NULL
  best_pval <- NULL
  best_condtruth <- NULL
  surrogate_class <- "polynomial"
  
  Xtmp <- matrix(1, nrow = p, ncol = 1)
  res <- MDproj2(delta, var, Xtmp) 
  bic <- res$MD + log(p) 
  best_fit <- res
  best_df <- ncol(Xtmp)
  best_obj <- res$obj
  best_J <- res$J
  best_pval <- res$pval
  best_bic <- bic
  mbhsf <- apply(abs(rd %*% best_J), 1, max)
  
  for (degree in 1:3) {
    Xtmp <- cbind(Xtmp, (1:p)^degree)
    res <- MDproj2(delta, var, Xtmp) # custom implementation
    bic <- res$MD + log(p) * ncol(Xtmp)
    if (!is.null(truth)) condtruth <- MDproj2(truth, var, Xtmp)$d else condtruth <- NULL
    if (bic < best_bic) {
      best_fit <- res
      best_df <- ncol(Xtmp)
      best_obj <- res$obj
      best_J <- res$J
      best_pval <- res$pval
      surrogate_class <- "polynomial"
      best_bic <- bic
      best_condtruth <- condtruth
    }
    mbhsf <- pmax(mbhsf, apply(abs(rd %*% res$J), 1, max))
  }

  # --- Unrestricted estimate (MATLAB: diag, sqrtm) ---
  bic_unrestricted <- log(p) * p
  if (bic_unrestricted < best_bic) {
    d <- delta
    Vd <- var
    obj <- 0
    pval <- 1
    stdVd <- sqrt(diag(var))
    cVdb <- diag(1 / stdVd) %*% var %*% diag(1 / stdVd)
    J <- tryCatch(chol(cVdb), error = function(e) diag(p))
    if (!is.null(truth)) condtruth <- truth else condtruth <- NULL
    best_fit <- list(d = d, Vd = Vd, obj = obj, J = J, pval = pval)
    best_df <- p
    surrogate_class <- "unrestricted"
    best_bic <- bic_unrestricted
    best_condtruth <- condtruth
  }
  mbhsf <- pmax(mbhsf, apply(abs(rd %*% best_fit$J), 1, max)) # do not understand this part

  # --- More flexible smoothers (trend filtering, MATLAB: MDprojl2tf, skipped unless p >= 6) ---
  if (p >= 6) {
    target_df <- p - 1
    lam_bounds <- find_lam_bounds(p, var, target_df) # custom implementation
    loglam1_range <- lam_bounds$lam1_range
    loglam2_range <- lam_bounds$lam2_range
    bestK <- -1
    bestlam1 <- -1
    bestlam2 <- -1
    for (K in 1:(p-1)) {
      print(K)
      loglambda_grid <- setup_grid(20, loglam1_range, loglam2_range, K, var, 4, target_df) # custom implementation
      for (jj in 1:nrow(loglambda_grid)) {
        res <- MDprojl2tf(delta, var, exp(loglambda_grid[jj, 1]), exp(loglambda_grid[jj, 2]), K) # custom implementation
        if (!is.null(truth)) condtruth <- MDprojl2tf(truth, var, exp(loglambda_grid[jj, 1]), exp(loglambda_grid[jj, 2]), K)$d else condtruth <- NULL
        bic <- res$MD + log(p) * res$df
        if (bic < best_bic) {
          bestlam1 <- exp(loglambda_grid[jj, 1])
          bestlam2 <- exp(loglambda_grid[jj, 2])
          bestK <- K
          best_fit <- res
          best_df <- res$df
          best_obj <- res$obj
          best_J <- res$J
          best_pval <- res$pval
          surrogate_class <- "M"
          best_bic <- bic
          best_condtruth <- condtruth
        }
        mbhsf <- pmax(mbhsf, apply(abs(rd %*% res$J), 1, max))
      }
    }
  }

  # --- POSI critical value (MATLAB: prctile) ---
  suptb <- as.numeric(quantile(mbhsf, 0.95))

  # --- Wald bounds on cumulative effect (MATLAB: Wald_bounds) ---
  wald_bounds <- Wald_bounds(delta, var, alpha = 0.05) # custom implementation
  lb <- sum(wald_bounds$LB)
  ub <- sum(wald_bounds$UB)

  # --- Diagnostics (MATLAB: diagnostics struct) ---
   diagnostics <- list()
  #   diagnostics$inpw <- all(truth < delta + 1.96 * sqrt(diag(var)) & truth > delta - 1.96 * sqrt(diag(var)))
  #   diagnostics$inwald <- sum(truth) < sum(wald_bounds$UB) & sum(truth) > sum(wald_bounds$LB)
  #   diagnostics$insupt <- all(truth < delta + supt * sqrt(diag(var)) & truth > delta - supt * sqrt(diag(var)))
  #   diagnostics$inrestricted <- all(truth < best_fit$d + suptb * sqrt(diag(best_fit$Vd)) & truth > best_fit$d - suptb * sqrt(diag(best_fit$Vd)))
  #   diagnostics$incond <- if (!is.null(best_condtruth)) all(best_condtruth < best_fit$d + suptb * sqrt(diag(best_fit$Vd)) & best_condtruth > best_fit$d - suptb * sqrt(diag(best_fit$Vd))) else NA
  #   diagnostics$mse_pw <- mean((delta - truth)^2)
  #   diagnostics$mse_smooth <- mean((best_fit$d - truth)^2)
     diagnostics$pw_width <- mean(2 * 1.96 * sqrt(diag(var)))
     diagnostics$supt_width <- mean(2 * supt_critval * sqrt(diag(var)))
     diagnostics$restricted_width <- mean(2 * suptb * sqrt(diag(best_fit$Vd)))
     diagnostics$suptb <- suptb
  #   diagnostics$wald_width <- (sum(wald_bounds$UB) - sum(wald_bounds$LB)) / p
     diagnostics$df <- best_df
     diagnostics$K <- if (exists("bestK")) bestK else NA
     diagnostics$lam1 <- if (exists("bestlam1")) bestlam1 else NA
     diagnostics$lam2 <- if (exists("bestlam2")) bestlam2 else NA
     diagnostics$suptbands <- cbind(delta + supt_critval * sqrt(diag(var)), delta - supt_critval * sqrt(diag(var)))
  #   diagnostics$surrogate <- best_condtruth
  #   diagnostics$surrogate_class <- surrogate_class
  

  # Calculate restricted bounds
  restricted_LB <- best_fit$d - suptb * sqrt(diag(best_fit$Vd))
  restricted_UB <- best_fit$d + suptb * sqrt(diag(best_fit$Vd))
  
  list(
    diagnostics = diagnostics,
    lb = lb,
    ub = ub,
    supt_LB = supt_LB,
    supt_UB = supt_UB,
    supt_critval = supt_critval,
    suptb = suptb,
    best_fit = best_fit,
    surrogate = best_fit$d,
    restricted_LB = restricted_LB,
    restricted_UB = restricted_UB,
    mbhsf = mbhsf
  )
}

# --- Helper function stubs (to be implemented) ---

bands_plugin <- function(delta, var, p, nsim, level) {
  # Simulate from N(0, var), calculate sup-t bands, return LB and UB
  rd <- MASS::mvrnorm(nsim, mu = rep(0, p), Sigma = var)
  std_devs <- sqrt(diag(var))
  rd_standardized <- rd / matrix(std_devs, nrow = nsim, ncol = length(std_devs), byrow = TRUE)
  sup_t <- as.numeric(quantile(apply(abs(rd_standardized), 1, max), level))
 # sup_t <- as.numeric(quantile(apply(abs(rd), 1, max), level))
  list(
    LB = delta - sup_t * sqrt(diag(var)),
    UB = delta + sup_t * sqrt(diag(var)),
    sup_t = sup_t
  )
}

# Polynomial projection taking as input beta_hat, V_beta, design matrix defining degrees of freedom
MDproj2 <- function(delta, V, X) {
  p <- length(delta)                       
  Vb <- solve(t(X) %*% solve(V) %*% X)      # Covariance of projected coefficients 
  beta <- Vb %*% (t(X) %*% solve(V, delta)) # Projected coefficients
  d <- X %*% beta                           # Projected restricted estimates
  Vd <- X %*% Vb %*% t(X)                   # Covariance of projected estimates (corresponds to V_M from eq 4)
  stdVd <- sqrt(diag(Vd))
  cVdb <- diag(1 / stdVd) %*% Vd %*% diag(1 / stdVd)
  
  eig <- eigen(cVdb)
  eigVec <- eig$vectors
  D <- eig$values
  J <- eigVec %*% diag(sqrt(D * (Re(D) > 1e-12))) %*% t(eigVec)
  
  MD <- t(delta - d) %*% solve(V, delta - d)
  pv <- 1 - pchisq(MD, p - length(beta))
  
  return(list(d = d, Vd = Vd, MD = MD, J = J, pv = pv))
}


find_lam_bounds <- function(p, V, target_df) {
  # with no penalty on third diff (lam2=0), p>df>K+1.
  # with no penalty on first diff (lam1=0), p>df>3
  
  lim <- 10
  V <- V / mean(diag(V))
  
  # Upper bound on lambda2 does not depend on K for the reason above
  f <- function(lam2) {
    diff_df(-lim, lam2, 1, V, target_df)
  }
  
  # Using optim with Nelder-Mead method (equivalent to fminsearch)
  result <- optim(par = 1, fn = f, method = "Brent", lower = -50, upper = 50)
  lam2_upper <- result$par
  
  # cur_df = df(-lim, lam2_upper, 1, V)  # check: should be 4
  
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
    arrange(Var1)
  
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



Wald_bounds <- function(dhat, Vhat, alpha, df = 1) {
  h <- length(dhat)
  critval <- qchisq(1 - alpha, df)
  
  lambda1 <- sqrt(sum(Vhat) / (4 * critval))
  lambda2 <- -sqrt(sum(Vhat) / (4 * critval))
  
  UB <- t(dhat) + (1 / (2 * lambda1)) * (t(rep(1, h)) %*% Vhat)
  LB <- t(dhat) + (1 / (2 * lambda2)) * (t(rep(1, h)) %*% Vhat)
  
  return(list(LB = t(LB), UB = t(UB)))
}


MDprojl2tf <- function(delta, V, lambda1, lambda2, K) {
  # Load required libraries
  library(Matrix)
  
  p <- length(delta)
  
  
  # First difference matrix
  D1 <- matrix(0, nrow = p, ncol = p)
  D1[cbind(2:p, 1:(p-1))] <- 1  # Put 1's on the subdiagonal
  D1 <- D1 - diag(p)  
  D1 <- D1[-1, ]  # Remove first row (equivalent to D1(2:end,:))
  
  # Second and third difference matrices
  D2 <- D1[-1, -1]  # Remove first row and first column
  D3 <- D2[-1, -1]  # Remove first row and first column again
  D3 <- D3 %*% D2 %*% D1  # Third difference
  
  # l2 trend filtering problem formulation
  # min (delta - b)'*inv(V)*(delta - b) + lambda_1*b'*D1'*W1*D1*b +
  # lambda_2*b'*D3'*W2*D3*b
  
  iV <- solve(V)
  scaledV <- V / mean(diag(V))
  
  # Weight matrices
  vD1 <- D1 %*% scaledV %*% t(D1)
  zeros_block <- matrix(0, nrow = K-1, ncol = K-1)
  partial_vD1 <- vD1[K:nrow(vD1), K:ncol(vD1)] %>% as.matrix()
  vD1_block <- partial_vD1 / mean(diag(partial_vD1))
  W1 <- as.matrix(bdiag(zeros_block, vD1_block))
  
  vD3 <- D3 %*% scaledV %*% t(D3)
  W3 <- vD3 / mean(diag(vD3))
  
  scalediV <- solve(scaledV)
  
  # Denominator calculation
  den <- solve(scalediV + lambda1 * t(D1) %*% diag(diag(W1)) %*% D1 + 
                 lambda2 * t(D3) %*% diag(diag(W3)) %*% D3)
  
  # Degrees of freedom
  df <- sum(diag(den %*% scalediV))
  
  # Fitted values
  d <- den %*% (scalediV %*% delta)
  
  # Covariance matrix
  Vd <- den %*% scalediV %*% t(den)
  Vd <- Vd * mean(diag(V))
  
  # Standardized covariance and eigendecomposition
  stdVd <- sqrt(diag(Vd))
  cVdb <- diag(1/stdVd) %*% Vd %*% diag(1/stdVd)
  
  eigen_result <- eigen(cVdb)
  eigVec <- eigen_result$vectors
  D_eigen <- diag(eigen_result$values)
  
  # Check if eigenvalues are real
  if (any(Im(eigen_result$values) > 0)) {
    cat("hi - perhaps issue with covariance?\n")
  }
  
  # J matrix calculation
  sqrt_D <- sqrt(pmax(D_eigen, 1e-12))
  sqrt_D[D_eigen <= 1e-12] <- 0
  J <- eigVec %*% sqrt_D %*% t(eigVec)
  
  # Mahalanobis distance
  residual <- delta - d
  MD <- as.numeric(t(residual) %*% solve(V) %*% residual)
  
  # P-value using chi-squared distribution
  pv <- 1 - pchisq(MD, df = p - df)
  
  # Return results as a list
  return(list(
    d = as.vector(d),
    Vd = Vd,
    MD = MD,
    J = J,
    pv = pv,
    df = df
  ))
}


diff_df <- function(loglam1, loglam2, K, V, targetdf) {
  p <- nrow(V)
  
  # First difference matrix
  D1 <- matrix(0, nrow = p, ncol = p)
  D1[cbind(2:p, 1:(p-1))] <- 1  # Put 1's on the subdiagonal
  D1 <- D1 - diag(p)  
  D1 <- D1[-1, ]  # Remove first row (equivalent to D1(2:end,:))
  
  # Second and third difference matrices
  D2 <- D1[-1, -1]  # Remove first row and first column
  D3 <- D2[-1, -1]  # Remove first row and first column again
  D3 <- D3 %*% D2 %*% D1  # Third difference
  
  # Compute vD1 and W1
  vD1 <- D1 %*% V %*% t(D1)
  
  # Create block diagonal matrix W1
  zeros_block <- matrix(0, nrow = K-1, ncol = K-1)
  vD1_subset <- vD1[K:nrow(vD1), K:ncol(vD1)]
  vD1_normalized <- vD1_subset / mean(diag(vD1_subset))
  
  # Create W1 using block diagonal structure
  W1 <- matrix(0, nrow = nrow(vD1), ncol = ncol(vD1))
  W1[K:nrow(W1), K:ncol(W1)] <- vD1_normalized
  
  # Compute vD3 and W2
  vD3 <- D3 %*% V %*% t(D3)
  W2 <- vD3 / mean(diag(vD3))
  
  # Inverse of V
  iV <- solve(V)
  
  # Compute denominator
  den <- solve(iV + exp(loglam1) * t(D1) %*% diag(diag(W1)) %*% D1 + 
                 exp(loglam2) * t(D3) %*% diag(diag(W2)) %*% D3)
  
  # Compute degrees of freedom
  df <- sum(diag(den %*% iV))  # trace is sum of diagonal elements
  
  # Return squared difference
  f <- (df - targetdf)^2
  
  return(f)
}


my_df <- function(loglam1, loglam2, K, V) {
    library(Matrix)  # For block diagonal operations
    
    p <- nrow(V)
    
    # First difference matrix
    # Create D1: diag(ones(p-1,1),-1) - eye(p) in MATLAB
    D1 <- matrix(0, p, p)
    diag(D1[-1, ]) <- 1  # Lower diagonal
    D1 <- D1 - diag(1, p)  # Subtract identity
    D1 <- D1[-1, ]  # Remove first row: (p-1) x p
    
    # D2 is D1 with first row and first column removed
    D2 <- D1[-1, -1]  # (p-2) x (p-1)
    
    # D3 is D2 with first row and first column removed  
    D3 <- D2[-1, -1]  # (p-3) x (p-2)
    
    # Third difference: D3 * D2 * D1
    D3 <- D3 %*% D2 %*% D1  # (p-3) x p
    
    # Calculate vD1 and W1
    vD1 <- D1 %*% V %*% t(D1)
    
    # Create block diagonal matrix W1
    zeros_block <- matrix(0, K-1, K-1)
    vD1_subset <- vD1[K:nrow(vD1), K:ncol(vD1)] %>% as.matrix()
    normalized_vD1 <- vD1_subset / mean(diag(vD1_subset))
    W1 <- bdiag(zeros_block, normalized_vD1)
    W1 <- as.matrix(W1)  # Convert to regular matrix
    
    # Calculate vD3 and W2
    vD3 <- D3 %*% V %*% t(D3)
    W2 <- vD3 / mean(diag(vD3))
    
    # Calculate inverse of V
    iV <- solve(V)
    
    # Calculate denominator (using diagonal of W1 and W2)
    den <- solve(iV + exp(loglam1) * t(D1) %*% diag(diag(W1)) %*% D1 + 
                   exp(loglam2) * t(D3) %*% diag(diag(W2)) %*% D3)
    
    # Calculate degrees of freedom
    df_result <- sum(diag(den %*% iV))
    
    return(df_result)
  }

