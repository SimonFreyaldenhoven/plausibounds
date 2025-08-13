# Internal utility functions

# Calculate sup-t bands
bands_plugin <- function(delta, var, p, nsim = 2000, level = 0.95) {
  # Simulate from N(0, var), calculate sup-t bands, return LB and UB
  rd <- MASS::mvrnorm(nsim, mu = rep(0, p), Sigma = var)
  std_devs <- sqrt(diag(var))
  rd_standardized <- rd / matrix(std_devs, nrow = nsim, ncol = length(std_devs), byrow = TRUE)
  sup_t <- as.numeric(quantile(apply(abs(rd_standardized), 1, max), level))
  
  list(
    LB = delta - sup_t * sqrt(diag(var)),
    UB = delta + sup_t * sqrt(diag(var)),
    sup_t = sup_t
  )
}

