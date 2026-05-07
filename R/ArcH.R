# ------------------------------------------------------------
# Angular Hellinger Compositional Normalization
# ------------------------------------------------------------
# Input:
#   X: matrix or data.frame
#      rows = samples, columns = taxa/features/components
#
# Output:
#   transformed matrix with same dimension as X
#
# Formula:
#   pi_i = X_i / sum_j X_ij
#   c(pi_i) = (1 / sqrt(k)) * sum_j sqrt(pi_ij)
#   s(pi_i) = sqrt(1 - c(pi_i)^2)
#   T(pi_i) = asin(s(pi_i)) / s(pi_i) *
#             (sqrt(pi_i) - c(pi_i) * sqrt(pi_0))
#
# where sqrt(pi_0) = (1/sqrt(k), ..., 1/sqrt(k)).
# ------------------------------------------------------------

angular_hellinger_normalization <- function(X) {
  
  # Convert to matrix
  X <- as.matrix(X)
  
  # Check input
  if (!is.numeric(X)) {
    stop("X must be numeric.")
  }
  
  if (any(X < 0, na.rm = TRUE)) {
    stop("X must contain non-negative values only.")
  }
  
  if (anyNA(X)) {
    stop("X contains NA values. Please handle missing values first.")
  }
  
  # Number of components
  k <- ncol(X)
  
  if (k < 2) {
    stop("X must contain at least two components.")
  }
  
  # Row sums
  rs <- rowSums(X)
  
  if (any(rs <= 0)) {
    stop("Each row must have positive total abundance.")
  }
  
  # Convert count/abundance data to compositional data
  Pi <- sweep(X, 1, rs, FUN = "/")
  
  # Square-root composition
  sqrt_Pi <- sqrt(Pi)
  
  # Uniform direction in square-root space
  sqrt_pi0 <- rep(1 / sqrt(k), k)
  
  # c(pi): cosine-type Hellinger similarity to the uniform composition
  c_pi <- as.vector(sqrt_Pi %*% sqrt_pi0)
  
  # Numerical protection: theoretically c(pi) is in [0, 1]
  c_pi <- pmin(pmax(c_pi, 0), 1)
  
  # s(pi) = sqrt(1 - c(pi)^2)
  s_pi <- sqrt(1 - c_pi^2)
  
  # Scaling factor
  scale_factor <- asin(s_pi) / s_pi
  
  # Handle the exact uniform-composition case where s(pi) = 0
  scale_factor[s_pi == 0] <- 1
  
  # Centered square-root composition
  centered <- sqrt_Pi - outer(c_pi, sqrt_pi0)
  
  # Final transformation
  T <- centered * scale_factor
  
  # Keep row and column names
  rownames(T) <- rownames(X)
  colnames(T) <- colnames(X)
  
  return(T)
}