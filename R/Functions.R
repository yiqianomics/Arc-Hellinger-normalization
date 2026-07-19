# =====================================================================
# Simplex Hellinger tools: coordinates and diversity
# =====================================================================
# Coordinates on the closed simplex:
#   HRIC()    centered Hellinger-Riemann intrinsic coordinates
#   RHRIC()   reference-based coordinates
#
# Diversity built from HRIC coordinates:
#   SHalpha() per-sample evenness (alpha)
#   SHgamma() evenness of the normalized mean Hellinger direction (gamma)
#   SHbeta()  non-negative spherical partition, gamma - mean(alpha) (beta)
#   SHdelta() between-group dispersion after removing the group-size
#             weighted within-group dispersion (delta)
#
# The file uses base R only (no package dependencies).
# =====================================================================


# ---------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------

#' Prepare Simplex Hellinger Quantities
#'
#' Internal helper for Simplex Hellinger transformations.
#'
#' @param X Numeric matrix or data frame.
#'
#' @return A list of intermediate quantities.
#' @keywords internal
.prepare_simplex_hellinger <- function(X) {
  
  X <- as.matrix(X)
  
  if (!is.numeric(X)) {
    stop("X must be numeric.", call. = FALSE)
  }
  
  if (nrow(X) < 1L) {
    stop("X must contain at least one sample.", call. = FALSE)
  }
  
  if (any(!is.finite(X))) {
    stop(
      "X must contain finite, non-missing values only.",
      call. = FALSE
    )
  }
  
  if (any(X < 0)) {
    stop(
      "X must contain non-negative values only.",
      call. = FALSE
    )
  }
  
  p <- ncol(X)
  
  if (p < 2L) {
    stop(
      "X must contain at least two components.",
      call. = FALSE
    )
  }
  
  row_totals <- rowSums(X)
  
  if (any(row_totals <= 0)) {
    stop(
      "Each row must have positive total abundance.",
      call. = FALSE
    )
  }
  
  Pi <- sweep(X, 1, row_totals, FUN = "/")
  sqrt_Pi <- sqrt(Pi)
  
  sqrt_pi0 <- rep(1 / sqrt(p), p)
  
  c_pi <- as.vector(sqrt_Pi %*% sqrt_pi0)
  c_pi <- pmin(pmax(c_pi, 0), 1)
  
  s_pi <- sqrt(pmax(0, 1 - c_pi^2))
  
  scale_factor <- rep(1, length(s_pi))
  non_uniform <- s_pi > 0
  scale_factor[non_uniform] <-
    asin(s_pi[non_uniform]) / s_pi[non_uniform]
  
  list(
    X = X,
    p = p,
    Pi = Pi,
    sqrt_Pi = sqrt_Pi,
    sqrt_pi0 = sqrt_pi0,
    c_pi = c_pi,
    s_pi = s_pi,
    scale_factor = scale_factor
  )
}

#' Match Reference Component
#'
#' Internal helper for reference-based coordinates.
#'
#' @param reference Column index or column name.
#' @param X Numeric matrix.
#'
#' @return Integer column index.
#' @keywords internal
.match_reference_component <- function(reference, X) {
  
  p <- ncol(X)
  
  if (length(reference) != 1) {
    stop("reference must be a single column index or column name.", call. = FALSE)
  }
  
  if (is.character(reference)) {
    if (is.null(colnames(X))) {
      stop("reference cannot be a column name when X has no column names.", call. = FALSE)
    }
    
    reference <- match(reference, colnames(X))
  }
  
  if (!is.numeric(reference) || is.na(reference) || reference != floor(reference)) {
    stop("reference must be a valid column index or column name.", call. = FALSE)
  }
  
  reference <- as.integer(reference)
  
  if (reference < 1 || reference > p) {
    stop("reference must be a valid column index or column name.", call. = FALSE)
  }
  
  reference
}


# ---------------------------------------------------------------------
# Coordinates
# ---------------------------------------------------------------------

#' Hellinger-Riemann Intrinsic Coordinates
#'
#' Computes centered Hellinger-Riemann intrinsic coordinates for compositional
#' data on the closed simplex.
#'
#' @param X Numeric matrix or data frame. Rows are samples and columns are
#'   components, taxa, or features.
#'
#' @return A numeric matrix with the same dimensions as `X`.
#'
#' @export
HRIC <- function(X) {
  
  prep <- .prepare_simplex_hellinger(X)
  
  centered <- prep$sqrt_Pi - outer(prep$c_pi, prep$sqrt_pi0)
  out <- sweep(centered, 1, prep$scale_factor, FUN = "*")
  
  rownames(out) <- rownames(prep$X)
  colnames(out) <- colnames(prep$X)
  
  out
}


#' Reference-Based Hellinger-Riemann Intrinsic Coordinates
#'
#' Computes reference-based Hellinger-Riemann intrinsic coordinates.
#'
#' For reference component `r`, the coordinates are
#' `a(pi) * (sqrt(pi_j) - sqrt(pi_r))` for all `j != r`.
#'
#' @param X Numeric matrix or data frame. Rows are samples and columns are
#'   components, taxa, or features.
#' @param reference Reference component, given as a column index or column name.
#'   Defaults to the last column.
#'
#' @return A numeric matrix with `ncol(X) - 1` columns.
#'
#' @export
RHRIC <- function(X, reference = NULL) {
  
  prep <- .prepare_simplex_hellinger(X)
  
  if (is.null(reference)) {
    reference <- prep$p
  }
  
  reference <- .match_reference_component(reference, prep$X)
  keep <- setdiff(seq_len(prep$p), reference)
  
  out <- sweep(
    prep$sqrt_Pi[, keep, drop = FALSE],
    1,
    prep$sqrt_Pi[, reference],
    FUN = "-"
  )
  
  out <- sweep(out, 1, prep$scale_factor, FUN = "*")
  
  component_names <- colnames(prep$X)
  if (is.null(component_names)) {
    component_names <- paste0("V", seq_len(prep$p))
  }
  
  rownames(out) <- rownames(prep$X)
  colnames(out) <- component_names[keep]
  
  attr(out, "reference") <- component_names[reference]
  attr(out, "reference_index") <- reference
  
  out
}


# ---------------------------------------------------------------------
# Diversity
# ---------------------------------------------------------------------

#' Simplex Hellinger Alpha Diversity
#'
#' Computes the per-sample Simplex Hellinger alpha diversity from the
#' Hellinger-Riemann intrinsic coordinates. The measure is an evenness score.
#'
#' For sample \eqn{i} with coordinate vector \eqn{Z_i = \mathrm{HRIC}(\pi_i)},
#' \deqn{\alpha_{\mathrm{SH}}(\pi_i) = 1 -
#'       \frac{\lVert Z_i \rVert_2^2}{A_p^2}, \qquad
#'       A_p = \arcsin\!\left(\sqrt{1 - 1/p}\right),}
#' where \eqn{A_p} is the maximum Simplex Hellinger distance from the uniform
#' composition and \eqn{p} is the number of components.
#'
#' Values near 1 indicate a composition close to uniform abundance across
#' components. Lower values indicate concentration toward fewer components.
#' Exact zeros are allowed because HRIC is defined on the closed simplex.
#'
#' @param X Numeric matrix or data frame. Rows are samples and columns are
#'   components, taxa, or features.
#'
#' @return A named numeric vector of length `nrow(X)`, one alpha value per
#'   sample, each in \eqn{[0, 1]}.
#'
#' @examples
#' X <- matrix(c(1, 1, 1, 1,
#'               10, 1, 1, 1,
#'               1, 0, 0, 0),
#'             nrow = 3, byrow = TRUE)
#' SHalpha(X)
#'
#' @export
SHalpha <- function(X) {
  
  Z <- HRIC(X)
  p <- ncol(Z)
  A_p <- asin(sqrt(1 - 1 / p))
  
  norm_sq <- rowSums(Z^2)
  alpha <- 1 - norm_sq / A_p^2
  alpha <- pmin(pmax(alpha, 0), 1)
  
  names(alpha) <- rownames(Z)
  alpha
}


#' Simplex Hellinger Gamma Diversity
#'
#' Computes the Simplex Hellinger gamma diversity from the normalized mean
#' Hellinger direction of the sample compositions.
#'
#' For sample \eqn{i}, let
#' \deqn{\pi_i =
#'       \frac{X_i}{\sum_{j=1}^p X_{ij}}}
#' denote its total-sum-scaled composition, and let
#' \deqn{q_i = \sqrt{\pi_i}}
#' denote its element-wise square-root representation. Each \eqn{q_i} lies on
#' the positive orthant of the unit sphere because
#' \eqn{\lVert q_i\rVert_2 = 1}.
#'
#' The cohort Hellinger direction is obtained by averaging the square-root
#' compositions and projecting the result back onto the unit sphere:
#' \deqn{\bar{q} = \frac{1}{n}\sum_{i=1}^n q_i, \qquad
#'       q_\gamma =
#'       \frac{\bar{q}}{\lVert\bar{q}\rVert_2}.}
#'
#' The corresponding cohort composition is
#' \deqn{\pi_\gamma = q_\gamma^{\circ 2},}
#' where the square is applied element-wise. The gamma diversity is then
#' \deqn{\gamma_{\mathrm{SH}} =
#'       1 -
#'       \frac{
#'       \lVert\mathrm{HRIC}(\pi_\gamma)\rVert_2^2
#'       }{A_p^2},
#'       \qquad
#'       A_p = \arcsin\!\left(\sqrt{1 - 1/p}\right).}
#'
#' Each sample contributes equally, regardless of its total abundance or
#' sequencing depth.
#'
#' @param X Numeric matrix or data frame. Rows are samples and columns are
#'   components, taxa, or features.
#'
#' @return A single numeric value in \eqn{[0, 1]}.
#'
#' @examples
#' X <- matrix(c(1, 1, 1, 1,
#'               10, 1, 1, 1,
#'               1, 0, 0, 0),
#'             nrow = 3, byrow = TRUE)
#' SHgamma(X)
#'
#' @export
SHgamma <- function(X) {
  
  prep <- .prepare_simplex_hellinger(X)
  
  ## Each row of prep$sqrt_Pi is sqrt(pi_i) and has Euclidean norm one.
  mean_sqrt_composition <- colMeans(prep$sqrt_Pi)
  
  ## Project the mean Hellinger direction back onto the unit sphere.
  mean_sqrt_norm <- sqrt(sum(mean_sqrt_composition^2))
  
  if (!is.finite(mean_sqrt_norm) || mean_sqrt_norm <= 0) {
    stop(
      "The mean square-root composition has invalid Euclidean norm.",
      call. = FALSE
    )
  }
  
  gamma_sqrt_composition <-
    mean_sqrt_composition / mean_sqrt_norm
  
  ## Convert the normalized Hellinger direction back to a composition.
  gamma_composition <- gamma_sqrt_composition^2
  
  ## Protect against floating-point deviation from a unit sum.
  gamma_composition <-
    gamma_composition / sum(gamma_composition)
  
  gamma_composition <- matrix(
    gamma_composition,
    nrow = 1,
    dimnames = list(NULL, colnames(prep$X))
  )
  
  Z_gamma <- HRIC(gamma_composition)
  
  p <- ncol(Z_gamma)
  A_p <- asin(sqrt(1 - 1 / p))
  
  gamma <- 1 - sum(Z_gamma^2) / A_p^2
  
  min(max(gamma, 0), 1)
}


#' Simplex Hellinger Beta Diversity
#'
#' Computes the non-negative additive beta component as the difference between
#' the gamma diversity of the normalized mean Hellinger direction and the mean
#' sample-level alpha diversity:
#' \deqn{\beta_{\mathrm{SH}} =
#'       \gamma_{\mathrm{SH}} -
#'       \frac{1}{n}\sum_{i=1}^n
#'       \alpha_{\mathrm{SH}}(\pi_i).}
#'
#' Let
#' \deqn{q_i = \sqrt{\pi_i}, \qquad
#'       q_\gamma =
#'       \frac{\sum_{i=1}^n q_i}
#'       {\left\lVert\sum_{i=1}^n q_i\right\rVert_2},}
#' and let
#' \deqn{\pi_\gamma = q_\gamma^{\circ 2}.}
#' Then
#' \deqn{\beta_{\mathrm{SH}} =
#'       \frac{
#'       n^{-1}\sum_{i=1}^n
#'       \lVert\mathrm{HRIC}(\pi_i)\rVert_2^2
#'       -
#'       \lVert\mathrm{HRIC}(\pi_\gamma)\rVert_2^2
#'       }{A_p^2}.}
#'
#' This quantity is non-negative because the normalized mean Hellinger
#' direction cannot be farther from the uniform direction, in squared
#' spherical distance, than the average squared distance of the individual
#' sample directions.
#'
#' Consequently,
#' \deqn{\gamma_{\mathrm{SH}} =
#'       \frac{1}{n}\sum_{i=1}^n
#'       \alpha_{\mathrm{SH}}(\pi_i)
#'       + \beta_{\mathrm{SH}}.}
#'
#' This beta component is based on the normalized Hellinger-spherical cohort
#' center. It is not, in general, identical to the Euclidean variance of the
#' HRIC coordinate vectors around their coordinate-wise mean.
#'
#' @param X Numeric matrix or data frame. Rows are samples and columns are
#'   components, taxa, or features.
#'
#' @return A single numeric value in \eqn{[0, 1]}.
#'
#' @seealso [SHalpha()], [SHgamma()]
#'
#' @examples
#' X <- matrix(c(1, 1, 1, 1,
#'               10, 1, 1, 1,
#'               1, 0, 0, 0),
#'             nrow = 3, byrow = TRUE)
#' SHbeta(X)
#'
#' @export
SHbeta <- function(X) {
  
  gamma <- SHgamma(X)
  mean_alpha <- mean(SHalpha(X))
  
  beta <- gamma - mean_alpha
  
  ## Theoretical beta is in [0, 1]; clipping removes only numerical error.
  min(max(beta, 0), 1)
}

#' Simplex Hellinger Between-Group Turnover
#'
#' Partitions the total dispersion of the Hellinger-Riemann intrinsic
#' coordinates into within-group and between-group components and returns the
#' between-group component \eqn{\delta}. This is the residual turnover that
#' remains after subtracting the group-size weighted within-group dispersion.
#'
#' Let \eqn{\bar{Z}} be the grand center and \eqn{\bar{Z}_k} the center of
#' group \eqn{k} with \eqn{n_k} samples. Define
#' \deqn{D_{\mathrm{cohort}} = \frac{1}{n} \sum_{i=1}^n
#'       \lVert Z_i - \bar{Z} \rVert_2^2, \qquad
#'       D_{g_k} = \frac{1}{n_k} \sum_{g_i = k}
#'       \lVert Z_i - \bar{Z}_k \rVert_2^2 .}
#' Then
#' \deqn{\delta = D_{\mathrm{cohort}} -
#'       \sum_{k} \frac{n_k}{n} \, D_{g_k} .}
#' A positive \eqn{\delta} means some turnover remains after removing the
#' within-group dispersion. The value is non-negative and generalizes to any
#' number of groups (two or more).
#'
#' @param X Numeric matrix or data frame. Rows are samples and columns are
#'   components, taxa, or features.
#' @param group A vector of group labels with one entry per row of `X`, in the
#'   same order as the rows. Coerced to a factor; unused levels are dropped.
#'   At least two distinct groups are required.
#'
#' @return A single non-negative numeric value.
#'
#' @examples
#' X <- matrix(c(1, 1, 1, 1,
#'               8, 1, 1, 1,
#'               1, 1, 1, 8,
#'               1, 1, 8, 1,
#'               2, 2, 1, 1,
#'               1, 1, 2, 2),
#'             nrow = 6, byrow = TRUE)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' SHdelta(X, grp)
#'
#' @export
SHdelta <- function(X, group) {
  
  if (missing(group)) {
    stop("group is required for SHdelta.", call. = FALSE)
  }
  
  Z <- HRIC(X)
  n <- nrow(Z)
  
  if (length(group) != n) {
    stop(
      "group must have one entry per row of X (length nrow(X)).",
      call. = FALSE
    )
  }
  
  if (anyNA(group)) {
    stop(
      "group contains NA values. Please handle missing values first.",
      call. = FALSE
    )
  }
  
  g <- droplevels(as.factor(group))
  
  if (nlevels(g) < 2L) {
    stop(
      "group must contain at least two distinct groups.",
      call. = FALSE
    )
  }
  
  Z_bar <- colMeans(Z)
  
  between <- 0
  
  for (k in levels(g)) {
    idx <- which(g == k)
    n_k <- length(idx)
    
    Z_k <- Z[idx, , drop = FALSE]
    Z_bar_k <- colMeans(Z_k)
    
    center_difference <- Z_bar_k - Z_bar
    
    between <- between +
      (n_k / n) * sum(center_difference^2)
  }
  
  max(between, 0)
}