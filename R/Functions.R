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
  
  if (anyNA(X)) {
    stop("X contains NA values. Please handle missing values first.", call. = FALSE)
  }
  
  if (any(X < 0)) {
    stop("X must contain non-negative values only.", call. = FALSE)
  }
  
  p <- ncol(X)
  
  if (p < 2) {
    stop("X must contain at least two components.", call. = FALSE)
  }
  
  row_totals <- rowSums(X)
  
  if (any(row_totals <= 0)) {
    stop("Each row must have positive total abundance.", call. = FALSE)
  }
  
  Pi <- sweep(X, 1, row_totals, FUN = "/")
  sqrt_Pi <- sqrt(Pi)
  
  sqrt_pi0 <- rep(1 / sqrt(p), p)
  
  c_pi <- as.vector(sqrt_Pi %*% sqrt_pi0)
  c_pi <- pmin(pmax(c_pi, 0), 1)
  
  s_pi <- sqrt(pmax(0, 1 - c_pi^2))
  
  scale_factor <- rep(1, length(s_pi))
  non_uniform <- s_pi > 0
  scale_factor[non_uniform] <- asin(s_pi[non_uniform]) / s_pi[non_uniform]
  
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


#' Simplex Hellinger Alpha Diversity
#'
#' Computes Simplex Hellinger alpha diversity from Hellinger-Riemann intrinsic
#' coordinates.
#'
#' The dominance score is the geodesic distance from the uniform composition:
#' `D(pi) = ||HRIC(pi)|| = asin(s(pi))`.
#'
#' The evenness score is:
#' `1 - D(pi) / acos(1 / sqrt(p))`.
#'
#' @param X Numeric matrix or data frame. Rows are samples and columns are
#'   components, taxa, or features.
#' @param measure Alpha diversity measure to return. Options are
#'   `"dominance"`, `"evenness"`, or `"both"`.
#'
#' @return A numeric vector if `measure` is `"dominance"` or `"evenness"`;
#'   otherwise a data frame with both measures.
#'
#' @export
Simplex_Hellinger_alpha <- function(
    X,
    measure = c("dominance", "evenness", "both")
) {
  
  measure <- match.arg(measure)
  
  prep <- .prepare_simplex_hellinger(X)
  
  dominance <- asin(prep$s_pi)
  
  D_max <- acos(1 / sqrt(prep$p))
  
  evenness <- 1 - dominance / D_max
  evenness <- pmin(pmax(evenness, 0), 1)
  
  if (measure == "both") {
    out <- data.frame(
      dominance = dominance,
      evenness = evenness
    )
    
    rownames(out) <- rownames(prep$X)
    
    return(out)
  }
  
  out <- switch(
    measure,
    dominance = dominance,
    evenness = evenness
  )
  
  names(out) <- rownames(prep$X)
  
  out
}


#' Simplex Hellinger Beta Diversity
#'
#' Computes Simplex Hellinger beta diversity as the Euclidean distance between
#' Hellinger-Riemann intrinsic coordinate vectors.
#'
#' @param X Numeric matrix or data frame. Rows are samples and columns are
#'   components, taxa, or features.
#' @param output Output type. Either `"dist"` or `"matrix"`.
#' @param diag Logical. Include the diagonal if returning a `dist` object.
#' @param upper Logical. Include the upper triangle if returning a `dist` object.
#'
#' @return A `dist` object or a square distance matrix.
#'
#' @export
Simplex_Hellinger_beta <- function(
    X,
    output = c("dist", "matrix"),
    diag = FALSE,
    upper = FALSE
) {
  
  output <- match.arg(output)
  
  coordinates <- HRIC(X)
  
  out <- stats::dist(
    coordinates,
    method = "euclidean",
    diag = diag,
    upper = upper
  )
  
  attr(out, "method") <- "Simplex Hellinger"
  
  if (output == "matrix") {
    return(as.matrix(out))
  }
  
  out
}