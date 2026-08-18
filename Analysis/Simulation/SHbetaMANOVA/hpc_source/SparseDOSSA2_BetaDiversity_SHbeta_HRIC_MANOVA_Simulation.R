.libPaths(c("/home/zhang.16383/RPackage_r451", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
rep_id <- as.integer(args[1])
if (is.na(rep_id) || rep_id < 1L) stop("rep_id must be a positive integer.")

source("polar_cluster.R")

suppressPackageStartupMessages({
  library(ape)
  library(cluster)
  library(SparseDOSSA2)
  library(phyloseq)
  library(GUniFrac)
  library(vegan)
  library(MiSPU)   # for throat.otu.tab / throat.tree data in your original workflow
  library(HRIC)    # for SHbeta(), HRIC(), and HRIC MANOVA decomposition
})

## -------------------- user / HPC settings --------------------
FIT_RDS      <- Sys.getenv("FIT_RDS", "SparseDOSSA2_fit_URT_lambda0.1.rds")
OUTDIR       <- Sys.getenv("OUTDIR", "beta_div_sparse12_out")
BASE_SEED    <- as.integer(Sys.getenv("BASE_SEED", "12345"))
N_SAMPLES    <- as.integer(Sys.getenv("N_SAMPLES", "100"))
ADONIS_PERM  <- as.integer(Sys.getenv("ADONIS_PERM", "999"))

## If K_COMM is supplied (the revised Slurm file supplies K_COMM=12), run only
## those comma-separated community counts. If it is absent, retain the original
## grid seq(2, 20, 2).
K_COMM_ENV <- trimws(Sys.getenv("K_COMM", ""))
if (nzchar(K_COMM_ENV)) {
  K_COMM_GRID <- suppressWarnings(as.integer(strsplit(K_COMM_ENV, ",", fixed = TRUE)[[1]]))
  if (!length(K_COMM_GRID) || anyNA(K_COMM_GRID) || any(K_COMM_GRID < 2L)) {
    stop("K_COMM must contain integer values >= 2, optionally comma-separated.")
  }
  K_COMM_GRID <- sort(unique(K_COMM_GRID))
} else {
  K_COMM_GRID <- seq(2, 20, 2)
}

if (!file.exists(FIT_RDS)) stop("FIT_RDS not found: ", FIT_RDS)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "raw"), recursive = TRUE, showWarnings = FALSE)

## beta settings: beta = 0 is null / type I setting
beta_grid <- seq(0, 1, 0.05)
signal_mode_grid <- c("abundance", "prevalence", "both")

## ============================================================
## Pseudocount settings for Aitchison and ILR + Euclidean
##
## Original method:
##   aitchison = vegan::vegdist(..., method = "aitchison", pseudocount = 0.5)
##
## Added Aitchison pseudocounts:
##   1, 0.1, 0.01
##
## Added ILR + Euclidean distances using the same pseudocounts:
##   0.5, 1, 0.1, 0.01
##
## Notes:
##   - ILR + Euclidean is theoretically equivalent to Aitchison distance
##     when the same positive data matrix and same pseudocount are used.
##   - Here ILR is implemented manually using an orthonormal Helmert basis,
##     so no additional package is required.
## ============================================================

AITCHISON_METHODS <- c(
  "aitchison",          # original: pseudocount = 0.5
  "aitchison_pc1",
  "aitchison_pc0.1",
  "aitchison_pc0.01"
)

ILR_EUCLIDEAN_METHODS <- c(
  "ILR_euclidean_pc0.5",
  "ILR_euclidean_pc1",
  "ILR_euclidean_pc0.1",
  "ILR_euclidean_pc0.01"
)

DISTANCE_PSEUDOCOUNT <- setNames(
  c(
    0.5,
    1,
    0.1,
    0.01,
    0.5,
    1,
    0.1,
    0.01
  ),
  c(
    AITCHISON_METHODS,
    ILR_EUCLIDEAN_METHODS
  )
)

## -------------------- helpers --------------------
row_tss <- function(x) {
  x <- as.matrix(x)
  rs <- rowSums(x)
  out <- matrix(0, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  keep <- is.finite(rs) & rs > 0
  if (any(keep)) {
    out[keep, ] <- x[keep, , drop = FALSE] / rs[keep]
  }
  out
}

## ============================================================
## AHC transformation:
##
## Given a composition pi = (pi_1, ..., pi_k),
##
##   c(pi) = (1 / sqrt(k)) * sum_i sqrt(pi_i)
##   s(pi) = sqrt(1 - c(pi)^2)
##   T(pi) = asin(s(pi)) / s(pi) * (sqrt(pi) - c(pi) * u0)
##
## where u0 = (1/sqrt(k), ..., 1/sqrt(k)) is the square-root
## uniform direction.
##
## This is the centered AHC version. Euclidean distance is then
## computed on transformed rows.
## ============================================================
AHC_transform <- function(counts_mat, zero_row_to_uniform = TRUE, eps_clip = 1e-12) {
  X <- as.matrix(counts_mat)
  n <- nrow(X)
  k <- ncol(X)
  
  if (k < 2L) {
    stop("AHC_transform requires at least two features.")
  }
  
  rs <- rowSums(X)
  P <- matrix(0, nrow = n, ncol = k, dimnames = dimnames(X))
  
  keep <- is.finite(rs) & rs > 0
  
  if (any(keep)) {
    P[keep, ] <- X[keep, , drop = FALSE] / rs[keep]
  }
  
  if (any(!keep)) {
    if (zero_row_to_uniform) {
      P[!keep, ] <- 1 / k
    } else {
      P[!keep, ] <- NA_real_
    }
  }
  
  ## Numerical guard: keep compositions inside [0, 1]
  P[P < 0] <- 0
  P[P > 1] <- 1
  
  Z <- sqrt(P)
  
  ## Uniform direction in square-root space
  u0 <- rep(1 / sqrt(k), k)
  
  ## c(pi): cosine similarity between sqrt(pi) and uniform direction
  c_val <- as.numeric(Z %*% u0)
  c_val <- pmin(pmax(c_val, -1), 1)
  
  ## s(pi)
  s_val <- sqrt(pmax(0, 1 - c_val^2))
  
  ## scale = asin(s) / s, with limit 1 when s -> 0
  scale <- rep(1, length(s_val))
  nonzero <- is.finite(s_val) & s_val > eps_clip
  scale[nonzero] <- asin(s_val[nonzero]) / s_val[nonzero]
  
  centered <- Z - tcrossprod(c_val, u0)
  
  Tmat <- centered * scale
  
  rownames(Tmat) <- rownames(X)
  colnames(Tmat) <- colnames(X)
  
  Tmat
}

safe_AHC_euclidean_distance <- function(counts_mat) {
  out <- tryCatch({
    Tmat <- AHC_transform(counts_mat)
    if (any(!is.finite(Tmat))) return(NULL)
    stats::dist(Tmat, method = "euclidean")
  }, error = function(e) {
    NULL
  })
  
  out
}

## ============================================================
## Aitchison distance with user-specified pseudocount
##
## This wraps vegan::vegdist(..., method = "aitchison").
## ============================================================
safe_aitchison_distance <- function(counts_mat, pseudocount = 0.5) {
  out <- tryCatch({
    X <- as.matrix(counts_mat)
    
    if (nrow(X) < 2L) {
      stop("safe_aitchison_distance requires at least two samples.")
    }
    
    if (ncol(X) < 2L) {
      stop("safe_aitchison_distance requires at least two features.")
    }
    
    if (any(!is.finite(X))) {
      stop("counts_mat contains NA, NaN, or Inf values.")
    }
    
    if (any(X < 0)) {
      stop("Aitchison distance requires non-negative counts.")
    }
    
    if (!is.finite(pseudocount) || pseudocount <= 0) {
      stop("pseudocount must be a positive finite number.")
    }
    
    vegan::vegdist(X, method = "aitchison", pseudocount = pseudocount)
  }, error = function(e) {
    NULL
  })
  
  out
}

## ============================================================
## ILR transformation and ILR + Euclidean distance
##
## For a positive vector x, ILR coordinates are obtained by projecting
## log(x) onto an orthonormal basis of the zero-sum contrast space.
##
## Here, x = counts + pseudocount.
##
## Since the Helmert basis is orthonormal and has column sums equal to 0,
## Euclidean distance on these ILR coordinates is equivalent to Aitchison
## distance on the same pseudocount-adjusted data.
## ============================================================
ILR_transform <- function(counts_mat, pseudocount = 0.5) {
  X <- as.matrix(counts_mat)
  
  if (nrow(X) < 1L) {
    stop("ILR_transform requires at least one sample.")
  }
  
  if (ncol(X) < 2L) {
    stop("ILR_transform requires at least two features.")
  }
  
  if (any(!is.finite(X))) {
    stop("counts_mat contains NA, NaN, or Inf values.")
  }
  
  if (any(X < 0)) {
    stop("ILR transformation requires non-negative counts before pseudocount addition.")
  }
  
  if (!is.finite(pseudocount) || pseudocount <= 0) {
    stop("pseudocount must be a positive finite number.")
  }
  
  X_pc <- X + pseudocount
  
  if (any(!is.finite(X_pc)) || any(X_pc <= 0)) {
    stop("pseudocount-adjusted matrix must be positive and finite.")
  }
  
  k <- ncol(X_pc)
  
  ## Helmert contrast basis: k x (k - 1)
  ## stats::contr.helmert(k) gives orthogonal zero-sum columns.
  ## We normalize each column to make the basis orthonormal.
  H <- stats::contr.helmert(k)
  V <- sweep(H, 2L, sqrt(colSums(H^2)), FUN = "/")
  
  Z <- log(X_pc) %*% V
  
  rownames(Z) <- rownames(X)
  colnames(Z) <- paste0("ilr", seq_len(ncol(Z)))
  
  Z
}

safe_ILR_euclidean_distance <- function(counts_mat, pseudocount = 0.5) {
  out <- tryCatch({
    Z <- ILR_transform(counts_mat, pseudocount = pseudocount)
    if (any(!is.finite(Z))) return(NULL)
    stats::dist(Z, method = "euclidean")
  }, error = function(e) {
    NULL
  })
  
  out
}

## ============================================================
## User-defined Jaccard distance
##
## binary = TRUE:
##   Standard presence/absence Jaccard distance:
##
##     d(i, j) = 1 - |A_i ∩ A_j| / |A_i ∪ A_j|
##
##   where A_i is the set of OTUs observed in sample i.
##
## binary = FALSE:
##   Quantitative Jaccard / Ruzicka distance:
##
##     d(i, j) = 1 - sum_k min(x_ik, x_jk) / sum_k max(x_ik, x_jk)
##
## Recommendation:
##   For microbiome beta-diversity, if you mean ordinary Jaccard,
##   use binary = TRUE.
## ============================================================

my_jaccard_distance <- function(counts_mat, binary = TRUE) {
  X <- as.matrix(counts_mat)
  
  if (nrow(X) < 2L) {
    stop("my_jaccard_distance requires at least two samples.")
  }
  
  if (ncol(X) < 1L) {
    stop("my_jaccard_distance requires at least one feature.")
  }
  
  if (any(!is.finite(X))) {
    stop("counts_mat contains NA, NaN, or Inf values.")
  }
  
  if (any(X < 0)) {
    stop("Jaccard distance requires non-negative values.")
  }
  
  n <- nrow(X)
  D <- matrix(0, nrow = n, ncol = n)
  
  rownames(D) <- rownames(X)
  colnames(D) <- rownames(X)
  
  if (binary) {
    ## Presence / absence Jaccard
    B <- X > 0
    
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        intersection_ij <- sum(B[i, ] & B[j, ])
        union_ij        <- sum(B[i, ] | B[j, ])
        
        if (union_ij == 0L) {
          d_ij <- 0
        } else {
          d_ij <- 1 - intersection_ij / union_ij
        }
        
        D[i, j] <- d_ij
        D[j, i] <- d_ij
      }
    }
    
  } else {
    ## Quantitative Jaccard / Ruzicka distance
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        numerator_ij   <- sum(pmin(X[i, ], X[j, ]))
        denominator_ij <- sum(pmax(X[i, ], X[j, ]))
        
        if (denominator_ij == 0) {
          d_ij <- 0
        } else {
          d_ij <- 1 - numerator_ij / denominator_ij
        }
        
        D[i, j] <- d_ij
        D[j, i] <- d_ij
      }
    }
  }
  
  stats::as.dist(D)
}

safe_my_jaccard_distance <- function(counts_mat, binary = TRUE) {
  out <- tryCatch({
    D <- my_jaccard_distance(counts_mat, binary = binary)
    if (any(!is.finite(as.numeric(D)))) return(NULL)
    D
  }, error = function(e) {
    NULL
  })
  
  out
}

get_est_mean_abs <- function(fit_obj) {
  f <- fit_obj$EM_fit$fit
  if (is.null(f)) f <- fit_obj$fit
  if (is.null(f)) stop("Cannot find EM_fit$fit or fit in FIT_RDS object.")
  
  pi0 <- f$pi0
  mu  <- f$mu
  
  if (!is.null(f$sigma)) {
    sigma <- f$sigma
  } else if (!is.null(f$sigma2)) {
    sigma <- sqrt(f$sigma2)
  } else {
    stop("Cannot find sigma or sigma2 in fit object.")
  }
  
  if (is.null(names(mu)) || is.null(names(pi0)) || is.null(names(sigma))) {
    stop("mu/pi0/sigma must be named by OTU IDs.")
  }
  
  common_ids <- Reduce(intersect, list(names(mu), names(pi0), names(sigma)))
  if (!length(common_ids)) stop("No overlapping OTU IDs across mu/pi0/sigma.")
  
  mu    <- mu[common_ids]
  pi0   <- pi0[common_ids]
  sigma <- sigma[common_ids]
  
  (1 - pi0) * exp(mu + 0.5 * sigma^2)
}

safe_distance <- function(physeq, method) {
  otu_mat <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
  
  out <- tryCatch({
    if (method == "AHC_euclidean") {
      safe_AHC_euclidean_distance(otu_mat)
      
    } else if (method == "GLaD_0.5_weighted") {
      if (!requireNamespace("GLaD", quietly = TRUE)) return(NULL)
      GLaD::GLaD(physeq)
      
    } else if (method == "GLaD_0.5_unweighted") {
      if (!requireNamespace("GLaD", quietly = TRUE)) return(NULL)
      GLaD::GLaD(physeq, weighted = FALSE)
      
    } else if (method == "GLaD_1_weighted") {
      if (!requireNamespace("GLaD", quietly = TRUE)) return(NULL)
      GLaD::GLaD_eigen(physeq)
      
    } else if (method == "GLaD_1_unweighted") {
      if (!requireNamespace("GLaD", quietly = TRUE)) return(NULL)
      GLaD::GLaD_eigen(physeq, weighted = FALSE)
      
    } else if (method == "GUniFrac") {
      tree <- phy_tree(physeq)
      unifrac_arr <- GUniFrac::GUniFrac(otu_mat, tree, alpha = c(0, 0.5, 1))$unifracs
      as.dist(unifrac_arr[, , "d_0.5"])
      
    } else if (method == "jaccard") {
      ## Use your own Jaccard implementation instead of vegan::vegdist().
      ##
      ## binary = TRUE gives ordinary presence/absence Jaccard.
      ## If you want quantitative Jaccard instead, change TRUE to FALSE.
      safe_my_jaccard_distance(otu_mat, binary = TRUE)
      
    } else if (method %in% AITCHISON_METHODS) {
      pc <- DISTANCE_PSEUDOCOUNT[[method]]
      safe_aitchison_distance(otu_mat, pseudocount = pc)
      
    } else if (method %in% ILR_EUCLIDEAN_METHODS) {
      pc <- DISTANCE_PSEUDOCOUNT[[method]]
      safe_ILR_euclidean_distance(otu_mat, pseudocount = pc)
      
    } else if (method == "unifrac_unweighted") {
      phyloseq::distance(physeq, method = "unifrac", weighted = FALSE)
      
    } else if (method == "unifrac_weighted_phyloseq") {
      phyloseq::distance(physeq, method = "unifrac", weighted = TRUE)
      
    } else {
      vegan::vegdist(otu_mat, method = method)
    }
  }, error = function(e) {
    NULL
  })
  
  if (is.matrix(out)) out <- as.dist(out)
  out
}

safe_adonis_p <- function(D, x, permutations = 999) {
  if (is.null(D)) return(NA_real_)
  if (any(!is.finite(as.numeric(D)))) return(NA_real_)
  
  out <- tryCatch({
    tab <- vegan::adonis2(D ~ x, permutations = permutations)
    as.numeric(tab$`Pr(>F)`[1])
  }, error = function(e) {
    NA_real_
  })
  out
}

## ============================================================
## Direct SHbeta group-difference permutation test
##
## HRIC::SHbeta(X) returns one scalar beta-diversity value for the
## collection of samples in X. For the binary simulated group variable x,
## compare the two group-specific SHbeta values using
##
##   T_SHbeta = | SHbeta(X_group1) - SHbeta(X_group0) |.
##
## The p-value is obtained by permuting the group labels while preserving
## the original group sizes. This test does not use adonis2. It tests whether
## the two groups differ in their within-group Simplex Hellinger beta
## diversity (dispersion), rather than testing a general centroid/location
## difference between the groups.
## ============================================================
SHBETA_TEST_METHOD <- "SHbeta"

validate_SHbeta_input <- function(counts_mat, group) {
  X <- as.matrix(counts_mat)
  
  if (nrow(X) < 2L) {
    stop("SHbeta permutation test requires at least two samples.")
  }
  if (ncol(X) < 2L) {
    stop("SHbeta permutation test requires at least two features.")
  }
  if (length(group) != nrow(X)) {
    stop("group must have one entry per row of counts_mat.")
  }
  if (any(!is.finite(X))) {
    stop("counts_mat contains NA, NaN, or Inf values.")
  }
  if (any(X < 0)) {
    stop("SHbeta requires non-negative values.")
  }
  if (any(rowSums(X) <= 0)) {
    stop("Every row must have a positive total for SHbeta.")
  }
  if (anyNA(group)) {
    stop("group contains NA values.")
  }
  
  g <- droplevels(as.factor(group))
  if (nlevels(g) != 2L) {
    stop("The direct SHbeta group-difference test requires exactly two groups.")
  }
  if (any(table(g) < 2L)) {
    stop("Each group must contain at least two samples.")
  }
  
  list(X = X, group = g)
}

## Direct package-based definition. This explicitly evaluates HRIC::SHbeta()
## separately within the two groups.
SHbeta_group_difference_stat_direct <- function(counts_mat, group) {
  checked <- validate_SHbeta_input(counts_mat, group)
  X <- checked$X
  g <- checked$group
  lev <- levels(g)
  
  beta_by_group <- vapply(lev, function(level_g) {
    idx <- which(g == level_g)
    HRIC::SHbeta(X[idx, , drop = FALSE])
  }, numeric(1))
  
  abs(beta_by_group[2L] - beta_by_group[1L])
}

## For computational efficiency, permutations use the algebraically identical
## SHbeta expression from HRIC coordinates computed once for the full dataset.
## HRIC is a row-wise transformation, so subsetting the precomputed coordinate
## matrix is exactly equivalent to applying HRIC() to each permuted group.
SHbeta_from_Z <- function(Z_group, A_p_sq) {
  Z_bar <- colMeans(Z_group)
  dev <- sweep(Z_group, 2L, Z_bar, FUN = "-")
  mean(rowSums(dev^2)) / A_p_sq
}

SHbeta_group_difference_stat_from_Z <- function(Z, group, A_p_sq) {
  g <- droplevels(as.factor(group))
  lev <- levels(g)
  
  beta_by_group <- vapply(lev, function(level_g) {
    idx <- which(g == level_g)
    SHbeta_from_Z(Z[idx, , drop = FALSE], A_p_sq = A_p_sq)
  }, numeric(1))
  
  abs(beta_by_group[2L] - beta_by_group[1L])
}

safe_SHbeta_permutation_p <- function(
    counts_mat,
    group,
    permutations = 999,
    seed = 1L,
    equality_tol = 1e-10,
    permutation_tol = 1e-12
) {
  out <- tryCatch({
    if (!is.finite(permutations) || permutations < 1L) {
      stop("permutations must be a positive integer.")
    }
    permutations <- as.integer(permutations)
    
    checked <- validate_SHbeta_input(counts_mat, group)
    X <- checked$X
    g <- checked$group
    
    Z <- HRIC::HRIC(X)
    if (any(!is.finite(Z))) {
      stop("HRIC coordinates contain non-finite values.")
    }
    
    p <- ncol(Z)
    A_p_sq <- asin(sqrt(1 - 1 / p))^2
    if (!is.finite(A_p_sq) || A_p_sq <= 0) {
      stop("Invalid SHbeta normalization constant.")
    }
    
    observed_fast <- SHbeta_group_difference_stat_from_Z(Z, g, A_p_sq)
    observed_direct <- SHbeta_group_difference_stat_direct(X, g)
    
    if (!isTRUE(all.equal(
      observed_fast,
      observed_direct,
      tolerance = equality_tol,
      check.attributes = FALSE
    ))) {
      stop(sprintf(
        paste0(
          "Fast and direct SHbeta group-difference statistics disagree: ",
          "fast=%.16g, direct=%.16g"
        ),
        observed_fast, observed_direct
      ))
    }
    
    set.seed(seed)
    perm_stats <- numeric(permutations)
    for (b in seq_len(permutations)) {
      g_perm <- sample(g, size = length(g), replace = FALSE)
      perm_stats[b] <- SHbeta_group_difference_stat_from_Z(
        Z = Z,
        group = g_perm,
        A_p_sq = A_p_sq
      )
    }
    
    (1 + sum(perm_stats >= observed_fast - permutation_tol)) /
      (permutations + 1)
  }, error = function(e) {
    NA_real_
  })
  
  out
}


## ============================================================
## HRIC one-way trace-MANOVA decomposition
##
## Let Z_i = HRIC(X_i), Z_bar be the grand coordinate mean, and Z_bar_g
## be the coordinate mean in group g. The trace sums of squares are
##
##   SS_between = sum_g n_g ||Z_bar_g - Z_bar||^2,
##   SS_within  = sum_g sum_{i in g} ||Z_i - Z_bar_g||^2,
##   SS_total   = sum_i ||Z_i - Z_bar||^2.
##
## The reported F statistic is the Euclidean trace / pseudo-F statistic
##
##   F_HRIC = (SS_between / df_between) / (SS_within / df_within).
##
## This is the MANOVA-style table naturally associated with Euclidean HRIC
## coordinates. It is not a Wilks/Pillai/Roy classical MANOVA statistic.
## The existing AHC_euclidean adonis2 test supplies the corresponding
## permutation p-value; the table below is stored descriptively for trends.
## ============================================================
HRIC_MANOVA_SUMMARY_COLUMNS <- c(
  "shbeta_total",
  "shbeta_group0",
  "shbeta_group1",
  "shbeta_abs_diff",
  "shbeta_weighted_within",
  "shbeta_between",
  "hric_ss_between",
  "hric_ss_within",
  "hric_ss_total",
  "hric_df_between",
  "hric_df_within",
  "hric_df_total",
  "hric_ms_between",
  "hric_ms_within",
  "hric_manova_F",
  "hric_r2_between",
  "hric_n_group0",
  "hric_n_group1",
  "hric_coordinate_dimension"
)

empty_HRIC_manova_summary <- function() {
  as.data.frame(
    as.list(setNames(
      rep(NA_real_, length(HRIC_MANOVA_SUMMARY_COLUMNS)),
      HRIC_MANOVA_SUMMARY_COLUMNS
    )),
    stringsAsFactors = FALSE
  )
}

safe_HRIC_manova_summary <- function(
    counts_mat,
    group,
    equality_tol = 1e-9
) {
  tryCatch({
    checked <- validate_SHbeta_input(counts_mat, group)
    X <- checked$X
    g <- checked$group
    lev <- levels(g)

    Z <- HRIC::HRIC(X)
    if (any(!is.finite(Z))) {
      stop("HRIC coordinates contain non-finite values.")
    }

    n <- nrow(Z)
    q <- ncol(Z)
    G <- nlevels(g)
    n_by_group <- as.numeric(table(g)[lev])

    Z_bar <- colMeans(Z)
    total_dev <- sweep(Z, 2L, Z_bar, FUN = "-")
    ss_total <- sum(total_dev^2)

    ss_within <- 0
    ss_between <- 0
    beta_by_group <- numeric(G)

    for (k in seq_along(lev)) {
      idx <- which(g == lev[k])
      Z_k <- Z[idx, , drop = FALSE]
      Z_bar_k <- colMeans(Z_k)

      within_dev_k <- sweep(Z_k, 2L, Z_bar_k, FUN = "-")
      ss_within <- ss_within + sum(within_dev_k^2)
      ss_between <- ss_between + length(idx) * sum((Z_bar_k - Z_bar)^2)

      ## Save the group beta diversities using the package function directly.
      beta_by_group[k] <- HRIC::SHbeta(X[idx, , drop = FALSE])
    }

    ## Remove only negligible negative roundoff.
    if (ss_between < 0 && abs(ss_between) <= equality_tol) ss_between <- 0
    if (ss_within < 0 && abs(ss_within) <= equality_tol) ss_within <- 0
    if (ss_total < 0 && abs(ss_total) <= equality_tol) ss_total <- 0

    scale_ref <- max(1, abs(ss_total), abs(ss_between) + abs(ss_within))
    if (abs(ss_total - ss_between - ss_within) > equality_tol * scale_ref) {
      stop(sprintf(
        paste0(
          "HRIC SS decomposition failed: total=%.16g, ",
          "between+within=%.16g"
        ),
        ss_total,
        ss_between + ss_within
      ))
    }

    df_between <- G - 1L
    df_within <- n - G
    df_total <- n - 1L

    ms_between <- if (df_between > 0L) ss_between / df_between else NA_real_
    ms_within <- if (df_within > 0L) ss_within / df_within else NA_real_

    manova_F <- if (
      is.finite(ms_between) && is.finite(ms_within) && ms_within > 0
    ) {
      ms_between / ms_within
    } else {
      NA_real_
    }

    r2_between <- if (is.finite(ss_total) && ss_total > 0) {
      ss_between / ss_total
    } else {
      NA_real_
    }

    A_p_sq <- asin(sqrt(1 - 1 / q))^2
    if (!is.finite(A_p_sq) || A_p_sq <= 0) {
      stop("Invalid SHbeta normalization constant.")
    }

    shbeta_total <- HRIC::SHbeta(X)
    shbeta_weighted_within <- sum((n_by_group / n) * beta_by_group)
    shbeta_between <- shbeta_total - shbeta_weighted_within

    ## Check the exact correspondence between beta diversity and HRIC SS:
    ## total beta          = SS_total   / (n A_p^2)
    ## weighted within beta= SS_within  / (n A_p^2)
    ## between beta        = SS_between / (n A_p^2)
    beta_formula <- c(
      total = ss_total / (n * A_p_sq),
      weighted_within = ss_within / (n * A_p_sq),
      between = ss_between / (n * A_p_sq)
    )
    beta_direct <- c(
      total = shbeta_total,
      weighted_within = shbeta_weighted_within,
      between = shbeta_between
    )

    if (!isTRUE(all.equal(
      beta_formula,
      beta_direct,
      tolerance = equality_tol,
      check.attributes = FALSE
    ))) {
      stop(
        "Direct SHbeta values disagree with the HRIC sum-of-squares partition."
      )
    }

    data.frame(
      shbeta_total = shbeta_total,
      shbeta_group0 = beta_by_group[1L],
      shbeta_group1 = beta_by_group[2L],
      shbeta_abs_diff = abs(beta_by_group[2L] - beta_by_group[1L]),
      shbeta_weighted_within = shbeta_weighted_within,
      shbeta_between = max(shbeta_between, 0),
      hric_ss_between = ss_between,
      hric_ss_within = ss_within,
      hric_ss_total = ss_total,
      hric_df_between = df_between,
      hric_df_within = df_within,
      hric_df_total = df_total,
      hric_ms_between = ms_between,
      hric_ms_within = ms_within,
      hric_manova_F = manova_F,
      hric_r2_between = r2_between,
      hric_n_group0 = n_by_group[1L],
      hric_n_group1 = n_by_group[2L],
      hric_coordinate_dimension = q,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    empty_HRIC_manova_summary()
  })
}

build_HRIC_manova_long_table <- function(results_df) {
  id_columns <- intersect(
    c(
      "rep_id", "setting_id", "k", "signal_cluster", "cluster_ordered",
      "signal_mode", "beta", "is_null", "n_samples", "y_mean",
      "seed_sim", "seed_test"
    ),
    colnames(results_df)
  )

  rows <- vector("list", nrow(results_df))

  for (i in seq_len(nrow(results_df))) {
    r <- results_df[i, , drop = FALSE]
    ids <- r[rep(1L, 3L), id_columns, drop = FALSE]
    rownames(ids) <- NULL

    rows[[i]] <- cbind(
      ids,
      data.frame(
        source = c("Between", "Within", "Total"),
        df = c(r$hric_df_between, r$hric_df_within, r$hric_df_total),
        SS = c(r$hric_ss_between, r$hric_ss_within, r$hric_ss_total),
        MS = c(r$hric_ms_between, r$hric_ms_within, NA_real_),
        F = c(r$hric_manova_F, NA_real_, NA_real_),
        R2 = c(r$hric_r2_between, NA_real_, 1),
        shbeta_total = rep(r$shbeta_total, 3L),
        shbeta_group0 = rep(r$shbeta_group0, 3L),
        shbeta_group1 = rep(r$shbeta_group1, 3L),
        shbeta_abs_diff = rep(r$shbeta_abs_diff, 3L),
        shbeta_weighted_within = rep(r$shbeta_weighted_within, 3L),
        shbeta_between = rep(r$shbeta_between, 3L),
        stringsAsFactors = FALSE
      )
    )
  }

  do.call(rbind, rows)
}

run_permanova <- function(physeq, distance_methods, permutations = 999, seed = 1L) {
  x <- as.numeric(sample_data(physeq)$x)
  test_methods <- c(distance_methods, SHBETA_TEST_METHOD)
  pvals <- setNames(rep(NA_real_, length(test_methods)), test_methods)
  
  ## Keep every existing distance/PERMANOVA method and its seed exactly as in
  ## the pasted script. The new direct SHbeta test is appended afterward.
  for (j in seq_along(distance_methods)) {
    method <- distance_methods[j]
    set.seed(seed + j)
    D <- safe_distance(physeq, method)
    pvals[method] <- safe_adonis_p(D, x, permutations = permutations)
  }
  
  otu_mat <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
  
  pvals[SHBETA_TEST_METHOD] <- safe_SHbeta_permutation_p(
    counts_mat = otu_mat,
    group = x,
    permutations = permutations,
    seed = seed + length(distance_methods) + 1L
  )
  
  pvals
}

make_phyloseq_full <- function(counts_mat, tree, y_vec) {
  counts_mat <- as.matrix(counts_mat)
  stopifnot(nrow(counts_mat) == length(y_vec))
  
  if (is.null(rownames(counts_mat))) {
    rownames(counts_mat) <- paste0("sample_", seq_len(nrow(counts_mat)))
  }
  
  phyloseq(
    otu_table(counts_mat, taxa_are_rows = FALSE),
    sample_data(data.frame(
      x = as.numeric(y_vec),
      row.names = rownames(counts_mat)
    )),
    phy_tree(tree)
  )
}

make_metadata_matrix <- function(y_vec) {
  y_vec <- as.integer(y_vec)
  out <- matrix(
    y_vec,
    ncol = 1,
    dimnames = list(paste0("sample_", seq_along(y_vec)), "Y")
  )
  out
}

make_spike_config <- function(feature_ids, signal_mode, effect_size) {
  props <- switch(
    signal_mode,
    abundance = "abundance",
    prevalence = "prevalence",
    both = c("abundance", "prevalence"),
    stop("Unknown signal_mode: ", signal_mode)
  )
  
  do.call(rbind, lapply(props, function(prop) {
    data.frame(
      metadata_datum = 1,
      feature_spiked = feature_ids,
      associated_property = prop,
      effect_size = effect_size,
      stringsAsFactors = FALSE
    )
  }))
}

safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (length(x) <= 1L || all(is.na(x))) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_quantile <- function(x, prob) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE, type = 7))
}

upper_tri_values <- function(mat) {
  if (nrow(mat) <= 1L) return(numeric(0))
  mat[upper.tri(mat, diag = FALSE)]
}

get_boxplot_summary <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(data.frame(
      within_cop_n_pair = 0L,
      within_cop_mean = NA_real_,
      within_cop_sd = NA_real_,
      within_cop_min = NA_real_,
      within_cop_q1 = NA_real_,
      within_cop_median = NA_real_,
      within_cop_q3 = NA_real_,
      within_cop_max = NA_real_,
      within_cop_iqr = NA_real_,
      within_cop_lower_whisker = NA_real_,
      within_cop_upper_whisker = NA_real_,
      within_cop_n_outlier = 0L,
      stringsAsFactors = FALSE
    ))
  }
  
  bp <- boxplot.stats(x, coef = 1.5, do.conf = FALSE, do.out = TRUE)
  
  data.frame(
    within_cop_n_pair = length(x),
    within_cop_mean = mean(x),
    within_cop_sd = if (length(x) > 1L) sd(x) else NA_real_,
    within_cop_min = min(x),
    within_cop_q1 = safe_quantile(x, 0.25),
    within_cop_median = safe_quantile(x, 0.50),
    within_cop_q3 = safe_quantile(x, 0.75),
    within_cop_max = max(x),
    within_cop_iqr = IQR(x),
    within_cop_lower_whisker = bp$stats[1],
    within_cop_upper_whisker = bp$stats[5],
    within_cop_n_outlier = length(bp$out),
    stringsAsFactors = FALSE
  )
}

## -------------------- fixed manuscript-style setup --------------------
data("throat.otu.tab")
data("throat.tree")

counts <- as.matrix(throat.otu.tab)
common <- intersect(colnames(counts), throat.tree$tip.label)
counts <- counts[, common, drop = FALSE]
tree0  <- ape::drop.tip(throat.tree, setdiff(throat.tree$tip.label, colnames(counts)))

keep <- colSums(counts > 0) > 1
counts_f <- counts[, keep, drop = FALSE]
tree0    <- ape::drop.tip(tree0, setdiff(tree0$tip.label, colnames(counts_f)))

fit <- readRDS(FIT_RDS)

est_abs_mean <- get_est_mean_abs(fit)
missing_otus <- setdiff(colnames(counts_f), names(est_abs_mean))
if (length(missing_otus) > 0) {
  stop(
    "These OTUs are in counts_f but missing from est_abs_mean: ",
    paste(head(missing_otus, 20), collapse = ", ")
  )
}
est_abs_mean <- est_abs_mean[colnames(counts_f)]

## template-level relative abundance and prevalence
counts_f_tss <- row_tss(counts_f)
otu_mean_rel_template <- colMeans(counts_f_tss, na.rm = TRUE)
otu_mean_prev_template <- colMeans(counts_f > 0, na.rm = TRUE)

## -------------------- community distance matrix reused across K --------------------
cop <- ape::cophenetic.phylo(tree0)
cop <- cop[colnames(counts_f), colnames(counts_f)]

distance_methods <- c(
  "GLaD_0.5_weighted",
  "GLaD_0.5_unweighted",
  "GLaD_1_weighted",
  "GLaD_1_unweighted",
  "AHC_euclidean",
  "euclidean",
  "bray",
  "jaccard",
  AITCHISON_METHODS,
  ILR_EUCLIDEAN_METHODS,
  "unifrac_unweighted",
  "unifrac_weighted_phyloseq",
  "GUniFrac"
)

## Existing distance methods above are unchanged. The direct SHbeta
## group-difference permutation test is added as p_SHbeta and does not use adonis2.
test_methods <- c(distance_methods, SHBETA_TEST_METHOD)

total_settings <- sum(K_COMM_GRID) * length(beta_grid) * length(signal_mode_grid)

pam_list              <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
rank_map_list         <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
abund_sum_list        <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
n_otus_list           <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
settings_list         <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))
cluster_summary_list  <- setNames(vector("list", length(K_COMM_GRID)), as.character(K_COMM_GRID))

cat(sprintf(
  "[START] rep_id=%d | total_settings=%d | K grid=%s | N=%d | adonis_perm=%d\n",
  rep_id, total_settings, paste(K_COMM_GRID, collapse = ","), N_SAMPLES, ADONIS_PERM
))

cat("[DISTANCE METHODS]\n")
cat(paste(distance_methods, collapse = ", "), "\n")
cat("[ADDITIONAL SHBETA TEST]\n")
cat(SHBETA_TEST_METHOD, "\n")

cat("[PSEUDOCOUNT METHODS]\n")
print(DISTANCE_PSEUDOCOUNT)

rows_out <- vector("list", total_settings)
row_idx <- 0L

for (K_COMM in K_COMM_GRID) {
  
  ## -------------------- community definition --------------------
  pam_k <- polar_cluster(
    dist_matrix = as.matrix(cop),
    k = K_COMM,
    max_cluster_size = Inf,
    first_pole_method = "median",
    subsequent_pole_quantile = 1,
    assignment_method = "balanced",
    outlier_filter = TRUE,
    outlier_threshold = 1.5
  )$clustering
  
  abund_sum <- tapply(est_abs_mean[names(pam_k)], pam_k, sum, na.rm = TRUE)
  n_otus    <- tapply(pam_k, pam_k, length)
  
  ord_labels <- names(sort(abund_sum, decreasing = TRUE))
  rank_map   <- setNames(seq_along(ord_labels), ord_labels)
  
  ## silhouette
  sil_width_by_otu <- rep(NA_real_, length(pam_k))
  names(sil_width_by_otu) <- names(pam_k)
  sil_width_by_otu[] <- tryCatch({
    sil_obj <- cluster::silhouette(as.integer(pam_k), dist = as.dist(cop))
    as.numeric(sil_obj[, "sil_width"])
  }, error = function(e) {
    rep(NA_real_, length(pam_k))
  })
  names(sil_width_by_otu) <- names(pam_k)
  
  ## cluster summary table for this K
  cluster_summary_k <- vector("list", K_COMM)
  
  for (g in seq_len(K_COMM)) {
    otus_g <- names(pam_k)[pam_k == g]
    cop_g  <- cop[otus_g, otus_g, drop = FALSE]
    within_vals <- upper_tri_values(cop_g)
    bp_sum <- get_boxplot_summary(within_vals)
    
    cluster_summary_k[[g]] <- cbind(
      data.frame(
        k = K_COMM,
        signal_cluster = g,
        cluster_ordered = as.integer(rank_map[as.character(g)]),
        cluster_abund_sum_template = as.numeric(abund_sum[as.character(g)]),
        cluster_n_otus_template = as.integer(n_otus[as.character(g)]),
        cluster_avg_silhouette = safe_mean(sil_width_by_otu[otus_g]),
        cluster_mean_relative_abundance = safe_mean(otu_mean_rel_template[otus_g]),
        cluster_mean_prevalence = safe_mean(otu_mean_prev_template[otus_g]),
        stringsAsFactors = FALSE
      ),
      bp_sum
    )
  }
  
  cluster_summary_k <- do.call(rbind, cluster_summary_k)
  
  ## settings: one signal community at a time, binary Y, no confounder
  settings <- expand.grid(
    cluster_idx = seq_len(K_COMM),
    beta = beta_grid,
    signal_mode = signal_mode_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  pam_list[[as.character(K_COMM)]]             <- pam_k
  rank_map_list[[as.character(K_COMM)]]        <- rank_map
  abund_sum_list[[as.character(K_COMM)]]       <- abund_sum
  n_otus_list[[as.character(K_COMM)]]          <- n_otus
  settings_list[[as.character(K_COMM)]]        <- settings
  cluster_summary_list[[as.character(K_COMM)]] <- cluster_summary_k
  
  cat(sprintf("[K=%d] settings=%d\n", K_COMM, nrow(settings)))
  
  for (s in seq_len(nrow(settings))) {
    row_idx <- row_idx + 1L
    job <- settings[s, , drop = FALSE]
    
    cluster_idx_chr <- as.character(job$cluster_idx)
    C_otus <- names(pam_k)[pam_k == job$cluster_idx]
    
    cluster_info <- cluster_summary_k[cluster_summary_k$signal_cluster == job$cluster_idx, , drop = FALSE]
    
    if (!length(C_otus)) {
      rows_out[[row_idx]] <- cbind(
        data.frame(
          rep_id = rep_id,
          setting_id = s,
          k = K_COMM,
          signal_cluster = job$cluster_idx,
          signal_mode = as.character(job$signal_mode),
          beta = job$beta,
          n_samples = NA_integer_,
          y_mean = NA_real_,
          seed_sim = NA_integer_,
          seed_test = NA_integer_,
          runtime_sec = NA_real_,
          empty_HRIC_manova_summary(),
          matrix(
            NA_real_,
            nrow = 1,
            ncol = length(test_methods),
            dimnames = list(NULL, paste0("p_", test_methods))
          ),
          stringsAsFactors = FALSE
        ),
        cluster_info[, setdiff(colnames(cluster_info), c("k", "signal_cluster")), drop = FALSE]
      )
      next
    }
    
    cluster_ordered <- as.integer(rank_map[cluster_idx_chr])
    cluster_abund_sum_template <- as.numeric(abund_sum[cluster_idx_chr])
    cluster_n_otus_template <- as.integer(n_otus[cluster_idx_chr])
    
    seed_sim  <- BASE_SEED + (rep_id - 1L) * 10000L + row_idx
    seed_test <- BASE_SEED + 2000000L + (rep_id - 1L) * 10000L + row_idx
    
    cat(sprintf(
      "[rep %03d | %3d/%3d] K=%02d cluster=%02d signal=%s beta=%.2f\n",
      rep_id, row_idx, total_settings, K_COMM, job$cluster_idx,
      as.character(job$signal_mode), job$beta
    ))
    
    ## ---------- create balanced-binomial Y first, then spike signal during simulation ----------
    set.seed(seed_sim + 500000L)
    Y <- rbinom(N_SAMPLES, 1, 0.5)
    metadata_matrix <- make_metadata_matrix(Y)
    spike_config <- make_spike_config(
      feature_ids = C_otus,
      signal_mode = as.character(job$signal_mode),
      effect_size = job$beta
    )
    
    set.seed(seed_sim)
    sim <- SparseDOSSA2::SparseDOSSA2(
      template = fit,
      n_sample = N_SAMPLES,
      new_features = FALSE,
      spike_metadata = spike_config,
      metadata_matrix = metadata_matrix,
      verbose = FALSE
    )
    
    A_abs      <- t(sim$simulated_matrices$a_spiked)
    counts_sim <- t(sim$simulated_data)
    
    tree <- ape::drop.tip(tree0, setdiff(tree0$tip.label, colnames(A_abs)))
    A_abs      <- A_abs[, tree$tip.label, drop = FALSE]
    counts_sim <- counts_sim[, tree$tip.label, drop = FALSE]
    
    libsize <- rowSums(counts_sim)
    keep_samp <- is.finite(libsize) & libsize > 0
    if (!all(keep_samp)) {
      A_abs      <- A_abs[keep_samp, , drop = FALSE]
      counts_sim <- counts_sim[keep_samp, , drop = FALSE]
      Y          <- Y[keep_samp]
    }
    
    counts_sim_tss_full <- row_tss(counts_sim)
    
    C_curr <- intersect(C_otus, colnames(A_abs))
    if (!length(C_curr)) {
      rows_out[[row_idx]] <- cbind(
        data.frame(
          rep_id = rep_id,
          setting_id = s,
          k = K_COMM,
          signal_cluster = job$cluster_idx,
          signal_mode = as.character(job$signal_mode),
          beta = job$beta,
          n_samples = 0L,
          y_mean = NA_real_,
          seed_sim = seed_sim,
          seed_test = seed_test,
          runtime_sec = NA_real_,
          empty_HRIC_manova_summary(),
          matrix(
            NA_real_,
            nrow = 1,
            ncol = length(test_methods),
            dimnames = list(NULL, paste0("p_", test_methods))
          ),
          stringsAsFactors = FALSE
        ),
        cluster_info[, setdiff(colnames(cluster_info), c("k", "signal_cluster")), drop = FALSE]
      )
      next
    }
    
    ## ---------- whole OTU table beta-diversity tests ----------
    physeq_full <- make_phyloseq_full(counts_mat = counts_sim, tree = tree, y_vec = Y)
    
    S <- rowSums(A_abs[, C_curr, drop = FALSE])
    
    t0 <- proc.time()[["elapsed"]]
    pvals <- run_permanova(
      physeq = physeq_full,
      distance_methods = distance_methods,
      permutations = ADONIS_PERM,
      seed = seed_test
    )

    ## Descriptive HRIC beta-diversity and trace-MANOVA quantities.
    ## This does not change any existing test or p-value.
    hric_manova_summary <- safe_HRIC_manova_summary(
      counts_mat = counts_sim,
      group = Y
    )

    t1 <- proc.time()[["elapsed"]]
    
    rows_out[[row_idx]] <- cbind(
      data.frame(
        rep_id = rep_id,
        setting_id = s,
        k = K_COMM,
        signal_cluster = job$cluster_idx,
        cluster_ordered = cluster_ordered,
        cluster_abund_sum_template = cluster_abund_sum_template,
        cluster_n_otus_template = cluster_n_otus_template,
        signal_mode = as.character(job$signal_mode),
        beta = job$beta,
        is_null = as.integer(job$beta == 0),
        n_samples = nrow(counts_sim),
        y_mean = mean(Y),
        signal_mean_abs = mean(S, na.rm = TRUE),
        signal_sd_abs = sd(S, na.rm = TRUE),
        seed_sim = seed_sim,
        seed_test = seed_test,
        runtime_sec = (t1 - t0),
        hric_manova_summary,
        as.data.frame(as.list(setNames(as.numeric(pvals), paste0("p_", names(pvals))))),
        stringsAsFactors = FALSE
      ),
      cluster_info[, c(
        "cluster_avg_silhouette",
        "cluster_mean_relative_abundance",
        "cluster_mean_prevalence",
        "within_cop_n_pair",
        "within_cop_mean",
        "within_cop_sd",
        "within_cop_min",
        "within_cop_q1",
        "within_cop_median",
        "within_cop_q3",
        "within_cop_max",
        "within_cop_iqr",
        "within_cop_lower_whisker",
        "within_cop_upper_whisker",
        "within_cop_n_outlier"
      ), drop = FALSE]
    )
  }
}

res_df <- do.call(rbind, rows_out[seq_len(row_idx)])
hric_manova_table_df <- build_HRIC_manova_long_table(res_df)

outfile_csv <- file.path(
  OUTDIR,
  "raw",
  sprintf("beta_div_sparse12_shbeta_hric_manova_rep%03d.csv", rep_id)
)
outfile_manova_csv <- file.path(
  OUTDIR,
  "raw",
  sprintf("beta_div_sparse12_hric_manova_table_rep%03d.csv", rep_id)
)
outfile_rds <- file.path(
  OUTDIR,
  "raw",
  sprintf("beta_div_sparse12_shbeta_hric_manova_rep%03d.rds", rep_id)
)

write.csv(res_df, outfile_csv, row.names = FALSE)
write.csv(hric_manova_table_df, outfile_manova_csv, row.names = FALSE)

saveRDS(
  list(
    results = res_df,
    hric_manova_table = hric_manova_table_df,
    settings = settings_list[["12"]],
    settings_by_k = settings_list,
    pam12 = pam_list[["12"]],
    pam_by_k = pam_list,
    rank_map = rank_map_list[["12"]],
    rank_map_by_k = rank_map_list,
    abund_sum = abund_sum_list[["12"]],
    abund_sum_by_k = abund_sum_list,
    n_otus = n_otus_list[["12"]],
    n_otus_by_k = n_otus_list,
    cluster_summary12 = cluster_summary_list[["12"]],
    cluster_summary_by_k = cluster_summary_list,
    K_COMM_grid = K_COMM_GRID,
    distance_methods = distance_methods,
    test_methods = test_methods,
    shbeta_test_method = SHBETA_TEST_METHOD,
    shbeta_test_definition = paste0(
      "abs[SHbeta(X_group1) - SHbeta(X_group0)], ",
      "tested by label permutation without adonis2"
    ),
    hric_manova_definition = paste0(
      "Trace MANOVA on HRIC coordinates: ",
      "F=(SS_between/df_between)/(SS_within/df_within); ",
      "stored descriptively without replacing existing p-values"
    ),
    hric_manova_summary_columns = HRIC_MANOVA_SUMMARY_COLUMNS,
    distance_pseudocounts = DISTANCE_PSEUDOCOUNT,
    aitchison_methods = AITCHISON_METHODS,
    ilr_euclidean_methods = ILR_EUCLIDEAN_METHODS,
    signal_mode_grid = signal_mode_grid,
    rep_id = rep_id,
    n_samples = N_SAMPLES,
    adonis_perm = ADONIS_PERM
  ),
  outfile_rds
)

cat(sprintf(
  "[DONE] rep_id=%d | wrote:\n  %s\n  %s\n  %s\n",
  rep_id,
  outfile_csv,
  outfile_manova_csv,
  outfile_rds
))