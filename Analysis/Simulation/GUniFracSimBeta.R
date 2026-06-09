############################################################
## MSeqU-style simulation for beta-diversity method comparison
##
## Main design:
## 1. Keep the first MSeqU-style simulation framework.
## 2. Replace OTU-level scattered signal by polar-cluster community signal.
## 3. Test beta diversity using PERMANOVA: vegan::adonis2(D ~ x).
## 4. Compare all methods from the second pasted beta-diversity code:
##      GLaD_0.5_weighted
##      GLaD_0.5_unweighted
##      GLaD_1_weighted
##      GLaD_1_unweighted
##      euclidean
##      bray
##      jaccard
##      aitchison
##      unifrac_unweighted
##      unifrac_weighted_phyloseq
##      GUniFrac
##    plus:
##      AHC_euclidean
##
## Primary recommendation:
##   No confounder in this run, because the beta-diversity test is D ~ x.
##
## Output:
##   1. replicate-level p-values
##   2. setting-level power/type-I summaries
##   3. method-level overall summaries
############################################################


## ============================================================
## 1. Package path and packages
## ============================================================

.libPaths(c("/home/zhang.16383/RPackage_r451", .libPaths()))

user_lib <- "/home/zhang.16383/RPackage_r451"

if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
}

options(repos = c(CRAN = "https://cloud.r-project.org"))

install_and_load <- function(pkg, lib) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = lib)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

packages <- c(
  "MASS",
  "dplyr",
  "tibble",
  "dirmult",
  "ape",
  "cluster",
  "phyloseq",
  "GUniFrac",
  "vegan",
  "MiSPU"
)

for (pkg in packages) {
  install_and_load(pkg, user_lib)
}

## GLaD is optional because it may be your local package.
## The code will return NA for GLaD methods if GLaD is unavailable.
has_GLaD <- requireNamespace("GLaD", quietly = TRUE)


## ============================================================
## 2. Project path and input loading
## ============================================================

project_dir <- Sys.getenv("PROJECT_DIR", "/home/zhang.16383/AHC")
outdir <- file.path(project_dir, "beta_div_msequ_out")

dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "raw"), recursive = TRUE, showWarnings = FALSE)

para_file <- file.path(project_dir, "para1.RData")

if (!file.exists(para_file)) {
  if (file.exists("para1.RData")) {
    para_file <- "para1.RData"
  } else {
    stop("Cannot find para1.RData. Put para1.RData under project_dir or current working directory.")
  }
}

load(para_file)

if (!exists("para1")) {
  stop("para1.RData was loaded, but object 'para1' was not found.")
}

## polar_cluster.R should contain your polar_cluster() function.
## This follows the second simulation workflow.
if (file.exists("polar_cluster.R")) {
  source("polar_cluster.R")
} else if (file.exists(file.path(project_dir, "polar_cluster.R"))) {
  source(file.path(project_dir, "polar_cluster.R"))
} else {
  stop("Cannot find polar_cluster.R. Put it in the current directory or project_dir.")
}

if (!exists("polar_cluster", mode = "function")) {
  stop("polar_cluster() was not found after sourcing polar_cluster.R.")
}

## Load throat tree as default phylogenetic tree.
## This is used by UniFrac/GUniFrac/GLaD and for polar clustering.
data("throat.tree", package = "MiSPU")


## ============================================================
## 3. General helpers
## ============================================================

safe_scale <- function(x) {
  x <- as.matrix(x)
  out <- scale(x)
  out[!is.finite(out)] <- 0
  out
}

safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (length(x) <= 1L || all(is.na(x))) return(NA_real_)
  sd(x, na.rm = TRUE)
}

row_tss <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  
  rs <- rowSums(x)
  out <- matrix(
    0,
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = dimnames(x)
  )
  
  keep <- is.finite(rs) & rs > 0
  
  if (any(keep)) {
    out[keep, ] <- sweep(x[keep, , drop = FALSE], 1, rs[keep], FUN = "/")
  }
  
  out
}

normalize_rows_allow_all_zero <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  
  n <- nrow(X)
  p <- ncol(X)
  
  rs <- rowSums(X)
  
  Pi <- matrix(NA_real_, nrow = n, ncol = p)
  
  positive_rows <- rs > 0
  zero_rows <- rs == 0
  
  Pi[positive_rows, ] <- sweep(
    X[positive_rows, , drop = FALSE],
    1,
    rs[positive_rows],
    FUN = "/"
  )
  
  ## If a sample has total count zero, map it to the uniform composition.
  if (any(zero_rows)) {
    Pi[zero_rows, ] <- matrix(
      1 / p,
      nrow = sum(zero_rows),
      ncol = p
    )
  }
  
  rownames(Pi) <- rownames(X)
  colnames(Pi) <- colnames(X)
  
  Pi
}

get_tree_for_otus <- function(tree, otu_ids) {
  
  if (is.null(tree) || !inherits(tree, "phylo")) {
    stop("tree must be a valid phylo object.")
  }
  
  missing_otus <- setdiff(otu_ids, tree$tip.label)
  
  if (length(missing_otus) > 0) {
    stop(
      "The tree is missing some OTUs used in the simulation. First missing OTUs: ",
      paste(head(missing_otus, 20), collapse = ", ")
    )
  }
  
  drop_tips <- setdiff(tree$tip.label, otu_ids)
  
  if (length(drop_tips) > 0) {
    tree_use <- ape::drop.tip(tree, drop_tips)
  } else {
    tree_use <- tree
  }
  
  ## Keep tree order synchronized later by subsetting count matrix to tree$tip.label.
  tree_use
}


## ============================================================
## 4. Original MSeqU-style helper functions
## ============================================================

rdirichlet.m <- function(alpha) {
  Gam <- matrix(
    rgamma(length(alpha), shape = alpha),
    nrow = nrow(alpha),
    ncol = ncol(alpha)
  )
  t(t(Gam) / colSums(Gam))
}

EstPara <- function(ref.otu.tab) {
  
  if (is.null(rownames(ref.otu.tab))) {
    rownames(ref.otu.tab) <- paste0("OTU", seq_len(nrow(ref.otu.tab)))
  }
  
  samplenames <- colnames(ref.otu.tab)
  taxnames <- rownames(ref.otu.tab)
  
  dirmult.paras <- dirmult::dirmult(t(ref.otu.tab))
  
  gamma <- dirmult.paras$gamma
  names(gamma) <- names(dirmult.paras$pi)
  
  ref.otu.tab <- sapply(
    seq_len(ncol(ref.otu.tab)),
    function(i) gamma + ref.otu.tab[, i]
  )
  
  ref.otu.tab.p <- rdirichlet.m(ref.otu.tab)
  
  colnames(ref.otu.tab.p) <- samplenames
  rownames(ref.otu.tab.p) <- taxnames
  
  ord <- order(rowMeans(ref.otu.tab.p), decreasing = TRUE)
  ref.otu.tab.p <- ref.otu.tab.p[ord, ]
  
  Si <- exp(rnorm(ncol(ref.otu.tab.p)))
  ref.otu.tab0 <- t(t(ref.otu.tab.p) * Si)
  
  colnames(ref.otu.tab0) <- colnames(ref.otu.tab.p)
  rownames(ref.otu.tab0) <- rownames(ref.otu.tab.p)
  
  list(mu = ref.otu.tab.p, ref.otu.tab = ref.otu.tab0)
}


## ============================================================
## 5. Revised SimulateMSeqU
##    Main revision:
##    true signal OTUs are selected by polar-cluster community.
## ============================================================

SimulateMSeqU <- function(
    para,
    nSam = 100,
    nOTU = 100,
    diff.otu.pct = 0.1,
    diff.otu.direct = c("balanced", "unbalanced"),
    diff.otu.mode = c("abundant", "rare", "mix", "user_specified", "cluster"),
    user_specified_otu = NULL,
    covariate.type = c("binary", "continuous"),
    grp.ratio = 1,
    covariate.eff.mean = 1,
    covariate.eff.sd = 0,
    confounder.type = c("none", "binary", "continuous", "both"),
    conf.cov.cor = 0.6,
    conf.diff.otu.pct = 0,
    conf.nondiff.otu.pct = 0,
    confounder.eff.mean = 0,
    confounder.eff.sd = 0,
    error.sd = 0,
    depth.mu = 10000,
    depth.theta = 5,
    depth.conf.factor = 0,
    cont.conf = 0,
    epsilon = NULL,
    
    ## Community/cluster-level signal arguments
    use.cluster.signal = TRUE,
    phylo.tree = NULL,
    K_COMM = NULL,
    signal.cluster = NULL,
    cluster.assignment = NULL,
    cluster.first_pole_method = "median",
    cluster.subsequent_pole_quantile = 1,
    cluster.assignment_method = "balanced",
    cluster.outlier_filter = TRUE,
    cluster.outlier_threshold = 1.5
) {
  
  diff.otu.direct <- match.arg(diff.otu.direct)
  diff.otu.mode <- match.arg(diff.otu.mode)
  covariate.type <- match.arg(covariate.type)
  confounder.type <- match.arg(confounder.type)
  
  model.paras <- para
  
  ## ----------------------------
  ## Select number of OTUs and samples
  ## ----------------------------
  
  ref.otu.tab <- model.paras$ref.otu.tab[seq_len(nOTU), , drop = FALSE]
  
  if (is.null(rownames(ref.otu.tab))) {
    rownames(ref.otu.tab) <- paste0("OTU", seq_len(nrow(ref.otu.tab)))
  }
  
  idx.otu <- rownames(ref.otu.tab)
  idx.sample <- colnames(model.paras$ref.otu.tab)[seq_len(nSam)]
  ref.otu.tab <- ref.otu.tab[, idx.sample, drop = FALSE]
  
  ## ----------------------------
  ## Confounder
  ## ----------------------------
  
  if (confounder.type == "none") {
    Z <- matrix(0, nrow = nSam, ncol = 1)
    colnames(Z) <- "Z0"
    confounder.eff.mean <- 0
    confounder.eff.sd <- 0
  }
  
  if (confounder.type == "continuous") {
    Z <- matrix(cont.conf, nrow = nSam, ncol = 1)
    colnames(Z) <- "Z1"
  }
  
  if (confounder.type == "binary") {
    Z <- matrix(
      c(rep(0, nSam %/% 2), rep(1, nSam - nSam %/% 2)),
      ncol = 1
    )
    colnames(Z) <- "Z1"
  }
  
  if (confounder.type == "both") {
    Z <- cbind(
      Z1 = rnorm(nSam),
      Z2 = c(rep(0, nSam %/% 2), rep(1, nSam - nSam %/% 2))
    )
  }
  
  rownames(Z) <- colnames(ref.otu.tab)
  
  ## ----------------------------
  ## Covariate of interest
  ##
  ## Important change:
  ## If there is no confounder, generate X independently.
  ## Otherwise, generate X correlated with Z.
  ## ----------------------------
  
  if (confounder.type == "none") {
    
    X_cont <- rnorm(nSam)
    
  } else {
    
    rho <- sqrt(conf.cov.cor^2 / (1 - conf.cov.cor^2))
    
    Z_summary <- safe_scale(Z) %*% rep(1, ncol(Z))
    
    if (is.null(epsilon)) {
      epsilon_vec <- rnorm(nSam)
    } else if (length(epsilon) == 1L) {
      epsilon_vec <- rnorm(nSam, mean = 0, sd = as.numeric(epsilon))
    } else if (length(epsilon) == nSam) {
      epsilon_vec <- as.numeric(epsilon)
    } else {
      stop("epsilon must be NULL, a scalar standard deviation, or a vector of length nSam.")
    }
    
    X_cont <- as.numeric(rho * safe_scale(Z_summary) + epsilon_vec)
  }
  
  if (covariate.type == "continuous") {
    X <- as.matrix(X_cont)
  }
  
  if (covariate.type == "binary") {
    cutoff <- quantile(X_cont, grp.ratio / (1 + grp.ratio), names = FALSE)
    X <- cbind(ifelse(X_cont <= cutoff, 0, 1))
  }
  
  colnames(X) <- "X"
  rownames(X) <- colnames(ref.otu.tab)
  
  group_tab <- table(as.numeric(X[, 1]))
  if (length(group_tab) != 2 || any(group_tab < 2)) {
    stop("Generated X does not have two valid groups. Check covariate generation.")
  }
  
  covariate.eff.mean1 <- covariate.eff.mean
  covariate.eff.mean2 <- covariate.eff.mean
  
  ## ----------------------------
  ## OTU-level covariate effects
  ## Original effect-generation logic is preserved.
  ## The only major change is which OTUs receive nonzero effects.
  ## ----------------------------
  
  if (diff.otu.direct == "balanced") {
    
    if (diff.otu.mode %in% c("abundant", "rare", "user_specified", "cluster")) {
      
      effect_vec <- sample(c(
        rnorm(floor(nOTU / 2), mean = -covariate.eff.mean2, sd = covariate.eff.sd),
        rnorm(nOTU - floor(nOTU / 2), mean = covariate.eff.mean2, sd = covariate.eff.sd)
      ))
      
    } else if (diff.otu.mode == "mix") {
      
      effect_vec <- c(
        sample(c(
          rnorm(floor(nOTU / 4), mean = -covariate.eff.mean1, sd = covariate.eff.sd),
          rnorm(floor(nOTU / 2) - floor(nOTU / 4), mean = covariate.eff.mean1, sd = covariate.eff.sd)
        )),
        sample(c(
          rnorm(floor((nOTU - floor(nOTU / 2)) / 2), mean = -covariate.eff.mean2, sd = covariate.eff.sd),
          rnorm(
            nOTU - floor(nOTU / 2) - floor((nOTU - floor(nOTU / 2)) / 2),
            mean = covariate.eff.mean2,
            sd = covariate.eff.sd
          )
        ))
      )
    }
  }
  
  if (diff.otu.direct == "unbalanced") {
    
    if (diff.otu.mode %in% c("abundant", "rare", "user_specified", "cluster")) {
      
      effect_vec <- rnorm(
        nOTU,
        mean = covariate.eff.mean2,
        sd = covariate.eff.sd
      )
      
    } else if (diff.otu.mode == "mix") {
      
      effect_vec <- c(
        sample(rnorm(
          floor(nOTU / 2),
          mean = covariate.eff.mean1,
          sd = covariate.eff.sd
        )),
        sample(rnorm(
          nOTU - floor(nOTU / 2),
          mean = covariate.eff.mean2,
          sd = covariate.eff.sd
        ))
      )
    }
  }
  
  eta.diff <- effect_vec %*% t(safe_scale(X))
  
  ## ----------------------------
  ## Confounder effects
  ## For the primary beta-diversity simulation, confounder.type = "none",
  ## so this becomes exactly zero.
  ## ----------------------------
  
  if (confounder.type == "none") {
    eta.conf <- matrix(0, nrow = nOTU, ncol = nSam)
    conf.otu.ind <- integer(0)
  } else {
    
    Z_summary <- safe_scale(Z) %*% rep(1, ncol(Z))
    
    conf_effect_vec <- sample(c(
      rnorm(floor(nOTU / 2), mean = -confounder.eff.mean, sd = confounder.eff.sd),
      rnorm(nOTU - floor(nOTU / 2), mean = confounder.eff.mean, sd = confounder.eff.sd)
    ))
    
    eta.conf <- conf_effect_vec %*% t(safe_scale(Z_summary))
  }
  
  ## ============================================================
  ## Choose true signal OTUs by polar-cluster community
  ## ============================================================
  
  otu.ord <- seq_len(nOTU)
  diff.otu.num <- round(diff.otu.pct * nOTU)
  
  pam_k <- NULL
  cluster_abund_sum <- NULL
  cluster_n_otus <- NULL
  cluster_rank_map <- NULL
  signal.cluster.used <- NA_integer_
  signal.cluster.ordered <- NA_integer_
  
  if (use.cluster.signal && diff.otu.mode != "user_specified") {
    
    if (is.null(K_COMM)) {
      if (!is.finite(diff.otu.pct) || diff.otu.pct <= 0 || diff.otu.pct > 1) {
        stop("When K_COMM is NULL, diff.otu.pct must be in (0, 1].")
      }
      K_COMM <- round(1 / diff.otu.pct)
      K_COMM <- max(1L, min(nOTU, as.integer(K_COMM)))
    }
    
    if (K_COMM < 1 || K_COMM > nOTU) {
      stop("K_COMM must be between 1 and nOTU.")
    }
    
    if (!is.null(cluster.assignment)) {
      
      pam_k <- cluster.assignment
      
      if (is.null(names(pam_k))) {
        stop("cluster.assignment must be a named vector with OTU IDs as names.")
      }
      
      missing_cluster_otus <- setdiff(idx.otu, names(pam_k))
      
      if (length(missing_cluster_otus) > 0) {
        stop(
          "cluster.assignment is missing some simulated OTUs: ",
          paste(head(missing_cluster_otus, 20), collapse = ", ")
        )
      }
      
      pam_k <- pam_k[idx.otu]
      
    } else {
      
      if (is.null(phylo.tree)) {
        if (!is.null(model.paras$tree)) {
          phylo.tree <- model.paras$tree
        } else {
          stop("use.cluster.signal = TRUE requires phylo.tree, cluster.assignment, or para$tree.")
        }
      }
      
      tree_use <- get_tree_for_otus(phylo.tree, idx.otu)
      cop <- ape::cophenetic.phylo(tree_use)
      cop <- cop[idx.otu, idx.otu, drop = FALSE]
      
      pam_k <- polar_cluster(
        dist_matrix = as.matrix(cop),
        k = K_COMM,
        max_cluster_size = Inf,
        first_pole_method = cluster.first_pole_method,
        subsequent_pole_quantile = cluster.subsequent_pole_quantile,
        assignment_method = cluster.assignment_method,
        outlier_filter = cluster.outlier_filter,
        outlier_threshold = cluster.outlier_threshold
      )$clustering
      
      pam_k <- pam_k[idx.otu]
    }
    
    otu_template_abund <- rowMeans(ref.otu.tab, na.rm = TRUE)
    names(otu_template_abund) <- idx.otu
    
    cluster_abund_sum <- tapply(
      otu_template_abund[names(pam_k)],
      pam_k,
      sum,
      na.rm = TRUE
    )
    
    cluster_n_otus <- tapply(pam_k, pam_k, length)
    
    ord_labels <- names(sort(cluster_abund_sum, decreasing = TRUE))
    cluster_rank_map <- setNames(seq_along(ord_labels), ord_labels)
    
    if (is.null(signal.cluster)) {
      
      if (diff.otu.mode == "abundant") {
        signal.cluster.used <- as.integer(names(which.max(cluster_abund_sum)))
      }
      
      if (diff.otu.mode == "rare") {
        signal.cluster.used <- as.integer(names(which.min(cluster_abund_sum)))
      }
      
      if (diff.otu.mode %in% c("mix", "cluster")) {
        signal.cluster.used <- sample(as.integer(names(cluster_abund_sum)), 1)
      }
      
    } else {
      
      signal.cluster.used <- as.integer(signal.cluster)
      
      if (!signal.cluster.used %in% as.integer(names(cluster_abund_sum))) {
        stop("signal.cluster is not one of the cluster labels produced by polar_cluster().")
      }
    }
    
    signal.cluster.ordered <- as.integer(
      cluster_rank_map[as.character(signal.cluster.used)]
    )
    
    diff.otu.ind <- which(pam_k == signal.cluster.used)
    
  } else {
    
    ## Original OTU-level selection fallback.
    if (diff.otu.mode == "user_specified") {
      if (is.null(user_specified_otu)) {
        stop("diff.otu.mode = 'user_specified' requires user_specified_otu.")
      }
      diff.otu.ind <- which(idx.otu %in% user_specified_otu)
    }
    
    if (diff.otu.mode == "mix") {
      diff.otu.ind <- sample(otu.ord, diff.otu.num)
    }
    
    if (diff.otu.mode == "abundant") {
      abundant_pool <- seq_len(round(length(otu.ord) / 4))
      diff.otu.ind <- sample(abundant_pool, diff.otu.num)
    }
    
    if (diff.otu.mode == "rare") {
      rare_start <- round(3 * length(otu.ord) / 4)
      rare_pool <- rare_start:length(otu.ord)
      diff.otu.ind <- sample(rare_pool, diff.otu.num)
    }
    
    if (diff.otu.mode == "cluster") {
      stop("diff.otu.mode = 'cluster' requires use.cluster.signal = TRUE.")
    }
  }
  
  if (length(diff.otu.ind) == 0) {
    stop("No signal OTUs were selected. Check cluster/tree/OTU names.")
  }
  
  ## ----------------------------
  ## Choose confounded OTUs if confounders are used.
  ## For primary beta-diversity simulation, this section is skipped.
  ## ----------------------------
  
  if (confounder.type != "none") {
    
    n_conf_diff <- round(nOTU * conf.diff.otu.pct)
    n_conf_nondiff <- round(nOTU * conf.nondiff.otu.pct)
    
    if (length(diff.otu.ind) >= n_conf_diff) {
      conf.otu.ind1 <- if (n_conf_diff > 0) {
        sample(diff.otu.ind, n_conf_diff)
      } else {
        integer(0)
      }
    } else {
      conf.otu.ind1 <- diff.otu.ind
    }
    
    nondiff_pool <- setdiff(seq_len(nOTU), diff.otu.ind)
    
    if (n_conf_nondiff > length(nondiff_pool)) {
      n_conf_nondiff <- length(nondiff_pool)
    }
    
    conf.otu.ind2 <- if (n_conf_nondiff > 0) {
      sample(nondiff_pool, n_conf_nondiff)
    } else {
      integer(0)
    }
    
    conf.otu.ind <- c(conf.otu.ind1, conf.otu.ind2)
    
  } else {
    
    conf.otu.ind <- integer(0)
  }
  
  ## ----------------------------
  ## Apply effects only to selected OTUs
  ## ----------------------------
  
  eta.diff[setdiff(seq_len(nOTU), diff.otu.ind), ] <- 0
  
  if (confounder.type != "none") {
    eta.conf[setdiff(seq_len(nOTU), conf.otu.ind), ] <- 0
  } else {
    eta.conf[,] <- 0
  }
  
  eta.error <- matrix(
    rnorm(nOTU * nSam, 0, error.sd),
    nrow = nOTU,
    ncol = nSam
  )
  
  eta.exp <- exp(t(eta.diff + eta.conf + eta.error))
  eta.exp <- eta.exp * t(ref.otu.tab)
  
  ref.otu.tab.prop <- eta.exp / rowSums(eta.exp)
  ref.otu.tab.prop <- t(ref.otu.tab.prop)
  
  ## ----------------------------
  ## Sequencing depth
  ## Primary beta-diversity setting uses depth.conf.factor = 0.
  ## ----------------------------
  
  nSeq <- MASS::rnegbin(
    nSam,
    mu = as.numeric(depth.mu * exp(safe_scale(X) * depth.conf.factor)),
    theta = depth.theta
  )
  
  nSeq[nSeq < 1] <- 1L
  
  otu.tab.sim <- sapply(
    seq_len(ncol(ref.otu.tab.prop)),
    function(i) {
      rmultinom(1, nSeq[i], ref.otu.tab.prop[, i])
    }
  )
  
  colnames(otu.tab.sim) <- rownames(eta.exp)
  rownames(otu.tab.sim) <- rownames(ref.otu.tab)
  
  diff.otu.ind.logical <- seq_len(nOTU) %in% diff.otu.ind
  conf.otu.ind.logical <- seq_len(nOTU) %in% conf.otu.ind
  
  list(
    otu.tab.sim = otu.tab.sim,
    covariate = X,
    confounder = Z,
    diff.otu.ind = diff.otu.ind.logical,
    otu.names = idx.otu,
    conf.otu.ind = conf.otu.ind.logical,
    
    use.cluster.signal = use.cluster.signal,
    K_COMM = K_COMM,
    cluster.assignment = pam_k,
    signal.cluster = signal.cluster.used,
    signal.cluster.ordered = signal.cluster.ordered,
    cluster.abund.sum = cluster_abund_sum,
    cluster.n.otus = cluster_n_otus,
    cluster.rank.map = cluster_rank_map,
    signal.otu.names = idx.otu[diff.otu.ind]
  )
}


## ============================================================
## 6. AHC transformation for beta-diversity distance
## ============================================================

ahc_matrix <- function(counts_mat) {
  
  counts_mat <- as.matrix(counts_mat)
  storage.mode(counts_mat) <- "numeric"
  
  Pi <- normalize_rows_allow_all_zero(counts_mat)
  
  k <- ncol(Pi)
  
  sqrt_Pi <- sqrt(Pi)
  sqrt_pi0 <- rep(1 / sqrt(k), k)
  
  c_pi <- as.vector(sqrt_Pi %*% sqrt_pi0)
  c_pi <- pmin(pmax(c_pi, 0), 1)
  
  s_pi <- sqrt(pmax(1 - c_pi^2, 0))
  
  scale_factor <- asin(s_pi) / s_pi
  scale_factor[!is.finite(scale_factor)] <- 1
  scale_factor[s_pi == 0] <- 1
  
  ahc_data <- (sqrt_Pi - outer(c_pi, sqrt_pi0)) * scale_factor
  
  rownames(ahc_data) <- rownames(counts_mat)
  colnames(ahc_data) <- colnames(counts_mat)
  
  ahc_data
}


## ============================================================
## 7. Phyloseq and beta-diversity test functions
## ============================================================

make_phyloseq_full <- function(counts_mat, tree, x_vec, Z = NULL) {
  
  counts_mat <- as.matrix(counts_mat)
  storage.mode(counts_mat) <- "numeric"
  
  if (is.null(rownames(counts_mat))) {
    rownames(counts_mat) <- paste0("sample_", seq_len(nrow(counts_mat)))
  }
  
  stopifnot(nrow(counts_mat) == length(x_vec))
  
  tree_use <- get_tree_for_otus(tree, colnames(counts_mat))
  
  ## Synchronize columns to tree tip order.
  counts_mat <- counts_mat[, tree_use$tip.label, drop = FALSE]
  
  meta <- data.frame(
    x = as.numeric(x_vec),
    row.names = rownames(counts_mat)
  )
  
  if (!is.null(Z)) {
    Z <- as.data.frame(Z)
    rownames(Z) <- rownames(counts_mat)
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    meta <- cbind(meta, Z)
  }
  
  phyloseq::phyloseq(
    phyloseq::otu_table(counts_mat, taxa_are_rows = FALSE),
    phyloseq::sample_data(meta),
    phyloseq::phy_tree(tree_use)
  )
}

safe_distance <- function(physeq, method) {
  
  otu_mat <- as(phyloseq::otu_table(physeq), "matrix")
  
  if (phyloseq::taxa_are_rows(physeq)) {
    otu_mat <- t(otu_mat)
  }
  
  out <- tryCatch({
    
    if (method == "AHC_euclidean") {
      
      ahc_dat <- ahc_matrix(otu_mat)
      stats::dist(ahc_dat, method = "euclidean")
      
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
      
      tree <- phyloseq::phy_tree(physeq)
      unifrac_arr <- GUniFrac::GUniFrac(
        otu_mat,
        tree,
        alpha = c(0, 0.5, 1)
      )$unifracs
      
      as.dist(unifrac_arr[, , "d_0.5"])
      
    } else if (method == "aitchison") {
      
      vegan::vegdist(otu_mat, method = "aitchison", pseudocount = 0.5)
      
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
  
  if (is.matrix(out)) {
    out <- as.dist(out)
  }
  
  out
}

safe_adonis_p <- function(D,
                          metadata,
                          permutations = 999,
                          adjust_confounders = FALSE) {
  
  if (is.null(D)) return(NA_real_)
  
  d_vec <- as.numeric(D)
  
  if (length(d_vec) == 0L || any(!is.finite(d_vec))) {
    return(NA_real_)
  }
  
  metadata <- as.data.frame(metadata)
  
  if (!("x" %in% colnames(metadata))) {
    stop("metadata must contain column x.")
  }
  
  out <- tryCatch({
    
    if (!adjust_confounders) {
      
      tab <- vegan::adonis2(D ~ x, data = metadata, permutations = permutations)
      as.numeric(tab$`Pr(>F)`[1])
      
    } else {
      
      z_cols <- grep("^Z", colnames(metadata), value = TRUE)
      
      ## Keep only nonconstant confounders.
      z_cols <- z_cols[sapply(metadata[, z_cols, drop = FALSE], function(v) {
        length(unique(v[is.finite(v)])) > 1
      })]
      
      if (length(z_cols) == 0L) {
        tab <- vegan::adonis2(D ~ x, data = metadata, permutations = permutations)
        as.numeric(tab$`Pr(>F)`[1])
      } else {
        rhs <- paste(c("x", z_cols), collapse = " + ")
        form <- as.formula(paste("D ~", rhs))
        tab <- vegan::adonis2(
          form,
          data = metadata,
          permutations = permutations,
          by = "margin"
        )
        as.numeric(tab["x", "Pr(>F)"])
      }
    }
    
  }, error = function(e) {
    NA_real_
  })
  
  out
}

run_permanova <- function(physeq,
                          distance_methods,
                          permutations = 999,
                          seed = 1L,
                          adjust_confounders = FALSE) {
  
  metadata <- as.data.frame(phyloseq::sample_data(physeq))
  
  pvals <- setNames(rep(NA_real_, length(distance_methods)), distance_methods)
  
  for (j in seq_along(distance_methods)) {
    
    method <- distance_methods[j]
    
    set.seed(seed + j)
    
    D <- safe_distance(physeq, method)
    
    pvals[method] <- safe_adonis_p(
      D = D,
      metadata = metadata,
      permutations = permutations,
      adjust_confounders = adjust_confounders
    )
  }
  
  pvals
}


## ============================================================
## 8. Prepare fixed tree and fixed community assignment
## ============================================================

base_seed <- as.integer(Sys.getenv("BASE_SEED", "12345"))
set.seed(base_seed)

n_simulation <- as.integer(Sys.getenv("N_SIMULATION", "100"))
nSam <- as.integer(Sys.getenv("N_SAMPLES", "100"))
nOTU <- as.integer(Sys.getenv("N_OTU", "100"))
significance_level <- as.numeric(Sys.getenv("ALPHA", "0.05"))
ADONIS_PERM <- as.integer(Sys.getenv("ADONIS_PERM", "999"))

## Primary beta-diversity setting:
## no confounder and no library-size confounding.
USE_CONFOUNDER <- FALSE
ADJUST_FOR_CONFOUNDERS <- FALSE

diff_otu_pct <- 0.1
K_COMM <- round(1 / diff_otu_pct)

idx.otu.global <- rownames(para1$ref.otu.tab)[seq_len(nOTU)]

if (is.null(idx.otu.global)) {
  stop("rownames(para1$ref.otu.tab) are required and must match the tree tip labels.")
}

tree_global <- get_tree_for_otus(throat.tree, idx.otu.global)

cop_global <- ape::cophenetic.phylo(tree_global)
cop_global <- cop_global[idx.otu.global, idx.otu.global, drop = FALSE]

pam_global <- polar_cluster(
  dist_matrix = as.matrix(cop_global),
  k = K_COMM,
  max_cluster_size = Inf,
  first_pole_method = "median",
  subsequent_pole_quantile = 1,
  assignment_method = "balanced",
  outlier_filter = TRUE,
  outlier_threshold = 1.5
)$clustering

pam_global <- pam_global[idx.otu.global]

cat("Fixed community assignment created.\n")
cat("K_COMM =", K_COMM, "\n")
cat("Cluster sizes:\n")
print(table(pam_global))


## ============================================================
## 9. Simulation settings
## ============================================================

distance_methods <- c(
  "AHC_euclidean",
  "GLaD_0.5_weighted",
  "GLaD_0.5_unweighted",
  "GLaD_1_weighted",
  "GLaD_1_unweighted",
  "euclidean",
  "bray",
  "jaccard",
  "aitchison",
  "unifrac_unweighted",
  "unifrac_weighted_phyloseq",
  "GUniFrac"
)

## Effect size 0 is the null setting for type I error.
## Effect size > 0 is the alternative setting for power.
effect_size_grid <- c(0, 1)

param_combinations <- expand.grid(
  diff_otu_direct = c("unbalanced", "balanced"),
  diff_otu_mode = c("rare", "mix", "abundant"),
  effect_size = effect_size_grid,
  depth_mu = c(10000, 1000, 100, 10),
  depth_theta = c(5, 10, 15),
  covariate_eff_sd = c(0, 0.5),
  stringsAsFactors = FALSE
)

total_param_sets <- nrow(param_combinations)

cat("Total parameter settings:", total_param_sets, "\n")
cat("Replications per setting:", n_simulation, "\n")
cat("Distance methods:", paste(distance_methods, collapse = ", "), "\n")
cat("ADONIS permutations:", ADONIS_PERM, "\n")
cat("Base seed:", base_seed, "\n")
cat("Use confounder:", USE_CONFOUNDER, "\n")
cat("Adjust for confounders:", ADJUST_FOR_CONFOUNDERS, "\n\n")


## ============================================================
## 10. Run one simulated dataset
## ============================================================

run_one_beta_dataset <- function(Simulated_data,
                                 rep_id,
                                 param_set,
                                 method_seed,
                                 distance_methods,
                                 tree,
                                 permutations = 999,
                                 alpha = 0.05,
                                 adjust_confounders = FALSE) {
  
  otu_tab <- Simulated_data$otu.tab.sim
  
  if (!is.matrix(otu_tab)) {
    otu_tab <- as.matrix(otu_tab)
  }
  
  ## SimulateMSeqU returns OTU x sample.
  ## Beta-diversity methods need sample x OTU.
  counts_mat <- t(otu_tab)
  
  storage.mode(counts_mat) <- "numeric"
  
  if (is.null(rownames(counts_mat))) {
    rownames(counts_mat) <- paste0("sample_", seq_len(nrow(counts_mat)))
  }
  
  ## Remove all-zero samples if they occur.
  libsize <- rowSums(counts_mat)
  keep_samp <- is.finite(libsize) & libsize > 0
  
  if (!all(keep_samp)) {
    counts_mat <- counts_mat[keep_samp, , drop = FALSE]
  }
  
  x_vec <- as.numeric(Simulated_data$covariate[, 1])
  x_vec <- x_vec[keep_samp]
  
  Z <- Simulated_data$confounder
  if (!is.null(Z)) {
    Z <- Z[keep_samp, , drop = FALSE]
  }
  
  if (length(unique(x_vec)) != 2L) {
    out <- tibble(
      rep_id = rep_id,
      param_set = param_set,
      method = distance_methods,
      p_value = NA_real_,
      reject = NA_integer_,
      n_samples = nrow(counts_mat),
      group0_n = sum(x_vec == 0),
      group1_n = sum(x_vec == 1)
    )
    return(out)
  }
  
  tree_use <- get_tree_for_otus(tree, colnames(counts_mat))
  counts_mat <- counts_mat[, tree_use$tip.label, drop = FALSE]
  
  physeq_full <- make_phyloseq_full(
    counts_mat = counts_mat,
    tree = tree_use,
    x_vec = x_vec,
    Z = Z
  )
  
  pvals <- run_permanova(
    physeq = physeq_full,
    distance_methods = distance_methods,
    permutations = permutations,
    seed = method_seed,
    adjust_confounders = adjust_confounders
  )
  
  tibble(
    rep_id = rep_id,
    param_set = param_set,
    method = names(pvals),
    p_value = as.numeric(pvals),
    reject = as.integer(!is.na(pvals) & pvals < alpha),
    n_samples = nrow(counts_mat),
    group0_n = sum(x_vec == 0),
    group1_n = sum(x_vec == 1)
  )
}


## ============================================================
## 11. Run simulation
## ============================================================

all_replicate_results_list <- vector("list", total_param_sets)

start_time <- Sys.time()

for (param_set in seq_len(total_param_sets)) {
  
  diff_otu_direct <- param_combinations$diff_otu_direct[param_set]
  diff_otu_mode <- param_combinations$diff_otu_mode[param_set]
  effect_size <- param_combinations$effect_size[param_set]
  depth_mu <- param_combinations$depth_mu[param_set]
  depth_theta <- param_combinations$depth_theta[param_set]
  covariate_eff_sd_original <- param_combinations$covariate_eff_sd[param_set]
  
  ## For a true null, both mean and SD of the covariate effect must be zero.
  ## Otherwise rnorm(mean = 0, sd = 0.5) would still create a non-null signal.
  covariate_eff_sd <- if (effect_size == 0) {
    0
  } else {
    covariate_eff_sd_original
  }
  
  cat(
    "Running parameter set", param_set, "of", total_param_sets,
    "| direct =", diff_otu_direct,
    "| mode =", diff_otu_mode,
    "| effect_size =", effect_size,
    "| depth_mu =", depth_mu,
    "| depth_theta =", depth_theta,
    "| cov_sd =", covariate_eff_sd,
    "\n"
  )
  
  replicate_results_this_setting <- vector("list", n_simulation)
  
  for (rep_id in seq_len(n_simulation)) {
    
    seed_sim <- base_seed + param_set * 100000L + rep_id
    seed_test <- base_seed + 50000000L + param_set * 100000L + rep_id
    
    set.seed(seed_sim)
    
    Simulated_data <- SimulateMSeqU(
      para = para1,
      nSam = nSam,
      nOTU = nOTU,
      diff.otu.pct = diff_otu_pct,
      diff.otu.direct = diff_otu_direct,
      diff.otu.mode = diff_otu_mode,
      user_specified_otu = NULL,
      covariate.type = "binary",
      grp.ratio = 1,
      covariate.eff.mean = effect_size,
      covariate.eff.sd = covariate_eff_sd,
      
      ## Primary beta-diversity run: no confounder.
      confounder.type = "none",
      conf.cov.cor = 0,
      conf.diff.otu.pct = 0,
      conf.nondiff.otu.pct = 0,
      confounder.eff.mean = 0,
      confounder.eff.sd = 0,
      
      error.sd = 0,
      depth.mu = depth_mu,
      depth.theta = depth_theta,
      
      ## No library-size confounding in primary beta-diversity test.
      depth.conf.factor = 0,
      cont.conf = 0,
      epsilon = NULL,
      
      ## Cluster/community signal.
      use.cluster.signal = TRUE,
      phylo.tree = tree_global,
      K_COMM = K_COMM,
      signal.cluster = NULL,
      cluster.assignment = pam_global,
      cluster.first_pole_method = "median",
      cluster.subsequent_pole_quantile = 1,
      cluster.assignment_method = "balanced",
      cluster.outlier_filter = TRUE,
      cluster.outlier_threshold = 1.5
    )
    
    beta_res <- run_one_beta_dataset(
      Simulated_data = Simulated_data,
      rep_id = rep_id,
      param_set = param_set,
      method_seed = seed_test,
      distance_methods = distance_methods,
      tree = tree_global,
      permutations = ADONIS_PERM,
      alpha = significance_level,
      adjust_confounders = ADJUST_FOR_CONFOUNDERS
    )
    
    beta_res <- beta_res %>%
      mutate(
        diff_otu_direct = diff_otu_direct,
        diff_otu_mode = diff_otu_mode,
        effect_size = effect_size,
        is_null = as.integer(effect_size == 0),
        depth_mu = depth_mu,
        depth_theta = depth_theta,
        covariate_eff_sd = covariate_eff_sd,
        covariate_eff_sd_original = covariate_eff_sd_original,
        K_COMM = K_COMM,
        signal_cluster = Simulated_data$signal.cluster,
        signal_cluster_ordered = Simulated_data$signal.cluster.ordered,
        signal_n_otus = sum(Simulated_data$diff.otu.ind),
        signal_otu_names = paste(Simulated_data$signal.otu.names, collapse = ";"),
        seed_sim = seed_sim,
        seed_test = seed_test
      )
    
    replicate_results_this_setting[[rep_id]] <- beta_res
    
    if (rep_id %% 10 == 0) {
      cat("  finished rep", rep_id, "of", n_simulation, "\n")
    }
  }
  
  replicate_results_this_setting <- bind_rows(replicate_results_this_setting)
  
  all_replicate_results_list[[param_set]] <- replicate_results_this_setting
  
  ## Save raw results for each parameter setting.
  raw_file <- file.path(
    outdir,
    "raw",
    sprintf("beta_div_msequ_param%03d.csv", param_set)
  )
  
  write.csv(replicate_results_this_setting, raw_file, row.names = FALSE)
}

all_replicate_results <- bind_rows(all_replicate_results_list)


## ============================================================
## 12. Summaries: power and type I error
## ============================================================

setting_summary <- all_replicate_results %>%
  group_by(
    method,
    diff_otu_direct,
    diff_otu_mode,
    effect_size,
    is_null,
    depth_mu,
    depth_theta,
    covariate_eff_sd_original,
    K_COMM
  ) %>%
  summarize(
    n_total = n(),
    n_valid = sum(!is.na(p_value)),
    valid_rate = n_valid / n_total,
    rejection_rate = mean(reject == 1, na.rm = TRUE),
    type1_error = ifelse(unique(is_null) == 1, rejection_rate, NA_real_),
    power = ifelse(unique(is_null) == 0, rejection_rate, NA_real_),
    mean_p = mean(p_value, na.rm = TRUE),
    median_p = median(p_value, na.rm = TRUE),
    min_p = min(p_value, na.rm = TRUE),
    max_p = max(p_value, na.rm = TRUE),
    mean_n_samples = mean(n_samples, na.rm = TRUE),
    mean_group0_n = mean(group0_n, na.rm = TRUE),
    mean_group1_n = mean(group1_n, na.rm = TRUE),
    mean_signal_n_otus = mean(signal_n_otus, na.rm = TRUE),
    .groups = "drop"
  )

method_summary <- all_replicate_results %>%
  group_by(method, effect_size, is_null) %>%
  summarize(
    n_total = n(),
    n_valid = sum(!is.na(p_value)),
    valid_rate = n_valid / n_total,
    rejection_rate = mean(reject == 1, na.rm = TRUE),
    type1_error = ifelse(unique(is_null) == 1, rejection_rate, NA_real_),
    power = ifelse(unique(is_null) == 0, rejection_rate, NA_real_),
    mean_p = mean(p_value, na.rm = TRUE),
    median_p = median(p_value, na.rm = TRUE),
    .groups = "drop"
  )

method_overall_summary <- all_replicate_results %>%
  group_by(method, is_null) %>%
  summarize(
    n_total = n(),
    n_valid = sum(!is.na(p_value)),
    valid_rate = n_valid / n_total,
    rejection_rate = mean(reject == 1, na.rm = TRUE),
    type1_error = ifelse(unique(is_null) == 1, rejection_rate, NA_real_),
    power = ifelse(unique(is_null) == 0, rejection_rate, NA_real_),
    mean_p = mean(p_value, na.rm = TRUE),
    median_p = median(p_value, na.rm = TRUE),
    .groups = "drop"
  )


## ============================================================
## 13. Save outputs
## ============================================================

replicate_csv <- file.path(outdir, "beta_div_msequ_replicate_results.csv")
setting_summary_csv <- file.path(outdir, "beta_div_msequ_setting_summary.csv")
method_summary_csv <- file.path(outdir, "beta_div_msequ_method_summary_by_effect.csv")
method_overall_summary_csv <- file.path(outdir, "beta_div_msequ_method_overall_summary.csv")
rds_file <- file.path(outdir, "beta_div_msequ_all_results.rds")

write.csv(all_replicate_results, replicate_csv, row.names = FALSE)
write.csv(setting_summary, setting_summary_csv, row.names = FALSE)
write.csv(method_summary, method_summary_csv, row.names = FALSE)
write.csv(method_overall_summary, method_overall_summary_csv, row.names = FALSE)

saveRDS(
  list(
    replicate_results = all_replicate_results,
    setting_summary = setting_summary,
    method_summary = method_summary,
    method_overall_summary = method_overall_summary,
    param_combinations = param_combinations,
    distance_methods = distance_methods,
    pam_global = pam_global,
    K_COMM = K_COMM,
    n_simulation = n_simulation,
    nSam = nSam,
    nOTU = nOTU,
    alpha = significance_level,
    ADONIS_PERM = ADONIS_PERM,
    base_seed = base_seed,
    use_confounder = USE_CONFOUNDER,
    adjust_for_confounders = ADJUST_FOR_CONFOUNDERS
  ),
  rds_file
)

end_time <- Sys.time()

cat("\n[DONE]\n")
cat("Start time:", as.character(start_time), "\n")
cat("End time:", as.character(end_time), "\n")
cat("Elapsed:", as.character(end_time - start_time), "\n\n")

cat("Saved files:\n")
cat("  ", replicate_csv, "\n")
cat("  ", setting_summary_csv, "\n")
cat("  ", method_summary_csv, "\n")
cat("  ", method_overall_summary_csv, "\n")
cat("  ", rds_file, "\n\n")

cat("Overall method summary:\n")
print(method_overall_summary, n = Inf)