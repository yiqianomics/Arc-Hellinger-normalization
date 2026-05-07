############################################################
## Corrected MSeqU simulation:
## Compare AHC / New Method vs ALR, CLR, ILR
##
## Main fixes:
## 1. Preserve sample IDs explicitly before merging metadata.
## 2. Check that merged data has nSam rows and two groups.
## 3. FDR = 0 when no discoveries.
## 4. Store raw and adjusted p-value diagnostics.
## 5. Keep current setting: 100 reps, same parameter grid.
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
  "compositions",
  "MASS",
  "dplyr",
  "tibble",
  "dirmult"
)

for (pkg in packages) {
  install_and_load(pkg, user_lib)
}


## ============================================================
## 2. Project path and para1 loading
## ============================================================

project_dir <- "/home/zhang.16383/AHC"
outdir <- file.path(project_dir, "out")

if (!dir.exists(project_dir)) {
  dir.create(project_dir, recursive = TRUE)
}

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

para_file <- file.path(project_dir, "para1.RData")

if (!file.exists(para_file)) {
  if (file.exists("para1.RData")) {
    para_file <- "para1.RData"
  } else {
    stop("Cannot find para1.RData. Put para1.RData under /home/zhang.16383/AHC or current working directory.")
  }
}

load(para_file)

if (!exists("para1")) {
  stop("para1.RData was loaded, but object 'para1' was not found.")
}


## ============================================================
## 3. Simulation functions
## Based on the pasted MSeqU-style simulation
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
  
  return(list(mu = ref.otu.tab.p, ref.otu.tab = ref.otu.tab0))
}


SimulateMSeqU <- function(
    para,
    nSam = 100,
    nOTU = 100,
    diff.otu.pct = 0.1,
    diff.otu.direct = c("balanced", "unbalanced"),
    diff.otu.mode = c("abundant", "rare", "mix", "user_specified"),
    user_specified_otu = NULL,
    covariate.type = c("binary", "continuous"),
    grp.ratio = 1,
    covariate.eff.mean = 1,
    covariate.eff.sd = 0,
    confounder.type = c("none", "binary", "continuous", "both"),
    conf.cov.cor = 0.6,
    conf.diff.otu.pct = 0,
    conf.nondiff.otu.pct = 0.1,
    confounder.eff.mean = 0,
    confounder.eff.sd = 0,
    error.sd = 0,
    depth.mu = 10000,
    depth.theta = 5,
    depth.conf.factor = 0,
    cont.conf = 0,
    epsilon = 0
) {
  
  diff.otu.direct <- match.arg(diff.otu.direct)
  diff.otu.mode <- match.arg(diff.otu.mode)
  covariate.type <- match.arg(covariate.type)
  confounder.type <- match.arg(confounder.type)
  
  model.paras <- para
  
  ## Select number of OTUs and samples
  ref.otu.tab <- model.paras$ref.otu.tab[seq_len(nOTU), , drop = FALSE]
  
  idx.otu <- rownames(ref.otu.tab)
  idx.sample <- colnames(model.paras$ref.otu.tab)[seq_len(nSam)]
  ref.otu.tab <- ref.otu.tab[, idx.sample, drop = FALSE]
  
  ## ----------------------------
  ## Confounder
  ## ----------------------------
  
  if (confounder.type == "none") {
    Z <- matrix(0, nrow = nSam, ncol = 1)
    confounder.eff.mean <- 0
    confounder.eff.sd <- 0
  }
  
  if (confounder.type == "continuous") {
    Z <- matrix(cont.conf, nrow = nSam, ncol = 1)
  }
  
  if (confounder.type == "binary") {
    Z <- matrix(
      c(rep(0, nSam %/% 2), rep(1, nSam - nSam %/% 2)),
      ncol = 1
    )
  }
  
  if (confounder.type == "both") {
    Z <- cbind(
      rnorm(nSam),
      c(rep(0, nSam %/% 2), rep(1, nSam - nSam %/% 2))
    )
  }
  
  rownames(Z) <- colnames(ref.otu.tab)
  
  ## ----------------------------
  ## Covariate of interest
  ## ----------------------------
  
  rho <- sqrt(conf.cov.cor^2 / (1 - conf.cov.cor^2))
  
  Z_summary <- scale(Z) %*% rep(1, ncol(Z))
  X_cont <- rho * scale(Z_summary) + epsilon
  
  if (covariate.type == "continuous") {
    X <- as.matrix(X_cont)
  }
  
  if (covariate.type == "binary") {
    cutoff <- quantile(X_cont, grp.ratio / (1 + grp.ratio))
    X <- cbind(ifelse(X_cont <= cutoff, 0, 1))
  }
  
  rownames(X) <- colnames(ref.otu.tab)
  
  covariate.eff.mean1 <- covariate.eff.mean
  covariate.eff.mean2 <- covariate.eff.mean
  
  ## ----------------------------
  ## OTU-level covariate effects
  ## ----------------------------
  
  if (diff.otu.direct == "balanced") {
    
    if (diff.otu.mode %in% c("abundant", "rare", "user_specified")) {
      
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
    
    if (diff.otu.mode %in% c("abundant", "rare", "user_specified")) {
      
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
  
  eta.diff <- effect_vec %*% t(scale(X))
  
  ## ----------------------------
  ## Confounder effects
  ## ----------------------------
  
  conf_effect_vec <- sample(c(
    rnorm(floor(nOTU / 2), mean = -confounder.eff.mean, sd = confounder.eff.sd),
    rnorm(nOTU - floor(nOTU / 2), mean = confounder.eff.mean, sd = confounder.eff.sd)
  ))
  
  eta.conf <- conf_effect_vec %*% t(scale(Z_summary))
  
  ## ----------------------------
  ## Choose true differential OTUs
  ## ----------------------------
  
  otu.ord <- seq_len(nOTU)
  diff.otu.num <- round(diff.otu.pct * nOTU)
  
  if (diff.otu.mode == "user_specified") {
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
  
  ## ----------------------------
  ## Choose confounded OTUs
  ## ----------------------------
  
  n_conf_diff <- round(nOTU * conf.diff.otu.pct)
  n_conf_nondiff <- round(nOTU * conf.nondiff.otu.pct)
  
  if (length(diff.otu.ind) >= n_conf_diff) {
    conf.otu.ind1 <- sample(diff.otu.ind, n_conf_diff)
  } else {
    conf.otu.ind1 <- diff.otu.ind
  }
  
  nondiff_pool <- setdiff(seq_len(nOTU), diff.otu.ind)
  
  conf.otu.ind <- c(
    conf.otu.ind1,
    sample(nondiff_pool, n_conf_nondiff)
  )
  
  ## ----------------------------
  ## Apply effects only to selected OTUs
  ## ----------------------------
  
  eta.diff[setdiff(seq_len(nOTU), diff.otu.ind), ] <- 0
  eta.conf[setdiff(seq_len(nOTU), conf.otu.ind), ] <- 0
  
  eta.error <- matrix(rnorm(nOTU * nSam, 0, error.sd), nrow = nOTU, ncol = nSam)
  
  eta.exp <- exp(t(eta.diff + eta.conf + eta.error))
  eta.exp <- eta.exp * t(ref.otu.tab)
  
  ref.otu.tab.prop <- eta.exp / rowSums(eta.exp)
  ref.otu.tab.prop <- t(ref.otu.tab.prop)
  
  ## ----------------------------
  ## Sequencing depth
  ## ----------------------------
  
  nSeq <- MASS::rnegbin(
    nSam,
    mu = as.numeric(depth.mu * exp(scale(X) * depth.conf.factor)),
    theta = depth.theta
  )
  
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
  
  return(list(
    otu.tab.sim = otu.tab.sim,
    covariate = X,
    confounder = Z,
    diff.otu.ind = diff.otu.ind.logical,
    otu.names = idx.otu,
    conf.otu.ind = conf.otu.ind.logical
  ))
}


## ============================================================
## 4. Correct data formatting
## This is the most important fix.
## ============================================================

make_analysis_data <- function(Simulated_data, nSam = 100, nOTU = 100) {
  
  otu_mat <- Simulated_data$otu.tab.sim
  
  if (!is.matrix(otu_mat)) {
    otu_mat <- as.matrix(otu_mat)
  }
  
  ## OTU table should be OTU x sample
  sample_names <- colnames(otu_mat)
  otu_names <- rownames(otu_mat)
  
  if (length(sample_names) != nSam) {
    stop("Sample number mismatch: colnames(otu.tab.sim) does not have length nSam.")
  }
  
  if (length(otu_names) != nOTU) {
    stop("OTU number mismatch: rownames(otu.tab.sim) does not have length nOTU.")
  }
  
  ## Explicitly preserve sample IDs.
  ## Do NOT rely on rownames after as.numeric().
  meta.dat <- data.frame(
    sample = sample_names,
    X = as.numeric(Simulated_data$covariate),
    Z1 = Simulated_data$confounder[, 1],
    Z2 = Simulated_data$confounder[, 2],
    stringsAsFactors = FALSE
  )
  
  ## Transpose OTU table to sample x OTU.
  t_otu_df <- as.data.frame(t(otu_mat), check.names = FALSE)
  t_otu_df$sample <- rownames(t_otu_df)
  
  ## Use left_join to preserve metadata/sample order.
  merged_data_current <- dplyr::left_join(
    meta.dat,
    t_otu_df,
    by = "sample"
  )
  
  if (nrow(merged_data_current) != nSam) {
    stop("Merged data does not have nSam rows. This indicates a sample-ID merge problem.")
  }
  
  if (anyNA(merged_data_current[, otu_names, drop = FALSE])) {
    stop("NA values found after merging OTU table with metadata. Check sample IDs.")
  }
  
  merged_data_current$Group <- ifelse(merged_data_current$X == 0, "A", "B")
  
  group_tab <- table(merged_data_current$Group)
  
  if (length(group_tab) != 2) {
    stop("Only one group found after creating Group. Check covariate generation.")
  }
  
  if (any(group_tab < 2)) {
    stop("At least one group has fewer than 2 samples. t-test is not valid.")
  }
  
  ## Keep only Group + OTU count columns.
  otu_df_processed <- merged_data_current[, c("Group", otu_names)]
  
  ## Construct truth vector in the same OTU order.
  diff_otu_ids <- Simulated_data$otu.names[Simulated_data$diff.otu.ind]
  diff_indices <- match(diff_otu_ids, otu_names)
  
  if (any(is.na(diff_indices))) {
    stop("Some differential OTU names could not be matched to the simulated OTU table.")
  }
  
  true_effects <- seq_len(nOTU) %in% diff_indices
  
  if (sum(true_effects) == 0) {
    stop("No true differential OTUs were found.")
  }
  
  ## Rename OTU columns after truth vector is constructed.
  colnames(otu_df_processed) <- c("Group", paste0("X", seq_len(nOTU)))
  
  return(list(
    data = otu_df_processed,
    true_effects = true_effects,
    group_table = group_tab
  ))
}


## ============================================================
## 5. Transformation helper functions
## ============================================================

replace_zeros_matrix <- function(X, constant = 0.5) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  X[X == 0] <- constant
  return(X)
}


normalize_rows <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  
  rs <- rowSums(X)
  
  if (any(rs <= 0)) {
    stop("Some rows have non-positive row sums.")
  }
  
  sweep(X, 1, rs, FUN = "/")
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
  
  ## If a sample has total count zero, map it to uniform composition.
  ## This prevents AHC from failing at very low sequencing depth.
  if (any(zero_rows)) {
    Pi[zero_rows, ] <- matrix(
      1 / p,
      nrow = sum(zero_rows),
      ncol = p
    )
  }
  
  rownames(Pi) <- rownames(X)
  colnames(Pi) <- colnames(X)
  
  return(Pi)
}


## ============================================================
## 6. Four transformations: AHC, ALR, CLR, ILR
## ============================================================

ahc_transformation <- function(data, group_factor = 1) {
  
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  
  ## AHC does not require ordinary zero replacement.
  Pi <- normalize_rows_allow_all_zero(X)
  
  k <- ncol(Pi)
  
  sqrt_Pi <- sqrt(Pi)
  sqrt_pi0 <- rep(1 / sqrt(k), k)
  
  c_pi <- as.vector(sqrt_Pi %*% sqrt_pi0)
  c_pi <- pmin(pmax(c_pi, 0), 1)
  
  s_pi <- sqrt(1 - c_pi^2)
  
  scale_factor <- asin(s_pi) / s_pi
  scale_factor[s_pi == 0] <- 1
  
  ahc_data <- (sqrt_Pi - outer(c_pi, sqrt_pi0)) * scale_factor
  
  colnames(ahc_data) <- colnames(Pi)
  
  out <- data.frame(Group = groups, ahc_data, check.names = FALSE)
  return(out)
}


alr_transformation <- function(data,
                               group_factor = 1,
                               ref_component,
                               zero_constant = 0.5) {
  
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  
  X <- replace_zeros_matrix(X, constant = zero_constant)
  Pi <- normalize_rows(X)
  
  p <- ncol(Pi)
  
  if (missing(ref_component) || is.null(ref_component)) {
    ref_component <- p
  }
  
  if (ref_component < 1 || ref_component > p) {
    stop("ref_component is out of range.")
  }
  
  alr_data <- log(Pi[, -ref_component, drop = FALSE] / Pi[, ref_component])
  
  colnames(alr_data) <- colnames(Pi)[-ref_component]
  
  out <- data.frame(Group = groups, alr_data, check.names = FALSE)
  return(out)
}


clr_transformation <- function(data,
                               group_factor = 1,
                               zero_constant = 0.5) {
  
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  
  X <- replace_zeros_matrix(X, constant = zero_constant)
  Pi <- normalize_rows(X)
  
  gm <- exp(rowMeans(log(Pi)))
  clr_data <- log(Pi / gm)
  
  colnames(clr_data) <- colnames(Pi)
  
  out <- data.frame(Group = groups, clr_data, check.names = FALSE)
  return(out)
}


ilr_transformation <- function(data,
                               group_factor = 1,
                               zero_constant = 0.5) {
  
  X <- data[, -group_factor, drop = FALSE]
  groups <- data[[group_factor]]
  
  X <- replace_zeros_matrix(X, constant = zero_constant)
  Pi <- normalize_rows(X)
  
  ilr_data <- compositions::ilr(compositions::acomp(Pi))
  ilr_data <- as.matrix(ilr_data)
  
  colnames(ilr_data) <- paste0("ILR", seq_len(ncol(ilr_data)))
  
  out <- data.frame(Group = groups, ilr_data, check.names = FALSE)
  return(out)
}


## ============================================================
## 7. Testing and metric functions
## ============================================================

run_feature_ttests <- function(transformed_data) {
  
  groups <- factor(transformed_data$Group)
  X <- transformed_data[, -1, drop = FALSE]
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  
  group_tab <- table(groups)
  
  if (length(group_tab) != 2) {
    stop("t-test requires exactly two groups.")
  }
  
  if (any(group_tab < 2)) {
    stop("Each group must have at least two samples.")
  }
  
  raw_p <- numeric(ncol(X))
  
  for (j in seq_len(ncol(X))) {
    
    y <- X[, j]
    
    raw_p[j] <- tryCatch({
      t.test(y ~ groups)$p.value
    }, error = function(e) {
      NA_real_
    })
  }
  
  adj_p <- p.adjust(raw_p, method = "bonferroni")
  
  out <- data.frame(
    feature = seq_along(raw_p),
    raw_p = raw_p,
    adj_p = adj_p
  )
  
  return(out)
}


choose_nondiff_reference <- function(true_effects) {
  
  nondiff_index <- which(!true_effects)
  
  if (length(nondiff_index) == 0) {
    stop("No non-differential feature is available as ALR reference.")
  }
  
  ## Current setting:
  ## choose the last known non-differential component as ALR reference.
  ref_component <- tail(nondiff_index, 1)
  
  return(ref_component)
}


safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  min(x, na.rm = TRUE)
}


calculate_power_fdr <- function(p_values,
                                true_effects,
                                significance_level = 0.05) {
  
  if (length(p_values) != length(true_effects)) {
    stop(
      "Length mismatch: length(p_values) = ",
      length(p_values),
      ", length(true_effects) = ",
      length(true_effects)
    )
  }
  
  rejected <- !is.na(p_values) & p_values < significance_level
  
  TP <- sum(rejected & true_effects)
  FP <- sum(rejected & !true_effects)
  FN <- sum(!rejected & true_effects)
  
  significant_count <- sum(true_effects)
  
  power <- if (significant_count == 0) {
    NA_real_
  } else {
    TP / significant_count
  }
  
  ## Fixed:
  ## If there are no discoveries, realized FDP/FDR is 0.
  fdr <- if ((TP + FP) == 0) {
    0
  } else {
    FP / (TP + FP)
  }
  
  return(list(
    Power = power,
    FDR = fdr,
    TP = TP,
    FP = FP,
    FN = FN,
    Rejections = sum(rejected)
  ))
}


make_method_result <- function(rep_id,
                               transformation,
                               p_table,
                               true_effects,
                               alr_ref,
                               significance_level = 0.05) {
  
  metrics <- calculate_power_fdr(
    p_values = p_table$adj_p,
    true_effects = true_effects,
    significance_level = significance_level
  )
  
  signal_raw <- p_table$raw_p[true_effects]
  null_raw <- p_table$raw_p[!true_effects]
  
  signal_adj <- p_table$adj_p[true_effects]
  null_adj <- p_table$adj_p[!true_effects]
  
  out <- tibble(
    rep_id = rep_id,
    transformation = transformation,
    power = metrics$Power,
    fdr = metrics$FDR,
    TP = metrics$TP,
    FP = metrics$FP,
    FN = metrics$FN,
    rejections = metrics$Rejections,
    n_tested = length(p_table$adj_p),
    n_na_p = sum(is.na(p_table$adj_p)),
    min_raw_p_signal = safe_min(signal_raw),
    min_raw_p_null = safe_min(null_raw),
    min_adj_p_signal = safe_min(signal_adj),
    min_adj_p_null = safe_min(null_adj),
    n_raw_p_signal_less_005 = sum(signal_raw < 0.05, na.rm = TRUE),
    n_adj_p_signal_less_005 = sum(signal_adj < 0.05, na.rm = TRUE),
    alr_ref = alr_ref
  )
  
  return(out)
}


run_one_dataset <- function(data,
                            true_effects,
                            rep_id,
                            significance_level = 0.05) {
  
  alr_ref <- choose_nondiff_reference(true_effects)
  
  method_results <- tibble()
  
  ## --------------------
  ## AHC
  ## --------------------
  
  transformed_AHC <- ahc_transformation(data)
  p_AHC <- run_feature_ttests(transformed_AHC)
  
  method_results <- bind_rows(
    method_results,
    make_method_result(
      rep_id = rep_id,
      transformation = "AHC",
      p_table = p_AHC,
      true_effects = true_effects,
      alr_ref = NA_integer_,
      significance_level = significance_level
    )
  )
  
  ## --------------------
  ## ALR
  ## --------------------
  
  transformed_ALR <- alr_transformation(data, ref_component = alr_ref)
  p_ALR <- run_feature_ttests(transformed_ALR)
  
  truth_ALR <- true_effects[-alr_ref]
  
  method_results <- bind_rows(
    method_results,
    make_method_result(
      rep_id = rep_id,
      transformation = "ALR",
      p_table = p_ALR,
      true_effects = truth_ALR,
      alr_ref = alr_ref,
      significance_level = significance_level
    )
  )
  
  ## --------------------
  ## CLR
  ## --------------------
  
  transformed_CLR <- clr_transformation(data)
  p_CLR <- run_feature_ttests(transformed_CLR)
  
  method_results <- bind_rows(
    method_results,
    make_method_result(
      rep_id = rep_id,
      transformation = "CLR",
      p_table = p_CLR,
      true_effects = true_effects,
      alr_ref = NA_integer_,
      significance_level = significance_level
    )
  )
  
  ## --------------------
  ## ILR
  ## --------------------
  
  transformed_ILR <- ilr_transformation(data)
  p_ILR <- run_feature_ttests(transformed_ILR)
  
  ## Important:
  ## ILR coordinates are not one-to-one OTU features.
  ## This keeps the same p-1 dimensional comparison structure as ALR.
  truth_ILR <- true_effects[-alr_ref]
  
  method_results <- bind_rows(
    method_results,
    make_method_result(
      rep_id = rep_id,
      transformation = "ILR",
      p_table = p_ILR,
      true_effects = truth_ILR,
      alr_ref = alr_ref,
      significance_level = significance_level
    )
  )
  
  return(method_results)
}


## ============================================================
## 8. Current simulation settings
## ============================================================

base_seed <- as.integer(Sys.getenv("BASE_SEED", "12345"))
set.seed(base_seed)

n_simulation <- 100
nSam <- 100
nOTU <- 100
significance_level <- 0.05

param_combinations <- expand.grid(
  diff_otu_direct = c("unbalanced", "balanced"),
  diff_otu_mode = c("rare", "mix", "abundant"),
  depth_mu = c(10000, 1000, 100, 10),
  depth_theta = c(5, 10, 15),
  covariate_eff_sd = c(0, 0.5),
  confounder_eff_sd = c(0, 0.5),
  depth_conf_factor = c(0, 0.5),
  stringsAsFactors = FALSE
)

total_param_sets <- nrow(param_combinations)

all_replicate_results_list <- list()
all_summary_results_list <- list()

start_time <- Sys.time()

cat("Total parameter settings:", total_param_sets, "\n")
cat("Replications per setting:", n_simulation, "\n")
cat("Base seed:", base_seed, "\n\n")


## ============================================================
## 9. Run simulation
## ============================================================

for (param_set in seq_len(total_param_sets)) {
  
  diff_otu_direct <- param_combinations$diff_otu_direct[param_set]
  diff_otu_mode <- param_combinations$diff_otu_mode[param_set]
  depth_mu <- param_combinations$depth_mu[param_set]
  depth_theta <- param_combinations$depth_theta[param_set]
  covariate_eff_sd <- param_combinations$covariate_eff_sd[param_set]
  confounder_eff_sd <- param_combinations$confounder_eff_sd[param_set]
  depth_conf_factor <- param_combinations$depth_conf_factor[param_set]
  
  cat(
    "Running parameter set", param_set, "of", total_param_sets,
    "| direct =", diff_otu_direct,
    "| mode =", diff_otu_mode,
    "| depth_mu =", depth_mu,
    "| depth_theta =", depth_theta,
    "| cov_sd =", covariate_eff_sd,
    "| conf_sd =", confounder_eff_sd,
    "| depth_conf =", depth_conf_factor,
    "\n"
  )
  
  replicate_results_this_setting <- tibble()
  
  for (rep_id in seq_len(n_simulation)) {
    
    ## Reproducible seed for each setting and replicate
    set.seed(base_seed + param_set * 100000 + rep_id)
    
    Simulated_data <- SimulateMSeqU(
      para = para1,
      nSam = nSam,
      nOTU = nOTU,
      diff.otu.pct = 0.1,
      diff.otu.direct = diff_otu_direct,
      diff.otu.mode = diff_otu_mode,
      user_specified_otu = NULL,
      covariate.type = "binary",
      grp.ratio = 1,
      covariate.eff.mean = 1,
      covariate.eff.sd = covariate_eff_sd,
      confounder.type = "both",
      conf.cov.cor = 0.6,
      conf.diff.otu.pct = 0,
      conf.nondiff.otu.pct = 0.1,
      confounder.eff.mean = 0,
      confounder.eff.sd = confounder_eff_sd,
      error.sd = 0,
      depth.mu = depth_mu,
      depth.theta = depth_theta,
      depth.conf.factor = depth_conf_factor,
      cont.conf = 0,
      epsilon = 0
    )
    
    analysis_obj <- make_analysis_data(
      Simulated_data = Simulated_data,
      nSam = nSam,
      nOTU = nOTU
    )
    
    one_result <- run_one_dataset(
      data = analysis_obj$data,
      true_effects = analysis_obj$true_effects,
      rep_id = rep_id,
      significance_level = significance_level
    )
    
    one_result <- one_result %>%
      mutate(
        diff_otu_direct = diff_otu_direct,
        diff_otu_mode = diff_otu_mode,
        depth_mu = depth_mu,
        depth_theta = depth_theta,
        covariate_eff_sd = covariate_eff_sd,
        confounder_eff_sd = confounder_eff_sd,
        depth_conf_factor = depth_conf_factor,
        n_group_A = as.integer(analysis_obj$group_table["A"]),
        n_group_B = as.integer(analysis_obj$group_table["B"])
      )
    
    replicate_results_this_setting <- bind_rows(
      replicate_results_this_setting,
      one_result
    )
  }
  
  summary_this_setting <- replicate_results_this_setting %>%
    group_by(
      diff_otu_direct,
      diff_otu_mode,
      depth_mu,
      depth_theta,
      covariate_eff_sd,
      confounder_eff_sd,
      depth_conf_factor,
      transformation
    ) %>%
    summarize(
      n_reps = n(),
      mean_power = mean(power, na.rm = TRUE),
      sd_power = sd(power, na.rm = TRUE),
      mean_fdr = mean(fdr, na.rm = TRUE),
      sd_fdr = sd(fdr, na.rm = TRUE),
      mean_TP = mean(TP, na.rm = TRUE),
      mean_FP = mean(FP, na.rm = TRUE),
      mean_FN = mean(FN, na.rm = TRUE),
      mean_rejections = mean(rejections, na.rm = TRUE),
      mean_n_na_p = mean(n_na_p, na.rm = TRUE),
      mean_min_raw_p_signal = mean(min_raw_p_signal, na.rm = TRUE),
      mean_min_adj_p_signal = mean(min_adj_p_signal, na.rm = TRUE),
      mean_n_raw_p_signal_less_005 = mean(n_raw_p_signal_less_005, na.rm = TRUE),
      mean_n_adj_p_signal_less_005 = mean(n_adj_p_signal_less_005, na.rm = TRUE),
      .groups = "drop"
    )
  
  all_replicate_results_list[[param_set]] <- replicate_results_this_setting
  all_summary_results_list[[param_set]] <- summary_this_setting
  
  cat("Completed parameter set", param_set, "of", total_param_sets, "\n\n")
}


## ============================================================
## 10. Combine and save results
## ============================================================

all_replicate_results <- bind_rows(all_replicate_results_list)
final_results <- bind_rows(all_summary_results_list)

replicate_file <- file.path(
  outdir,
  "replicate_results_AHC_ALR_CLR_ILR_MSeqU_100rep_corrected.csv"
)

summary_file <- file.path(
  outdir,
  "summary_results_AHC_ALR_CLR_ILR_MSeqU_100rep_corrected.csv"
)

write.csv(
  all_replicate_results,
  replicate_file,
  row.names = FALSE
)

write.csv(
  final_results,
  summary_file,
  row.names = FALSE
)

end_time <- Sys.time()

cat("\nSimulation finished.\n")
cat("Start time:", as.character(start_time), "\n")
cat("End time:", as.character(end_time), "\n")
cat("Total runtime:", as.character(end_time - start_time), "\n\n")

cat("Replicate-level results saved to:\n")
cat(replicate_file, "\n\n")

cat("Summary results saved to:\n")
cat(summary_file, "\n\n")


## ============================================================
## 11. Quick sanity checks printed at the end
## ============================================================

cat("Quick check: overall method summary\n")

overall_check <- all_replicate_results %>%
  group_by(transformation) %>%
  summarize(
    mean_power = mean(power, na.rm = TRUE),
    max_power = max(power, na.rm = TRUE),
    mean_fdr = mean(fdr, na.rm = TRUE),
    mean_rejections = mean(rejections, na.rm = TRUE),
    max_rejections = max(rejections, na.rm = TRUE),
    mean_n_na_p = mean(n_na_p, na.rm = TRUE),
    mean_raw_signal_hits = mean(n_raw_p_signal_less_005, na.rm = TRUE),
    mean_adj_signal_hits = mean(n_adj_p_signal_less_005, na.rm = TRUE),
    .groups = "drop"
  )

print(overall_check)

cat("\nQuick check: settings with largest power\n")

top_power_settings <- final_results %>%
  arrange(desc(mean_power), desc(mean_rejections)) %>%
  head(20)

print(top_power_settings)