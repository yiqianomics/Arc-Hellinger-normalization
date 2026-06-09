############################################################
## Real gut data analysis for Arc-Hellinger normalization
##
## Outcome:
##   arm
##
## Confounders:
##   age, day
##
## Notes:
##   day < 0 is before surgery, day = 0 is surgery day,
##   and day > 0 is after surgery. The adjusted models use
##   day as a numeric covariate.
##
## Outputs:
##   Analysis/RealDataGut/real_data_gut_results/tables
##   Analysis/RealDataGut/real_data_gut_results/figures
##
## The AHC implementation is sourced from R/ArcH.R and is not
## modified here.
############################################################

options(stringsAsFactors = FALSE, width = 140)
set.seed(20260608)

get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE))
  }
  NA_character_
}

script_path <- get_script_path()
analysis_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()
project_dir <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = TRUE)

out_dir <- file.path(analysis_dir, "real_data_gut_results")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
if (dir.exists(out_dir)) {
  unlink(file.path(out_dir, c("figures", "tables", "run_summary.txt")), recursive = TRUE, force = TRUE)
}
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

core_packages <- c(
  "phyloseq", "vegan", "ggplot2", "dplyr", "tidyr",
  "tibble", "purrr", "patchwork", "lefser", "SummarizedExperiment",
  "limma", "UpSetR", "cowplot"
)
missing_core <- core_packages[!vapply(core_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_core) > 0) {
  stop("Missing required packages: ", paste(missing_core, collapse = ", "))
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(patchwork)
})

theme_set(
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold", color = "grey15"),
      plot.title.position = "plot",
      legend.position = "bottom"
    )
)

save_plot <- function(plot, filename, width, height) {
  png_file <- file.path(fig_dir, paste0(filename, ".png"))
  pdf_file <- file.path(fig_dir, paste0(filename, ".pdf"))
  ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = 320)
  ggplot2::ggsave(pdf_file, plot, width = width, height = height)
  invisible(c(png_file, pdf_file))
}

clean_method_label <- function(x) {
  dplyr::recode(
    x,
    "AHC_D_norm" = "AHC (dominance)",
    "AHC_E" = "AHC (evenness)",
    "Shannon_drop0_evenness" = "Shannon (drop 0s)",
    "Shannon_pc0.01_evenness" = "Shannon (ps = 0.01)",
    "Shannon_pc0.5_evenness" = "Shannon (ps = 0.5)",
    "Shannon_pc1_evenness" = "Shannon (ps = 1)",
    "AHC_beta" = "AHC beta diversity",
    "AHC" = "AHC",
    "Aitchison_pc0.01" = "Aitchison (pc = 0.01)",
    "Aitchison_pc0.5" = "Aitchison (pc = 0.5)",
    "Aitchison_pc1" = "Aitchison (pc = 1)",
    "CLR_pc0.01" = "CLR (pc = 0.01)",
    "CLR_pc0.5" = "CLR (pc = 0.5)",
    "CLR_pc1" = "CLR (pc = 1)",
    "DESeq2" = "DESeq2",
    "DESeq2_LRT" = "DESeq2",
    "LEfSe" = "LEfSe",
    .default = x
  )
}

format_contrast <- function(g1, g2) paste0(g1, "_vs_", g2)

format_p <- function(p) {
  ifelse(
    is.na(p), "p = NA",
    ifelse(p < 1e-4, "p < 1e-4", paste0("p = ", formatC(p, format = "fg", digits = 2)))
  )
}

arm_levels <- c("control", "treatment")
arm_labels <- c(control = "Control", treatment = "Treatment")
arm_palette <- c(control = "#0072B2", treatment = "#D55E00")
method_palette <- c(
  AHC = "#AB2428",
  CLR_pc0.01 = "#2F74B8",
  CLR_pc0.5 = "#79B0D7",
  CLR_pc1 = "#A0C8E8",
  DESeq2 = "#009E73",
  LEfSe = "#E69F00"
)

pseudo_counts <- c(0.01, 0.5, 1)
names(pseudo_counts) <- c("0.01", "0.5", "1")

alpha_level <- 0.05
da_min_prevalence_fraction <- 0.10
da_p_adjust_method <- "fdr"
lefse_p_adjust_method <- "fdr"
da_min_total_count <- 50
lefse_lda_threshold <- 2

############################################################
## Load data and AHC transform
############################################################

ahc_file <- file.path(project_dir, "R", "ArcH.R")
if (!file.exists(ahc_file)) {
  stop("Cannot find AHC source file: ", ahc_file)
}
source(ahc_file)
if (!exists("angular_hellinger_normalization")) {
  stop("R/ArcH.R did not define angular_hellinger_normalization().")
}

gut_file <- file.path(analysis_dir, "gut.RData")
if (!file.exists(gut_file)) {
  stop("Cannot find gut.RData in: ", analysis_dir)
}
load(gut_file)
if (!exists("gut") || !inherits(gut, "phyloseq")) {
  stop("gut.RData must contain a phyloseq object named gut.")
}

ps_gut <- gut
ps_gut <- phyloseq::prune_samples(phyloseq::sample_sums(ps_gut) > 0, ps_gut)
ps_gut <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_gut) > 0, ps_gut)

meta_raw <- as(phyloseq::sample_data(ps_gut), "data.frame")
meta_raw$SampleID <- rownames(meta_raw)
required_columns <- c("arm", "age", "day")
missing_columns <- setdiff(required_columns, names(meta_raw))
if (length(missing_columns) > 0) {
  stop("Missing required sample_data columns: ", paste(missing_columns, collapse = ", "))
}

meta_raw <- meta_raw %>%
  dplyr::mutate(
    arm = factor(as.character(arm), levels = arm_levels),
    age = suppressWarnings(as.numeric(age)),
    day = suppressWarnings(as.numeric(day)),
    DayPhase = factor(
      dplyr::case_when(
        day < 0 ~ "Before surgery",
        day == 0 ~ "Surgery day",
        day > 0 ~ "After surgery",
        TRUE ~ NA_character_
      ),
      levels = c("Before surgery", "Surgery day", "After surgery")
    )
  )

complete_sample <- stats::complete.cases(meta_raw[, c("arm", "age", "day")])
excluded_samples <- meta_raw %>%
  dplyr::filter(!complete_sample) %>%
  dplyr::select(SampleID, arm, age, day)
write.csv(excluded_samples, file.path(tab_dir, "samples_excluded_missing_arm_age_day.csv"), row.names = FALSE)

ps_gut <- phyloseq::prune_samples(meta_raw$SampleID[complete_sample], ps_gut)
ps_gut <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_gut) > 0, ps_gut)

meta <- meta_raw[match(phyloseq::sample_names(ps_gut), meta_raw$SampleID), , drop = FALSE]
rownames(meta) <- meta$SampleID
meta$arm <- droplevels(meta$arm)
if (length(levels(meta$arm)) < 2) {
  stop("At least two arm levels are required after complete-case filtering.")
}
if (!all(levels(meta$arm) %in% names(arm_labels))) {
  arm_labels <- stats::setNames(tools::toTitleCase(levels(meta$arm)), levels(meta$arm))
  arm_palette <- stats::setNames(grDevices::hcl.colors(length(levels(meta$arm)), "Dark 3"), levels(meta$arm))
}

arm_n <- table(meta$arm)
arm_axis_labels <- setNames(
  paste0(arm_labels[levels(meta$arm)], "\n(n = ", as.integer(arm_n[levels(meta$arm)]), ")"),
  levels(meta$arm)
)
meta$ArmLabel <- factor(arm_labels[as.character(meta$arm)], levels = unname(arm_labels[levels(meta$arm)]))
meta$ArmLabelN <- factor(arm_axis_labels[as.character(meta$arm)], levels = arm_axis_labels[levels(meta$arm)])

counts_all <- as(phyloseq::otu_table(ps_gut), "matrix")
if (phyloseq::taxa_are_rows(ps_gut)) {
  counts_all <- t(counts_all)
}
storage.mode(counts_all) <- "numeric"
counts_all <- counts_all[meta$SampleID, , drop = FALSE]
counts_all <- counts_all[, colSums(counts_all) > 0, drop = FALSE]

tax_df <- as(phyloseq::tax_table(ps_gut), "matrix")
tax_df <- as.data.frame(tax_df, stringsAsFactors = FALSE, check.names = FALSE)
tax_df$Taxon <- rownames(tax_df)
tax_df <- tax_df %>% dplyr::relocate(Taxon)

lib_size <- rowSums(counts_all)
rel_all <- sweep(counts_all, 1, lib_size, FUN = "/")

data_summary <- tibble::tibble(
  metric = c(
    "samples_original",
    "samples_nonzero",
    "samples_complete_case_arm_age_day",
    "samples_excluded_missing_arm_age_day",
    "taxa_original",
    "taxa_used_nonzero_in_complete_cases",
    "zero_fraction_selected",
    "median_library_size",
    "min_library_size",
    "max_library_size",
    paste0(levels(meta$arm), "_n"),
    "before_surgery_n",
    "surgery_day_n",
    "after_surgery_n"
  ),
  value = c(
    phyloseq::nsamples(gut),
    nrow(meta_raw),
    nrow(counts_all),
    nrow(excluded_samples),
    phyloseq::ntaxa(gut),
    ncol(counts_all),
    mean(counts_all == 0),
    median(lib_size),
    min(lib_size),
    max(lib_size),
    as.integer(arm_n[levels(meta$arm)]),
    sum(meta$day < 0),
    sum(meta$day == 0),
    sum(meta$day > 0)
  )
)
write.csv(data_summary, file.path(tab_dir, "data_summary.csv"), row.names = FALSE)

covariate_summary <- meta %>%
  dplyr::group_by(arm, ArmLabel) %>%
  dplyr::summarise(
    n = dplyr::n(),
    age_mean = mean(age),
    age_sd = sd(age),
    age_median = median(age),
    day_mean = mean(day),
    day_sd = sd(day),
    day_median = median(day),
    before_surgery_n = sum(day < 0),
    surgery_day_n = sum(day == 0),
    after_surgery_n = sum(day > 0),
    .groups = "drop"
  )
write.csv(covariate_summary, file.path(tab_dir, "sample_covariate_summary_by_arm.csv"), row.names = FALSE)

############################################################
## AHC checks
############################################################

ahc_all <- angular_hellinger_normalization(counts_all)
p_all <- ncol(counts_all)
ahc_norm <- sqrt(rowSums(ahc_all^2))
ahc_bound <- acos(1 / sqrt(p_all))

c_pi <- as.vector(sqrt(rel_all) %*% rep(1 / sqrt(p_all), p_all))
c_pi <- pmin(pmax(c_pi, 0), 1)
theta <- acos(c_pi)

ahc_check <- tibble::tibble(
  check = c(
    "max_abs_coordinate_row_sum",
    "max_abs_norm_minus_angle",
    "min_norm",
    "max_norm",
    "theoretical_bound",
    "n_bound_violations"
  ),
  value = c(
    max(abs(rowSums(ahc_all))),
    max(abs(ahc_norm - theta)),
    min(ahc_norm),
    max(ahc_norm),
    ahc_bound,
    sum(ahc_norm > ahc_bound + 1e-10)
  )
)
write.csv(ahc_check, file.path(tab_dir, "ahc_function_checks.csv"), row.names = FALSE)

############################################################
## 1. Alpha diversity
############################################################

shannon_from_counts <- function(x) {
  vegan::diversity(x, index = "shannon")
}

shannon_with_pseudocount <- function(x, pc) {
  y <- x + pc
  p <- sweep(y, 1, rowSums(y), FUN = "/")
  -rowSums(p * log(p))
}

alpha_wide <- tibble::tibble(
  SampleID = meta$SampleID,
  arm = meta$arm,
  ArmLabel = meta$ArmLabel,
  ArmLabelN = meta$ArmLabelN,
  age = meta$age,
  day = meta$day,
  DayPhase = meta$DayPhase,
  LibrarySize = lib_size,
  Richness = vegan::specnumber(counts_all),
  AHC_D = ahc_norm,
  AHC_D_norm = ahc_norm / ahc_bound,
  AHC_E = 1 - ahc_norm / ahc_bound,
  Shannon_drop0 = shannon_from_counts(counts_all)
)
alpha_wide$Shannon_drop0_evenness <- alpha_wide$Shannon_drop0 / log(p_all)

for (pc_name in names(pseudo_counts)) {
  pc <- pseudo_counts[[pc_name]]
  shannon_name <- paste0("Shannon_pc", pc_name)
  evenness_name <- paste0(shannon_name, "_evenness")
  alpha_wide[[shannon_name]] <- shannon_with_pseudocount(counts_all, pc)
  alpha_wide[[evenness_name]] <- alpha_wide[[shannon_name]] / log(p_all)
}

write.csv(alpha_wide, file.path(tab_dir, "alpha_diversity_by_sample.csv"), row.names = FALSE)

alpha_norm_cols <- c(
  "AHC_D_norm",
  "AHC_E",
  "Shannon_drop0_evenness",
  "Shannon_pc0.01_evenness",
  "Shannon_pc0.5_evenness",
  "Shannon_pc1_evenness"
)

alpha_long <- alpha_wide %>%
  dplyr::select(SampleID, arm, ArmLabel, ArmLabelN, age, day, DayPhase, dplyr::all_of(alpha_norm_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(alpha_norm_cols),
    names_to = "method",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    MethodLabel = factor(clean_method_label(method), levels = clean_method_label(alpha_norm_cols))
  )

alpha_summary <- alpha_long %>%
  dplyr::group_by(method, MethodLabel, arm, ArmLabel) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(alpha_summary, file.path(tab_dir, "alpha_arm_summary.csv"), row.names = FALSE)

alpha_adjusted_lm_df <- function(dat) {
  fit <- stats::lm(value ~ arm + age + day, data = dat)
  drop_tab <- stats::drop1(fit, test = "F")
  arm_row <- drop_tab["arm", , drop = FALSE]
  coef_tab <- stats::coef(summary(fit))
  arm_coef <- grep("^arm", rownames(coef_tab), value = TRUE)
  tibble::tibble(
    df_arm = arm_row$Df,
    statistic = arm_row$F,
    p_value = arm_row$`Pr(>F)`,
    coefficient_name = if (length(arm_coef) == 1) arm_coef else NA_character_,
    coefficient_estimate = if (length(arm_coef) == 1) coef_tab[arm_coef, "Estimate"] else NA_real_,
    coefficient_se = if (length(arm_coef) == 1) coef_tab[arm_coef, "Std. Error"] else NA_real_,
    coefficient_t = if (length(arm_coef) == 1) coef_tab[arm_coef, "t value"] else NA_real_,
    coefficient_p_value = if (length(arm_coef) == 1) coef_tab[arm_coef, "Pr(>|t|)"] else NA_real_
  )
}

alpha_lm_adjusted <- alpha_long %>%
  dplyr::group_by(method, MethodLabel) %>%
  dplyr::group_modify(~ alpha_adjusted_lm_df(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(q_value = p.adjust(p_value, method = "BH"))
write.csv(alpha_lm_adjusted, file.path(tab_dir, "alpha_arm_lm_adjusted_tests.csv"), row.names = FALSE)

pairwise_wilcox_df <- function(dat) {
  pairs <- utils::combn(levels(dat$arm), 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pair) {
    sub <- dat %>% dplyr::filter(arm %in% pair)
    wt <- wilcox.test(value ~ droplevels(arm), data = sub, exact = FALSE)
    tibble::tibble(
      group1 = pair[1],
      group2 = pair[2],
      contrast = format_contrast(pair[1], pair[2]),
      statistic = unname(wt$statistic),
      p_value = wt$p.value
    )
  }) %>%
    dplyr::mutate(q_value = p.adjust(p_value, method = "BH"))
}

alpha_pairwise <- alpha_long %>%
  dplyr::group_by(method, MethodLabel) %>%
  dplyr::group_modify(~ pairwise_wilcox_df(.x)) %>%
  dplyr::ungroup()
write.csv(alpha_pairwise, file.path(tab_dir, "alpha_pairwise_wilcox_tests_unadjusted.csv"), row.names = FALSE)

alpha_lm_lab <- alpha_lm_adjusted %>%
  dplyr::mutate(label = paste0("Adjusted F = ", formatC(statistic, format = "fg", digits = 3), "\n", format_p(p_value)))

p_alpha <- ggplot(alpha_long, aes(x = ArmLabelN, y = value, color = arm, fill = arm)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.18, linewidth = 0.35) +
  geom_jitter(width = 0.13, height = 0, alpha = 0.42, size = 1.15) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 2.2,
    fill = "white",
    color = "grey15",
    stroke = 0.35
  ) +
  geom_label(
    data = alpha_lm_lab,
    aes(label = label),
    x = Inf,
    y = Inf,
    hjust = 1.03,
    vjust = 1.08,
    inherit.aes = FALSE,
    size = 2.65,
    linewidth = 0.15,
    label.padding = grid::unit(0.12, "lines"),
    fill = "white",
    alpha = 0.88,
    color = "grey15"
  ) +
  facet_wrap(~ MethodLabel, ncol = 3) +
  scale_color_manual(values = arm_palette[levels(meta$arm)], guide = "none") +
  scale_fill_manual(values = arm_palette[levels(meta$arm)], guide = "none") +
  labs(x = NULL, y = "Normalized alpha index") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text.x = element_text(size = 9.5),
    panel.spacing = grid::unit(0.85, "lines")
  )
save_plot(p_alpha, "fig1_alpha_diversity_by_arm", width = 10.5, height = 6.7)

make_scalar_config <- function(v) {
  d <- stats::dist(scale(v))
  pts <- stats::cmdscale(d, k = 2, eig = FALSE)
  pts <- as.matrix(pts)
  if (ncol(pts) < 2) {
    pts <- cbind(pts[, 1], 0)
  }
  colnames(pts) <- c("Axis1", "Axis2")
  pts
}

pairwise_procrustes <- function(configs, permutations = 999) {
  methods <- names(configs)
  purrr::map_dfr(methods, function(m1) {
    purrr::map_dfr(methods, function(m2) {
      if (m1 == m2) {
        return(tibble::tibble(method1 = m1, method2 = m2, procrustes_r = 1, p_value = NA_real_))
      }
      out <- tryCatch(
        vegan::protest(configs[[m1]], configs[[m2]], permutations = permutations),
        error = function(e) NULL
      )
      tibble::tibble(
        method1 = m1,
        method2 = m2,
        procrustes_r = if (is.null(out)) NA_real_ else unname(out$t0),
        p_value = if (is.null(out)) NA_real_ else out$signif
      )
    })
  })
}

alpha_configs <- lapply(alpha_norm_cols, function(nm) make_scalar_config(alpha_wide[[nm]]))
names(alpha_configs) <- alpha_norm_cols
alpha_procrustes <- pairwise_procrustes(alpha_configs) %>%
  dplyr::mutate(
    Method1 = factor(clean_method_label(method1), levels = clean_method_label(alpha_norm_cols)),
    Method2 = factor(clean_method_label(method2), levels = rev(clean_method_label(alpha_norm_cols))),
    method1_index = match(method1, alpha_norm_cols),
    method2_index = match(method2, alpha_norm_cols),
    label = sprintf("%.2f", procrustes_r)
  )
write.csv(alpha_procrustes, file.path(tab_dir, "alpha_method_procrustes.csv"), row.names = FALSE)

p_alpha_procrustes_data <- alpha_procrustes %>%
  dplyr::filter(method1_index >= method2_index)

p_alpha_procrustes <- ggplot(p_alpha_procrustes_data, aes(x = Method1, y = Method2, fill = procrustes_r)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 3.15, color = "grey10") +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC", limits = c(0, 1), na.value = "grey96") +
  labs(
    x = NULL,
    y = NULL,
    fill = "Procrustes r",
    caption = "Triangle shown; permutation p-values are saved in alpha_method_procrustes.csv."
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    plot.caption = element_text(color = "grey35", hjust = 0)
  )
save_plot(p_alpha_procrustes, "fig1b_alpha_method_procrustes", width = 7.4, height = 5.9)

############################################################
## 2. Beta diversity
############################################################

clr_transform <- function(x, pc) {
  y <- x + pc
  p <- sweep(y, 1, rowSums(y), FUN = "/")
  log_p <- log(p)
  sweep(log_p, 1, rowMeans(log_p), FUN = "-")
}

beta_dists <- list(AHC_beta = stats::dist(ahc_all, method = "euclidean"))
for (pc_name in names(pseudo_counts)) {
  beta_dists[[paste0("Aitchison_pc", pc_name)]] <- stats::dist(clr_transform(counts_all, pseudo_counts[[pc_name]]))
}

run_adonis_adjusted <- function(d, metadata, permutations = 999) {
  tab <- vegan::adonis2(d ~ arm + age + day, data = metadata, permutations = permutations, by = "margin")
  out <- as.data.frame(tab)
  out$term <- rownames(out)
  out <- out %>%
    dplyr::filter(term %in% c("arm", "age", "day")) %>%
    dplyr::transmute(
      term,
      df = Df,
      sumsq = SumOfSqs,
      r2 = R2,
      statistic = F,
      p_value = `Pr(>F)`
    )
  tibble::as_tibble(out)
}

beta_permanova <- purrr::imap_dfr(beta_dists, function(d, method_name) {
  run_adonis_adjusted(d, meta) %>%
    dplyr::mutate(method = method_name, MethodLabel = clean_method_label(method_name), .before = 1)
}) %>%
  dplyr::group_by(term) %>%
  dplyr::mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup()
write.csv(beta_permanova, file.path(tab_dir, "beta_permanova_adjusted.csv"), row.names = FALSE)

pcoa_one <- function(d, method_name) {
  fit <- cmdscale(d, k = 2, eig = TRUE)
  eig <- fit$eig
  pos_sum <- sum(eig[eig > 0], na.rm = TRUE)
  var1 <- if (pos_sum > 0) 100 * eig[1] / pos_sum else NA_real_
  var2 <- if (pos_sum > 0) 100 * eig[2] / pos_sum else NA_real_
  tibble::tibble(
    SampleID = rownames(fit$points),
    Axis1 = fit$points[, 1],
    Axis2 = fit$points[, 2],
    method = method_name,
    MethodLabel = clean_method_label(method_name),
    AxisLabel = sprintf("PCoA1 %.1f%%\nPCoA2 %.1f%%", var1, var2)
  )
}

pcoa_df <- purrr::imap_dfr(beta_dists, pcoa_one) %>%
  dplyr::left_join(meta %>% dplyr::select(SampleID, arm, ArmLabel, age, day, DayPhase), by = "SampleID") %>%
  dplyr::mutate(MethodLabel = factor(MethodLabel, levels = clean_method_label(names(beta_dists))))
write.csv(pcoa_df, file.path(tab_dir, "beta_pcoa_coordinates.csv"), row.names = FALSE)

pcoa_lab <- pcoa_df %>%
  dplyr::group_by(MethodLabel, AxisLabel) %>%
  dplyr::summarise(
    x = min(Axis1, na.rm = TRUE),
    y = max(Axis2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    beta_permanova %>%
      dplyr::filter(term == "arm") %>%
      dplyr::select(MethodLabel, r2, statistic, p_value),
    by = "MethodLabel"
  ) %>%
  dplyr::mutate(
    label = paste0(
      AxisLabel,
      "\nAdjusted PERMANOVA R2 = ", sprintf("%.3f", r2),
      "\n", format_p(p_value)
    )
  )

p_beta_pcoa <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2)) +
  stat_ellipse(aes(color = arm, group = arm), linewidth = 0.45, alpha = 0.78, show.legend = FALSE) +
  geom_point(aes(fill = arm), shape = 21, size = 2.0, alpha = 0.86, color = "white", stroke = 0.22) +
  geom_label(
    data = pcoa_lab,
    aes(label = label),
    x = -Inf,
    y = Inf,
    inherit.aes = FALSE,
    hjust = -0.03,
    vjust = 1.08,
    size = 2.55,
    lineheight = 0.92,
    linewidth = 0.15,
    label.padding = grid::unit(0.13, "lines"),
    fill = "white",
    alpha = 0.88,
    color = "grey20"
  ) +
  facet_wrap(~ MethodLabel, scales = "free", nrow = 1) +
  scale_color_manual(values = arm_palette[levels(meta$arm)], labels = arm_labels[levels(meta$arm)]) +
  scale_fill_manual(values = arm_palette[levels(meta$arm)], labels = arm_labels[levels(meta$arm)]) +
  labs(x = "PCoA coordinate 1", y = "PCoA coordinate 2", fill = "Arm") +
  theme(panel.spacing = grid::unit(0.9, "lines"))
save_plot(p_beta_pcoa, "fig2_beta_pcoa_by_method", width = 12.8, height = 4.6)

beta_configs <- lapply(names(beta_dists), function(nm) {
  pcoa_df %>%
    dplyr::filter(method == nm) %>%
    dplyr::arrange(match(SampleID, meta$SampleID)) %>%
    dplyr::select(Axis1, Axis2) %>%
    as.matrix()
})
names(beta_configs) <- names(beta_dists)

beta_procrustes <- pairwise_procrustes(beta_configs) %>%
  dplyr::mutate(
    Method1 = factor(clean_method_label(method1), levels = clean_method_label(names(beta_dists))),
    Method2 = factor(clean_method_label(method2), levels = rev(clean_method_label(names(beta_dists)))),
    method1_index = match(method1, names(beta_dists)),
    method2_index = match(method2, names(beta_dists)),
    label = sprintf("%.2f", procrustes_r)
  )
write.csv(beta_procrustes, file.path(tab_dir, "beta_pcoa_procrustes.csv"), row.names = FALSE)

p_beta_procrustes_data <- beta_procrustes %>%
  dplyr::filter(method1_index >= method2_index)

p_beta_procrustes <- ggplot(p_beta_procrustes_data, aes(x = Method1, y = Method2, fill = procrustes_r)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 3.4, color = "grey10") +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC", limits = c(0, 1), na.value = "grey96") +
  labs(
    x = NULL,
    y = NULL,
    fill = "Procrustes r",
    caption = "Triangle shown; permutation p-values are saved in beta_pcoa_procrustes.csv."
  ) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    plot.caption = element_text(color = "grey35", hjust = 0)
  )
save_plot(p_beta_procrustes, "fig2b_beta_pcoa_procrustes", width = 6.3, height = 5.2)

############################################################
## 3. Differential abundance
############################################################

tax_prev <- colSums(counts_all > 0)
tax_total <- colSums(counts_all)
da_keep <- tax_prev >= ceiling(da_min_prevalence_fraction * nrow(counts_all)) &
  tax_total >= da_min_total_count

counts_da <- counts_all[, da_keep, drop = FALSE]
rel_da <- sweep(counts_da, 1, rowSums(counts_da), FUN = "/")

feature_info <- tibble::tibble(
  Taxon = colnames(counts_da),
  total_count = colSums(counts_da),
  prevalence_n = colSums(counts_da > 0),
  prevalence_fraction = colMeans(counts_da > 0),
  mean_relative_abundance = colMeans(rel_da)
) %>%
  dplyr::left_join(tax_df, by = "Taxon")

feature_filter_summary <- tibble::tibble(
  metric = c(
    "samples",
    "taxa_before_filter",
    "taxa_after_filter",
    "min_prevalence_fraction",
    "min_prevalence_n",
    "min_total_count",
    "zero_fraction_after_filter"
  ),
  value = c(
    nrow(counts_all),
    ncol(counts_all),
    ncol(counts_da),
    da_min_prevalence_fraction,
    ceiling(da_min_prevalence_fraction * nrow(counts_all)),
    da_min_total_count,
    mean(counts_da == 0)
  )
)
write.csv(feature_filter_summary, file.path(tab_dir, "da_feature_filter_summary.csv"), row.names = FALSE)
write.csv(feature_info, file.path(tab_dir, "da_feature_info_filtered_taxa.csv"), row.names = FALSE)

adjusted_lm_matrix <- function(mat, metadata) {
  metadata <- metadata[rownames(mat), , drop = FALSE]
  full_design <- stats::model.matrix(~ arm + age + day, data = metadata)
  reduced_design <- stats::model.matrix(~ age + day, data = metadata)

  full_fit <- stats::lm.fit(full_design, mat)
  reduced_fit <- stats::lm.fit(reduced_design, mat)

  rss_full <- colSums(full_fit$residuals^2)
  rss_reduced <- colSums(reduced_fit$residuals^2)
  df_arm <- full_fit$rank - reduced_fit$rank
  df_residual <- nrow(mat) - full_fit$rank

  rss_drop <- pmax(rss_reduced - rss_full, 0)
  ms_arm <- rss_drop / df_arm
  ms_residual <- rss_full / df_residual
  statistic <- ms_arm / ms_residual
  p_value <- stats::pf(statistic, df_arm, df_residual, lower.tail = FALSE)

  no_variance <- !is.finite(statistic) | ms_residual <= 0
  statistic[no_variance] <- NA_real_
  p_value[no_variance] <- NA_real_

  out <- tibble::tibble(
    Taxon = colnames(mat),
    df_arm = df_arm,
    df_residual = df_residual,
    statistic = statistic,
    p_value = p_value
  )

  groups <- levels(metadata$arm)
  for (lev in groups) {
    out[[paste0("mean_", lev)]] <- colMeans(mat[metadata$arm == lev, , drop = FALSE])
  }

  if (length(groups) == 2) {
    out$unadjusted_mean_difference <- out[[paste0("mean_", groups[2])]] - out[[paste0("mean_", groups[1])]]
  }

  coef_names <- colnames(full_design)
  arm_coef <- grep("^arm", coef_names, value = TRUE)
  if (length(arm_coef) == 1) {
    coef_idx <- match(arm_coef, coef_names)
    out$coefficient_name <- arm_coef
    out$coefficient_estimate <- full_fit$coefficients[coef_idx, ]
  } else {
    out$coefficient_name <- NA_character_
    out$coefficient_estimate <- NA_real_
  }

  out
}

da_transforms <- list(AHC = angular_hellinger_normalization(counts_da))
for (pc_name in names(pseudo_counts)) {
  da_transforms[[paste0("CLR_pc", pc_name)]] <- clr_transform(counts_da, pseudo_counts[[pc_name]])
}

da_anova <- purrr::imap_dfr(da_transforms, function(mat, method_name) {
  adjusted_lm_matrix(mat, meta) %>%
    dplyr::mutate(method = method_name, MethodLabel = clean_method_label(method_name), .before = 1)
}) %>%
  dplyr::group_by(method) %>%
  dplyr::mutate(
    padj = p.adjust(p_value, method = da_p_adjust_method),
    significant = !is.na(padj) & padj < alpha_level
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(feature_info, by = "Taxon")

write.csv(da_anova, file.path(tab_dir, "da_omnibus_lm_adjusted_AHC_CLR.csv"), row.names = FALSE)

top_ahc_anova <- da_anova %>%
  dplyr::filter(method == "AHC") %>%
  dplyr::arrange(padj, p_value) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::ungroup()
write.csv(top_ahc_anova, file.path(tab_dir, "da_top30_AHC_adjusted.csv"), row.names = FALSE)

run_deseq2_omnibus <- function(counts, metadata) {
  empty <- tibble::tibble(
    method = character(),
    MethodLabel = character(),
    Taxon = character(),
    baseMean = numeric(),
    statistic = numeric(),
    p_value = numeric(),
    padj = numeric(),
    significant = logical()
  )
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    message("DESeq2 is not installed; skipping DESeq2 comparison.")
    return(empty)
  }

  count_mat <- t(round(counts))
  storage.mode(count_mat) <- "integer"
  coldata <- metadata[rownames(counts), , drop = FALSE]
  coldata$arm <- droplevels(coldata$arm)
  coldata$age <- as.numeric(coldata$age)
  coldata$day <- as.numeric(coldata$day)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = ~ age + day + arm
  )
  dds <- dds[rowSums(DESeq2::counts(dds)) >= da_min_total_count, ]

  dds <- tryCatch(
    DESeq2::DESeq(dds, test = "LRT", reduced = ~ age + day, sfType = "poscounts", fitType = "parametric", quiet = TRUE),
    error = function(e) {
      message("DESeq2 parametric fit failed; retrying with local fit. Original error: ", conditionMessage(e))
      DESeq2::DESeq(dds, test = "LRT", reduced = ~ age + day, sfType = "poscounts", fitType = "local", quiet = TRUE)
    }
  )

  res <- DESeq2::results(dds, alpha = alpha_level, pAdjustMethod = da_p_adjust_method)
  res_df <- as.data.frame(res)
  res_df$Taxon <- rownames(res_df)
  tibble::as_tibble(res_df) %>%
    dplyr::transmute(
      method = "DESeq2_LRT",
      MethodLabel = "DESeq2",
      Taxon,
      baseMean,
      statistic = stat,
      p_value = pvalue,
      padj,
      significant = !is.na(padj) & padj < alpha_level
    ) %>%
    dplyr::left_join(feature_info, by = "Taxon")
}

da_deseq2 <- run_deseq2_omnibus(counts_da, meta)
write.csv(da_deseq2, file.path(tab_dir, "da_deseq2_omnibus_lrt_adjusted.csv"), row.names = FALSE)

run_lefse_one_vs_rest <- function(counts, metadata) {
  empty <- tibble::tibble(
    method = character(),
    MethodLabel = character(),
    Taxon = character(),
    target_arm = character(),
    lda_score = numeric(),
    significant = logical()
  )
  if (!requireNamespace("lefser", quietly = TRUE)) {
    message("The lefser package is not installed; skipping LEfSe analysis.")
    return(empty)
  }

  counts_feature_sample <- t(counts)
  storage.mode(counts_feature_sample) <- "numeric"
  coldata <- metadata[rownames(counts), , drop = FALSE]

  purrr::map_dfr(levels(coldata$arm), function(target_arm) {
    coldata_run <- coldata
    coldata_run$ArmBinary <- factor(
      ifelse(coldata_run$arm == target_arm, target_arm, "Other"),
      levels = c("Other", target_arm)
    )
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_feature_sample),
      colData = coldata_run
    )
    se_ra <- lefser::relativeAb(se)
    set.seed(20260608)
    res <- tryCatch(
      lefser::lefser(
        se_ra,
        classCol = "ArmBinary",
        kruskal.threshold = alpha_level,
        wilcox.threshold = alpha_level,
        lda.threshold = lefse_lda_threshold,
        method = lefse_p_adjust_method
      ),
      error = function(e) {
        message("lefser failed for ", target_arm, " vs rest: ", conditionMessage(e))
        data.frame(features = character(0), scores = numeric(0))
      }
    )
    if (nrow(res) == 0) {
      return(empty)
    }
    tibble::as_tibble(res) %>%
      dplyr::transmute(
        method = "LEfSe",
        MethodLabel = "LEfSe",
        Taxon = features,
        target_arm = target_arm,
        lda_score = scores,
        significant = TRUE
      )
  }) %>%
    dplyr::left_join(feature_info, by = "Taxon")
}

da_lefse <- run_lefse_one_vs_rest(counts_da, meta)
write.csv(da_lefse, file.path(tab_dir, "da_lefse_one_vs_rest_unadjusted.csv"), row.names = FALSE)

da_omnibus_counts <- dplyr::bind_rows(
  da_anova %>% dplyr::select(method, MethodLabel, Taxon, padj, significant),
  da_deseq2 %>% dplyr::select(method, MethodLabel, Taxon, padj, significant),
  da_lefse %>%
    dplyr::distinct(method, MethodLabel, Taxon, significant) %>%
    dplyr::mutate(padj = NA_real_)
) %>%
  dplyr::group_by(method, MethodLabel) %>%
  dplyr::summarise(
    n_tested = dplyr::n_distinct(Taxon),
    n_significant = dplyr::n_distinct(Taxon[significant]),
    .groups = "drop"
  )
write.csv(da_omnibus_counts, file.path(tab_dir, "da_omnibus_discovery_counts.csv"), row.names = FALSE)

da_sets <- tibble::tibble(Taxon = colnames(counts_da)) %>%
  dplyr::mutate(
    AHC = Taxon %in% (da_anova %>% dplyr::filter(method == "AHC", significant) %>% dplyr::pull(Taxon)),
    CLR_pc0.01 = Taxon %in% (da_anova %>% dplyr::filter(method == "CLR_pc0.01", significant) %>% dplyr::pull(Taxon)),
    CLR_pc0.5 = Taxon %in% (da_anova %>% dplyr::filter(method == "CLR_pc0.5", significant) %>% dplyr::pull(Taxon)),
    CLR_pc1 = Taxon %in% (da_anova %>% dplyr::filter(method == "CLR_pc1", significant) %>% dplyr::pull(Taxon)),
    DESeq2 = Taxon %in% (da_deseq2 %>% dplyr::filter(significant) %>% dplyr::pull(Taxon)),
    LEfSe = Taxon %in% unique(da_lefse$Taxon)
  )
write.csv(da_sets, file.path(tab_dir, "da_omnibus_overlap_membership.csv"), row.names = FALSE)

save_limma_venn <- function(membership, set_cols, labels, filename, main) {
  mat <- as.matrix(membership[, set_cols, drop = FALSE])
  storage.mode(mat) <- "logical"
  colnames(mat) <- set_cols
  vc <- limma::vennCounts(mat)
  vc_df <- as.data.frame(unclass(vc))
  write.csv(vc_df, file.path(tab_dir, paste0(filename, "_counts.csv")), row.names = FALSE)

  png_file <- file.path(fig_dir, paste0(filename, ".png"))
  pdf_file <- file.path(fig_dir, paste0(filename, ".pdf"))
  cols <- method_palette[set_cols]
  cols[is.na(cols)] <- grDevices::hcl.colors(sum(is.na(cols)), "Dark 3")

  grDevices::png(png_file, width = 7.5, height = 7.0, units = "in", res = 320)
  limma::vennDiagram(vc, names = labels, circle.col = cols, counts.col = "black", cex = c(1.15, 0.95, 0.72), main = main)
  grDevices::dev.off()

  grDevices::pdf(pdf_file, width = 7.5, height = 7.0)
  limma::vennDiagram(vc, names = labels, circle.col = cols, counts.col = "black", cex = c(1.15, 0.95, 0.72), main = main)
  grDevices::dev.off()

  invisible(vc_df)
}

save_three_set_venn <- function(membership, set_cols, labels, filename) {
  if (length(set_cols) != 3L || length(labels) != 3L) {
    stop("save_three_set_venn() expects exactly three sets and three labels.")
  }

  dat <- as.data.frame(membership[, set_cols, drop = FALSE], check.names = FALSE)
  dat[] <- lapply(dat, as.logical)
  a <- dat[[set_cols[1]]]
  b <- dat[[set_cols[2]]]
  c <- dat[[set_cols[3]]]

  vc_df <- tibble::tibble(
    !!set_cols[1] := as.integer(c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)),
    !!set_cols[2] := as.integer(c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)),
    !!set_cols[3] := as.integer(c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)),
    Counts = c(
      sum(!a & !b & !c),
      sum(!a & !b & c),
      sum(!a & b & !c),
      sum(!a & b & c),
      sum(a & !b & !c),
      sum(a & !b & c),
      sum(a & b & !c),
      sum(a & b & c)
    )
  )
  write.csv(vc_df, file.path(tab_dir, paste0(filename, "_counts.csv")), row.names = FALSE)

  region_counts <- c(
    a_only = sum(a & !b & !c),
    b_only = sum(!a & b & !c),
    c_only = sum(!a & !b & c),
    ab_only = sum(a & b & !c),
    ac_only = sum(a & !b & c),
    bc_only = sum(!a & b & c),
    abc = sum(a & b & c)
  )

  draw_venn <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(0.05, 0.05, 0.05, 0.05), xpd = NA)
    centers <- rbind(c(3.45, 5.35), c(6.55, 5.35), c(5.00, 3.15))
    radius <- 2.65
    circle_points <- function(center, n = 720) {
      theta <- seq(0, 2 * pi, length.out = n)
      cbind(center[1] + radius * cos(theta), center[2] + radius * sin(theta))
    }
    inside_all <- function(points) {
      inside <- vapply(seq_len(nrow(centers)), function(i) {
        (points[, 1] - centers[i, 1])^2 + (points[, 2] - centers[i, 2])^2 <= radius^2 + 1e-8
      }, logical(nrow(points)))
      rowSums(inside) == nrow(centers)
    }

    graphics::plot(
      NA,
      xlim = c(0, 10),
      ylim = c(0, 9),
      asp = 1,
      axes = FALSE,
      ann = FALSE,
      bty = "n"
    )

    for (i in seq_len(nrow(centers))) {
      pts <- circle_points(centers[i, ])
      graphics::polygon(pts[, 1], pts[, 2], col = "#EAF1F8", border = NA)
    }
    triple_pts <- do.call(rbind, lapply(seq_len(nrow(centers)), function(i) {
      pts <- circle_points(centers[i, ])
      pts[inside_all(pts), , drop = FALSE]
    }))
    triple_center <- colMeans(triple_pts)
    triple_order <- order(atan2(triple_pts[, 2] - triple_center[2], triple_pts[, 1] - triple_center[1]))
    triple_pts <- triple_pts[triple_order, , drop = FALSE]
    graphics::polygon(triple_pts[, 1], triple_pts[, 2], col = "#3479B9", border = NA)
    for (i in seq_len(nrow(centers))) {
      pts <- circle_points(centers[i, ])
      graphics::lines(pts[, 1], pts[, 2], col = "black", lwd = 2.2)
    }

    graphics::text(1.25, 8.35, labels[1], cex = 1.45, adj = c(0, 0.5))
    graphics::text(8.75, 8.35, labels[2], cex = 1.45, adj = c(1, 0.5))
    graphics::text(5.00, 0.32, labels[3], cex = 1.45, adj = c(0.5, 0.5))

    graphics::text(2.55, 5.55, region_counts["a_only"], cex = 1.35)
    graphics::text(7.45, 5.55, region_counts["b_only"], cex = 1.35)
    graphics::text(5.00, 1.80, region_counts["c_only"], cex = 1.35)
    graphics::text(5.00, 5.88, region_counts["ab_only"], cex = 1.35)
    graphics::text(3.85, 4.05, region_counts["ac_only"], cex = 1.35)
    graphics::text(6.15, 4.05, region_counts["bc_only"], cex = 1.35)
    graphics::text(5.00, 4.55, region_counts["abc"], cex = 1.35)
  }

  png_file <- file.path(fig_dir, paste0(filename, ".png"))
  pdf_file <- file.path(fig_dir, paste0(filename, ".pdf"))

  grDevices::png(png_file, width = 7.0, height = 6.1, units = "in", res = 320)
  draw_venn()
  grDevices::dev.off()

  grDevices::pdf(pdf_file, width = 7.0, height = 6.1)
  draw_venn()
  grDevices::dev.off()

  invisible(c(png_file, pdf_file))
}

save_upset_plot <- function(membership, set_cols, labels, filename, title, n_intersections = 18) {
  dat <- as.data.frame(membership[, set_cols, drop = FALSE], check.names = FALSE)
  dat[] <- lapply(dat, as.logical)

  combo_mat <- as.matrix(dat)
  storage.mode(combo_mat) <- "integer"
  detected <- rowSums(combo_mat) > 0

  if (!any(detected)) {
    empty_intersections <- tibble::tibble(Intersection = character(), size = integer(), methods = character())
    write.csv(empty_intersections, file.path(tab_dir, paste0(filename, "_top_intersections.csv")), row.names = FALSE)
    p_empty <- ggplot() +
      annotate("text", x = 0, y = 0, label = paste0(title, "\nNo detected taxa"), size = 4) +
      theme_void()
    save_plot(p_empty, filename, width = 8.0, height = 4.5)
    return(invisible(NULL))
  }

  combo_key <- apply(combo_mat[detected, , drop = FALSE], 1, paste0, collapse = "")
  combo_counts <- sort(table(combo_key), decreasing = TRUE)
  top_keys <- names(combo_counts)[seq_len(min(n_intersections, length(combo_counts)))]

  intersections <- tibble::tibble(
    Intersection = factor(seq_along(top_keys), levels = seq_along(top_keys)),
    combo_key = top_keys,
    size = as.integer(combo_counts[top_keys])
  )

  bit_df <- as.data.frame(do.call(rbind, strsplit(top_keys, split = "")), stringsAsFactors = FALSE)
  colnames(bit_df) <- set_cols
  matrix_df <- bit_df %>%
    dplyr::mutate(Intersection = factor(seq_along(top_keys), levels = seq_along(top_keys))) %>%
    tidyr::pivot_longer(dplyr::all_of(set_cols), names_to = "set_col", values_to = "active") %>%
    dplyr::mutate(
      active = active == "1",
      set_index = match(set_col, set_cols),
      SetLabel = labels[set_index]
    )

  set_sizes <- tibble::tibble(
    set_col = set_cols,
    SetLabel = labels,
    set_index = seq_along(set_cols),
    size = colSums(dat)
  )

  segment_df <- matrix_df %>%
    dplyr::filter(active) %>%
    dplyr::group_by(Intersection) %>%
    dplyr::summarise(
      ymin = min(set_index),
      ymax = max(set_index),
      .groups = "drop"
    ) %>%
    dplyr::filter(ymin < ymax)

  top_intersections <- matrix_df %>%
    dplyr::filter(active) %>%
    dplyr::left_join(intersections, by = "Intersection") %>%
    dplyr::group_by(Intersection, size) %>%
    dplyr::summarise(methods = paste(SetLabel, collapse = " + "), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(size))
  write.csv(top_intersections, file.path(tab_dir, paste0(filename, "_top_intersections.csv")), row.names = FALSE)

  band_df <- tibble::tibble(
    ymin = seq_along(set_cols) - 0.5,
    ymax = seq_along(set_cols) + 0.5,
    shade = rep(c("#FFFFFF", "#F7F7F7"), length.out = length(set_cols))
  )

  p_main <- ggplot(intersections, aes(x = Intersection, y = size)) +
    geom_col(width = 0.62, fill = "#303030") +
    geom_text(aes(label = size), vjust = -0.28, size = 2.65, color = "grey20") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Intersection size") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = grid::unit(c(4, 5, 1, 5), "pt")
    )

  p_matrix <- ggplot(matrix_df, aes(x = Intersection, y = set_index)) +
    geom_rect(
      data = band_df,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = shade),
      inherit.aes = FALSE,
      color = NA
    ) +
    geom_segment(
      data = segment_df,
      aes(x = Intersection, xend = Intersection, y = ymin, yend = ymax),
      inherit.aes = FALSE,
      color = "#252525",
      linewidth = 0.45
    ) +
    geom_point(color = "#DADADA", size = 2.7) +
    geom_point(data = dplyr::filter(matrix_df, active), color = "#252525", size = 2.75) +
    scale_y_continuous(
      breaks = seq_along(set_cols),
      labels = NULL,
      limits = c(0.5, length(set_cols) + 0.5),
      expand = c(0, 0)
    ) +
    scale_fill_identity() +
    labs(x = "Largest intersections sorted by size", y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = grid::unit(c(1, 5, 4, 5), "pt")
    )

  p_sets <- ggplot(set_sizes) +
    geom_rect(
      aes(
        xmin = 0,
        xmax = size,
        ymin = set_index - 0.23,
        ymax = set_index + 0.23,
        fill = set_col
      )
    ) +
    geom_text(aes(x = size, y = set_index, label = size), hjust = 1.08, size = 2.55, color = "white") +
    scale_fill_manual(values = method_palette[set_cols], guide = "none") +
    scale_y_continuous(
      breaks = seq_along(set_cols),
      labels = labels,
      limits = c(0.5, length(set_cols) + 0.5),
      expand = c(0, 0)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = "Detected taxa per method", y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y = element_text(color = "grey15"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = grid::unit(c(1, 2, 4, 5), "pt")
    )

  spacer <- cowplot::ggdraw() +
    cowplot::draw_label(title, x = 0.04, y = 0.92, hjust = 0, vjust = 1, size = 10, fontface = "bold", color = "grey15") +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  right_panels <- cowplot::align_plots(p_main, p_matrix, align = "v", axis = "l")
  top_row <- cowplot::plot_grid(spacer, right_panels[[1]], ncol = 2, rel_widths = c(0.32, 1), align = "h")
  bottom_row <- cowplot::plot_grid(p_sets, right_panels[[2]], ncol = 2, rel_widths = c(0.32, 1), align = "h")
  combined <- cowplot::plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(0.62, 0.38), align = "v") +
    theme(plot.background = element_rect(fill = "white", color = NA))

  png_file <- file.path(fig_dir, paste0(filename, ".png"))
  pdf_file <- file.path(fig_dir, paste0(filename, ".pdf"))
  ggplot2::ggsave(png_file, combined, width = 10.0, height = 5.8, dpi = 320)
  ggplot2::ggsave(pdf_file, combined, width = 10.0, height = 5.8)

  invisible(c(png_file, pdf_file))
}

save_limma_venn(
  da_sets,
  set_cols = c("AHC", "CLR_pc0.01", "CLR_pc0.5", "CLR_pc1", "LEfSe"),
  labels = c("AHC", "CLR\npc = 0.01", "CLR\npc = 0.5", "CLR\npc = 1", "LEfSe"),
  filename = "fig3_da_venn_AHC_CLR_LEfSe",
  main = "Adjusted arm DA overlap: AHC, CLR pseudo-counts, LEfSe"
)

save_upset_plot(
  da_sets,
  set_cols = c("AHC", "CLR_pc0.01", "CLR_pc0.5", "CLR_pc1", "LEfSe"),
  labels = c("AHC", "CLR pc = 0.01", "CLR pc = 0.5", "CLR pc = 1", "LEfSe"),
  filename = "fig3b_da_upset_AHC_CLR_LEfSe",
  title = "AHC + CLR + LEfSe; FDR < 0.05"
)

save_limma_venn(
  da_sets,
  set_cols = c("AHC", "CLR_pc0.01", "CLR_pc0.5", "CLR_pc1", "DESeq2"),
  labels = c("AHC", "CLR\npc = 0.01", "CLR\npc = 0.5", "CLR\npc = 1", "DESeq2"),
  filename = "fig4_da_venn_AHC_CLR_DESeq2",
  main = "Adjusted arm DA overlap: AHC, CLR pseudo-counts, DESeq2"
)

save_upset_plot(
  da_sets,
  set_cols = c("AHC", "CLR_pc0.01", "CLR_pc0.5", "CLR_pc1", "DESeq2"),
  labels = c("AHC", "CLR pc = 0.01", "CLR pc = 0.5", "CLR pc = 1", "DESeq2"),
  filename = "fig4b_da_upset_AHC_CLR_DESeq2",
  title = "AHC + CLR + DESeq2; FDR < 0.05"
)

save_three_set_venn(
  da_sets,
  set_cols = c("CLR_pc0.01", "CLR_pc0.5", "CLR_pc1"),
  labels = c("CLR (pc = 0.01)", "CLR (pc = 0.5)", "CLR (pc = 1)"),
  filename = "fig5_da_venn_CLR_pseudocounts"
)

############################################################
## Run summary
############################################################

summary_lines <- c(
  "RealDataGut AHC analysis completed.",
  paste("Run time:", as.character(Sys.time())),
  paste("Analysis directory:", analysis_dir),
  paste("Output directory:", out_dir),
  paste("Outcome:", "arm"),
  paste("Confounders:", "age, day"),
  paste("Day coding:", "negative = before surgery; zero = surgery day; positive = after surgery"),
  paste("Arm counts:", paste(names(table(meta$arm)), as.integer(table(meta$arm)), sep = "=", collapse = ", ")),
  paste("Complete-case samples:", nrow(counts_all)),
  paste("Excluded samples missing arm/age/day:", nrow(excluded_samples)),
  paste("All selected taxa:", ncol(counts_all)),
  paste("DA filtered taxa:", ncol(counts_da)),
  paste("Pseudo-counts:", paste(pseudo_counts, collapse = ", ")),
  paste("AHC/CLR DA model:", "feature ~ arm + age + day compared to feature ~ age + day"),
  paste("DESeq2 DA model:", "design = ~ age + day + arm; LRT reduced = ~ age + day"),
  paste("AHC/CLR and DESeq2 DA adjustment:", da_p_adjust_method, "adjusted p-values; threshold =", alpha_level),
  paste("LEfSe:", "lefser package, one-vs-rest arm runs; no age/day adjustment available in this LEfSe step; p adjustment =", lefse_p_adjust_method, "; LDA threshold =", lefse_lda_threshold),
  "",
  "Key tables:",
  file.path(tab_dir, "data_summary.csv"),
  file.path(tab_dir, "alpha_arm_summary.csv"),
  file.path(tab_dir, "alpha_arm_lm_adjusted_tests.csv"),
  file.path(tab_dir, "beta_permanova_adjusted.csv"),
  file.path(tab_dir, "beta_pcoa_procrustes.csv"),
  file.path(tab_dir, "da_omnibus_lm_adjusted_AHC_CLR.csv"),
  file.path(tab_dir, "da_deseq2_omnibus_lrt_adjusted.csv"),
  file.path(tab_dir, "da_lefse_one_vs_rest_unadjusted.csv"),
  "",
  "Key figures:",
  file.path(fig_dir, "fig1_alpha_diversity_by_arm.png"),
  file.path(fig_dir, "fig1b_alpha_method_procrustes.png"),
  file.path(fig_dir, "fig2_beta_pcoa_by_method.png"),
  file.path(fig_dir, "fig2b_beta_pcoa_procrustes.png"),
  file.path(fig_dir, "fig3_da_venn_AHC_CLR_LEfSe.png"),
  file.path(fig_dir, "fig3b_da_upset_AHC_CLR_LEfSe.png"),
  file.path(fig_dir, "fig4_da_venn_AHC_CLR_DESeq2.png"),
  file.path(fig_dir, "fig4b_da_upset_AHC_CLR_DESeq2.png")
)

writeLines(summary_lines, con = file.path(out_dir, "run_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
