############################################################
## Real ocean data analysis for Arc-Hellinger normalization
##
## Outputs:
##   Analysis/RealDataOcean/real_data_ocean_results/tables
##   Analysis/RealDataOcean/real_data_ocean_results/figures
##
## The AHC implementation is sourced from R/ArcH.R and is not
## modified here.
############################################################

options(stringsAsFactors = FALSE, width = 140)
set.seed(20260604)

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

out_dir <- file.path(analysis_dir, "real_data_ocean_results")
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

layer_levels <- c("SRF", "DCM", "MES")
layer_labels <- c(SRF = "Surface", DCM = "DCM", MES = "Mesopelagic")
layer_palette <- c(SRF = "#0072B2", DCM = "#D55E00", MES = "#009E73")
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

ocean_file <- file.path(analysis_dir, "Ocean.RData")
if (!file.exists(ocean_file)) {
  stop("Cannot find Ocean.RData in: ", analysis_dir)
}
load(ocean_file)
if (!exists("ps") || !inherits(ps, "phyloseq")) {
  stop("Ocean.RData must contain a phyloseq object named ps.")
}

ps_layer <- phyloseq::subset_samples(ps, Layer %in% layer_levels)
ps_layer <- phyloseq::prune_samples(phyloseq::sample_sums(ps_layer) > 0, ps_layer)
ps_layer <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_layer) > 0, ps_layer)

meta <- as(phyloseq::sample_data(ps_layer), "data.frame")
meta$SampleID <- rownames(meta)
meta$Layer <- factor(meta$Layer, levels = layer_levels)
meta$LayerLabel <- factor(layer_labels[as.character(meta$Layer)], levels = unname(layer_labels[layer_levels]))
layer_n <- table(meta$Layer)
layer_axis_labels <- setNames(
  paste0(layer_labels[layer_levels], "\n(n = ", as.integer(layer_n[layer_levels]), ")"),
  layer_levels
)
meta$LayerLabelN <- factor(layer_axis_labels[as.character(meta$Layer)], levels = layer_axis_labels[layer_levels])

counts_all <- as(phyloseq::otu_table(ps_layer), "matrix")
if (phyloseq::taxa_are_rows(ps_layer)) {
  counts_all <- t(counts_all)
}
storage.mode(counts_all) <- "numeric"
counts_all <- counts_all[meta$SampleID, , drop = FALSE]
counts_all <- counts_all[, colSums(counts_all) > 0, drop = FALSE]

tax_df <- as(phyloseq::tax_table(ps_layer), "matrix")
tax_df <- as.data.frame(tax_df, stringsAsFactors = FALSE, check.names = FALSE)
tax_df$Taxon <- rownames(tax_df)
tax_df <- tax_df %>% dplyr::relocate(Taxon)

lib_size <- rowSums(counts_all)
rel_all <- sweep(counts_all, 1, lib_size, FUN = "/")

data_summary <- tibble::tibble(
  metric = c(
    "samples_original",
    "samples_used_SRF_DCM_MES",
    "taxa_original",
    "taxa_used_nonzero_in_selected_layers",
    "zero_fraction_selected",
    "median_library_size",
    "min_library_size",
    "max_library_size",
    "SRF_n",
    "DCM_n",
    "MES_n"
  ),
  value = c(
    phyloseq::nsamples(ps),
    nrow(counts_all),
    phyloseq::ntaxa(ps),
    ncol(counts_all),
    mean(counts_all == 0),
    median(lib_size),
    min(lib_size),
    max(lib_size),
    sum(meta$Layer == "SRF"),
    sum(meta$Layer == "DCM"),
    sum(meta$Layer == "MES")
  )
)
write.csv(data_summary, file.path(tab_dir, "data_summary.csv"), row.names = FALSE)

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
  Layer = meta$Layer,
  LayerLabel = meta$LayerLabel,
  LayerLabelN = meta$LayerLabelN,
  LibrarySize = lib_size,
  Richness = vegan::specnumber(counts_all),
  AHC_D = ahc_norm,
  AHC_D_norm = ahc_norm / ahc_bound,
  AHC_E = 1 - ahc_norm / ahc_bound,
  Shannon_drop0 = shannon_from_counts(counts_all),
  Shannon_drop0_evenness = Shannon_drop0 / log(p_all)
)

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
  dplyr::select(SampleID, Layer, LayerLabel, LayerLabelN, dplyr::all_of(alpha_norm_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(alpha_norm_cols),
    names_to = "method",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    MethodLabel = factor(clean_method_label(method), levels = clean_method_label(alpha_norm_cols))
  )

alpha_summary <- alpha_long %>%
  dplyr::group_by(method, MethodLabel, Layer, LayerLabel) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q1 = quantile(value, 0.25, na.rm = TRUE),
    q3 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(alpha_summary, file.path(tab_dir, "alpha_layer_summary.csv"), row.names = FALSE)

alpha_kw <- alpha_long %>%
  dplyr::group_by(method, MethodLabel) %>%
  dplyr::summarise(
    statistic = unname(kruskal.test(value ~ Layer)$statistic),
    p_value = kruskal.test(value ~ Layer)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(q_value = p.adjust(p_value, method = "BH"))
write.csv(alpha_kw, file.path(tab_dir, "alpha_layer_kruskal_tests.csv"), row.names = FALSE)

pairwise_wilcox_df <- function(dat) {
  pairs <- utils::combn(levels(dat$Layer), 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pair) {
    sub <- dat %>% dplyr::filter(Layer %in% pair)
    wt <- wilcox.test(value ~ droplevels(Layer), data = sub, exact = FALSE)
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
write.csv(alpha_pairwise, file.path(tab_dir, "alpha_pairwise_wilcox_tests.csv"), row.names = FALSE)

anova_layer_df <- function(dat) {
  fit <- stats::lm(value ~ Layer, data = dat)
  tab <- stats::anova(fit)
  tibble::tibble(
    df_layer = tab$Df[1],
    df_residual = tab$Df[2],
    statistic = tab$`F value`[1],
    p_value = tab$`Pr(>F)`[1]
  )
}

alpha_anova <- alpha_long %>%
  dplyr::group_by(method, MethodLabel) %>%
  dplyr::group_modify(~ anova_layer_df(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(q_value = p.adjust(p_value, method = "BH"))
write.csv(alpha_anova, file.path(tab_dir, "alpha_layer_anova_tests.csv"), row.names = FALSE)

alpha_anova_lab <- alpha_anova %>%
  dplyr::mutate(label = paste0("F = ", formatC(statistic, format = "fg", digits = 3), "\n", format_p(p_value)))

p_alpha <- ggplot(alpha_long, aes(x = LayerLabelN, y = value, color = Layer, fill = Layer)) +
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
    data = alpha_anova_lab,
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
  scale_color_manual(values = layer_palette, guide = "none") +
  scale_fill_manual(values = layer_palette, guide = "none") +
  labs(x = NULL, y = "Normalized alpha index") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text.x = element_text(size = 9.5),
    panel.spacing = grid::unit(0.85, "lines")
  )
save_plot(p_alpha, "fig1_alpha_diversity_by_layer", width = 10.5, height = 6.7)

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
    label = ifelse(
      is.na(p_value),
      sprintf("%.2f", procrustes_r),
      sprintf("%.2f", procrustes_r)
    )
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

run_adonis_layer <- function(d, metadata, permutations = 999) {
  tab <- vegan::adonis2(d ~ Layer, data = metadata, permutations = permutations)
  tibble::tibble(
    df = tab$Df[1],
    sumsq = tab$SumOfSqs[1],
    r2 = tab$R2[1],
    statistic = tab$F[1],
    p_value = tab$`Pr(>F)`[1]
  )
}

beta_permanova <- purrr::imap_dfr(beta_dists, function(d, method_name) {
  run_adonis_layer(d, meta) %>%
    dplyr::mutate(method = method_name, MethodLabel = clean_method_label(method_name), .before = 1)
}) %>%
  dplyr::mutate(q_value = p.adjust(p_value, method = "BH"))
write.csv(beta_permanova, file.path(tab_dir, "beta_permanova_layer.csv"), row.names = FALSE)

pairwise_adonis <- function(d, metadata, method_name) {
  pairs <- utils::combn(layer_levels, 2, simplify = FALSE)
  dmat <- as.matrix(d)
  purrr::map_dfr(pairs, function(pair) {
    keep <- metadata$Layer %in% pair
    sub_meta <- metadata[keep, , drop = FALSE]
    sub_meta$Layer <- droplevels(sub_meta$Layer)
    sub_d <- stats::as.dist(dmat[keep, keep, drop = FALSE])
    run_adonis_layer(sub_d, sub_meta) %>%
      dplyr::mutate(
        method = method_name,
        MethodLabel = clean_method_label(method_name),
        group1 = pair[1],
        group2 = pair[2],
        contrast = format_contrast(pair[1], pair[2]),
        .before = 1
      )
  })
}

beta_pairwise_permanova <- purrr::imap_dfr(beta_dists, ~ pairwise_adonis(.x, meta, .y)) %>%
  dplyr::group_by(method) %>%
  dplyr::mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup()
write.csv(beta_pairwise_permanova, file.path(tab_dir, "beta_pairwise_permanova.csv"), row.names = FALSE)

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
  dplyr::left_join(meta %>% dplyr::select(SampleID, Layer, LayerLabel), by = "SampleID") %>%
  dplyr::mutate(MethodLabel = factor(MethodLabel, levels = clean_method_label(names(beta_dists))))

pcoa_lab <- pcoa_df %>%
  dplyr::group_by(MethodLabel, AxisLabel) %>%
  dplyr::summarise(
    x = min(Axis1, na.rm = TRUE),
    y = max(Axis2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    beta_permanova %>%
      dplyr::select(MethodLabel, r2, statistic, p_value),
    by = "MethodLabel"
  ) %>%
  dplyr::mutate(
    label = paste0(
      AxisLabel,
      "\nPERMANOVA R2 = ", sprintf("%.3f", r2),
      "\n", format_p(p_value)
    )
  )

p_beta_pcoa <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2)) +
  stat_ellipse(aes(color = Layer, group = Layer), linewidth = 0.45, alpha = 0.78, show.legend = FALSE) +
  geom_point(aes(fill = Layer), shape = 21, size = 2.0, alpha = 0.86, color = "white", stroke = 0.22) +
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
  scale_color_manual(values = layer_palette, labels = layer_labels[layer_levels]) +
  scale_fill_manual(values = layer_palette, labels = layer_labels[layer_levels]) +
  labs(x = "PCoA coordinate 1", y = "PCoA coordinate 2", fill = "Layer") +
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
    label = ifelse(
      is.na(p_value),
      sprintf("%.2f", procrustes_r),
      sprintf("%.2f", procrustes_r)
    )
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

oneway_anova_matrix <- function(mat, group) {
  group <- factor(group, levels = layer_levels)
  keep <- !is.na(group)
  mat <- mat[keep, , drop = FALSE]
  group <- droplevels(group[keep])

  n <- nrow(mat)
  groups <- levels(group)
  g <- length(groups)
  grand <- colMeans(mat)

  ss_between <- rep(0, ncol(mat))
  ss_within <- rep(0, ncol(mat))
  group_means <- list()

  for (lev in groups) {
    xg <- mat[group == lev, , drop = FALSE]
    ng <- nrow(xg)
    mg <- colMeans(xg)
    group_means[[lev]] <- mg
    ss_between <- ss_between + ng * (mg - grand)^2
    ss_within <- ss_within + colSums(sweep(xg, 2, mg, FUN = "-")^2)
  }

  df_between <- g - 1
  df_within <- n - g
  ms_between <- ss_between / df_between
  ms_within <- ss_within / df_within
  statistic <- ms_between / ms_within
  p_value <- stats::pf(statistic, df_between, df_within, lower.tail = FALSE)

  no_variance <- !is.finite(statistic) | ms_within <= 0
  statistic[no_variance] <- NA_real_
  p_value[no_variance] <- NA_real_

  out <- tibble::tibble(
    Taxon = colnames(mat),
    df_layer = df_between,
    df_residual = df_within,
    statistic = statistic,
    p_value = p_value
  )

  for (lev in groups) {
    out[[paste0("mean_", lev)]] <- group_means[[lev]]
  }

  out
}

da_transforms <- list(AHC = angular_hellinger_normalization(counts_da))
for (pc_name in names(pseudo_counts)) {
  da_transforms[[paste0("CLR_pc", pc_name)]] <- clr_transform(counts_da, pseudo_counts[[pc_name]])
}

da_anova <- purrr::imap_dfr(da_transforms, function(mat, method_name) {
  oneway_anova_matrix(mat, meta$Layer) %>%
    dplyr::mutate(method = method_name, MethodLabel = clean_method_label(method_name), .before = 1)
}) %>%
  dplyr::group_by(method) %>%
  dplyr::mutate(
    padj = p.adjust(p_value, method = da_p_adjust_method),
    significant = !is.na(padj) & padj < alpha_level
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(feature_info, by = "Taxon")

write.csv(da_anova, file.path(tab_dir, "da_omnibus_anova_AHC_CLR.csv"), row.names = FALSE)

top_ahc_anova <- da_anova %>%
  dplyr::filter(method == "AHC") %>%
  dplyr::arrange(padj, p_value) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::ungroup()
write.csv(top_ahc_anova, file.path(tab_dir, "da_top30_AHC_omnibus.csv"), row.names = FALSE)

run_deseq2_omnibus <- function(counts, metadata) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    message("DESeq2 is not installed; skipping DESeq2 comparison.")
    return(tibble::tibble())
  }
  count_mat <- t(round(counts))
  storage.mode(count_mat) <- "integer"
  coldata <- metadata[rownames(counts), , drop = FALSE]
  coldata$Layer <- factor(coldata$Layer, levels = layer_levels)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = ~ Layer
  )
  dds <- dds[rowSums(DESeq2::counts(dds)) >= da_min_total_count, ]

  dds <- tryCatch(
    DESeq2::DESeq(dds, test = "LRT", reduced = ~ 1, sfType = "poscounts", fitType = "parametric", quiet = TRUE),
    error = function(e) {
      message("DESeq2 parametric fit failed; retrying with local fit. Original error: ", conditionMessage(e))
      DESeq2::DESeq(dds, test = "LRT", reduced = ~ 1, sfType = "poscounts", fitType = "local", quiet = TRUE)
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
write.csv(da_deseq2, file.path(tab_dir, "da_deseq2_omnibus_lrt.csv"), row.names = FALSE)

run_lefse_one_vs_rest <- function(counts, metadata) {
  if (!requireNamespace("lefser", quietly = TRUE)) {
    stop("The lefser package is required for the LEfSe analysis.")
  }

  counts_feature_sample <- t(counts)
  storage.mode(counts_feature_sample) <- "numeric"
  coldata <- metadata[rownames(counts), , drop = FALSE]

  purrr::map_dfr(layer_levels, function(target_layer) {
    coldata_run <- coldata
    coldata_run$LayerBinary <- factor(
      ifelse(coldata_run$Layer == target_layer, target_layer, "Other"),
      levels = c("Other", target_layer)
    )
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts_feature_sample),
      colData = coldata_run
    )
    se_ra <- lefser::relativeAb(se)
    set.seed(20260604)
    res <- tryCatch(
      lefser::lefser(
        se_ra,
        classCol = "LayerBinary",
        kruskal.threshold = alpha_level,
        wilcox.threshold = alpha_level,
        lda.threshold = lefse_lda_threshold,
        method = lefse_p_adjust_method
      ),
      error = function(e) {
        message("lefser failed for ", target_layer, " vs rest: ", conditionMessage(e))
        data.frame(features = character(0), scores = numeric(0))
      }
    )
    if (nrow(res) == 0) {
      return(tibble::tibble())
    }
    tibble::as_tibble(res) %>%
      dplyr::transmute(
        method = "LEfSe",
        MethodLabel = "LEfSe",
        Taxon = features,
        target_layer = target_layer,
        lda_score = scores,
        significant = TRUE
      )
  }) %>%
    dplyr::left_join(feature_info, by = "Taxon")
}

da_lefse <- run_lefse_one_vs_rest(counts_da, meta)
write.csv(da_lefse, file.path(tab_dir, "da_lefse_one_vs_rest.csv"), row.names = FALSE)

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
  cols <- c("#AB2428", "#2F74B8", "#79B0D7", "#A0C8E8", "#009E73")

  grDevices::png(png_file, width = 7.5, height = 7.0, units = "in", res = 320)
  limma::vennDiagram(vc, names = labels, circle.col = cols, counts.col = "black", cex = c(1.15, 0.95, 0.72))
  grDevices::dev.off()

  grDevices::pdf(pdf_file, width = 7.5, height = 7.0)
  limma::vennDiagram(vc, names = labels, circle.col = cols, counts.col = "black", cex = c(1.15, 0.95, 0.72))
  grDevices::dev.off()

  invisible(vc_df)
}

save_upset_plot <- function(membership, set_cols, labels, filename, title, n_intersections = 18) {
  dat <- as.data.frame(membership[, set_cols, drop = FALSE], check.names = FALSE)
  dat[] <- lapply(dat, as.logical)

  combo_mat <- as.matrix(dat)
  storage.mode(combo_mat) <- "integer"
  detected <- rowSums(combo_mat) > 0
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
  main = "Omnibus DA overlap: AHC, CLR pseudo-counts, LEfSe"
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
  main = "Omnibus DA overlap: AHC, CLR pseudo-counts, DESeq2"
)

save_upset_plot(
  da_sets,
  set_cols = c("AHC", "CLR_pc0.01", "CLR_pc0.5", "CLR_pc1", "DESeq2"),
  labels = c("AHC", "CLR pc = 0.01", "CLR pc = 0.5", "CLR pc = 1", "DESeq2"),
  filename = "fig4b_da_upset_AHC_CLR_DESeq2",
  title = "AHC + CLR + DESeq2; FDR < 0.05"
)

save_limma_venn(
  da_sets,
  set_cols = c("CLR_pc0.01", "CLR_pc0.5", "CLR_pc1"),
  labels = c("CLR\npc = 0.01", "CLR\npc = 0.5", "CLR\npc = 1"),
  filename = "fig5_da_venn_CLR_pseudocounts",
  main = "Omnibus DA overlap: CLR pseudo-counts"
)
############################################################
## Run summary
############################################################

summary_lines <- c(
  "RealDataOcean AHC analysis completed.",
  paste("Run time:", as.character(Sys.time())),
  paste("Analysis directory:", analysis_dir),
  paste("Output directory:", out_dir),
  paste("Layers used:", paste(layer_levels, collapse = ", ")),
  paste("Layer counts:", paste(names(table(meta$Layer)), as.integer(table(meta$Layer)), sep = "=", collapse = ", ")),
  paste("All selected taxa:", ncol(counts_all)),
  paste("DA filtered taxa:", ncol(counts_da)),
  paste("Pseudo-counts:", paste(pseudo_counts, collapse = ", ")),
  paste("AHC/CLR and DESeq2 DA adjustment:", da_p_adjust_method, "adjusted p-values; threshold =", alpha_level),
  paste("LEfSe:", "lefser package, one-vs-rest layer runs unioned into one LEfSe set; p adjustment =", lefse_p_adjust_method, "; LDA threshold =", lefse_lda_threshold),
  "",
  "Key tables:",
  file.path(tab_dir, "alpha_layer_summary.csv"),
  file.path(tab_dir, "alpha_layer_anova_tests.csv"),
  file.path(tab_dir, "beta_permanova_layer.csv"),
  file.path(tab_dir, "beta_pcoa_procrustes.csv"),
  file.path(tab_dir, "da_omnibus_anova_AHC_CLR.csv"),
  file.path(tab_dir, "da_deseq2_omnibus_lrt.csv"),
  file.path(tab_dir, "da_lefse_one_vs_rest.csv"),
  "",
  "Key figures:",
  file.path(fig_dir, "fig1_alpha_diversity_by_layer.png"),
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
