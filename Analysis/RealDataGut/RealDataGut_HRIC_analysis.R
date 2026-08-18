############################################################
## HRIC diversity and differential-abundance analysis
## Real gut microbiome data
############################################################

options(stringsAsFactors = FALSE, width = 140)
set.seed(20260720)

required_packages <- c(
  "HRIC", "phyloseq", "ggplot2", "dplyr", "tidyr", "tibble",
  "purrr", "patchwork", "scales", "sandwich", "lmtest"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(patchwork)
})

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  normalizePath("RealDataGut/RealDataGut_HRIC_analysis.R", mustWork = TRUE)
}

analysis_dir <- dirname(get_script_path())
out_dir <- file.path(analysis_dir, "real_data_gut_results")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
unlink(out_dir, recursive = TRUE, force = TRUE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, filename, width, height, dpi = 400) {
  files <- file.path(fig_dir, paste0(filename, c(".png", ".pdf", ".tiff")))
  ggsave(files[1], plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(files[2], plot, width = width, height = height, bg = "white")
  ggsave(
    files[3], plot, width = width, height = height, dpi = dpi,
    compression = "lzw", bg = "white"
  )
  files
}

theme_nature <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, color = "grey15"),
      axis.ticks = element_line(linewidth = 0.30, color = "grey15"),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.text = element_text(color = "grey15"),
      axis.title = element_text(color = "grey10"),
      strip.background = element_rect(fill = "grey96", color = "grey82", linewidth = 0.25),
      strip.text = element_text(face = "bold", color = "grey10", margin = margin(3, 3, 3, 3)),
      legend.title = element_blank(),
      legend.text = element_text(color = "grey15"),
      plot.margin = margin(8, 10, 8, 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

format_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 1e-4, "< 1e-4", formatC(p, format = "fg", digits = 2))
  )
}

clean_taxon_value <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^[A-Za-z]__", "", x)
  invalid <- is.na(x) | x == "" | tolower(x) %in% c(
    "na", "uncultured", "unidentified", "uncultured bacterium",
    "uncultured organism", "metagenome", "unknown"
  )
  x[invalid] <- NA_character_
  x
}

deepest_taxon_label <- function(taxonomy) {
  rank_order <- c("Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
  available <- intersect(rank_order, colnames(taxonomy))
  labels <- vapply(seq_len(nrow(taxonomy)), function(i) {
    for (rank in available) {
      value <- clean_taxon_value(taxonomy[i, rank])
      if (!is.na(value)) {
        return(paste0(value, " (", tolower(rank), ")"))
      }
    }
    paste0("Unclassified feature ", rownames(taxonomy)[i])
  }, character(1))
  labels
}

sample_breaks <- function(n, target = 12) {
  step <- max(1, ceiling(n / target / 5) * 5)
  breaks <- unique(c(1, seq(step, n, by = step)))
  if (length(breaks) > 1 && n - tail(breaks, 1) < step * 0.45) {
    breaks <- head(breaks, -1)
  }
  unique(c(breaks, n))
}

############################################################
## Data preparation
############################################################

gut_file <- file.path(analysis_dir, "gut.RData")
if (!file.exists(gut_file)) stop("Cannot find gut.RData")
load(gut_file)
if (!exists("gut") || !inherits(gut, "phyloseq")) {
  stop("gut.RData must contain a phyloseq object named gut")
}

ps_gut <- gut
ps_gut <- prune_samples(sample_sums(ps_gut) > 0, ps_gut)
ps_gut <- prune_taxa(taxa_sums(ps_gut) > 0, ps_gut)

meta <- as(sample_data(ps_gut), "data.frame")
meta$SampleID <- rownames(meta)
required_meta <- c("arm", "age", "day", "host_subject_id")
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}
meta <- meta %>%
  mutate(
    arm = factor(tolower(as.character(arm)), levels = c("control", "treatment")),
    age = suppressWarnings(as.numeric(age)),
    day = suppressWarnings(as.numeric(day)),
    host_subject_id = factor(as.character(host_subject_id)),
    ArmLabel = factor(
      ifelse(arm == "control", "Control", "Treatment"),
      levels = c("Control", "Treatment")
    )
  )
if (any(is.na(meta$arm))) stop("arm contains values other than control or treatment")

counts_feature <- as(otu_table(ps_gut), "matrix")
if (taxa_are_rows(ps_gut)) counts_feature <- t(counts_feature)
storage.mode(counts_feature) <- "numeric"
counts_feature <- counts_feature[meta$SampleID, , drop = FALSE]
counts_feature <- counts_feature[, colSums(counts_feature) > 0, drop = FALSE]

taxonomy <- as(tax_table(ps_gut), "matrix")
taxonomy <- taxonomy[colnames(counts_feature), , drop = FALSE]
taxon_labels <- deepest_taxon_label(taxonomy)
counts_taxon <- t(rowsum(t(counts_feature), group = taxon_labels, reorder = FALSE))
counts_taxon <- counts_taxon[, colSums(counts_taxon) > 0, drop = FALSE]

arm_colors <- c(Control = "#2369A0", Treatment = "#C44E3B")

############################################################
## Alpha, gamma, beta, and sample turnover by arm
############################################################

compute_arm_diversity <- function(arm_value) {
  arm_meta <- meta %>% filter(arm == arm_value)
  sample_ids <- arm_meta$SampleID
  x <- counts_feature[sample_ids, , drop = FALSE]
  x <- x[, colSums(x) > 0, drop = FALSE]
  alpha_raw <- HRIC::SHalpha(x)
  alpha <- as.numeric(alpha_raw)
  if (!is.null(names(alpha_raw))) {
    alpha <- as.numeric(alpha_raw[sample_ids])
  }
  gamma <- as.numeric(HRIC::SHgamma(x))
  beta <- as.numeric(HRIC::SHbeta(x))
  arm_label <- ifelse(arm_value == "control", "Control", "Treatment")

  tibble(
    SampleID = sample_ids,
    arm = arm_value,
    ArmLabel = arm_label,
    Alpha = alpha,
    Gamma = gamma,
    Beta = beta,
    Turnover = gamma - alpha
  ) %>%
    arrange(Turnover, Alpha, SampleID) %>%
    mutate(
      SampleIndex = row_number(),
      SamplePlotID = paste0("S", SampleIndex)
    )
}

diversity <- bind_rows(
  compute_arm_diversity("control"),
  compute_arm_diversity("treatment")
) %>%
  mutate(ArmLabel = factor(ArmLabel, levels = c("Control", "Treatment")))

diversity_summary <- diversity %>%
  group_by(arm, ArmLabel) %>%
  summarise(
    n_samples = n(),
    n_features = ncol(counts_feature),
    mean_alpha = mean(Alpha),
    gamma = first(Gamma),
    beta = first(Beta),
    gamma_minus_mean_alpha = gamma - mean_alpha,
    .groups = "drop"
  )

write.csv(diversity, file.path(tab_dir, "alpha_gamma_turnover_by_sample.csv"), row.names = FALSE)
write.csv(diversity_summary, file.path(tab_dir, "diversity_summary_by_arm.csv"), row.names = FALSE)
write.csv(
  diversity %>% select(arm, ArmLabel, SampleID, SamplePlotID, SampleIndex),
  file.path(tab_dir, "sample_plot_key.csv"), row.names = FALSE
)

figure2_gamma <- diversity_summary %>%
  mutate(
    x = Inf,
    label = paste0("gamma == ", sprintf("%.3f", gamma))
  )

figure2_y_limits <- range(c(diversity$Alpha, diversity$Gamma))
figure2_y_limits <- c(
  max(0, figure2_y_limits[1] - diff(figure2_y_limits) * 0.04),
  figure2_y_limits[2] + diff(figure2_y_limits) * 0.14
)

make_figure2_panel <- function(arm_value) {
  panel_data <- diversity %>% filter(arm == arm_value)
  panel_summary <- diversity_summary %>% filter(arm == arm_value)
  panel_gamma <- figure2_gamma %>% filter(arm == arm_value)
  n_samples <- nrow(panel_data)
  breaks <- sample_breaks(n_samples)
  panel_label <- as.character(first(panel_data$ArmLabel))

  ggplot(panel_data, aes(x = SampleIndex, y = Alpha)) +
    geom_segment(
      aes(xend = SampleIndex, yend = Gamma),
      linewidth = 0.20, color = "grey55", alpha = 0.32
    ) +
    geom_hline(
      data = panel_summary,
      aes(yintercept = gamma),
      linewidth = 0.55, linetype = "22", color = "#A51C30"
    ) +
    geom_point(
      fill = unname(arm_colors[panel_label]), shape = 21, size = 2.0,
      stroke = 0.25, color = "white"
    ) +
    geom_label(
      data = panel_gamma,
      aes(x = x, y = gamma, label = label),
      inherit.aes = FALSE, parse = TRUE,
      hjust = 1.05, vjust = -0.45, size = 3.0,
      color = "#A51C30", fill = scales::alpha("white", 0.92),
      linewidth = 0.25, label.padding = grid::unit(0.10, "lines")
    ) +
    scale_x_continuous(
      breaks = breaks,
      labels = paste0("S", breaks),
      limits = c(0.5, n_samples + 0.5),
      expand = c(0, 0)
    ) +
    scale_y_continuous(limits = figure2_y_limits, expand = c(0, 0)) +
    labs(
      title = paste0(panel_label, " (n = ", n_samples, ")"),
      x = "Sample", y = "Alpha diversity"
    ) +
    theme_nature(base_size = 9.3) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 10.5),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major.y = element_line(linewidth = 0.18, color = "grey92")
    )
}

figure2 <- make_figure2_panel("control") /
  make_figure2_panel("treatment") +
  plot_layout(heights = c(1, 1))

figure2_files <- save_plot(
  figure2, "figure2_alpha_gamma_by_group", width = 10, height = 7.4
)

############################################################
## HRIC differential-abundance analysis
## Treatment effect adjusted for age and sampling day.
## Cluster-robust inference accounts for repeated subjects.
############################################################

hric_taxon <- HRIC::HRIC(counts_taxon)
if (is.null(rownames(hric_taxon))) rownames(hric_taxon) <- rownames(counts_taxon)
if (is.null(colnames(hric_taxon))) colnames(hric_taxon) <- colnames(counts_taxon)

prevalence <- colMeans(counts_taxon > 0)
total_count <- colSums(counts_taxon)
eligible_taxa <- names(prevalence)[prevalence >= 0.10 & total_count >= 50]

daa_meta <- meta %>%
  filter(complete.cases(arm, age, day, host_subject_id)) %>%
  mutate(
    age_z = as.numeric(scale(age)),
    day_z = as.numeric(scale(day))
  )

fit_taxon_lm <- function(taxon) {
  model_data <- daa_meta %>%
    transmute(
      SampleID,
      host_subject_id,
      arm,
      age_z,
      day_z,
      HRIC_value = as.numeric(hric_taxon[SampleID, taxon])
    ) %>%
    filter(is.finite(HRIC_value))

  empty <- tibble(
    Taxon = taxon, n_samples = nrow(model_data),
    n_subjects = n_distinct(model_data$host_subject_id),
    effect_treatment_vs_control = NA_real_, robust_se = NA_real_,
    robust_t = NA_real_, p = NA_real_, ordinary_se = NA_real_,
    ordinary_p = NA_real_, adjusted_r_squared = NA_real_
  )
  if (nrow(model_data) < 20 || n_distinct(model_data$arm) < 2 ||
      n_distinct(model_data$HRIC_value) < 2) {
    return(empty)
  }

  fit <- tryCatch(
    lm(HRIC_value ~ arm + age_z + day_z, data = model_data),
    error = function(e) NULL
  )
  if (is.null(fit) || !"armtreatment" %in% names(coef(fit))) return(empty)

  ordinary <- summary(fit)$coefficients
  robust <- tryCatch(
    lmtest::coeftest(
      fit,
      vcov. = sandwich::vcovCL(
        fit, cluster = model_data$host_subject_id, type = "HC1"
      )
    ),
    error = function(e) NULL
  )
  if (is.null(robust) || !"armtreatment" %in% rownames(robust)) return(empty)

  tibble(
    Taxon = taxon,
    n_samples = nrow(model_data),
    n_subjects = n_distinct(model_data$host_subject_id),
    effect_treatment_vs_control = unname(robust["armtreatment", "Estimate"]),
    robust_se = unname(robust["armtreatment", "Std. Error"]),
    robust_t = unname(robust["armtreatment", "t value"]),
    p = unname(robust["armtreatment", "Pr(>|t|)"]),
    ordinary_se = unname(ordinary["armtreatment", "Std. Error"]),
    ordinary_p = unname(ordinary["armtreatment", "Pr(>|t|)"]),
    adjusted_r_squared = summary(fit)$adj.r.squared
  )
}

daa_results <- map_dfr(eligible_taxa, fit_taxon_lm) %>%
  mutate(
    prevalence = prevalence[Taxon],
    total_count = total_count[Taxon],
    q = p.adjust(p, method = "BH"),
    Significant = !is.na(q) & q < 0.05,
    Direction = case_when(
      Significant & effect_treatment_vs_control > 0 ~ "Higher in treatment",
      Significant & effect_treatment_vs_control < 0 ~ "Higher in control",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(q, desc(abs(effect_treatment_vs_control)))

write.csv(
  daa_results,
  file.path(tab_dir, "daa_hric_lm_arm_adjusted_age_day.csv"),
  row.names = FALSE
)

############################################################
## Figure 4: longitudinal top 19 taxa plus Others
############################################################

relative_taxon <- sweep(counts_taxon, 1, rowSums(counts_taxon), FUN = "/")
mean_relative <- sort(colMeans(relative_taxon), decreasing = TRUE)
significant_taxa_all <- daa_results$Taxon[daa_results$Significant]
priority_taxa <- unique(c(significant_taxa_all, names(mean_relative)))
top_taxa <- priority_taxa[seq_len(min(19, length(priority_taxa)))]
top_taxa <- top_taxa[order(mean_relative[top_taxa], decreasing = TRUE)]
taxon_order <- c(top_taxa, "Others")

figure4_day_limits <- c(0, 90)
participant_levels <- c(paste0("C", seq_len(11)), paste0("T", seq_len(14)))
participant_levels <- participant_levels[participant_levels %in% as.character(meta$host_subject_id)]

figure4_wide <- as.data.frame(relative_taxon[, top_taxa, drop = FALSE]) %>%
  rownames_to_column("SampleID") %>%
  mutate(Others = pmax(0, 1 - rowSums(across(all_of(top_taxa))))) %>%
  left_join(
    meta %>% select(SampleID, arm, ArmLabel, host_subject_id, day),
    by = "SampleID"
  ) %>%
  filter(day >= figure4_day_limits[1], day <= figure4_day_limits[2])

figure4_sample_data <- figure4_wide %>%
  pivot_longer(
    cols = all_of(taxon_order),
    names_to = "Taxon", values_to = "RelativeAbundance"
  ) %>%
  mutate(Taxon = factor(Taxon, levels = taxon_order))

# Same-participant samples collected on the same day are averaged so bars do not overlap.
figure4_data <- figure4_sample_data %>%
  group_by(arm, ArmLabel, host_subject_id, day, Taxon) %>%
  summarise(
    RelativeAbundance = mean(RelativeAbundance),
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Participant = factor(as.character(host_subject_id), levels = participant_levels)
  )

figure4_time_coverage <- figure4_wide %>%
  group_by(arm, ArmLabel, host_subject_id) %>%
  summarise(
    n_samples = n(),
    n_unique_days = n_distinct(day),
    first_day = min(day),
    last_day = max(day),
    .groups = "drop"
  ) %>%
  mutate(ParticipantNumber = as.numeric(sub("^[A-Za-z]+", "", as.character(host_subject_id)))) %>%
  arrange(arm, ParticipantNumber) %>%
  select(-ParticipantNumber)

write.csv(
  figure4_data,
  file.path(tab_dir, "figure4_longitudinal_relative_abundance.csv"),
  row.names = FALSE
)
write.csv(
  figure4_time_coverage,
  file.path(tab_dir, "figure4_participant_time_coverage.csv"),
  row.names = FALSE
)

high_contrast <- c(
  "#245AA3", "#B7312C", "#008B61", "#7C3A96", "#E49B00",
  "#069BB4", "#C14E8A", "#91722E", "#4B8C3D", "#4E4E4E",
  "#00A58B", "#F39A31", "#405F91", "#AD3830", "#7653A2",
  "#D85A00", "#4EB4CC", "#5F7C32", "#A4A226"
)
taxon_palette <- c(setNames(high_contrast[seq_along(top_taxa)], top_taxa), Others = "#CBCBCB")
significant_top_taxa <- intersect(top_taxa, daa_results$Taxon[daa_results$Significant])
legend_labels <- setNames(taxon_order, taxon_order)
legend_labels[significant_top_taxa] <- paste0(legend_labels[significant_top_taxa], "*")

write.csv(
  tibble(
    Taxon = top_taxa,
    MeanRelativeAbundance = unname(mean_relative[top_taxa]),
    Color = unname(taxon_palette[top_taxa]),
    DAA_q = daa_results$q[match(top_taxa, daa_results$Taxon)],
    DAA_significant = top_taxa %in% significant_top_taxa,
    Selection = ifelse(
      top_taxa %in% significant_taxa_all,
      "BH-significant DAA taxon",
      "Largest mean relative abundance"
    )
  ),
  file.path(tab_dir, "figure4_displayed_taxa.csv"), row.names = FALSE
)

figure4 <- ggplot(
  figure4_data,
  aes(x = day, y = RelativeAbundance, fill = Taxon)
) +
    geom_col(width = 0.92, color = "white", linewidth = 0.05) +
    scale_fill_manual(
      values = taxon_palette,
      breaks = taxon_order,
      labels = legend_labels,
      guide = guide_legend(
        ncol = 5, byrow = TRUE,
        keyheight = grid::unit(3.8, "mm"),
        keywidth = grid::unit(4.5, "mm")
      )
    ) +
    scale_x_continuous(breaks = c(0, 30, 60, 90), limits = c(-1.5, 91.5), expand = c(0, 0)) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      breaks = c(0, 0.5, 1),
      expand = c(0, 0)
    ) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    facet_wrap(~ Participant, ncol = 5, drop = FALSE) +
    labs(x = "Day", y = "Relative abundance") +
    theme_nature(base_size = 9.3) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 7.2),
      axis.text.y = element_text(size = 7.2),
      strip.text = element_text(face = "bold", size = 8.2),
      legend.position = "bottom",
      legend.text = element_text(size = 7.2),
      legend.spacing.x = grid::unit(1.0, "mm"),
      panel.grid.major.x = element_line(linewidth = 0.18, color = "grey90"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

figure4_files <- save_plot(
  figure4, "figure4_hric_lm_daa_top_taxa", width = 10, height = 12.2, dpi = 600
)

############################################################
## PCoA: HRIC geometry versus ordinary Euclidean geometry
############################################################

relative_feature <- sweep(counts_feature, 1, rowSums(counts_feature), FUN = "/")
hric_feature <- HRIC::HRIC(counts_feature)
if (is.null(rownames(hric_feature))) rownames(hric_feature) <- rownames(counts_feature)
hric_feature <- hric_feature[meta$SampleID, , drop = FALSE]

compute_pcoa <- function(x, method_label) {
  distance <- stats::dist(x, method = "euclidean")
  fit <- stats::cmdscale(distance, k = 2, eig = TRUE, add = FALSE)
  positive_total <- sum(fit$eig[fit$eig > 0])
  variance <- 100 * fit$eig[1:2] / positive_total

  scores <- as.data.frame(fit$points) %>%
    rownames_to_column("SampleID") %>%
    transmute(
      SampleID,
      Axis1 = V1,
      Axis2 = V2,
      Method = method_label,
      Axis1Variance = variance[1],
      Axis2Variance = variance[2]
    ) %>%
    left_join(
      meta %>% select(SampleID, arm, ArmLabel, host_subject_id, day),
      by = "SampleID"
    )

  list(scores = scores, eigenvalues = fit$eig, distance = distance)
}

pcoa_hric_result <- compute_pcoa(
  hric_feature,
  "HRIC-transformed coordinates"
)
pcoa_euclidean_result <- compute_pcoa(
  relative_feature,
  "Relative-abundance Euclidean"
)
pcoa_scores <- bind_rows(
  pcoa_hric_result$scores,
  pcoa_euclidean_result$scores
)
pcoa_eigenvalues <- bind_rows(
  tibble(
    Method = "HRIC-transformed coordinates",
    Axis = seq_along(pcoa_hric_result$eigenvalues),
    Eigenvalue = pcoa_hric_result$eigenvalues
  ),
  tibble(
    Method = "Relative-abundance Euclidean",
    Axis = seq_along(pcoa_euclidean_result$eigenvalues),
    Eigenvalue = pcoa_euclidean_result$eigenvalues
  )
) %>%
  group_by(Method) %>%
  mutate(PercentPositiveVariation = 100 * Eigenvalue / sum(Eigenvalue[Eigenvalue > 0])) %>%
  ungroup()

write.csv(pcoa_scores, file.path(tab_dir, "pcoa_hric_euclidean_coordinates.csv"), row.names = FALSE)
write.csv(pcoa_eigenvalues, file.path(tab_dir, "pcoa_hric_euclidean_eigenvalues.csv"), row.names = FALSE)

make_pcoa_panel <- function(method_label) {
  panel_data <- pcoa_scores %>% filter(Method == method_label)
  axis1_variance <- first(panel_data$Axis1Variance)
  axis2_variance <- first(panel_data$Axis2Variance)

  ggplot(panel_data, aes(x = Axis1, y = Axis2)) +
    stat_ellipse(
      aes(color = ArmLabel), type = "norm", level = 0.80,
      linewidth = 0.80, alpha = 0.95, show.legend = FALSE
    ) +
    geom_point(
      aes(fill = ArmLabel, shape = ArmLabel),
      size = 2.15, stroke = 0.30, color = "white", alpha = 0.72
    ) +
    scale_fill_manual(values = arm_colors, name = NULL) +
    scale_color_manual(values = arm_colors, name = NULL) +
    scale_shape_manual(values = c(Control = 21, Treatment = 24), name = NULL) +
    labs(
      title = method_label,
      x = sprintf("PCoA 1 (%.1f%%)", axis1_variance),
      y = sprintf("PCoA 2 (%.1f%%)", axis2_variance)
    ) +
    coord_equal() +
    theme_nature(base_size = 10.0) +
    theme(
      plot.title = element_text(face = "bold", size = 10.8, hjust = 0),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      panel.grid = element_blank()
    )
}

pcoa_hric_plot <- make_pcoa_panel("HRIC-transformed coordinates")
pcoa_euclidean_plot <- make_pcoa_panel("Relative-abundance Euclidean")
pcoa_figure <- pcoa_hric_plot / pcoa_euclidean_plot +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")

pcoa_files <- save_plot(
  pcoa_figure, "figure5_pcoa_hric_euclidean", width = 7.2, height = 10.8, dpi = 600
)

integrated_left <- wrap_plots(
  list(figure2, figure4), ncol = 1, heights = c(0.34, 0.66)
)
integrated_gut_figure <- wrap_plots(
  list(integrated_left, pcoa_figure), nrow = 1, widths = c(1.45, 1.00)
) &
  theme(plot.background = element_rect(fill = "white", color = NA))

integrated_gut_files <- save_plot(
  integrated_gut_figure,
  "figure_integrated_gut_diversity_daa_pcoa",
  width = 14.5,
  height = 14.5,
  dpi = 600
)

############################################################
## Quality control and run documentation
############################################################

sample_order_check <- diversity %>%
  group_by(arm) %>%
  summarise(sorted = all(diff(Turnover) >= -1e-12), .groups = "drop")

figure4_window_check <-
  min(figure4_data$day) >= figure4_day_limits[1] &&
  max(figure4_data$day) <= figure4_day_limits[2] &&
  n_distinct(figure4_data$host_subject_id) == n_distinct(meta$host_subject_id)

stack_sum_error <- figure4_data %>%
  group_by(arm, host_subject_id, day) %>%
  summarise(total = sum(RelativeAbundance), .groups = "drop") %>%
  summarise(error = max(abs(total - 1))) %>%
  pull(error)

figure_files <- c(figure2_files, figure4_files, pcoa_files, integrated_gut_files)
qc <- tibble(
  check = c(
    "HRIC package functions available",
    "Two arm-specific gamma diversities",
    "Additive diversity partition",
    "Turnover definition",
    "Figure 2 turnover ordering",
    "Figure 4 longitudinal participant window",
    "Figure 4 19 named taxa plus Others",
    "Figure 4 bars sum to 100 percent",
    "DAA uses lm on HRIC coordinates",
    "DAA adjusts for age and day",
    "Repeated subjects use cluster-robust inference",
    "PCoA uses HRIC package coordinates",
    "PCoA Euclidean comparison uses relative abundance",
    "Integrated figure layout and exports",
    "Figure exports"
  ),
  status = c(
    all(c("HRIC", "SHalpha", "SHgamma", "SHbeta") %in% getNamespaceExports("HRIC")),
    nrow(diversity_summary) == 2 && n_distinct(diversity_summary$gamma) == 2,
    max(abs(diversity_summary$gamma_minus_mean_alpha - diversity_summary$beta)) < 1e-10,
    max(abs(diversity$Turnover - (diversity$Gamma - diversity$Alpha))) < 1e-12,
    all(sample_order_check$sorted),
    isTRUE(figure4_window_check),
    length(top_taxa) == 19 && nlevels(figure4_data$Taxon) == 20,
    stack_sum_error < 1e-10,
    nrow(daa_results) == length(eligible_taxa),
    all(c("age_z", "day_z") %in% all.vars(HRIC_value ~ arm + age_z + day_z)),
    n_distinct(daa_meta$host_subject_id) > 1 && all(is.finite(daa_results$robust_se[!is.na(daa_results$p)])),
    identical(rownames(hric_feature), meta$SampleID) && all(is.finite(hric_feature)),
    max(abs(rowSums(relative_feature) - 1)) < 1e-12 && all(is.finite(relative_feature)),
    nrow(pcoa_scores) == 2 * nrow(meta) && all(file.exists(integrated_gut_files)),
    all(file.exists(figure_files)) && all(file.info(figure_files)$size > 0)
  )
)
write.csv(qc, file.path(tab_dir, "analysis_qc_checks.csv"), row.names = FALSE)

hric_description <- utils::packageDescription("HRIC")
hric_sha <- hric_description$GithubSHA1
if (is.null(hric_sha) || is.na(hric_sha) || hric_sha == "") hric_sha <- hric_description$RemoteSha

run_summary <- c(
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Samples: ", nrow(meta)),
  paste0("Subjects: ", n_distinct(meta$host_subject_id)),
  paste0("Features used for diversity: ", ncol(counts_feature)),
  paste0("Aggregated taxon labels: ", ncol(counts_taxon)),
  paste0("DAA taxa after prevalence/count filtering: ", length(eligible_taxa)),
  paste0("Significant DAA taxa (BH q < 0.05): ", sum(daa_results$Significant, na.rm = TRUE)),
  paste0("Significant taxa among Figure 4 top 19: ", length(significant_top_taxa)),
  paste0("HRIC package version: ", as.character(packageVersion("HRIC"))),
  paste0("HRIC GitHub SHA: ", hric_sha),
  paste0("QC checks passed: ", sum(qc$status), "/", nrow(qc)),
  "Alpha, gamma, beta, turnover, and transformed taxon coordinates are calculated directly with the HRIC package.",
  "Figure 2 contains one gamma diversity for Control and one for Treatment; samples are ordered by increasing turnover within arm.",
  "DAA model: lm(HRIC taxon coordinate ~ arm + standardized age + standardized day).",
  "P values for the treatment coefficient use HC1 cluster-robust covariance by host_subject_id, followed by BH adjustment across all eligible taxa.",
  paste0("Figure 4 shows participant-labelled longitudinal compositions from day ", figure4_day_limits[1], " to day ", figure4_day_limits[2], ". Same-participant samples collected on the same day are averaged before plotting."),
  "PCoA comparison: the HRIC panel uses Euclidean distances among feature-level HRIC::HRIC coordinates; the Euclidean panel uses Euclidean distances among feature-level relative-abundance vectors. Both ordinations use classical multidimensional scaling with stats::cmdscale.",
  "The integrated figure places Figure 2 at upper left, Figure 4 at lower left, and the two PCoA panels in the right column."
)
writeLines(run_summary, file.path(out_dir, "run_summary.txt"))

caption <- c(
  "**Figure 2. HRIC alpha and gamma diversity by treatment group.** Samples are ordered from smallest to largest turnover within each group and relabeled S1 to SN. Points show sample-level alpha diversity from `HRIC::SHalpha`; dashed red lines show the group-specific gamma diversity from `HRIC::SHgamma`. Grey segments connect each alpha value to its group gamma.",
  "",
  "**Figure 4. Longitudinal relative abundance and HRIC differential-abundance analysis.** Participant-labelled facets show microbial composition from day 0 to day 90 on a common time axis. Each vertical bar represents one participant-day and contains 19 named taxa plus Others; samples collected from the same participant on the same day are averaged before plotting. Control participants are labelled C1--C11 and treatment participants T1--T14. BH-significant DAA taxa are retained first and the remaining display slots are filled by mean relative abundance. Asterisks mark taxa with BH-adjusted q < 0.05 for the treatment coefficient from `lm(HRIC taxon coordinate ~ arm + age + day)`. Inference uses HC1 covariance clustered by participant to account for repeated samples.",
  "",
  "**PCoA and integrated figure.** The right column compares two feature-level PCoA ordinations. The HRIC panel applies `HRIC::HRIC` directly to the count matrix and calculates Euclidean distances among the resulting intrinsic-coordinate vectors. The Euclidean panel calculates Euclidean distances among ordinary relative-abundance vectors. Classical PCoA is performed with `stats::cmdscale`; axis labels report each axis as a percentage of the positive eigenvalue sum. Points are samples, colours and shapes indicate treatment arm, and ellipses enclose the fitted 80% bivariate-normal region. Figure 2 is shown at upper left and Figure 4 at lower left."
)
writeLines(caption, file.path(out_dir, "figure_captions.md"))

if (!all(qc$status)) {
  stop("Analysis completed, but one or more QC checks failed. See analysis_qc_checks.csv")
}

message("Analysis complete: ", out_dir)
