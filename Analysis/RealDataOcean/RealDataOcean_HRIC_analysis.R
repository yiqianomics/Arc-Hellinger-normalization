############################################################
## RealDataOcean HRIC analysis
##
## Diversity and HRIC coordinates are computed directly from
## the yiqianomics/HRIC package.
############################################################

options(stringsAsFactors = FALSE, width = 140)
set.seed(20260717)

get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

script_path <- get_script_path()
analysis_dir <- if (file.exists(script_path)) dirname(script_path) else script_path
analysis_dir <- normalizePath(analysis_dir, mustWork = TRUE)

out_dir <- file.path(analysis_dir, "real_data_ocean_results")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
if (dir.exists(out_dir)) {
  unlink(list.files(out_dir, full.names = TRUE), recursive = TRUE, force = TRUE)
}
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "phyloseq", "HRIC", "ggplot2", "dplyr", "tidyr", "tibble",
  "purrr", "ggrepel", "maps", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(HRIC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggrepel)
  library(maps)
  library(scales)
})

theme_nature <- function(base_size = 7.8) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.24, color = "grey12"),
      axis.ticks = element_line(linewidth = 0.22, color = "grey12"),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.title = element_text(size = rel(0.95)),
      legend.text = element_text(size = rel(0.88)),
      strip.background = element_rect(fill = "grey96", color = "grey72", linewidth = 0.22),
      strip.text = element_text(face = "bold", color = "black"),
      plot.tag = element_text(face = "bold", size = rel(1.35), color = "black", margin = margin(r = 3, b = 2)),
      plot.tag.position = "topleft",
      plot.margin = margin(7, 5.5, 5.5, 7)
    )
}

save_plot <- function(plot, filename, width, height, dpi = 600) {
  pdf_file <- file.path(fig_dir, paste0(filename, ".pdf"))
  png_file <- file.path(fig_dir, paste0(filename, ".png"))
  tiff_file <- file.path(fig_dir, paste0(filename, ".tiff"))
  ggsave(pdf_file, plot, width = width, height = height, units = "in", limitsize = FALSE)
  ggsave(png_file, plot, width = width, height = height, units = "in", dpi = dpi, limitsize = FALSE)
  ggsave(
    tiff_file, plot, width = width, height = height, units = "in",
    dpi = dpi, compression = "lzw", limitsize = FALSE
  )
  invisible(c(pdf_file, png_file, tiff_file))
}

format_p <- function(p) {
  ifelse(
    is.na(p), "p = NA",
    ifelse(p < 1e-4, "p < 1e-4", paste0("p = ", formatC(p, format = "fg", digits = 2)))
  )
}

format_p_plain <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 1e-4, "<1e-4", formatC(p, format = "fg", digits = 2))
  )
}

sample_breaks <- function(n) {
  if (n <= 45) {
    seq_len(n)
  } else if (n <= 90) {
    unique(c(1, seq(10, n, by = 10), n))
  } else {
    unique(c(1, seq(20, n, by = 20), n))
  }
}

safe_id <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

as_numeric_safe <- function(x) suppressWarnings(as.numeric(as.character(x)))

short_region <- function(x) {
  y <- as.character(x)
  y <- sub("^\\[[^]]+\\][[:space:]]*", "", y)
  y <- sub("\\s*\\(MRGID:.*$", "", y)
  trimws(y)
}

clean_taxon_value <- function(x) {
  y <- trimws(as.character(x))
  y[y == "" | is.na(y)] <- NA_character_
  generic <- tolower(y) %in% c("uncultured", "unclassified", "unknown", "metagenome")
  y[generic] <- NA_character_
  y
}

make_taxon_labels <- function(tax_df) {
  ranks <- c("Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
  rank_label <- c(
    Genus = "genus", Family = "family", Order = "order",
    Class = "class", Phylum = "phylum", Kingdom = "kingdom"
  )
  out <- rep(NA_character_, nrow(tax_df))
  for (rank in ranks) {
    if (!rank %in% names(tax_df)) next
    values <- clean_taxon_value(tax_df[[rank]])
    fill <- is.na(out) & !is.na(values)
    out[fill] <- paste0(values[fill], " (", rank_label[[rank]], ")")
  }
  out[is.na(out)] <- rownames(tax_df)[is.na(out)]
  out
}

ocean_file <- file.path(analysis_dir, "Ocean.RData")
if (!file.exists(ocean_file)) {
  stop("Cannot find Ocean.RData in ", analysis_dir)
}
load(ocean_file)
if (!exists("ps") || !inherits(ps, "phyloseq")) {
  stop("Ocean.RData must contain a phyloseq object named ps.")
}

meta <- as(phyloseq::sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
meta$RegionShort <- short_region(meta$Ocean.region)
meta$Latitude <- as_numeric_safe(meta$Latitude)
meta$Longitude <- as_numeric_safe(meta$Longitude)
meta$Depth.nominal <- as_numeric_safe(meta$Depth.nominal)
meta$Temperature <- as_numeric_safe(meta$Temperature)
meta$Salinity <- as_numeric_safe(meta$Salinity)
meta$ChlorophyllA <- as_numeric_safe(meta$ChlorophyllA)
meta$Layer <- factor(as.character(meta$Layer), levels = c("SRF", "DCM", "MES", "MIX"))

counts_feature <- as(phyloseq::otu_table(ps), "matrix")
if (phyloseq::taxa_are_rows(ps)) {
  counts_feature <- t(counts_feature)
}
storage.mode(counts_feature) <- "numeric"
counts_feature <- counts_feature[meta$SampleID, , drop = FALSE]

keep_samples <- rowSums(counts_feature, na.rm = TRUE) > 0
counts_feature <- counts_feature[keep_samples, , drop = FALSE]
meta <- meta[match(rownames(counts_feature), meta$SampleID), , drop = FALSE]
keep_taxa <- colSums(counts_feature, na.rm = TRUE) > 0
counts_feature <- counts_feature[, keep_taxa, drop = FALSE]

tax_df <- as.data.frame(as(phyloseq::tax_table(ps), "matrix"), stringsAsFactors = FALSE)
tax_df <- tax_df[colnames(counts_feature), , drop = FALSE]
tax_df$FeatureID <- rownames(tax_df)
tax_df$TaxonLabel <- make_taxon_labels(tax_df)

counts_taxon <- rowsum(t(counts_feature), group = tax_df$TaxonLabel, reorder = FALSE)
counts_taxon <- t(counts_taxon)
counts_taxon <- counts_taxon[, colSums(counts_taxon) > 0, drop = FALSE]

brunt_col <- grep("^Brunt", names(meta), value = TRUE)
if (length(brunt_col) != 1) {
  stop("Cannot identify the Brunt-Vaisala covariate column.")
}

covariate_cols <- c(
  "Latitude", "Longitude",
  "Depth.nominal", "Temperature", "Gradient.Surface.temp.SST.",
  "Salinity", "Density", "Oxygen",
  "ChlorophyllA", "Fluorescence",
  "PAR.PC",
  "NO3", "NO2", "NO2NO3", "PO4", "Si",
  "Ammonium.5m", "Iron.5m", "Nitracline",
  "Carbon.total", "CO3", "HCO3", "Alkalinity.total",
  "Depth.Mixed.Layer", "Depth.Min.O2", "Depth.Max.O2",
  brunt_col, "Lyapunov", "Okubo.Weiss", "Residence.time"
)
covariate_labels <- c(
  Latitude = "Latitude (deg)",
  Longitude = "Longitude (deg)",
  Depth.nominal = "Depth (m)",
  Temperature = "Temperature (deg C)",
  Gradient.Surface.temp.SST. = "Surface temp. gradient",
  Salinity = "Salinity (PSU)",
  Density = "Density",
  Oxygen = "Oxygen",
  ChlorophyllA = "Chlorophyll a",
  Fluorescence = "Fluorescence",
  PAR.PC = "PAR",
  NO3 = "Nitrate",
  NO2 = "Nitrite",
  NO2NO3 = "Nitrite + nitrate",
  PO4 = "Phosphate",
  Si = "Silicate",
  Ammonium.5m = "Ammonium at 5 m",
  Iron.5m = "Iron at 5 m",
  Nitracline = "Nitracline",
  Carbon.total = "Total carbon",
  CO3 = "Carbonate",
  HCO3 = "Bicarbonate",
  Alkalinity.total = "Total alkalinity",
  Depth.Mixed.Layer = "Mixed-layer depth",
  Depth.Min.O2 = "Minimum-O2 depth",
  Depth.Max.O2 = "Maximum-O2 depth",
  Lyapunov = "Lyapunov exponent",
  Okubo.Weiss = "Okubo-Weiss",
  Residence.time = "Residence time"
)
covariate_labels[brunt_col] <- "Brunt-Vaisala frequency"
missing_covariates <- setdiff(covariate_cols, names(meta))
if (length(missing_covariates) > 0) {
  stop("Cannot find requested covariate(s): ", paste(missing_covariates, collapse = ", "))
}
covariate_labels <- covariate_labels[covariate_cols]
meta[covariate_cols] <- lapply(meta[covariate_cols], as_numeric_safe)

region_counts <- sort(table(meta$Ocean.region), decreasing = TRUE)
region_lookup <- setNames(names(region_counts), short_region(names(region_counts)))
focus_region_labels <- c("Arctic Ocean", "North Atlantic Ocean")
missing_focus_regions <- setdiff(focus_region_labels, names(region_lookup))
if (length(missing_focus_regions) > 0) {
  stop("Cannot find requested focus region(s): ", paste(missing_focus_regions, collapse = ", "))
}
selected_region_values <- unname(region_lookup[focus_region_labels])
selected_region_labels <- focus_region_labels
selected_region_palette <- setNames(
  c("#2B6EA6", "#16877A")[seq_along(selected_region_values)],
  selected_region_labels
)

panel_defs <- tibble::tibble(
  panel = letters[seq_along(selected_region_values)],
  panel_label = selected_region_labels,
  region_value = selected_region_values,
  file_id = safe_id(selected_region_labels)
)

write.csv(
  tibble::tibble(
    region = names(region_counts),
    region_short = short_region(names(region_counts)),
    n_samples = as.integer(region_counts),
    selected_for_regional_figures = names(region_counts) %in% selected_region_values
  ),
  file.path(tab_dir, "selected_region_counts.csv"),
  row.names = FALSE
)

sparse_covariates <- intersect(
  c("Carbon.total", "CO3", "HCO3", "Alkalinity.total", "PAR.PC"),
  names(meta)
)
availability_cols <- unique(c(covariate_cols, sparse_covariates))
availability_labels <- c(
  covariate_labels,
  Carbon.total = "Total carbon",
  CO3 = "Carbonate",
  HCO3 = "Bicarbonate",
  Alkalinity.total = "Total alkalinity",
  PAR.PC = "PAR"
)
covariate_availability <- purrr::map_dfr(seq_along(selected_region_values), function(i) {
  panel_meta <- meta[meta$Ocean.region == selected_region_values[i], , drop = FALSE]
  purrr::map_dfr(availability_cols, function(covar) {
    values <- as_numeric_safe(panel_meta[[covar]])
    complete <- is.finite(values)
    tibble::tibble(
      region = selected_region_labels[i],
      covariate = covar,
      covariate_label = unname(availability_labels[covar]),
      selected_for_figure3 = covar %in% covariate_cols,
      n_total = nrow(panel_meta),
      n_complete = sum(complete),
      n_unique = length(unique(values[complete])),
      complete_fraction = ifelse(nrow(panel_meta) == 0, NA_real_, sum(complete) / nrow(panel_meta))
    )
  })
})
write.csv(covariate_availability, file.path(tab_dir, "figure3_covariate_availability.csv"), row.names = FALSE)

order_samples <- function(meta_panel, panel_region_value = NA_character_) {
  if (is.na(panel_region_value)) {
    meta_panel %>%
      arrange(RegionShort, Longitude, Latitude, Depth.nominal, SampleID)
  } else {
    meta_panel %>%
      arrange(Longitude, Latitude, Depth.nominal, SampleID)
  }
}

get_panel_meta <- function(region_value) {
  if (is.na(region_value)) {
    panel_meta <- meta
  } else {
    panel_meta <- meta[meta$Ocean.region == region_value, , drop = FALSE]
  }
  order_samples(panel_meta, region_value)
}

compute_diversity <- function(panel_meta, panel, panel_label) {
  sample_ids <- panel_meta$SampleID
  x <- counts_feature[sample_ids, , drop = FALSE]
  x <- x[, colSums(x) > 0, drop = FALSE]
  alpha <- HRIC::SHalpha(x)
  gamma <- as.numeric(HRIC::SHgamma(x))
  beta <- as.numeric(HRIC::SHbeta(x))
  diversity_df <- tibble::tibble(
    Panel = panel,
    PanelLabel = panel_label,
    SampleID = sample_ids,
    Alpha = as.numeric(alpha[sample_ids]),
    Gamma = gamma,
    Beta = beta,
    Turnover0 = gamma - as.numeric(alpha[sample_ids])
  ) %>%
    left_join(
      panel_meta %>%
        select(SampleID, RegionShort, Ocean.region, Layer, all_of(covariate_cols)),
      by = "SampleID"
    )

  diversity_df %>%
    arrange(Turnover0, Alpha, SampleID) %>%
    mutate(
      SamplePlotID = paste0("sample", row_number()),
      SampleIndex = row_number()
    ) %>%
    select(Panel, PanelLabel, SampleID, SamplePlotID, SampleIndex, everything())
}

all_diversity <- purrr::pmap_dfr(
  panel_defs,
  function(panel, panel_label, region_value, file_id) {
    compute_diversity(get_panel_meta(region_value), panel, panel_label)
  }
)

diversity_summary <- all_diversity %>%
  group_by(Panel, PanelLabel) %>%
  summarise(
    n_samples = n(),
    n_features = ncol(counts_feature),
    mean_alpha = mean(Alpha),
    gamma = first(Gamma),
    beta = first(Beta),
    mean_turnover0 = mean(Turnover0),
    .groups = "drop"
  )

write.csv(all_diversity, file.path(tab_dir, "figure2_alpha_gamma_turnover0_values.csv"), row.names = FALSE)
write.csv(diversity_summary, file.path(tab_dir, "diversity_summary_HRIC.csv"), row.names = FALSE)

sample_order_check <- all_diversity %>%
  group_by(Panel, PanelLabel) %>%
  summarise(
    sorted_by_turnover0 = all(diff(Turnover0) >= -1e-12),
    first_turnover0 = first(Turnover0),
    last_turnover0 = last(Turnover0),
    .groups = "drop"
  )

map_data <- ggplot2::map_data("world")
map_points <- meta %>%
  mutate(
    RegionGroup = ifelse(Ocean.region %in% selected_region_values, RegionShort, "Other regions"),
    RegionGroup = factor(RegionGroup, levels = c(selected_region_labels, "Other regions"))
  )

region_centroids <- map_points %>%
  filter(Ocean.region %in% selected_region_values) %>%
  group_by(RegionShort) %>%
  summarise(
    Longitude = mean(Longitude, na.rm = TRUE),
    Latitude = mean(Latitude, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0(RegionShort, "\n", "n = ", n))

map_palette <- c(selected_region_palette, "Other regions" = "#B6B6B6")

figure1 <- ggplot() +
  geom_polygon(
    data = map_data,
    aes(x = long, y = lat, group = group),
    fill = "#EFEFEB", color = "#D0D0CC", linewidth = 0.11
  ) +
  geom_point(
    data = map_points,
    aes(x = Longitude, y = Latitude, fill = RegionGroup),
    shape = 21, color = "grey15", stroke = 0.16, size = 1.85, alpha = 0.92
  ) +
  ggrepel::geom_label_repel(
    data = region_centroids,
    aes(x = Longitude, y = Latitude, label = label),
    size = 2.25, label.size = 0.16, label.padding = unit(0.13, "lines"),
    fill = "#FFFFFF", color = "grey8", min.segment.length = 0,
    segment.color = "grey35", segment.size = 0.18, seed = 20260717
  ) +
  scale_fill_manual(values = map_palette, name = "Sample region") +
  guides(fill = "none") +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-78, 85), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", tag = "a") +
  theme_nature(base_size = 7.9) +
  theme(
    axis.line = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "#FBFDFF", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.tag = element_blank()
  )

save_plot(figure1, "figure1_ocean_sample_map", width = 7.05, height = 3.8)

make_figure2 <- function(div_df, panel_letter, panel_label) {
  n_samples <- nrow(div_df)
  div_df <- div_df %>%
    mutate(
      SampleIndex = seq_len(n_samples),
      SamplePlotID = factor(SamplePlotID, levels = SamplePlotID)
    )
  gamma <- unique(div_df$Gamma)
  x_breaks <- sample_breaks(n_samples)
  x_labels <- paste0("sample", x_breaks)
  axis_size <- ifelse(n_samples > 120, 5.8, ifelse(n_samples > 55, 6.4, 6.9))
  turnover_label_size <- ifelse(n_samples > 120, 1.35, ifelse(n_samples > 55, 1.75, 2.15))
  point_size <- ifelse(n_samples > 120, 1.35, 1.75)

  ggplot(div_df, aes(x = SampleIndex)) +
    geom_segment(
      aes(xend = SampleIndex, y = Alpha, yend = Gamma),
      linewidth = 0.18, color = "#B7B7B7", alpha = 0.78
    ) +
    geom_tile(
      aes(y = -0.065, fill = Turnover0),
      width = 0.92, height = 0.045, alpha = 0.92
    ) +
    geom_text(
      aes(y = -0.065, label = sprintf("%.2f", Turnover0)),
      angle = 90, size = turnover_label_size, color = "grey8",
      vjust = 0.5, hjust = 0.5
    ) +
    geom_hline(yintercept = gamma, linetype = "22", linewidth = 0.42, color = "#9D1F2E") +
    geom_point(
      aes(y = Alpha),
      color = "#2B6EA6", size = point_size, alpha = 0.95
    ) +
    annotate(
      "label", x = max(1, round(n_samples * 0.055)), y = gamma + 0.048,
      label = paste0("gamma = ", sprintf("%.3f", gamma)),
      size = 2.25, fill = "white", color = "#9D1F2E"
    ) +
    scale_fill_gradientn(
      colors = c("#F4F7FB", "#BBD1E8", "#4779B5", "#183A67"),
      limits = range(all_diversity$Turnover0, na.rm = TRUE),
      name = "turnover0"
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      labels = function(x) ifelse(x < 0, "", sprintf("%.1f", x)),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = expansion(mult = c(0.003, 0.01))) +
    coord_cartesian(ylim = c(-0.13, 1), clip = "off") +
    labs(
      x = "Sample",
      y = "Alpha diversity (HRIC evenness)",
      tag = panel_letter
    ) +
    guides(fill = guide_colorbar(order = 1, barwidth = 3.3, barheight = 0.28, title.position = "top")) +
    theme_nature(base_size = 7.9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = axis_size),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.major.y = element_line(linewidth = 0.12, color = "grey90"),
      plot.title = element_blank()
    )
}

figure2_plots <- list()
for (i in seq_len(nrow(panel_defs))) {
  panel <- panel_defs$panel[i]
  panel_label <- panel_defs$panel_label[i]
  file_id <- panel_defs$file_id[i]
  div_df <- all_diversity %>% filter(Panel == panel)
  p <- make_figure2(div_df, panel, panel_label)
  figure2_plots[[panel]] <- p
  width <- ifelse(nrow(div_df) > 120, 14.5, max(6.9, 0.18 * nrow(div_df)))
  save_plot(p, paste0("figure2", panel, "_", file_id, "_alpha_gamma_turnover0"), width = width, height = 4.55)
}

cor_test_safe <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 4 || length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) {
    return(tibble::tibble(n = sum(keep), rho = NA_real_, p = NA_real_))
  }
  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  tibble::tibble(n = sum(keep), rho = unname(test$estimate), p = test$p.value)
}

figure3_stats <- all_diversity %>%
  select(Panel, PanelLabel, SampleID, SamplePlotID, Turnover0, all_of(covariate_cols)) %>%
  tidyr::pivot_longer(all_of(covariate_cols), names_to = "Covariate", values_to = "CovariateValue") %>%
  group_by(Panel, PanelLabel, Covariate) %>%
  group_modify(~ cor_test_safe(.x$CovariateValue, .x$Turnover0)) %>%
  ungroup() %>%
  mutate(
    CovariateLabel = covariate_labels[Covariate],
    label = paste0("rho = ", sprintf("%.2f", rho), "\n", format_p(p), "\n", "n = ", n)
  )

write.csv(figure3_stats, file.path(tab_dir, "figure3_spearman_turnover0_covariates.csv"), row.names = FALSE)

make_figure3 <- function(div_df, stats_df, panel_letter) {
  facet_columns <- 5
  point_size <- ifelse(length(covariate_cols) > 20, 1.15, 1.55)
  label_size <- ifelse(length(covariate_cols) > 20, 1.68, 2.25)
  turnover_limits <- range(div_df$Turnover0, na.rm = TRUE)
  turnover_pad <- diff(turnover_limits) * 0.18
  if (!is.finite(turnover_pad) || turnover_pad == 0) turnover_pad <- 0.05
  long_df <- div_df %>%
    select(SampleID, SamplePlotID, Turnover0, all_of(covariate_cols)) %>%
    tidyr::pivot_longer(all_of(covariate_cols), names_to = "Covariate", values_to = "CovariateValue") %>%
    mutate(CovariateLabel = factor(covariate_labels[Covariate], levels = covariate_labels))
  smooth_df <- long_df %>%
    group_by(Covariate, CovariateLabel) %>%
    filter(
      sum(is.finite(CovariateValue) & is.finite(Turnover0)) >= 4,
      length(unique(CovariateValue[is.finite(CovariateValue) & is.finite(Turnover0)])) >= 2
    ) %>%
    ungroup()
  stats_df <- stats_df %>%
    mutate(CovariateLabel = factor(covariate_labels[Covariate], levels = covariate_labels))

  ggplot(long_df, aes(x = CovariateValue, y = Turnover0)) +
    geom_point(shape = 21, size = point_size, stroke = 0.12, color = "grey12", fill = "#2B6EA6", alpha = 0.78, na.rm = TRUE) +
    geom_smooth(
      data = smooth_df,
      method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.34,
      color = "#9D1F2E", fill = "#9D1F2E", alpha = 0.10, na.rm = TRUE
    ) +
    geom_text(
      data = stats_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE, hjust = -0.08, vjust = 1.18, size = label_size, lineheight = 0.93, color = "grey12"
    ) +
    facet_wrap(~ CovariateLabel, scales = "free_x", ncol = facet_columns) +
    coord_cartesian(ylim = c(turnover_limits[1] - turnover_pad, turnover_limits[2] + turnover_pad)) +
    labs(x = NULL, y = "turnover0 (gamma - alpha)", tag = panel_letter) +
    theme_nature(base_size = 7.9) +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(linewidth = 0.11, color = "grey91")
    )
}

figure3_plots <- list()
for (i in seq_len(nrow(panel_defs))) {
  panel <- panel_defs$panel[i]
  file_id <- panel_defs$file_id[i]
  div_df <- all_diversity %>% filter(Panel == panel)
  stats_df <- figure3_stats %>% filter(Panel == panel)
  p <- make_figure3(div_df, stats_df, panel)
  figure3_plots[[panel]] <- p
  save_plot(
    p,
    paste0("figure3", panel, "_", file_id, "_turnover0_covariates"),
    width = 12.2,
    height = 1.75 * ceiling(length(covariate_cols) / 5) + 0.75
  )
}

taxa_palette <- function(taxa) {
  values <- grDevices::hcl.colors(length(taxa), palette = "Dynamic")
  names(values) <- taxa
  if ("Others" %in% taxa) values["Others"] <- "#D4D4D4"
  values
}

ttest_taxa <- function(counts_taxon_panel, panel_meta, top_taxa, covariates) {
  x <- counts_taxon_panel[, colSums(counts_taxon_panel) > 0, drop = FALSE]
  hric_x <- HRIC::HRIC(x)
  top_taxa <- intersect(top_taxa, colnames(hric_x))
  out <- purrr::map_dfr(top_taxa, function(taxon) {
    purrr::map_dfr(covariates, function(covar) {
      values <- as_numeric_safe(panel_meta[[covar]])
      y <- hric_x[panel_meta$SampleID, taxon]
      keep <- is.finite(values) & is.finite(y)
      if (sum(keep) < 6 || length(unique(values[keep])) < 2) {
        return(tibble::tibble(
          Taxon = taxon, Covariate = covar, n = sum(keep), median_cut = NA_real_,
          mean_low = NA_real_, mean_high = NA_real_, p = NA_real_
        ))
      }
      cut <- median(values[keep], na.rm = TRUE)
      group <- ifelse(values[keep] <= cut, "low", "high")
      if (length(unique(group)) < 2 || min(table(group)) < 3 || length(unique(y[keep])) < 2) {
        return(tibble::tibble(
          Taxon = taxon, Covariate = covar, n = sum(keep), median_cut = cut,
          mean_low = mean(y[keep][group == "low"]), mean_high = mean(y[keep][group == "high"]), p = NA_real_
        ))
      }
      test <- tryCatch(t.test(y[keep] ~ group), error = function(e) NULL)
      tibble::tibble(
        Taxon = taxon,
        Covariate = covar,
        n = sum(keep),
        median_cut = cut,
        mean_low = mean(y[keep][group == "low"]),
        mean_high = mean(y[keep][group == "high"]),
        p = if (is.null(test)) NA_real_ else test$p.value
      )
    })
  })
  out %>%
    mutate(
      q = p.adjust(p, method = "BH"),
      CovariateLabel = covariate_labels[Covariate],
      Significant = !is.na(q) & q < 0.05
    )
}

make_figure4 <- function(panel_meta, panel_letter, panel_label) {
  sample_ids <- panel_meta$SampleID
  counts_panel <- counts_taxon[sample_ids, , drop = FALSE]
  counts_panel <- counts_panel[, colSums(counts_panel) > 0, drop = FALSE]
  rel_panel <- sweep(counts_panel, 1, rowSums(counts_panel), FUN = "/")
  mean_rel <- sort(colMeans(rel_panel, na.rm = TRUE), decreasing = TRUE)
  top_taxa <- names(mean_rel)[seq_len(min(19, length(mean_rel)))]

  tests <- ttest_taxa(counts_panel, panel_meta, top_taxa, covariate_cols)
  significant_taxa <- tests %>%
    filter(Significant) %>%
    pull(Taxon) %>%
    unique()

  stack_counts <- counts_panel[, top_taxa, drop = FALSE]
  others <- rowSums(counts_panel[, setdiff(colnames(counts_panel), top_taxa), drop = FALSE])
  stack_counts <- cbind(stack_counts, Others = others)
  stack_rel <- sweep(stack_counts, 1, rowSums(stack_counts), FUN = "/") * 100

  stack_df <- as.data.frame(stack_rel, check.names = FALSE) %>%
    tibble::rownames_to_column("SampleID") %>%
    tidyr::pivot_longer(-SampleID, names_to = "Taxon", values_to = "RelativeAbundance") %>%
    left_join(panel_meta %>% select(SampleID), by = "SampleID") %>%
    mutate(
      SamplePlotID = paste0("sample", match(SampleID, sample_ids)),
      SamplePlotID = factor(SamplePlotID, levels = paste0("sample", seq_along(sample_ids))),
      Taxon = factor(Taxon, levels = c(top_taxa, "Others"))
    )

  legend_labels <- setNames(levels(stack_df$Taxon), levels(stack_df$Taxon))
  legend_labels[significant_taxa] <- paste0(legend_labels[significant_taxa], "*")
  pal <- taxa_palette(levels(stack_df$Taxon))
  n_samples <- length(sample_ids)
  x_breaks <- paste0("sample", sample_breaks(n_samples))
  axis_size <- ifelse(n_samples > 120, 5.8, ifelse(n_samples > 55, 6.4, 6.9))

  plot <- ggplot(stack_df, aes(x = SamplePlotID, y = RelativeAbundance, fill = Taxon)) +
    geom_col(width = 0.98, linewidth = 0, color = NA) +
    scale_fill_manual(values = pal, labels = legend_labels, drop = FALSE, name = "Taxon (* q < 0.05)") +
    scale_x_discrete(breaks = x_breaks, labels = x_breaks, drop = FALSE) +
    scale_y_continuous(
      breaks = seq(0, 100, 25),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(ylim = c(0, 100), expand = FALSE) +
    labs(x = "Sample", y = "Relative abundance", tag = panel_letter) +
    guides(fill = guide_legend(ncol = 4, byrow = TRUE, override.aes = list(linewidth = 0))) +
    theme_nature(base_size = 7.9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = axis_size),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 6.5),
      legend.text = element_text(size = 6.1),
      legend.key.size = unit(3.3, "mm"),
      panel.grid.major.y = element_line(linewidth = 0.11, color = "grey91")
    )

  list(
    plot = plot,
    top_taxa = tibble::tibble(
      Panel = panel_letter,
      PanelLabel = panel_label,
      Taxon = top_taxa,
      mean_relative_abundance = as.numeric(mean_rel[top_taxa]),
      significant_by_HRIC_ttest = top_taxa %in% significant_taxa
    ),
    tests = tests %>% mutate(Panel = panel_letter, PanelLabel = panel_label, .before = 1),
    stack = stack_df %>% mutate(Panel = panel_letter, PanelLabel = panel_label, .before = 1)
  )
}

figure4_plots <- list()
figure4_top_taxa <- list()
figure4_tests <- list()
figure4_stack_data <- list()
for (i in seq_len(nrow(panel_defs))) {
  panel <- panel_defs$panel[i]
  panel_label <- panel_defs$panel_label[i]
  region_value <- panel_defs$region_value[i]
  file_id <- panel_defs$file_id[i]
  panel_order <- all_diversity %>%
    filter(Panel == panel) %>%
    arrange(SampleIndex) %>%
    pull(SampleID)
  panel_meta <- get_panel_meta(region_value)
  panel_meta <- panel_meta[match(panel_order, panel_meta$SampleID), , drop = FALSE]
  if (!identical(panel_meta$SampleID, panel_order)) {
    stop("Figure 4 sample order does not match Figure 2 for panel ", panel)
  }
  result <- make_figure4(panel_meta, panel, panel_label)
  figure4_plots[[panel]] <- result$plot
  figure4_top_taxa[[panel]] <- result$top_taxa
  figure4_tests[[panel]] <- result$tests
  figure4_stack_data[[panel]] <- result$stack
  width <- ifelse(nrow(panel_meta) > 120, 14.5, max(6.9, 0.18 * nrow(panel_meta)))
  save_plot(result$plot, paste0("figure4", panel, "_", file_id, "_stacked_taxa"), width = width, height = 5.15)
}

write.csv(bind_rows(figure4_top_taxa), file.path(tab_dir, "figure4_top19_taxa_by_panel.csv"), row.names = FALSE)
write.csv(bind_rows(figure4_tests), file.path(tab_dir, "figure4_HRIC_coordinate_ttests_by_panel.csv"), row.names = FALSE)
write.csv(bind_rows(figure4_stack_data), file.path(tab_dir, "figure4_stacked_bar_plot_data.csv"), row.names = FALSE)

figure4_order_check <- bind_rows(figure4_stack_data) %>%
  distinct(Panel, SamplePlotID, SampleID) %>%
  left_join(
    all_diversity %>% select(Panel, SamplePlotID, Figure2SampleID = SampleID),
    by = c("Panel", "SamplePlotID")
  ) %>%
  summarise(matches_figure2_order = all(SampleID == Figure2SampleID), .groups = "drop") %>%
  pull(matches_figure2_order)

combined_figure_files <- file.path(
  fig_dir,
  paste0(
    c(
      "figure2_alpha_gamma_turnover0_all_panels",
      "figure3_turnover0_covariates_all_panels",
      "figure4_stacked_taxa_all_panels"
    ),
    rep(c(".pdf", ".png", ".tiff"), each = 3)
  )
)

sample_panel_key <- all_diversity %>%
  select(
    Panel, PanelLabel, SamplePlotID, SampleIndex, SampleID, RegionShort, Ocean.region,
    Layer, Turnover0, Alpha, Gamma, all_of(covariate_cols)
  )
write.csv(sample_panel_key, file.path(tab_dir, "sample_panel_key.csv"), row.names = FALSE)

qc_checks <- tibble::tibble(
  check = c(
    "HRIC package functions available",
    "Additive diversity partition",
    "turnover0 definition",
    "Figure 2 and Figure 4 sample order",
    "Regional-only figure organization",
    "Figure 3 panel count",
    "Figure 4 taxon grouping"
  ),
  status = c(
    all(c("HRIC", "SHalpha", "SHgamma", "SHbeta") %in% getNamespaceExports("HRIC")),
    all(abs(diversity_summary$gamma - diversity_summary$mean_alpha - diversity_summary$beta) < 1e-10),
    max(abs(all_diversity$Turnover0 - (all_diversity$Gamma - all_diversity$Alpha))) < 1e-12,
    all(sample_order_check$sorted_by_turnover0) && isTRUE(figure4_order_check),
    nrow(panel_defs) == length(selected_region_labels) && !any(file.exists(combined_figure_files)),
    nrow(figure3_stats) == nrow(panel_defs) * length(covariate_cols),
    all(bind_rows(figure4_top_taxa) %>% count(Panel) %>% pull(n) <= 19)
  ),
  detail = c(
    "Required exports available: HRIC, SHalpha, SHgamma, SHbeta",
    paste0("max abs(gamma - mean_alpha - beta) = ", formatC(max(abs(diversity_summary$gamma - diversity_summary$mean_alpha - diversity_summary$beta)), format = "e", digits = 2)),
    paste0("max abs(turnover0 - (gamma - alpha)) = ", formatC(max(abs(all_diversity$Turnover0 - (all_diversity$Gamma - all_diversity$Alpha))), format = "e", digits = 2)),
    "Figure 2 samples are sorted from smallest to largest turnover0, and Figure 4 uses the same sample1...sampleN order.",
    "Figures 2, 3, and 4 are generated only as separate Arctic Ocean and North Atlantic Ocean files; no combined versions are written.",
    paste0(nrow(figure3_stats), " Spearman tests for ", nrow(panel_defs), " panels x ", length(covariate_cols), " covariates"),
    "Each Figure 4 panel uses the top 19 aggregated taxa plus Others in the stacked plot."
  )
)
write.csv(qc_checks, file.path(tab_dir, "analysis_qc_checks.csv"), row.names = FALSE)

hric_description <- utils::packageDescription("HRIC")
hric_github_sha <- hric_description$GithubSHA1
if (is.null(hric_github_sha) || is.na(hric_github_sha) || identical(hric_github_sha, "")) {
  hric_github_sha <- hric_description$RemoteSha
}
if (is.null(hric_github_sha) || is.na(hric_github_sha) || identical(hric_github_sha, "")) {
  hric_github_sha <- "not recorded"
}
covariate_list_text <- paste(unname(covariate_labels), collapse = ", ")

run_summary <- c(
  "RealDataOcean HRIC analysis",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Samples used: ", nrow(meta)),
  paste0("Original feature count used for HRIC diversity: ", ncol(counts_feature)),
  paste0("Aggregated taxon labels used for stacked bars: ", ncol(counts_taxon)),
  paste0("HRIC package version: ", as.character(utils::packageVersion("HRIC"))),
  paste0("HRIC GitHub SHA: ", hric_github_sha),
  paste0("Selected regional panels: ", paste(selected_region_labels, collapse = "; ")),
  paste0("QC checks passed: ", sum(qc_checks$status), "/", nrow(qc_checks)),
  "Figures 2, 3, and 4 are saved only as separate regional files, with no combined all-panel versions.",
  "Figure 2 turnover0 is gamma - alpha, with alpha/gamma/beta from HRIC::SHalpha, HRIC::SHgamma, and HRIC::SHbeta.",
  "Figures do not separate or encode samples by layer.",
  "Figure 2 samples are ordered from smallest to largest turnover0; Figure 4 uses the same sample order for each corresponding panel.",
  paste0("Figure 3 uses Spearman cor.test between turnover0 and these covariates: ", covariate_list_text, "."),
  "Figure 4 asterisks in the legend mark top taxa with at least one BH-adjusted Welch t-test q < 0.05 after median-splitting the Figure 3 covariates and testing HRIC::HRIC coordinates."
)
writeLines(run_summary, file.path(out_dir, "run_summary.txt"))

caption_text <- c(
  "# RealDataOcean figure captions",
  "",
  "**Figure 1. Ocean sample locations.** Points show samples with available latitude and longitude. Highlighted colors denote the Arctic Ocean and North Atlantic Ocean focus regions selected for regional panels; other samples are grey.",
  "",
  "**Figure 2. HRIC alpha diversity, gamma diversity, and sample-level turnover0.** Figure 2a shows Arctic Ocean and Figure 2b shows North Atlantic Ocean. Each panel orders samples from the smallest to the largest turnover0, relabeled as sample1, sample2, and so on. Points are HRIC alpha diversity values from `HRIC::SHalpha`; dashed red lines are panel-specific gamma diversity from `HRIC::SHgamma`; the aligned lower strip labels turnover0, computed as gamma minus alpha for each sample.",
  "",
  "**Figure 3. Correlations between turnover0 and environmental covariates.** Figure 3a shows Arctic Ocean and Figure 3b shows North Atlantic Ocean. Facets show Spearman correlations between turnover0 and the selected location, hydrographic, oxygen, productivity, nutrient, mixed-layer, stability, circulation, and iron covariates. Text in each facet reports Spearman rho, p value from `cor.test(method = \"spearman\")`, and complete-case sample size. Covariate availability is listed in `tables/figure3_covariate_availability.csv`.",
  "",
  "**Figure 4. Relative abundance of the top taxa.** Figure 4a shows Arctic Ocean and Figure 4b shows North Atlantic Ocean. Stacked bars show the top 19 aggregated taxa by mean relative abundance within each panel plus all remaining taxa as Others, using the same sample order as the corresponding Figure 2 panel. Taxa marked with an asterisk in the legend have at least one BH-adjusted Welch t-test q < 0.05 when `HRIC::HRIC` coordinates are compared across median-split Figure 3 covariates."
)
writeLines(caption_text, file.path(out_dir, "figure_captions.md"))

sink(file.path(out_dir, "session_info.txt"))
print(sessionInfo())
sink()

message("Analysis complete. Results written to: ", out_dir)
