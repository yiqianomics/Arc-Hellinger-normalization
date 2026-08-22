# Ocean microbiome analysis for the Arctic and North Atlantic regions.

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

out_dir <- file.path(analysis_dir, "results")
fig_dir <- out_dir
tab_dir <- file.path(out_dir, "tables")
unlink(out_dir, recursive = TRUE, force = TRUE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "phyloseq", "HRIC", "ggplot2", "dplyr", "tidyr", "tibble",
  "purrr", "ggrepel", "maps", "scales", "patchwork"
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
  library(patchwork)
})

theme_nature <- function(base_size = 9.0) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.34, color = "grey10"),
      axis.ticks = element_line(linewidth = 0.30, color = "grey10"),
      axis.title = element_text(color = "grey8", size = rel(1.00), face = "bold"),
      axis.text = element_text(color = "grey10", size = rel(0.90)),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.title = element_text(size = rel(0.90), face = "bold"),
      legend.text = element_text(size = rel(0.82)),
      strip.background = element_rect(fill = "#F7F7F7", color = "#D6D6D6", linewidth = 0.16),
      strip.text = element_text(face = "bold", color = "grey8", size = rel(0.82), margin = margin(1.5, 2, 1.5, 2)),
      plot.tag = element_text(face = "bold", size = rel(1.25), color = "grey8", margin = margin(r = 3, b = 2)),
      plot.tag.position = "topleft",
      plot.margin = margin(4.5, 4.5, 4.5, 5.5)
    )
}

save_plot <- function(plot, filename, width, height) {
  pdf_file <- file.path(fig_dir, paste0(filename, ".pdf"))
  ggsave(
    pdf_file,
    plot,
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
  invisible(pdf_file)
}

format_p_plotmath <- function(p) {
  ifelse(
    is.na(p), 'italic(p) == "NA"',
    ifelse(p < 1e-4, "italic(p) < 1e-4", paste0("italic(p) == ", formatC(p, format = "fg", digits = 2)))
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
  "Depth.nominal", "Temperature", "Oxygen",
  "ChlorophyllA", "PO4", "NO3"
)
covariate_labels <- c(
  Latitude = "Latitude (deg)",
  Longitude = "Longitude (deg)",
  Depth.nominal = "Depth (m)",
  Temperature = "Temperature (°C)",
  Gradient.Surface.temp.SST. = "Surface temp. gradient",
  Salinity = "Salinity (PSU)",
  Density = "Density",
  Oxygen = "Oxygen (µmol kg-1)",
  ChlorophyllA = "Chlorophyll a (mg m-3)",
  Fluorescence = "Fluorescence",
  PAR.PC = "PAR",
  NO3 = "Nitrate (µmol l-1)",
  NO2 = "Nitrite",
  NO2NO3 = "Nitrite + nitrate",
  PO4 = "Phosphate (µmol l-1)",
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
selected_region_ocean_fill <- setNames(
  c("#DCECF7", "#DFF2EC")[seq_along(selected_region_values)],
  selected_region_labels
)
selected_region_mark <- setNames(
  c("#1A5D91", "#047C70")[seq_along(selected_region_values)],
  selected_region_labels
)

panel_defs <- tibble::tibble(
  panel = letters[seq_along(selected_region_values)],
  panel_label = selected_region_labels,
  region_value = selected_region_values
)

region_map_settings <- list(
  "Arctic Ocean" = list(
    xlim = c(-180, 180), ylim = c(58, 85),
    shade_xlim = c(-180, 180), shade_ylim = c(58, 85),
    background_xlim = c(-180, 180), background_ylim = c(-55, 85),
    show_circulation = FALSE
  ),
  "North Atlantic Ocean" = list(
    xlim = c(-92, 18), ylim = c(8, 66),
    shade_xlim = c(-92, 18), shade_ylim = c(8, 66),
    background_xlim = c(-126, 52), background_ylim = c(8, 66),
    show_circulation = TRUE
  )
)

integrated_size <- c(width = 210 / 25.4, height = 280 / 25.4)

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
    Turnover = gamma - as.numeric(alpha[sample_ids])
  ) %>%
    left_join(
      panel_meta %>%
        select(SampleID, RegionShort, Ocean.region, Layer, all_of(covariate_cols)),
      by = "SampleID"
    )

  diversity_df %>%
    arrange(Turnover, Alpha, SampleID) %>%
    mutate(
      SamplePlotID = paste0("S", row_number()),
      SampleIndex = row_number()
    ) %>%
    select(Panel, PanelLabel, SampleID, SamplePlotID, SampleIndex, everything())
}

all_diversity <- purrr::pmap_dfr(
  panel_defs,
  function(panel, panel_label, region_value) {
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
    mean_turnover = mean(Turnover),
    .groups = "drop"
  )

write.csv(all_diversity, file.path(tab_dir, "figure2_alpha_gamma_turnover_values.csv"), row.names = FALSE)
write.csv(diversity_summary, file.path(tab_dir, "diversity_summary_HRIC.csv"), row.names = FALSE)

sample_order_check <- all_diversity %>%
  group_by(Panel, PanelLabel) %>%
  summarise(
    sorted_by_turnover = all(diff(Turnover) >= -1e-12),
    first_turnover = first(Turnover),
    last_turnover = last(Turnover),
    .groups = "drop"
  )

map_data <- ggplot2::map_data("world")

make_shade_rect <- function(region_label) {
  settings <- region_map_settings[[region_label]]
  shade_xlim <- if (!is.null(settings$shade_xlim)) settings$shade_xlim else settings$xlim
  shade_ylim <- if (!is.null(settings$shade_ylim)) settings$shade_ylim else settings$ylim
  tibble::tibble(
    xmin = shade_xlim[1],
    xmax = shade_xlim[2],
    ymin = shade_ylim[1],
    ymax = shade_ylim[2]
  )
}

make_hatch_segments <- function(region_label) {
  rect <- make_shade_rect(region_label)
  shade_width <- rect$xmax - rect$xmin
  shade_height <- rect$ymax - rect$ymin
  if (!is.finite(shade_width) || !is.finite(shade_height) || shade_width <= 0 || shade_height <= 0) {
    return(tibble::tibble(x = numeric(), y = numeric(), xend = numeric(), yend = numeric()))
  }
  hatch_dx <- shade_width * ifelse(region_label == "Arctic Ocean", 0.11, 0.22)
  spacing <- shade_width / ifelse(region_label == "Arctic Ocean", 34, 30)
  starts <- seq(rect$xmin - hatch_dx, rect$xmax, by = spacing)
  purrr::map_dfr(starts, function(x_bottom) {
    y_low <- max(rect$ymin, rect$ymin + shade_height * (rect$xmin - x_bottom) / hatch_dx)
    y_high <- min(rect$ymax, rect$ymin + shade_height * (rect$xmax - x_bottom) / hatch_dx)
    if (!is.finite(y_low) || !is.finite(y_high) || y_low >= y_high) {
      return(tibble::tibble(x = numeric(), y = numeric(), xend = numeric(), yend = numeric()))
    }
    tibble::tibble(
      x = x_bottom + hatch_dx * (y_low - rect$ymin) / shade_height,
      y = y_low,
      xend = x_bottom + hatch_dx * (y_high - rect$ymin) / shade_height,
      yend = y_high
    )
  })
}

north_atlantic_currents <- tibble::tribble(
  ~Current, ~CurrentType, ~Longitude, ~Latitude, ~PointOrder,
  "Gulf Stream", "Warm current", -81, 25, 1,
  "Gulf Stream", "Warm current", -80, 29, 2,
  "Gulf Stream", "Warm current", -76, 33, 3,
  "Gulf Stream", "Warm current", -70, 37, 4,
  "Gulf Stream", "Warm current", -63, 40, 5,
  "Gulf Stream", "Warm current", -57, 42, 6,
  "North Atlantic Current", "Warm current", -57, 42, 1,
  "North Atlantic Current", "Warm current", -49, 45, 2,
  "North Atlantic Current", "Warm current", -41, 48, 3,
  "North Atlantic Current", "Warm current", -33, 51, 4,
  "North Atlantic Current", "Warm current", -25, 54, 5,
  "North Atlantic Current", "Warm current", -17, 57, 6,
  "North Atlantic Current", "Warm current", -10, 59, 7,
  "North Equatorial Current", "Warm current", -15, 17, 1,
  "North Equatorial Current", "Warm current", -24, 16.5, 2,
  "North Equatorial Current", "Warm current", -34, 16.2, 3,
  "North Equatorial Current", "Warm current", -44, 16.5, 4,
  "North Equatorial Current", "Warm current", -54, 17, 5,
  "North Equatorial Current", "Warm current", -64, 18, 6,
  "Labrador Current", "Cold current", -56, 62, 1,
  "Labrador Current", "Cold current", -57, 59, 2,
  "Labrador Current", "Cold current", -56, 56, 3,
  "Labrador Current", "Cold current", -55, 53, 4,
  "Labrador Current", "Cold current", -54, 50, 5,
  "Labrador Current", "Cold current", -52, 47, 6,
  "Labrador Current", "Cold current", -49, 44, 7,
  "Canary Current", "Cold current", -12, 40, 1,
  "Canary Current", "Cold current", -13, 36, 2,
  "Canary Current", "Cold current", -14.5, 32, 3,
  "Canary Current", "Cold current", -16, 28, 4,
  "Canary Current", "Cold current", -17, 24, 5,
  "Canary Current", "Cold current", -16.5, 20, 6,
  "Canary Current", "Cold current", -15, 16, 7
)

north_atlantic_current_labels <- tibble::tribble(
  ~Current, ~CurrentType, ~LabelLongitude, ~LabelLatitude,
  "Gulf Stream", "Warm current", -82.0, 36.5,
  "North Atlantic Current", "Warm current", -35.5, 55.5,
  "North Equatorial Current", "Warm current", -43.0, 13.4,
  "Labrador Current", "Cold current", -65.0, 53.5,
  "Canary Current", "Cold current", -24.0, 29.5
)
make_region_map <- function(region_label, region_value, xlim, ylim, show_circulation = FALSE) {
  region_points <- meta %>%
    filter(Ocean.region == region_value)
  shade_fill <- unname(selected_region_ocean_fill[region_label])
  hatch_color <- unname(selected_region_mark[region_label])
  shade_rect <- make_shade_rect(region_label)
  hatch_segments <- if (region_label == "Arctic Ocean") {
    make_hatch_segments(region_label)
  } else {
    tibble::tibble(x = numeric(), y = numeric(), xend = numeric(), yend = numeric())
  }
  region_label_df <- tibble::tibble(
    Longitude = xlim[1] + 0.08 * diff(xlim),
    Latitude = ylim[2] - 0.12 * diff(ylim),
    label = paste0(region_label, "\n", "n = ", nrow(region_points))
  )

  p <- ggplot() +
    geom_rect(
      data = shade_rect,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = scales::alpha(shade_fill, 0.18), color = NA
    ) +
    geom_segment(
      data = hatch_segments,
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE, color = scales::alpha(hatch_color, 0.36), linewidth = 0.22
    ) +
    geom_polygon(
      data = map_data,
      aes(x = long, y = lat, group = group),
      fill = "#EFEFEB", color = "#D0D0CC", linewidth = 0.11
    )

  if (show_circulation) {
    p <- p +
      geom_path(
        data = north_atlantic_currents,
        aes(x = Longitude, y = Latitude, group = Current, color = CurrentType),
        linewidth = 0.72, lineend = "round", linejoin = "round", alpha = 0.94,
        arrow = grid::arrow(type = "closed", length = grid::unit(0.075, "in"))
      ) +
      geom_text(
        data = north_atlantic_current_labels,
        aes(x = LabelLongitude, y = LabelLatitude, label = Current, color = CurrentType),
        inherit.aes = FALSE, size = 2.15, fontface = "bold"
      ) +
      scale_color_manual(
        values = c("Warm current" = "#B6423A", "Cold current" = "#1F5A99"),
        name = NULL
      )
  }

  p +
    geom_point(
      data = region_points,
      aes(x = Longitude, y = Latitude),
      shape = 4, color = "#111111", stroke = 0.88, size = 2.20, alpha = 0.98
    ) +
    geom_label(
      data = region_label_df,
      aes(x = Longitude, y = Latitude, label = label),
      size = 2.85, linewidth = 0.18, label.padding = grid::unit(0.12, "lines"),
      fill = scales::alpha("white", 0.88), color = "grey8", hjust = 0
    ) +
    coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    theme_nature(base_size = 9.0) +
    theme(
      axis.line = element_blank(),
      legend.position = if (show_circulation) "bottom" else "none",
      legend.direction = "horizontal",
      legend.text = element_text(size = 8.0),
      panel.background = element_rect(fill = "#F4F8F8", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

rescale_to_range <- function(x, from, to) {
  to[1] + (x - from[1]) / diff(from) * diff(to)
}

make_figure2_background_data <- function(region_label, region_value, n_samples, y_upper) {
  settings <- region_map_settings[[region_label]]
  region_points <- meta %>%
    filter(Ocean.region == region_value)
  background_xlim <- if (!is.null(settings$background_xlim)) settings$background_xlim else settings$xlim
  background_ylim <- if (!is.null(settings$background_ylim)) settings$background_ylim else settings$ylim
  # Fill the complete Figure 2 data rectangle for every regional map.
  target_xlim <- c(0.5, n_samples + 0.5)
  target_ylim <- c(0, y_upper)
  shade_fill <- unname(selected_region_ocean_fill[region_label])
  mark_color <- scales::alpha("#111111", 0.56)
  shade_rect <- make_shade_rect(region_label)
  hatch_segments <- make_hatch_segments(region_label)

  transform_lon <- function(x) rescale_to_range(x, background_xlim, target_xlim)
  transform_lat <- function(y) rescale_to_range(y, background_ylim, target_ylim)

  map_background <- map_data %>%
    mutate(FigureX = transform_lon(long), FigureY = transform_lat(lat))

  shade_rect <- shade_rect %>%
    transmute(
      xmin = transform_lon(pmax(xmin, background_xlim[1])),
      xmax = transform_lon(pmin(xmax, background_xlim[2])),
      ymin = transform_lat(pmax(ymin, background_ylim[1])),
      ymax = transform_lat(pmin(ymax, background_ylim[2]))
    ) %>%
    filter(xmin < xmax, ymin < ymax)

  hatch_segments <- hatch_segments %>%
    mutate(
      FigureX = transform_lon(x),
      FigureY = transform_lat(y),
      FigureXEnd = transform_lon(xend),
      FigureYEnd = transform_lat(yend)
    )

  figure_points <- region_points %>%
    filter(
      Longitude >= background_xlim[1], Longitude <= background_xlim[2],
      Latitude >= background_ylim[1], Latitude <= background_ylim[2]
    ) %>%
    mutate(FigureX = transform_lon(Longitude), FigureY = transform_lat(Latitude))

  currents <- tibble::tibble(
    Current = character(), CurrentType = character(),
    FigureX = numeric(), FigureY = numeric(), PointOrder = integer()
  )
  if (isTRUE(settings$show_circulation)) {
    currents <- north_atlantic_currents %>%
      mutate(
        FigureX = transform_lon(Longitude),
        FigureY = transform_lat(Latitude)
      ) %>%
      arrange(Current, PointOrder)
  }

  list(
    ocean = tibble::tibble(
      xmin = target_xlim[1], xmax = target_xlim[2],
      ymin = target_ylim[1], ymax = target_ylim[2]
    ),
    shade = shade_rect,
    hatch = hatch_segments,
    map = map_background,
    points = figure_points,
    currents = currents,
    shade_fill = shade_fill,
    hatch_color = unname(selected_region_mark[region_label]),
    point_color = mark_color
  )
}

add_figure2_background_layers <- function(plot, bg) {
  plot +
    geom_rect(
      data = bg$ocean,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = scales::alpha("#F4F8F8", 0.14), color = NA
    ) +
    geom_rect(
      data = bg$shade,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = scales::alpha(bg$shade_fill, 0.13), color = NA
    ) +
    geom_polygon(
      data = bg$map,
      aes(x = FigureX, y = FigureY, group = group),
      inherit.aes = FALSE,
      fill = scales::alpha("#ECEBE6", 0.34),
      color = scales::alpha("#B9B9B3", 0.20),
      linewidth = 0.06
    ) +
    geom_segment(
      data = bg$hatch,
      aes(x = FigureX, y = FigureY, xend = FigureXEnd, yend = FigureYEnd),
      inherit.aes = FALSE, color = scales::alpha(bg$hatch_color, 0.25), linewidth = 0.16
    ) +
    geom_path(
      data = bg$currents,
      aes(x = FigureX, y = FigureY, group = Current, color = CurrentType),
      inherit.aes = FALSE, linewidth = 0.56, lineend = "round", alpha = 0.44,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.068, "in"))
    ) +
    geom_point(
      data = bg$points,
      aes(x = FigureX, y = FigureY),
      inherit.aes = FALSE, shape = 4, color = bg$point_color, stroke = 0.72, size = 1.70, alpha = 0.88
    ) +
    scale_color_manual(
      values = c("Warm current" = "#B6423A", "Cold current" = "#1F5A99"),
      guide = "none"
    )
}

figure1a_settings <- region_map_settings[["Arctic Ocean"]]
figure1b_settings <- region_map_settings[["North Atlantic Ocean"]]

figure1a <- make_region_map(
  "Arctic Ocean", selected_region_values[match("Arctic Ocean", selected_region_labels)],
  xlim = figure1a_settings$background_xlim, ylim = figure1a_settings$background_ylim,
  show_circulation = figure1a_settings$show_circulation
)
figure1b <- make_region_map(
  "North Atlantic Ocean", selected_region_values[match("North Atlantic Ocean", selected_region_labels)],
  xlim = figure1b_settings$background_xlim, ylim = figure1b_settings$background_ylim,
  show_circulation = figure1b_settings$show_circulation
)
figure1_plots <- list(a = figure1a, b = figure1b)

make_figure2 <- function(div_df, panel_label) {
  n_samples <- nrow(div_df)
  div_df <- div_df %>%
    mutate(
      SampleIndex = seq_len(n_samples),
      SamplePlotID = factor(SamplePlotID, levels = SamplePlotID)
    )
  gamma <- unique(div_df$Gamma)
  x_breaks <- sample_breaks(n_samples)
  x_labels <- paste0("S", x_breaks)
  axis_size <- ifelse(n_samples > 120, 7.2, ifelse(n_samples > 55, 7.7, 8.2))
  point_size <- ifelse(n_samples > 120, 1.55, 2.05)
  y_data <- c(div_df$Alpha, gamma)
  y_span <- diff(range(y_data, na.rm = TRUE))
  if (!is.finite(y_span) || y_span == 0) y_span <- 0.05
  y_lower <- max(0, min(y_data, na.rm = TRUE) - max(0.015, y_span * 0.08))
  y_upper <- min(1, max(y_data, na.rm = TRUE) + max(0.025, y_span * 0.16))
  y_breaks <- pretty(c(y_lower, y_upper), n = 5)
  y_breaks <- y_breaks[y_breaks >= y_lower & y_breaks <= y_upper]
  gamma_label_x <- n_samples - 0.25
  gamma_label_y <- min(y_upper - 0.012, gamma + max(0.012, y_span * 0.10))
  mark_color <- unname(selected_region_mark[panel_label])

  ggplot(div_df, aes(x = SampleIndex)) +
    geom_segment(
      aes(xend = SampleIndex, y = Alpha, yend = Gamma),
      linewidth = 0.26, color = "#B7B7B7", alpha = 0.72
    ) +
    geom_hline(yintercept = gamma, linetype = "22", linewidth = 0.62, color = "#9D1F2E") +
    geom_point(
      aes(y = Alpha),
      color = mark_color, size = point_size, alpha = 0.95
    ) +
    annotate(
      "label", x = gamma_label_x, y = gamma_label_y,
      label = paste0("gamma == ", sprintf("%.3f", gamma)),
      parse = TRUE, size = 2.85, fill = scales::alpha("white", 0.96), color = "#9D1F2E",
      hjust = 1
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = function(x) sprintf("%.1f", x),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = expansion(mult = c(0.003, 0.01))) +
    coord_cartesian(xlim = c(0.5, n_samples + 0.5), ylim = c(y_lower, y_upper), clip = "on") +
    labs(
      x = "Sample",
      y = "Alpha diversity (HRIC evenness)"
    ) +
    theme_nature(base_size = 9.0) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = axis_size),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.title = element_blank()
    )
}

figure2_plots <- list()
for (i in seq_len(nrow(panel_defs))) {
  panel <- panel_defs$panel[i]
  panel_label <- panel_defs$panel_label[i]
  div_df <- all_diversity %>% filter(Panel == panel)
  p <- make_figure2(div_df, panel_label)
  figure2_plots[[panel]] <- p
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
  select(Panel, PanelLabel, SampleID, SamplePlotID, Turnover, all_of(covariate_cols)) %>%
  tidyr::pivot_longer(all_of(covariate_cols), names_to = "Covariate", values_to = "CovariateValue") %>%
  group_by(Panel, PanelLabel, Covariate) %>%
  group_modify(~ cor_test_safe(.x$CovariateValue, .x$Turnover)) %>%
  ungroup() %>%
  mutate(
    CovariateLabel = covariate_labels[Covariate],
    label = paste0(
      "rho == ", sprintf("%.2f", rho),
      "*','~~", format_p_plotmath(p),
      "*','~~italic(n) == ", n
    )
  )

write.csv(figure3_stats, file.path(tab_dir, "figure3_spearman_turnover_covariates.csv"), row.names = FALSE)

make_figure3 <- function(div_df, stats_df) {
  facet_columns <- ifelse(length(covariate_cols) <= 6, 3, 5)
  label_size <- ifelse(length(covariate_cols) > 20, 2.05, 2.70)
  turnover_limits <- range(div_df$Turnover, na.rm = TRUE)
  turnover_pad <- diff(turnover_limits) * 0.18
  if (!is.finite(turnover_pad) || turnover_pad == 0) turnover_pad <- 0.05
  long_df <- div_df %>%
    select(SampleID, SamplePlotID, Turnover, all_of(covariate_cols)) %>%
    tidyr::pivot_longer(all_of(covariate_cols), names_to = "Covariate", values_to = "CovariateValue") %>%
    mutate(CovariateLabel = factor(covariate_labels[Covariate], levels = covariate_labels))
  smooth_df <- long_df %>%
    group_by(Covariate, CovariateLabel) %>%
    filter(
      sum(is.finite(CovariateValue) & is.finite(Turnover)) >= 4,
      length(unique(CovariateValue[is.finite(CovariateValue) & is.finite(Turnover)])) >= 2
    ) %>%
    ungroup()
  stats_df <- stats_df %>%
    mutate(CovariateLabel = factor(covariate_labels[Covariate], levels = covariate_labels)) %>%
    left_join(
      long_df %>%
        filter(is.finite(CovariateValue), is.finite(Turnover)) %>%
        group_by(Covariate) %>%
        summarise(
          x_mid = mean(range(CovariateValue)),
          y_mid = mean(range(Turnover)),
          top_left = sum(CovariateValue <= x_mid & Turnover >= y_mid),
          top_right = sum(CovariateValue > x_mid & Turnover >= y_mid),
          .groups = "drop"
        ),
      by = "Covariate"
    ) %>%
    mutate(
      LabelX = ifelse(top_left <= top_right, -Inf, Inf),
      LabelHjust = ifelse(top_left <= top_right, -0.04, 1.04)
    )

  ggplot(long_df, aes(x = CovariateValue, y = Turnover)) +
    geom_smooth(
      data = smooth_df,
      method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.82,
      color = "#9D1F2E", fill = "#9D1F2E", alpha = 0.12, na.rm = TRUE
    ) +
    geom_label(
      data = stats_df,
      aes(x = LabelX, y = Inf, label = label, hjust = LabelHjust),
      inherit.aes = FALSE, vjust = 1.10, size = label_size, color = "grey12",
      fill = scales::alpha("white", 0.80), linewidth = 0,
      label.padding = grid::unit(0.045, "lines"),
      label.r = grid::unit(0, "lines"), parse = TRUE
    ) +
    facet_wrap(~ CovariateLabel, scales = "free_x", ncol = facet_columns) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
    coord_cartesian(ylim = c(turnover_limits[1] - turnover_pad, turnover_limits[2] + turnover_pad)) +
    labs(x = NULL, y = expression(atop("local-regional diversity gap", gamma - alpha))) +
    theme_nature(base_size = 9.0) +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    )
}

figure3_plots <- list()
for (i in seq_len(nrow(panel_defs))) {
  panel <- panel_defs$panel[i]
  div_df <- all_diversity %>% filter(Panel == panel)
  stats_df <- figure3_stats %>% filter(Panel == panel)
  p <- make_figure3(div_df, stats_df)
  figure3_plots[[panel]] <- p
}

taxa_colors <- c(
  "#1B4F9C", "#B72E2B", "#008A5B", "#7B3294", "#E69F00",
  "#0099B4", "#C44E8B", "#8C6D31", "#4B8B3B", "#4B4B4B",
  "#00A087", "#F39B30", "#3C5488", "#A73030", "#6E4D9B",
  "#D55E00", "#4DBBD5", "#55752F", "#9E9E2D"
)

get_top_taxa <- function(sample_ids) {
  counts_panel <- counts_taxon[sample_ids, , drop = FALSE]
  counts_panel <- counts_panel[, colSums(counts_panel) > 0, drop = FALSE]
  rel_panel <- sweep(counts_panel, 1, rowSums(counts_panel), FUN = "/")
  mean_rel <- sort(colMeans(rel_panel, na.rm = TRUE), decreasing = TRUE)
  names(mean_rel)[seq_len(min(19, length(mean_rel)))]
}

figure4_panel_top_taxa <- setNames(
  lapply(panel_defs$panel, function(panel) {
    sample_ids <- all_diversity %>%
      filter(Panel == panel) %>%
      arrange(SampleIndex) %>%
      pull(SampleID)
    get_top_taxa(sample_ids)
  }),
  panel_defs$panel
)
figure4_shared_taxa <- Reduce(intersect, figure4_panel_top_taxa)
figure4_shared_taxa <- figure4_panel_top_taxa[[1]][figure4_panel_top_taxa[[1]] %in% figure4_shared_taxa]
figure4_panel_palettes <- setNames(
  lapply(panel_defs$panel, function(panel) {
    taxon_order <- c(figure4_shared_taxa, setdiff(figure4_panel_top_taxa[[panel]], figure4_shared_taxa))
    setNames(taxa_colors[seq_along(taxon_order)], taxon_order)
  }),
  panel_defs$panel
)
figure4_color_key <- purrr::imap_dfr(figure4_panel_palettes, function(palette, panel) {
  tibble::tibble(
    Panel = panel,
    PanelLabel = panel_defs$panel_label[match(panel, panel_defs$panel)],
    Taxon = names(palette),
    Color = unname(palette),
    SharedAcrossPanels = names(palette) %in% figure4_shared_taxa
  )
})
write.csv(figure4_color_key, file.path(tab_dir, "figure4_taxon_color_key.csv"), row.names = FALSE)

lm_taxa <- function(counts_taxon_panel, panel_meta, top_taxa, covariate) {
  x <- counts_taxon_panel[, colSums(counts_taxon_panel) > 0, drop = FALSE]
  hric_x <- HRIC::HRIC(x)
  top_taxa <- intersect(top_taxa, colnames(hric_x))
  out <- purrr::map_dfr(top_taxa, function(taxon) {
    values <- as_numeric_safe(panel_meta[[covariate]])
    y <- hric_x[panel_meta$SampleID, taxon]
    keep <- is.finite(values) & is.finite(y)
    if (sum(keep) < 4 || length(unique(values[keep])) < 2 || length(unique(y[keep])) < 2) {
      return(tibble::tibble(
        Taxon = taxon,
        Covariate = covariate,
        n = sum(keep),
        slope = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p = NA_real_,
        r_squared = NA_real_
      ))
    }
    fit <- tryCatch(lm(y[keep] ~ values[keep]), error = function(e) NULL)
    if (is.null(fit)) {
      return(tibble::tibble(
        Taxon = taxon,
        Covariate = covariate,
        n = sum(keep),
        slope = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p = NA_real_,
        r_squared = NA_real_
      ))
    }
    fit_summary <- summary(fit)
    coef_table <- fit_summary$coefficients
    tibble::tibble(
      Taxon = taxon,
      Covariate = covariate,
      n = sum(keep),
      slope = unname(coef_table[2, "Estimate"]),
      std_error = unname(coef_table[2, "Std. Error"]),
      statistic = unname(coef_table[2, "t value"]),
      p = unname(coef_table[2, "Pr(>|t|)"]),
      r_squared = unname(fit_summary$r.squared)
    )
  })
  out %>%
    mutate(
      q = p.adjust(p, method = "BH"),
      CovariateLabel = unname(covariate_labels[Covariate]),
      Significant = !is.na(q) & q < 0.05
    )
}

make_figure4 <- function(panel_meta, panel_letter, panel_label, daa_covariate) {
  sample_ids <- panel_meta$SampleID
  counts_panel <- counts_taxon[sample_ids, , drop = FALSE]
  counts_panel <- counts_panel[, colSums(counts_panel) > 0, drop = FALSE]
  rel_panel <- sweep(counts_panel, 1, rowSums(counts_panel), FUN = "/")
  mean_rel <- sort(colMeans(rel_panel, na.rm = TRUE), decreasing = TRUE)
  top_taxa <- names(mean_rel)[seq_len(min(19, length(mean_rel)))]
  if (!identical(top_taxa, figure4_panel_top_taxa[[panel_letter]])) {
    stop("Figure 4 top-taxon order changed unexpectedly for panel ", panel_letter)
  }

  tests <- lm_taxa(counts_panel, panel_meta, top_taxa, daa_covariate)
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
      SamplePlotID = paste0("S", match(SampleID, sample_ids)),
      SamplePlotID = factor(SamplePlotID, levels = paste0("S", seq_along(sample_ids))),
      Taxon = factor(Taxon, levels = c(top_taxa, "Others"))
    )

  legend_labels <- setNames(levels(stack_df$Taxon), levels(stack_df$Taxon))
  legend_labels[significant_taxa] <- paste0(legend_labels[significant_taxa], "*")
  pal <- c(figure4_panel_palettes[[panel_letter]], Others = "#CFCFCF")
  pal <- pal[levels(stack_df$Taxon)]
  n_samples <- length(sample_ids)
  x_breaks <- paste0("S", sample_breaks(n_samples))
  axis_size <- ifelse(n_samples > 120, 7.0, ifelse(n_samples > 55, 7.4, 7.8))

  plot <- ggplot(stack_df, aes(x = SamplePlotID, y = RelativeAbundance, fill = Taxon)) +
    geom_col(width = 0.98, linewidth = 0.055, color = "white") +
    scale_fill_manual(values = pal, labels = legend_labels, drop = FALSE, name = NULL) +
    scale_x_discrete(breaks = x_breaks, labels = x_breaks, drop = FALSE) +
    scale_y_continuous(
      breaks = seq(0, 100, 25),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(ylim = c(0, 100), expand = FALSE) +
    labs(x = "Sample", y = "Relative abundance") +
    guides(fill = guide_legend(
      ncol = 5, byrow = TRUE, title.position = "top", title.hjust = 0,
      override.aes = list(linewidth = 0)
    )) +
    theme_nature(base_size = 9.0) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = axis_size),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 7.4),
      legend.text = element_text(size = 7.0),
      legend.key.size = grid::unit(3.7, "mm"),
      legend.spacing.x = grid::unit(1.3, "mm"),
      panel.grid = element_blank()
    )

  list(
    plot = plot,
    top_taxa = tibble::tibble(
      Panel = panel_letter,
      PanelLabel = panel_label,
      Covariate = daa_covariate,
      CovariateLabel = unname(covariate_labels[daa_covariate]),
      Taxon = top_taxa,
      mean_relative_abundance = as.numeric(mean_rel[top_taxa]),
      significant_by_HRIC_lm = top_taxa %in% significant_taxa
    ),
    tests = tests %>% mutate(Panel = panel_letter, PanelLabel = panel_label, .before = 1),
    stack = stack_df %>% mutate(
      Panel = panel_letter,
      PanelLabel = panel_label,
      Covariate = daa_covariate,
      CovariateLabel = unname(covariate_labels[daa_covariate]),
      .before = 1
    )
  )
}

figure4_plots_by_panel_covariate <- setNames(vector("list", nrow(panel_defs)), panel_defs$panel)
figure4_top_taxa <- list()
figure4_lm_results <- list()
figure4_stack_data <- list()
for (i in seq_len(nrow(panel_defs))) {
  panel <- panel_defs$panel[i]
  panel_label <- panel_defs$panel_label[i]
  region_value <- panel_defs$region_value[i]
  panel_order <- all_diversity %>%
    filter(Panel == panel) %>%
    arrange(SampleIndex) %>%
    pull(SampleID)
  panel_meta <- get_panel_meta(region_value)
  panel_meta <- panel_meta[match(panel_order, panel_meta$SampleID), , drop = FALSE]
  if (!identical(panel_meta$SampleID, panel_order)) {
    stop("Figure 4 sample order does not match Figure 2 for panel ", panel)
  }
  figure4_plots_by_panel_covariate[[panel]] <- list()
  for (covar in covariate_cols) {
    result <- make_figure4(panel_meta, panel, panel_label, covar)
    figure4_plots_by_panel_covariate[[panel]][[covar]] <- result$plot
    result_key <- paste(panel, covar, sep = "__")
    figure4_top_taxa[[result_key]] <- result$top_taxa
    figure4_lm_results[[result_key]] <- result$tests
    figure4_stack_data[[result_key]] <- result$stack
  }
}

write.csv(bind_rows(figure4_top_taxa), file.path(tab_dir, "figure4_top19_taxa_by_panel_covariate.csv"), row.names = FALSE)
figure4_lm_table <- bind_rows(figure4_lm_results)
write.csv(figure4_lm_table, file.path(tab_dir, "figure4_HRIC_lm_by_panel_covariate.csv"), row.names = FALSE)
write.csv(bind_rows(figure4_stack_data), file.path(tab_dir, "figure4_stacked_bar_plot_data_by_covariate.csv"), row.names = FALSE)

figure4_selection_stats <- figure4_lm_table %>%
  group_by(Covariate, CovariateLabel) %>%
  summarise(
    n_significant_tests = sum(Significant, na.rm = TRUE),
    n_significant_taxa = n_distinct(Taxon[Significant]),
    min_q = ifelse(all(is.na(q)), NA_real_, min(q, na.rm = TRUE)),
    median_q = ifelse(all(is.na(q)), NA_real_, median(q, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_significant_tests), desc(n_significant_taxa), is.na(min_q), min_q, CovariateLabel)
main_figure4_covariate <- figure4_selection_stats$Covariate[1]
main_figure4_covariate_label <- unname(covariate_labels[main_figure4_covariate])
figure4_selection_stats <- figure4_selection_stats %>%
  mutate(selected_for_integrated_figure = Covariate == main_figure4_covariate)
write.csv(figure4_selection_stats, file.path(tab_dir, "figure4_lm_covariate_selection.csv"), row.names = FALSE)

figure4_plots <- setNames(
  lapply(panel_defs$panel, function(panel) {
    figure4_plots_by_panel_covariate[[panel]][[main_figure4_covariate]]
  }),
  panel_defs$panel
)

figure4_order_check <- bind_rows(figure4_stack_data) %>%
  distinct(Panel, Covariate, SamplePlotID, SampleID) %>%
  left_join(
    all_diversity %>% select(Panel, SamplePlotID, Figure2SampleID = SampleID),
    by = c("Panel", "SamplePlotID")
  ) %>%
  summarise(matches_figure2_order = all(SampleID == Figure2SampleID), .groups = "drop") %>%
  pull(matches_figure2_order)

make_integrated_region <- function(panel) {
  panel_label <- panel_defs$panel_label[match(panel, panel_defs$panel)]
  panel_n <- sum(all_diversity$Panel == panel)
  panel_tags <- if (panel == "a") letters[1:4] else letters[5:8]
  compact_theme <- theme(
    axis.title = element_text(size = 7.0, face = "bold"),
    axis.text = element_text(size = 5.8),
    legend.title = element_text(size = 6.0, face = "bold"),
    legend.text = element_text(size = 5.2),
    strip.text = element_text(size = 5.8, face = "bold"),
    plot.title = element_text(size = 7.0, face = "bold"),
    plot.subtitle = element_text(size = 6.2, face = "bold"),
    plot.tag = element_text(size = 9.0, face = "bold", color = "grey5"),
    plot.tag.position = "topleft"
  )
  map_plot <- figure1_plots[[panel]] +
    labs(tag = panel_tags[1]) +
    compact_theme +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 5.2),
      plot.margin = margin(1, 0, 1, 2)
    )
  map_plot$layers <- Filter(
    function(layer) !inherits(layer$geom, "GeomLabel"),
    map_plot$layers
  )
  for (i in seq_along(map_plot$layers)) {
    layer <- map_plot$layers[[i]]
    if (inherits(layer$geom, "GeomText")) {
      map_plot$layers[[i]]$aes_params$size <- 1.42
    } else if (inherits(layer$geom, "GeomPath")) {
      map_plot$layers[[i]]$aes_params$linewidth <- 0.55
    } else if (inherits(layer$geom, "GeomPoint")) {
      map_plot$layers[[i]]$aes_params$size <- 1.55
      map_plot$layers[[i]]$aes_params$stroke <- 0.70
    }
  }
  alpha_plot <- figure2_plots[[panel]] +
    labs(y = "Alpha diversity", tag = panel_tags[2]) +
    compact_theme +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.2),
      axis.text.y = element_text(size = 5.4),
      plot.margin = margin(1, 2, 1, 0)
    )
  gamma_value <- unique(all_diversity$Gamma[all_diversity$Panel == panel])
  for (i in seq_along(alpha_plot$layers)) {
    if (inherits(alpha_plot$layers[[i]]$geom, "GeomLabel")) {
      alpha_plot$layers[[i]]$data$y <- gamma_value + 0.004
    }
  }
  taxa_plot <- figure4_plots[[panel]] +
    labs(tag = panel_tags[3]) +
    compact_theme +
    guides(fill = guide_legend(
      ncol = 5, byrow = TRUE, title.position = "top", title.hjust = 0,
      override.aes = list(linewidth = 0)
    )) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.0),
      axis.text.y = element_text(size = 5.2),
      legend.key.size = grid::unit(2.4, "mm"),
      legend.spacing.x = grid::unit(0.5, "mm"),
      legend.spacing.y = grid::unit(0, "mm"),
      legend.box.spacing = grid::unit(0.5, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(2, 2, 2, 2)
    )
  covariate_plot <- figure3_plots[[panel]] +
    labs(tag = panel_tags[4]) +
    facet_wrap(~ CovariateLabel, scales = "free_x", nrow = 1) +
    compact_theme +
    theme(
      axis.text = element_text(size = 5.0),
      axis.title.y = element_text(size = 6.5, face = "bold"),
      strip.text = element_text(size = 5.2, face = "bold"),
      plot.margin = margin(2, 2, 2, 2)
    )
  for (i in seq_along(covariate_plot$layers)) {
    if (inherits(covariate_plot$layers[[i]]$geom, "GeomLabel")) {
      covariate_plot$layers[[i]]$aes_params$size <- 1.35
      covariate_plot$layers[[i]]$data$LabelX <- covariate_plot$layers[[i]]$data$x_mid
      covariate_plot$layers[[i]]$data$LabelHjust <- 0.5
    }
  }

  first_row <- wrap_plots(
    map_plot,
    alpha_plot,
    nrow = 1,
    widths = c(0.50, 0.50)
  )

  region_header <- ggplot() +
    annotate(
      "label", x = 0, y = 0.5,
      label = paste0(panel_label, "   n = ", panel_n),
      hjust = 0, vjust = 0.5, size = 3.8, fontface = "bold",
      color = "grey8", fill = "white", linewidth = 0.45,
      label.padding = grid::unit(0.16, "lines")
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(1, 2, 1, 7))

  wrap_plots(
    region_header,
    first_row,
    taxa_plot,
    covariate_plot,
    ncol = 1,
    heights = c(0.22, 1.45, 1.35, 1.30)
  )
}

integrated_figure <- wrap_plots(
  make_integrated_region("a"),
  make_integrated_region("b"),
  ncol = 1,
  heights = c(1, 1)
) &
  theme(plot.background = element_rect(fill = "white", color = NA))

save_plot(
  integrated_figure,
  "figure_integrated_ocean_regions",
  width = integrated_size["width"],
  height = integrated_size["height"]
)

sample_panel_key <- all_diversity %>%
  select(
    Panel, PanelLabel, SamplePlotID, SampleIndex, SampleID, RegionShort, Ocean.region,
    Layer, Turnover, Alpha, Gamma, all_of(covariate_cols)
  )
write.csv(sample_panel_key, file.path(tab_dir, "sample_panel_key.csv"), row.names = FALSE)

stopifnot(
  all(c("HRIC", "SHalpha", "SHgamma", "SHbeta") %in% getNamespaceExports("HRIC")),
  all(abs(diversity_summary$gamma - diversity_summary$mean_alpha - diversity_summary$beta) < 1e-10),
  max(abs(all_diversity$Turnover - (all_diversity$Gamma - all_diversity$Alpha))) < 1e-12,
  all(sample_order_check$sorted_by_turnover),
  isTRUE(figure4_order_check),
  file.exists(file.path(fig_dir, "figure_integrated_ocean_regions.pdf"))
)

message("Analysis complete. Results written to: ", out_dir)
