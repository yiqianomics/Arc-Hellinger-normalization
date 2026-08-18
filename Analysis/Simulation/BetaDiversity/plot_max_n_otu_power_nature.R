#!/usr/bin/env Rscript

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(ragg)
})

if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Package 'mgcv' is required to fit the power curves.")
}

script_dir_override <- Sys.getenv("SIM_POWER_DIR", unset = "")
script_dir <- if (nzchar(script_dir_override)) {
  normalizePath(script_dir_override)
} else {
  get_script_dir()
}
outdir <- file.path(script_dir, "weighted_power_output")
input_file <- file.path(outdir, "max_n_otu_cluster_power.csv")

if (!file.exists(input_file)) {
  stop("Power table not found: ", input_file)
}

power_raw <- read.csv(input_file, stringsAsFactors = FALSE)

required_columns <- c(
  "k", "signal_mode", "beta", "method", "power", "n_nonmissing"
)
missing_columns <- setdiff(required_columns, names(power_raw))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

# The ordering creates three coherent legend rows: focal method, conventional
# ecological distances, Aitchison distances, and ILR-Euclidean distances.
method_order <- c(
  "AHC_euclidean",
  "bray",
  "euclidean",
  "jaccard",
  "aitchison_pc0.01",
  "aitchison_pc0.1",
  "aitchison",
  "aitchison_pc1",
  "ILR_euclidean_pc0.01",
  "ILR_euclidean_pc0.1",
  "ILR_euclidean_pc0.5",
  "ILR_euclidean_pc1"
)

method_labels <- c(
  "AHC_euclidean" = "HRIC distance",
  "bray" = "Bray-Curtis",
  "euclidean" = "Euclidean",
  "jaccard" = "Jaccard",
  "aitchison_pc0.01" = "Aitchison (0.01)",
  "aitchison_pc0.1" = "Aitchison (0.1)",
  "aitchison" = "Aitchison (0.5)",
  "aitchison_pc1" = "Aitchison (1)",
  "ILR_euclidean_pc0.01" = "ILR-Euclidean (0.01)",
  "ILR_euclidean_pc0.1" = "ILR-Euclidean (0.1)",
  "ILR_euclidean_pc0.5" = "ILR-Euclidean (0.5)",
  "ILR_euclidean_pc1" = "ILR-Euclidean (1)"
)

method_colors <- c(
  "AHC_euclidean" = "#111111",
  "bray" = "#D55E00",
  "euclidean" = "#A6508A",
  "jaccard" = "#7A7A7A",
  "aitchison_pc0.01" = "#73BFE2",
  "aitchison_pc0.1" = "#3B8BC2",
  "aitchison" = "#1768AC",
  "aitchison_pc1" = "#244A73",
  "ILR_euclidean_pc0.01" = "#57B894",
  "ILR_euclidean_pc0.1" = "#2A9D8F",
  "ILR_euclidean_pc0.5" = "#4E8B57",
  "ILR_euclidean_pc1" = "#6F7F2A"
)

method_linetypes <- c(
  "AHC_euclidean" = "solid",
  "bray" = "solid",
  "euclidean" = "longdash",
  "jaccard" = "dotted",
  "aitchison_pc0.01" = "solid",
  "aitchison_pc0.1" = "longdash",
  "aitchison" = "dotdash",
  "aitchison_pc1" = "dotted",
  "ILR_euclidean_pc0.01" = "solid",
  "ILR_euclidean_pc0.1" = "longdash",
  "ILR_euclidean_pc0.5" = "dotdash",
  "ILR_euclidean_pc1" = "dotted"
)

signal_labels <- c(
  "abundance" = "Abundance",
  "both" = "Joint",
  "prevalence" = "Prevalence"
)

power_plot_data <- power_raw |>
  filter(
    k %in% c(2, 4, 6, 8),
    beta > 0,
    method %in% method_order
  ) |>
  mutate(
    method = factor(method, levels = method_order),
    signal_mode = factor(
      signal_mode,
      levels = c("abundance", "both", "prevalence")
    ),
    signal_label = factor(
      unname(signal_labels[as.character(signal_mode)]),
      levels = unname(signal_labels[c("abundance", "both", "prevalence")])
    ),
    k_label = factor(
      paste0("k = ", k),
      levels = paste0("k = ", c(2, 4, 6, 8))
    )
  )

# A binomial generalized additive model smooths Monte Carlo variation in each
# empirical power series. cummax removes only negligible fitted reversals so
# that increasing effect sizes cannot produce a visually decreasing guide.
smooth_power_curve <- function(dat, n_grid = 401L) {
  dat <- dat[order(dat$beta), , drop = FALSE]
  grid <- seq(min(dat$beta), max(dat$beta), length.out = n_grid)

  successes <- round(dat$power * dat$n_nonmissing)
  model <- mgcv::gam(
    cbind(successes, dat$n_nonmissing - successes) ~
      s(beta, k = 5, bs = "cs"),
    data = dat,
    family = stats::binomial(),
    method = "REML"
  )

  fitted_power <- stats::predict(
    model,
    newdata = data.frame(beta = grid),
    type = "response"
  )

  data.frame(
    beta = grid,
    power_curve = cummax(pmin(pmax(fitted_power, 0), 1))
  )
}

power_curves <- power_plot_data |>
  group_by(k, k_label, signal_mode, signal_label, method) |>
  group_modify(~ smooth_power_curve(.x)) |>
  ungroup()

plot_base <- ggplot(
  power_curves,
  aes(
    x = beta,
    y = power_curve,
    colour = method,
    linetype = method,
    group = method
  )
) +
  geom_hline(
    yintercept = 0.8,
    colour = "#D3D3D3",
    linewidth = 0.20,
    linetype = "22"
  ) +
  geom_line(linewidth = 0.30, lineend = "round") +
  geom_line(
    data = filter(power_curves, method == "AHC_euclidean"),
    linewidth = 0.42,
    lineend = "round",
    show.legend = FALSE
  ) +
  facet_grid(
    rows = vars(signal_label),
    cols = vars(k_label),
    switch = "y"
  ) +
  scale_colour_manual(
    values = method_colors,
    breaks = method_order,
    labels = unname(method_labels[method_order]),
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = method_linetypes,
    breaks = method_order,
    labels = unname(method_labels[method_order]),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0.05, 1),
    breaks = c(0.25, 0.50, 0.75, 1.00),
    labels = c("0.25", "0.50", "0.75", "1.00"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = c("0", "0.25", "0.50", "0.75", "1"),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(
    x = expression("Effect size, " * beta),
    y = "Power",
    colour = NULL,
    linetype = NULL
  ) +
  guides(
    colour = guide_legend(
      nrow = 3,
      byrow = TRUE,
      override.aes = list(linewidth = 0.45, alpha = 1)
    )
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    panel.border = element_rect(
      colour = "#333333",
      fill = NA,
      linewidth = 0.28
    ),
    panel.grid.major.y = element_line(
      colour = "#E7E7E7",
      linewidth = 0.18
    ),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.30),
    axis.ticks.length = unit(1.8, "pt"),
    axis.title = element_text(size = 7, colour = "#111111"),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4)),
    axis.text = element_text(size = 6.2, colour = "#222222"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(
      size = 7,
      face = "bold",
      colour = "#111111",
      margin = margin(t = 2.5, b = 2.5)
    ),
    strip.text.y.left = element_text(
      size = 7,
      face = "bold",
      angle = 90,
      colour = "#111111",
      margin = margin(l = 2.5, r = 2.5)
    ),
    panel.spacing.x = unit(4.2, "pt"),
    panel.spacing.y = unit(4.2, "pt"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 6.0, colour = "#222222"),
    legend.key = element_blank(),
    legend.key.width = unit(18, "pt"),
    legend.key.height = unit(8.5, "pt"),
    legend.spacing.x = unit(4, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.box.margin = margin(t = 3),
    plot.margin = margin(t = 2, r = 5, b = 1, l = 2)
  )

png_file <- file.path(
  outdir,
  "max_n_otu_cluster_power_lineplot_smooth.png"
)
pdf_file <- file.path(
  outdir,
  "max_n_otu_cluster_power_lineplot_smooth.pdf"
)
tiff_file <- file.path(
  outdir,
  "max_n_otu_cluster_power_lineplot_smooth.tiff"
)
source_data_file <- file.path(
  outdir,
  "max_n_otu_cluster_power_source_data.csv"
)

ggsave(
  png_file,
  plot = plot_base,
  device = ragg::agg_png,
  width = 180,
  height = 160,
  units = "mm",
  dpi = 600,
  bg = "white"
)

ggsave(
  pdf_file,
  plot = plot_base,
  device = grDevices::cairo_pdf,
  width = 180,
  height = 160,
  units = "mm",
  bg = "white"
)

ggsave(
  tiff_file,
  plot = plot_base,
  device = ragg::agg_tiff,
  width = 180,
  height = 160,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

source_data <- power_plot_data |>
  transmute(
    panel = as.character(signal_label),
    k,
    beta,
    method = unname(method_labels[as.character(method)]),
    empirical_power = power,
    n_replicates = n_nonmissing
  )
write.csv(source_data, source_data_file, row.names = FALSE)

message("Saved:")
message("  ", png_file)
message("  ", pdf_file)
message("  ", tiff_file)
message("  ", source_data_file)
