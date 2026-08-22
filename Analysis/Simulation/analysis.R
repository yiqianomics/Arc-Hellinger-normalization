#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

required_packages <- c("dplyr", "ggplot2", "mgcv", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(scales)
})

script_path <- grep("^--file=", commandArgs(FALSE), value = TRUE)
analysis_dir <- if (length(script_path) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_path[1])))
} else {
  normalizePath(getwd())
}

data_dir <- file.path(analysis_dir, "data")
result_dir <- file.path(analysis_dir, "results")
unlink(result_dir, recursive = TRUE, force = TRUE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

power_file <- file.path(data_dir, "max_n_otu_cluster_power.csv")
manova_file <- file.path(data_dir, "table_replication_level_hric_manova.csv")
theory_file <- file.path(data_dir, "table_hric_shbeta_theory_check.csv")
input_files <- c(power_file, manova_file, theory_file)
if (any(!file.exists(input_files))) {
  stop("Missing input files: ", paste(input_files[!file.exists(input_files)], collapse = ", "))
}

method_order <- c(
  "AHC_euclidean", "bray", "euclidean", "jaccard",
  "aitchison_pc0.01", "aitchison_pc0.1", "aitchison", "aitchison_pc1",
  "ILR_euclidean_pc0.01", "ILR_euclidean_pc0.1",
  "ILR_euclidean_pc0.5", "ILR_euclidean_pc1"
)
method_labels <- c(
  AHC_euclidean = "HRIC distance",
  bray = "Bray-Curtis",
  euclidean = "Euclidean",
  jaccard = "Jaccard",
  aitchison_pc0.01 = "Aitchison (0.01)",
  aitchison_pc0.1 = "Aitchison (0.1)",
  aitchison = "Aitchison (0.5)",
  aitchison_pc1 = "Aitchison (1)",
  ILR_euclidean_pc0.01 = "ILR-Euclidean (0.01)",
  ILR_euclidean_pc0.1 = "ILR-Euclidean (0.1)",
  ILR_euclidean_pc0.5 = "ILR-Euclidean (0.5)",
  ILR_euclidean_pc1 = "ILR-Euclidean (1)"
)
method_colors <- c(
  AHC_euclidean = "#111111",
  bray = "#D55E00",
  euclidean = "#A6508A",
  jaccard = "#7A7A7A",
  aitchison_pc0.01 = "#73BFE2",
  aitchison_pc0.1 = "#3B8BC2",
  aitchison = "#1768AC",
  aitchison_pc1 = "#244A73",
  ILR_euclidean_pc0.01 = "#57B894",
  ILR_euclidean_pc0.1 = "#2A9D8F",
  ILR_euclidean_pc0.5 = "#4E8B57",
  ILR_euclidean_pc1 = "#6F7F2A"
)
method_linetypes <- c(
  AHC_euclidean = "solid",
  bray = "solid",
  euclidean = "longdash",
  jaccard = "dotted",
  aitchison_pc0.01 = "solid",
  aitchison_pc0.1 = "longdash",
  aitchison = "dotdash",
  aitchison_pc1 = "dotted",
  ILR_euclidean_pc0.01 = "solid",
  ILR_euclidean_pc0.1 = "longdash",
  ILR_euclidean_pc0.5 = "dotdash",
  ILR_euclidean_pc1 = "dotted"
)
signal_labels <- c(
  abundance = "Abundance",
  both = "Joint",
  prevalence = "Prevalence"
)

power_raw <- read.csv(power_file)
required_power_columns <- c(
  "k", "signal_mode", "beta", "method", "power", "n_nonmissing"
)
if (length(setdiff(required_power_columns, names(power_raw))) > 0L) {
  stop("The power input does not contain the required columns.")
}

power_data <- power_raw |>
  filter(k %in% c(2, 4, 6, 8), beta > 0, method %in% method_order) |>
  mutate(
    method = factor(method, levels = method_order),
    signal_mode = factor(signal_mode, levels = c("abundance", "both", "prevalence")),
    signal_label = factor(
      unname(signal_labels[as.character(signal_mode)]),
      levels = unname(signal_labels[c("abundance", "both", "prevalence")])
    ),
    k_label = factor(paste0("k = ", k), levels = paste0("k = ", c(2, 4, 6, 8)))
  )

# The smooths summarize Monte Carlo rejection rates without replacing the
# empirical values used in the simulation analysis.
smooth_power_curve <- function(data, n_grid = 401L) {
  data <- data[order(data$beta), , drop = FALSE]
  grid <- seq(min(data$beta), max(data$beta), length.out = n_grid)
  successes <- round(data$power * data$n_nonmissing)
  model <- mgcv::gam(
    cbind(successes, data$n_nonmissing - successes) ~ s(beta, k = 5, bs = "cs"),
    data = data,
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

power_curves <- power_data |>
  group_by(k, k_label, signal_mode, signal_label, method) |>
  group_modify(~ smooth_power_curve(.x)) |>
  ungroup()

power_plot <- ggplot(
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
  facet_grid(rows = vars(signal_label), cols = vars(k_label), switch = "y") +
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
    panel.border = element_rect(colour = "#333333", fill = NA, linewidth = 0.28),
    panel.grid.major.y = element_line(colour = "#E7E7E7", linewidth = 0.18),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.30),
    axis.ticks.length = unit(1.8, "pt"),
    axis.title = element_text(size = 6.5, colour = "#111111"),
    axis.title.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(margin = margin(r = 4)),
    axis.text = element_text(size = 5.6, colour = "#222222"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(
      size = 6.4,
      face = "bold",
      colour = "#111111",
      margin = margin(t = 2.5, b = 2.5)
    ),
    strip.text.y.left = element_text(
      size = 6.4,
      face = "bold",
      angle = 90,
      colour = "#111111",
      margin = margin(l = 2.5, r = 2.5)
    ),
    panel.spacing.x = unit(3.5, "pt"),
    panel.spacing.y = unit(3.5, "pt"),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 5.25, colour = "#222222"),
    legend.key = element_blank(),
    legend.key.width = unit(14, "pt"),
    legend.key.height = unit(6.7, "pt"),
    legend.spacing.x = unit(2.5, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.box.margin = margin(t = 1),
    plot.margin = margin(4, 4, 1, 3)
  )

rep_data <- read.csv(manova_file) |>
  filter(
    is.finite(shbeta_between),
    is.finite(hric_manova_F),
    is.finite(hric_ss_between)
  ) |>
  mutate(
    signal_mode = factor(
      signal_mode,
      levels = c("Abundance", "Both", "Prevalence"),
      labels = c("Abundance", "Joint", "Prevalence")
    )
  )
required_manova_columns <- c(
  "rep_id", "k", "signal_mode", "beta", "shbeta_between",
  "hric_manova_F", "hric_ss_between"
)
if (length(setdiff(required_manova_columns, names(rep_data))) > 0L) {
  stop("The MANOVA input does not contain the required columns.")
}

mode_colors <- c(
  Abundance = "#0072B2",
  Joint = "#009E73",
  Prevalence = "#D55E00"
)

fit_predictions <- function(data, response, n_grid = 301L) {
  split(data, data$signal_mode) |>
    lapply(function(mode_data) {
      model <- stats::lm(
        stats::reformulate("shbeta_between", response = response),
        data = mode_data
      )
      grid <- data.frame(
        shbeta_between = seq(
          min(mode_data$shbeta_between),
          max(mode_data$shbeta_between),
          length.out = n_grid
        )
      )
      fitted <- stats::predict(model, newdata = grid, se.fit = TRUE)
      grid$fit <- as.numeric(fitted$fit)
      grid$lower <- grid$fit - 1.96 * as.numeric(fitted$se.fit)
      grid$upper <- grid$fit + 1.96 * as.numeric(fitted$se.fit)
      grid$signal_mode <- unique(mode_data$signal_mode)
      grid$r_squared <- summary(model)$r.squared
      grid
    }) |>
    bind_rows() |>
    mutate(signal_mode = factor(signal_mode, levels = levels(rep_data$signal_mode)))
}

f_predictions <- fit_predictions(rep_data, "hric_manova_F")
r2_labels <- f_predictions |>
  group_by(signal_mode) |>
  summarise(
    label = sprintf("italic(R)^2 == %.3f", first(r_squared)),
    .groups = "drop"
  )

base_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.18),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "#333333", linewidth = 0.28),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.24),
    axis.ticks.length = unit(1.8, "pt"),
    axis.title = element_text(size = 7, colour = "#111111"),
    axis.text = element_text(size = 6.2, colour = "#222222"),
    legend.position = "none",
    plot.margin = margin(3, 4, 2, 3)
  )

f_plot <- ggplot() +
  geom_ribbon(
    data = f_predictions,
    aes(x = shbeta_between, ymin = lower, ymax = upper, fill = signal_mode),
    alpha = 0.12,
    linewidth = 0
  ) +
  geom_point(
    data = rep_data,
    aes(x = shbeta_between, y = hric_manova_F, colour = signal_mode),
    alpha = 0.24,
    size = 0.46,
    stroke = 0
  ) +
  geom_line(
    data = f_predictions,
    aes(x = shbeta_between, y = fit, colour = signal_mode),
    linewidth = 0.40,
    lineend = "round"
  ) +
  geom_text(
    data = r2_labels,
    aes(x = Inf, y = Inf, label = label),
    parse = TRUE,
    hjust = 1.12,
    vjust = 1.35,
    size = 2.05,
    family = "Arial",
    colour = "#333333"
  ) +
  facet_wrap(~signal_mode, nrow = 1) +
  scale_colour_manual(values = mode_colors) +
  scale_fill_manual(values = mode_colors) +
  scale_x_continuous(
    breaks = c(0.004, 0.007, 0.010),
    labels = c("0.004", "0.007", "0.010"),
    expand = expansion(mult = c(0.03, 0.04))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.2),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  labs(x = "Normalized between-group turnover", y = "HRIC MANOVA pseudo-F") +
  base_theme +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 6.1,
      face = "bold",
      colour = "#222222",
      margin = margin(2.2, 1.5, 2.2, 1.5)
    ),
    axis.title = element_text(size = 6.6),
    axis.text = element_text(size = 5.6),
    panel.spacing.x = unit(6.5, "pt"),
    plot.margin = margin(3, 5, 3, 4)
  )

# The identity in panel b is checked against replication-level values rather
# than estimated from rendered graphics.
theory_check <- read.csv(theory_file)
max_ss_error <- theory_check$value[
  theory_check$check == "Maximum |SS_between error|"
]
if (length(max_ss_error) != 1L || !is.finite(max_ss_error)) {
  stop("The theoretical SS-between check is unavailable.")
}

theory_slope <- stats::median(
  rep_data$hric_ss_between / rep_data$shbeta_between
)
ss_range <- range(rep_data$shbeta_between)
ss_line <- data.frame(
  shbeta_between = ss_range,
  hric_ss_between = theory_slope * ss_range
)

ss_plot <- ggplot() +
  geom_point(
    data = rep_data,
    aes(x = shbeta_between, y = hric_ss_between),
    colour = "#3E7396",
    alpha = 0.32,
    size = 0.50,
    stroke = 0
  ) +
  geom_line(
    data = ss_line,
    aes(x = shbeta_between, y = hric_ss_between),
    colour = "#111111",
    linewidth = 0.28,
    linetype = "22",
    lineend = "round"
  ) +
  annotate(
    "text",
    x = ss_range[1] + 0.05 * diff(ss_range),
    y = max(rep_data$hric_ss_between) - 0.06 * diff(range(rep_data$hric_ss_between)),
    label = sprintf("SS[between] == %.1f %%.%% delta[SH]^\"*\"", theory_slope),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 2.35,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = ss_range[1] + 0.05 * diff(ss_range),
    y = max(rep_data$hric_ss_between) - 0.20 * diff(range(rep_data$hric_ss_between)),
    label = sprintf("max |error| = %.1e", max_ss_error),
    hjust = 0,
    vjust = 1,
    size = 1.95,
    family = "Arial"
  ) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.001),
    expand = expansion(mult = c(0.03, 0.04))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.2),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  labs(x = "Normalized between-group turnover", y = expression(SS[between])) +
  base_theme +
  theme(
    axis.title = element_text(size = 6.6),
    axis.text = element_text(size = 5.8),
    plot.margin = margin(3, 4, 3, 5)
  )

equivalence_row <- (f_plot | ss_plot) + plot_layout(widths = c(1.8, 1))
integrated_plot <- equivalence_row /
  power_plot +
  plot_layout(heights = c(1, 1.95)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(
      family = "Arial",
      face = "bold",
      size = 8,
      colour = "#111111"
    ),
    plot.tag.position = c(0, 1)
  )

output_file <- file.path(
  result_dir,
  "figure_integrated_simulation_hric_turnover.pdf"
)
ggsave(
  output_file,
  plot = integrated_plot,
  device = grDevices::cairo_pdf,
  width = 180,
  height = 170,
  units = "mm",
  bg = "white"
)

write.csv(
  power_data |>
    transmute(
      panel = as.character(signal_label),
      k,
      beta,
      method = unname(method_labels[as.character(method)]),
      empirical_power = power,
      n_replicates = n_nonmissing
    ),
  file.path(result_dir, "power_values.csv"),
  row.names = FALSE
)
write.csv(
  rep_data |>
    transmute(
      replicate = rep_id,
      k,
      signal_mode,
      beta,
      normalized_between_group_turnover = shbeta_between,
      manova_pseudo_F = hric_manova_F,
      ss_between = hric_ss_between
    ),
  file.path(result_dir, "turnover_manova_values.csv"),
  row.names = FALSE
)

stopifnot(
  nrow(power_data) == 2880L,
  nrow(rep_data) == 7875L,
  max_ss_error < 1e-10,
  all(rep_data |>
    count(k, signal_mode, beta) |>
    pull(n) == 125L),
  file.exists(output_file)
)

message("Saved ", output_file)
