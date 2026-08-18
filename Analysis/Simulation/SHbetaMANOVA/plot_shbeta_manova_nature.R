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
  library(scales)
})

script_dir_override <- Sys.getenv("SIM_MANOVA_DIR", unset = "")
script_dir <- if (nzchar(script_dir_override)) {
  normalizePath(script_dir_override)
} else {
  get_script_dir()
}

source_file <- file.path(
  script_dir,
  "source_data",
  "table_replication_level_hric_manova.csv"
)
theory_file <- file.path(
  script_dir,
  "source_data",
  "table_hric_shbeta_theory_check.csv"
)
derived_dir <- file.path(script_dir, "derived_data")
figure_dir <- file.path(script_dir, "figures")

dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_file)) {
  stop("Replication-level source table not found: ", source_file)
}
if (!file.exists(theory_file)) {
  stop("Theory-check source table not found: ", theory_file)
}

rep_data <- read.csv(source_file, stringsAsFactors = FALSE) |>
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

required_columns <- c(
  "rep_id",
  "k",
  "signal_mode",
  "beta",
  "shbeta_between",
  "hric_manova_F",
  "hric_ss_between"
)
missing_columns <- setdiff(required_columns, names(rep_data))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

mode_colors <- c(
  "Abundance" = "#0072B2",
  "Joint" = "#009E73",
  "Prevalence" = "#D55E00"
)

fit_predictions <- function(dat, response, n_grid = 301L) {
  split(dat, dat$signal_mode) |>
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
    mutate(
      signal_mode = factor(
        signal_mode,
        levels = levels(rep_data$signal_mode)
      )
    )
}

f_predictions <- fit_predictions(rep_data, "hric_manova_F")
ss_predictions <- fit_predictions(rep_data, "hric_ss_between")
r2_labels <- f_predictions |>
  group_by(signal_mode) |>
  summarise(
    label = sprintf("italic(R)^2 == %.3f", first(r_squared)),
    .groups = "drop"
  )

write.csv(
  f_predictions,
  file.path(derived_dir, "manova_F_vs_shbeta_fitted_curves.csv"),
  row.names = FALSE
)
write.csv(
  ss_predictions,
  file.path(derived_dir, "ss_between_vs_shbeta_fitted_curves.csv"),
  row.names = FALSE
)

base_theme <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid.major.y = element_line(
      colour = "#E6E6E6",
      linewidth = 0.18
    ),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "#333333", linewidth = 0.28),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.24),
    axis.ticks.length = unit(1.8, "pt"),
    axis.title = element_text(size = 7, colour = "#111111"),
    axis.text = element_text(size = 6.2, colour = "#222222"),
    legend.position = "none",
    plot.margin = margin(3, 4, 2, 3)
  )

plot_f_vs_shbeta <- ggplot() +
  geom_ribbon(
    data = f_predictions,
    aes(
      x = shbeta_between,
      ymin = lower,
      ymax = upper,
      fill = signal_mode
    ),
    alpha = 0.12,
    linewidth = 0
  ) +
  geom_point(
    data = rep_data,
    aes(
      x = shbeta_between,
      y = hric_manova_F,
      colour = signal_mode
    ),
    alpha = 0.24,
    size = 0.46,
    stroke = 0
  ) +
  geom_line(
    data = f_predictions,
    aes(
      x = shbeta_between,
      y = fit,
      colour = signal_mode
    ),
    linewidth = 0.40,
    lineend = "round"
  ) +
  geom_text(
    data = r2_labels,
    aes(
      x = Inf,
      y = Inf,
      label = label
    ),
    parse = TRUE,
    hjust = 1.12,
    vjust = 1.35,
    size = 2.05,
    family = "Arial",
    colour = "#333333"
  ) +
  facet_wrap(
    ~signal_mode,
    nrow = 1
  ) +
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
  labs(
    x = "Normalized between-group turnover",
    y = "HRIC MANOVA pseudo-F"
  ) +
  base_theme +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 6.2,
      face = "bold",
      colour = "#222222",
      margin = margin(2.2, 1.5, 2.2, 1.5)
    ),
    panel.spacing.x = unit(6.5, "pt")
  )

theory_check <- read.csv(theory_file, stringsAsFactors = FALSE)
max_ss_error <- theory_check$value[
  theory_check$check == "Maximum |SS_between error|"
]
if (length(max_ss_error) != 1L || !is.finite(max_ss_error)) {
  stop("Could not recover the maximum SS_between error.")
}

theory_slope <- stats::median(
  rep_data$hric_ss_between / rep_data$shbeta_between
)
ss_range <- range(rep_data$shbeta_between)
ss_line <- data.frame(
  shbeta_between = ss_range,
  hric_ss_between = theory_slope * ss_range
)

plot_ss_vs_shbeta <- ggplot() +
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
    y = max(rep_data$hric_ss_between) - 0.06 *
      diff(range(rep_data$hric_ss_between)),
    label = sprintf(
      "SS[between] == %.1f %%.%% delta[SH]^\"*\"",
      theory_slope
    ),
    parse = TRUE,
    hjust = 0,
    vjust = 1,
    size = 2.35,
    family = "Arial"
  ) +
  annotate(
    "text",
    x = ss_range[1] + 0.05 * diff(ss_range),
    y = max(rep_data$hric_ss_between) - 0.20 *
      diff(range(rep_data$hric_ss_between)),
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
  labs(
    x = "Normalized between-group turnover",
    y = expression(SS[between])
  ) +
  base_theme

save_plot <- function(plot, stem, width_mm = 90, height_mm = 62) {
  ggsave(
    file.path(figure_dir, paste0(stem, ".pdf")),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width_mm,
    height = height_mm,
    units = "mm",
    bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".png")),
    plot = plot,
    device = ragg::agg_png,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    file.path(figure_dir, paste0(stem, ".tiff")),
    plot = plot,
    device = ragg::agg_tiff,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
}

save_plot(
  plot_f_vs_shbeta,
  "fig_manova_F_vs_shbeta_between_nature",
  width_mm = 125
)
save_plot(
  plot_ss_vs_shbeta,
  "fig_ss_between_vs_shbeta_between_nature",
  width_mm = 78
)

message("Saved SHbeta-MANOVA figures to: ", figure_dir)
