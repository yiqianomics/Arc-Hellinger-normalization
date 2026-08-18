#!/usr/bin/env Rscript

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(ragg)
})

script_dir <- get_script_dir()
simulation_dir <- dirname(script_dir)
power_dir <- file.path(simulation_dir, "BetaDiversity")
manova_dir <- file.path(simulation_dir, "SHbetaMANOVA")

power_script <- file.path(
  power_dir,
  "plot_max_n_otu_power_nature.R"
)
manova_script <- file.path(
  manova_dir,
  "plot_shbeta_manova_nature.R"
)

if (!file.exists(power_script)) {
  stop("Power plotting script not found: ", power_script)
}
if (!file.exists(manova_script)) {
  stop("SHbeta-MANOVA plotting script not found: ", manova_script)
}

old_power_dir <- Sys.getenv("SIM_POWER_DIR", unset = NA_character_)
old_manova_dir <- Sys.getenv("SIM_MANOVA_DIR", unset = NA_character_)

on.exit({
  if (is.na(old_power_dir)) {
    Sys.unsetenv("SIM_POWER_DIR")
  } else {
    Sys.setenv(SIM_POWER_DIR = old_power_dir)
  }
  if (is.na(old_manova_dir)) {
    Sys.unsetenv("SIM_MANOVA_DIR")
  } else {
    Sys.setenv(SIM_MANOVA_DIR = old_manova_dir)
  }
}, add = TRUE)

power_env <- new.env(parent = globalenv())
Sys.setenv(SIM_POWER_DIR = power_dir)
sys.source(power_script, envir = power_env)

manova_env <- new.env(parent = globalenv())
Sys.setenv(SIM_MANOVA_DIR = manova_dir)
sys.source(manova_script, envir = manova_env)

power_panel <- power_env$plot_base +
  guides(
    colour = guide_legend(
      nrow = 3,
      byrow = TRUE,
      override.aes = list(linewidth = 0.45, alpha = 1)
    )
  ) +
  theme(
    axis.title = element_text(size = 6.5),
    axis.text = element_text(size = 5.6),
    strip.text.x = element_text(size = 6.4),
    strip.text.y.left = element_text(size = 6.4),
    legend.text = element_text(size = 5.25),
    legend.key.width = unit(14, "pt"),
    legend.key.height = unit(6.7, "pt"),
    legend.spacing.x = unit(2.5, "pt"),
    legend.box.margin = margin(t = 1),
    panel.spacing.x = unit(3.5, "pt"),
    panel.spacing.y = unit(3.5, "pt"),
    plot.margin = margin(4, 4, 1, 3)
  )

f_panel <- manova_env$plot_f_vs_shbeta +
  theme(
    axis.title = element_text(size = 6.6),
    axis.text = element_text(size = 5.6),
    strip.text = element_text(size = 6.1),
    panel.spacing.x = unit(6.5, "pt"),
    plot.margin = margin(3, 5, 3, 4)
  )

ss_panel <- manova_env$plot_ss_vs_shbeta +
  theme(
    axis.title = element_text(size = 6.6),
    axis.text = element_text(size = 5.8),
    plot.margin = margin(3, 4, 3, 5)
  )

equivalence_row <- f_panel | ss_panel
equivalence_row <- equivalence_row + plot_layout(widths = c(1.8, 1))

integrated_plot <- equivalence_row /
  power_panel +
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

output_stem <- file.path(
  script_dir,
  "figure_integrated_simulation_hric_turnover"
)

ggsave(
  paste0(output_stem, ".pdf"),
  plot = integrated_plot,
  device = grDevices::cairo_pdf,
  width = 180,
  height = 170,
  units = "mm",
  bg = "white"
)
ggsave(
  paste0(output_stem, ".png"),
  plot = integrated_plot,
  device = ragg::agg_png,
  width = 180,
  height = 170,
  units = "mm",
  dpi = 600,
  bg = "white"
)
ggsave(
  paste0(output_stem, ".tiff"),
  plot = integrated_plot,
  device = ragg::agg_tiff,
  width = 180,
  height = 170,
  units = "mm",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

message("Saved integrated simulation figure: ", output_stem)
