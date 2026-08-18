#!/usr/bin/env Rscript

## ============================================================
## Summarize SparseDOSSA2 SHbeta and HRIC MANOVA simulations
##
## Main quantities:
##   1. SHbeta_between
##   2. HRIC SS_between
##   3. HRIC MANOVA pseudo-F
##   4. HRIC between-group R^2
##   5. SHbeta total / within / between decomposition
##   6. Exact relationship:
##
##        SS_between = n * A_p^2 * SHbeta_between
##
##      where
##
##        A_p = asin(sqrt(1 - 1/p)).
##
## Summary strategy:
##   - First average across signal clusters within each replication.
##   - Then summarize across independent replications.
##
## Output:
##   - PDF figures only
##   - CSV summary tables
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## ------------------------------------------------------------
## User/HPC settings
## ------------------------------------------------------------

OUTDIR <- Sys.getenv(
  "OUTDIR",
  "/home/zhang.16383/AHC/beta_div_sparse12_shbeta_hric_manova_out"
)

RAW_DIR <- file.path(OUTDIR, "raw")
SUMMARY_DIR <- file.path(OUTDIR, "summary_hric_manova")

K_TARGET <- as.integer(Sys.getenv("K_COMM", "12"))

dir.create(
  SUMMARY_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Raw-data directory:\n  ", RAW_DIR, "\n", sep = "")
cat("Summary directory:\n  ", SUMMARY_DIR, "\n", sep = "")
cat("Target K:\n  ", K_TARGET, "\n\n", sep = "")

## ------------------------------------------------------------
## Locate result files
## ------------------------------------------------------------

csv_files <- list.files(
  RAW_DIR,
  pattern = "^beta_div_sparse12_shbeta_hric_manova_rep[0-9]+\\.csv$",
  full.names = TRUE
)

## Fallback in case the filename prefix differs slightly
if (length(csv_files) == 0L) {
  csv_files <- list.files(
    RAW_DIR,
    pattern = "rep[0-9]+\\.csv$",
    full.names = TRUE
  )

  ## Do not accidentally read the long-format MANOVA tables
  csv_files <- csv_files[
    !grepl("manova_table", basename(csv_files), ignore.case = TRUE)
  ]
}

if (length(csv_files) == 0L) {
  stop(
    "No simulation CSV files were found in:\n  ",
    RAW_DIR
  )
}

cat("Number of raw CSV files found: ", length(csv_files), "\n", sep = "")
cat("First file:\n  ", csv_files[1L], "\n\n", sep = "")

## ------------------------------------------------------------
## Read and combine
## ------------------------------------------------------------

raw_list <- lapply(csv_files, function(f) {
  out <- tryCatch(
    fread(f, showProgress = FALSE),
    error = function(e) {
      warning(
        "Could not read file: ",
        f,
        "\nReason: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (!is.null(out)) {
    out[, source_file := basename(f)]
  }

  out
})

raw_list <- Filter(Negate(is.null), raw_list)

if (length(raw_list) == 0L) {
  stop("None of the simulation CSV files could be read.")
}

dat <- rbindlist(
  raw_list,
  use.names = TRUE,
  fill = TRUE
)

cat("Combined rows: ", nrow(dat), "\n", sep = "")
cat("Combined columns: ", ncol(dat), "\n", sep = "")

## ------------------------------------------------------------
## Required columns
## ------------------------------------------------------------

required_columns <- c(
  "rep_id",
  "k",
  "signal_mode",
  "beta",
  "n_samples",
  "shbeta_total",
  "shbeta_weighted_within",
  "shbeta_between",
  "hric_ss_between",
  "hric_ss_within",
  "hric_ss_total",
  "hric_df_between",
  "hric_df_within",
  "hric_ms_between",
  "hric_ms_within",
  "hric_manova_F",
  "hric_r2_between",
  "hric_coordinate_dimension"
)

missing_columns <- setdiff(required_columns, names(dat))

if (length(missing_columns) > 0L) {
  stop(
    "The combined result is missing the following required columns:\n  ",
    paste(missing_columns, collapse = "\n  ")
  )
}

## ------------------------------------------------------------
## Basic cleaning
## ------------------------------------------------------------

numeric_columns <- intersect(
  c(
    "rep_id",
    "k",
    "beta",
    "n_samples",
    "shbeta_total",
    "shbeta_group0",
    "shbeta_group1",
    "shbeta_abs_diff",
    "shbeta_weighted_within",
    "shbeta_between",
    "hric_ss_between",
    "hric_ss_within",
    "hric_ss_total",
    "hric_df_between",
    "hric_df_within",
    "hric_df_total",
    "hric_ms_between",
    "hric_ms_within",
    "hric_manova_F",
    "hric_r2_between",
    "hric_coordinate_dimension",
    "p_SHbeta"
  ),
  names(dat)
)

for (v in numeric_columns) {
  set(dat, j = v, value = as.numeric(dat[[v]]))
}

dat[, signal_mode := factor(
  signal_mode,
  levels = c("abundance", "prevalence", "both")
)]

dat <- dat[
  k == K_TARGET &
    is.finite(beta) &
    !is.na(signal_mode)
]

if (nrow(dat) == 0L) {
  stop(
    "No valid rows remain after filtering to k = ",
    K_TARGET,
    "."
  )
}

cat("\nReplications represented:\n")
print(sort(unique(dat$rep_id)))

cat("\nNumber of unique replications: ",
    uniqueN(dat$rep_id),
    "\n",
    sep = "")

cat("\nRows by signal mode:\n")
print(dat[, .N, by = signal_mode])

cat("\nRows by beta:\n")
print(dat[, .N, by = beta][order(beta)])

## ------------------------------------------------------------
## Check duplicate setting rows
## ------------------------------------------------------------

setting_id_columns <- intersect(
  c(
    "rep_id",
    "k",
    "signal_cluster",
    "signal_mode",
    "beta"
  ),
  names(dat)
)

duplicate_check <- dat[
  ,
  .N,
  by = setting_id_columns
][N > 1L]

if (nrow(duplicate_check) > 0L) {
  warning(
    "Duplicate setting rows were detected. ",
    "The first-stage aggregation will average them."
  )
  print(head(duplicate_check, 20L))
}

## ============================================================
## Theory checks
## ============================================================

## A_p = asin(sqrt(1 - 1/p))
dat[
  ,
  hric_Ap := asin(
    sqrt(
      1 - 1 / hric_coordinate_dimension
    )
  )
]

## Reconstruct SS components from SHbeta decomposition
dat[
  ,
  ss_between_from_shbeta :=
    n_samples *
    hric_Ap^2 *
    shbeta_between
]

dat[
  ,
  ss_total_from_shbeta :=
    n_samples *
    hric_Ap^2 *
    shbeta_total
]

dat[
  ,
  ss_within_from_shbeta :=
    n_samples *
    hric_Ap^2 *
    shbeta_weighted_within
]

dat[
  ,
  error_ss_between :=
    hric_ss_between -
    ss_between_from_shbeta
]

dat[
  ,
  error_ss_total :=
    hric_ss_total -
    ss_total_from_shbeta
]

dat[
  ,
  error_ss_within :=
    hric_ss_within -
    ss_within_from_shbeta
]

## Reconstruct F from SHbeta between/within components
dat[
  ,
  f_from_shbeta :=
    (shbeta_between / hric_df_between) /
    (shbeta_weighted_within / hric_df_within)
]

dat[
  ,
  error_f_from_shbeta :=
    hric_manova_F -
    f_from_shbeta
]

safe_max_abs <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  max(abs(x))
}

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)

  if (sum(keep) < 3L) {
    return(NA_real_)
  }

  if (sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(NA_real_)
  }

  cor(x[keep], y[keep])
}

theory_check <- data.table(
  quantity = c(
    "SS_between",
    "SS_within",
    "SS_total",
    "MANOVA_F",
    "cor(SHbeta_between, SS_between)"
  ),
  value = c(
    safe_max_abs(dat$error_ss_between),
    safe_max_abs(dat$error_ss_within),
    safe_max_abs(dat$error_ss_total),
    safe_max_abs(dat$error_f_from_shbeta),
    safe_cor(
      dat$shbeta_between,
      dat$hric_ss_between
    )
  )
)

cat("\nTheory consistency checks:\n")
print(theory_check)

fwrite(
  theory_check,
  file.path(
    SUMMARY_DIR,
    "table_hric_shbeta_theory_check.csv"
  )
)

## ============================================================
## Stage 1:
## Average over signal clusters within each replication
##
## This makes the independent unit the simulation replication,
## rather than treating all clusters as independent observations.
## ============================================================

metrics <- c(
  "shbeta_total",
  "shbeta_weighted_within",
  "shbeta_between",
  "hric_ss_between",
  "hric_ss_within",
  "hric_ss_total",
  "hric_manova_F",
  "hric_r2_between"
)

if ("shbeta_abs_diff" %in% names(dat)) {
  metrics <- c(metrics, "shbeta_abs_diff")
}

if ("p_SHbeta" %in% names(dat)) {
  metrics <- c(metrics, "p_SHbeta")
}

mean_finite <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(x)
}

rep_level <- dat[
  ,
  c(
    list(
      n_cluster_rows = .N,
      n_distinct_clusters = if (
        "signal_cluster" %in% names(dat)
      ) {
        uniqueN(signal_cluster)
      } else {
        NA_integer_
      }
    ),
    lapply(.SD, mean_finite)
  ),
  by = .(
    rep_id,
    k,
    signal_mode,
    beta
  ),
  .SDcols = metrics
]

setorder(
  rep_level,
  signal_mode,
  beta,
  rep_id
)

fwrite(
  rep_level,
  file.path(
    SUMMARY_DIR,
    "table_replication_level_hric_manova.csv"
  )
)

## ============================================================
## Stage 2:
## Summarize across independent replications
## ============================================================

summarize_vector <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (n == 0L) {
    return(list(
      n_rep = 0L,
      mean = NA_real_,
      sd = NA_real_,
      mcse = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      median = NA_real_,
      q25 = NA_real_,
      q75 = NA_real_
    ))
  }

  x_mean <- mean(x)
  x_sd <- if (n > 1L) sd(x) else NA_real_
  x_mcse <- if (n > 1L) x_sd / sqrt(n) else NA_real_

  list(
    n_rep = n,
    mean = x_mean,
    sd = x_sd,
    mcse = x_mcse,
    lower = if (is.finite(x_mcse)) {
      x_mean - 1.96 * x_mcse
    } else {
      NA_real_
    },
    upper = if (is.finite(x_mcse)) {
      x_mean + 1.96 * x_mcse
    } else {
      NA_real_
    },
    median = median(x),
    q25 = as.numeric(quantile(
      x,
      0.25,
      names = FALSE,
      type = 7
    )),
    q75 = as.numeric(quantile(
      x,
      0.75,
      names = FALSE,
      type = 7
    ))
  )
}

summary_long <- rbindlist(
  lapply(metrics, function(metric_name) {
    rep_level[
      ,
      {
        s <- summarize_vector(get(metric_name))

        data.table(
          metric = metric_name,
          n_rep = s$n_rep,
          mean = s$mean,
          sd = s$sd,
          mcse = s$mcse,
          lower = s$lower,
          upper = s$upper,
          median = s$median,
          q25 = s$q25,
          q75 = s$q75
        )
      },
      by = .(
        k,
        signal_mode,
        beta
      )
    ]
  }),
  use.names = TRUE,
  fill = TRUE
)

setorder(
  summary_long,
  metric,
  signal_mode,
  beta
)

fwrite(
  summary_long,
  file.path(
    SUMMARY_DIR,
    "table_summary_hric_manova_long.csv"
  )
)

summary_wide <- dcast(
  summary_long,
  k + signal_mode + beta ~ metric,
  value.var = c(
    "n_rep",
    "mean",
    "sd",
    "mcse",
    "lower",
    "upper",
    "median",
    "q25",
    "q75"
  )
)

fwrite(
  summary_wide,
  file.path(
    SUMMARY_DIR,
    "table_summary_hric_manova_wide.csv"
  )
)

cat("\nMain summary preview:\n")
print(
  summary_long[
    metric %in% c(
      "shbeta_between",
      "hric_ss_between",
      "hric_manova_F",
      "hric_r2_between"
    )
  ][
    order(metric, signal_mode, beta)
  ]
)

## ============================================================
## Plotting helpers
## ============================================================

base_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    axis.title = element_text(
      face = "plain"
    )
  )

get_metric_summary <- function(metric_name) {
  out <- summary_long[
    metric == metric_name
  ]

  if (nrow(out) == 0L) {
    stop(
      "No summary rows found for metric: ",
      metric_name
    )
  }

  out
}

make_trend_plot <- function(
    metric_name,
    y_label,
    filename,
    force_zero = FALSE,
    log_y = FALSE
) {
  plot_dat <- get_metric_summary(metric_name)

  p <- ggplot(
    plot_dat,
    aes(
      x = beta,
      y = mean,
      group = 1
    )
  ) +
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper
      ),
      alpha = 0.18,
      na.rm = TRUE
    ) +
    geom_line(
      linewidth = 0.75,
      na.rm = TRUE
    ) +
    geom_point(
      size = 1.8,
      na.rm = TRUE
    ) +
    facet_wrap(
      ~ signal_mode,
      nrow = 1,
      scales = "free_y"
    ) +
    labs(
      x = "SparseDOSSA2 effect size",
      y = y_label
    ) +
    base_theme

  if (force_zero && !log_y) {
    p <- p + expand_limits(y = 0)
  }

  if (log_y) {
    p <- p +
      scale_y_log10() +
      annotation_logticks(sides = "l")
  }

  ggsave(
    filename = file.path(
      SUMMARY_DIR,
      filename
    ),
    plot = p,
    width = 10,
    height = 4,
    device = cairo_pdf
  )

  invisible(p)
}

## ============================================================
## Figure 1: between-group SHbeta trend
## ============================================================

make_trend_plot(
  metric_name = "shbeta_between",
  y_label = expression(
    "Between-group " * beta[SH]
  ),
  filename = "fig_shbeta_between_trend.pdf",
  force_zero = TRUE
)

## ============================================================
## Figure 2: HRIC SS_between trend
## ============================================================

make_trend_plot(
  metric_name = "hric_ss_between",
  y_label = expression(
    SS[between]
  ),
  filename = "fig_hric_ss_between_trend.pdf",
  force_zero = TRUE
)

## ============================================================
## Figure 3: HRIC MANOVA pseudo-F trend
## ============================================================

make_trend_plot(
  metric_name = "hric_manova_F",
  y_label = "HRIC MANOVA pseudo-F",
  filename = "fig_hric_manova_F_trend.pdf",
  force_zero = TRUE
)

## ============================================================
## Figure 4: HRIC between-group R-squared trend
## ============================================================

make_trend_plot(
  metric_name = "hric_r2_between",
  y_label = expression(
    R[between]^2
  ),
  filename = "fig_hric_r2_between_trend.pdf",
  force_zero = TRUE
)

## ============================================================
## Figure 5:
## SHbeta decomposition:
## total = weighted within + between
## ============================================================

decomposition_metrics <- c(
  shbeta_total = "Total",
  shbeta_weighted_within = "Weighted within",
  shbeta_between = "Between"
)

decomposition_dat <- summary_long[
  metric %in% names(decomposition_metrics)
]

decomposition_dat[
  ,
  component := factor(
    decomposition_metrics[metric],
    levels = c(
      "Total",
      "Weighted within",
      "Between"
    )
  )
]

p_decomposition <- ggplot(
  decomposition_dat,
  aes(
    x = beta,
    y = mean,
    linetype = component,
    shape = component,
    group = component
  )
) +
  geom_line(
    linewidth = 0.75,
    na.rm = TRUE
  ) +
  geom_point(
    size = 1.7,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ signal_mode,
    nrow = 1,
    scales = "free_y"
  ) +
  labs(
    x = "SparseDOSSA2 effect size",
    y = expression(beta[SH]),
    linetype = NULL,
    shape = NULL
  ) +
  base_theme

ggsave(
  filename = file.path(
    SUMMARY_DIR,
    "fig_shbeta_decomposition_trend.pdf"
  ),
  plot = p_decomposition,
  width = 10,
  height = 4,
  device = cairo_pdf
)

## ============================================================
## Figure 6:
## MANOVA F versus SHbeta_between
##
## This shows how F standardizes between-group beta diversity
## by within-group dispersion.
## ============================================================

relationship_rep <- rep_level[
  is.finite(shbeta_between) &
    is.finite(hric_manova_F)
]

p_f_vs_beta <- ggplot(
  relationship_rep,
  aes(
    x = shbeta_between,
    y = hric_manova_F,
    group = interaction(
      signal_mode,
      beta
    )
  )
) +
  geom_point(
    alpha = 0.25,
    size = 1.1
  ) +
  geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ signal_mode,
    nrow = 1,
    scales = "free"
  ) +
  labs(
    x = expression(
      "Between-group " * beta[SH]
    ),
    y = "HRIC MANOVA pseudo-F"
  ) +
  base_theme

ggsave(
  filename = file.path(
    SUMMARY_DIR,
    "fig_manova_F_vs_shbeta_between.pdf"
  ),
  plot = p_f_vs_beta,
  width = 10,
  height = 4,
  device = cairo_pdf
)

## ============================================================
## Figure 7:
## Exact row-level relationship:
##
## SS_between = n * A_p^2 * SHbeta_between
##
## The theoretical relationship is shown as y = x after
## converting SHbeta_between to the SS scale.
## ============================================================

relationship_raw <- dat[
  is.finite(ss_between_from_shbeta) &
    is.finite(hric_ss_between)
]

p_ss_identity <- ggplot(
  relationship_raw,
  aes(
    x = ss_between_from_shbeta,
    y = hric_ss_between
  )
) +
  geom_point(
    alpha = 0.20,
    size = 0.9
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2,
    linewidth = 0.7
  ) +
  facet_wrap(
    ~ signal_mode,
    nrow = 1,
    scales = "free"
  ) +
  coord_equal() +
  labs(
    x = expression(
      n * A[p]^2 * beta[SH][", between"]
    ),
    y = expression(
      SS[between]
    )
  ) +
  base_theme

ggsave(
  filename = file.path(
    SUMMARY_DIR,
    "fig_ss_between_vs_scaled_shbeta_between.pdf"
  ),
  plot = p_ss_identity,
  width = 10,
  height = 4,
  device = cairo_pdf
)

## ============================================================
## Figure 8:
## Direct unscaled relationship between SHbeta_between and
## SS_between.
##
## Since p and n are fixed, this should also be exactly linear.
## ============================================================

p_ss_unscaled <- ggplot(
  relationship_raw,
  aes(
    x = shbeta_between,
    y = hric_ss_between
  )
) +
  geom_point(
    alpha = 0.20,
    size = 0.9
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ signal_mode,
    nrow = 1,
    scales = "free"
  ) +
  labs(
    x = expression(
      "Between-group " * beta[SH]
    ),
    y = expression(
      SS[between]
    )
  ) +
  base_theme

ggsave(
  filename = file.path(
    SUMMARY_DIR,
    "fig_ss_between_vs_shbeta_between.pdf"
  ),
  plot = p_ss_unscaled,
  width = 10,
  height = 4,
  device = cairo_pdf
)

## ============================================================
## Optional figure:
## group-specific SHbeta values, if available
## ============================================================

if (
  all(
    c(
      "shbeta_group0",
      "shbeta_group1"
    ) %in% names(dat)
  )
) {
  group_rep <- dat[
    ,
    .(
      shbeta_group0 = mean_finite(
        shbeta_group0
      ),
      shbeta_group1 = mean_finite(
        shbeta_group1
      )
    ),
    by = .(
      rep_id,
      k,
      signal_mode,
      beta
    )
  ]

  group_long <- melt(
    group_rep,
    id.vars = c(
      "rep_id",
      "k",
      "signal_mode",
      "beta"
    ),
    measure.vars = c(
      "shbeta_group0",
      "shbeta_group1"
    ),
    variable.name = "group",
    value.name = "shbeta"
  )

  group_long[
    ,
    group := factor(
      group,
      levels = c(
        "shbeta_group0",
        "shbeta_group1"
      ),
      labels = c(
        "Group 0",
        "Group 1"
      )
    )
  ]

  group_summary <- group_long[
    ,
    {
      s <- summarize_vector(shbeta)

      .(
        n_rep = s$n_rep,
        mean = s$mean,
        sd = s$sd,
        mcse = s$mcse,
        lower = s$lower,
        upper = s$upper
      )
    },
    by = .(
      signal_mode,
      beta,
      group
    )
  ]

  fwrite(
    group_summary,
    file.path(
      SUMMARY_DIR,
      "table_group_specific_shbeta_summary.csv"
    )
  )

  p_group_shbeta <- ggplot(
    group_summary,
    aes(
      x = beta,
      y = mean,
      linetype = group,
      shape = group,
      group = group
    )
  ) +
    geom_line(
      linewidth = 0.75,
      na.rm = TRUE
    ) +
    geom_point(
      size = 1.7,
      na.rm = TRUE
    ) +
    facet_wrap(
      ~ signal_mode,
      nrow = 1,
      scales = "free_y"
    ) +
    labs(
      x = "SparseDOSSA2 effect size",
      y = expression(beta[SH]),
      linetype = NULL,
      shape = NULL
    ) +
    base_theme

  ggsave(
    filename = file.path(
      SUMMARY_DIR,
      "fig_group_specific_shbeta_trend.pdf"
    ),
    plot = p_group_shbeta,
    width = 10,
    height = 4,
    device = cairo_pdf
  )
}

## ============================================================
## Compact comparison table
## ============================================================

main_metrics <- c(
  "shbeta_between",
  "hric_ss_between",
  "hric_manova_F",
  "hric_r2_between"
)

main_summary <- summary_long[
  metric %in% main_metrics
]

main_summary_wide <- dcast(
  main_summary,
  k + signal_mode + beta + n_rep ~ metric,
  value.var = c(
    "mean",
    "mcse",
    "lower",
    "upper"
  )
)

setorder(
  main_summary_wide,
  signal_mode,
  beta
)

fwrite(
  main_summary_wide,
  file.path(
    SUMMARY_DIR,
    "table_main_shbeta_manova_trends.csv"
  )
)

cat("\nMain trend table:\n")
print(main_summary_wide)

## ============================================================
## Output listing
## ============================================================

cat("\n[DONE] Summary outputs written to:\n  ",
    SUMMARY_DIR,
    "\n\n",
    sep = "")

output_files <- list.files(
  SUMMARY_DIR,
  full.names = FALSE
)

cat(paste0("  ", output_files, collapse = "\n"))
cat("\n")