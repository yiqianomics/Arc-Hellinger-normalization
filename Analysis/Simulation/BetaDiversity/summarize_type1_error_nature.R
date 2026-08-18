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
  library(tidyr)
})

script_dir <- get_script_dir()
raw_file <- file.path(script_dir, "beta_div_sparse12_all_reps.csv")
output_dir <- file.path(script_dir, "weighted_power_output")

if (!file.exists(raw_file)) {
  stop("Raw simulation table not found: ", raw_file)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
  "bray" = "Bray--Curtis",
  "euclidean" = "Euclidean",
  "jaccard" = "Jaccard",
  "aitchison_pc0.01" = "Aitchison (0.01)",
  "aitchison_pc0.1" = "Aitchison (0.1)",
  "aitchison" = "Aitchison (0.5)",
  "aitchison_pc1" = "Aitchison (1)",
  "ILR_euclidean_pc0.01" = "\\(\\operatorname{ILR}\\)--Euclidean (0.01)",
  "ILR_euclidean_pc0.1" = "\\(\\operatorname{ILR}\\)--Euclidean (0.1)",
  "ILR_euclidean_pc0.5" = "\\(\\operatorname{ILR}\\)--Euclidean (0.5)",
  "ILR_euclidean_pc1" = "\\(\\operatorname{ILR}\\)--Euclidean (1)"
)

p_columns <- paste0("p_", method_order)
required_columns <- c(
  "k",
  "signal_mode",
  "beta",
  "signal_cluster",
  "cluster_n_otus_template",
  "cluster_abund_sum_template",
  "cluster_ordered",
  p_columns
)

raw_data <- read.csv(raw_file, stringsAsFactors = FALSE)
missing_columns <- setdiff(required_columns, names(raw_data))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

null_data <- raw_data |>
  filter(k %in% c(2, 4, 6, 8), beta == 0)

# Match panel c: choose the largest taxon cluster within every k and nominal
# signal mode, breaking size ties by fitted abundance rank.
cluster_choice <- null_data |>
  distinct(
    k,
    signal_mode,
    signal_cluster,
    cluster_n_otus_template,
    cluster_abund_sum_template,
    cluster_ordered
  ) |>
  arrange(
    k,
    signal_mode,
    desc(cluster_n_otus_template),
    cluster_ordered,
    signal_cluster
  ) |>
  group_by(k, signal_mode) |>
  slice(1L) |>
  ungroup()

selected_null <- null_data |>
  inner_join(
    select(cluster_choice, k, signal_mode, signal_cluster),
    by = c("k", "signal_mode", "signal_cluster")
  )

type1_by_setting <- selected_null |>
  select(k, signal_mode, all_of(p_columns)) |>
  pivot_longer(
    cols = all_of(p_columns),
    names_to = "method",
    values_to = "p_value"
  ) |>
  mutate(
    method = sub("^p_", "", method),
    false_rejection = p_value < 0.05
  ) |>
  group_by(k, signal_mode, method) |>
  summarise(
    type1_error = mean(false_rejection, na.rm = TRUE),
    n_replicates = sum(is.finite(p_value)),
    .groups = "drop"
  ) |>
  mutate(
    signal_mode = recode(
      signal_mode,
      abundance = "Abundance",
      both = "Joint",
      prevalence = "Prevalence"
    ),
    method = factor(method, levels = method_order),
    method_label = unname(method_labels[as.character(method)])
  ) |>
  arrange(method, signal_mode, k)

if (any(type1_by_setting$n_replicates != 100L)) {
  stop("Every null setting must contain exactly 100 finite P values.")
}

type1_summary <- type1_by_setting |>
  group_by(method, method_label) |>
  summarise(
    Abundance = mean(type1_error[signal_mode == "Abundance"]),
    Joint = mean(type1_error[signal_mode == "Joint"]),
    Prevalence = mean(type1_error[signal_mode == "Prevalence"]),
    Overall = mean(type1_error),
    Minimum = min(type1_error),
    Maximum = max(type1_error),
    .groups = "drop"
  ) |>
  arrange(method) |>
  mutate(
    overall_range = sprintf(
      "%.4f (%.4f--%.4f)",
      Overall,
      Minimum,
      Maximum
    )
  )

write.csv(
  type1_by_setting,
  file.path(output_dir, "type1_error_by_k_and_signal_mode.csv"),
  row.names = FALSE
)
write.csv(
  type1_summary,
  file.path(output_dir, "type1_error_summary_nature.csv"),
  row.names = FALSE
)
write.csv(
  cluster_choice,
  file.path(output_dir, "type1_error_selected_clusters.csv"),
  row.names = FALSE
)

table_rows <- sprintf(
  "%s & %.4f & %.4f & %.4f & %s \\\\",
  type1_summary$method_label,
  type1_summary$Abundance,
  type1_summary$Joint,
  type1_summary$Prevalence,
  type1_summary$overall_range
)

table_lines <- c(
  "% Requires \\usepackage{booktabs}.",
  "\\begin{table*}[t]",
  "    \\centering",
  "    \\caption{\\textbf{Empirical type I error of distance-based PERMANOVA in null simulations.}",
  "    Values under each nominal signal mode are mean rejection proportions across",
  "    \\(k=2,4,6\\) and 8 at effect size zero. The final column reports the",
  "    overall mean and the range across all 12 combinations of mode and \\(k\\).",
  "    Each combination contained 100 independent simulated datasets, and tests",
  "    used a nominal level of 0.05. Parenthetical values in method names",
  "    denote the pseudocount added before log-ratio transformation.}",
  "    \\label{tab:simulation-type1}",
  "    \\begin{tabular}{lcccc}",
  "        \\toprule",
  "        Method & Abundance & Joint & Prevalence & Overall mean (range) \\\\",
  "        \\midrule",
  paste0("        ", table_rows),
  "        \\bottomrule",
  "    \\end{tabular}",
  "\\end{table*}"
)

writeLines(
  table_lines,
  file.path(output_dir, "table_type1_error_nature.tex")
)

print(type1_summary)
message("Saved type I error summaries and LaTeX table to: ", output_dir)
