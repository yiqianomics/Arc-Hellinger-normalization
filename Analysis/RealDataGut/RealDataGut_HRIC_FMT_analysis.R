############################################################
## HRIC analysis of the randomized auto-FMT trial
## Taur et al., Science Translational Medicine (2018)
############################################################

options(stringsAsFactors = FALSE, width = 140)
set.seed(20260723)

required_packages <- c(
  "HRIC", "phyloseq", "ggplot2", "dplyr", "tidyr", "tibble",
  "patchwork", "scales", "splines", "lme4", "sandwich", "reformulas"
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
  library(patchwork)
  library(splines)
  library(lme4)
})

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  normalizePath(
    "RealDataGut/RealDataGut_HRIC_FMT_analysis.R",
    mustWork = TRUE
  )
}

analysis_dir <- dirname(get_script_path())
out_dir <- file.path(analysis_dir, "hric_fmt_diversity_results")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, filename, width, height, dpi = 600) {
  paths <- file.path(fig_dir, paste0(filename, c(".png", ".pdf", ".tiff")))
  ggsave(paths[1], plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(paths[2], plot, width = width, height = height, bg = "white")
  ggsave(
    paths[3], plot, width = width, height = height, dpi = dpi,
    compression = "lzw", bg = "white"
  )
  paths
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

theme_nature <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.45, color = "grey10"),
      axis.ticks = element_line(linewidth = 0.40, color = "grey10"),
      axis.ticks.length = grid::unit(1.6, "mm"),
      axis.text = element_text(color = "grey10", size = base_size - 0.5),
      axis.title = element_text(
        color = "grey5", face = "bold", size = base_size + 0.5
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold", color = "grey5", size = base_size + 0.5,
        hjust = 0
      ),
      legend.title = element_text(face = "bold", color = "grey5"),
      legend.text = element_text(color = "grey10"),
      plot.title = element_text(
        face = "bold", size = base_size + 1.5, color = "grey5", hjust = 0
      ),
      plot.subtitle = element_text(
        size = base_size - 0.2, color = "grey25", hjust = 0,
        margin = margin(b = 5)
      ),
      plot.tag = element_text(face = "bold", size = base_size + 3),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(7, 9, 7, 9)
    )
}

mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  se <- if (n > 1) stats::sd(x) / sqrt(n) else 0
  tibble(
    y = mean(x),
    ymin = mean(x) - stats::qt(0.975, max(1, n - 1)) * se,
    ymax = mean(x) + stats::qt(0.975, max(1, n - 1)) * se
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
  vapply(seq_len(nrow(taxonomy)), function(i) {
    for (rank in available) {
      value <- clean_taxon_value(taxonomy[i, rank])
      if (!is.na(value)) {
        return(paste0(value, " (", tolower(rank), ")"))
      }
    }
    "Unclassified bacteria"
  }, character(1))
}

make_alluvial_ribbons <- function(data, taxon_levels, points_per_link = 18L) {
  ordered <- data %>%
    mutate(Taxon = factor(Taxon, levels = taxon_levels)) %>%
    arrange(PatientLabel, SampleOrder, Taxon) %>%
    group_by(PatientLabel, SampleOrder) %>%
    mutate(
      Lower = lag(cumsum(Contribution), default = 0),
      Upper = cumsum(Contribution)
    ) %>%
    ungroup()

  ribbons <- bind_rows(lapply(split(ordered, ordered$PatientLabel), function(d) {
    sample_orders <- sort(unique(d$SampleOrder))
    if (length(sample_orders) < 2) return(tibble())
    bind_rows(lapply(seq_len(length(sample_orders) - 1L), function(i) {
      left <- d %>% filter(SampleOrder == sample_orders[i]) %>% arrange(Taxon)
      right <- d %>% filter(SampleOrder == sample_orders[i + 1L]) %>% arrange(Taxon)
      bind_rows(lapply(seq_along(taxon_levels), function(j) {
        s <- seq(0, 1, length.out = points_per_link)
        smooth <- 3 * s^2 - 2 * s^3
        x <- left$SampleOrder[j] +
          (right$SampleOrder[j] - left$SampleOrder[j]) * s
        lower <- left$Lower[j] + (right$Lower[j] - left$Lower[j]) * smooth
        upper <- left$Upper[j] + (right$Upper[j] - left$Upper[j]) * smooth
        tibble(
          PatientLabel = left$PatientLabel[j],
          Taxon = taxon_levels[j],
          RibbonID = paste(left$PatientLabel[j], i, j, sep = "_"),
          x = c(x, rev(x)),
          y = c(lower, rev(upper))
        )
      }))
    }))
  }))

  list(ordered = ordered, ribbons = ribbons)
}

predict_fixed_curve <- function(model, newdata, covariance = NULL) {
  is_mixed <- inherits(model, "merMod")
  fixed_formula <- if (is_mixed) {
    reformulas::nobars(formula(model))
  } else {
    formula(model)
  }
  fixed_terms <- delete.response(terms(fixed_formula))
  design <- model.matrix(fixed_terms, newdata)
  coefficients <- if (is_mixed) lme4::fixef(model) else coef(model)
  coefficient_names <- names(coefficients)
  design <- design[, coefficient_names, drop = FALSE]
  if (is.null(covariance)) covariance <- vcov(model)
  covariance <- as.matrix(covariance)
  fit <- as.numeric(design %*% coefficients)
  se <- sqrt(pmax(0, rowSums((design %*% covariance) * design)))
  bind_cols(
    newdata,
    tibble(
      Estimate = fit,
      Lower = fit - 1.96 * se,
      Upper = fit + 1.96 * se
    )
  )
}

fixed_contrast <- function(
  model,
  newdata_1,
  newdata_0,
  covariance = NULL
) {
  is_mixed <- inherits(model, "merMod")
  fixed_formula <- if (is_mixed) {
    reformulas::nobars(formula(model))
  } else {
    formula(model)
  }
  fixed_terms <- delete.response(terms(fixed_formula))
  coefficients <- if (is_mixed) lme4::fixef(model) else coef(model)
  coefficient_names <- names(coefficients)
  x1 <- model.matrix(fixed_terms, newdata_1)[, coefficient_names, drop = FALSE]
  x0 <- model.matrix(fixed_terms, newdata_0)[, coefficient_names, drop = FALSE]
  contrast <- x1 - x0
  if (is.null(covariance)) covariance <- vcov(model)
  covariance <- as.matrix(covariance)
  estimate <- as.numeric(contrast %*% coefficients)
  se <- sqrt(pmax(0, rowSums((contrast %*% covariance) * contrast)))
  tibble(
    Estimate = estimate,
    SE = se,
    Lower = estimate - 1.96 * se,
    Upper = estimate + 1.96 * se,
    P = 2 * pnorm(abs(estimate / se), lower.tail = FALSE)
  )
}

fixed_linear_contrast <- function(
  model,
  contrast,
  covariance = NULL
) {
  coefficients <- if (inherits(model, "merMod")) {
    lme4::fixef(model)
  } else {
    coef(model)
  }
  contrast <- as.numeric(contrast[names(coefficients)])
  if (is.null(covariance)) covariance <- vcov(model)
  covariance <- as.matrix(covariance)[
    names(coefficients), names(coefficients), drop = FALSE
  ]
  estimate <- sum(contrast * coefficients)
  se <- sqrt(max(0, as.numeric(
    t(contrast) %*% covariance %*% contrast
  )))
  tibble(
    Estimate = estimate,
    SE = se,
    Lower = estimate - 1.96 * se,
    Upper = estimate + 1.96 * se,
    P = 2 * pnorm(abs(estimate / se), lower.tail = FALSE)
  )
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

ps_gut <- prune_samples(sample_sums(gut) > 0, gut)
ps_gut <- prune_taxa(taxa_sums(ps_gut) > 0, ps_gut)

meta <- as(sample_data(ps_gut), "data.frame")
meta$SampleID <- rownames(meta)
required_meta <- c(
  "sample_name", "arm", "day", "engraftmentday",
  "host_subject_id", "randomizationday", "inversesimpson"
)
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "))
}

meta <- meta %>%
  mutate(
    arm = factor(tolower(as.character(arm)), levels = c("control", "treatment")),
    Arm = factor(
      ifelse(arm == "control", "Control", "Auto-FMT"),
      levels = c("Control", "Auto-FMT")
    ),
    day = as.numeric(day),
    engraftmentday = as.numeric(engraftmentday),
    randomizationday = as.numeric(randomizationday),
    InverseSimpson = as.numeric(inversesimpson),
    Patient = factor(as.character(host_subject_id)),
    EventTime = day - randomizationday,
    PostIndex = as.integer(day >= randomizationday),
    TreatmentPost = as.integer(arm == "treatment" & PostIndex == 1)
  )

if (
  any(!is.finite(meta$InverseSimpson)) ||
    any(meta$InverseSimpson < 1)
) {
  stop("Metadata inversesimpson values must be finite and at least 1.")
}

counts <- as(otu_table(ps_gut), "matrix")
if (taxa_are_rows(ps_gut)) counts <- t(counts)
storage.mode(counts) <- "numeric"
counts <- counts[meta$SampleID, , drop = FALSE]
counts <- counts[, colSums(counts) > 0, drop = FALSE]

# Taur et al. quantified alpha diversity as
# log(1 / sum_j p_ij^2), where p_ij is the relative abundance of OTU j.
relative_counts <- counts / rowSums(counts)
meta$InverseSimpsonCalculated <- 1 / rowSums(relative_counts^2)
meta$LogInverseSimpsonCalculated <- log(
  meta$InverseSimpsonCalculated
)
observed_richness <- rowSums(relative_counts > 0)
shannon_entropy <- -rowSums(
  ifelse(
    relative_counts > 0,
    relative_counts * log(relative_counts),
    0
  )
)
meta$PielouEvenness <- shannon_entropy / log(observed_richness)
if (
  any(!is.finite(meta$InverseSimpsonCalculated)) ||
    any(meta$InverseSimpsonCalculated < 1) ||
    any(!is.finite(meta$LogInverseSimpsonCalculated))
) {
  stop("Calculated inverse Simpson values are invalid.")
}
if (
  any(observed_richness < 2) ||
    any(!is.finite(meta$PielouEvenness)) ||
    any(meta$PielouEvenness < 0) ||
    any(meta$PielouEvenness > 1)
) {
  stop("Calculated Pielou evenness values are invalid.")
}

# Every HRIC calculation uses this same feature dictionary. This matters because
# the package's alpha scale depends on the number of simplex components.
stopifnot(identical(rownames(counts), meta$SampleID))
alpha_all <- HRIC::SHalpha(counts)
meta$Alpha <- as.numeric(alpha_all[meta$SampleID])

# Squared HRIC coordinates provide an additive decomposition of each sample's
# alpha deficit: sum_j Z_ij^2 / A_p^2 = 1 - alpha_i.
Z <- HRIC::HRIC(counts)
hric_radius <- asin(sqrt(1 - 1 / ncol(counts)))
feature_alpha_deficit <- Z^2 / hric_radius^2
stopifnot(
  max(abs(rowSums(feature_alpha_deficit) - (1 - meta$Alpha))) < 1e-10
)

taxonomy <- as(tax_table(ps_gut), "matrix")
taxonomy <- taxonomy[colnames(counts), , drop = FALSE]
bacteroidetes_features <- rownames(taxonomy)[which(
  clean_taxon_value(taxonomy[, "Phylum"]) == "Bacteroidetes"
)]
meta$BacteroidetesPercent <- 100 * rowSums(
  counts[, bacteroidetes_features, drop = FALSE]
) / rowSums(counts)

taxon_labels <- deepest_taxon_label(taxonomy)
taxon_alpha_deficit <- t(rowsum(
  t(feature_alpha_deficit),
  group = taxon_labels,
  reorder = FALSE
))
rownames(taxon_alpha_deficit) <- meta$SampleID

arm_colors <- c(Control = "#2F6FA3", `Auto-FMT` = "#C94F3D")
component_colors <- c(`Mean alpha` = "#4D6F8C", Beta = "#D9A441")
phase_levels <- c(
  "Pre-HSCT",
  "Post-HSCT / pre-index",
  "Post-index"
)
meta <- meta %>%
  mutate(
    Phase = factor(
      case_when(
        day < 0 ~ "Pre-HSCT",
        day <= randomizationday ~ "Post-HSCT / pre-index",
        TRUE ~ "Post-index"
      ),
      levels = phase_levels
    )
  )

############################################################
## Figure 1a-d: HRIC alpha trajectories for all patients
############################################################

extract_participant_id <- function(sample_name) {
  sample_name <- as.character(sample_name)
  is_fmt_name <- grepl("\\.FMT\\.", sample_name)
  participant_id <- rep(NA_character_, length(sample_name))
  participant_id[is_fmt_name] <- sub(
    "^.*\\.FMT\\.([0-9]{4}).*$", "\\1", sample_name[is_fmt_name]
  )
  participant_id[!is_fmt_name] <- sub(
    "^.*\\.([0-9]{3,4})[A-Za-z].*$", "\\1", sample_name[!is_fmt_name]
  )
  sprintf("%04d", as.integer(participant_id))
}

participant_mapping <- meta %>%
  mutate(
    Patient = as.character(Patient),
    CandidateID = extract_participant_id(sample_name),
    UsesFMTIdentifier = grepl(".FMT.", sample_name, fixed = TRUE)
  ) %>%
  group_by(Patient, Arm) %>%
  summarise(
    ParticipantID = if (any(UsesFMTIdentifier)) {
      first(CandidateID[UsesFMTIdentifier])
    } else {
      first(CandidateID)
    },
    .groups = "drop"
  )

if (
  anyNA(participant_mapping$ParticipantID) ||
    anyDuplicated(participant_mapping$ParticipantID)
) {
  stop("Each internal patient ID must map to one participant identifier.")
}

meta$ParticipantID <- participant_mapping$ParticipantID[
  match(as.character(meta$Patient), participant_mapping$Patient)
]

all_patient_key <- meta %>%
  mutate(Patient = as.character(Patient)) %>%
  group_by(Patient, Arm, ParticipantID) %>%
  summarise(
    Samples = n(),
    EngraftmentDay = first(engraftmentday),
    IndexDay = first(randomizationday),
    .groups = "drop"
  ) %>%
  arrange(Arm, as.integer(ParticipantID)) %>%
  group_by(Arm) %>%
  mutate(
    ArmOrder = row_number(),
    AtlasPage = case_when(
      Arm == "Control" ~ ceiling(ArmOrder / 6),
      TRUE ~ ceiling(ArmOrder / 7)
    ),
    PatientLabel = sprintf("Patient %s (n=%d)", ParticipantID, Samples)
  ) %>%
  ungroup()

all_patient_baseline <- meta %>%
  mutate(Patient = as.character(Patient)) %>%
  filter(day < 0) %>%
  group_by(Patient) %>%
  arrange(day, SampleID, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    Patient,
    BaselineAlpha = Alpha,
    BaselineDay = day,
    Recovery80 = 0.8 * Alpha
  )

all_patient_recovery <- bind_rows(lapply(
  seq_len(nrow(all_patient_key)),
  function(i) {
    patient_row <- all_patient_key[i, ]
    recovery_threshold <- all_patient_baseline$Recovery80[
      all_patient_baseline$Patient == patient_row$Patient
    ]
    post_days <- meta %>%
      filter(
        as.character(Patient) == patient_row$Patient,
        EventTime >= 1,
        EventTime <= 60
      ) %>%
      group_by(EventTime) %>%
      summarise(Alpha = mean(Alpha), .groups = "drop") %>%
      arrange(EventTime)

    sustained_day <- NA_real_
    if (nrow(post_days) >= 2) {
      crossings <- which(post_days$Alpha >= recovery_threshold)
      for (crossing in crossings) {
        later <- which(post_days$EventTime > post_days$EventTime[crossing])
        if (
          length(later) > 0 &&
            post_days$Alpha[later[1]] >= recovery_threshold
        ) {
          sustained_day <- post_days$EventTime[crossing]
          break
        }
      }
    }

    tibble(
      Patient = patient_row$Patient,
      SustainedRecoveryDay = sustained_day,
      RecoveryLabel = if (is.na(sustained_day)) {
        "80%: not sustained by +60 d"
      } else {
        paste0("80%: +", sustained_day, " d")
      }
    )
  }
))

all_patient_trajectory <- meta %>%
  mutate(Patient = as.character(Patient)) %>%
  left_join(all_patient_key, by = c("Patient", "Arm", "ParticipantID")) %>%
  left_join(all_patient_baseline, by = "Patient") %>%
  left_join(all_patient_recovery, by = "Patient") %>%
  arrange(Arm, ArmOrder, day, SampleID)

all_patient_events <- all_patient_key %>%
  select(Patient, Arm, PatientLabel, AtlasPage, EngraftmentDay, IndexDay) %>%
  pivot_longer(
    cols = c(EngraftmentDay, IndexDay),
    names_to = "EventType",
    values_to = "EventDay"
  ) %>%
  mutate(
    Event = case_when(
      EventType == "EngraftmentDay" ~ "Engraftment",
      Arm == "Control" ~ "Randomization (no FMT)",
      TRUE ~ "Auto-FMT / randomization"
    ),
    Event = factor(
      Event,
      levels = c(
        "Engraftment", "Randomization (no FMT)",
        "Auto-FMT / randomization"
      )
    )
  )

atlas_x_limits <- range(meta$day, na.rm = TRUE) + c(-2, 2)
atlas_y_upper <- max(meta$Alpha, na.rm = TRUE) * 1.06
atlas_event_colors <- c(
  Engraftment = "#686868",
  `Randomization (no FMT)` = "#7057A3",
  `Auto-FMT / randomization` = "#D97706"
)

make_patient_alpha_atlas <- function(arm_name, page_number, tag) {
  plot_data <- all_patient_trajectory %>%
    filter(Arm == arm_name, AtlasPage == page_number)
  patient_levels <- all_patient_key %>%
    filter(Arm == arm_name, AtlasPage == page_number) %>%
    arrange(ArmOrder) %>%
    pull(PatientLabel)
  plot_data <- plot_data %>%
    mutate(PatientLabel = factor(PatientLabel, levels = patient_levels))
  baseline_data <- plot_data %>%
    distinct(PatientLabel, BaselineAlpha, Recovery80)
  recovery_data <- plot_data %>%
    distinct(PatientLabel, RecoveryLabel)
  event_data <- all_patient_events %>%
    filter(Arm == arm_name, AtlasPage == page_number) %>%
    mutate(PatientLabel = factor(PatientLabel, levels = patient_levels))

  ggplot(
    plot_data,
    aes(x = day, y = Alpha, group = PatientLabel)
  ) +
    annotate(
      "rect",
      xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
      fill = "grey70", alpha = 0.10
    ) +
    geom_vline(
      xintercept = 0, linewidth = 0.52, linetype = "22", color = "grey15"
    ) +
    geom_vline(
      data = event_data,
      aes(xintercept = EventDay, color = Event),
      inherit.aes = FALSE, linewidth = 0.58, linetype = "13"
    ) +
    geom_hline(
      data = baseline_data,
      aes(yintercept = BaselineAlpha),
      inherit.aes = FALSE, linewidth = 0.50,
      linetype = "longdash", color = "grey30"
    ) +
    geom_hline(
      data = baseline_data,
      aes(yintercept = Recovery80),
      inherit.aes = FALSE, linewidth = 0.45,
      linetype = "dotted", color = "grey55"
    ) +
    geom_line(
      linewidth = 0.62, alpha = 0.82,
      color = unname(arm_colors[[arm_name]])
    ) +
    geom_point(
      shape = 21, fill = "white", stroke = 0.68, size = 2.0,
      color = unname(arm_colors[[arm_name]])
    ) +
    geom_point(
      data = plot_data %>% filter(day == IndexDay),
      shape = 23, fill = "white", stroke = 0.82, size = 2.7,
      color = unname(arm_colors[[arm_name]])
    ) +
    geom_label(
      data = recovery_data,
      aes(x = Inf, y = Inf, label = RecoveryLabel),
      inherit.aes = FALSE,
      hjust = 1.03, vjust = 1.20,
      size = 2.15, linewidth = 0.16,
      label.padding = grid::unit(0.65, "mm"),
      fill = scales::alpha("white", 0.90), color = "grey20"
    ) +
    facet_wrap(vars(PatientLabel), ncol = 3) +
    scale_color_manual(
      values = atlas_event_colors,
      breaks = c(
        "Engraftment", "Randomization (no FMT)",
        "Auto-FMT / randomization"
      )
    ) +
    scale_x_continuous(
      limits = atlas_x_limits,
      breaks = c(-20, 0, 30, 60, 90),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = c(0, atlas_y_upper),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      tag = tag,
      title = paste0(
        arm_name, " patients, page ", page_number, " of 2"
      ),
      subtitle = paste0(
        "Shared axes; black dashed line marks allo-HSCT day 0. ",
        "n in each strip is the longitudinal stool-sample count."
      ),
      x = "Day relative to allo-HSCT infusion",
      y = "HRIC alpha diversity",
      color = NULL
    ) +
    theme_nature(9.5) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.width = grid::unit(5.2, "mm"),
      panel.spacing = grid::unit(4.2, "mm"),
      strip.text = element_text(size = 9.2, face = "bold", hjust = 0),
      plot.title = element_text(size = 12.0),
      plot.subtitle = element_text(size = 8.8),
      plot.tag = element_text(size = 13.5),
      axis.text = element_text(size = 8.0),
      axis.title = element_text(size = 9.4)
    )
}

patient_alpha_atlas_specs <- tribble(
  ~Arm, ~Page, ~Tag, ~Filename, ~Height,
  "Control", 1, "a", "figure1a_all_patient_alpha_control_1", 6.0,
  "Control", 2, "b", "figure1b_all_patient_alpha_control_2", 6.0,
  "Auto-FMT", 1, "c", "figure1c_all_patient_alpha_autofmt_1", 7.9,
  "Auto-FMT", 2, "d", "figure1d_all_patient_alpha_autofmt_2", 7.9
)

patient_alpha_atlas_files <- unlist(lapply(
  seq_len(nrow(patient_alpha_atlas_specs)),
  function(i) {
    spec <- patient_alpha_atlas_specs[i, ]
    atlas_plot <- make_patient_alpha_atlas(
      spec$Arm, spec$Page, spec$Tag
    )
    save_plot(
      atlas_plot, spec$Filename,
      width = 7.2, height = spec$Height
    )
  }
))

write.csv(
  all_patient_trajectory %>%
    select(
      SampleID, sample_name, Patient, ParticipantID, PatientLabel,
      Arm, AtlasPage, day, engraftmentday, randomizationday,
      EventTime, Alpha, BaselineAlpha, Recovery80
    ),
  file.path(tab_dir, "figure1_all_patient_alpha_trajectories.csv"),
  row.names = FALSE
)
write.csv(
  all_patient_key %>%
    left_join(all_patient_baseline, by = "Patient") %>%
    left_join(all_patient_recovery, by = "Patient"),
  file.path(tab_dir, "figure1_all_patient_recovery_summary.csv"),
  row.names = FALSE
)

############################################################
## Figure 1 alternative: original metadata inverse Simpson
############################################################

all_patient_is_baseline <- meta %>%
  mutate(Patient = as.character(Patient)) %>%
  filter(day < 0) %>%
  group_by(Patient) %>%
  arrange(day, SampleID, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    Patient,
    BaselineInverseSimpson = InverseSimpson,
    BaselineDay = day,
    Recovery80InverseSimpson = 0.8 * InverseSimpson
  )

all_patient_is_recovery <- bind_rows(lapply(
  seq_len(nrow(all_patient_key)),
  function(i) {
    patient_row <- all_patient_key[i, ]
    recovery_threshold <-
      all_patient_is_baseline$Recovery80InverseSimpson[
        all_patient_is_baseline$Patient == patient_row$Patient
      ]
    post_days <- meta %>%
      filter(
        as.character(Patient) == patient_row$Patient,
        EventTime >= 1,
        EventTime <= 60
      ) %>%
      group_by(EventTime) %>%
      summarise(
        InverseSimpson = mean(InverseSimpson),
        .groups = "drop"
      ) %>%
      arrange(EventTime)

    sustained_day <- NA_real_
    if (nrow(post_days) >= 2) {
      crossings <- which(
        post_days$InverseSimpson >= recovery_threshold
      )
      for (crossing in crossings) {
        later <- which(post_days$EventTime > post_days$EventTime[crossing])
        if (
          length(later) > 0 &&
            post_days$InverseSimpson[later[1]] >= recovery_threshold
        ) {
          sustained_day <- post_days$EventTime[crossing]
          break
        }
      }
    }

    tibble(
      Patient = patient_row$Patient,
      SustainedRecoveryDay = sustained_day,
      RecoveryLabel = if (is.na(sustained_day)) {
        "80% threshold not sustained"
      } else {
        paste0("80% threshold: +", sustained_day, " d")
      }
    )
  }
))

all_patient_is_trajectory <- meta %>%
  mutate(Patient = as.character(Patient)) %>%
  left_join(all_patient_key, by = c("Patient", "Arm", "ParticipantID")) %>%
  left_join(all_patient_is_baseline, by = "Patient") %>%
  left_join(all_patient_is_recovery, by = "Patient") %>%
  arrange(Arm, ArmOrder, day, SampleID)

atlas_is_y_upper <- max(meta$InverseSimpson, na.rm = TRUE) * 1.06

make_patient_inverse_simpson_atlas <- function(
  arm_name,
  page_number,
  tag
) {
  plot_data <- all_patient_is_trajectory %>%
    filter(Arm == arm_name, AtlasPage == page_number)
  patient_levels <- all_patient_key %>%
    filter(Arm == arm_name, AtlasPage == page_number) %>%
    arrange(ArmOrder) %>%
    pull(PatientLabel)
  plot_data <- plot_data %>%
    mutate(PatientLabel = factor(PatientLabel, levels = patient_levels))
  line_data <- plot_data %>%
    mutate(
      TimeSegment = case_when(
        day < IndexDay ~ "Pre-index",
        day > IndexDay ~ "Post-index",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(TimeSegment))
  baseline_data <- plot_data %>%
    distinct(
      PatientLabel, BaselineInverseSimpson,
      Recovery80InverseSimpson
    )
  recovery_data <- plot_data %>%
    distinct(PatientLabel, RecoveryLabel)
  event_data <- all_patient_events %>%
    filter(Arm == arm_name, AtlasPage == page_number) %>%
    mutate(PatientLabel = factor(PatientLabel, levels = patient_levels))

  ggplot(
    plot_data,
    aes(x = day, y = InverseSimpson, group = PatientLabel)
  ) +
    annotate(
      "rect",
      xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
      fill = "grey70", alpha = 0.10
    ) +
    geom_vline(
      xintercept = 0, linewidth = 0.52,
      linetype = "22", color = "grey15"
    ) +
    geom_vline(
      data = event_data,
      aes(xintercept = EventDay, color = Event),
      inherit.aes = FALSE, linewidth = 0.58, linetype = "13"
    ) +
    geom_hline(
      data = baseline_data,
      aes(yintercept = BaselineInverseSimpson),
      inherit.aes = FALSE, linewidth = 0.50,
      linetype = "longdash", color = "grey30"
    ) +
    geom_hline(
      data = baseline_data,
      aes(yintercept = Recovery80InverseSimpson),
      inherit.aes = FALSE, linewidth = 0.45,
      linetype = "dotted", color = "grey55"
    ) +
    geom_line(
      data = line_data,
      aes(
        group = interaction(PatientLabel, TimeSegment, drop = TRUE)
      ),
      linewidth = 0.62, alpha = 0.82,
      color = unname(arm_colors[[arm_name]])
    ) +
    geom_point(
      data = plot_data %>% filter(day != IndexDay),
      shape = 21, fill = "white", stroke = 0.68, size = 2.0,
      color = unname(arm_colors[[arm_name]])
    ) +
    geom_point(
      data = plot_data %>% filter(day == IndexDay),
      shape = 23, fill = "grey88", stroke = 0.82, size = 2.7,
      color = "grey25"
    ) +
    geom_label(
      data = recovery_data,
      aes(x = Inf, y = Inf, label = RecoveryLabel),
      inherit.aes = FALSE,
      hjust = 1.03, vjust = 1.20,
      size = 2.10, linewidth = 0.16,
      label.padding = grid::unit(0.65, "mm"),
      fill = scales::alpha("white", 0.90), color = "grey20"
    ) +
    facet_wrap(vars(PatientLabel), ncol = 3) +
    scale_color_manual(
      values = atlas_event_colors,
      breaks = c(
        "Engraftment", "Randomization (no FMT)",
        "Auto-FMT / randomization"
      )
    ) +
    scale_x_continuous(
      limits = atlas_x_limits,
      breaks = c(-20, 0, 30, 60, 90),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = c(0, atlas_is_y_upper),
      breaks = seq(0, 25, by = 5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      tag = tag,
      title = paste0(
        arm_name, " patients, page ", page_number, " of 2"
      ),
      subtitle = paste0(
        "Original metadata inverse Simpson; grey diamonds are index-day ",
        "stools with unknown within-day order."
      ),
      x = "Day relative to allo-HSCT infusion",
      y = "Inverse Simpson diversity",
      color = NULL
    ) +
    theme_nature(9.5) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.width = grid::unit(5.2, "mm"),
      panel.spacing = grid::unit(4.2, "mm"),
      strip.text = element_text(size = 9.2, face = "bold", hjust = 0),
      plot.title = element_text(size = 12.0),
      plot.subtitle = element_text(size = 8.6),
      plot.tag = element_text(size = 13.5),
      axis.text = element_text(size = 8.0),
      axis.title = element_text(size = 9.4)
    )
}

patient_inverse_simpson_atlas_specs <- tribble(
  ~Arm, ~Page, ~Tag, ~Filename, ~Height,
  "Control", 1, "a",
  "figure1a_all_patient_inverse_simpson_control_1", 6.0,
  "Control", 2, "b",
  "figure1b_all_patient_inverse_simpson_control_2", 6.0,
  "Auto-FMT", 1, "c",
  "figure1c_all_patient_inverse_simpson_autofmt_1", 7.9,
  "Auto-FMT", 2, "d",
  "figure1d_all_patient_inverse_simpson_autofmt_2", 7.9
)

patient_inverse_simpson_atlas_files <- unlist(lapply(
  seq_len(nrow(patient_inverse_simpson_atlas_specs)),
  function(i) {
    spec <- patient_inverse_simpson_atlas_specs[i, ]
    atlas_plot <- make_patient_inverse_simpson_atlas(
      spec$Arm, spec$Page, spec$Tag
    )
    save_plot(
      atlas_plot, spec$Filename,
      width = 7.2, height = spec$Height
    )
  }
))

write.csv(
  all_patient_is_trajectory %>%
    select(
      SampleID, sample_name, Patient, ParticipantID, PatientLabel,
      Arm, AtlasPage, day, engraftmentday, randomizationday,
      EventTime, InverseSimpson, BaselineInverseSimpson,
      Recovery80InverseSimpson
    ),
  file.path(
    tab_dir,
    "figure1_all_patient_inverse_simpson_trajectories.csv"
  ),
  row.names = FALSE
)
write.csv(
  all_patient_key %>%
    left_join(all_patient_is_baseline, by = "Patient") %>%
    left_join(all_patient_is_recovery, by = "Patient"),
  file.path(
    tab_dir,
    "figure1_all_patient_inverse_simpson_recovery_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## Figure 1 alternative: inverse Simpson calculated by formula
############################################################

all_patient_formula_is_baseline <- meta %>%
  mutate(Patient = as.character(Patient)) %>%
  filter(day < 0) %>%
  group_by(Patient) %>%
  arrange(day, SampleID, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    Patient,
    BaselineInverseSimpsonCalculated = InverseSimpsonCalculated,
    BaselineLogInverseSimpsonCalculated =
      LogInverseSimpsonCalculated,
    BaselineDay = day,
    Recovery80InverseSimpsonCalculated =
      0.8 * InverseSimpsonCalculated,
    Recovery80LogInverseSimpsonCalculated =
      log(0.8 * InverseSimpsonCalculated)
  )

all_patient_formula_is_recovery <- bind_rows(lapply(
  seq_len(nrow(all_patient_key)),
  function(i) {
    patient_row <- all_patient_key[i, ]
    recovery_threshold <-
      all_patient_formula_is_baseline$
        Recovery80InverseSimpsonCalculated[
          all_patient_formula_is_baseline$Patient ==
            patient_row$Patient
        ]
    post_days <- meta %>%
      filter(
        as.character(Patient) == patient_row$Patient,
        EventTime >= 1,
        EventTime <= 60
      ) %>%
      group_by(EventTime) %>%
      summarise(
        InverseSimpsonCalculated =
          mean(InverseSimpsonCalculated),
        .groups = "drop"
      ) %>%
      arrange(EventTime)

    sustained_day <- NA_real_
    if (nrow(post_days) >= 2) {
      crossings <- which(
        post_days$InverseSimpsonCalculated >= recovery_threshold
      )
      for (crossing in crossings) {
        later <- which(
          post_days$EventTime > post_days$EventTime[crossing]
        )
        if (
          length(later) > 0 &&
            post_days$InverseSimpsonCalculated[later[1]] >=
              recovery_threshold
        ) {
          sustained_day <- post_days$EventTime[crossing]
          break
        }
      }
    }

    tibble(
      Patient = patient_row$Patient,
      SustainedRecoveryDay = sustained_day,
      RecoveryLabel = if (is.na(sustained_day)) {
        "80% threshold not sustained"
      } else {
        paste0("80% threshold: +", sustained_day, " d")
      }
    )
  }
))

all_patient_formula_is_trajectory <- meta %>%
  mutate(Patient = as.character(Patient)) %>%
  left_join(all_patient_key, by = c("Patient", "Arm", "ParticipantID")) %>%
  left_join(all_patient_formula_is_baseline, by = "Patient") %>%
  left_join(all_patient_formula_is_recovery, by = "Patient") %>%
  arrange(Arm, ArmOrder, day, SampleID)

atlas_formula_log_is_y_upper <-
  max(meta$LogInverseSimpsonCalculated, na.rm = TRUE) * 1.04

make_patient_formula_log_is_atlas <- function(
  arm_name,
  page_number,
  tag
) {
  plot_data <- all_patient_formula_is_trajectory %>%
    filter(Arm == arm_name, AtlasPage == page_number)
  patient_levels <- all_patient_key %>%
    filter(Arm == arm_name, AtlasPage == page_number) %>%
    arrange(ArmOrder) %>%
    pull(PatientLabel)
  plot_data <- plot_data %>%
    mutate(PatientLabel = factor(PatientLabel, levels = patient_levels))
  line_data <- plot_data %>%
    mutate(
      TimeSegment = case_when(
        day < IndexDay ~ "Pre-index",
        day > IndexDay ~ "Post-index",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(TimeSegment))
  baseline_data <- plot_data %>%
    distinct(
      PatientLabel, BaselineLogInverseSimpsonCalculated,
      Recovery80LogInverseSimpsonCalculated
    )
  recovery_data <- plot_data %>%
    distinct(PatientLabel, RecoveryLabel)
  event_data <- all_patient_events %>%
    filter(Arm == arm_name, AtlasPage == page_number) %>%
    mutate(PatientLabel = factor(PatientLabel, levels = patient_levels))

  ggplot(
    plot_data,
    aes(
      x = day,
      y = LogInverseSimpsonCalculated,
      group = PatientLabel
    )
  ) +
    annotate(
      "rect",
      xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
      fill = "grey70", alpha = 0.10
    ) +
    geom_vline(
      xintercept = 0, linewidth = 0.52,
      linetype = "22", color = "grey15"
    ) +
    geom_vline(
      data = event_data,
      aes(xintercept = EventDay, color = Event),
      inherit.aes = FALSE, linewidth = 0.58, linetype = "13"
    ) +
    geom_hline(
      data = baseline_data,
      aes(yintercept = BaselineLogInverseSimpsonCalculated),
      inherit.aes = FALSE, linewidth = 0.50,
      linetype = "longdash", color = "grey30"
    ) +
    geom_hline(
      data = baseline_data,
      aes(yintercept = Recovery80LogInverseSimpsonCalculated),
      inherit.aes = FALSE, linewidth = 0.45,
      linetype = "dotted", color = "grey55"
    ) +
    geom_line(
      data = line_data,
      aes(
        group = interaction(PatientLabel, TimeSegment, drop = TRUE)
      ),
      linewidth = 0.62, alpha = 0.82,
      color = unname(arm_colors[[arm_name]])
    ) +
    geom_point(
      data = plot_data %>% filter(day != IndexDay),
      shape = 21, fill = "white", stroke = 0.68, size = 2.0,
      color = unname(arm_colors[[arm_name]])
    ) +
    geom_point(
      data = plot_data %>% filter(day == IndexDay),
      shape = 23, fill = "grey88", stroke = 0.82, size = 2.7,
      color = "grey25"
    ) +
    geom_label(
      data = recovery_data,
      aes(x = Inf, y = Inf, label = RecoveryLabel),
      inherit.aes = FALSE,
      hjust = 1.03, vjust = 1.20,
      size = 2.10, linewidth = 0.16,
      label.padding = grid::unit(0.65, "mm"),
      fill = scales::alpha("white", 0.90), color = "grey20"
    ) +
    facet_wrap(vars(PatientLabel), ncol = 3) +
    scale_color_manual(
      values = atlas_event_colors,
      breaks = c(
        "Engraftment", "Randomization (no FMT)",
        "Auto-FMT / randomization"
      )
    ) +
    scale_x_continuous(
      limits = atlas_x_limits,
      breaks = c(-20, 0, 30, 60, 90),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = c(0, atlas_formula_log_is_y_upper),
      breaks = seq(0, 3.5, by = 0.5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      tag = tag,
      title = paste0(
        arm_name, " patients, page ", page_number, " of 2"
      ),
      subtitle = paste0(
        "Calculated from OTU counts as log[1 / sum(p_j^2)]; ",
        "grey diamonds mark index-day stools."
      ),
      x = "Day relative to allo-HSCT infusion",
      y = "Log inverse Simpson diversity",
      color = NULL
    ) +
    theme_nature(9.5) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.width = grid::unit(5.2, "mm"),
      panel.spacing = grid::unit(4.2, "mm"),
      strip.text = element_text(size = 9.2, face = "bold", hjust = 0),
      plot.title = element_text(size = 12.0),
      plot.subtitle = element_text(size = 8.6),
      plot.tag = element_text(size = 13.5),
      axis.text = element_text(size = 8.0),
      axis.title = element_text(size = 9.4)
    )
}

patient_formula_log_is_atlas_specs <- tribble(
  ~Arm, ~Page, ~Tag, ~Filename, ~Height,
  "Control", 1, "a",
  "figure1a_all_patient_log_inverse_simpson_calculated_control_1", 6.0,
  "Control", 2, "b",
  "figure1b_all_patient_log_inverse_simpson_calculated_control_2", 6.0,
  "Auto-FMT", 1, "c",
  "figure1c_all_patient_log_inverse_simpson_calculated_autofmt_1", 7.9,
  "Auto-FMT", 2, "d",
  "figure1d_all_patient_log_inverse_simpson_calculated_autofmt_2", 7.9
)

patient_formula_log_is_atlas_files <- unlist(lapply(
  seq_len(nrow(patient_formula_log_is_atlas_specs)),
  function(i) {
    spec <- patient_formula_log_is_atlas_specs[i, ]
    atlas_plot <- make_patient_formula_log_is_atlas(
      spec$Arm, spec$Page, spec$Tag
    )
    save_plot(
      atlas_plot, spec$Filename,
      width = 7.2, height = spec$Height
    )
  }
))

write.csv(
  all_patient_formula_is_trajectory %>%
    select(
      SampleID, sample_name, Patient, ParticipantID, PatientLabel,
      Arm, AtlasPage, day, engraftmentday, randomizationday,
      EventTime, InverseSimpson, InverseSimpsonCalculated,
      LogInverseSimpsonCalculated,
      BaselineInverseSimpsonCalculated,
      BaselineLogInverseSimpsonCalculated,
      Recovery80InverseSimpsonCalculated,
      Recovery80LogInverseSimpsonCalculated
    ),
  file.path(
    tab_dir,
    "figure1_all_patient_log_inverse_simpson_calculated_trajectories.csv"
  ),
  row.names = FALSE
)
write.csv(
  all_patient_key %>%
    left_join(all_patient_formula_is_baseline, by = "Patient") %>%
    left_join(all_patient_formula_is_recovery, by = "Patient"),
  file.path(
    tab_dir,
    "figure1_all_patient_log_inverse_simpson_calculated_recovery_summary.csv"
  ),
  row.names = FALSE
)
write.csv(
  meta %>%
    transmute(
      SampleID,
      MetadataInverseSimpson = InverseSimpson,
      CalculatedInverseSimpson = InverseSimpsonCalculated,
      CalculatedLogInverseSimpson = LogInverseSimpsonCalculated,
      RawDifference =
        InverseSimpsonCalculated - InverseSimpson
    ),
  file.path(
    tab_dir,
    "inverse_simpson_formula_validation_by_sample.csv"
  ),
  row.names = FALSE
)

############################################################
## Figure 1 alternative: Pielou evenness calculated by formula
############################################################

all_patient_pielou_baseline <- meta %>%
  mutate(Patient = as.character(Patient)) %>%
  filter(day < 0) %>%
  group_by(Patient) %>%
  arrange(day, SampleID, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    Patient,
    BaselinePielouEvenness = PielouEvenness,
    BaselineDay = day,
    Recovery80PielouEvenness = 0.8 * PielouEvenness
  )

all_patient_pielou_recovery <- bind_rows(lapply(
  seq_len(nrow(all_patient_key)),
  function(i) {
    patient_row <- all_patient_key[i, ]
    recovery_threshold <-
      all_patient_pielou_baseline$Recovery80PielouEvenness[
        all_patient_pielou_baseline$Patient == patient_row$Patient
      ]
    post_days <- meta %>%
      filter(
        as.character(Patient) == patient_row$Patient,
        EventTime >= 1,
        EventTime <= 60
      ) %>%
      group_by(EventTime) %>%
      summarise(
        PielouEvenness = mean(PielouEvenness),
        .groups = "drop"
      ) %>%
      arrange(EventTime)

    sustained_day <- NA_real_
    if (nrow(post_days) >= 2) {
      crossings <- which(
        post_days$PielouEvenness >= recovery_threshold
      )
      for (crossing in crossings) {
        later <- which(
          post_days$EventTime > post_days$EventTime[crossing]
        )
        if (
          length(later) > 0 &&
            post_days$PielouEvenness[later[1]] >= recovery_threshold
        ) {
          sustained_day <- post_days$EventTime[crossing]
          break
        }
      }
    }

    tibble(
      Patient = patient_row$Patient,
      SustainedRecoveryDay = sustained_day,
      RecoveryLabel = if (is.na(sustained_day)) {
        "80% threshold not sustained"
      } else {
        paste0("80% threshold: +", sustained_day, " d")
      }
    )
  }
))

all_patient_pielou_trajectory <- meta %>%
  mutate(Patient = as.character(Patient)) %>%
  left_join(all_patient_key, by = c("Patient", "Arm", "ParticipantID")) %>%
  left_join(all_patient_pielou_baseline, by = "Patient") %>%
  left_join(all_patient_pielou_recovery, by = "Patient") %>%
  arrange(Arm, ArmOrder, day, SampleID)

make_patient_pielou_atlas <- function(
  arm_name,
  page_number,
  tag
) {
  plot_data <- all_patient_pielou_trajectory %>%
    filter(Arm == arm_name, AtlasPage == page_number)
  patient_levels <- all_patient_key %>%
    filter(Arm == arm_name, AtlasPage == page_number) %>%
    arrange(ArmOrder) %>%
    pull(PatientLabel)
  plot_data <- plot_data %>%
    mutate(PatientLabel = factor(PatientLabel, levels = patient_levels))
  line_data <- plot_data %>%
    mutate(
      TimeSegment = case_when(
        day < IndexDay ~ "Pre-index",
        day > IndexDay ~ "Post-index",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(TimeSegment))
  baseline_data <- plot_data %>%
    distinct(
      PatientLabel, BaselinePielouEvenness,
      Recovery80PielouEvenness
    )
  recovery_data <- plot_data %>%
    distinct(PatientLabel, RecoveryLabel)
  event_data <- all_patient_events %>%
    filter(Arm == arm_name, AtlasPage == page_number) %>%
    mutate(PatientLabel = factor(PatientLabel, levels = patient_levels))

  ggplot(
    plot_data,
    aes(
      x = day,
      y = PielouEvenness,
      group = PatientLabel
    )
  ) +
    annotate(
      "rect",
      xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
      fill = "grey70", alpha = 0.10
    ) +
    geom_vline(
      xintercept = 0, linewidth = 0.52,
      linetype = "22", color = "grey15"
    ) +
    geom_vline(
      data = event_data,
      aes(xintercept = EventDay, color = Event),
      inherit.aes = FALSE, linewidth = 0.58, linetype = "13"
    ) +
    geom_hline(
      data = baseline_data,
      aes(yintercept = BaselinePielouEvenness),
      inherit.aes = FALSE, linewidth = 0.50,
      linetype = "longdash", color = "grey30"
    ) +
    geom_hline(
      data = baseline_data,
      aes(yintercept = Recovery80PielouEvenness),
      inherit.aes = FALSE, linewidth = 0.45,
      linetype = "dotted", color = "grey55"
    ) +
    geom_line(
      data = line_data,
      aes(
        group = interaction(PatientLabel, TimeSegment, drop = TRUE)
      ),
      linewidth = 0.62, alpha = 0.82,
      color = unname(arm_colors[[arm_name]])
    ) +
    geom_point(
      data = plot_data %>% filter(day != IndexDay),
      shape = 21, fill = "white", stroke = 0.68, size = 2.0,
      color = unname(arm_colors[[arm_name]])
    ) +
    geom_point(
      data = plot_data %>% filter(day == IndexDay),
      shape = 23, fill = "grey88", stroke = 0.82, size = 2.7,
      color = "grey25"
    ) +
    geom_label(
      data = recovery_data,
      aes(
        x = atlas_x_limits[2] - 1,
        y = Inf,
        label = RecoveryLabel
      ),
      inherit.aes = FALSE,
      hjust = 1.0, vjust = 1.20,
      size = 2.10, linewidth = 0.16,
      label.padding = grid::unit(0.65, "mm"),
      fill = scales::alpha("white", 0.90), color = "grey20"
    ) +
    facet_wrap(vars(PatientLabel), ncol = 3) +
    scale_color_manual(
      values = atlas_event_colors,
      breaks = c(
        "Engraftment", "Randomization (no FMT)",
        "Auto-FMT / randomization"
      )
    ) +
    scale_x_continuous(
      limits = atlas_x_limits,
      breaks = c(-20, 0, 30, 60, 90),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = c(0, 0.8),
      breaks = seq(0, 0.8, by = 0.1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      tag = tag,
      title = paste0(
        arm_name, " patients, page ", page_number, " of 2"
      ),
      subtitle = paste0(
        "Pielou's J = Shannon entropy / log(observed richness); ",
        "grey diamonds mark index-day stools."
      ),
      x = "Day relative to allo-HSCT infusion",
      y = "Pielou's evenness",
      color = NULL
    ) +
    theme_nature(9.5) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.width = grid::unit(5.2, "mm"),
      panel.spacing = grid::unit(4.2, "mm"),
      strip.text = element_text(size = 9.2, face = "bold", hjust = 0),
      plot.title = element_text(size = 12.0),
      plot.subtitle = element_text(size = 8.6),
      plot.tag = element_text(size = 13.5),
      axis.text = element_text(size = 8.0),
      axis.title = element_text(size = 9.4)
    )
}

patient_pielou_atlas_specs <- tribble(
  ~Arm, ~Page, ~Tag, ~Filename, ~Height,
  "Control", 1, "a",
  "figure1a_all_patient_pielou_evenness_control_1", 6.0,
  "Control", 2, "b",
  "figure1b_all_patient_pielou_evenness_control_2", 6.0,
  "Auto-FMT", 1, "c",
  "figure1c_all_patient_pielou_evenness_autofmt_1", 7.9,
  "Auto-FMT", 2, "d",
  "figure1d_all_patient_pielou_evenness_autofmt_2", 7.9
)

patient_pielou_atlas_files <- unlist(lapply(
  seq_len(nrow(patient_pielou_atlas_specs)),
  function(i) {
    spec <- patient_pielou_atlas_specs[i, ]
    atlas_plot <- make_patient_pielou_atlas(
      spec$Arm, spec$Page, spec$Tag
    )
    save_plot(
      atlas_plot, spec$Filename,
      width = 7.2, height = spec$Height
    )
  }
))

write.csv(
  all_patient_pielou_trajectory %>%
    select(
      SampleID, sample_name, Patient, ParticipantID, PatientLabel,
      Arm, AtlasPage, day, engraftmentday, randomizationday,
      EventTime, PielouEvenness, BaselinePielouEvenness,
      Recovery80PielouEvenness
    ),
  file.path(
    tab_dir,
    "figure1_all_patient_pielou_evenness_trajectories.csv"
  ),
  row.names = FALSE
)
write.csv(
  all_patient_key %>%
    left_join(all_patient_pielou_baseline, by = "Patient") %>%
    left_join(all_patient_pielou_recovery, by = "Patient"),
  file.path(
    tab_dir,
    "figure1_all_patient_pielou_evenness_recovery_summary.csv"
  ),
  row.names = FALSE
)

############################################################
## Figure 1: longitudinal HRIC for patients 0009 and 0028
############################################################

focus_key <- tribble(
  ~Patient, ~PatientLabel, ~ClinicalLabel,
  "C1", "Patient 0009 (control)", "Randomization (no FMT)",
  "T4", "Patient 0028 (auto-FMT)", "Auto-FMT / randomization"
)

focus_meta <- meta %>%
  filter(as.character(Patient) %in% focus_key$Patient) %>%
  mutate(Patient = as.character(Patient)) %>%
  left_join(focus_key, by = "Patient") %>%
  mutate(
    PatientLabel = factor(PatientLabel, levels = focus_key$PatientLabel),
    Patient = factor(Patient, levels = focus_key$Patient)
  ) %>%
  arrange(Patient, day)

focus_baseline <- focus_meta %>%
  filter(day < 0) %>%
  group_by(PatientLabel) %>%
  arrange(day, SampleID, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    PatientLabel,
    BaselineAlpha = Alpha,
    BaselineDay = day,
    Recovery80 = 0.8 * Alpha
  )

focus_recovery_labels <- bind_rows(lapply(
  levels(focus_meta$PatientLabel),
  function(patient_label) {
    threshold <- focus_baseline$Recovery80[
      focus_baseline$PatientLabel == patient_label
    ]
    post_days <- focus_meta %>%
      filter(
        PatientLabel == patient_label,
        EventTime >= 1,
        EventTime <= 60
      ) %>%
      group_by(EventTime) %>%
      summarise(Alpha = mean(Alpha), .groups = "drop") %>%
      arrange(EventTime)

    sustained_day <- NA_real_
    if (nrow(post_days) >= 2) {
      hits <- which(post_days$Alpha >= threshold)
      for (i in hits) {
        later <- which(post_days$EventTime > post_days$EventTime[i])
        if (length(later) > 0 && post_days$Alpha[later[1]] >= threshold) {
          sustained_day <- post_days$EventTime[i]
          break
        }
      }
    }
    tibble(
      PatientLabel = factor(
        patient_label,
        levels = levels(focus_meta$PatientLabel)
      ),
      RecoveryLabel = if (is.na(sustained_day)) {
        "Sustained 80% recovery not observed by day 60"
      } else {
        paste0("Sustained 80% recovery: day +", sustained_day)
      }
    )
  }
))

compute_phase_hric <- function(patient_id, phase_name) {
  ids <- focus_meta %>%
    filter(as.character(Patient) == patient_id, as.character(Phase) == phase_name) %>%
    pull(SampleID)
  x <- counts[ids, , drop = FALSE]
  alpha <- HRIC::SHalpha(x)
  tibble(
    Patient = patient_id,
    Phase = phase_name,
    n = nrow(x),
    MeanAlpha = mean(alpha),
    Gamma = as.numeric(HRIC::SHgamma(x)),
    Beta = as.numeric(HRIC::SHbeta(x))
  )
}

focus_decomp <- bind_rows(lapply(focus_key$Patient, function(patient_id) {
  bind_rows(lapply(phase_levels, function(phase_name) {
    compute_phase_hric(patient_id, phase_name)
  }))
})) %>%
  left_join(focus_key, by = "Patient") %>%
  mutate(
    PatientLabel = factor(PatientLabel, levels = focus_key$PatientLabel),
    Phase = factor(Phase, levels = phase_levels)
  )

focus_decomp_long <- focus_decomp %>%
  select(PatientLabel, Phase, n, MeanAlpha, Beta, Gamma) %>%
  pivot_longer(
    cols = c(MeanAlpha, Beta),
    names_to = "Component",
    values_to = "Value"
  ) %>%
  mutate(
    Component = factor(
      recode(Component, MeanAlpha = "Mean alpha"),
      levels = c("Beta", "Mean alpha")
    )
  )

event_data <- focus_meta %>%
  distinct(PatientLabel, Patient, engraftmentday, randomizationday, ClinicalLabel) %>%
  pivot_longer(
    cols = c(engraftmentday, randomizationday),
    names_to = "EventType",
    values_to = "EventDay"
  ) %>%
  mutate(
    Event = case_when(
      EventType == "engraftmentday" ~ "Engraftment",
      Patient == "C1" ~ "Randomization (no FMT)",
      TRUE ~ "Auto-FMT / randomization"
    ),
    Event = factor(
      Event,
      levels = c(
        "Engraftment", "Randomization (no FMT)", "Auto-FMT / randomization"
      )
    )
  )

p1_trajectory <- ggplot(
  focus_meta,
  aes(x = day, y = Alpha, group = PatientLabel, color = Arm)
) +
  annotate(
    "rect",
    xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.10
  ) +
  geom_vline(
    xintercept = 0, linewidth = 0.55, linetype = "22", color = "grey20"
  ) +
  geom_vline(
    data = event_data,
    aes(xintercept = EventDay, color = Event),
    inherit.aes = FALSE, linewidth = 0.65, linetype = "13"
  ) +
  geom_hline(
    data = focus_baseline,
    aes(yintercept = BaselineAlpha),
    inherit.aes = FALSE, linewidth = 0.55,
    linetype = "longdash", color = "grey35"
  ) +
  geom_hline(
    data = focus_baseline,
    aes(yintercept = Recovery80),
    inherit.aes = FALSE, linewidth = 0.48,
    linetype = "dotted", color = "grey55"
  ) +
  geom_line(linewidth = 0.65, alpha = 0.78) +
  geom_point(shape = 21, fill = "white", stroke = 0.75, size = 2.4) +
  geom_point(
    data = focus_meta %>% filter(day == randomizationday),
    shape = 23, fill = "white", stroke = 0.9, size = 3.2,
    show.legend = FALSE
  ) +
  geom_label(
    data = focus_recovery_labels,
    aes(x = Inf, y = Inf, label = RecoveryLabel),
    inherit.aes = FALSE,
    hjust = 1.03, vjust = 1.18,
    size = 2.55, linewidth = 0.20,
    fill = scales::alpha("white", 0.88), color = "grey20"
  ) +
  facet_wrap(
    vars(PatientLabel), ncol = 1, scales = "free_x",
    strip.position = "top"
  ) +
  scale_color_manual(
    values = c(
      arm_colors,
      Engraftment = "#666666",
      `Randomization (no FMT)` = "#6A51A3",
      `Auto-FMT / randomization` = "#D97706"
    ),
    breaks = c(
      "Control", "Auto-FMT", "Engraftment",
      "Randomization (no FMT)", "Auto-FMT / randomization"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    x = "Day relative to allo-HSCT infusion",
    y = "HRIC alpha diversity",
    color = NULL
  ) +
  theme_nature(10) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.spacing.y = grid::unit(6, "mm")
  )

p1_bacteroidetes <- ggplot(
  focus_meta,
  aes(x = day, y = BacteroidetesPercent, group = PatientLabel, color = Arm)
) +
  annotate(
    "rect",
    xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.10
  ) +
  geom_vline(
    xintercept = 0, linewidth = 0.55, linetype = "22", color = "grey20"
  ) +
  geom_vline(
    data = event_data,
    aes(xintercept = EventDay, color = Event),
    inherit.aes = FALSE, linewidth = 0.65, linetype = "13"
  ) +
  geom_hline(
    yintercept = 0.1, linewidth = 0.65, linetype = "longdash",
    color = "grey25"
  ) +
  geom_line(linewidth = 0.65, alpha = 0.78) +
  geom_point(shape = 21, fill = "white", stroke = 0.75, size = 2.2) +
  facet_wrap(
    vars(PatientLabel), ncol = 1, scales = "free_x",
    strip.position = "top"
  ) +
  scale_color_manual(
    values = c(
      arm_colors,
      Engraftment = "#666666",
      `Randomization (no FMT)` = "#6A51A3",
      `Auto-FMT / randomization` = "#D97706"
    )
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.001),
    breaks = c(0, 0.01, 0.1, 1, 10, 100),
    labels = function(x) ifelse(x == 0, "0", scales::number(x, accuracy = 0.01)),
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  labs(
    title = "Bacteroidetes abundance around trial eligibility",
    subtitle = paste(
      "Dashed line: trial threshold of 0.1% by qPCR;",
      "points: descriptive 16S relative abundance."
    ),
    x = "Day relative to allo-HSCT infusion",
    y = "Bacteroidetes relative abundance (%)"
  ) +
  theme_nature(9.5) +
  theme(
    legend.position = "none",
    panel.spacing.y = grid::unit(5, "mm")
  ) +
  guides(color = "none")

p1_decomp <- ggplot(
focus_decomp_long,
  aes(x = Phase, y = Value, fill = Component)
) +
  geom_col(width = 0.66, color = "white", linewidth = 0.30) +
  geom_point(
    data = focus_decomp,
    aes(x = Phase, y = Gamma),
    inherit.aes = FALSE, shape = 23, size = 2.2,
    fill = "white", color = "grey10", stroke = 0.55
  ) +
  geom_text(
    data = focus_decomp,
    aes(
      x = Phase, y = Gamma,
      label = sprintf("atop(gamma == %.3f, italic(n) == %d)", Gamma, n)
    ),
    inherit.aes = FALSE, vjust = -0.55, size = 2.8,
    parse = TRUE
  ) +
  facet_wrap(vars(PatientLabel), ncol = 1, strip.position = "top") +
  scale_fill_manual(
    values = component_colors,
    breaks = c("Mean alpha", "Beta"),
    labels = c("Mean alpha", "Regional beta")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.20)),
    limits = c(0, max(focus_decomp$Gamma) * 1.25)
  ) +
  labs(
    x = NULL,
    y = "HRIC diversity",
    fill = NULL
  ) +
  theme_nature(10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 22, hjust = 1),
    panel.spacing.y = grid::unit(6, "mm")
  )

focus_taxon_long <- as_tibble(
  taxon_alpha_deficit[focus_meta$SampleID, , drop = FALSE],
  rownames = "SampleID",
  .name_repair = "minimal"
) %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "Taxon",
    values_to = "Contribution"
  ) %>%
  left_join(
    focus_meta %>%
      group_by(PatientLabel) %>%
      arrange(day, SampleID, .by_group = TRUE) %>%
      mutate(SampleOrder = row_number()) %>%
      ungroup() %>%
      select(SampleID, PatientLabel, day, SampleOrder),
    by = "SampleID"
  )

focus_sample_top_drivers <- focus_taxon_long %>%
  group_by(SampleID, PatientLabel, day, SampleOrder) %>%
  mutate(
    AlphaDeficit = sum(Contribution),
    FractionOfAlphaDeficit = Contribution / AlphaDeficit
  ) %>%
  arrange(desc(Contribution), Taxon, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  mutate(DriverRank = row_number()) %>%
  ungroup() %>%
  select(
    SampleID, PatientLabel, day, SampleOrder, DriverRank, Taxon,
    Contribution, FractionOfAlphaDeficit, AlphaDeficit
  )

focus_top_taxa <- focus_taxon_long %>%
  group_by(Taxon) %>%
  summarise(TotalContribution = sum(Contribution), .groups = "drop") %>%
  arrange(desc(TotalContribution), Taxon) %>%
  slice_head(n = 8)

focus_taxon_display <- focus_taxon_long %>%
  mutate(
    Taxon = ifelse(Taxon %in% focus_top_taxa$Taxon, Taxon, "Other taxa")
  ) %>%
  group_by(SampleID, PatientLabel, day, SampleOrder, Taxon) %>%
  summarise(Contribution = sum(Contribution), .groups = "drop") %>%
  complete(
    nesting(SampleID, PatientLabel, day, SampleOrder),
    Taxon = c(focus_top_taxa$Taxon, "Other taxa"),
    fill = list(Contribution = 0)
  )

focus_taxon_levels <- c(focus_top_taxa$Taxon, "Other taxa")
focus_alluvial <- make_alluvial_ribbons(
  focus_taxon_display,
  focus_taxon_levels
)

taxon_colors <- setNames(
  c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#6A3D9A", "#8C6D31",
    "#B8B8B8"
  ),
  focus_taxon_levels
)

p1_taxa <- ggplot() +
  geom_polygon(
    data = focus_alluvial$ribbons,
    aes(x = x, y = y, group = RibbonID, fill = Taxon),
    color = NA, alpha = 0.92
  ) +
  geom_col(
    data = focus_alluvial$ordered,
    aes(x = SampleOrder, y = Contribution, fill = Taxon),
    width = 0.16, color = "white", linewidth = 0.08
  ) +
  facet_wrap(
    vars(PatientLabel), ncol = 1, scales = "free_x",
    strip.position = "top"
  ) +
  scale_fill_manual(values = taxon_colors, breaks = focus_taxon_levels) +
  scale_x_continuous(
    breaks = seq_len(max(focus_taxon_display$SampleOrder)),
    labels = function(x) paste0("S", x),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(
    x = "Sequential stool sample",
    y = expression("Taxon contribution to " * (1 - alpha)),
    fill = NULL
  ) +
  theme_nature(9.5) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "bottom",
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(5, "mm"),
    panel.spacing.y = grid::unit(6, "mm")
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

figure1 <- ((p1_trajectory + p1_decomp) /
  p1_bacteroidetes /
  p1_taxa) +
  plot_layout(
    widths = c(1.65, 1),
    heights = c(1.18, 0.74, 1.08),
    guides = "keep"
  ) +
  plot_annotation(
    title = "Patient-level HRIC diversity and taxon contributions across allo-HSCT",
    subtitle = paste0(
      "Horizontal lines mark personal baseline (dashed) and 80% recovery (dotted);",
      " taxon ribbons decompose 1 - alpha.\n",
      "Bacteroidetes provides eligibility context; ",
      "vertical dotted lines mark engraftment and randomization/auto-FMT; ",
      "qPCR and 16S abundance are shown as distinct measurements."
    ),
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(size = 10.2, color = "grey25", hjust = 0),
      plot.tag = element_text(face = "bold", size = 14)
    )
  ) &
  theme(legend.position = "bottom")

figure1_files <- save_plot(
  figure1, "figure1_fmt_hric_patients_0009_0028",
  width = 13.8, height = 14.8
)

write.csv(
  focus_meta %>%
    select(
      SampleID, sample_name, Patient, PatientLabel, Arm, day,
      engraftmentday, randomizationday, EventTime, Phase, Alpha,
      BacteroidetesPercent
    ),
  file.path(tab_dir, "figure1_patient_alpha_trajectories.csv"),
  row.names = FALSE
)
write.csv(
  focus_meta %>%
    select(
      SampleID, sample_name, Patient, PatientLabel, Arm, day,
      engraftmentday, randomizationday, EventTime,
      BacteroidetesPercent
    ),
  file.path(tab_dir, "figure1_patient_bacteroidetes_trajectory.csv"),
  row.names = FALSE
)
write.csv(
  focus_decomp,
  file.path(tab_dir, "figure1_patient_phase_alpha_beta_gamma.csv"),
  row.names = FALSE
)
write.csv(
  focus_taxon_display,
  file.path(tab_dir, "figure1_patient_taxon_alpha_deficit.csv"),
  row.names = FALSE
)
write.csv(
  focus_sample_top_drivers,
  file.path(tab_dir, "figure1_sample_top_taxon_drivers.csv"),
  row.names = FALSE
)
write.csv(
  focus_baseline %>%
    left_join(focus_recovery_labels, by = "PatientLabel"),
  file.path(tab_dir, "figure1_personal_recovery_thresholds.csv"),
  row.names = FALSE
)

############################################################
## Figure 2: randomized-arm comparison around the index day
############################################################

# No stool was collected exactly on the intervention day for every patient.
# One sample per patient is selected in each transparent, prespecified window.
landmark_specs <- tribble(
  ~Landmark, ~WindowLabel, ~Lower, ~Upper, ~Target,
  "Pre-index", "-28 to -1 d", -28, -1, -1,
  "Early post-index", "1 to 14 d", 1, 14, 1,
  "Later post-index", "15 to 60 d", 15, 60, 30
) %>%
  mutate(
    Landmark = factor(
      Landmark,
      levels = c("Pre-index", "Early post-index", "Later post-index")
    )
  )

select_landmark <- function(spec_row) {
  candidates <- meta %>%
    filter(
      EventTime >= spec_row$Lower,
      EventTime <= spec_row$Upper
    ) %>%
    mutate(
      Landmark = spec_row$Landmark,
      TargetDistance = abs(EventTime - spec_row$Target)
    ) %>%
    group_by(Patient) %>%
    arrange(TargetDistance, abs(EventTime), day, SampleID, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup()
  candidates
}

landmarks <- bind_rows(lapply(seq_len(nrow(landmark_specs)), function(i) {
  select_landmark(landmark_specs[i, ])
})) %>%
  mutate(
    Landmark = factor(
      as.character(Landmark),
      levels = levels(landmark_specs$Landmark)
    )
  )

landmark_summary <- landmarks %>%
  group_by(Landmark, Arm) %>%
  summarise(
    n = n(),
    mean_alpha = mean(Alpha),
    alpha_sd = sd(Alpha),
    .groups = "drop"
  ) %>%
  group_by(Landmark, Arm) %>%
  group_modify(~ bind_cols(.x, mean_ci(
    landmarks$Alpha[
      landmarks$Landmark == .y$Landmark &
        landmarks$Arm == .y$Arm
    ]
  ))) %>%
  ungroup()

compute_group_decomp <- function(landmark_name, arm_name) {
  d <- landmarks %>%
    filter(
      as.character(Landmark) == landmark_name,
      as.character(Arm) == arm_name
    )
  x <- counts[d$SampleID, , drop = FALSE]
  tibble(
    Landmark = landmark_name,
    Arm = arm_name,
    n = nrow(x),
    MeanAlpha = mean(HRIC::SHalpha(x)),
    Gamma = as.numeric(HRIC::SHgamma(x)),
    Beta = as.numeric(HRIC::SHbeta(x))
  )
}

group_decomp <- bind_rows(lapply(levels(landmarks$Landmark), function(landmark_name) {
  bind_rows(lapply(levels(meta$Arm), function(arm_name) {
    compute_group_decomp(landmark_name, arm_name)
  }))
})) %>%
  mutate(
    Landmark = factor(Landmark, levels = levels(landmarks$Landmark)),
    Arm = factor(Arm, levels = levels(meta$Arm))
  )

group_decomp_long <- group_decomp %>%
  pivot_longer(
    cols = c(MeanAlpha, Beta),
    names_to = "Component",
    values_to = "Value"
  ) %>%
  mutate(
    Component = factor(
      recode(Component, MeanAlpha = "Mean alpha"),
      levels = c("Beta", "Mean alpha")
    )
  )

delta_by_landmark <- bind_rows(lapply(levels(landmarks$Landmark), function(landmark_name) {
  d <- landmarks %>% filter(as.character(Landmark) == landmark_name)
  tibble(
    Landmark = landmark_name,
    Delta = as.numeric(HRIC::SHdelta(counts[d$SampleID, , drop = FALSE], d$Arm))
  )
})) %>%
  mutate(Landmark = factor(Landmark, levels = levels(landmarks$Landmark)))

############################################################
## Time-adjusted inference
############################################################

# Alpha model: a common smooth trend over allo-HSCT day, a post-index change,
# and the additional post-index change attributable to randomized auto-FMT.
alpha_model <- lmer(
  Alpha ~ Arm + ns(day, df = 4) + PostIndex + TreatmentPost + (1 | Patient),
  data = meta,
  REML = FALSE
)
alpha_coef <- summary(alpha_model)$coefficients["TreatmentPost", ]
alpha_effect <- tibble(
  Metric = "Treatment-specific post-index effect",
  Estimate = unname(alpha_coef["Estimate"]),
  SE = unname(alpha_coef["Std. Error"]),
  Lower = Estimate - 1.96 * SE,
  Upper = Estimate + 1.96 * SE,
  P = 2 * stats::pnorm(abs(Estimate / SE), lower.tail = FALSE)
)

# HRIC-coordinate comparison: remove the common nonlinear allo-HSCT-day trend,
# calculate each patient's early-post minus pre-index coordinate change, and
# compare those change vectors between randomized arms. The permutation keeps
# every patient's complete longitudinal change vector intact.
time_design <- model.matrix(~ ns(day, df = 4), data = meta)
Z_time_residual <- qr.resid(qr(time_design), Z)
rownames(Z_time_residual) <- meta$SampleID

paired_landmarks <- landmarks %>%
  filter(
    Landmark %in% c("Pre-index", "Early post-index")
  ) %>%
  select(Patient, Arm, Landmark, SampleID) %>%
  pivot_wider(names_from = Landmark, values_from = SampleID) %>%
  filter(
    !is.na(`Pre-index`),
    !is.na(`Early post-index`)
  )

change_matrix <- t(vapply(seq_len(nrow(paired_landmarks)), function(i) {
  post_id <- paired_landmarks$`Early post-index`[i]
  pre_id <- paired_landmarks$`Pre-index`[i]
  Z_time_residual[post_id, ] - Z_time_residual[pre_id, ]
}, numeric(ncol(Z_time_residual))))
rownames(change_matrix) <- as.character(paired_landmarks$Patient)

between_change_stat <- function(group, gram_matrix) {
  group <- factor(group, levels = c("Control", "Auto-FMT"))
  n_total <- length(group)
  n_control <- sum(group == "Control")
  n_treatment <- sum(group == "Auto-FMT")
  weights <- ifelse(
    group == "Auto-FMT",
    1 / n_treatment,
    -1 / n_control
  )
  mean_difference_sq <- as.numeric(t(weights) %*% gram_matrix %*% weights)
  (n_control * n_treatment / n_total^2) * mean_difference_sq
}

change_gram <- tcrossprod(change_matrix)
observed_change_stat <- between_change_stat(paired_landmarks$Arm, change_gram)
n_permutations <- 4999L
permuted_change_stats <- replicate(
  n_permutations,
  between_change_stat(sample(paired_landmarks$Arm), change_gram)
)
change_p <- (1 + sum(permuted_change_stats >= observed_change_stat)) /
  (n_permutations + 1)

change_test <- tibble(
  n_patients = nrow(paired_landmarks),
  statistic = observed_change_stat,
  permutations = n_permutations,
  p_value = change_p,
  definition = paste(
    "Between-arm separation of patient-level early-post minus pre-index",
    "HRIC-coordinate changes after removing a nonlinear allo-HSCT-day trend"
  )
)

############################################################
## Baseline-adjusted recovery speed
############################################################

# The earliest observed pre-HSCT stool is the closest available proxy for the
# stored donor/healthy baseline. The nearest pre-index stool captures damage at
# the moment immediately preceding the randomized intervention.
baseline_covariates <- meta %>%
  filter(day < 0) %>%
  group_by(Patient) %>%
  arrange(day, SampleID, .by_group = TRUE) %>%
  summarise(
    BaselineSampleID = first(SampleID),
    BaselineDay = first(day),
    BaselineAlpha = first(Alpha),
    BaselineAlphaMean = mean(Alpha),
    .groups = "drop"
  )

preindex_covariates <- landmarks %>%
  filter(Landmark == "Pre-index") %>%
  transmute(
    Patient,
    PreIndexSampleID = SampleID,
    PreIndexAlpha = Alpha
  )

patient_covariates <- meta %>%
  distinct(Patient, Arm, randomizationday) %>%
  left_join(baseline_covariates, by = "Patient") %>%
  left_join(preindex_covariates, by = "Patient")

baseline_ids_by_sample <- patient_covariates$BaselineSampleID[
  match(meta$Patient, patient_covariates$Patient)
]
meta$HRICBaselineDistance <- sqrt(rowSums(
  (Z - Z[baseline_ids_by_sample, , drop = FALSE])^2
))

preindex_distance <- meta %>%
  select(Patient, SampleID, HRICBaselineDistance) %>%
  inner_join(
    patient_covariates %>% select(Patient, PreIndexSampleID),
    by = "Patient"
  ) %>%
  filter(SampleID == PreIndexSampleID) %>%
  transmute(
    Patient,
    PreIndexDistance = HRICBaselineDistance
  )

patient_covariates <- patient_covariates %>%
  left_join(preindex_distance, by = "Patient")

post_recovery <- meta %>%
  filter(EventTime >= 1, EventTime <= 60) %>%
  left_join(
    patient_covariates %>%
      select(
        Patient, BaselineAlpha, BaselineAlphaMean,
        PreIndexAlpha, PreIndexDistance
      ),
    by = "Patient"
  ) %>%
  mutate(
    BaselineAlphaC = BaselineAlpha - mean(BaselineAlpha, na.rm = TRUE),
    BaselineAlphaMeanC = BaselineAlphaMean -
      mean(BaselineAlphaMean, na.rm = TRUE),
    PreIndexAlphaC = PreIndexAlpha - mean(PreIndexAlpha, na.rm = TRUE),
    PreIndexDistanceC = PreIndexDistance -
      mean(PreIndexDistance, na.rm = TRUE),
    IndexDayC = randomizationday - mean(randomizationday, na.rm = TRUE),
    EarlyTime = pmin(pmax(EventTime - 1, 0), 13),
    LateTime = pmax(EventTime - 14, 0)
  )

# Sensitivity data include index-day stools. Their within-day order relative to
# auto-FMT is unknown, so they are excluded from the primary estimand.
post_recovery_day0 <- meta %>%
  filter(EventTime >= 0, EventTime <= 60) %>%
  left_join(
    patient_covariates %>%
      select(
        Patient, BaselineAlpha, BaselineAlphaMean,
        PreIndexAlpha, PreIndexDistance
      ),
    by = "Patient"
  ) %>%
  mutate(
    BaselineAlphaC = BaselineAlpha - mean(BaselineAlpha, na.rm = TRUE),
    BaselineAlphaMeanC = BaselineAlphaMean -
      mean(BaselineAlphaMean, na.rm = TRUE),
    PreIndexAlphaC = PreIndexAlpha - mean(PreIndexAlpha, na.rm = TRUE),
    PreIndexDistanceC = PreIndexDistance -
      mean(PreIndexDistance, na.rm = TRUE),
    IndexDayC = randomizationday - mean(randomizationday, na.rm = TRUE),
    EarlyTime = pmin(EventTime, 14),
    LateTime = pmax(EventTime - 14, 0)
  )

alpha_recovery_model <- lmer(
  Alpha ~ BaselineAlphaC + PreIndexAlphaC + IndexDayC +
    Arm * (EarlyTime + LateTime) + (1 | Patient),
  data = post_recovery,
  REML = FALSE
)

distance_recovery_model <- lm(
  HRICBaselineDistance ~ BaselineAlphaC + PreIndexDistanceC + IndexDayC +
    Arm * (EarlyTime + LateTime),
  data = post_recovery
)
distance_recovery_covariance <- sandwich::vcovCL(
  distance_recovery_model,
  cluster = post_recovery$Patient,
  type = "HC2"
)

# Sensitivity analysis: replace the earliest baseline alpha with the mean of all
# available pre-HSCT alpha measurements.
alpha_recovery_model_mean_baseline <- lmer(
  Alpha ~ BaselineAlphaMeanC + PreIndexAlphaC + IndexDayC +
    Arm * (EarlyTime + LateTime) + (1 | Patient),
  data = post_recovery,
  REML = FALSE
)

alpha_recovery_model_day0 <- lmer(
  Alpha ~ BaselineAlphaC + PreIndexAlphaC + IndexDayC +
    Arm * (EarlyTime + LateTime) + (1 | Patient),
  data = post_recovery_day0,
  REML = FALSE
)

distance_recovery_model_day0 <- lm(
  HRICBaselineDistance ~ BaselineAlphaC + PreIndexDistanceC + IndexDayC +
    Arm * (EarlyTime + LateTime),
  data = post_recovery_day0
)
distance_recovery_covariance_day0 <- sandwich::vcovCL(
  distance_recovery_model_day0,
  cluster = post_recovery_day0$Patient,
  type = "HC2"
)

recovery_grid <- expand_grid(
  EventTime = seq(1, 60, by = 1),
  Arm = factor(c("Control", "Auto-FMT"), levels = levels(meta$Arm))
) %>%
  mutate(
    BaselineAlphaC = 0,
    BaselineAlphaMeanC = 0,
    PreIndexAlphaC = 0,
    PreIndexDistanceC = 0,
    IndexDayC = 0,
    EarlyTime = pmin(pmax(EventTime - 1, 0), 13),
    LateTime = pmax(EventTime - 14, 0)
  )

alpha_recovery_curve <- predict_fixed_curve(
  alpha_recovery_model,
  recovery_grid
)
distance_recovery_curve <- predict_fixed_curve(
  distance_recovery_model,
  recovery_grid,
  covariance = distance_recovery_covariance
)

extract_speed_terms <- function(model, outcome, covariance = NULL) {
  is_mixed <- inherits(model, "merMod")
  coefficients <- if (is_mixed) lme4::fixef(model) else coef(model)
  if (is.null(covariance)) covariance <- vcov(model)
  standard_errors <- sqrt(diag(covariance))
  coefficient_names <- names(coefficients)
  terms <- c("EarlyTime", "LateTime")
  bind_rows(lapply(terms, function(time_term) {
    name <- coefficient_names[
      grepl("ArmAuto-FMT", coefficient_names, fixed = TRUE) &
        grepl(time_term, coefficient_names, fixed = TRUE)
    ]
    if (length(name) != 1) stop("Cannot identify ", time_term, " interaction")
    estimate <- coefficients[name]
    se <- standard_errors[name]
    tibble(
      Outcome = outcome,
      Period = ifelse(time_term == "EarlyTime", "Day 1-14", "Day 15-60"),
      EstimatePerDay = estimate,
      SE = se,
      Lower = estimate - 1.96 * se,
      Upper = estimate + 1.96 * se,
      P = 2 * pnorm(abs(estimate / se), lower.tail = FALSE)
    )
  }))
}

recovery_speed_effects <- bind_rows(
  extract_speed_terms(alpha_recovery_model, "HRIC alpha diversity"),
  extract_speed_terms(
    distance_recovery_model,
    "HRIC distance to baseline",
    covariance = distance_recovery_covariance
  ),
  extract_speed_terms(
    alpha_recovery_model_mean_baseline,
    "HRIC alpha diversity (mean-baseline sensitivity)"
  )
)

make_arm_contrast <- function(
  model,
  day_value,
  outcome,
  covariance = NULL,
  include_day0 = FALSE
) {
  common <- tibble(
    EventTime = day_value,
    BaselineAlphaC = 0,
    BaselineAlphaMeanC = 0,
    PreIndexAlphaC = 0,
    PreIndexDistanceC = 0,
    IndexDayC = 0,
    EarlyTime = if (include_day0) {
      pmin(pmax(day_value, 0), 14)
    } else {
      pmin(pmax(day_value - 1, 0), 13)
    },
    LateTime = pmax(day_value - 14, 0)
  )
  treated <- common %>%
    mutate(Arm = factor("Auto-FMT", levels = levels(meta$Arm)))
  control <- common %>%
    mutate(Arm = factor("Control", levels = levels(meta$Arm)))
  fixed_contrast(
    model,
    treated,
    control,
    covariance = covariance
  ) %>%
    mutate(Outcome = outcome, EventTime = day_value, .before = 1)
}

recovery_arm_contrasts <- bind_rows(lapply(c(14, 30, 60), function(day_value) {
  bind_rows(
    make_arm_contrast(
      alpha_recovery_model, day_value, "HRIC alpha diversity"
    ),
    make_arm_contrast(
      distance_recovery_model,
      day_value,
      "HRIC distance to baseline",
      covariance = distance_recovery_covariance
    )
  )
}))

recovery_sensitivity_contrasts <- bind_rows(lapply(c(14, 30, 60), function(day_value) {
  bind_rows(
    make_arm_contrast(
      alpha_recovery_model,
      day_value,
      "HRIC alpha diversity"
    ) %>%
      mutate(Analysis = "Primary: earliest pre-HSCT baseline; exclude day 0"),
    make_arm_contrast(
      alpha_recovery_model_mean_baseline,
      day_value,
      "HRIC alpha diversity"
    ) %>%
      mutate(Analysis = "Sensitivity: mean pre-HSCT baseline; exclude day 0"),
    make_arm_contrast(
      alpha_recovery_model_day0,
      day_value,
      "HRIC alpha diversity",
      include_day0 = TRUE
    ) %>%
      mutate(Analysis = "Sensitivity: earliest pre-HSCT baseline; include day 0"),
    make_arm_contrast(
      distance_recovery_model,
      day_value,
      "HRIC distance to baseline",
      covariance = distance_recovery_covariance
    ) %>%
      mutate(Analysis = "Primary: earliest pre-HSCT baseline; exclude day 0"),
    make_arm_contrast(
      distance_recovery_model_day0,
      day_value,
      "HRIC distance to baseline",
      covariance = distance_recovery_covariance_day0,
      include_day0 = TRUE
    ) %>%
      mutate(Analysis = "Sensitivity: earliest pre-HSCT baseline; include day 0")
  )
})) %>%
  relocate(Analysis, .after = Outcome)

# The primary recovery estimand is the time-averaged adjusted arm difference
# over days 1-30. This directly targets early recovery and avoids selecting a
# single favorable sampling day. Days 1-60 is retained as a secondary window.
recovery_design_matrix <- function(model, arm_name, days) {
  is_mixed <- inherits(model, "merMod")
  fixed_formula <- if (is_mixed) {
    reformulas::nobars(formula(model))
  } else {
    formula(model)
  }
  fixed_terms <- delete.response(terms(fixed_formula))
  newdata <- tibble(
    EventTime = days,
    Arm = factor(arm_name, levels = levels(meta$Arm)),
    BaselineAlphaC = 0,
    BaselineAlphaMeanC = 0,
    PreIndexAlphaC = 0,
    PreIndexDistanceC = 0,
    IndexDayC = 0,
    EarlyTime = pmin(pmax(days - 1, 0), 13),
    LateTime = pmax(days - 14, 0)
  )
  coefficient_names <- names(if (is_mixed) lme4::fixef(model) else coef(model))
  model.matrix(fixed_terms, newdata)[
    , coefficient_names, drop = FALSE
  ]
}

make_average_arm_contrast <- function(
  model,
  end_day,
  outcome,
  covariance = NULL
) {
  days <- seq(1, end_day, by = 1)
  weights <- if (length(days) == 1) {
    1
  } else {
    c(0.5, rep(1, length(days) - 2), 0.5) / (end_day - 1)
  }
  treatment_design <- recovery_design_matrix(model, "Auto-FMT", days)
  control_design <- recovery_design_matrix(model, "Control", days)
  contrast <- colSums(
    (treatment_design - control_design) * weights
  )
  names(contrast) <- colnames(treatment_design)
  fixed_linear_contrast(
    model,
    contrast,
    covariance = covariance
  ) %>%
    mutate(
      Outcome = outcome,
      WindowEnd = end_day,
      Window = paste0("Day 1-", end_day),
      Estimand = "Time-averaged adjusted Auto-FMT - control difference",
      .before = 1
    )
}

recovery_average_effects_analytic <- bind_rows(lapply(c(30, 60), function(end_day) {
  bind_rows(
    make_average_arm_contrast(
      alpha_recovery_model,
      end_day,
      "HRIC alpha diversity"
    ),
    make_average_arm_contrast(
      distance_recovery_model,
      end_day,
      "HRIC distance to personal baseline",
      covariance = distance_recovery_covariance
    )
  )
}))

resample_patients_within_arm <- function(data) {
  bind_rows(lapply(levels(data$Arm), function(arm_name) {
    arm_data <- data %>% filter(Arm == arm_name)
    patient_ids <- unique(as.character(arm_data$Patient))
    sampled_ids <- sample(patient_ids, length(patient_ids), replace = TRUE)
    bind_rows(lapply(seq_along(sampled_ids), function(i) {
      arm_data %>%
        filter(as.character(Patient) == sampled_ids[i]) %>%
        mutate(Patient = paste0(arm_name, "_bootstrap_", i))
    }))
  })) %>%
    mutate(
      Patient = factor(Patient),
      Arm = factor(as.character(Arm), levels = levels(meta$Arm))
    )
}

bootstrap_recovery_effects <- function(data, iterations) {
  draws <- vector("list", iterations)
  completed <- 0L
  attempts <- 0L
  maximum_attempts <- max(iterations * 3L, iterations + 20L)

  while (completed < iterations && attempts < maximum_attempts) {
    attempts <- attempts + 1L
    bootstrap_data <- resample_patients_within_arm(data)
    fitted <- tryCatch(
      suppressMessages(suppressWarnings({
        alpha_fit <- lmer(
          Alpha ~ BaselineAlphaC + PreIndexAlphaC + IndexDayC +
            Arm * (EarlyTime + LateTime) + (1 | Patient),
          data = bootstrap_data,
          REML = FALSE
        )
        distance_fit <- lm(
          HRICBaselineDistance ~ BaselineAlphaC + PreIndexDistanceC +
            IndexDayC + Arm * (EarlyTime + LateTime),
          data = bootstrap_data
        )
        list(alpha = alpha_fit, distance = distance_fit)
      })),
      error = function(e) NULL
    )
    if (is.null(fitted)) next

    completed <- completed + 1L
    draws[[completed]] <- bind_rows(lapply(c(30, 60), function(end_day) {
      bind_rows(
        make_average_arm_contrast(
          fitted$alpha,
          end_day,
          "HRIC alpha diversity"
        ),
        make_average_arm_contrast(
          fitted$distance,
          end_day,
          "HRIC distance to personal baseline"
        )
      )
    })) %>%
      transmute(
        BootstrapIteration = completed,
        Outcome,
        WindowEnd,
        Window,
        Estimate
      )
  }

  if (completed < iterations) {
    warning(
      "Completed ", completed, " of ", iterations,
      " requested patient bootstrap iterations."
    )
  }
  bind_rows(draws[seq_len(completed)])
}

bootstrap_iterations <- suppressWarnings(as.integer(
  Sys.getenv("HRIC_BOOTSTRAPS", unset = "999")
))
if (!is.finite(bootstrap_iterations) || bootstrap_iterations < 19) {
  bootstrap_iterations <- 999L
}
set.seed(20260723)
recovery_bootstrap_draws <- bootstrap_recovery_effects(
  post_recovery,
  bootstrap_iterations
)

recovery_bootstrap_summary <- recovery_bootstrap_draws %>%
  group_by(Outcome, WindowEnd, Window) %>%
  summarise(
    BootstrapIterations = n(),
    BootstrapSE = sd(Estimate),
    BootstrapLower = quantile(Estimate, 0.025, names = FALSE),
    BootstrapUpper = quantile(Estimate, 0.975, names = FALSE),
    BootstrapP = min(
      1,
      2 * min(
        (1 + sum(Estimate <= 0)) / (n() + 1),
        (1 + sum(Estimate >= 0)) / (n() + 1)
      )
    ),
    .groups = "drop"
  )

recovery_average_effects <- recovery_average_effects_analytic %>%
  left_join(
    recovery_bootstrap_summary,
    by = c("Outcome", "WindowEnd", "Window")
  )

# Sensitivity analysis: estimate the pre-index arm-specific trend and then
# subtract the adjusted arm difference at day -1 from the post-index contrast.
# This checks whether chance imbalance in the pre-FMT trajectory explains the
# primary early-window result.
event_recovery <- meta %>%
  filter(
    (EventTime >= -28 & EventTime <= -1) |
      (EventTime >= 1 & EventTime <= 60)
  ) %>%
  left_join(
    patient_covariates %>%
      select(Patient, BaselineAlpha),
    by = "Patient"
  ) %>%
  mutate(
    BaselineAlphaC = BaselineAlpha - mean(BaselineAlpha, na.rm = TRUE),
    IndexDayC = randomizationday - mean(randomizationday, na.rm = TRUE),
    PreTime = ifelse(EventTime <= -1, EventTime + 1, 0),
    Post = as.integer(EventTime >= 1),
    EarlyTime = ifelse(
      EventTime >= 1,
      pmin(pmax(EventTime - 1, 0), 13),
      0
    ),
    LateTime = ifelse(EventTime >= 1, pmax(EventTime - 14, 0), 0)
  )

alpha_event_model <- lmer(
  Alpha ~ BaselineAlphaC + IndexDayC +
    Arm * (PreTime + Post + EarlyTime + LateTime) + (1 | Patient),
  data = event_recovery,
  REML = FALSE
)
distance_event_model <- lm(
  HRICBaselineDistance ~ BaselineAlphaC + IndexDayC +
    Arm * (PreTime + Post + EarlyTime + LateTime),
  data = event_recovery
)
distance_event_covariance <- sandwich::vcovCL(
  distance_event_model,
  cluster = event_recovery$Patient,
  type = "HC2"
)

event_design_matrix <- function(model, arm_name, event_times) {
  is_mixed <- inherits(model, "merMod")
  fixed_formula <- if (is_mixed) {
    reformulas::nobars(formula(model))
  } else {
    formula(model)
  }
  fixed_terms <- delete.response(terms(fixed_formula))
  newdata <- tibble(
    EventTime = event_times,
    Arm = factor(arm_name, levels = levels(meta$Arm)),
    BaselineAlphaC = 0,
    IndexDayC = 0,
    PreTime = ifelse(event_times <= -1, event_times + 1, 0),
    Post = as.integer(event_times >= 1),
    EarlyTime = ifelse(
      event_times >= 1,
      pmin(pmax(event_times - 1, 0), 13),
      0
    ),
    LateTime = ifelse(
      event_times >= 1,
      pmax(event_times - 14, 0),
      0
    )
  )
  coefficient_names <- names(if (is_mixed) lme4::fixef(model) else coef(model))
  model.matrix(fixed_terms, newdata)[
    , coefficient_names, drop = FALSE
  ]
}

make_event_did_average <- function(
  model,
  end_day,
  outcome,
  covariance = NULL
) {
  post_days <- seq(1, end_day)
  weights <- c(0.5, rep(1, length(post_days) - 2), 0.5) /
    (end_day - 1)
  post_difference <- event_design_matrix(
    model, "Auto-FMT", post_days
  ) - event_design_matrix(model, "Control", post_days)
  preindex_difference <- event_design_matrix(
    model, "Auto-FMT", -1
  ) - event_design_matrix(model, "Control", -1)
  contrast <- colSums(post_difference * weights) -
    as.numeric(preindex_difference)
  names(contrast) <- colnames(post_difference)
  fixed_linear_contrast(
    model,
    contrast,
    covariance = covariance
  ) %>%
    mutate(
      Outcome = outcome,
      WindowEnd = end_day,
      Window = paste0("Day 1-", end_day),
      Estimand = paste(
        "Time-averaged adjusted arm difference minus",
        "the adjusted arm difference at day -1"
      ),
      .before = 1
    )
}

event_did_average_effects <- bind_rows(lapply(c(30, 60), function(end_day) {
  bind_rows(
    make_event_did_average(
      alpha_event_model,
      end_day,
      "HRIC alpha diversity"
    ),
    make_event_did_average(
      distance_event_model,
      end_day,
      "HRIC distance to personal baseline",
      covariance = distance_event_covariance
    )
  )
}))

extract_named_interaction <- function(model, pattern, outcome, covariance = NULL) {
  coefficients <- if (inherits(model, "merMod")) {
    lme4::fixef(model)
  } else {
    coef(model)
  }
  if (is.null(covariance)) covariance <- vcov(model)
  covariance <- as.matrix(covariance)
  coefficient_name <- names(coefficients)[
    grepl("ArmAuto-FMT", names(coefficients), fixed = TRUE) &
      grepl(pattern, names(coefficients), fixed = TRUE)
  ]
  if (length(coefficient_name) != 1) {
    stop("Cannot identify interaction term for ", pattern)
  }
  estimate <- coefficients[coefficient_name]
  se <- sqrt(covariance[coefficient_name, coefficient_name])
  tibble(
    Outcome = outcome,
    Term = pattern,
    Estimate = unname(estimate),
    SE = unname(se),
    Lower = unname(estimate - 1.96 * se),
    Upper = unname(estimate + 1.96 * se),
    P = unname(2 * pnorm(abs(estimate / se), lower.tail = FALSE))
  )
}

event_pretrend_tests <- bind_rows(
  extract_named_interaction(
    alpha_event_model,
    "PreTime",
    "HRIC alpha diversity"
  ),
  extract_named_interaction(
    distance_event_model,
    "PreTime",
    "HRIC distance to personal baseline",
    covariance = distance_event_covariance
  )
)

# A clinically intuitive threshold is retained as a descriptive secondary
# endpoint because stool sampling is irregular. A sustained recovery requires
# the first crossing to be confirmed at the next available sampling day.
patient_day_alpha <- meta %>%
  filter(EventTime >= 1, EventTime <= 60) %>%
  group_by(Patient, Arm, EventTime) %>%
  summarise(Alpha = mean(Alpha), .groups = "drop")

threshold_recovery <- bind_rows(lapply(
  seq_len(nrow(patient_covariates)),
  function(i) {
    patient_row <- patient_covariates[i, ]
    patient_days <- patient_day_alpha %>%
      filter(Patient == patient_row$Patient) %>%
      arrange(EventTime)
    threshold <- 0.8 * patient_row$BaselineAlpha
    crossings <- which(patient_days$Alpha >= threshold)
    first_crossing <- if (length(crossings) > 0) {
      patient_days$EventTime[crossings[1]]
    } else {
      NA_real_
    }
    sustained_crossing <- NA_real_
    if (length(crossings) > 0) {
      for (crossing in crossings) {
        later <- which(
          patient_days$EventTime >
            patient_days$EventTime[crossing]
        )
        if (
          length(later) > 0 &&
            patient_days$Alpha[later[1]] >= threshold
        ) {
          sustained_crossing <- patient_days$EventTime[crossing]
          break
        }
      }
    }
    tibble(
      Patient = patient_row$Patient,
      Arm = patient_row$Arm,
      BaselineAlpha = patient_row$BaselineAlpha,
      PreIndexAlpha = patient_row$PreIndexAlpha,
      RecoveryThreshold80 = threshold,
      AtRiskAtIndex = patient_row$PreIndexAlpha < threshold,
      FirstCrossingDay = first_crossing,
      SustainedCrossingDay = sustained_crossing,
      FirstCrossingBy30 = !is.na(first_crossing) & first_crossing <= 30,
      SustainedCrossingBy30 =
        !is.na(sustained_crossing) & sustained_crossing <= 30
    )
  }
))

threshold_recovery_summary <- threshold_recovery %>%
  filter(AtRiskAtIndex) %>%
  group_by(Arm) %>%
  summarise(
    AtRisk = n(),
    FirstCrossingBy30 = sum(FirstCrossingBy30),
    SustainedCrossingBy30 = sum(SustainedCrossingBy30),
    .groups = "drop"
  )

threshold_fisher_test <- function(endpoint) {
  summary_data <- threshold_recovery_summary
  cases <- summary_data[[endpoint]]
  totals <- summary_data$AtRisk
  test <- fisher.test(cbind(cases, totals - cases))
  tibble(
    Endpoint = endpoint,
    ControlCases = cases[summary_data$Arm == "Control"],
    ControlTotal = totals[summary_data$Arm == "Control"],
    AutoFMTCases = cases[summary_data$Arm == "Auto-FMT"],
    AutoFMTTotal = totals[summary_data$Arm == "Auto-FMT"],
    FisherP = test$p.value
  )
}

threshold_recovery_tests <- bind_rows(
  threshold_fisher_test("FirstCrossingBy30"),
  threshold_fisher_test("SustainedCrossingBy30")
)

############################################################
## Taxa contributing to arm-level alpha recovery
############################################################

compute_taxon_recovery <- function(post_landmark) {
  paired <- landmarks %>%
    filter(Landmark %in% c("Pre-index", post_landmark)) %>%
    select(Patient, Arm, Landmark, SampleID) %>%
    pivot_wider(names_from = Landmark, values_from = SampleID) %>%
    filter(!is.na(`Pre-index`), !is.na(.data[[post_landmark]]))

  recovery_matrix <- t(vapply(seq_len(nrow(paired)), function(i) {
    taxon_alpha_deficit[paired$`Pre-index`[i], ] -
      taxon_alpha_deficit[paired[[post_landmark]][i], ]
  }, numeric(ncol(taxon_alpha_deficit))))

  as_tibble(
    recovery_matrix,
    .name_repair = "minimal"
  ) %>%
    mutate(
      Patient = paired$Patient,
      Arm = paired$Arm,
      Landmark = post_landmark,
      .before = 1
    ) %>%
    pivot_longer(
      cols = -c(Patient, Arm, Landmark),
      names_to = "Taxon",
      values_to = "RecoveryContribution"
    )
}

taxon_recovery <- bind_rows(
  compute_taxon_recovery("Early post-index"),
  compute_taxon_recovery("Later post-index")
)

taxon_recovery_summary <- taxon_recovery %>%
  group_by(Landmark, Arm, Taxon) %>%
  summarise(
    MeanRecoveryContribution = mean(RecoveryContribution),
    SE = sd(RecoveryContribution) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

top_recovery_taxa <- taxon_recovery_summary %>%
  group_by(Taxon) %>%
  summarise(
    DriverScore = max(abs(MeanRecoveryContribution)),
    .groups = "drop"
  ) %>%
  arrange(desc(DriverScore), Taxon) %>%
  slice_head(n = 10)

taxon_profile_similarity <- bind_rows(lapply(
  c("Early post-index", "Later post-index"),
  function(landmark_name) {
    wide <- taxon_recovery_summary %>%
      filter(Landmark == landmark_name) %>%
      select(Taxon, Arm, MeanRecoveryContribution) %>%
      pivot_wider(names_from = Arm, values_from = MeanRecoveryContribution)
    top_control <- wide %>%
      arrange(desc(abs(Control))) %>%
      slice_head(n = 10) %>%
      pull(Taxon)
    top_treatment <- wide %>%
      arrange(desc(abs(`Auto-FMT`))) %>%
      slice_head(n = 10) %>%
      pull(Taxon)
    comparison_taxa <- union(
      wide %>%
        arrange(desc(abs(Control))) %>%
        slice_head(n = 50) %>%
        pull(Taxon),
      wide %>%
        arrange(desc(abs(`Auto-FMT`))) %>%
        slice_head(n = 50) %>%
        pull(Taxon)
    )
    comparison <- wide %>% filter(Taxon %in% comparison_taxa)
    tibble(
      Landmark = landmark_name,
      SpearmanRho = cor(
        comparison$Control,
        comparison$`Auto-FMT`,
        method = "spearman"
      ),
      Top10Overlap = length(intersect(top_control, top_treatment)),
      Top10Union = length(union(top_control, top_treatment))
    )
  }
))

############################################################
## Figure 2 panels
############################################################

p2_alpha <- ggplot(
  landmarks,
  aes(x = Landmark, y = Alpha, color = Arm)
) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.11, dodge.width = 0.55, seed = 17
    ),
    alpha = 0.58, size = 1.7
  ) +
  geom_errorbar(
    data = landmark_summary,
    aes(y = y, ymin = ymin, ymax = ymax, group = Arm),
    position = position_dodge(width = 0.55),
    width = 0.08, linewidth = 0.70
  ) +
  geom_point(
    data = landmark_summary,
    aes(y = y, group = Arm),
    position = position_dodge(width = 0.55),
    shape = 21, fill = "white", stroke = 0.85, size = 3.0
  ) +
  scale_color_manual(values = arm_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
  labs(
    title = "Patient-level alpha diversity",
    x = NULL,
    y = "HRIC alpha diversity",
    color = NULL
  ) +
  theme_nature(10) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.margin = margin(7, 9, 18, 9)
  )

p2_decomp <- ggplot(
  group_decomp_long,
  aes(x = Arm, y = Value, fill = Component)
) +
  geom_col(width = 0.64, color = "white", linewidth = 0.30) +
  geom_point(
    data = group_decomp,
    aes(x = Arm, y = Gamma),
    inherit.aes = FALSE, shape = 23, size = 2.1,
    fill = "white", color = "grey10", stroke = 0.55
  ) +
  geom_text(
    data = group_decomp,
    aes(
      x = Arm, y = Gamma,
      label = sprintf("atop(gamma == %.3f, italic(n) == %d)", Gamma, n)
    ),
    inherit.aes = FALSE, vjust = -0.65, size = 2.65,
    parse = TRUE
  ) +
  geom_text(
    data = delta_by_landmark,
    aes(
      x = 1.5,
      y = max(group_decomp$Gamma) * 1.22,
      label = sprintf("delta == %.4f", Delta)
    ),
    inherit.aes = FALSE, size = 2.8, color = "grey20",
    parse = TRUE
  ) +
  facet_wrap(~ Landmark, nrow = 1) +
  scale_fill_manual(
    values = component_colors,
    breaks = c("Mean alpha", "Beta"),
    labels = c("Mean alpha", "Regional beta")
  ) +
  scale_y_continuous(
    limits = c(0, max(group_decomp$Gamma) * 1.30),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Arm-level local-regional diversity decomposition",
    x = NULL,
    y = "HRIC diversity",
    fill = NULL
  ) +
  theme_nature(9.5) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.spacing.x = grid::unit(4, "mm"),
    plot.margin = margin(7, 9, 18, 9)
  )

p2_alpha_effect <- ggplot(alpha_effect, aes(y = Metric, x = Estimate)) +
  geom_vline(xintercept = 0, color = "grey45", linewidth = 0.55) +
  geom_errorbar(
    aes(xmin = Lower, xmax = Upper),
    orientation = "y", width = 0,
    linewidth = 0.85, color = arm_colors[["Auto-FMT"]]
  ) +
  geom_point(
    shape = 21, size = 3.4, stroke = 0.85, fill = "white",
    color = arm_colors[["Auto-FMT"]]
  ) +
  annotate(
    "text",
    x = Inf, y = 1,
    label = paste0(
      "estimate = ", sprintf("%.3f", alpha_effect$Estimate),
      "\n95% CI ", sprintf("%.3f", alpha_effect$Lower),
      " to ", sprintf("%.3f", alpha_effect$Upper),
      "\nP ", ifelse(alpha_effect$P < 0.001, "< 0.001", paste0("= ", format_p(alpha_effect$P)))
    ),
    hjust = 1, vjust = -0.15, size = 3.0, lineheight = 1.05
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.38))) +
  labs(
    title = "Time-adjusted alpha effect",
    subtitle = "Mixed model controls allo-HSCT day and repeated samples",
    x = "Additional post-index effect of auto-FMT",
    y = NULL
  ) +
  theme_nature(10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  )

null_df <- tibble(Statistic = permuted_change_stats)
p2_change_test <- ggplot(null_df, aes(x = Statistic)) +
  geom_histogram(
    bins = 34, fill = "grey78", color = "white", linewidth = 0.25
  ) +
  geom_vline(
    xintercept = observed_change_stat,
    color = arm_colors[["Auto-FMT"]], linewidth = 1.0
  ) +
  annotate(
    "text",
    x = observed_change_stat, y = Inf,
    label = paste0(
      "observed\nP = ", format_p(change_p)
    ),
    hjust = -0.08, vjust = 1.15, size = 3.0,
    color = arm_colors[["Auto-FMT"]], fontface = "bold"
  ) +
  labs(
    title = "Time-adjusted HRIC coordinate change",
    subtitle = paste0(
      "Patient-level permutation null (", n_permutations, " permutations)"
    ),
    x = "Between-arm change statistic",
    y = "Permutation count"
  ) +
  theme_nature(10)

alpha_day14 <- recovery_arm_contrasts %>%
  filter(Outcome == "HRIC alpha diversity", EventTime == 14)
alpha_day60 <- recovery_arm_contrasts %>%
  filter(Outcome == "HRIC alpha diversity", EventTime == 60)
distance_day14 <- recovery_arm_contrasts %>%
  filter(Outcome == "HRIC distance to baseline", EventTime == 14)
distance_day60 <- recovery_arm_contrasts %>%
  filter(Outcome == "HRIC distance to baseline", EventTime == 60)

alpha_early_average <- recovery_average_effects %>%
  filter(Outcome == "HRIC alpha diversity", WindowEnd == 30)
distance_early_average <- recovery_average_effects %>%
  filter(
    Outcome == "HRIC distance to personal baseline",
    WindowEnd == 30
  )
threshold_sustained <- threshold_recovery_tests %>%
  filter(Endpoint == "SustainedCrossingBy30")

format_signed <- function(x, digits = 3) {
  sprintf(paste0("%+.", digits, "f"), x)
}

p2_alpha_recovery <- ggplot(
  alpha_recovery_curve,
  aes(x = EventTime, y = Estimate, color = Arm, fill = Arm)
) +
  annotate(
    "rect",
    xmin = 1, xmax = 30, ymin = -Inf, ymax = Inf,
    fill = "grey35", alpha = 0.055
  ) +
  geom_hline(
    yintercept = mean(patient_covariates$BaselineAlpha),
    linewidth = 0.55, linetype = "22", color = "grey35"
  ) +
  geom_ribbon(
    aes(ymin = Lower, ymax = Upper),
    color = NA, alpha = 0.16
  ) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = arm_colors) +
  scale_fill_manual(values = arm_colors) +
  scale_x_continuous(breaks = c(1, 14, 30, 45, 60)) +
  labs(
    title = "Alpha-diversity recovery",
    subtitle = paste0(
      "Days 1-30 adjusted mean: ",
      format_signed(alpha_early_average$Estimate),
      " (95% bootstrap CI ",
      format_signed(alpha_early_average$BootstrapLower),
      " to ", format_signed(alpha_early_average$BootstrapUpper), ")\n",
      "Day 60 contrast: ", format_signed(alpha_day60$Estimate),
      " (95% model CI ", format_signed(alpha_day60$Lower),
      " to ", format_signed(alpha_day60$Upper), ")"
    ),
    x = "Days after randomization/FMT index",
    y = "Adjusted HRIC alpha diversity",
    color = NULL,
    fill = NULL
  ) +
  theme_nature(10) +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(lineheight = 1.05)
  )

p2_distance_recovery <- ggplot(
  distance_recovery_curve,
  aes(x = EventTime, y = Estimate, color = Arm, fill = Arm)
) +
  annotate(
    "rect",
    xmin = 1, xmax = 30, ymin = -Inf, ymax = Inf,
    fill = "grey35", alpha = 0.055
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.55, linetype = "22", color = "grey35"
  ) +
  geom_ribbon(
    aes(ymin = Lower, ymax = Upper),
    color = NA, alpha = 0.16
  ) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = arm_colors) +
  scale_fill_manual(values = arm_colors) +
  scale_x_continuous(breaks = c(1, 14, 30, 45, 60)) +
  labs(
    title = "Restoration of personal composition",
    subtitle = paste0(
      "Days 1-30 adjusted mean: ",
      format_signed(distance_early_average$Estimate),
      " (95% bootstrap CI ",
      format_signed(distance_early_average$BootstrapLower),
      " to ", format_signed(distance_early_average$BootstrapUpper), ")\n",
      "Day 60 contrast: ", format_signed(distance_day60$Estimate),
      " (95% model CI ", format_signed(distance_day60$Lower),
      " to ", format_signed(distance_day60$Upper), "); lower is closer"
    ),
    x = "Days after randomization/FMT index",
    y = "Adjusted HRIC distance to baseline",
    color = NULL,
    fill = NULL
  ) +
  theme_nature(10) +
  theme(
    legend.position = "top",
    plot.subtitle = element_text(lineheight = 1.05)
  )

taxon_plot_data <- taxon_recovery_summary %>%
  filter(Taxon %in% top_recovery_taxa$Taxon) %>%
  mutate(
    Taxon = factor(Taxon, levels = rev(top_recovery_taxa$Taxon)),
    MeanContributionScaled = 1000 * MeanRecoveryContribution,
    Landmark = factor(
      Landmark,
      levels = c("Early post-index", "Later post-index")
    )
  )

taxon_plot_segments <- taxon_plot_data %>%
  select(Landmark, Taxon, Arm, MeanContributionScaled) %>%
  pivot_wider(names_from = Arm, values_from = MeanContributionScaled)

similarity_label <- paste(vapply(
  seq_len(nrow(taxon_profile_similarity)),
  function(i) {
    row <- taxon_profile_similarity[i, ]
    paste0(
      ifelse(row$Landmark == "Early post-index", "Early", "Later"),
      ": rho = ", sprintf("%.2f", row$SpearmanRho),
      ", top-10 overlap ", row$Top10Overlap, "/", row$Top10Union
    )
  },
  character(1)
), collapse = "; ")

p2_taxa <- ggplot(
  taxon_plot_data,
  aes(x = MeanContributionScaled, y = Taxon, color = Arm)
) +
  geom_vline(xintercept = 0, linewidth = 0.45, color = "grey55") +
  geom_segment(
    data = taxon_plot_segments,
    aes(
      x = Control, xend = `Auto-FMT`,
      y = Taxon, yend = Taxon
    ),
    inherit.aes = FALSE, linewidth = 0.65, color = "grey72"
  ) +
  geom_point(size = 2.5, alpha = 0.92) +
  facet_wrap(~ Landmark, nrow = 1) +
  scale_color_manual(values = arm_colors) +
  labs(
    title = "Taxa contributing to alpha recovery",
    subtitle = paste0(
      "Positive values reduce contributions to 1 - alpha; negative values increase them.\n",
      "Profile concordance: ", similarity_label
    ),
    x = expression("Mean contribution to alpha recovery (" * 10^-3 * ")"),
    y = NULL,
    color = NULL
  ) +
  theme_nature(9.2) +
  theme(
    legend.position = "none",
    panel.spacing.x = grid::unit(5, "mm"),
    panel.grid.major.y = element_line(linewidth = 0.25, color = "grey90")
  )

figure2 <- (p2_alpha_recovery | p2_distance_recovery) /
  (p2_decomp | p2_taxa) +
  plot_layout(
    heights = c(0.88, 1.25),
    widths = c(1.12, 1),
    guides = "keep"
  ) +
  plot_annotation(
    title = "Early gut microbiota recovery after auto-FMT in allo-HSCT recipients",
    subtitle = paste0(
      "Primary window: days 1-30; models adjust personal pre-HSCT alpha, ",
      "pre-index damage and randomization day, with patient-level resampling.\n",
      "Descriptive sustained recovery to at least 80% of personal alpha baseline by day 30: ",
      threshold_sustained$AutoFMTCases, "/", threshold_sustained$AutoFMTTotal,
      " auto-FMT versus ", threshold_sustained$ControlCases, "/",
      threshold_sustained$ControlTotal, " control (Fisher's exact P = ",
      format_p(threshold_sustained$FisherP),
      "); index-day stools were excluded."
    ),
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = element_text(
        size = 9.8, color = "grey25", hjust = 0, lineheight = 1.05
      ),
      plot.tag = element_text(face = "bold", size = 14)
    )
  )

figure2_files <- save_plot(
  figure2, "figure2_fmt_hric_randomized_comparison",
  width = 13.8, height = 10.8
)

write.csv(
  landmarks %>%
    select(
      SampleID, sample_name, Patient, Arm, day, randomizationday,
      EventTime, Landmark, Alpha
    ),
  file.path(tab_dir, "figure2_selected_landmark_samples.csv"),
  row.names = FALSE
)
write.csv(
  group_decomp %>% left_join(delta_by_landmark, by = "Landmark"),
  file.path(tab_dir, "figure2_group_alpha_beta_gamma_delta.csv"),
  row.names = FALSE
)
write.csv(
  alpha_effect,
  file.path(tab_dir, "figure2_time_adjusted_alpha_model.csv"),
  row.names = FALSE
)
write.csv(
  change_test,
  file.path(tab_dir, "figure2_time_adjusted_hric_change_test.csv"),
  row.names = FALSE
)
write.csv(
  null_df,
  file.path(tab_dir, "figure2_time_adjusted_hric_change_null.csv"),
  row.names = FALSE
)
write.csv(
  patient_covariates,
  file.path(tab_dir, "figure2_patient_baseline_covariates.csv"),
  row.names = FALSE
)
write.csv(
  recovery_speed_effects,
  file.path(tab_dir, "figure2_baseline_adjusted_recovery_speed.csv"),
  row.names = FALSE
)
write.csv(
  recovery_arm_contrasts,
  file.path(tab_dir, "figure2_baseline_adjusted_arm_contrasts.csv"),
  row.names = FALSE
)
write.csv(
  recovery_sensitivity_contrasts,
  file.path(tab_dir, "figure2_recovery_sensitivity_contrasts.csv"),
  row.names = FALSE
)
write.csv(
  recovery_average_effects,
  file.path(tab_dir, "figure2_time_averaged_recovery_effects.csv"),
  row.names = FALSE
)
write.csv(
  recovery_bootstrap_draws,
  file.path(tab_dir, "figure2_patient_bootstrap_recovery_draws.csv"),
  row.names = FALSE
)
write.csv(
  event_pretrend_tests,
  file.path(tab_dir, "figure2_preindex_trend_tests.csv"),
  row.names = FALSE
)
write.csv(
  event_did_average_effects,
  file.path(tab_dir, "figure2_event_time_did_sensitivity.csv"),
  row.names = FALSE
)
write.csv(
  threshold_recovery,
  file.path(tab_dir, "figure2_patient_recovery_thresholds.csv"),
  row.names = FALSE
)
write.csv(
  threshold_recovery_summary %>%
    tidyr::crossing(
      Endpoint = threshold_recovery_tests$Endpoint
    ) %>%
    left_join(
      threshold_recovery_tests %>%
        select(Endpoint, FisherP),
      by = "Endpoint"
    ),
  file.path(tab_dir, "figure2_recovery_threshold_summary.csv"),
  row.names = FALSE
)
write.csv(
  taxon_recovery,
  file.path(tab_dir, "figure2_patient_taxon_recovery_contributions.csv"),
  row.names = FALSE
)
write.csv(
  taxon_recovery_summary,
  file.path(tab_dir, "figure2_arm_taxon_recovery_summary.csv"),
  row.names = FALSE
)
write.csv(
  taxon_profile_similarity,
  file.path(tab_dir, "figure2_taxon_profile_similarity.csv"),
  row.names = FALSE
)

capture.output(
  summary(alpha_model),
  file = file.path(out_dir, "time_adjusted_alpha_model_summary.txt")
)
capture.output(
  summary(alpha_recovery_model),
  file = file.path(out_dir, "baseline_adjusted_alpha_recovery_model.txt")
)
capture.output(
  summary(distance_recovery_model),
  file = file.path(out_dir, "baseline_adjusted_distance_recovery_model.txt")
)
capture.output(
  summary(alpha_recovery_model_day0),
  file = file.path(out_dir, "sensitivity_alpha_recovery_model_including_day0.txt")
)
capture.output(
  summary(distance_recovery_model_day0),
  file = file.path(out_dir, "sensitivity_distance_recovery_model_including_day0.txt")
)
capture.output(
  summary(alpha_event_model),
  file = file.path(out_dir, "sensitivity_alpha_event_time_model.txt")
)
capture.output(
  summary(distance_event_model),
  file = file.path(out_dir, "sensitivity_distance_event_time_model.txt")
)

run_summary <- c(
  "HRIC auto-FMT analysis completed.",
  paste("HRIC package version:", as.character(packageVersion("HRIC"))),
  paste("Samples:", nrow(meta)),
  paste("Patients:", nlevels(meta$Patient)),
  paste("Features in fixed dictionary:", ncol(counts)),
  paste("Control patients:", n_distinct(meta$Patient[meta$Arm == "Control"])),
  paste("Auto-FMT patients:", n_distinct(meta$Patient[meta$Arm == "Auto-FMT"])),
  paste("Primary recovery samples (days 1-60):", nrow(post_recovery)),
  paste(
    "Patient bootstrap iterations:",
    min(recovery_bootstrap_summary$BootstrapIterations)
  ),
  paste("Index-day stools excluded from primary recovery models:", sum(meta$EventTime == 0)),
  paste("Paired patients in adjusted HRIC change test:", nrow(paired_landmarks)),
  paste(
    "Figure 1 all-patient alpha atlas:",
    paste(patient_alpha_atlas_files, collapse = ", ")
  ),
  paste(
    "Figure 1 all-patient inverse Simpson atlas:",
    paste(patient_inverse_simpson_atlas_files, collapse = ", ")
  ),
  paste(
    "Figure 1 calculated log inverse Simpson atlas:",
    paste(patient_formula_log_is_atlas_files, collapse = ", ")
  ),
  paste(
    "Figure 1 calculated Pielou evenness atlas:",
    paste(patient_pielou_atlas_files, collapse = ", ")
  ),
  paste("Figure 1:", paste(figure1_files, collapse = ", ")),
  paste("Figure 2:", paste(figure2_files, collapse = ", "))
)
writeLines(run_summary, file.path(out_dir, "run_summary.txt"))
message(paste(run_summary, collapse = "\n"))
