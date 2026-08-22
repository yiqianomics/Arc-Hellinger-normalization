# Longitudinal gut microbiome analysis of the randomized auto-FMT trial.

options(stringsAsFactors = FALSE, width = 140)
set.seed(20260723)

required_packages <- c(
  "HRIC", "phyloseq", "ggplot2", "dplyr", "tidyr", "tibble",
  "patchwork", "scales", "lme4", "sandwich", "reformulas"
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
  library(lme4)
})

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

script_path <- get_script_path()
analysis_dir <- if (file.exists(script_path)) dirname(script_path) else script_path
out_dir <- file.path(analysis_dir, "results")
fig_dir <- out_dir
tab_dir <- file.path(out_dir, "tables")
unlink(out_dir, recursive = TRUE, force = TRUE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, filename, width, height) {
  path <- file.path(fig_dir, paste0(filename, ".pdf"))
  ggsave(path, plot, width = width, height = height, bg = "white")
  invisible(path)
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

# Data preparation

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
  "host_subject_id", "randomizationday"
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
    Patient = factor(as.character(host_subject_id)),
    EventTime = day - randomizationday
  )

counts <- as(otu_table(ps_gut), "matrix")
if (taxa_are_rows(ps_gut)) counts <- t(counts)
storage.mode(counts) <- "numeric"
counts <- counts[meta$SampleID, , drop = FALSE]
counts <- counts[, colSums(counts) > 0, drop = FALSE]

# Every HRIC calculation uses this same feature dictionary. This matters because
# the package's alpha scale depends on the number of simplex components.
stopifnot(identical(rownames(counts), meta$SampleID))
alpha_all <- HRIC::SHalpha(counts)
meta$Alpha <- as.numeric(alpha_all[meta$SampleID])

# Squared HRIC coordinates decompose 1 - alpha into taxon contributions:
# sum_j Z_ij^2 / A_p^2 = 1 - alpha_i.
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

# Figure 1: longitudinal HRIC for patients 0009 and 0028

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

# Figure 2: randomized-arm comparison around the index day

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

# Baseline-adjusted recovery

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
        Patient, BaselineAlpha, PreIndexAlpha, PreIndexDistance
      ),
    by = "Patient"
  ) %>%
  mutate(
    BaselineAlphaC = BaselineAlpha - mean(BaselineAlpha, na.rm = TRUE),
    PreIndexAlphaC = PreIndexAlpha - mean(PreIndexAlpha, na.rm = TRUE),
    PreIndexDistanceC = PreIndexDistance -
      mean(PreIndexDistance, na.rm = TRUE),
    IndexDayC = randomizationday - mean(randomizationday, na.rm = TRUE),
    EarlyTime = pmin(pmax(EventTime - 1, 0), 13),
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

recovery_grid <- expand_grid(
  EventTime = seq(1, 60, by = 1),
  Arm = factor(c("Control", "Auto-FMT"), levels = levels(meta$Arm))
) %>%
  mutate(
    BaselineAlphaC = 0,
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

make_arm_contrast <- function(
  model,
  day_value,
  outcome,
  covariance = NULL
) {
  common <- tibble(
    EventTime = day_value,
    BaselineAlphaC = 0,
    PreIndexAlphaC = 0,
    PreIndexDistanceC = 0,
    IndexDayC = 0,
    EarlyTime = pmin(pmax(day_value - 1, 0), 13),
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

bootstrap_iterations <- 4999L
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

# Taxa contributing to arm-level alpha recovery

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

# Figure 2 panels

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

alpha_day60 <- recovery_arm_contrasts %>%
  filter(Outcome == "HRIC alpha diversity", EventTime == 60)
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
  patient_covariates,
  file.path(tab_dir, "figure2_patient_baseline_covariates.csv"),
  row.names = FALSE
)
write.csv(
  bind_rows(
    alpha_recovery_curve %>% mutate(Outcome = "HRIC alpha diversity"),
    distance_recovery_curve %>%
      mutate(Outcome = "HRIC distance to personal baseline")
  ) %>%
    select(Outcome, EventTime, Arm, Estimate, Lower, Upper),
  file.path(tab_dir, "figure2_adjusted_recovery_curves.csv"),
  row.names = FALSE
)
write.csv(
  recovery_arm_contrasts,
  file.path(tab_dir, "figure2_baseline_adjusted_arm_contrasts.csv"),
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
  taxon_recovery_summary,
  file.path(tab_dir, "figure2_arm_taxon_recovery_summary.csv"),
  row.names = FALSE
)
write.csv(
  taxon_profile_similarity,
  file.path(tab_dir, "figure2_taxon_profile_similarity.csv"),
  row.names = FALSE
)

bootstrap_check <- recovery_bootstrap_draws %>%
  group_by(Outcome, WindowEnd) %>%
  summarise(
    n = n(),
    Lower = quantile(Estimate, 0.025, names = FALSE),
    Upper = quantile(Estimate, 0.975, names = FALSE),
    .groups = "drop"
  ) %>%
  left_join(
    recovery_bootstrap_summary %>%
      select(Outcome, WindowEnd, BootstrapIterations, BootstrapLower, BootstrapUpper),
    by = c("Outcome", "WindowEnd")
  )

stopifnot(
  all(c("HRIC", "SHalpha", "SHbeta", "SHgamma", "SHdelta") %in%
    getNamespaceExports("HRIC")),
  max(abs(rowSums(feature_alpha_deficit) - (1 - meta$Alpha))) < 1e-10,
  all(abs(group_decomp$Gamma - group_decomp$MeanAlpha - group_decomp$Beta) < 1e-10),
  all(bootstrap_check$n == bootstrap_iterations),
  all(abs(bootstrap_check$Lower - bootstrap_check$BootstrapLower) < 1e-12),
  all(abs(bootstrap_check$Upper - bootstrap_check$BootstrapUpper) < 1e-12),
  file.exists(file.path(fig_dir, "figure1_fmt_hric_patients_0009_0028.pdf")),
  file.exists(file.path(fig_dir, "figure2_fmt_hric_randomized_comparison.pdf"))
)

message("Analysis complete. Results written to: ", out_dir)
