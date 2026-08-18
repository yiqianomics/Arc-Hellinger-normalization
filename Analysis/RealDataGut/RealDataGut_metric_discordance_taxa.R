suppressPackageStartupMessages({
  library(phyloseq)
  library(HRIC)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

analysis_dir <- normalizePath(
  "RealDataGut",
  winslash = "/",
  mustWork = TRUE
)
result_dir <- file.path(
  analysis_dir,
  "hric_fmt_diversity_results"
)
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(analysis_dir, "gut.RData"))

counts <- as(otu_table(gut), "matrix")
if (taxa_are_rows(gut)) counts <- t(counts)
storage.mode(counts) <- "numeric"
counts <- counts[, colSums(counts) > 0, drop = FALSE]
relative_counts <- counts / rowSums(counts)

alpha_data <- read.csv(
  file.path(
    table_dir,
    "figure1_all_patient_alpha_trajectories.csv"
  ),
  stringsAsFactors = FALSE
)
inverse_simpson_data <- read.csv(
  file.path(
    table_dir,
    "figure1_all_patient_log_inverse_simpson_calculated_trajectories.csv"
  ),
  stringsAsFactors = FALSE
)
evenness_data <- read.csv(
  file.path(
    table_dir,
    "figure1_all_patient_pielou_evenness_trajectories.csv"
  ),
  stringsAsFactors = FALSE
)

metric_data <- alpha_data %>%
  select(
    SampleID, ParticipantID, day, EventTime, Alpha
  ) %>%
  left_join(
    inverse_simpson_data %>%
      select(SampleID, LogInverseSimpsonCalculated),
    by = "SampleID"
  ) %>%
  left_join(
    evenness_data %>%
      select(SampleID, PielouEvenness),
    by = "SampleID"
  ) %>%
  mutate(
    ParticipantID = sprintf("%04d", as.integer(ParticipantID))
  ) %>%
  filter(ParticipantID %in% c("0020", "0061")) %>%
  arrange(ParticipantID, day, SampleID)

transition_data <- metric_data %>%
  group_by(ParticipantID) %>%
  mutate(
    PreviousSampleID = lag(SampleID),
    PreviousDay = lag(day),
    DeltaHRIC = Alpha - lag(Alpha),
    DeltaLogInverseSimpson =
      LogInverseSimpsonCalculated -
      lag(LogInverseSimpsonCalculated),
    DeltaPielou = PielouEvenness - lag(PielouEvenness)
  ) %>%
  ungroup() %>%
  filter(
    !is.na(PreviousSampleID),
    DeltaHRIC * DeltaLogInverseSimpson < 0 |
      DeltaHRIC * DeltaPielou < 0
  ) %>%
  mutate(
    Transition = paste0(
      ParticipantID, ": day ", PreviousDay, " to ", day
    ),
    TransitionOrder = row_number()
  )

taxonomy <- as.data.frame(
  as(tax_table(gut), "matrix"),
  stringsAsFactors = FALSE
)
taxonomy$FeatureID <- rownames(taxonomy)

clean_rank <- function(x) {
  x <- sub("^[a-z]__", "", x)
  x[is.na(x) | x == ""] <- NA_character_
  x
}

taxonomy <- taxonomy %>%
  mutate(
    GenusClean = clean_rank(Genus),
    FamilyClean = clean_rank(Family),
    OrderClean = clean_rank(Order),
    PhylumClean = clean_rank(Phylum),
    Taxon = case_when(
      !is.na(GenusClean) ~ paste0(GenusClean, " (genus)"),
      !is.na(FamilyClean) ~ paste0(FamilyClean, " (family)"),
      !is.na(OrderClean) ~ paste0(OrderClean, " (order)"),
      !is.na(PhylumClean) ~ paste0(PhylumClean, " (phylum)"),
      TRUE ~ "Unclassified"
    )
  )

taxonomy <- taxonomy[match(colnames(counts), taxonomy$FeatureID), ]
stopifnot(identical(taxonomy$FeatureID, colnames(counts)))

hric_coordinates <- HRIC::HRIC(counts)
hric_alpha <- HRIC::SHalpha(counts)
hric_radius <- asin(sqrt(1 - 1 / ncol(counts)))
hric_deficit <- hric_coordinates^2 / hric_radius^2

alpha_check <- max(abs(
  hric_alpha[metric_data$SampleID] - metric_data$Alpha
))
if (!is.finite(alpha_check) || alpha_check > 1e-12) {
  stop("Saved HRIC alpha values do not match HRIC::SHalpha.")
}

entropy_components <- ifelse(
  relative_counts > 0,
  -relative_counts * log(relative_counts),
  0
)

decompose_transition <- function(i) {
  transition <- transition_data[i, ]
  before_id <- transition$PreviousSampleID
  after_id <- transition$SampleID

  before_p <- relative_counts[before_id, ]
  after_p <- relative_counts[after_id, ]
  before_h <- entropy_components[before_id, ]
  after_h <- entropy_components[after_id, ]
  before_richness <- sum(before_p > 0)
  after_richness <- sum(after_p > 0)
  before_entropy <- sum(before_h)
  after_entropy <- sum(after_h)
  before_inverse_log_richness <- 1 / log(before_richness)
  after_inverse_log_richness <- 1 / log(after_richness)

  pielou_taxon_contribution <-
    0.5 *
    (
      before_inverse_log_richness +
        after_inverse_log_richness
    ) *
    (after_h - before_h)
  pielou_richness_contribution <-
    0.5 *
    (before_entropy + after_entropy) *
    (
      after_inverse_log_richness -
        before_inverse_log_richness
    )

  feature_result <- tibble(
    Transition = transition$Transition,
    TransitionOrder = transition$TransitionOrder,
    ParticipantID = transition$ParticipantID,
    BeforeSampleID = before_id,
    AfterSampleID = after_id,
    BeforeDay = transition$PreviousDay,
    AfterDay = transition$day,
    FeatureID = colnames(counts),
    Taxon = taxonomy$Taxon,
    BeforeRelativeAbundance = as.numeric(before_p),
    AfterRelativeAbundance = as.numeric(after_p),
    HRICContribution = as.numeric(
      hric_deficit[before_id, ] -
        hric_deficit[after_id, ]
    ),
    SimpsonDirectionContribution = as.numeric(
      before_p^2 - after_p^2
    ),
    PielouEntropyContribution = as.numeric(
      pielou_taxon_contribution
    )
  )

  check_values <- tibble(
    Transition = transition$Transition,
    TransitionOrder = transition$TransitionOrder,
    ParticipantID = transition$ParticipantID,
    BeforeSampleID = before_id,
    AfterSampleID = after_id,
    BeforeDay = transition$PreviousDay,
    AfterDay = transition$day,
    DeltaHRIC = transition$DeltaHRIC,
    ReconstructedDeltaHRIC =
      sum(feature_result$HRICContribution),
    DeltaLogInverseSimpson =
      transition$DeltaLogInverseSimpson,
    DeltaSimpsonConcentration =
      sum(feature_result$SimpsonDirectionContribution),
    DeltaPielou = transition$DeltaPielou,
    PielouTaxonEntropyComponent =
      sum(feature_result$PielouEntropyContribution),
    PielouRichnessNormalizationComponent =
      pielou_richness_contribution,
    ReconstructedDeltaPielou =
      sum(feature_result$PielouEntropyContribution) +
      pielou_richness_contribution,
    BeforeRichness = before_richness,
    AfterRichness = after_richness
  )

  list(features = feature_result, checks = check_values)
}

decompositions <- lapply(
  seq_len(nrow(transition_data)),
  decompose_transition
)
feature_contributions <- bind_rows(lapply(
  decompositions,
  `[[`,
  "features"
))
decomposition_checks <- bind_rows(lapply(
  decompositions,
  `[[`,
  "checks"
))

if (
  max(abs(
    decomposition_checks$DeltaHRIC -
      decomposition_checks$ReconstructedDeltaHRIC
  )) > 1e-12 ||
    max(abs(
      decomposition_checks$DeltaPielou -
        decomposition_checks$ReconstructedDeltaPielou
    )) > 1e-12
) {
  stop("Taxon decomposition did not reconstruct metric changes.")
}

taxon_contributions <- feature_contributions %>%
  group_by(
    Transition, TransitionOrder, ParticipantID,
    BeforeSampleID, AfterSampleID, BeforeDay, AfterDay,
    Taxon
  ) %>%
  summarise(
    BeforeRelativeAbundance = sum(BeforeRelativeAbundance),
    AfterRelativeAbundance = sum(AfterRelativeAbundance),
    HRICContribution = sum(HRICContribution),
    SimpsonDirectionContribution =
      sum(SimpsonDirectionContribution),
    PielouEntropyContribution =
      sum(PielouEntropyContribution),
    .groups = "drop"
  )

plot_contributions <- taxon_contributions %>%
  select(
    Transition, TransitionOrder, Taxon,
    HRICContribution,
    SimpsonDirectionContribution,
    PielouEntropyContribution
  ) %>%
  pivot_longer(
    cols = c(
      HRICContribution,
      SimpsonDirectionContribution,
      PielouEntropyContribution
    ),
    names_to = "Metric",
    values_to = "Contribution"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      HRICContribution = "HRIC alpha",
      SimpsonDirectionContribution = "Inverse Simpson",
      PielouEntropyContribution = "Pielou evenness"
    )
  ) %>%
  group_by(Transition, Metric) %>%
  mutate(
    SignedShare = 100 * Contribution / sum(abs(Contribution))
  ) %>%
  ungroup()

pielou_richness_plot <- decomposition_checks %>%
  transmute(
    Transition,
    TransitionOrder,
    Taxon = "Richness normalization",
    Metric = "Pielou evenness",
    Contribution = PielouRichnessNormalizationComponent
  )

pielou_denominators <- bind_rows(
  taxon_contributions %>%
    group_by(Transition) %>%
    summarise(
      AbsoluteTaxonContribution =
        sum(abs(PielouEntropyContribution)),
      .groups = "drop"
    ),
  decomposition_checks %>%
    transmute(
      Transition,
      AbsoluteTaxonContribution =
        abs(PielouRichnessNormalizationComponent)
    )
) %>%
  group_by(Transition) %>%
  summarise(
    PielouAbsoluteTotal = sum(AbsoluteTaxonContribution),
    .groups = "drop"
  )

plot_contributions <- plot_contributions %>%
  filter(Metric != "Pielou evenness") %>%
  bind_rows(
    plot_contributions %>%
      filter(Metric == "Pielou evenness") %>%
      select(
        Transition, TransitionOrder, Taxon,
        Metric, Contribution
      ) %>%
      bind_rows(pielou_richness_plot) %>%
      left_join(pielou_denominators, by = "Transition") %>%
      mutate(
        SignedShare =
          100 * Contribution / PielouAbsoluteTotal
      ) %>%
      select(
        Transition, TransitionOrder, Taxon,
        Metric, Contribution, SignedShare
      )
  )

taxon_ranking <- plot_contributions %>%
  filter(Taxon != "Richness normalization") %>%
  select(Transition, Taxon, Metric, SignedShare) %>%
  pivot_wider(
    names_from = Metric,
    values_from = SignedShare,
    values_fill = 0
  ) %>%
  mutate(
    DivergenceScore = pmax(
      abs(`HRIC alpha` - `Inverse Simpson`),
      abs(`HRIC alpha` - `Pielou evenness`)
    )
  ) %>%
  group_by(Transition) %>%
  slice_max(DivergenceScore, n = 8, with_ties = FALSE) %>%
  ungroup()

recurrent_taxa <- taxon_ranking %>%
  group_by(Taxon) %>%
  summarise(
    TotalDivergenceScore = sum(DivergenceScore),
    MeanDivergenceScore = mean(DivergenceScore),
    DiscordantTransitions = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalDivergenceScore))

selected_plot_data <- plot_contributions %>%
  semi_join(
    taxon_ranking %>% select(Transition, Taxon),
    by = c("Transition", "Taxon")
  ) %>%
  bind_rows(
    plot_contributions %>%
      filter(Taxon == "Richness normalization")
  ) %>%
  left_join(
    taxon_ranking %>%
      select(Transition, Taxon, DivergenceScore),
    by = c("Transition", "Taxon")
  ) %>%
  mutate(
    DivergenceScore = ifelse(
      Taxon == "Richness normalization",
      Inf,
      DivergenceScore
    )
  ) %>%
  group_by(Transition) %>%
  mutate(
    TaxonOrder = dense_rank(desc(DivergenceScore)),
    TaxonPanel = factor(
      paste(Transition, Taxon, sep = "___"),
      levels = rev(unique(
        paste(
          Transition[order(TaxonOrder, Taxon)],
          Taxon[order(TaxonOrder, Taxon)],
          sep = "___"
        )
      ))
    )
  ) %>%
  ungroup() %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c(
        "HRIC alpha",
        "Inverse Simpson",
        "Pielou evenness"
      )
    ),
    Transition = factor(
      Transition,
      levels = transition_data$Transition
    )
  )

metric_colors <- c(
  "HRIC alpha" = "#222222",
  "Inverse Simpson" = "#2F6DA3",
  "Pielou evenness" = "#D65A45"
)

discordance_plot <- ggplot(
  selected_plot_data,
  aes(
    x = SignedShare,
    y = TaxonPanel,
    color = Metric
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    color = "grey50"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = SignedShare,
      yend = TaxonPanel
    ),
    linewidth = 0.62,
    alpha = 0.76
  ) +
  geom_point(size = 2.15, stroke = 0) +
  facet_wrap(
    vars(Transition),
    ncol = 2,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = function(x) sub("^.*___", "", x)
  ) +
  scale_color_manual(values = metric_colors, drop = FALSE) +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.06, 0.08))
  ) +
  labs(
    title = "Taxa underlying discordant alpha-diversity trajectories",
    subtitle = paste0(
      "Positive values drive an increase in the named metric; ",
      "negative values drive a decrease."
    ),
    x = "Signed share of absolute metric contributions",
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
    axis.text = element_text(color = "grey15"),
    strip.background = element_rect(
      fill = "grey96",
      color = "grey25",
      linewidth = 0.45
    ),
    strip.text = element_text(
      face = "bold",
      size = 10.5,
      margin = margin(5, 5, 5, 5)
    ),
    panel.spacing = grid::unit(8, "mm"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "grey25"),
    legend.position = "top",
    legend.justification = "left",
    legend.key.width = grid::unit(5, "mm"),
    legend.margin = margin(0, 0, 5, 0)
  )

figure_base <- file.path(
  figure_dir,
  "figure_metric_discordance_taxon_contributions_0020_0061"
)
ggsave(
  paste0(figure_base, ".pdf"),
  discordance_plot,
  width = 10.0,
  height = 11.0,
  units = "in",
  device = cairo_pdf
)
ggsave(
  paste0(figure_base, ".png"),
  discordance_plot,
  width = 10.0,
  height = 11.0,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  paste0(figure_base, ".tiff"),
  discordance_plot,
  width = 10.0,
  height = 11.0,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

write.csv(
  transition_data,
  file.path(
    table_dir,
    "metric_discordance_transitions_0020_0061.csv"
  ),
  row.names = FALSE
)
write.csv(
  feature_contributions,
  file.path(
    table_dir,
    "metric_discordance_feature_contributions_0020_0061.csv"
  ),
  row.names = FALSE
)
write.csv(
  taxon_contributions,
  file.path(
    table_dir,
    "metric_discordance_taxon_contributions_0020_0061.csv"
  ),
  row.names = FALSE
)
write.csv(
  decomposition_checks,
  file.path(
    table_dir,
    "metric_discordance_decomposition_checks_0020_0061.csv"
  ),
  row.names = FALSE
)
write.csv(
  taxon_ranking %>%
    arrange(Transition, desc(DivergenceScore)),
  file.path(
    table_dir,
    "metric_discordance_top_taxa_0020_0061.csv"
  ),
  row.names = FALSE
)
write.csv(
  recurrent_taxa,
  file.path(
    table_dir,
    "metric_discordance_recurrent_taxa_0020_0061.csv"
  ),
  row.names = FALSE
)

message(
  "Metric-discordance taxon analysis completed.\n",
  "Discordant transitions: ", nrow(transition_data), "\n",
  "Maximum HRIC reconstruction error: ",
  format(
    max(abs(
      decomposition_checks$DeltaHRIC -
        decomposition_checks$ReconstructedDeltaHRIC
    )),
    scientific = TRUE
  ),
  "\nMaximum Pielou reconstruction error: ",
  format(
    max(abs(
      decomposition_checks$DeltaPielou -
        decomposition_checks$ReconstructedDeltaPielou
    )),
    scientific = TRUE
  ),
  "\nFigure: ", paste0(figure_base, ".pdf")
)
