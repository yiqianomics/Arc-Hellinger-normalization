# HRIC auto-FMT analysis

## Study frame

- **Cohort:** Adults undergoing allogeneic hematopoietic stem-cell
  transplantation (allo-HSCT) at Memorial Sloan Kettering Cancer Center.
- **Design:** Phase 2, single-center, open-label, 1:1 randomized controlled
  trial, stratified by cord-blood versus non-cord-blood graft source.
- **Auto-FMT arm:** 14 evaluatable patients assigned a single rectal retention
  enema prepared from their own stool collected before conditioning.
- **Control arm:** 11 evaluatable patients assigned no FMT.
- **Evaluatable population:** The first 31 randomized participants comprised
  16 treatment and 15 control patients. Six lacked a post-randomization stool
  within 14 days, leaving 14 treatment and 11 control patients in the paper's
  microbiome analysis and in this local dataset.

## Timeline

- **Allo-HSCT day 0** is the stem-cell infusion day, not the auto-FMT day.
- A negative allo-HSCT day precedes the stem-cell infusion.
- `randomizationday` is a post-engraftment trial landmark on the allo-HSCT-day
  scale. Treatment participants received auto-FMT at this intervention
  landmark; controls received no procedure.
- The local metadata contain no separate `fmt_day` or within-day collection
  time. The analysis therefore calls randomization day the
  randomization/FMT **index day** and excludes exact index-day stools from the
  primary post-intervention model. A sensitivity analysis includes them.

## Baseline adjustment

Baseline diversity should be adjusted because pre-existing diversity can
persist into the post-intervention trajectory. Two distinct baselines are used:

1. **Personal pre-HSCT target:** the earliest available sample before allo-HSCT
   day 0. This is the closest local-data proxy for the stored pre-conditioning
   community that the patient is intended to recover toward.
2. **Pre-index damage state:** the sample nearest the randomization/FMT index
   within days -28 to -1. This captures how depleted the microbiota was
   immediately before the randomized intervention.

The primary alpha model adjusts both alpha values, randomization day, and
within-patient repeated measurements:

`alpha ~ baseline alpha + pre-index alpha + randomization day + arm * (early time + late time) + (1 | patient)`

The composition model replaces the outcome with Euclidean distance in the
package-defined HRIC coordinates from the patient's personal pre-HSCT sample
and adjusts the corresponding pre-index distance. Patient-clustered HC2
standard errors are used because the random-intercept variance for this outcome
collapsed to zero.

## Figure 1

Each stool is a local community. Package functions `HRIC::SHalpha`,
`HRIC::SHgamma`, and `HRIC::SHbeta` are applied directly using one fixed
6,647-feature dictionary.

- Patient 0009 is metadata patient C1, a control patient with 24 stools.
- Patient 0028 is metadata patient T4, an auto-FMT patient with 15 stools.
- Dashed horizontal lines mark each patient's earliest available pre-HSCT
  alpha diversity.
- Any index-day stool is classified with the pre-index period because the
  metadata do not establish its within-day order relative to auto-FMT. The
  selected treatment example has no index-day stool.
- Patient 0028 was selected as the treatment example because sampling closely
  brackets the intervention. Auto-FMT occurred on allo-HSCT day 28, 8 days
  after engraftment. HRIC alpha rose from 0.0467 at day 26 (2 days before
  auto-FMT) to approximately 0.066 at day 29 and 0.0816 at day 36. HRIC
  distance to the personal baseline was 1.989 before auto-FMT; the two day-29
  stools had distances of 0.408 and 1.451. Thus, this example is selected for
  its tightly sampled alpha recovery, not as proof of uniform compositional
  recovery.
- Trial eligibility after engraftment required Bacteroidetes below 0.1% of
  total 16S by quantitative PCR. Figure 1 therefore includes an aligned
  Bacteroidetes trajectory. Its points are relative abundance calculated from
  the local 16S sequencing table, whereas the horizontal 0.1% line is the
  protocol's qPCR threshold shown only as clinical context. The assays are not
  treated as interchangeable.
- The alpha-beta-gamma bars satisfy the package identity
  `gamma = mean(alpha) + beta`.
- The Sankey-style ribbons decompose every sample's **alpha deficit**. With
  `Z = HRIC::HRIC(counts)` and
  `A_p = asin(sqrt(1 - 1 / p))`, a feature contributes `Z^2 / A_p^2`; the
  feature contributions sum numerically to `1 - alpha` within `1e-10`.
  Ribbons are aggregated to the deepest available genus or family. They are
  not relative abundance and should not be interpreted as abundance increases
  or causal effects.

The dominant displayed contributors across the two example patients are
Listeriaceae, Streptococcus, Actinomyces, Enterobacteriaceae, Lactobacillus,
Blautia, [Eubacterium], and Peptostreptococcaceae; all remaining taxa are
retained as "Other taxa." The accompanying CSV preserves the sample-specific
values.

## Figure 2

The primary longitudinal analysis uses post-index days 1-60 and separates
days 1-14 from days 15-60. Exact day-0 stools are excluded because the metadata
do not establish whether they were collected before or after auto-FMT.

The primary early-recovery estimand is the time-averaged adjusted arm
difference over days 1-30. It is calculated from the complete fitted treatment
contrast curve by trapezoidal integration, rather than selecting one favorable
day. Uncertainty is evaluated by 4,999 nonparametric patient-level bootstrap
resamples, stratified by randomized arm.

The day 1-30 alpha-diversity difference was +0.0235. Its model-based 95% CI was
0.0099 to 0.0371, whereas the patient-bootstrap interval was -0.0015 to 0.0373.
The corresponding HRIC distance-to-baseline difference was -0.496, with a
patient-bootstrap interval of -0.617 to -0.230. Thus, the strongest robust
early-window evidence is for recovery of personal composition; the alpha
signal is concordant but less precise under patient resampling.

After adjustment for both baselines and randomization day, the auto-FMT minus
control alpha contrast was +0.0219 at day 14 (95% CI 0.0024 to 0.0414,
P = 0.0278), +0.0136 at day 30 (95% CI 0.0005 to 0.0268, P = 0.0422), and
-0.0018 at day 60 (P = 0.867). The HRIC distance-to-baseline contrast was
-0.471 at day 14 (95% CI -0.786 to -0.156, P = 0.00334), -0.349 at day 30
(95% CI -0.530 to -0.167, P = 0.000165), and -0.119 at day 60 (P = 0.591);
negative values mean closer to the personal pre-HSCT composition.

These results support **earlier recovery/separation** in the auto-FMT arm,
followed by convergence by day 60. They do not establish a significantly
steeper continuous recovery slope: the arm-by-early-time interactions were
P = 0.208 for alpha and P = 0.269 for HRIC distance.

The taxon panel compares how much each taxon reduces the alpha deficit from
pre-index to early and later post-index samples. Positive values indicate a
reduced taxon-specific alpha deficit; negative values indicate an increased
deficit. Among taxa in the union of each arm's 50 largest absolute
contributors, the control and auto-FMT profiles are related but not identical:
Spearman rho is 0.58 early and 0.49 later. Five taxa overlap between the
arm-specific top-10 sets at each period, whose union contains 15 taxa.

## Sensitivity analyses

- Replacing the earliest pre-HSCT alpha with the mean of all available
  pre-HSCT alpha values gives nearly identical day-14 and day-30 contrasts.
- Including exact day-0 stools gives an alpha contrast of +0.0233 at day 14
  (P = 0.0217) and +0.0139 at day 30 (P = 0.0399), with no day-60 difference.
  The corresponding HRIC-distance contrasts remain negative at day 14
  (-0.361, P = 0.0283) and day 30 (-0.295, P = 0.000627).
- In an event-time model spanning days -28 to -1 and days 1-60, the pre-index
  arm difference in alpha slope was borderline (P = 0.058), whereas the
  distance pre-trend was not evident (P = 0.373). A difference-in-differences
  sensitivity that subtracts the adjusted arm difference at day -1 retains
  an early alpha advantage of +0.0580 (model-based 95% CI 0.0350 to 0.0809)
  and a distance advantage of -0.603 (-0.858 to -0.348).
- Among patients below 80% of their personal alpha baseline at the index day,
  sustained recovery by day 30 was observed for 6 of 13 auto-FMT patients and
  1 of 8 controls. This irregularly sampled endpoint is descriptive and the
  two-sided Fisher exact test is not significant (P = 0.174).

## Remaining limitations

The local analysis uses the 25-patient evaluatable microbiome subset rather
than all 31 randomized participants. Selection based on post-randomization
stool availability can weaken a strict intention-to-treat interpretation.
Stool sampling is also irregular and may be informative, and the metadata do
not identify the within-day order of auto-FMT and stool collection.

The local metadata also lack time-varying antibiotics, diet, and G-CSF. These
are post-randomization variables and should not automatically enter the primary
model for the **total** effect of treatment, because they may be mediators.
They are useful for clearly defined sensitivity or direct-effect estimands.
For a definitive statement about recovery **rate**, the strongest additions
would be all randomized participants, exact intervention and sample times, a
prospectively specified recovery threshold, and a model that handles irregular
or informative stool collection.

## Source

Taur Y, et al. *Reconstitution of the gut microbiota of antibiotic-treated
patients by autologous fecal microbiota transplant*. Science Translational
Medicine. 2018;10(460):eaap9489.
https://doi.org/10.1126/scitranslmed.aap9489
