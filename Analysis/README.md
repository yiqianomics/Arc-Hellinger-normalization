# HRIC analyses

This directory contains the complete analysis inputs, code, supporting tables,
and four manuscript figures for the ocean, auto-FMT, and simulation analyses.
Each analysis has one R script and recreates its own `results` directory from
the retained source data.

## Reproduction

Run the analyses from this directory:

```sh
Rscript RealDataOcean/analysis.R
Rscript RealDataGut/analysis.R
Rscript Simulation/analysis.R
```

The auto-FMT analysis performs 4,999 patient-level bootstrap iterations and
therefore takes longer than the other two analyses.

The analyses require R, the current HRIC package, and the packages checked at
the beginning of each script. HRIC can be installed directly from its source
repository:

```r
remotes::install_github("yiqianomics/HRIC")
```

The final verification run used R 4.5.0 and HRIC 0.1.0.

## Ocean analysis

`RealDataOcean/Ocean.RData` contains the phyloseq object `ps`, comprising Tara
Oceans microbial profiles and matched environmental metadata. The script uses
`HRIC::SHalpha`, `HRIC::SHgamma`, and `HRIC::SHbeta` directly for the diversity
partition and `HRIC::HRIC` for taxon-coordinate models. It writes:

- `RealDataOcean/results/figure_integrated_ocean_regions.pdf`
- the plotted diversity, environmental-association, taxon-model, composition,
  colour-key, and sample-order values in `RealDataOcean/results/tables`

## Auto-FMT analysis

`RealDataGut/gut.RData` contains the phyloseq object `gut` from the randomized
auto-FMT study reported by Taur et al. (Science Translational Medicine, 2018;
doi:10.1126/scitranslmed.aap9489). The script uses `HRIC::SHalpha`,
`HRIC::SHgamma`, `HRIC::SHbeta`, `HRIC::SHdelta`, and `HRIC::HRIC` directly. It
writes:

- `RealDataGut/results/figure1_fmt_hric_patients_0009_0028.pdf`
- `RealDataGut/results/figure2_fmt_hric_randomized_comparison.pdf`
- the plotted trajectories, diversity partitions, model contrasts, bootstrap
  draws, recovery summaries, and taxon contributions in
  `RealDataGut/results/tables`

Exact index-day stools are excluded from the primary post-index recovery
models because the metadata do not establish whether collection occurred
before or after auto-FMT. Bootstrap resampling is stratified by randomized arm
and retains all longitudinal observations for each sampled patient.

## Simulation analysis

`Simulation/data` contains the three numerical inputs used by the integrated
simulation figure. The replication-level MANOVA table and theoretical identity
checks are direct HPC exports; no values were digitized or estimated from a
rendered figure. The legacy method key `AHC_euclidean` in the power table is
displayed as `HRIC distance`. The script writes:

- `Simulation/results/figure_integrated_simulation_hric_turnover.pdf`
- the plotted power and replication-level turnover/MANOVA values in
  `Simulation/results`

The script verifies the expected 7,875 independent simulation rows, 125
replicates per parameter combination, and the HRIC turnover--sum-of-squares
identity to numerical precision before completing.
