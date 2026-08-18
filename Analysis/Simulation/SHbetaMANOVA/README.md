# SHbeta-MANOVA simulation figure sources

This directory contains a local, reproducible copy of the results needed to
redraw the SHbeta-MANOVA relationship figures.

- `hpc_source/` contains unmodified copies of the HPC simulation and summary
  scripts.
- `source_data/` contains unmodified HPC output tables.
- `derived_data/` contains fitted curves generated locally.
- `figures/` contains the locally redrawn figures.

HPC source:

`/home/zhang.16383/AHC/beta_div_sparse12_shbeta_hric_manova_out/summary_hric_manova`

The two local relationship figures use
`table_replication_level_hric_manova.csv`, in which the independent simulation
replicate is the plotting and modelling unit. The source table contains 7,875
rows: 125 replicates for each of 21 effect sizes and three signal modes.

Panels b and c read these HPC-exported numerical values directly; no values
were digitized or estimated from the original HPC PDF figures. Within each
replication, the HPC summary script first averaged the 12 signal-cluster
results. The locally fitted lines and confidence intervals are derived display
layers and are not used as source observations.

The SHA-256 digest of the downloaded replication-level table is
`7a84bc24b2b0c0980eabbf664dda18160283fbae17e5d052384cc94594cadcc9`.
As a consistency check, recomputing all 63 signal-mode-by-effect-size means
from this table reproduced the independently downloaded HPC trend table to
floating-point precision (maximum absolute discrepancy below
`5.2e-15`). The exact SHbeta--sum-of-squares check had a maximum absolute
error of `6.1e-14`.
