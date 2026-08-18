# Integrated simulation figure

`plot_simulation_integrated_nature.R` regenerates all three displays from local
source data and combines them into a 180 x 170 mm figure:

- **a**, HRIC MANOVA pseudo-F versus normalized between-group turnover;
- **b**, HRIC between-group sum of squares versus normalized between-group
  turnover;
- **c**, power under abundance, joint and prevalence perturbations.

The power source data are in `../BetaDiversity/weighted_power_output/`. The
SHbeta-MANOVA source data and unmodified HPC scripts are in
`../SHbetaMANOVA/`.

The manuscript-ready outputs are named
`figure_integrated_simulation_hric_turnover` in PDF, PNG and TIFF formats.
