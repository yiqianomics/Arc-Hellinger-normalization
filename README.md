# Hellinger–Riemann Intrinsic Coordinates (HRIC)

Intrinsic coordinates and diversity summaries for sparse compositional data,
built on the Hellinger–Riemann geometry of the closed simplex.

HRIC maps count, abundance, or compositional data — where rows are samples and
columns are components, taxa, or features — into an intrinsic coordinate system
in which ordinary Euclidean operations (distances, means, dispersion) respect the
geometry of the simplex. Exact zeros are handled directly, so no pseudo-count
replacement is required.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yiqianomics/HRIC")
```

## Functions

### Coordinates

- `HRIC(X)` — centered Hellinger–Riemann intrinsic coordinates. Returns an
  `n x p` matrix.
- `RHRIC(X, reference)` — reference-based coordinates relative to a single
  component. Returns an `n x (p - 1)` matrix; the reference defaults to the last
  column.

### Diversity

All diversity summaries are computed from the HRIC coordinates and share a common
evenness scale, where 1 is uniform abundance across components and 0 is
concentration on a single component.

- `SHalpha(X)` — per-sample alpha diversity (evenness). Returns one value per
  sample in `[0, 1]`.
- `SHgamma(X)` — gamma diversity: the evenness of the cohort center. Returns a
  single value.
- `SHbeta(X)` — beta diversity: the additive partition, gamma minus the mean
  sample-level alpha. Returns a single non-negative value.
- `SHdelta(X, group)` — between-group turnover: the coordinate dispersion that
  remains after removing the group-size weighted within-group dispersion.
  Requires a grouping vector and generalizes to any number of groups.

## Example

```r
library(HRIC)

# counts: rows = samples, columns = taxa (exact zeros allowed)
X <- matrix(
  c(10, 20, 30,
     5, 15, 80,
     0, 10, 90,
    40,  5,  5),
  nrow = 4, byrow = TRUE,
  dimnames = list(paste0("sample", 1:4), paste0("taxon", 1:3))
)

# Intrinsic coordinates
HRIC(X)
RHRIC(X)

# Diversity summaries
SHalpha(X)                 # alpha diversity (per-sample evenness)
SHgamma(X)                 # gamma diversity (evenness of the cohort center)
SHbeta(X)                  # beta diversity (diversity among samples)

group <- c("A", "A", "B", "B")
SHdelta(X, group)          # between-group turnover
```

## Diversity framework

The alpha, gamma, and beta summaries form a Whittaker-style additive partition on
the Simplex Hellinger evenness scale:

```
gamma = mean(alpha) + beta
```

Here alpha is within-sample evenness, gamma is the evenness of the pooled cohort
center, and beta is the diversity among samples. Because beta equals the mean
squared coordinate deviation from the cohort center (rescaled to the evenness
range), it is always non-negative.

Given a grouping, `SHdelta` provides a separate variance decomposition of
coordinate dispersion and returns the between-group component. A positive value
means turnover remains after accounting for within-group variation.

## Notes

- Input may be a numeric matrix or a data frame. Rows are samples; columns are
  components, taxa, or features.
- Values must be non-negative and each row must have a positive total. Exact
  zeros are allowed; missing values are not and should be handled beforehand.
- The package depends only on base R.

## Citation

A manuscript describing the method is in preparation. Please check back for
citation details.