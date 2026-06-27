# ArcHellinger

ArcHellinger provides Hellinger-Riemann intrinsic coordinates and Simplex
Hellinger diversity summaries for sparse compositional data.

The package is designed for count, abundance, or compositional matrices where
rows are samples and columns are components, taxa, or features. Exact zeros are
allowed; no pseudo-count replacement is required.

## Main functions

- `HRIC()`: centered Hellinger-Riemann intrinsic coordinates.
- `RHRIC()`: reference-based Hellinger-Riemann intrinsic coordinates.
- `Simplex_Hellinger_alpha()`: normalized Simplex Hellinger dominance and
  evenness.
- `Simplex_Hellinger_beta()`: Simplex Hellinger beta diversity distances.

## Example

```r
X <- matrix(
  c(10, 20, 30,
    5, 15, 80,
    0, 10, 90),
  nrow = 3,
  byrow = TRUE
)

HRIC(X)
RHRIC(X)
Simplex_Hellinger_alpha(X, measure = "both")
Simplex_Hellinger_beta(X)
```
