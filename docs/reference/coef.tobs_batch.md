# Per-species coefficients from a batched occu_cover fit.

Per-species coefficients from a batched occu_cover fit.

## Usage

``` r
# S3 method for class 'tobs_batch'
coef(object, ...)
```

## Arguments

- object:

  A `tobs_batch` returned by
  [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) on
  multi-response `y`.

- ...:

  Passed to each per-species
  [`coef()`](https://rdrr.io/r/stats/coef.html).

## Value

A named list of per-species coefficient vectors.

## Examples

``` r
# \donttest{
sims <- lapply(1:2, function(s) simulate_occu_cover(N = 40, J = 3, seed = s))
y <- lapply(sims, function(s) s$y)
y_pos <- lapply(sims, function(s) { p <- s$y_pos; p[is.na(p)] <- 0; p })
fit <- tobs(~ occ_cov1, data = sims[[1]]$data,
            family = occu_cover("lognormal"), detection = ~ 1,
            y = y, y.pos = y_pos, species = c("a", "b"), method = "laplace")
coef(fit)
# }
```
