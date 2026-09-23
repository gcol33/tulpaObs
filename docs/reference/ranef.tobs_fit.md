# Random-effect estimates (BLUPs) for a tobs_fit

Returns the per-group random-effect posterior summaries. Under NUTS the
non-centred draws are reconstructed to the natural BLUP scale
(`b_{g,c} = sigma_c * (L z_g)_c`) and summarised; the deterministic
Laplace path returns the variance-component EM modes and their
Schur-complement standard errors. When the fit carries no formula random
effects this falls back to the generic flat random-effect table.

## Usage

``` r
# S3 method for class 'tobs_fit'
ranef(object, ...)
```

## Arguments

- object:

  A `tobs_fit` object.

- ...:

  Ignored.

## Value

A data frame with one row per group level and coefficient (`group`,
`level`, `term`, `estimate`, `std.error`), or the generic `ranef` table
when no `re_effects` are present.
