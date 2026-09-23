# Encode cover-hurdle data for the two-Laplace fit

Splits `y` into a binomial occurrence indicator and a positive-cover
subset, builds design matrices for each arm using the fixed-effects part
of `formula`, and extracts any structured terms it carried (an areal
[`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)/[`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
spatial field, plus
[`temporal()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
/ [`re()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
blocks). The formula is parsed against the NA-dropped observations so
the structured index codes align with both arms.

## Usage

``` r
encode_cover_hurdle(
  formula,
  data,
  y,
  positive = c("lognormal", "lognormal_trunc", "beta", "beta_oi", "ordinal", "gaussian"),
  breaks = NULL,
  autoscale = TRUE,
  presence_formula = NULL,
  positive_formula = NULL
)
```

## Arguments

- formula:

  State-process formula (no LHS); used for both occurrence and
  positive-cover arms.

- data:

  Data frame with `nrow(data) == length(y)`.

- y:

  Length-N numeric vector of cover in `[0, 1]`. NAs are dropped from
  both arms (treated as missing, not as zero cover).

- positive:

  Positive-arm likelihood: one of `"lognormal"`, `"lognormal_trunc"`,
  `"beta"`, `"beta_oi"`, `"ordinal"`, `"gaussian"`.

- breaks:

  Ordinal cut points, required when `positive = "ordinal"` and ignored
  otherwise.

- autoscale:

  Logical (default `TRUE`); rescale the positive arm away from the
  boundary so the Laplace engine never sees an exact 0 or 1.

- presence_formula:

  Optional formula overriding `formula` on the occurrence arm. `NULL`
  uses `formula` for both arms.

- positive_formula:

  Optional formula overriding `formula` on the positive-cover arm.
  `NULL` uses `formula` for both arms.

## Value

A list with: `occ_data`, `pos_data`, `spatial_spec` (a `tulpa_spatial`
built from the unweighted areal formula term, or NULL), `trend`
(per-observation weight + label from a weighted areal term, or NULL),
`temporal` and `re` (structured terms from the formula, or NULL), `N`,
`idx_pos` (row indices of the positive subset within `data`), `formula`
(the fixed-effects formula), `positive`.

## Details

For `positive = "lognormal"` the positive arm's response is
`log(y[occur == 1])`. For `positive = "beta"` it is `y[occur == 1]` on
the natural (0, 1) scale; an additional eps-clip is applied so the
Laplace engine does not see exact 1's introduced upstream.
