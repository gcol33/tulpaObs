# Parse a process formula into fixed effects and structured terms

Parse a process formula into fixed effects and structured terms

## Usage

``` r
.tobs_parse_formula(formula, data = NULL, env = environment(formula))
```

## Arguments

- formula:

  a one- or two-sided formula for a single process (e.g. the occupancy
  or detection linear predictor).

- data:

  the model data frame; columns are made available when evaluating
  structured-term calls.

- env:

  environment in which to resolve non-data symbols (e.g. an adjacency
  matrix passed as `graph = adj`). Defaults to the formula's
  environment.

## Value

A list with `fe_formula` (the fixed-effects formula for `model.matrix`)
and `terms` (a list of `tobs_*` specs).
