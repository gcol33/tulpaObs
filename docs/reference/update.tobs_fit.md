# Update and refit a tobs model

Re-enters
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) with
the arguments the fit was originally built from (`$tobs_call`, recorded
by every family), so the refit keeps the fit's family dispatcher,
`method` and `control` instead of dropping to a bare Laplace refit of
the occupancy state alone. A name in `...` that is one of
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)'s own
arguments (`formula`, `data`, `family`, `occurrence`, `detection`,
`positive`, `y`, `visits`, `method`, `priors`, `control`) replaces that
argument outright; any other name (`n.iter`, `seed`, `verbose`, ...) is
merged into `control`, matching where the recorded fit's own value for
it lives.

## Usage

``` r
# S3 method for class 'tobs_fit'
update(object, ..., evaluate = TRUE)
```

## Arguments

- object:

  A `tobs_fit` object (built via
  [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)).

- ...:

  Named arguments to override; see Details on where each lands.

- evaluate:

  If TRUE (default), refit the model. If FALSE, return the unevaluated
  [`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md) call.

## Value

Updated `tobs_fit` object (or an unevaluated call if
`evaluate = FALSE`).
