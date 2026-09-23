# Fit the two arms of a cover hurdle

Two independent
[`tulpa::tulpa_laplace()`](https://gillescolling.com/tulpa/reference/tulpa_laplace.html)
calls. For `positive = "lognormal"` the positive arm is a Gaussian fit
on `log(cover)` with sigma estimated post-hoc as the residual standard
error. For `positive = "beta"` the positive arm uses
[`tulpa::tulpa_laplace_beta()`](https://gillescolling.com/tulpa/reference/tulpa_laplace_beta.html)
which estimates the precision `phi` via an outer 1-D optimisation and
weights the Hessian accordingly.

## Usage

``` r
fit_cover_hurdle(
  enc,
  positive = enc$positive,
  engine = "laplace",
  priors = NULL,
  control = list()
)
```

## Arguments

- enc:

  Output of
  [`encode_cover_hurdle()`](https://gillescolling.com/tulpaObs/reference/encode_cover_hurdle.md).

- positive:

  `"lognormal"` or `"beta"` (taken from `enc$positive`).

- engine:

  `"laplace"` (default) or `"nested_laplace"`. The latter is routed
  through
  [`fit_cover_hurdle_joint_nested()`](https://gillescolling.com/tulpaObs/reference/fit_cover_hurdle_joint_nested.md).

- priors:

  Optional
  [`cover_priors()`](https://gillescolling.com/tulpaObs/reference/cover_priors.md)
  object (or a coercible list / `FALSE`). Adds a weakly-informative
  fixed-effect penalty on both arms (occurrence + positive, beta or
  lognormal); `NULL` / `FALSE` fit unpenalised. Rejected with a spatial
  formula (the spatial solver carries its own prior).

- control:

  List with optional `max.iter`, `tol`, `n.threads`.

## Value

List with `m_occ`, `m_pos`, `positive`, `pos_fit_n`, `pos_fit_p`, plus
one of `sigma_pos` (lognormal) or `phi_pos` (beta).
