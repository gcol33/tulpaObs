# Random-effect estimates (BLUPs) for a tobs_fit

Returns the per-group random-effect posterior summaries. Under NUTS the
non-centred draws are reconstructed to the natural BLUP scale
(`b_{g,c} = sigma_c * (L z_g)_c`) and summarised; the deterministic
Laplace path returns the variance-component EM modes and their
Schur-complement standard errors. A community fit reports its
per-species deviations from the community mean. Every family and engine
returns the same columns.

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

A data frame with one row per random-effect level and coefficient and
the columns `arm` (the linear predictor the effect enters: `psi`, `p`,
`lambda`, ...; `psi+p` for an effect shared across arms; `NA` when the
fit does not distinguish arms), `group` (the grouping variable;
`"species"` for the per-species deviations of a community fit), `level`
(the group level, or the species), `term` (the coefficient the effect
shifts, `"(Intercept)"` for a random intercept), `estimate` and
`std.error` (`NA` where the fitter reports no posterior SD). Zero rows
when the fit carries no random effects.

On a community fit, `std.error` is the per-species posterior SD of the
deviation, the square root of the matching diagonal entry of the
species' posterior covariance `fit$ms_community$Cinv[[s]]` (the analogue
of lme4's conditional variance). On the Laplace routes that covariance
is conditional on the community means; on a spatial fit integrated over
an outer hyperparameter grid it combines the per-node covariances and
the spread of the per-node BLUPs (law of total covariance); on the NUTS
and Polya-Gamma Gibbs routes it is the covariance of the sampled
deviations. The
[`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md)
`log_r` (negative binomial) and `logit_omega` (zero-inflated) deviations
carry an SD as well: every route that fits them (joint AGHQ, and NUTS
for `log_r`) keeps them in the covariance. `std.error` is `NA` for a
species whose per-species posterior solve failed.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 150, J = 4, n.occ.covs = 1, n.det.covs = 1,
                     n.visit.groups = 8, sigma.visit = 1, seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ (1 | visit_group), visits = sim$visits,
            y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
ranef(fit)
# }
```
