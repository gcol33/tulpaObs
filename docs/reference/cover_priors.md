# Weakly-informative priors for cover-hurdle Laplace fits

Builds an opt-in fixed-effect prior for
[`cover()`](https://gillescolling.com/tulpaObs/reference/cover.md)
models fit with a Laplace method. The cover hurdle has two arms with
their own coefficients – an occurrence arm (Bernoulli presence) and a
positive-cover arm (beta or lognormal) – each with an intercept and a
slope bucket. The penalty is the quadratic
`+ sum((beta_j - mu_j)^2 / (2 sd_j^2))` added to each arm's negative
log-likelihood, the same mechanism
[`occu_priors()`](https://gillescolling.com/tulpaObs/reference/occu_priors.md)
uses.

## Usage

``` r
cover_priors(
  occ_intercept = list(mean = 0, sd = 2),
  occ_slope = list(mean = 0, sd = 2.5),
  pos_intercept = list(mean = 0, sd = 3),
  pos_slope = list(mean = 0, sd = 2.5)
)
```

## Arguments

- occ_intercept:

  Prior on the occurrence (presence) intercept, logit scale.
  `list(mean, sd)`. Default `list(mean = 0, sd = 2)`.

- occ_slope:

  Prior on occurrence slopes. Default `list(mean = 0, sd = 2.5)`.

- pos_intercept:

  Prior on the positive-cover intercept (logit scale for `"beta"`, log
  scale for `"lognormal"`). Default `list(mean = 0, sd = 3)`.

- pos_slope:

  Prior on positive-cover slopes. Default `list(mean = 0, sd = 2.5)`.

## Value

A `cover_priors` object, ready to pass to `tobs(..., priors = ...)`.

## Details

Cover priors are opt-in: without `priors`, a cover fit is unpenalised.
The main use is regularising perfect separation in the occurrence arm at
small `N`. Set any `sd = Inf` to drop that component.

## Coverage

Applied on the separate-Laplace path (`method = "laplace"` /
`"laplace_sla"`) when the formula has no spatial term. Both arms are
penalised: the occurrence and lognormal-positive arms through
[`tulpa::tulpa_laplace()`](https://gillescolling.com/tulpa/reference/tulpa_laplace.html),
and the beta-positive arm through
[`tulpa::tulpa_laplace_beta()`](https://gillescolling.com/tulpa/reference/tulpa_laplace_beta.html).
The joint nested-Laplace path (`method = "nested_laplace"` /
`"nested_laplace_sla"`, used for spatial cover formulas) also threads
the priors, as a per-arm `beta_prior_mean` / `beta_prior_prec` on the
joint engine's occurrence and positive responses. The separate-Laplace
path keeps rejecting a prior alongside a spatial term (that solver
carries its own); add the spatial term through the nested-Laplace path
to combine the two.

## See also

[`occu_priors()`](https://gillescolling.com/tulpaObs/reference/occu_priors.md)

## Examples

``` r
# regularise the occurrence arm, leave the positive arm unpenalised
priors <- cover_priors(pos_intercept = list(mean = 0, sd = Inf),
                       pos_slope     = list(mean = 0, sd = Inf))
priors

# \donttest{
sim <- simulate_cover(N = 150, seed = 1)
fit <- tobs(cover ~ x, data = sim$data,
            family = cover(response = "lognormal"),
            method = "laplace", priors = priors)
coef(fit)
# }
```
