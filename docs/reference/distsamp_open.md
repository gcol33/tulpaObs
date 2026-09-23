# Open-population distance-sampling family

The `unmarked` `distsampOpen` model: a Dail-Madsen open N-mixture (as
[`dyn_abun()`](https://gillescolling.com/tulpaObs/reference/dyn_abun.md))
observed by distance sampling at each primary period – the
open-population counterpart of the single-season
[`gdistremoval()`](https://gillescolling.com/tulpaObs/reference/gdistremoval.md).
Initial abundance `N_1 ~ Poisson(lambda)`, then
`N_t = Binomial(N_{t-1}, omega) + Poisson(gamma)` (apparent survival
`omega` + recruitment `gamma`); at each primary period the detected
birds are distance-sampled into bins.

## Usage

``` r
distsamp_open(
  transect = c("line", "point"),
  cutpoints = NULL,
  K_max = NULL,
  mixture = c("poisson", "negbin", "zip", "zinb"),
  dynamics = c("constant", "notrend", "trend", "autoreg", "ricker", "gompertz")
)
```

## Arguments

- transect:

  Transect geometry: `"line"` (default) or `"point"`.

- cutpoints:

  Distance-bin edges, length `dim(y)[2] + 1`, strictly increasing and
  starting at `>= 0`.

- K_max:

  Truncation for the latent abundance HMM. Defaults to
  `3 * max(period total) + 40`.

- mixture:

  Initial-abundance mixing distribution. `"poisson"` (default),
  `"negbin"` (negative binomial, `Var(N_1) = lambda + lambda^2 / r`,
  with an overdispersion `r` estimated jointly and reported as `log_r`),
  `"zip"` (zero-inflated Poisson) or `"zinb"` (zero-inflated negative
  binomial). A structural-zero site is never occupied across any primary
  period, so all its band counts are zero; the observed per-site
  marginal is the two-component mixture
  `omega * 1{all zero} + (1 - omega) * L_open`, with `L_open` the exact
  open-population distance marginal and `omega` an intercept-only
  structural-zero probability reported as `zi_logit` (distinct from
  `omega`, the survival arm). The negative-binomial size /
  zero-inflation is layered over the same forward-HMM marginal; the
  Poisson path is unchanged.

- dynamics:

  Population-dynamics form for the transition `N_t | N_{t-1}` (following
  [`unmarked::distsampOpen()`](https://ecoverseR.github.io/unmarked/reference/distsampOpen.html)):
  `"constant"` (default) is
  `N_t = Binomial(N_{t-1}, omega) + Poisson(gamma)`; `"notrend"` ties
  recruitment so the expected abundance is stationary
  (`gamma = (1 - omega) * lambda`, no free `gamma` arm); `"trend"` is an
  exponential-growth process `N_t ~ Poisson(N_{t-1} * gamma)` with no
  survival term (`omega` unused); `"autoreg"` adds density-dependent
  recruitment `Poisson(N_{t-1} * gamma)` on top of survival; `"ricker"`
  and `"gompertz"` are the density-regulated recruitment forms with an
  estimated carrying capacity `K` (reported as `K`, on the log scale,
  modelled through the `omega = ~ ...` formula slot) and a growth rate
  reported as `r` (modelled through the `gamma = ~ ...` slot).
  `"constant"` / `"notrend"` fit with the exact analytic gradient; the
  density-dependent forms (`trend` / `autoreg` / `ricker` / `gompertz`)
  use the exact forward-HMM marginal with a numeric gradient.

## Value

A `tobs_family` object.

## Details

The distance-band allocation is conditional on the period total
detected, so it factors out of the abundance HMM (the
[`gdistremoval()`](https://gillescolling.com/tulpaObs/reference/gdistremoval.md)
trick): the marginal is the
[`dyn_abun()`](https://gillescolling.com/tulpaObs/reference/dyn_abun.md)
forward recursion with the detection probability set to the overall
distance detection `pdist`, plus the per-period band multinomials.

Four site-level arms: log abundance `lambda` (the
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
`formula`), log distance scale `sigma` (`detection`), logit survival
`omega` (`omega = ~ ...`), and log recruitment `gamma`
(`gamma = ~ ...`), the last two intercept-only by default.

## Inputs

`y` is a 3D array `[n_sites x n_bins x n_seasons]` of per-distance-band
counts at each primary period (secondary occasions absorbed into the
period total).

## References

Dail, D., Madsen, L. (2011). Models for estimating abundance from
repeated counts of an open metapopulation. *Biometrics* 67, 577-587.
Sollmann, R., Gardner, B., Chandler, R. B., Royle, J. A., Sillett, T. S.
(2015). An open-population hierarchical distance sampling model.
*Ecology* 96, 325-331.

## See also

[`gdistremoval()`](https://gillescolling.com/tulpaObs/reference/gdistremoval.md)
(single-season),
[`dyn_abun()`](https://gillescolling.com/tulpaObs/reference/dyn_abun.md)
(open N-mixture),
[`distance()`](https://gillescolling.com/tulpaObs/reference/distance.md).

## Examples

``` r
f <- distsamp_open(transect = "line", cutpoints = c(0, 10, 20, 30, 40),
                   dynamics = "notrend")
f
```
