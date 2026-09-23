# Open-population (Dail-Madsen) N-mixture family

Latent abundance evolves across primary seasons via apparent survival
and recruitment (Dail & Madsen 2011): `N_1 ~ Poisson(lambda)`; for
`t >= 2`, `N_t = S_t + G_t` with survivors
`S_t ~ Binomial(N_{t-1}, omega)` and recruits `G_t ~ Poisson(gamma)`;
observed via `Binomial(N_t, p)` over secondary visits. The latent
abundance sequence is summed out by an exact HMM forward recursion (it
is not closed form, unlike the static
[`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md));
analytic gradients come from forward-mode differentiation of the scaled
forward algorithm, so the fit is a direct maximum-likelihood / Laplace
fit with a NUTS path over the same marginal.

## Usage

``` r
dyn_abun(K_max = NULL, mixture = c("poisson", "negbin", "zip", "zinb"))
```

## Arguments

- K_max:

  abundance-state truncation for the forward recursion (states
  `0..K_max`). `NULL` (default) uses `max(count) + 40`; raise it if
  abundance may exceed that (the forward cost is roughly cubic in
  `K_max`).

- mixture:

  initial-abundance distribution: `"poisson"` (default), `"negbin"`
  (negative-binomial `N_1 ~ NB(mean = lambda, size = r)`), or their
  zero-inflated counterparts `"zip"` / `"zinb"` (a structural-zero share
  `omega` of sites is never occupied across any season; the remaining
  sites follow the Dail-Madsen open-population process). Zero-inflation
  is non-spatial `laplace` with an intercept-only `omega`; a field / RE
  / NUTS stay Poisson / negbin.

## Value

A `tobs_family` object.

## Details

Four arms: initial abundance `lambda` (the
[`tobs()`](https://gillescolling.com/tulpaObs/reference/tobs.md)
`formula`) and detection `p` (`detection`) are site-level; apparent
survival `omega` (the `omega` argument, default `~ 1`) and recruitment
`gamma` (the `gamma` argument, default `~ 1`) span the `T - 1`
transition intervals. A constant `omega` / `gamma` is shared across a
site's seasons; supplying a season-varying covariate (a
`[n_sites x (T - 1)]` matrix column of `data`, one column per transition
interval) on `omega` / `gamma` gives interval-specific vital rates. The
response `y` is a 3D array `[n_sites x max_visits x n_seasons]` (or a
list of per-season count matrices); missing visits are `NA`.

## References

Dail, D., Madsen, L. (2011). Models for estimating abundance from
repeated counts of an open metapopulation. *Biometrics* 67, 577-587.

## Examples

``` r
f <- dyn_abun(mixture = "negbin")
f
```
