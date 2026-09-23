# Multispecies N-mixture family

Per-species N-mixture with shared community-level hyperparameters.

## Usage

``` r
ms_abun(K_max = NULL, mixture = c("poisson", "negbin", "zip", "zinb"))
```

## Arguments

- K_max:

  upper bound for the latent-abundance marginal sum (the exact
  integration over `N` is truncated at `K_max`). `NULL` (default) lets
  the engine pick `max(y) + 100` (matching
  [`unmarked::pcount()`](https://ecoverseR.github.io/unmarked/reference/pcount.html));
  raise it if a fit warns that the posterior over `N` puts mass on the
  boundary.

- mixture:

  Abundance mixing distribution: `"poisson"` (default), `"negbin"`, or
  their zero-inflated counterparts `"zip"` / `"zinb"`. The
  negative-binomial dispersion and the zero-inflation structural-zero
  share are each a per-species random effect
  (`log_r_s ~ N(mu_log_r, sigma_log_r)`,
  `logit_omega_s ~ N(mu_omega, sigma_omega)`) integrated by the
  community AGHQ path alongside the abundance / detection coefficients.
  Zero-inflation is a non-spatial `laplace` fit with an intercept-only
  structural-zero logit; a shared field stays Poisson / negbin.

## Value

A `tobs_family` object.

## Reading the negbin community dispersion

Under `mixture = "negbin"` the community mean `mu_log_r` reaches
[`coef()`](https://rdrr.io/r/stats/coef.html) /
[`vcov()`](https://rdrr.io/r/stats/vcov.html) /
[`confint()`](https://rdrr.io/r/stats/confint.html) with a marginal Wald
SE, and the variance component `sigma_log_r` is reported on
`fit$ms_dispersion`. That interval is calibrated **conditional on the
variance component being recovered**, and the two are not independent:
`sigma_log_r` is a scalar variance over species and at few species it
can settle near its lower boundary, which shifts `mu_log_r` and narrows
its SE at the same time.

Measured on simulated data (39 `laplace` fits at 8 and 36 species,
`sigma_logr = 0.5`): where `sigma_log_r` came back at least 0.30, the
nominal 95% interval covered 33 of 34 with a `sqrt(mean(z^2))` of 0.88;
where it came back below, it covered 2 of 5, the point estimate was 2.2x
further from the truth and the SE 28% narrower. So read
`fit$ms_dispersion$sigma_log_r` before quoting a `mu_log_r` interval: a
value near zero, well under what the data should support, means the
interval is not trustworthy even though the fit converged and reported
no warning. Those thresholds are the ones this fixture separated on and
are not a general rule; what transfers is the conditioning, not the
number.

Away from that boundary the interval is calibrated, and calibrated at
every group count measured. Over 97 fits at 8 / 18 / 36 species the
scale `sd(mu_log_r) / mean(se)` reads 1.066 / 1.265 / 0.981, and the
middle figure is the seed block rather than the estimator: `mu_log_r` is
a POPULATION mean, each seed draws `S` log-dispersions around it, and
that draw supplies about two thirds of the across-seed spread. The
18-species blocks drew theirs 18-21% wider than `sigma_logr / sqrt(S)`.
Put the draw at its expectation and the scale is 1.077 / 1.101 / 0.977;
rebuild the SE at the simulated sigma as well and it is 0.990 / 1.035 /
0.963 (NOTES_measurements.md). The residual there is the attenuation of
`sigma_log_r` itself (0.418 / 0.448 / 0.487 against 0.5), which shrinks
as species are added and which the SE inherits, since
`se^2 = (sigma_log_r^2 + c) / S`.

A simulation study measuring this needs
[`simulate_ms_abun()`](https://gillescolling.com/tulpaObs/reference/simulate_ms_abun.md)'s
`truth$mu_log_r_real` – the mean of the log-dispersions that seed
actually drew – beside `truth$mu_log_r`. Score coverage of a community
mean against the constant; score point recovery and interval SCALE
against the realized mean, or the draw is charged to the estimator.

`control$logr.sigma.prior` puts a Penalized-Complexity prior on that
variance, which adds curvature at `sigma_log_r -> 0`. Measured on the
same fixture at `c(1, 0.05)`, it is worth reaching for when a fit
reports a `sigma_log_r` near zero or fails outright with a singular
marginal Hessian – one seed that errored under pure maximum likelihood
converged and covered under the prior, and variance components near
0.01-0.06 lifted by an order of magnitude. It does **not** repair
coverage: on 19 paired seeds the nominal 95% interval covered 17 either
way, with the same two seeds missing, because the prior bites near the
boundary while those misses sit at `sigma_log_r` around 0.2. It is off
by default because the fits that were already calibrated pick up a small
systematic shift under it (`mu_log_r` by -0.006, p = 0.001).

## Examples

``` r
f <- ms_abun(mixture = "negbin")
f
```
