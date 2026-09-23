# Decode the two-arm fit into a cover_fit object

Extracts beta vectors and SEs for each arm. SEs are scaled to match each
arm's dispersion convention:

## Usage

``` r
decode_cover_hurdle(fits, enc, family, approx = "gaussian_laplace")
```

## Details

- lognormal arm: `tulpa_laplace(family = "gaussian")` computes the
  Hessian assuming phi = 1, so SEs are rescaled by `sigma_pos^2`.

- beta arm: `tulpa_laplace_beta()` already weights the Hessian by phi
  (Fisher information), so SEs are returned at scale 1.

Under an SLA method (`method = "laplace_sla"` / `"nested_laplace_sla"`),
the cover-hurdle SLA gamma is computed via
[`.sla_compute_cover_hurdle()`](https://gillescolling.com/tulpaObs/reference/dot-sla_compute_cover_hurdle.md):
a per-arm 5-point FD of the *original* Bernoulli / Beta / Lognormal
log-likelihood against the arm's `solve(H_beta)` Sigma (raw Hessian – no
Louis correction needed here because both arms run real likelihoods at
the mode, not the pseudo-binomial M-step encoding). Per-arm pseudo-draws
are then resampled from skew-normals fit by moment-matching
`(beta_arm, se_arm, gamma_arm)`.
