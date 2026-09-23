# Compute SLA gamma for cover hurdle arms

Returns per-coefficient SLA gamma for each of the two arms of a cover
hurdle fit. Each arm's gamma is computed via
[`.sla_gamma_fd()`](https://gillescolling.com/tulpaObs/reference/dot-sla_gamma_fd.md)
using the arm's `solve(H_beta)` as Sigma (no Louis correction – both
arms are real likelihoods at the mode under the single-Laplace path).

## Usage

``` r
.sla_compute_cover_hurdle(fits, enc, positive)
```

## Arguments

- fits:

  The list returned by
  [`fit_cover_hurdle()`](https://gillescolling.com/tulpaObs/reference/fit_cover_hurdle.md)
  (with `m_occ`, `m_pos`, `positive`, and one of `phi_pos` /
  `sigma_pos`).

- enc:

  The encoded data from
  [`encode_cover_hurdle()`](https://gillescolling.com/tulpaObs/reference/encode_cover_hurdle.md).

- positive:

  One of `"beta"`, `"lognormal"`.

## Value

List(gamma_occ, gamma_pos, valid, reason).

## Details

phi (beta arm) / sigma_pos (lognormal arm) are treated as fixed nuisance
parameters: SLA gamma is computed only for the fixed-effect
coefficients.
