# Decode a joint-nested-Laplace cover-hurdle fit into a `cover_fit`.

Lighter-weight than the single-Laplace decode: the joint engine has
already produced posterior moments for beta and the spatial
hyperparameters, so we just shape them into the existing `cover_fit`
structure.

## Usage

``` r
decode_cover_hurdle_joint(fits, enc, family, approx = "gaussian_laplace")
```

## Details

Under an SLA method (`method = "nested_laplace_sla"`), the per-arm
marginal skewness is computed via `.sla_compute_cover_hurdle_joint()`
(mixture third-moment over the outer grid; per-grid FD of the joint
inner log-lik along the constraint-corrected Sigma columns) and stored
with the grid-combined marginals it corrects, which
`.tobs_cover_eta_draws()` reshapes the joint posterior draw bundle
against.
