# Check model identifiability

Diagnostics for potential identifiability issues. Checks for confounded
covariates, low detection or detection-free sites, and sparse
replication.

## Usage

``` r
tobs_check_id(model, fit = NULL)
```

## Arguments

- model:

  A `tobs_model` object (before fitting).

- fit:

  Optional `tobs_fit` object (for post-fit diagnostics).

## Value

A list with `identifiable`, the `issues` messages, and `prefit_checked`
(`FALSE` where the family carries no site-by-visit grid).

## Details

Applies to the families whose response is a site-by-visit grid:
[`occu()`](https://gillescolling.com/tulpaObs/reference/occu.md),
[`abun()`](https://gillescolling.com/tulpaObs/reference/abun.md),
[`removal()`](https://gillescolling.com/tulpaObs/reference/removal.md)
and
[`fp_occu()`](https://gillescolling.com/tulpaObs/reference/fp_occu.md).
On any other family the pre-fit checks are reported as not applicable
rather than silently skipped – an empty `issues` on a family the checks
never ran on is not a clean bill of health. The post-fit sampler checks
apply to any NUTS fit.

## Examples

``` r
# \donttest{
sim <- simulate_occu(N = 100, J = 3, n_occ_covs = 1, n_det_covs = 1,
                     seed = 1)
fit <- tobs(~ occ_cov1, data = sim$data, family = occu(),
            detection = ~ det_cov1, y = sim$y, method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
chk <- tobs_check_id(fit$model, fit)
chk$identifiable
# }
```
