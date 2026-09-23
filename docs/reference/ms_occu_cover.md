# Community (multispecies) joint occupancy-detection + cover family

The community version of
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md):
a per-species joint occupancy-cover model with Gaussian community
hyperpriors on the per-species coefficients of all three arms (occupancy
`psi`, detection `p`, and positive cover). Rare species borrow strength
from common ones through the shared community means and covariances, the
same pooling that stabilises
[`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md)
and
[`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md)
but on the joint occupancy + vegetation-cover response.

## Usage

``` r
ms_occu_cover(response = c("beta", "lognormal", "gaussian"))
```

## Arguments

- response:

  likelihood for the positive cover arm. `"beta"` (cover in (0, 1)),
  `"lognormal"` (log-cover Gaussian), or `"gaussian"` (an identity-link
  Gaussian magnitude, the delta-normal hurdle; for a pre-transformed /
  unbounded positive response, not raw cover fractions).

## Value

A `tobs_family` object.

## Details

Per species `s`, cell `i`, visit `j`:

    z_{s,i}        ~ Bernoulli(psi_{s,i})
    y_{s,i,j} | z  ~ Bernoulli(p_{s,i,j})
    c_{s,i,j} | y  ~ f_pos(eta_pos_{s,i,j}, disp)
    logit psi_{s,i}  = X_occ_i . (mu_occ + b_occ_s)
    logit p_{s,i,j}  = X_p_{ij} . (mu_p   + b_p_s)
    g(cover)         = X_pos_{ij} . (mu_pos + b_pos_s)
    b_occ_s ~ N(0, Sigma_occ), b_p_s ~ N(0, Sigma_p),
    b_pos_s ~ N(0, Sigma_pos)

The latent presence `z` integrates out per species-cell in closed form
(the same two-state mixture as
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md));
the per-species coefficient deviations are the random effects,
integrated by a Laplace-EM. The positive-arm dispersion is a shared
community parameter.

## Inputs

`y` and `y_pos` are 3D arrays `[n_sites x max_visits x n_species]` (or
named lists of `n_sites x max_visits` matrices, one per species);
`species =` is required. The occupancy `formula`, the `detection`
formula, and the cover `positive` formula carry community covariates
shared across species. [`coef()`](https://rdrr.io/r/stats/coef.html)
returns the community means;
[`ranef()`](https://gillescolling.com/tulpa/reference/ranef.html) the
per-species coefficient deviations.

## Reduced-rank spatial factors

A single areal field term (`icar(graph = adj)`,
`car_proper(graph = adj)`, or `bym2(graph = adj)`) on the occupancy
`formula` fits the reduced-rank (HMSC / spatial-gllvm) community model:
`K` shared latent fields with per-species loadings on the occupancy
predictor, `logit psi_sc = X_c (mu_occ + b_occ_s) + sum_k L_sk w_kc`. A
single shared field with only a species intercept can shift each map's
level but not its shape; the loadings give each species its own spatial
shape as a combination of the shared factors, borrowing strength so rare
taxa get a calibrated map. The number of factors is set by
`control = list(n.factors = K)`, or chosen automatically by
`control = list(n.factors = "auto", n.factors.max = M)`, which fits the
identified (lower-triangular loading) ladder and picks the rank that
maximises the empirical-Bayes Laplace marginal likelihood (the field is
integrated out, so its prior supplies the Occam penalty – latent-level
criteria such as WAIC track the field's effective dimension, not the
rank, and under-select). The per-rank evidence is returned in
`fit$spatial$K_selection`. `sd.load` sets the loading prior scale.

Adding the SAME `icar(graph = adj)` term to the cover `positive` formula
shares the latent fields across the two processes: the fields then also
load on the cover predictor through a free loading matrix,
`g(cover)_sc = X_c (mu_pos + b_pos_s) + sum_k Lpos_sk w_kc`, so a
species' spatial occupancy pattern and its spatial cover pattern are
linked through one set of factors. The cover loadings are returned in
`fit$spatial$loadings_cover`. The field is shared, so the cover-arm term
must match the occupancy arm (same type and graph), and a cover-arm
field without a matching occupancy field errors.

The field structure is set by the term:
[`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
(improper intrinsic CAR, the default), `car_proper(graph = adj)` (a
proper CAR whose per-factor correlation `rho` is returned in
`fit$spatial$rho_w`), or `bym2(graph = adj)` (the Riebler 2016
convolution whose per-factor spatial-variance fraction `phi` is returned
in `fit$spatial$phi_w`); the field hyperparameter is estimated by EM.
All three share one engine; `fit$spatial$field_type` records the choice.
Structured terms on the detection arm, or an unsupported field term,
error from the dispatcher.

The shared fields imply a residual species-association matrix (the
spatial-JSDM / HMSC output):
[`tobs_associations()`](https://gillescolling.com/tulpaObs/reference/tobs_associations.md)
returns the occupancy association `corr(L L')` and, with a cover-arm
factor, the cover association and the cross-arm occupancy-vs-cover
association, each marginalised over the loading posterior so an interval
accompanies the estimate.
[`predict()`](https://rdrr.io/r/stats/predict.html) returns the other
JSDM output – the per-species per-cell occupancy maps (`psi` with
`psi_lower` / `psi_upper`), marginalised over the loading + field
posterior so a rare species borrows strength across the shared factors
for a calibrated map.

## Scope

The non-spatial fit is Laplace-EM. A community spatial occu_cover – a
shared latent field coupled across the occupancy and cover arms with
per-species RE on all three arms – is the reduced-rank spatial-factor
fit: an
[`icar()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md) /
[`car_proper()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
/ [`bym2()`](https://gillescolling.com/tulpaObs/reference/tobs_terms.md)
term on the occupancy formula (and, with a matching term on the cover
formula, a cover-arm factor) fits per-species loadings on the shared
field via Laplace-EM and NUTS, with the per-species community
covariances `Sigma_occ` / `Sigma_p` / `Sigma_pos` on top. The free
per-species loading form generalises a single common field amplitude, so
it subsumes the common-amplitude coupling of
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)'s
joint-coupled engine. A community model whose per-arm variance
components are integrated on an outer grid (the engine route of
[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)'s
joint-coupled path) is not used: the joint nested-Laplace engine
integrates every variance component on its outer grid, so per-arm
community RE variances plus the field hyperparameters exceed the
engine's grid cap – the closed-form covariance M-step of the Laplace-EM
is the scaling route for community variance components. Structured terms
beyond the shared field error from the dispatcher rather than being
silently dropped.

## Community variance debias

The community-MEAN estimates
([`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html)) are unbiased. The
community-VARIANCE components – the per-arm covariance matrices in
`fit$ms_community$Sigma_occ` / `Sigma_p` / `Sigma_pos` and their `sd_*`,
which set the spread of the per-species deviations – pick up Laplace
small-cluster attenuation at small per-species n under the EM M-step.
They are debiased by default with an adaptive Gauss-Hermite quadrature
of the exact per-species RE posterior (`control$re.aghq = TRUE`, the
same correction the single-arm AGHQ path applies; disable with
`re.aghq = FALSE`, set nodes with `control$n.quad`). The cover hurdle
ties a species' psi / p / cover coefficients through the data, so the
per-species RE posterior is not separable across arms and the AGHQ is a
tensor product over the joint RE vector – `n.quad^P` nodes per species
in the total RE dimension `P`. That cost is exponential in `P`, so the
debias is a hard scope limit: it runs only up to a small total RE
dimension (`control$re.aghq.maxdim`, default 4, e.g. intercept-only on
all three arms), and larger RE designs keep the EM variance. The EM
(Laplace) variance is a documented lower bound on the true component –
attenuated toward zero at small per-species n, monotonically less
attenuated as n grows – not a bias of unknown sign.
`fit$ms_community$var_attenuation$debias` records `"aghq"` or `"none"`,
`$affects` names the lower-bounded components, and
[`print()`](https://rdrr.io/r/base/print.html) flags the EM case.

## See also

[`occu_cover()`](https://gillescolling.com/tulpaObs/reference/occu_cover.md)
(single species),
[`ms_occu()`](https://gillescolling.com/tulpaObs/reference/ms_occu.md)
(community occupancy, no cover),
[`ms_abun()`](https://gillescolling.com/tulpaObs/reference/ms_abun.md)
(community N-mixture).

## Examples

``` r
f <- ms_occu_cover(response = "lognormal")
f
```
