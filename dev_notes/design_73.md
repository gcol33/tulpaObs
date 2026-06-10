# tulpaObs#73 - ms_abun spatial NUTS (shared areal field on the community N-mixture)

## Goal

Combine the non-centered per-species community N-mixture sampler
(`src/ms_abun_nuts.cpp`, #14) with a SHARED fixed-hyper non-centered areal field
on the abundance arm (the `abun` NUTS+areal block in `src/marginal_count_nuts.h`
/ `R/abun_nuts.R::.tobs_fit_abun_nuts_spatial`, #51, the tulpa#87 pattern). The
result samples the exact joint posterior of the community means, the per-species
deviations, the community covariances, AND the shared spatial field jointly,
giving calibrated community + field intervals and per-(species, site) WAIC/LOO.

## Model

    log lambda_{s,i} = X_lambda_i . (mu_lambda + b_lambda_s) + f_{u(i)}
    logit p_{s,i,j}  = X_p_{ij}   . (mu_p      + b_p_s)
    b_{s,arm} ~ N(0, Sigma_arm),  f ~ N(0, [tau Q(rho)]^{-1})  shared across species

with N_{s,i} integrated out per species-site in closed form (the Royle marginal).

## Sampler parameterisation (both layers non-centered)

- Per-species block (unchanged from #14): whitened z_s ~ N(0, I),
  b_{s,arm} = C_arm z_{s,arm}, so the community covariances enter only the data
  term -> breaks the centered b<->Sigma funnel.
- Field block (mirrors #51): FIXED field precision tau Q(rho) at the
  nested-Laplace (#12 sfMsNMix) posterior mean; whitened raw ~ N(0, I),
  f = Linv %*% raw with Linv = (chol(tau Q(rho)))^{-1}. The hyperparameter is
  fixed, so there is no field-hyperparameter funnel and no log|Q(rho)| gradient
  (the tulpa#87 result: a fixed-hyper non-centered field calibrates the beta SDs
  to the nested-Laplace SEs). proper-CAR first: full-rank Q(rho) -> well-
  conditioned geometry; intrinsic ICAR has a flat field-mean direction that maxes
  treedepth without a sum-to-zero reparam (gate icar/bym2 to nested_laplace).

Packed vector:
    theta = ( mu, {z_s} species-major, chol_lambda, chol_p [, chol_logr], raw )
with `raw` length n_field_units appended after the chol blocks.

## Gradient

The field is a SHARED offset on every species' eta_lambda, so:

- log-posterior gains  sum_{s,i} [contribution through eta_lambda += f_{u(i)}]
  and the whitened-raw prior  -0.5 ||raw||^2.
- d log p / d f_u = sum_s sum_{i : u(i)=u} grad_eta_lambda_{s,i}  (the per-site
  abundance score already returned by compute_nmix_site, summed over species and
  over the sites mapping to unit u). Then
      d log p / d raw = Linv' %*% (d log p / d f) - raw .
- The per-species coefficient / chol gradients are UNCHANGED in form: each
  species sees its eta_lambda shifted by the current f, and grad_eta_lambda then
  flows into (mu_lambda, b_lambda_s) exactly as before. So no new likelihood
  math -- the field is one extra additive offset + one extra gradient channel.

## Implementation path (in-tree, single source of truth)

1. `MsNmixNutsData` (ms_abun_nuts.cpp): add the optional field block fields
   mirroring CountNutsData -- `n_field_units`, `o_raw`, `field_map` (0-based unit
   per site), `Linv` (row-major). Layout: `o_raw = total_no_field`, new total
   `+= n_field_units`.
2. `ms_abun_nuts_eval`: reconstruct `f = Linv %*% raw` once; in the per-species
   site loop add `f[field_map[site]]` to `eta_lambda`; accumulate a per-unit
   field-score vector `g_f[u] += grad_eta_lambda` across species (reduced
   serially in unit order for byte-exactness); after the species reduction,
   `g_raw = Linv' %*% g_f - raw` and add `-0.5||raw||^2` to lp. The chol / z / mu
   gradients are untouched.
3. R side (`.tobs_fit_ms_abun_nuts_spatial`, R/ms_abun_nuts.R): warm the field
   precision (tau[, rho]) + community means / covariances from the #12 sfMsNMix
   nested-Laplace fit (`.tobs_fit_ms_nmix_spatial`); build Linv = backsolve of
   chol(tau Q(rho) + ridge); pack `raw0 = chol(...) %*% f_warm` (so f = Linv raw
   returns the warm field); inv-metric 1 on the raw block. Reuse the #14 chains /
   rhat-ess / reconstruction wholesale, adding the field posterior mean
   `f = Linv %*% colMeans(raw)` to `fit$spatial_field`.
4. Dispatcher (`.dispatch_ms_abun`): a spatial term on the abundance formula
   under method="nuts" routes here (proper-CAR only; icar/bym2 -> nested_laplace
   pointer, matching the single-species abun gate).

## Why this route over the alternatives

- vs the tulpa#67 spatial-FACTOR route (reduced-rank loadings, as in
  ms_occu_cover spatial): the factor model is a DIFFERENT model (per-species
  loadings on a shared low-rank field, a JSDM association structure), not a
  single shared field. #73 asks for the shared-field sfMsNMix analogue, so the
  fixed-hyper field is the faithful sampler. The factor route stays the path for
  community ASSOCIATIONS (it is already shipped for occu_cover).
- vs sampling the field hyperparameter: the fixed-hyper choice is the tulpa#87
  result -- it removes the worst funnel at negligible cost to calibration because
  the field shape is well identified by the data at the family's sample sizes.

## Status

Design + core C++ field block implemented in this branch (the eval gains the
shared-field offset + the Linv-projected raw gradient, byte-checked against an R
oracle); the R front-door warm-starts from the #12 sfMsNMix nested-Laplace fit.
proper-CAR first; icar/bym2 NUTS+spatial gated to nested_laplace.
