# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Resolve a response on the top formula LHS. When `formula` is two-sided, the
# LHS is the response: it is the clean alternative to a separate `y =` for a
# single-vector-response family (`family$response == "vector"`, the cover
# hurdle). The LHS expression is evaluated against `data` first, then the
# calling environment, so it may be a bare column (`cover.flat`) or an
# expression (`log(cover + 1)`). The formula is stripped to one-sided so every
# downstream consumer (dispatchers, family encoders, the structured-term parser,
# which all assume a one-sided process formula with the response in `y`) is
# unchanged; the RHS -- spatial() / bars / other terms -- still flows through
# the existing parser untouched.
#
# Returns list(formula =, y =). A one-sided formula passes through unchanged
# (current interface). A two-sided formula errors for a matrix-response family
# (its response is a matrix / array, not a formula LHS) or when `y =` is also
# supplied (the response would be given twice).
.tobs_resolve_response_lhs <- function(formula, y, family,
                                       data, env = NULL) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  # The formula carries the user's calling environment; resolve LHS symbols not
  # found in `data` against it (matching how model.matrix resolves the RHS).
  if (is.null(env)) env <- environment(formula) %||% parent.frame()
  # One-sided formula (length 2: `~ rhs`) is the current interface, unchanged.
  if (length(formula) < 3L) {
    return(list(formula = formula, y = y))
  }

  # Two-sided formula: the LHS is a response. Only single-vector-response
  # families can take it there.
  if (!identical(family$response %||% "matrix", "vector")) {
    stop(
      sprintf(
        "%s()'s response is a matrix / array supplied via `y =`, not on the ",
        family$name),
      "formula left-hand side. Use a one-sided `formula = ~ predictors` and ",
      "pass the response as `y =` (see `?tobs`, the `y` argument).",
      call. = FALSE
    )
  }

  if (!is.null(y)) {
    stop(
      "the response is given twice: once on the formula left-hand side ",
      "(`", deparse(formula[[2L]]), " ~ ...`) and once via `y =`. Supply it ",
      "one way only -- either move it to the LHS and drop `y =`, or keep ",
      "`y =` and make `formula` one-sided.",
      call. = FALSE
    )
  }

  lhs <- formula[[2L]]
  data_env <- if (!missing(data) && !is.null(data) &&
                  (is.data.frame(data) || is.list(data))) {
    list2env(as.list(data), parent = env)
  } else env
  y <- tryCatch(
    eval(lhs, envir = data_env),
    error = function(e) stop(sprintf(
      "Could not evaluate the response `%s` on the formula left-hand side: %s",
      deparse(lhs), conditionMessage(e)), call. = FALSE)
  )

  # Strip to one-sided so downstream code sees the unchanged interface. Build
  # `~ rhs` from the RHS call object directly (not deparse-then-reparse), so a
  # long multi-line RHS -- e.g. a wide spatial() term -- survives intact.
  rhs_formula <- stats::as.formula(call("~", formula[[3L]]),
                                   env = environment(formula))
  list(formula = rhs_formula, y = y)
}

# Public `method` names are sugar over the orthogonal internal triple
# (engine, approx, correction). Each method names one fully-specified route;
# invalid cross-products (e.g. NUTS with an SLA marginal) simply have no name.
.tobs_method_table <- list(
  laplace            = list(engine = "laplace",        approx = "gaussian_laplace",   correction = "none"),
  laplace_sla        = list(engine = "laplace",        approx = "simplified_laplace", correction = "none"),
  laplace_gibbs      = list(engine = "laplace",        approx = "gaussian_laplace",   correction = "gibbs"),
  laplace_mi         = list(engine = "laplace",        approx = "gaussian_laplace",   correction = "mi"),
  pg_gibbs           = list(engine = "pg_gibbs",       approx = "gaussian_laplace",   correction = "none"),
  nested_laplace     = list(engine = "nested_laplace", approx = "gaussian_laplace",   correction = "none"),
  nested_laplace_sla = list(engine = "nested_laplace", approx = "simplified_laplace", correction = "none"),
  nuts               = list(engine = "nuts",           approx = "gaussian_laplace",   correction = "none")
)

# Resolve a public method name to the internal (engine, approx, correction)
# triple. `"auto"` maps the family's default engine to its base method.
.tobs_resolve_method <- function(method, family) {
  if (identical(method, "auto")) {
    method <- switch(
      family$default_engine,
      laplace        = "laplace",
      nested_laplace = "nested_laplace",
      nuts           = "nuts",
      pg_gibbs       = "pg_gibbs",
      stop(sprintf("Family '%s' has an unknown default_engine '%s'.",
                   family$name, family$default_engine), call. = FALSE)
    )
  }
  route <- .tobs_method_table[[method]]
  if (is.null(route)) {
    stop(sprintf("Unknown method '%s'. See `?tobs` for the route list.",
                 method), call. = FALSE)
  }
  route$method <- method   # resolved public name (auto -> concrete) for messages
  route
}

# ---------------------------------------------------------------------------
# Per-family backend coverage
#
# Single source of truth for which `method` each working family supports.
# `tobs()` validates the resolved (auto -> concrete) method against this set
# and errors with a pointer to the supported methods, rather than silently
# downgrading the engine (the old `.map_engine()` nested_laplace -> single-
# Laplace fall-back, which also mislabelled `fit$method`) or scattering the
# rejection across each family's dispatcher (the cover hurdle's bespoke
# `stop()`s). Gating mirrors the (engine, approx, correction) architecture: a
# method is listed iff its engine has a real execution path for the family.
#
#   * nested_laplace -- the nested-Laplace engine assembles a multi-block latent
#     prior (spatial / temporal / iid) and routes the state ("occ") M-step block
#     through `tulpa::tulpa_nested_laplace()`. Wired for single-season,
#     integrated, and dynamic occupancy (`.tobs_em_nested_laplace()`,
#     which shares the per-model-type callbacks with the Laplace path) and for
#     the cover hurdle's joint path (`tulpa_nested_laplace_joint()`). It also
#     supports INLA-style NA-response prediction (held-out sites), so the latent
#     field interpolates occupancy at unsurveyed sites.
#   * nested_laplace_sla -- the skew correction on the nested path is wired for
#     single-season occupancy and the cover hurdle only.
#   * laplace / laplace_sla / laplace_gibbs / laplace_mi -- run on tulpa's
#     EM+Laplace engine, which has callbacks for every occupancy family. The
#     cover hurdle is fit by a separate two-Laplace dispatcher with no EM
#     correction engine, so it takes no laplace_gibbs / laplace_mi.
#   * nuts -- sampled either by the unified C++ entry `cpp_occu_fit` (the
#     occupancy families) or by an in-tree FullGradFn over the family's own
#     closed-form marginal. The per-family comments below name the file each
#     path lives in.
#
# The list below is the roster, and every family constructor the package exports
# has an entry. A family object with no entry is a no-op for the validator and
# is rejected at dispatch instead.
.tobs_family_methods <- list(
  occu     = c("laplace", "laplace_sla", "laplace_gibbs", "laplace_mi",
               "pg_gibbs", "nested_laplace", "nested_laplace_sla", "nuts"),
  dyn_occu = c("laplace", "laplace_sla", "laplace_gibbs", "laplace_mi",
               "nested_laplace", "nuts"),
  # ms_occu: community single-season occupancy via the shared community
  # Laplace-EM (R/community_em.R) -- per-species occupancy / detection
  # coefficient RE with independent per-arm Gaussian community covariances. The
  # latent state marginalizes in closed form. nuts: the non-spatial community
  # sampler over the closed-form occupancy two-state per-(species, site) marginal
  # via the in-tree C++ FullGradFn (R/ms_occu_nuts.R, src/ms_occu_nuts.cpp) --
  # samples the community means, per-species deviations, AND the two independent
  # per-arm community covariances jointly, non-centered, warm-started at the
  # Laplace-EM mode. nested_laplace: a shared areal field (icar/bym2/car_proper)
  # on the occupancy arm via the in-tree community-spatial nested Laplace-EM
  # (R/ms_occu_spatial.R, src/ms_occu_spatial.cpp) -- the occupancy analogue of
  # sfMsNMix.
  ms_occu  = c("laplace", "nuts", "pg_gibbs", "nested_laplace"),
  int_occu = c("laplace", "laplace_sla", "laplace_gibbs", "laplace_mi",
               "nested_laplace", "nuts"),
  # jsdm: the joint species distribution model observes presence/absence directly
  # (no detection process), with per-species coefficients under a Gaussian
  # community covariance -- the spOccupancy lfJSDM / sfJSDM model class, i.e. the
  # community GLMM of ms_count() with a logit link, so it shares that binder,
  # community Laplace-EM, latent driver and NUTS target. nested_laplace: a shared
  # areal field (icar/car_proper/bym2), optionally with latent() factors
  # (sfJSDM), via block coordinate ascent. laplace: the non-spatial community EM,
  # or latent() factors alone (lfJSDM). nuts: the exact joint community posterior
  # over the Bernoulli response. The single-block correction routes (laplace_sla
  # / laplace_gibbs / laplace_mi) belonged to the former shared-FE +
  # scalar-species-intercept model and do not apply to the community EM.
  jsdm     = c("laplace", "nuts", "pg_gibbs",
               "nested_laplace"),
  # count: GLMM on the observed count / continuous response directly (no
  # detection, no latent state) -- the relative-abundance model of spAbundance
  # (abund). Non-spatial Laplace only for the first ship: a single tulpa GLMM
  # block (Poisson / neg_binomial_2 / gaussian). The negbin size / gaussian
  # residual variance is estimated by an outer dispersion loop in .dispatch_count
  # (tulpa_laplace takes a fixed phi). nested_laplace: a plain areal field (icar
  # / bym2 / car_proper) on the abundance formula (the spAbund analogue) -- the
  # field is a latent GMRF prior on the count block, integrated over its
  # hyperparameters via the shared nested-Laplace EM machinery. Community
  # (msAbund) / NUTS are the documented follow-ups.
  count    = c("laplace", "nested_laplace"),
  # ms_count: community relative-abundance GLMM (msAbund) via the shared community
  # Laplace-EM (R/community_em.R) -- per-species coefficient RE with a Gaussian
  # community covariance; Poisson / negbin (per-species dispersion RE) / gaussian.
  # nested_laplace: a shared areal field icar() on the abundance formula (the
  # sfMsAbund analogue, Poisson) via block coordinate ascent -- the community EM
  # (field as a per-site offset) alternated with a self-contained Poisson-ICAR
  # field update (R/ms_count_spatial.R), no C++. nuts: the exact joint community
  # count posterior via the in-tree C++ FullGradFn (R/ms_count_nuts.R,
  # src/ms_count_nuts.cpp), warm-started at the Laplace-EM mode; Poisson.
  ms_count = c("laplace", "nested_laplace", "nuts", "pg_gibbs"),
  # abun: non-spatial N-mixture (laplace; Poisson or negbin) + areal-spatial offset
  # (nested_laplace: icar / bym2 / car_proper on the abundance arm). tulpa's
  # spatial fitters return the grid-integrated coefficient covariance, so the
  # spatial SEs are calibrated (law-of-total-covariance over the hyperparameter
  # grid). nuts: the non-spatial sampler over the closed-form marginal via the
  # in-tree C++ FullGradFn (R/abun_nuts.R, src/abun_nuts.cpp); Poisson or negbin
  # (log_r sampled), warm-started at the Laplace mode. NUTS + areal: a fixed-hyper
  # non-centered field on the abundance arm -- car_proper via the square inverse
  # Cholesky, icar / bym2 via the sum-to-zero reparameterisation (the intrinsic
  # field's null-space direction dropped). A single intercept RE samples
  # non-spatially; spatial XOR RE.
  abun     = c("laplace", "nested_laplace", "nuts"),
  # royle_nichols: Royle-Nichols occupancy (abundance-induced detection
  # heterogeneity). Latent N ~ Poisson marginalised in closed form; the exact
  # marginal is maximised (optim BFGS) with an observed-information vcov. Site-
  # level detection; non-spatial laplace only for the first ship (visit-level
  # detection, areal fields, and NUTS are the documented follow-ups, #116).
  royle_nichols = c("laplace"),
  # occu_ttd: time-to-detection occupancy (exponential TTD, unmarked occuTTD).
  # Two-state occupancy marginal with a continuous censored-exponential emission,
  # maximised in closed form (optim BFGS) with an observed-information vcov.
  # Site-level rate; non-spatial laplace only for the first ship (Weibull shape,
  # visit-varying rate, areal fields, NUTS are the documented follow-ups, #116).
  occu_ttd = c("laplace"),
  # occu_multi: multi-species co-occurrence occupancy (Rota 2016, unmarked
  # occuMulti). Joint 2^S-state log-linear occupancy with first + second order
  # natural parameters; the exact marginal is enumerated and maximised (optim
  # BFGS) with an observed-information vcov. Shared covariate design, site-level
  # detection, non-spatial laplace only for the first ship (per-parameter
  # formulas, higher-order terms, visit-level detection are follow-ups, #116).
  occu_multi = c("laplace"),
  # double_observer: double-observer abundance (unmarked multinomPois pi). Poisson-
  # multinomial thinning makes the observable cell counts independent Poissons, so
  # the marginal is closed form (no latent-N sum); maximised by optim BFGS with an
  # observed-information vcov. type = "independent" (3 cells) or "dependent"
  # (removal-style, 2 cells, role-swapping `primary` for identifiability, #116).
  # Site-level detection, non-spatial laplace only.
  double_observer = c("laplace"),
  # gdistremoval: joint distance + removal sampling (unmarked gdistremoval,
  # Amundson et al. 2014). Single-season; the detected birds are cross-classified
  # by distance band and removal period. Binomial thinning of a Poisson N is
  # closed under thinning, so the total-detected count is Poisson and the band /
  # period allocations are two conditional multinomials -- a closed-form marginal
  # (no latent-N sum), maximised by optim BFGS with an observed-information vcov.
  # Site-level arms, half-normal key, non-spatial laplace only for the first ship.
  gdistremoval = c("laplace"),
  # distsamp_open: open-population distance sampling (unmarked distsampOpen). A
  # Dail-Madsen open N-mixture (dyn_abun) with a distance-bin multinomial emission
  # per primary period. The band allocation is conditional on the period total, so
  # it factors out of the abundance HMM: the marginal reuses the dyn_abun forward
  # kernel (eta_p = logit(pdist)) + the per-period band multinomials, maximised by
  # optim BFGS with an observed-information vcov. Site-level arms, half-normal key,
  # constant dynamics, non-spatial laplace only for the first ship.
  distsamp_open = c("laplace"),
  # dyn_int_occu: multi-season integrated occupancy (spOccupancy tIntPGOcc). A
  # dynamic-occupancy HMM whose per-season emission pools multiple detection
  # sources; the latent state integrates out by the two-state forward recursion,
  # maximised with analytic (forward-backward Fisher-identity) gradients + an
  # observed-information vcov. A shared areal icar() field on the first-season
  # occupancy formula fits stIntPGOcc under nested_laplace via the shared
  # areal-BFGS driver (the field gradient is the psi1 score w1 - psi1, #122).
  # v1 = full site / season overlap, constant transitions, site-level detection
  # (partial overlap, season-varying rates, bym2 / car_proper, NUTS are
  # documented follow-ups, #122).
  dyn_int_occu = c("laplace", "nested_laplace"),
  # t_occu: multi-season occupancy with an AR1 year random effect on the state
  # (spOccupancy tPGOcc). NOT colext -- a per-(site, season) Bernoulli GLMM with a
  # shared AR1 year effect; the seasons factorise given the year effects, so the
  # Polya-Gamma Gibbs sampler is exact (the engine spOccupancy uses). pg_gibbs only.
  t_occu = c("pg_gibbs"),
  # ms_abun: community / multispecies N-mixture via the in-tree C++ Laplace-EM
  # (per-species coefficient RE with Gaussian community covariances). A shared areal
  # field (icar / bym2 / car_proper) on the abundance arm fits under nested_laplace;
  # Poisson or grid-integrated negbin size. nuts: the non-spatial community sampler
  # over the closed-form per-(species, site) marginal via the in-tree C++ FullGradFn
  # (R/ms_abun_nuts.R, src/ms_abun_nuts.cpp) -- samples the community means,
  # per-species deviations, AND community covariances jointly; Poisson or
  # per-species negbin (log_r_s sampled), warm-started at the Laplace-EM mode. nuts
  # + a shared areal field (car_proper, Poisson) joins a fixed-hyper non-centered
  # proper-CAR field on the abundance arm; icar/bym2 + temporal / RE NUTS not yet
  # wired.
  ms_abun  = c("laplace", "nested_laplace", "nuts"),
  # removal: sequential-depletion removal sampling. Non-spatial closed-form marginal
  # Laplace (Poisson or negbin; the depleting-binomial product summed over latent
  # N), grouped-RE AGHQ Laplace, the in-tree C++ FullGradFn NUTS over the same
  # marginal, and an areal icar()/car_proper() field on the abundance arm via
  # nested_laplace (the shared count-marginal spatial driver). A temporal()
  # AR1/RW1/RW2/iid field composes WITH the areal field on the abundance arm under
  # nested_laplace via the shared areal-BFGS driver (a second latent block, both
  # grid-integrated), OR rides the family NUTS field block on its own
  # (temporal-only, no simultaneous areal field). spde not yet wired (R/removal.R,
  # R/removal_spatial.R).
  removal  = c("laplace", "nested_laplace", "nuts"),
  # distance: binned distance sampling (half-normal / hazard-rate key, line / point
  # transect). Non-spatial closed-form marginal Laplace (Poisson or negbin),
  # grouped-RE AGHQ Laplace (abundance arm), the in-tree C++ FullGradFn NUTS, and an
  # areal icar()/car_proper()/bym2() field on the abundance arm via nested_laplace
  # (the per-site var_N rank-1 cross-arm from distance_kernel.h). The hazard-rate
  # scalar log-shape is threaded into the areal-BFGS fixed block, so a hazard-key
  # areal fit recovers both the abundance field and the shape. A temporal()
  # AR1/RW1/RW2/iid field composes WITH the areal field on the abundance arm under
  # nested_laplace via the shared areal-BFGS driver (a second latent block), OR
  # rides the family NUTS field block on its own (temporal-only). The hazard-rate
  # key's global log-shape rides alongside the NUTS field block (areal or temporal),
  # so a hazard-key NUTS+field fit recovers both the abundance field and the shape.
  # Grouped-RE hazard and hazard-key DETECTION-arm areal not yet wired
  # (R/distance.R, R/distance_spatial.R, src/distance_*.cpp).
  distance = c("laplace", "nested_laplace", "nuts"),
  # ms_distance: community binned distance sampling (the spAbundance msDS
  # analogue) -- per-species distance sampling with Gaussian community
  # hyperpriors on the abundance / detection-scale coefficients, over the shared
  # community Laplace-EM (R/ms_distance.R). The latent N still integrates out in
  # closed form per species-site, so every fit is driven by the existing distance
  # kernel (cpp_distance_site_sweep); no new C++. laplace: the plain community
  # fit, or latent() factors alone (lfMsDS). nested_laplace: a shared field
  # (icar/car_proper/bym2/spde) with or without factors (sfMsDS), via the same
  # block coordinate ascent as every other community family
  # (R/community_latent.R). Poisson only (the negbin size is not yet a
  # per-species RE); NUTS not wired.
  ms_distance = c("laplace", "nested_laplace"),
  # fp_occu: multistate false-positive occupancy (Miller et al. 2011). Latent
  # occupancy z summed out in closed form (two states); four site-level logit arms
  # (psi, true detection p11, false-positive p10, certain-classification b).
  # Non-spatial analytic-gradient BFGS over the exact marginal with an
  # observed-information vcov (laplace), grouped-RE AGHQ Laplace (psi or p11 arm),
  # the in-tree C++ FullGradFn NUTS, and an areal icar()/car_proper() field on the
  # occupancy arm via nested_laplace (BFGS over the marginal + CAR prior, FD-Hessian
  # observed info). A temporal() AR1/RW1/RW2/iid field composes WITH the areal field
  # on the psi arm under nested_laplace via the shared areal-BFGS driver (a second
  # latent block), OR rides the family NUTS field block on its own (temporal-only, no
  # simultaneous areal field). (R/fp_occu.R, R/fp_occu_spatial.R, src/fp_occu_*.cpp).
  fp_occu  = c("laplace", "nested_laplace", "nuts"),
  # dyn_abun: Dail-Madsen open-population N-mixture (Poisson initial abundance,
  # binomial survival, Poisson recruitment, binomial detection). The latent
  # abundance sequence is summed out by an exact HMM forward recursion (not closed
  # form); analytic gradients by forward-mode differentiation. Non-spatial
  # analytic-gradient BFGS over the forward marginal with an observed-information
  # vcov (laplace), grouped-RE AGHQ Laplace (initial-abundance arm), the in-tree
  # C++ FullGradFn NUTS, and an areal icar()/car_proper() field on the initial-
  # abundance arm via nested_laplace (BFGS over the forward marginal + CAR prior,
  # FD-Hessian observed info). A temporal() AR1/RW1/RW2/iid field composes WITH the
  # areal field on the initial-abundance arm under nested_laplace via the shared
  # areal-BFGS driver (a second latent block). Season-varying survival /
  # recruitment: a covariate on omega_formula / gamma_formula carried as a [n_sites
  # x (T-1)] matrix column gives interval- indexed vital rates in the forward
  # kernel, on every backend (laplace, NUTS, nested_laplace areal on the
  # initial-abundance arm; R/dyn_abun.R, src/dyn_abun_kernel.h). NUTS+temporal not
  # yet wired (R/dyn_abun.R, R/dyn_abun_spatial.R, src/dyn_abun_*.cpp).
  dyn_abun = c("laplace", "nested_laplace", "nuts"),
  # cover: standalone vegetation-cover hurdle (presence Bernoulli + beta /
  # lognormal positive arm). Non-spatial Laplace via two independent
  # tulpa_laplace() calls (laplace / laplace_sla); a shared areal field across
  # the arms fits under nested_laplace / nested_laplace_sla. nuts: the
  # non-spatial sampler over the exact two-arm coefficient marginal
  # c(beta_presence, beta_positive, log_disp) via the in-tree C++ FullGradFn
  # (R/cover_nuts.R, src/cover_nuts.cpp), warm-started at the Laplace mode --
  # calibrated (non-Gaussian) intervals and a per-draw pointwise likelihood for
  # WAIC / LOO, beta or lognormal cover. Structured-term NUTS is gated (the
  # shared field is grid-integrated under nested_laplace).
  cover    = c("laplace", "laplace_sla", "nested_laplace", "nested_laplace_sla",
               "nuts"),
  # occu_cover: non-spatial Laplace via direct optim on the exact two-state
  # marginal; nested-Laplace adds a cell-level ICAR field shared across the psi
  # and cover arms with scaling alpha, in the INLA `copy =` idiom. A bym2() term
  # is read as ICAR (rho fixed to 1).
  # nuts: the non-spatial sampler over the exact two-state coefficient marginal
  # via the in-tree C++ FullGradFn (R/occu_cover_nuts.R, src/occu_cover_nuts.cpp),
  # warm-started at the Laplace mode -- calibrated (non-Gaussian) intervals and a
  # per-draw pointwise likelihood for WAIC / LOO, beta or lognormal cover. A
  # spatial occu_cover NUTS path is not yet wired (the shared coupled field is
  # grid-integrated under nested_laplace; the spatial-factor community sampler
  # ms_occu_cover() + icar() samples a shared field).
  occu_cover = c("laplace", "nested_laplace", "nuts"),
  # occu_multiscale_cover: three-level cell / plot / visit occupancy + cover.
  # "nested_laplace" carries the shared areal field (the four-arm cell-coupling
  # spec); "laplace" / "nuts" are the non-spatial path (iid cells, no field) --
  # the exact three-level marginal optimised directly (Laplace) or sampled
  # (NUTS, the exact coefficient posterior + calibrated WAIC / LOO). Cells are
  # declared the same way on every path, via icar(group_var = "<cell>") (the
  # graph is ignored under "laplace" / "nuts"). Both marginalize z (cells) and
  # a (plots) in closed form.
  occu_multiscale_cover = c("laplace", "nested_laplace", "nuts"),
  # ms_occu_cover: community joint occupancy-detection + cover. Per-species
  # coefficient RE with Gaussian community covariances across the psi / p / pos
  # arms; the latent presence z integrates out in closed form (the occu_cover
  # marginal) and the per-species deviations are integrated by a Laplace-EM.
  # Non-spatial only -- the community analogue of the joint-coupled spatial
  # engine (per-species RE layered on the shared coupled field) needs upstream
  # tulpa support, so nested_laplace is not offered. nuts: the reduced-rank
  # spatial-factor path (a shared icar/car/bym2 field with per-species loadings)
  # samples the exact joint posterior via tulpa's NUTS + the in-tree C++
  # FullGradFn. Non-spatial ms_occu_cover has no NUTS path (gated in the
  # dispatcher).
  ms_occu_cover = c("laplace", "nuts"),
  # ms_dyn_occu / ms_int_occu: community dynamic / integrated occupancy. Per-
  # species coefficient RE with per-arm Gaussian community covariances, fit by the
  # shared community Laplace-EM (R/community_em.R). The latent occupancy path (HMM
  # forward for dynamic, two-state mixture for integrated) marginalizes in closed
  # form. nuts (ms_dyn_occu): the non-spatial community sampler over the exact
  # HMM-forward per-species marginal via the in-tree C++ FullGradFn
  # (R/ms_dyn_occu_nuts.R, src/ms_dyn_occu_nuts.cpp) -- samples the community
  # means, per-species first-season / detection deviations, the two independent
  # per-arm community covariances, AND the shared colonisation / extinction globals
  # jointly, non-centered, warm-started at the Laplace-EM mode. Non-spatial only
  # (NUTS + areal field -> nested_laplace).
  ms_dyn_occu = c("laplace", "pg_gibbs", "nested_laplace", "nuts"),
  ms_int_occu = c("laplace", "pg_gibbs", "nuts"),
  # occu_categorical: presence + nominal K-class hurdle. A Bernoulli presence arm
  # and a baseline-category multinomial logit on the class given present (the
  # FD-validated tulpa multinomial kernel; the non-spatial fit is the vectorised R
  # Newton over the same closed forms). Non-spatial Laplace only for the first
  # ship; the native multi-process LikelihoodSpec path (spatial fields / NUTS) and
  # the latent-class confusion variant (the K-class generalisation of fp_occu) are
  # the documented follow-ups.
  occu_categorical = c("laplace")
)

# Validate a resolved public method name against the family's supported set.
# `method` is the concrete name ("auto" already resolved upstream). No-op for a
# family with no entry: every exported family constructor has one, so a family
# object that reaches here without an entry is a hand-built `obs_family()`, and
# `tobs()` rejects it at dispatch.
.tobs_validate_family_method <- function(method, family) {
  supported <- .tobs_family_methods[[family$name]]
  if (is.null(supported) || method %in% supported) return(invisible(NULL))
  stop(
    sprintf(
      "method = \"%s\" is not available for %s() (%s). Supported: %s.",
      method, family$name, family$class_long,
      paste0("\"", supported, "\"", collapse = ", ")
    ),
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# Control-option validation
#
# `control` is splatted as named args onto `.tobs_fit_model()`, whose formals
# are the union of every knob across all methods (plus a trailing `...`). That
# means a control that does not apply to the chosen method is silently ignored
# (e.g. `n.chains` under `"laplace"`) and a typo (`niter`) vanishes into `...`.
# We validate names up front against a per-route allowlist so misapplied or
# misspelled controls error instead.
#
# Single source of truth: control names are grouped by capability, and each
# engine/correction route admits a set of groups. `sigma.beta` is shared by the
# Laplace and NUTS paths; `seed` and `n.seeds` by the stochastic-correction and
# NUTS paths (the deterministic Laplace routes reject `n.seeds` here, since
# seed-variant fits would be identical -- see the ensemble branch in tobs()).
# ---------------------------------------------------------------------------
.tobs_control_groups <- list(
  laplace_em = c("max.iter", "tol", "damping", "sigma.beta",
                 "re.aghq", "n.quad", "n.quad.scalar", "re.lkj", "optimizer",
                 "omega.sigma.prior", "logr.sigma.prior",
                 "hessian", "inner.solver", "integration"),
  # The community latent routes -- a shared areal field and/or latent() factors
  # on a community family -- fit by block coordinate ascent between the
  # community Laplace-EM and the field / factor updates (R/community_latent.R).
  # `max.outer` caps that outer alternation. Admitted on both Laplace routes: a
  # factor-only model is method = "laplace", a shared field "nested_laplace".
  # `factor.starts` sets how many candidate starting directions the first factor
  # pass selects over; each costs a full loading-EM run against the family's own
  # oracle, so it is the dominant cost on families whose oracle marginalises a
  # latent state. Opted into per family via `obs_family(control_groups=)`, NOT
  # admitted route-wide: only the community families fit by
  # `.tobs_community_latent_ascent()` have an outer alternation to cap, and a
  # route-wide allowance let `max.outer` be passed to every Laplace family and
  # dropped. Still route-gated, so it is rejected under "nuts" as a wrong-method
  # control rather than an unknown one.
  block_coordinate = "max.outer",
  # Split from `block_coordinate` because the two are not co-extensive:
  # `ms_dyn_occu` reaches the driver with `latent = NULL`, so it has an outer
  # alternation but no candidate set of starting directions to widen.
  block_coordinate_factor = "factor.starts",
  # Outer-grid knobs for the standalone occu() varying-coefficient (SVC) bar,
  # which reroutes from the EM fixed-point path onto the joint direct-grid engine
  # under method = "nested_laplace". They are no-ops on the EM path (a plain
  # intercept field, temporal / re structure), which the SVC reroute predicate
  # skips; admitted on the nested_laplace route so a SVC fit can tune the grid /
  # threading without a separate method name.
  nested_laplace_joint = c("sigma.grid",
                           # The other three outer-grid axes a latent block
                           # can carry: the mixing / correlation parameter,
                           # the precision, and an SPDE range.
                           "rho.grid", "tau.grid", "range.grid",
                           "n.threads", "n.threads.outer",
                           "adaptive.grid", "adaptive.grid.edge.thresh",
                           "adaptive.grid.max.passes",
                           # Tuning for `integration = "grid_adaptive"` (the
                           # subset-lattice builder), distinct from the three
                           # post-integration refinement knobs above.
                           "adaptive.grid.cutoff", "adaptive.grid.stride",
                           "adaptive.grid.max.frac", "adaptive.grid.min.cells",
                           "var.of.means.consistency",
                           "var.of.means.min.ess", "diagnose.k", "diagnose.draws",
                           "k.samples", "k.bootstrap", "k.tail.points", "k.conf.bands",
                           "force.sparse", "inner.refresh", "checkpoint",
                           # Regularizing hyperpriors on the outer grid axes,
                           # forwarded to the joint driver's prior_sigma / _alpha
                           # / _phi (e.g. a PC prior on the spatial field SD).
                           "prior.sigma", "prior.alpha", "prior.phi"),
  correction = c("n.gibbs", "n.imputations", "seed", "n.seeds"),
  sampler    = c("n.iter", "n.warmup", "n.thin", "n.chains", "n.threads",
                 # OpenMP threads inside ONE gradient evaluation of the
                 # community NUTS targets, whose per-species loop is
                 # parallel. Distinct from `n.threads`, which spreads whole
                 # chains. 0 leaves the count to OpenMP.
                 "n.threads.grad",
                 "adapt.delta", "max.treedepth", "seed", "sigma.beta",
                 "n.seeds",
                 # Community-mean prior SD on the log-dispersion mu_log_r, for
                 # the negative-binomial NUTS paths that carry one
                 # (ms_abun(), ms_count(), jsdm()). Ignored by a family or
                 # mixture with no log-dispersion arm.
                 "sigma.logr",
                 # ms_occu_cover() NUTS per-species dispersion RE (#115 B7): opt
                 # into a fourth 1-D community arm on the cover log-dispersion.
                 "dispersion.re", "sigma.ld.init"),
  universal  = c("verbose",
                 "progress", "progress.every", "progress.throttle",
                 "progress.file")
)

# Capability groups admitted by a resolved (engine, correction) route.
.tobs_control_allow <- function(engine, correction) {
  switch(
    engine,
    laplace        = c("laplace_em", if (correction != "none") "correction"),
    nested_laplace = c("laplace_em", "nested_laplace_joint"),
    nuts           = "sampler",
    pg_gibbs       = "sampler",
    character(0)
  )
}

# Family-opted capability groups, gated by route. A family declares these via
# `obs_family(control_groups=)`; they are admitted only on the engines that host
# them, so a Laplace-only group stays rejected under a sampler route.
.tobs_family_group_hosts <- list(
  block_coordinate        = c("laplace", "nested_laplace"),
  block_coordinate_factor = c("laplace", "nested_laplace"))

.tobs_family_groups <- function(family, engine) {
  groups <- family$control_groups %||% character(0)
  if (!length(groups)) return(character(0))
  keep <- vapply(groups, function(g) {
    hosts <- .tobs_family_group_hosts[[g]]
    is.null(hosts) || engine %in% hosts
  }, logical(1))
  groups[keep]
}

# Public method names that accept a given control key (for "wrong method"
# hints). Derived from the route table + allowlist so it stays in sync.
.tobs_methods_for_control <- function(key) {
  in_group <- vapply(.tobs_control_groups, function(g) key %in% g, logical(1))
  groups   <- names(.tobs_control_groups)[in_group]
  if (!length(groups)) return(character(0))
  methods <- names(.tobs_method_table)
  keep <- vapply(methods, function(m) {
    r <- .tobs_method_table[[m]]
    # A family-opted group belongs to the methods whose engine hosts it, even
    # though no route admits it unconditionally -- otherwise a key in such a group
    # would report as applying to no method at all.
    hosted <- names(Filter(function(h) r$engine %in% h, .tobs_family_group_hosts))
    any(c(.tobs_control_allow(r$engine, r$correction), hosted) %in% groups)
  }, logical(1))
  methods[keep]
}

# Validate `control` names against the resolved route. Errors on (a) a known
# control that the chosen method does not use, or (b) an unrecognized name
# (with a fuzzy "did you mean" suggestion). Collects all offenders into one
# message. No-op for a valid (or empty) control list.
.tobs_validate_control <- function(control, route, family = NULL) {
  if (length(control) == 0L) return(invisible(NULL))
  nms <- names(control)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("`control` must be a fully named list, e.g. ",
         "control = list(n.iter = 4000).", call. = FALSE)
  }

  # Family-specific dispatchers (e.g. the cover hurdle) declare extra control
  # names via family$control_keys; admit those alongside the engine controls.
  family_keys <- if (!is.null(family)) family$control_keys %||% character(0)
                 else character(0)

  route_groups  <- .tobs_control_allow(route$engine, route$correction)
  family_groups <- .tobs_family_groups(family, route$engine)
  allowed_groups <- c("universal", route_groups, family_groups)
  allowed_keys <- unique(c(unlist(.tobs_control_groups[allowed_groups],
                                  use.names = FALSE),
                           family_keys))
  vocabulary   <- unique(c(unlist(.tobs_control_groups, use.names = FALSE),
                           family_keys))

  bad <- setdiff(nms, allowed_keys)
  if (!length(bad)) return(invisible(NULL))

  method <- route$method %||% route$engine
  # Keys whose group this ROUTE hosts but this FAMILY did not opt into. Pointing
  # at another method would be wrong -- no method makes them apply here.
  fam_only <- names(Filter(function(g) route$engine %in% (g %||% character(0)),
                           .tobs_family_group_hosts))
  fam_only_keys <- unlist(.tobs_control_groups[setdiff(fam_only, family_groups)],
                          use.names = FALSE)
  fam_label <- if (!is.null(family$name)) paste0(family$name, "()") else "this family"
  msgs <- vapply(bad, function(key) {
    if (key %in% fam_only_keys) {
      sprintf("  - '%s' is not used by %s.", key, fam_label)
    } else if (key %in% vocabulary) {
      uses <- .tobs_methods_for_control(key)
      sprintf("  - '%s' is not used by method = \"%s\"; it applies to %s.",
              key, method,
              paste0("method = \"", uses, "\"", collapse = " / "))
    } else {
      near <- agrep(key, vocabulary, value = TRUE, max.distance = 0.34)
      hint <- if (length(near))
        sprintf(" Did you mean %s?",
                paste0("'", near, "'", collapse = " / ")) else ""
      sprintf("  - '%s' is not a known control option.%s", key, hint)
    }
  }, character(1))

  stop("Invalid `control` option(s) for method = \"", method, "\":\n",
       paste(msgs, collapse = "\n"),
       "\nSee `?tobs` for the controls each method uses.", call. = FALSE)
}

.map_engine <- function(engine, family = NULL) {
  # Engine name translation between the tobs vocabulary and what the underlying
  # fitter currently understands.
  #
  # Which families the nested-Laplace engine serves is `.tobs_family_methods`'s
  # answer, read here rather than restated: the two lists had already drifted
  # apart in both directions, this one carrying two families that dispatch
  # directly and never reach it, and missing seven the registry does list.
  # `tobs()` runs `.tobs_validate_family_method()` before dispatch, so an
  # unsupported family arriving here is an internal mis-wire rather than a user
  # error to downgrade silently.
  if (engine == "nested_laplace") {
    if (isTRUE("nested_laplace" %in% .tobs_family_methods[[family]])) {
      return("nested_laplace")
    }
    stop(sprintf(
      "Internal error: nested_laplace reached .map_engine for family '%s'; the method registry should have rejected it.",
      family %||% "(unspecified)"), call. = FALSE)
  }
  switch(
    engine,
    laplace = "laplace",
    nuts    = "nuts",
    engine
  )
}

# Normalize the user's `visits` argument to the shape `.tobs_build_single`
# expects: a long data frame with `n_sites * max_visits` rows in site-major
# order (row `r` corresponds to site `(r-1) %/% max_visits + 1`, visit
# `(r-1) %% max_visits + 1`), plus the formula to apply to it.
#
# Accepts:
#   * NULL                                  -> NULL, no visit-level path
#   * named list of [n_sites, max_visits]   -> flatten to long DF, treat
#       matrices (the `tobs_data()` shape)     `detection` as visit-level
#                                              (intercept dropped, site-level
#                                              X_det is intercept-only)
#   * data.frame with N*J rows + a          -> existing dual-formula behavior:
#       "formula" attribute                    site-level `detection` against
#                                              `data`, attr formula against
#                                              `visits`
#   * data.frame with N*J rows, no formula  -> treat `detection` as visit-level
#       attribute                              (intercept dropped); site-level
#                                              X_det is intercept-only
#
# Returns a list with:
#   visits             — long data frame (or NULL)
#   det_visit_formula  — formula applied to visits (or NULL)
#   det_formula        — formula applied to site-level `data`
.normalize_visits <- function(visits, detection,
                              n_sites, max_visits) {
  if (is.null(visits)) {
    return(list(visits = NULL,
                det_visit_formula = NULL,
                det_formula = detection))
  }

  expected_rows <- n_sites * max_visits

  # Case 1: list of [n_sites, max_visits] matrices (tobs_data() output)
  if (is.list(visits) && !is.data.frame(visits)) {
    nms <- names(visits)
    if (is.null(nms) || any(!nzchar(nms))) {
      stop("`visits` (list of matrices) must be a named list; ",
           "names become the column names of the flattened frame.",
           call. = FALSE)
    }
    bad <- vapply(visits, function(m) {
      !is.matrix(m) || nrow(m) != n_sites || ncol(m) != max_visits
    }, logical(1))
    if (any(bad)) {
      stop(sprintf(
        "`visits` elements must be [%d x %d] matrices matching y; ",
        n_sites, max_visits),
        sprintf("element(s) %s have wrong shape.",
                paste(nms[bad], collapse = ", ")),
        call. = FALSE)
    }
    flat <- as.data.frame(
      lapply(visits, function(m) {
        col <- as.vector(t(m))
        if (isTRUE(attr(m, "tobs_factor"))) {
          factor(col, levels = attr(m, "tobs_levels"))
        } else {
          col
        }
      }),
      stringsAsFactors = FALSE
    )
    return(list(
      visits = flat,
      det_visit_formula = .drop_intercept(detection),
      det_formula = ~ 1
    ))
  }

  # Case 2 / 3: data frame
  if (is.data.frame(visits)) {
    if (nrow(visits) != expected_rows) {
      stop(sprintf(
        "`visits` (data frame) must have %d rows (n_sites * max_visits); ",
        expected_rows),
        sprintf("got %d.", nrow(visits)),
        call. = FALSE)
    }
    attached <- attr(visits, "formula")
    if (!is.null(attached)) {
      # Dual-formula power-user mode: detection stays site-level
      return(list(
        visits = visits,
        det_visit_formula = attached,
        det_formula = detection
      ))
    }
    return(list(
      visits = visits,
      det_visit_formula = .drop_intercept(detection),
      det_formula = ~ 1
    ))
  }

  stop("`visits` must be NULL, a named list of [n_sites x max_visits] ",
       "matrices, or a long data frame with n_sites * max_visits rows; ",
       "got ", paste(class(visits), collapse = "/"), ".",
       call. = FALSE)
}

# Compact (ragged) counterpart of .normalize_visits Case 1: the visit covariates
# arrive as a named list of length-V vectors (one entry per VALID visit, in the
# canonical order(site, visit)), not [n_sites x max_visits] matrices. Build the
# V-row visit frame directly -- no padding, no flatten -- and return the same
# (visits, det_visit_formula, det_formula) shape the dense path returns, so the
# downstream design build and arm assembly are identical. The visit-level design
# is then built with max_per_unit = NULL (the compact signal in
# .tobs_build_visit_X).
.normalize_visits_ragged <- function(visits, detection, n_visits_valid) {
  if (is.null(visits)) {
    return(list(visits = NULL, det_visit_formula = NULL, det_formula = detection))
  }
  if (!is.list(visits) || is.data.frame(visits)) {
    stop("compact `visits` must be a named list of length-V vectors ",
         "(tobs_data(compact = TRUE) output).", call. = FALSE)
  }
  nms <- names(visits)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("compact `visits` must be a named list; names become column names.",
         call. = FALSE)
  }
  bad <- vapply(visits, function(v) length(v) != n_visits_valid, logical(1))
  if (any(bad)) {
    stop(sprintf("compact `visits` element(s) %s have length != %d valid visits.",
                 paste(nms[bad], collapse = ", "), n_visits_valid), call. = FALSE)
  }
  flat <- as.data.frame(
    lapply(visits, function(v) {
      if (isTRUE(attr(v, "tobs_factor"))) factor(as.character(v),
                                                 levels = attr(v, "tobs_levels"))
      else as.numeric(v)
    }),
    stringsAsFactors = FALSE)
  names(flat) <- nms
  list(visits = flat,
       det_visit_formula = .drop_intercept(detection),
       det_formula = ~ 1)
}

# Drop the intercept term from a formula (returns `~ . - 1`-style update).
# Preserves the LHS if any (none of our detection formulas have one).
.drop_intercept <- function(f) {
  stats::update(f, ~ . - 1)
}

# ---------------------------------------------------------------------------
# Print method for fits
# ---------------------------------------------------------------------------

#' Print method for tobs_fit
#' @param x a `tobs_fit` object.
#' @param ... forwarded to underlying print methods.
#' @return `x`, invisibly.
#' @export
print.tobs_fit <- function(x, ...) {
  fam <- attr(x, "tobs_family")
  if (!is.null(fam)) {
    cat(sprintf("<tobs_fit: %s>\n", fam$class_long))
    cat(sprintf("  family         : %s (status: %s)\n", fam$name, fam$status))
    cat(sprintf("  default method : %s\n", fam$default_engine))
    cat("\n")
  }
  model <- x$model
  if (!is.null(model)) {
    if (model$model_type == "single" || model$model_type == "nmix" ||
        model$model_type == "removal") {
      lab <- if (identical(model$model_type, "removal")) "Passes" else "Max visits"
      cat(sprintf("  Sites: %d, %s: %d\n", model$n_sites, lab, model$max_visits))
    } else if (model$model_type == "dynamic") {
      cat(sprintf("  Sites: %d, Seasons: %d, Max visits: %d\n",
                  model$n_sites, model$n_seasons, model$max_visits))
    } else if (model$model_type == "ms_occu" ||
               model$model_type == "ms_nmix" ||
               model$model_type == "ms_occu_cover") {
      cat(sprintf("  Sites: %d, Species: %d\n", model$n_sites, model$n_species))
    } else if (model$model_type == "ms_dyn_occu") {
      cat(sprintf("  Sites: %d, Seasons: %d, Species: %d\n",
                  model$n_sites, model$n_seasons, model$n_species))
    } else if (model$model_type == "ms_int_occu") {
      cat(sprintf("  Sites: %d, Sources: %d, Species: %d\n",
                  model$n_sites, model$n_sources, model$n_species))
    } else if (model$model_type == "occu_multiscale_cover") {
      cat(sprintf("  Cells: %d, Plots: %d, Max visits: %d\n",
                  model$n_cells, model$n_plots, model$max_visits))
    } else if (model$model_type == "integrated") {
      cat(sprintf("  Sites: %d, Sources: %d\n", model$n_sites, model$n_sources))
    } else if (model$model_type == "jsdm") {
      cat(sprintf("  Sites: %d, Species: %d\n", model$n_sites, model$n_species))
    }
  }
  if (!is.null(x$n_samples)) {
    cat(sprintf("  Samples: %d", x$n_samples))
    if (!is.null(x$n_chains)) {
      cat(sprintf(" (%d chain%s%s)", x$n_chains,
                  if (x$n_chains > 1L) "s" else "",
                  if (!is.null(x$n_thin) && x$n_thin > 1L)
                    sprintf(", thin %d", x$n_thin) else ""))
    }
    if (!is.null(x$epsilon) && !is.na(x$epsilon)) {
      cat(sprintf(", step size: %.4f", x$epsilon))
    }
    cat("\n")
  }
  # Reproducibility: seeds for stochastic routes (NUTS chains / MI / Gibbs).
  if (!is.null(x$seeds)) {
    cat(sprintf("  Seeds: %s\n", paste(x$seeds, collapse = ", ")))
  } else if (!is.null(x$seed)) {
    cat(sprintf("  Seed: %d\n", x$seed))
  }
  conv_shown <- FALSE
  if (!is.null(x$convergence)) {
    rh <- x$convergence$rhat
    eb <- x$convergence$ess_bulk
    if (any(is.finite(rh)) || any(is.finite(eb))) {
      cat(sprintf("  Convergence: max Rhat %.3f, min bulk ESS %.0f\n",
                  max(rh, na.rm = TRUE), min(eb, na.rm = TRUE)))
      conv_shown <- TRUE
      if (any(rh > 1.01, na.rm = TRUE)) {
        cat("    WARNING: Rhat > 1.01 for some parameters; chains may not have ",
            "mixed. Increase n.iter / n.chains.\n", sep = "")
      }
    }
  }
  if (!is.null(x$sla_status) && !identical(x$sla_status, "off")) {
    cat(sprintf("  Marginals: %s\n", x$sla_status))
    clipped <- attr(x$draws, "sla_clipped")
    fallback <- attr(x$draws, "sla_fallback")
    if (length(clipped) > 0) {
      cat(sprintf("    skew clipped to +/-0.95 for: %s\n",
                  paste(clipped, collapse = ", ")))
    }
    if (length(fallback) > 0) {
      cat(sprintf("    fell back to Gaussian for: %s\n",
                  paste(fallback, collapse = ", ")))
    }
  }
  if (!is.null(x$divergent) && isTRUE(sum(x$divergent) > 0)) {
    cat(sprintf("  WARNING: %d divergent transitions\n", sum(x$divergent)))
  }
  # Scalar fallback for a fit carrying only the summary numbers (no per-parameter
  # record); the per-parameter branch above already prints the same line.
  if (!conv_shown && !is.null(x$max_rhat) && is.finite(x$max_rhat)) {
    cat(sprintf("  Convergence: max R-hat %.3f, min ESS %.0f\n",
                x$max_rhat, x$min_ess))
  }
  if (!is.null(x$intercepts)) {
    cat("\n")
    for (nm in names(x$intercepts)) {
      label <- switch(nm,
        psi  = "Mean occupancy (intercept)",
        psi1 = "Mean initial occupancy (intercept)",
        p    = "Mean detection (intercept)",
        gamma   = "Mean colonization (intercept)",
        epsilon = "Mean extinction (intercept)",
        lambda  = "Mean abundance (intercept)"
      )
      if (!is.null(label)) {
        cat(sprintf("%s: %.3f\n", label, x$intercepts[[nm]]))
      }
    }
  }
  if (!is.null(x$nmix_dispersion)) {
    d <- x$nmix_dispersion
    if (isTRUE(is.finite(d$r_sd))) {
      cat(sprintf("NB dispersion (size r): %.3f (SE %.3f)\n", d$r, d$r_sd))
    } else {
      cat(sprintf("NB dispersion (size r): %.3f\n", d$r))
    }
  }
  # Surface attenuated community variance components so the reported between-
  # species spread is not read as unbiased. Means are unaffected.
  va <- x$ms_community$var_attenuation
  if (!is.null(va) && !identical(va$debias, "aghq")) {
    cat(sprintf("  Note: %s\n", va$note))
  }
  invisible(x)
}
