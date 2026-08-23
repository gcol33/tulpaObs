# Pack the outer-grid control knobs into the named list every latent block
# reads. Absent knobs stay absent, so a block keeps its own default grid.
.tobs_outer_grids <- function(sigma = NULL, rho = NULL, tau = NULL,
                              range = NULL) {
  Filter(Negate(is.null),
         list(sigma = sigma, rho = rho, tau = tau, range = range))
}

#' Internal engine entry point
#'
#' Dispatches to Laplace (default) or NUTS for a built `tobs_model`. Not
#' user-facing; called from `tobs()` via the per-family `.dispatch_*` helpers.
#' Spatial / temporal / random-effect / SVC / latent structure is read from
#' the structured terms the formula carried (`model$structured_terms`), not
#' from arguments — there is a single user-facing specification path.
#'
#' @keywords internal
.tobs_fit_model <- function(model,
                            method = c("laplace", "nested_laplace", "nuts",
                                       "pg_gibbs"),
                            priors = NULL,
                            sigma.beta = NULL, sigma.re.scale = 1,
                            max.iter = 100L, tol = 1e-4, damping = 0.7,
                            n.iter = NULL, n.warmup = NULL, n.thin = NULL,
                            n.chains = NULL, n.threads = NULL,
                            sigma.logr = NULL,
                            max.treedepth = NULL, adapt.delta = NULL,
                            seed = NULL,
                            approx = c("gaussian_laplace", "simplified_laplace"),
                            correction = "none",
                            n.gibbs = 10L, n.imputations = 20L,
                            re.aghq = TRUE, n.quad = NULL, re.lkj = NULL,
                            K.max = NULL, mixture = "poisson",
                            integration = c("grid", "ccd"),
                            sigma.grid = NULL, rho.grid = NULL,
                            tau.grid = NULL, range.grid = NULL,
                            verbose = TRUE, ...) {

  method <- match.arg(method)
  approx <- match.arg(approx)
  integration <- match.arg(integration)

  if (!inherits(model, "tobs_model")) {
    stop("model must be a tobs_model object (from `.tobs_build_model()`)")
  }

  # Sampler knobs come from the one engine table. They arrive as NULL sentinels
  # rather than literal formals because this entry serves several engines and
  # the profile is only known once `method` is: a literal formal here would be a
  # second answer to "what is the default n.iter", and it was -- the whole
  # family fitter roster below is handed explicit values from this frame, so
  # their own formals never applied and this was the live answer for every one
  # of them, at 2000 draws against the table's 1000.
  #
  # That 2000 was not a decision. Commit 8975470 fixed `n.iter` from meaning the
  # TOTAL run to meaning kept post-warmup draws on exactly these paths, and left
  # the literal alone -- so a default that had always kept 2000 - 1000 = 1000
  # draws silently began keeping 2000 (that commit's own message records the
  # behaviour change). Reading the table restores the count these paths were
  # calibrated at. The knobs that ARE deliberate -- a wider coefficient prior, a
  # looser adaptation target, a different stream seed -- are on the record as a
  # family row instead (`.TOBS_SINGLE_SPECIES_NUTS`).
  # The `%||%` tails cover a non-sampling engine (`laplace` / `nested_laplace`
  # carry a ridge and no chain knobs) whose branches never read them anyway.
  prof <- .tobs_single_species_defaults(method)
  if (is.null(sigma.beta))    sigma.beta    <- prof$sigma.beta    %||% 10
  if (is.null(n.iter))        n.iter        <- prof$n.iter        %||% 1000L
  if (is.null(n.warmup))      n.warmup      <- prof$n.warmup      %||% 1000L
  if (is.null(n.thin))        n.thin        <- prof$n.thin        %||% 1L
  if (is.null(n.chains))      n.chains      <- prof$n.chains      %||% 1L
  if (is.null(max.treedepth)) max.treedepth <- prof$max.treedepth %||% 10L
  if (is.null(adapt.delta))   adapt.delta   <- prof$adapt.delta   %||% 0.9
  if (is.null(seed))          seed          <- prof$seed          %||% 1L
  # A regularization strength is an engine-level value like the ridge beside
  # it; it is read on the Laplace routes, so the profile is the Laplace one
  # regardless of `method`.
  if (is.null(re.lkj)) re.lkj <- .tobs_default("laplace", "re.lkj")
  # `n.quad` names one control across several marginals; this entry's route is
  # the formula-RE AGHQ debias.
  if (is.null(n.quad))  n.quad  <- .tobs_n_quad("re_aghq")

  # Engine-shaped structure specs derived from the formula's structured terms.
  structs  <- .tobs_structures_from_model(model)
  spatial  <- structs$spatial
  temporal <- structs$temporal
  re       <- structs$re
  svc      <- structs$svc
  latent   <- structs$latent

  # svc() (the continuous NNGP spatially-varying coefficient) is a latent field
  # block on the state arm, so every family whose marginal exposes a per-site eta
  # gradient to the shared areal-BFGS nested-Laplace driver (R/areal_bfgs.R) fits
  # it: single-season occu(), and the observation families removal / distance /
  # fp_occu / dyn_abun. Single-season occu() also samples it on the NUTS path, in
  # the compiled tulpa engine (populate_svc, src/occu_fit.cpp). The remaining
  # fitters -- the N-mixture families, whose areal path is the C++ count-spatial
  # driver rather than the areal-BFGS one, and every family NUTS target -- would
  # silently drop the term and fit a model missing what the user asked for, so
  # those error with a pointer. The areal analogue of a spatially-varying
  # coefficient is the weighted areal bar (spatial(~ 1 + w || cell, graph)), which
  # arrives as a `spatial` term, not `svc`.
  svc_bfgs <- c("removal", "distance", "fp_occu", "dyn_abun")
  svc_wired <-
    (identical(model$model_type, "single") &&
       method %in% c("laplace", "nested_laplace", "nuts")) ||
    (model$model_type %in% svc_bfgs &&
       method %in% c("laplace", "nested_laplace"))
  if (!is.null(svc) && !svc_wired) {
    stop(sprintf(paste0(
      "svc() (a spatially-varying coefficient) is wired for single-season ",
      "occu() under method = \"laplace\" / \"nested_laplace\" / \"nuts\", and for ",
      "removal() / distance() / fp_occu() / dyn_abun() under method = ",
      "\"laplace\" / \"nested_laplace\"; it is silently unsupported for ",
      "model_type = \"%s\"%s. For an areal spatially-varying coefficient use a ",
      "weighted areal bar -- spatial(~ 1 + w || cell, graph = adj) with method = ",
      "\"nested_laplace\" -- which is recovery-tested."),
      model$model_type,
      if (!method %in% c("laplace", "nested_laplace", "nuts"))
        sprintf(" under method = \"%s\"", method) else ""),
      call. = FALSE)
  }

  # Polya-Gamma Gibbs (spOccupancy PGOcc): a REAL MCMC chain over the exact
  # single-season occupancy posterior via PG data augmentation (distinct from
  # method = "laplace_gibbs", the stochastic-EM variance correction). Runs on the
  # natural-scale design (the conjugate Gaussian update needs no autoscaling). v1:
  # single-season, site-level detection, no structured terms (the PG-spatial
  # extensions -- tulpa's pg_binomial_{icar,bym2,...} -- are the documented
  # follow-up).
  if (identical(method, "pg_gibbs")) {
    if (!identical(model$model_type, "single"))
      stop("method = \"pg_gibbs\" is currently wired for single-season occu() ",
           "only.", call. = FALSE)
    if (!is.null(temporal) || !is.null(re) || !is.null(latent) || !is.null(svc))
      stop("method = \"pg_gibbs\" does not yet support temporal / RE / latent / ",
           "svc terms (an areal icar() field IS supported -- spPGOcc, ",
           ").", call. = FALSE)
    if (!is.null(spatial)) {
      # spPGOcc: an intrinsic areal (ICAR) field on the occupancy logit, jointly
      # updated with the coefficients as a Gaussian Markov random field.
      return(.tobs_fit_occu_pg_gibbs_spatial(
        model, spatial, priors = priors, sigma.beta = sigma.beta,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = max(n.chains, 2L),
        n.thin = n.thin, seed = seed, verbose = verbose))
    }
    return(.tobs_fit_occu_pg_gibbs(
      model, priors = priors, sigma.beta = sigma.beta,
      n.iter = n.iter, n.warmup = n.warmup, n.chains = max(n.chains, 2L),
      n.thin = n.thin, seed = seed, verbose = verbose))
  }

  # Autoscale every per-process design matrix before the engine sees it.
  # The engine optimizes on the centered+scaled design; per-process
  # betas / SEs / draws are transformed back to the user-facing natural
  # scale below. `model` (natural-scale) is restored on the returned fit
  # so `fitted()`, `residuals()`, `predict()`, and diagnostics see the
  # same X they would have without this hook.
  scale_info   <- .autoscale_model_X(model)
  fit_model    <- scale_info$model
  scales       <- scale_info$scales
  process_info <- model$process_info

  # Shared tail for the observation-family dispatch branches below. Each fitter
  # ran against the autoscaled `fit_model`, so the coefficients / covariance are
  # transformed back to the natural scale and the unscaled `model` is restored --
  # which drops every slot the fitter set on its own copy, so the fitted field's
  # per-arm eta offset is (re)attached against the restored model here, where a
  # single call reaches every family.
  .tobs_finalize_family_fit <- function(fit) {
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$vcov  <- .unscale_vcov(fit$vcov, scales, process_info)
    fit$model <- model
    fit <- .tobs_attach_model_eta_offset(fit)
    fit <- .tobs_nuts_field_loglik(fit)
    fit$intercepts <- compute_intercepts(model, fit$means)
    fit
  }

  # N-mixture abundance: a closed-form marginal Laplace fit (tulpa owns the
  # likelihood). No EM, no NUTS yet. `method` is "laplace" (non-spatial) or
  # "nested_laplace" (areal spatial offset). Reuses the per-process autoscaler;
  # the coefficient covariance is transformed back to natural scale alongside
  # the means / draws.
  if (identical(model$model_type, "nmix")) {
    # NUTS: sample the exact coefficient posterior of the non-spatial N-mixture
    # via the in-tree C++ FullGradFn over the closed-form marginal (R/abun_nuts.R).
    # Spatial / RE / temporal terms are not yet wired on the sampler (#51).
    if (identical(method, "nuts")) {
      if (!is.null(temporal)) {
        stop("method = \"nuts\" for abun() does not yet support temporal terms ",
             "(#51); use method = \"laplace\".", call. = FALSE)
      }
      if (!is.null(spatial)) {
        # Areal field on the abundance arm sampled by NUTS: a fixed-hyper
        # non-centered icar()/car_proper() field (the field precision fixed at
        # the nested-Laplace estimate), jointly with the coefficients.
        fit <- .tobs_fit_abun_nuts_spatial(
          fit_model, spatial, mixture = mixture, K_max = K.max,
          sigma.beta = sigma.beta, sigma.logr = sigma.logr,
          n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
          n.thin = n.thin, n.threads = n.threads,
          max.treedepth = max.treedepth, adapt.delta = adapt.delta,
          seed = seed, verbose = verbose)
      } else {
      # Random effects: a single intercept RE on one arm samples under NUTS
      # (non-centered per-site offset + log_sigma hyperparameter). Slopes /
      # multi-term / both-arm RE stay on the AGHQ Laplace path.
      fit <- .tobs_fit_abun_nuts(
        fit_model, mixture = mixture, K_max = K.max, sigma.beta = sigma.beta,
        sigma.logr = sigma.logr,
        re = re,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
        n.thin = n.thin, n.threads = n.threads,
        max.treedepth = max.treedepth, adapt.delta = adapt.delta,
        seed = seed, verbose = verbose)
      }
      return(.tobs_finalize_family_fit(fit))
    }
    nmix_method <- if (is.null(spatial)) "laplace" else "nested_laplace"
    fit <- .tobs_fit_nmix(fit_model, method = nmix_method, spatial = spatial,
                          temporal = temporal, re = re, priors = priors,
                          mixture = mixture, K_max = K.max,
                          max_iter = max.iter, tol = tol,
                          n_quad = n.quad, lkj_eta = re.lkj,
                          sigma_beta = sigma.beta,
                          verbose = verbose)
    return(.tobs_finalize_family_fit(fit))
  }

  # Removal sampling: the sequential-depletion abundance marginal (its latent N
  # summed out in closed form, like the N-mixture). Non-spatial fixed effects
  # only this round; "laplace" or "nuts".
  if (identical(model$model_type, "removal")) {
    if (!is.null(temporal))
      .tobs_check_count_temporal(temporal, spatial, method, "removal", "abundance",
                                 allow_temporal_only = TRUE,
                                 allow_nuts_temporal = TRUE)
    if (!is.null(spatial) || !is.null(temporal) || !is.null(svc)) {
      # Areal field on the abundance arm: icar() / car_proper() / bym2() under
      # the nested-Laplace driver, optionally composed with a temporal() block
      # or continuous varying-coefficient surfaces via the shared areal-BFGS
      # driver, or a fixed-hyper non-centered car_proper() field sampled jointly
      # with the coefficients under NUTS. A temporal() or svc() term on its own
      # runs the same areal-BFGS driver with just that block.
      if (identical(method, "nuts")) {
        fit <- .tobs_fit_removal_nuts_spatial(
          fit_model, spatial = spatial, temporal = temporal,
          mixture = mixture, K_max = K.max,
          sigma.beta = sigma.beta, sigma.logr = sigma.logr,
          n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
          n.thin = n.thin, n.threads = n.threads,
          max.treedepth = max.treedepth, adapt.delta = adapt.delta,
          seed = seed, verbose = verbose)
      } else if (is.null(spatial)) {
        fit <- .tobs_fit_removal_spatial_bfgs(fit_model, spatial = NULL,
                                              temporal = temporal, svc = svc,
                                              mixture = mixture,
                                              K_max = K.max, max_iter = max.iter,
                                              tol = tol, verbose = verbose)
      } else {
        fit <- .tobs_fit_removal_spatial(fit_model, spatial, temporal = temporal,
                                         svc = svc, mixture = mixture,
                                         K_max = K.max, max_iter = max.iter,
                                         tol = tol, verbose = verbose)
      }
    } else if (identical(method, "nuts")) {
      fit <- .tobs_fit_removal_nuts(
        fit_model, mixture = mixture, K_max = K.max, sigma.beta = sigma.beta,
        sigma.logr = sigma.logr, re = re,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
        n.thin = n.thin, n.threads = n.threads,
        max.treedepth = max.treedepth, adapt.delta = adapt.delta,
        seed = seed, verbose = verbose)
    } else if (!is.null(re)) {
      # Site-level grouped RE on the abundance OR detection arm via the shared
      # count-model AGHQ path. n_quad = 1 is the joint Laplace (the
      # small-cluster sigma attenuation regime); n_quad > 1 debiases it.
      fit <- .tobs_fit_removal_re(fit_model, re = re, mixture = mixture,
                                  K_max = K.max, max_iter = max.iter, tol = tol,
                                  n_quad = n.quad, lkj_eta = re.lkj,
                                  theta_prior_sd = sigma.beta, verbose = verbose)
    } else {
      fit <- .tobs_fit_removal(fit_model, mixture = mixture, K_max = K.max,
                               max_iter = max.iter, tol = tol, verbose = verbose)
    }
    return(.tobs_finalize_family_fit(fit))
  }

  # Distance sampling: the binned multinomial-over-N marginal (its latent N
  # summed out in closed form, like the N-mixture). Non-spatial fixed effects
  # only this round; "laplace" or "nuts".
  if (identical(model$model_type, "distance")) {
    if (!is.null(temporal))
      .tobs_check_count_temporal(temporal, spatial, method, "distance", "abundance",
                                 allow_temporal_only = TRUE,
                                 allow_nuts_temporal = TRUE)
    if (!is.null(spatial) || !is.null(temporal) || !is.null(svc)) {
      # Areal field on the abundance arm: icar() / car_proper() (half-normal or
      # hazard key) under the nested-Laplace driver, optionally composed with a
      # temporal() block or continuous varying-coefficient surfaces via the shared
      # areal-BFGS driver, or a fixed-hyper non-centered car_proper() field
      # sampled under NUTS (half-normal key). A temporal() or svc() term on its
      # own runs the same areal-BFGS driver with just that block.
      if (identical(method, "nuts")) {
        fit <- .tobs_fit_distance_nuts_spatial(
          fit_model, spatial = spatial, temporal = temporal,
          mixture = mixture, K_max = K.max,
          sigma.beta = sigma.beta,
          n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
          n.thin = n.thin, n.threads = n.threads,
          max.treedepth = max.treedepth, adapt.delta = adapt.delta,
          seed = seed, verbose = verbose)
      } else {
        fit <- .tobs_fit_distance_spatial(fit_model, spatial, temporal = temporal,
                                          svc = svc, mixture = mixture,
                                          K_max = K.max, max_iter = max.iter,
                                          tol = tol, verbose = verbose,
                                          integration = integration)
      }
    } else if (identical(method, "nuts")) {
      fit <- .tobs_fit_distance_nuts(
        fit_model, mixture = mixture, K_max = K.max, sigma.beta = sigma.beta,
        sigma.logr = sigma.logr, re = re,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
        n.thin = n.thin, n.threads = n.threads,
        max.treedepth = max.treedepth, adapt.delta = adapt.delta,
        seed = seed, verbose = verbose)
    } else if (!is.null(re)) {
      # Site-level grouped RE on the abundance arm via the shared count-model
      # AGHQ path (half-normal key, abundance arm only). n_quad = 1 is the joint
      # Laplace, n_quad > 1 debiases the small-cluster attenuation.
      fit <- .tobs_fit_distance_re(fit_model, re = re, mixture = mixture,
                                   K_max = K.max, max_iter = max.iter, tol = tol,
                                   n_quad = n.quad, lkj_eta = re.lkj,
                                   theta_prior_sd = sigma.beta, verbose = verbose)
    } else {
      fit <- .tobs_fit_distance(fit_model, mixture = mixture, K_max = K.max,
                                max_iter = max.iter, tol = tol, verbose = verbose)
    }
    return(.tobs_finalize_family_fit(fit))
  }

  # Open-population (Dail-Madsen) N-mixture: the latent abundance sequence summed
  # out by an exact HMM forward recursion (not closed form). Non-spatial fixed
  # effects only this round; "laplace" or "nuts".
  if (identical(model$model_type, "dyn_abun")) {
    # Zero-inflated open N-mixture (zip / zinb): a pure-R structural-zero layer
    # over the Dail-Madsen marginal. v1 is non-spatial laplace with an
    # intercept-only structural-zero probability; a field, an RE, or NUTS stay
    # Poisson / negbin.
    if (model$mixture %in% c("zip", "zinb")) {
      if (!is.null(spatial) || !is.null(temporal) || !is.null(re) ||
          identical(method, "nuts")) {
        stop("Zero-inflated open N-mixture (zip / zinb) does not yet compose ",
             "with a spatial field, a temporal term, a random effect, or NUTS; ",
             "use mixture = \"poisson\" / \"negbin\" for those, or drop the term.",
             call. = FALSE)
      }
      fit <- .tobs_fit_dyn_abun_zip(fit_model, max_iter = 300L, verbose = verbose)
      return(.tobs_finalize_family_fit(fit))
    }
    if (!is.null(temporal))
      .tobs_check_count_temporal(temporal, spatial, method, "dyn_abun",
                                 "initial-abundance", allow_temporal_only = TRUE,
                                 allow_nuts_temporal = TRUE)
    if (!is.null(spatial) || !is.null(temporal) || !is.null(svc)) {
      # Areal field on the initial-abundance arm: icar() / car_proper() under the
      # nested-Laplace forward-HMM driver, optionally composed with a temporal()
      # block or continuous varying-coefficient surfaces via the shared areal-BFGS
      # driver, or a fixed-hyper non-centered car_proper() field sampled jointly
      # with the coefficients under NUTS. A temporal() or svc() term on its own
      # runs the same areal-BFGS driver with just that block (#114, #144), or --
      # for temporal -- a fixed-hyper non-centered temporal field on the NUTS
      # field block (#114).
      if (identical(method, "nuts")) {
        fit <- if (is.null(spatial))
          .tobs_fit_dyn_abun_nuts_temporal(
            fit_model, temporal, mixture = model$mixture %||% "poisson",
            K_max = K.max, sigma.beta = sigma.beta,
            n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
            n.thin = n.thin, n.threads = n.threads,
            max.treedepth = max.treedepth, adapt.delta = adapt.delta,
            seed = seed, verbose = verbose)
        else .tobs_fit_dyn_abun_nuts_spatial(
          fit_model, spatial, mixture = model$mixture %||% "poisson",
          K_max = K.max, sigma.beta = sigma.beta,
          n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
          n.thin = n.thin, n.threads = n.threads,
          max.treedepth = max.treedepth, adapt.delta = adapt.delta,
          seed = seed, verbose = verbose)
      } else {
        fit <- .tobs_fit_dyn_abun_spatial(fit_model, spatial, temporal = temporal,
                                          svc = svc,
                                          mixture = model$mixture %||% "poisson",
                                          K_max = K.max, max_iter = 300L, tol = 1e-8,
                                          verbose = verbose, integration = integration)
      }
    } else if (identical(method, "nuts")) {
      fit <- .tobs_fit_dyn_abun_nuts(
        fit_model, sigma.beta = sigma.beta, re = re,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
        n.thin = n.thin, n.threads = n.threads,
        max.treedepth = max.treedepth, adapt.delta = adapt.delta,
        seed = seed, verbose = verbose)
    } else if (!is.null(re)) {
      # Site-level grouped RE on the initial-abundance (lambda) arm via the
      # exact HMM-forward AGHQ path. n_quad = 1 is the joint Laplace (the
      # small-cluster sigma attenuation regime); n_quad > 1 debiases it.
      fit <- .tobs_fit_dyn_abun_re(fit_model, re = re, max_iter = 300L,
                                   tol = 1e-8, verbose = verbose,
                                   n_quad = n.quad, lkj_eta = re.lkj,
                                   theta_prior_sd = sigma.beta)
    } else {
      fit <- .tobs_fit_dyn_abun(fit_model, max_iter = 300L, tol = 1e-8,
                                verbose = verbose)
    }
    return(.tobs_finalize_family_fit(fit))
  }

  # False-positive occupancy: the Miller et al. (2011) multistate marginal (its
  # latent occupancy z summed out in closed form). Non-spatial fixed effects only
  # this round; "laplace" (analytic-gradient BFGS over the exact marginal) or
  # "nuts".
  if (identical(model$model_type, "fp_occu")) {
    if (!is.null(temporal))
      .tobs_check_count_temporal(temporal, spatial, method, "fp_occu", "occupancy",
                                 allow_temporal_only = TRUE,
                                 allow_nuts_temporal = TRUE)
    if (!is.null(spatial) || !is.null(temporal) || !is.null(svc)) {
      # Areal field on the occupancy (psi) arm: icar() / car_proper() under the
      # nested-Laplace two-state driver, optionally composed with a temporal()
      # block or continuous varying-coefficient surfaces via the shared areal-BFGS
      # driver, or a fixed-hyper non-centered car_proper() field sampled jointly
      # with the coefficients under NUTS. A temporal() or svc() term on its own
      # runs the same areal-BFGS driver with just that block (#114, #144).
      if (identical(method, "nuts")) {
        fit <- .tobs_fit_fp_occu_nuts_spatial(
          fit_model, spatial = spatial, temporal = temporal, sigma.beta = sigma.beta,
          n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
          n.thin = n.thin, n.threads = n.threads,
          max.treedepth = max.treedepth, adapt.delta = adapt.delta,
          seed = seed, verbose = verbose)
      } else {
        fit <- .tobs_fit_fp_occu_spatial(fit_model, spatial, temporal = temporal,
                                         svc = svc, max_iter = max.iter,
                                         tol = 1e-8, verbose = verbose,
                                         integration = integration)
      }
    } else if (identical(method, "nuts")) {
      fit <- .tobs_fit_fp_occu_nuts(
        fit_model, sigma.beta = sigma.beta, re = re,
        n.iter = n.iter, n.warmup = n.warmup, n.chains = n.chains,
        n.thin = n.thin, n.threads = n.threads,
        max.treedepth = max.treedepth, adapt.delta = adapt.delta,
        seed = seed, verbose = verbose)
    } else if (!is.null(re)) {
      # Site-level grouped RE on the occupancy (psi) arm via the pure-R make_site
      # AGHQ path. n_quad = 1 is the joint Laplace, n_quad > 1 debiases the
      # small-cluster variance-component attenuation.
      fit <- .tobs_fit_fp_occu_re(fit_model, re = re, max_iter = max.iter,
                                  tol = tol, n_quad = n.quad, lkj_eta = re.lkj,
                                  sigma.beta = sigma.beta, verbose = verbose)
    } else {
      fit <- .tobs_fit_fp_occu(fit_model, max_iter = 500L, tol = 1e-8,
                               sigma.beta = NULL, verbose = verbose)
    }
    return(.tobs_finalize_family_fit(fit))
  }

  # Continuous NNGP varying coefficient(s) on single-season occupancy under the
  # deterministic backends. The surfaces are latent field blocks on the
  # occupancy logit, so the fit rides the shared areal-BFGS nested-Laplace
  # driver: the field hyperparameters (marginal SD, range) are integrated on an
  # outer grid either way, so `laplace` and `nested_laplace` land on the same
  # fitter and the fit reports `nested_laplace`. A structured term alongside
  # svc() would need a second field family in the same block list and is not
  # wired, so it errors instead of being dropped.
  if (!is.null(svc) && method %in% c("laplace", "nested_laplace")) {
    if (!is.null(spatial) || !is.null(temporal) || !is.null(re)) {
      stop("occu() + svc() on the Laplace backends fits the varying-coefficient ",
           "surface(s) alone; a spatial() / temporal() / re() term alongside it ",
           "is not wired. Use method = \"nuts\", or drop the extra term. ",
           "", call. = FALSE)
    }
    fit <- do.call(.tobs_fit_occu_svc, c(
      list(model = fit_model, svc = svc, priors = priors,
           max_iter = as.integer(max.iter), tol = 1e-8, verbose = verbose),
      list(...)))
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$model      <- model
    # The fitted surfaces contribute a per-site offset on the occupancy logit
    # (sum_k X[, index_k] * z_k), carried on the natural-scale model so the
    # in-sample fitted() reads it -- the surfaces are latent, so the offset is
    # scale-invariant. Attached after the model swap, which drops any slot the
    # fitter set on its own copy.
    fit$model$occ_eta_offset <- fit$svc_eta_offset
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  if (method == "laplace") {
    # .tobs_laplace() consumes the spatial and re terms; it has no temporal
    # channel, so a temporal() term reaching here would be dropped and the fit
    # would silently omit the field the user asked for. The temporal field is
    # assembled as a latent block by the nested-Laplace path
    # (.tobs_block_from_temporal, R/em_nested_laplace.R), where it composes with
    # the areal field and the iid blocks, and is populated directly on the
    # sampler (populate_temporal, src/occu_fit.cpp).
    if (!is.null(temporal)) {
      stop(paste0(
        "temporal() is not consumed by the Laplace engine (method = \"laplace\" ",
        "/ \"laplace_sla\" / \"laplace_gibbs\" / \"laplace_mi\"); the temporal ",
        "field is grid-integrated under method = \"nested_laplace\" (where it ",
        "composes with an areal field and re() blocks) and sampled under ",
        "method = \"nuts\". Re-fit with one of those, or drop the term."),
        call. = FALSE)
    }
    fit <- .tobs_laplace(fit_model, spatial = spatial, re = re,
                         priors = priors,
                         sigma_beta = sigma.beta,
                         max_iter = max.iter, tol = tol, damping = damping,
                         approx = approx,
                         correction = correction,
                         n_imputations = n.imputations, n_gibbs = n.gibbs,
                         seed = seed,
                         re_aghq = re.aghq, n_quad = n.quad, lkj_eta = re.lkj,
                         verbose = verbose)
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$model      <- model
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  if (method == "nested_laplace") {
    # Areal count / relative-abundance GLMM: the response is observed directly
    # (no latent state), so a count areal fit is a single tulpa_nested_laplace()
    # call over the count block with the areal field as its prior -- not the
    # occupancy EM. Grid-integrated fixed effects + field.
    if (identical(model$model_type, "count")) {
      fit <- .tobs_fit_count_spatial(fit_model, spatial,
                                     max_iter = as.integer(max.iter), tol = tol,
                                     sigma.grid = sigma.grid,
                                     rho.grid = rho.grid, tau.grid = tau.grid,
                                     range.grid = range.grid,
                                     verbose = verbose)
      fit <- .unscale_fit_per_process(fit, scales, process_info)
      fit$model      <- model
      # The latent field is a per-site eta offset (scale-invariant), carried on the
      # model so the field-aware fitted() / WAIC read it. Swapping in the unscaled
      # model drops any slot the fitter attached, so the offset is (re)built here
      # against that model. With a varying-coefficient bar the contribution is the
      # WEIGHTED sum over the intercept and trend fields, sum_k W[i,k] f_k[i] --
      # not the intercept field alone; a plain intercept field reduces to
      # spatial_field.
      fit$model$count_field_offset <-
        .tobs_spatial_field_offset(fit, spatial, model)
      fit$intercepts <- compute_intercepts(model, fit$means)
      return(fit)
    }
    # Standalone occu() varying-coefficient (SVC) spatial bar: route through the
    # joint direct-grid engine, single-arm (occupancy + detection, no cover arm),
    # which integrates the field hyperparameters on a direct outer grid. The EM
    # fixed-point path oscillates / does not converge on this case at EVA scale;
    # the reroute is scoped to the SVC occupancy fit the single-arm joint engine
    # covers (`.tobs_occu_reroute_to_joint`). The joint route consumes the
    # autoscaled `fit_model`; the per-process betas / SEs / draws are transformed
    # back to natural scale by the shared unscale below, exactly as the EM path's
    # output is.
    if (.tobs_occu_reroute_to_joint(fit_model, spatial, temporal, re)) {
      fields <- .tobs_resolve_occu_spatial_fields(spatial, fit_model)
      fit <- do.call(.tobs_fit_occu_joint, c(
        list(model = fit_model, fields = fields, priors = priors,
             max.iter = as.integer(max.iter), tol = tol, verbose = verbose),
        list(...)
      ))
      fit <- .unscale_fit_per_process(fit, scales, process_info)
      # Restore the natural-scale model for fitted()/predict()/diagnostics, but
      # carry the field geometry (the graph -> cell-node count and site -> node
      # map) the joint route resolved, which field-aware predict reads.
      model$n_cells   <- fit$model$n_cells
      model$site_cell <- fit$model$site_cell
      fit$model       <- model
      fit$intercepts  <- compute_intercepts(model, fit$means)
      return(fit)
    }
    # Nested-Laplace path: single-season, integrated, or dynamic
    # occupancy. The driver builds a multi-block latent prior from spatial +
    # temporal + re and attaches it to the state ("occ") M-step block, which
    # tulpa::tulpa_em_laplace() routes through tulpa::tulpa_nested_laplace().
    nl_max_iter <- min(as.integer(max.iter), 25L)
    # INLA NA-response prediction: single-season sites with an all-missing
    # detection history are held out of the likelihood (n_trials = 0) but stay
    # in the design so the latent field interpolates their occupancy from
    # neighbours / shared structure (`.tobs_heldout_sites()`).
    heldout_state <- .tobs_heldout_sites(model)
    fit <- .tobs_em_nested_laplace(
      model    = fit_model,
      spatial  = spatial,
      temporal = temporal,
      re       = re,
      priors   = priors,
      sigma_beta = sigma.beta,
      max_iter = nl_max_iter,
      tol      = tol,
      damping  = damping,
      heldout_state = heldout_state,
      # Outer-grid overrides, the same `control` names the joint cover route
      # takes; the validator already admits them on a nested-Laplace route.
      grids    = .tobs_outer_grids(sigma.grid, rho.grid, tau.grid, range.grid),
      verbose  = verbose
    )
    fit <- .unscale_fit_per_process(fit, scales, process_info)
    fit$model      <- model
    fit$intercepts <- compute_intercepts(model, fit$means)
    return(fit)
  }

  # Only the NUTS path reaches here (laplace / nested_laplace returned above).
  # The occupancy NUTS spec threads `sigma_beta` but not a user `priors` object,
  # so warn instead of silently ignoring a supplied prior. `priors = FALSE`
  # ("disable the penalty") and the default `NULL` are not user priors.
  if (!is.null(priors) && !isFALSE(priors)) {
    warning("`priors` is not applied on the occupancy NUTS path and was ",
            "ignored; the sampler uses the weakly-informative N(0, sigma.beta) ",
            "coefficient prior (set via control$sigma.beta). Use ",
            "method = \"laplace\" to apply an occu_priors() penalty.",
            call. = FALSE)
  }

  # spatial / temporal / re / svc / latent are produced by
  # .tobs_structures_from_model(), which guarantees their classes; no
  # user-input validation is needed here.
  model_type <- model$model_type

  # ---- Build spec list for C++ ----
  spec <- list(
    model_type = model_type,
    sigma_beta = sigma.beta,
    n_iter = as.integer(n.iter + n.warmup),    # engine total = warmup + kept samples
    n_warmup = as.integer(n.warmup),
    max_treedepth = as.integer(max.treedepth),
    adapt_delta = adapt.delta,
    seed = as.integer(seed),
    verbose = verbose
  )

  # ---- Process design matrices (autoscaled) ----
  spec$X_processes <- fit_model$X_processes

  # ---- Process names for column labels ----
  spec$process_names <- lapply(model$process_info, function(pi) {
    paste0(pi$name, "_", pi$coef_names)
  })

  # ---- Model-type-specific fields ----
  if (model_type == "single") {
    spec$y <- model$y
    if (!is.null(model$X_det_visit)) {
      spec$X_det_visit <- model$X_det_visit
      spec$extra_param_names <- paste0("p_visit_", model$det_visit_names)
    }

  } else if (model_type == "dynamic") {
    spec$y_flat <- model$y_flat
    spec$n_visits <- model$n_visits
    spec$any_detected <- model$any_detected
    spec$n_sites <- model$n_sites
    spec$n_seasons <- model$n_seasons
    spec$max_visits <- model$max_visits

  } else if (model_type == "integrated") {
    spec$y_sources <- model$y_sources
    spec$site_maps <- model$site_maps
    spec$n_sources <- model$n_sources
    spec$n_sites <- model$n_sites

  }

  # ---- Spatial ----
  if (!is.null(spatial)) {
    spatial_params <- build_spatial_params(spatial, model$n_sites)
    spec$spatial_params <- spatial_params
  }

  # ---- Temporal ----
  # Index codes were resolved when the temporal() term was constructed.
  if (!is.null(temporal)) {
    temp_spec <- list(type = temporal$type, shared = temporal$shared,
                      cyclic = temporal$cyclic,
                      tau_shape = temporal$tau_shape,
                      tau_rate = temporal$tau_rate,
                      time_idx = temporal$time_idx,
                      n_times  = temporal$n_times)
    if (!is.null(temporal$group_idx)) {
      temp_spec$group_idx <- temporal$group_idx
      temp_spec$n_groups  <- temporal$n_groups
    }
    spec$temporal_spec <- temp_spec
  }

  # ---- Random effects ----
  if (!is.null(re)) {
    # Accept single tobs_re or list of tobs_re
    if (inherits(re, "tobs_re")) re <- list(re)
    re_spec <- build_re_spec(re, model)
    spec$re_spec <- re_spec
  }

  # ---- SVC ----
  if (!is.null(svc)) {
    # Build X_svc from the design matrix columns. Pull from the
    # (autoscaled) `fit_model` so the SVC base column values match the
    # global beta's parameterization seen by the optimizer. For svc on
    # the intercept the values are 1.0 in both natural and scaled spaces;
    # for svc on a non-intercept numeric column the per-location offsets
    # land on the scaled-column scale.
    X_occ <- fit_model$X_processes[[1]]
    # `coefficients = ` names the columns, `indices = ` gives their positions;
    # both resolve here against the same design.
    svc_cols <- .tobs_svc_columns(svc, X_occ, "occu")
    svc_indices_0based <- svc_cols - 1L  # C++ 0-based
    X_svc_flat <- numeric(nrow(X_occ) * svc$n_svc)
    for (j in seq_along(svc_cols)) {
      col <- svc_cols[j]
      for (i in seq_len(nrow(X_occ))) {
        X_svc_flat[(i - 1) * svc$n_svc + j] <- X_occ[i, col]
      }
    }
    spec$svc_spec <- list(
      n_obs = svc$n_obs, n_svc = svc$n_svc, nn = svc$nn,
      coords = svc$coords, svc_indices = svc_indices_0based,
      X_svc = X_svc_flat,
      nn_idx = svc$nn_idx, nn_dist = svc$nn_dist,
      nn_order = svc$nn_order, nn_order_inv = svc$nn_order_inv,
      cov_type = svc$cov_type, shared = svc$shared,
      sigma2_prior_scale = svc$sigma2_prior_scale,
      phi_prior_U = svc$prior_range[1],
      phi_prior_alpha = svc$prior_range[2]
    )
  }

  # ---- Latent factors ----
  if (!is.null(latent)) {
    spec$latent_spec <- list(
      n_factors = latent$n_factors,
      shared = latent$shared,
      constraint = latent$constraint,
      sigma_prior_rate = latent$sigma_prior_rate
    )
  }

  # ---- Run NUTS: one or more chains, pooled in R (see nuts_chains.R) ----
  fit <- .tobs_run_chains(spec, n.chains = n.chains, n.thin = n.thin,
                          n.threads = n.threads, verbose = verbose)

  # Unscale per-process beta slices in means / sds / draws (the engine
  # optimized on the centered+scaled design).
  fit <- .unscale_fit_per_process(fit, scales, process_info)

  # ---- Build R parameter names ----
  param_names <- unlist(spec$process_names)
  if (!is.null(model$det_visit_names) && length(model$det_visit_names) > 0) {
    param_names <- c(param_names, paste0("p_visit_", model$det_visit_names))
  }

  # Name the random-effect block (log_sigma / chol / z, type-blocked per
  # tulpa's layout) and reconstruct per-group BLUPs into `re_effects` so
  # summary() / ranef() label them instead of showing param[i]. Counts
  # and positions are unchanged.
  if (!is.null(re)) {
    re_design <- .tobs_re_design(if (inherits(re, "tobs_re")) list(re) else re,
                                 model)
    n_lead <- length(param_names)
    re_nms <- .tobs_re_nuts_param_names(re_design)
    if (n_lead + length(re_nms) <= length(fit$means)) {
      param_names <- c(param_names, re_nms)
      names(fit$means)[seq_along(param_names)] <- param_names
      if (!is.null(fit$draws) && ncol(fit$draws) >= length(param_names)) {
        colnames(fit$draws)[seq_along(param_names)] <- param_names
      }
      fit$re_effects <- tryCatch(
        .tobs_re_nuts_effects(fit$draws, re_design, n_lead),
        error = function(e) NULL)
    }
  }
  fit$param_names <- param_names

  # ---- Compute probability-scale intercepts (on natural-scale means) ----
  fit$intercepts <- compute_intercepts(model, fit$means)

  fit$model <- model
  fit$spatial <- spatial
  fit$temporal <- temporal
  fit$re <- re
  fit$svc <- svc
  if (!is.null(svc)) {
    # `svc_cols` was resolved above against the design the surfaces were packed
    # against; reporting it rather than re-resolving keeps the reported columns
    # and the sampled ones the same object.
    fit$svc_indices <- svc_cols
    fit$svc_coefficients <- svc$coefficients
  }
  # The fitted per-location SVC surface. `fit$svc` is the term the user passed
  # in; `fit$svc_field` is what was estimated, mirroring `fit$spatial_field` on
  # the areal path. Without it the surface was in the draws but unreadable, so
  # the term could only ever be smoke-tested.
  fit$svc_field <- .tobs_svc_field(fit)
  # The fitted areal (icar/car_proper) field, mirroring the nested-Laplace
  # `fit$spatial_field`. NULL when there is no areal term.
  fit$spatial_field <- .tobs_areal_field(fit)
  # The fitted continuous GP field. Same slot as the areal surface -- they are
  # alternative spatial terms, never both present.
  if (is.null(fit$spatial_field)) fit$spatial_field <- .tobs_gp_field(fit)
  # Every sampled spatial surface enters the psi logit in the sampler, so it has
  # to enter the linear predictor `fitted()` rebuilds as well. Before this the
  # slot was set only on the Laplace svc route, so a NUTS fit with ANY spatial
  # term reported one psi per site -- the field was influencing the fit and
  # silently absent from everything read off it.
  fit$model$occ_eta_offset <- .tobs_nuts_occ_offset(fit, model)
  fit$latent <- latent
  # Resolved per-chain seeds (chain c used seed + c - 1) for reproducibility.
  fit$seeds <- as.integer(seed) + seq_len(as.integer(n.chains)) - 1L
  # Expose process_info at top level for tulpa generic S3 methods
  fit$process_info <- model$process_info
  # Cross-chain convergence diagnostics (Rhat / bulk + tail ESS), through the
  # writer every sampled path shares, so the record has one shape across
  # families. Computed on the named, natural-scale draws (Rhat / ESS are
  # scale-invariant).
  class(fit) <- c("tobs_fit", "tulpa_fit")
  fit <- .tobs_nuts_attach_convergence(
    fit, .tobs_nuts_chains_from_ids(fit$draws, fit$chain_id),
    par_names = colnames(fit$draws), n_iter = as.integer(n.iter))
  fit
}

# Fitted per-location SVC surface from a NUTS fit.
#
# The engine exports the SVC block's offsets on `ParamLayout` and the fitter
# returns them as `svc_layout`, so the surface is sliced by position rather than
# by parsing column names. Weights are stored `w_flat[j * n_obs + i]` (j indexes
# the SVC term, i the location), which is the transpose of the `X_svc` stride --
# hence the byrow = FALSE fill below.
#
# Returns an `n_obs x n_svc` matrix of posterior means (a bare vector when a
# single coefficient varies, matching `fit$spatial_field` on the areal path),
# carrying the per-draw surface as attribute "draws" for interval work. NULL
# when the fit has no SVC term or the backend did not report a layout (only the
# single-season occu NUTS path does).
.tobs_svc_field <- function(fit) {
  lay <- fit$svc_layout
  if (is.null(lay) || is.null(fit$means)) return(NULL)
  n_svc <- as.integer(lay$n_svc)
  n_obs <- as.integer(lay$n_obs)
  cols  <- seq.int(as.integer(lay$w_start), as.integer(lay$w_end))
  if (length(cols) != n_svc * n_obs) return(NULL)   # HSGP basis, not a surface
  if (length(fit$means) < max(cols)) return(NULL)

  surf <- matrix(as.numeric(fit$means[cols]), nrow = n_obs, ncol = n_svc)
  if (!is.null(fit$draws) && ncol(fit$draws) >= max(cols)) {
    attr(surf, "draws") <- fit$draws[, cols, drop = FALSE]
  }
  if (n_svc == 1L) {
    out <- as.numeric(surf)
    attr(out, "draws") <- attr(surf, "draws")
    return(out)
  }
  surf
}

# Fitted areal field from a single-season occupancy NUTS fit.
#
# The engine exports the field block's offsets on `ParamLayout` (emitted as
# `spatial_layout` in occu_fit.cpp), so the field is sliced by position rather
# than parsed from names. For icar / car_proper the field node enters the logit
# linear predictor directly, so the posterior-mean node IS the per-cell surface;
# the intrinsic level is confounded with the intercept, so it is centred (as the
# nested-Laplace areal summary is). Returns a length-`n_units` vector with the
# per-draw centred field on attribute "draws", or NULL when there is no areal
# term. bym2 names its blocks (`spatial_field` / `spatial_theta`) but its field
# is the Riebler rho-mix of both blocks scaled by the graph scale factor, so it
# is left to the named draws here rather than reconstructed on a partial block.
.tobs_areal_field <- function(fit) {
  lay <- fit$spatial_layout
  if (is.null(lay) || is.null(fit$means)) return(NULL)
  if (!identical(lay$type, "icar") && !identical(lay$type, "car_proper")) return(NULL)
  cols <- seq.int(as.integer(lay$field_start), as.integer(lay$field_end))
  if (length(fit$means) < max(cols)) return(NULL)
  field <- as.numeric(fit$means[cols])
  field <- field - mean(field)
  if (!is.null(fit$draws) && ncol(fit$draws) >= max(cols)) {
    d <- fit$draws[, cols, drop = FALSE]
    attr(field, "draws") <- d - rowMeans(d)
  }
  field
}

# Fitted continuous GP field from a NUTS fit.
#
# Sliced by position off `gp_layout`, the same way the areal and SVC surfaces
# are. NOT centred, unlike the intrinsic areal field: a GP's prior is
# N(0, sigma2 K), which is proper, so its level is identified rather than
# confounded with the intercept. NULL when there is no GP term, or when the fit
# collapsed the field (marginalised out, so there are no draws to report).
.tobs_gp_field <- function(fit) {
  lay <- fit$gp_layout
  if (is.null(lay) || is.null(fit$means)) return(NULL)
  if (isTRUE(lay$collapsed) || is.na(lay$field_start)) return(NULL)
  cols <- seq.int(as.integer(lay$field_start), as.integer(lay$field_end))
  if (length(fit$means) < max(cols)) return(NULL)
  field <- as.numeric(fit$means[cols])
  if (!is.null(fit$draws) && ncol(fit$draws) >= max(cols)) {
    attr(field, "draws") <- fit$draws[, cols, drop = FALSE]
  }
  field
}

# The per-site occupancy-logit offset a sampled spatial surface contributes, for
# the in-sample linear predictor `fitted()` rebuilds.
#
# The RAW field, not the centred one reported as `fit$spatial_field`: the
# sampler's intercept was drawn against the uncentred surface, so subtracting
# its mean here would shift every psi by that constant. Both the areal and GP
# blocks map observations to units 1:1 on this backend (`spatial_group[i] = i+1`
# / `obs_to_loc[i] = i` in populate_helpers.h), so the field vector IS the
# per-site offset and no map is applied.
#
# NULL when there is nothing to add, which leaves eta exactly as it was:
# no spatial term, a collapsed GP, bym2 (whose surface is the rho-mix of two
# blocks and is left to the named draws, as `.tobs_areal_field` documents), or a
# unit count that does not match the design's rows.
.tobs_nuts_occ_offset <- function(fit, model) {
  lay <- fit$spatial_layout %||% fit$gp_layout
  if (is.null(lay) || is.null(fit$means)) return(NULL)
  if (!is.null(fit$spatial_layout) &&
      !identical(fit$spatial_layout$type, "icar") &&
      !identical(fit$spatial_layout$type, "car_proper")) return(NULL)
  if (isTRUE(lay$collapsed) || is.na(lay$field_start)) return(NULL)
  cols <- seq.int(as.integer(lay$field_start), as.integer(lay$field_end))
  if (length(fit$means) < max(cols)) return(NULL)
  off <- as.numeric(fit$means[cols])
  X_occ <- model$X_processes[[1L]]
  if (is.null(X_occ) || nrow(X_occ) != length(off)) return(NULL)
  off
}

# Translate the structured terms a formula carried (`model$structured_terms`)
# into the engine-shaped specs the fitter consumes. A term's process
# membership (`$processes`, set by the formula parser) becomes the length-2
# `shared = c(occ, det)` vector the C++ engine and Laplace paths read.
#
# The engine carries two sharing booleans (occupancy/state + detection), so a
# structured effect on a later process (colonization / extinction / extra
# integrated sources) is rejected rather than silently dropped. At most one
# spatial / temporal / svc / latent term is allowed; random effects may be
# multiple.
.tobs_structures_from_model <- function(model) {
  out <- list(spatial = NULL, temporal = NULL, re = NULL,
              svc = NULL, latent = NULL)
  terms <- model$structured_terms
  if (is.null(terms) || length(terms) == 0L) return(out)

  proc_names <- vapply(model$process_info, function(pi) pi$name, character(1))
  re_list <- list()
  spatial_list <- list()

  one <- function(slot, label) {
    if (!is.null(out[[slot]])) {
      stop(sprintf("Only one %s term is supported per model.", label),
           call. = FALSE)
    }
  }

  for (t in terms) {
    spec  <- t$spec
    procs <- t$processes
    if (any(procs > 2L)) {
      late <- proc_names[procs[procs > 2L]]
      stop(sprintf(
        "Structured term `%s` enters the '%s' predictor; structured effects ",
        spec$label %||% class(spec)[1], paste(late, collapse = "', '")),
        "are supported on the occupancy/state and detection predictors only.",
        call. = FALSE)
    }
    shared <- c(1L %in% procs, 2L %in% procs)

    if (inherits(spec, "tobs_spatial")) {
      # A varying-coefficient bar (one `tobs_spatial` with `is_bar`) or several
      # areal terms (an unweighted intercept field plus weighted SVC fields)
      # describe one multi-field spatial structure. Collect every spatial spec;
      # `.tobs_collect_spatial()` returns either the single plain spec (back-
      # compat: laplace / abun / removal / etc. read `$type` / `$shared`
      # directly) or a combined `is_multifield` container the occu() nested-
      # Laplace path expands into one block per field. A weighted SVC field on a
      # non-nested consumer is still rejected there via
      # `.tobs_reject_weighted_spatial()`.
      spec$shared <- shared
      spatial_list[[length(spatial_list) + 1L]] <- spec
    } else if (inherits(spec, "tobs_temporal")) {
      one("temporal", "temporal"); spec$shared <- shared; out$temporal <- spec
    } else if (inherits(spec, "tobs_svc")) {
      one("svc", "svc"); spec$shared <- shared; out$svc <- spec
    } else if (inherits(spec, "tobs_latent")) {
      one("latent", "latent"); spec$shared <- any(shared); out$latent <- spec
    } else if (inherits(spec, "tobs_re")) {
      spec$shared <- shared
      re_list[[length(re_list) + 1L]] <- spec
    }
  }
  if (length(re_list)) out$re <- re_list
  if (length(spatial_list)) out$spatial <- .tobs_collect_spatial(spatial_list)
  out
}

# Collapse the spatial term(s) on a formula into the single `spatial` slot the
# fitters read. One plain unweighted areal / continuous term passes through
# unchanged (the common case; every fitter reads its `$type` / `$shared`). A
# varying-coefficient bar (`is_bar`) or several areal terms (an intercept field
# plus weighted SVC fields) describe one multi-field spatial structure: they are
# wrapped in a combined `tobs_spatial` carrying `is_multifield = TRUE` and the
# ordered `fields` list, with the intercept field's `type` / `shared` / graph at
# the top level so a consumer that does not understand multi-field spatial still
# dispatches on `$type` and rejects the SVC field through
# `.tobs_reject_weighted_spatial()`. The occu() nested-Laplace path is the one
# consumer that expands `fields` into one latent block each.
.tobs_collect_spatial <- function(spatial_list) {
  if (length(spatial_list) == 1L) {
    s <- spatial_list[[1L]]
    if (!isTRUE(s$is_bar)) return(s)
  }
  # Multiple terms must agree on which arms they enter (one shared field on a
  # given arm); take the first field's sharing as the structure's sharing.
  base <- spatial_list[[1L]]
  combined <- base
  combined$is_multifield <- TRUE
  combined$fields <- spatial_list
  combined
}

# Build multi-term RE spec for C++ from list of tobs_re objects
build_re_spec <- function(re_list, model) {
  n_terms <- length(re_list)
  N <- if (!is.null(model$N)) model$N else model$n_sites

  groups <- vector("list", n_terms)
  n_groups <- integer(n_terms)
  has_slopes <- FALSE
  n_coefs <- integer(n_terms)
  has_intercept <- rep(TRUE, n_terms)
  slope_matrices <- vector("list", n_terms)

  # Aggregate shared across all terms (union)
  max_proc <- length(model$process_info)
  shared <- rep(FALSE, max_proc)

  for (t in seq_len(n_terms)) {
    re <- re_list[[t]]

    # Group codes were resolved when the re() term was constructed.
    grp <- as.integer(re$group_idx)
    if (length(grp) != N) {
      stop(sprintf("RE group vector has %d elements but model has %d observations",
                   length(grp), N))
    }
    groups[[t]] <- grp
    n_groups[t] <- if (!is.null(re$n_groups)) re$n_groups else max(grp)

    # Sharing
    for (k in seq_along(re$shared)) {
      if (re$shared[k]) shared[k] <- TRUE
    }

    # Slopes. The covariate was resolved at construction: a numeric column or
    # a multi-column matrix (a bare symbol / cbind() in the formula), or column
    # names for a direct re() call. Each column is one slope; n_coefs is the
    # slope count plus the implicit intercept unless the block is slope-only.
    if (re$type == "slope" && !is.null(re$covariate)) {
      has_slopes <- TRUE
      Xs <- .tobs_re_slope_matrix(re$covariate, model$data)
      has_intercept[t] <- isTRUE(re$intercept)
      if (ncol(Xs) == 0L) {
        stop("re(): random slope resolved to zero covariate columns.",
             call. = FALSE)
      }
      n_coefs[t] <- (if (has_intercept[t]) 1L else 0L) + ncol(Xs)
      slope_matrices[[t]] <- Xs
    } else {
      n_coefs[t] <- 1L  # Intercept only
    }
  }

  spec <- list(
    n_terms = n_terms,
    groups = groups,
    n_groups = n_groups,
    shared = shared,
    re_has_intercept = as.integer(has_intercept),
    sigma_re_scale = re_list[[1]]$sigma_scale
  )

  if (has_slopes) {
    spec$has_slopes <- TRUE
    spec$n_coefs <- n_coefs
    spec$slope_matrices <- slope_matrices
    # Per-term correlation flag (0/1): a term is correlated only if it asked
    # for it and has more than one coefficient. The engine reads this per term
    # (re_correlated[t] / re_n_chol[t]), so mixed `|` / `||` blocks are honoured.
    spec$correlated <- as.integer(
      vapply(re_list, function(r) isTRUE(r$correlated), logical(1)) &
      (n_coefs > 1L))
  }

  spec
}

# Resolve an re() slope covariate to a numeric [N x n_slopes] design matrix.
# Accepts column name(s) (resolved via model.matrix, intercept dropped), a
# numeric matrix (e.g. cbind(x, z) from bar desugaring), or a single numeric
# vector. Column names are carried through for ranef() labelling.
.tobs_re_slope_matrix <- function(cov, data) {
  if (is.character(cov)) {
    X <- stats::model.matrix(stats::reformulate(cov), data)
    icpt <- match("(Intercept)", colnames(X))
    if (!is.na(icpt)) X <- X[, -icpt, drop = FALSE]
    return(X)
  }
  if (is.matrix(cov)) {
    storage.mode(cov) <- "double"
    if (is.null(colnames(cov))) {
      colnames(cov) <- paste0("slope", seq_len(ncol(cov)))
    }
    return(cov)
  }
  m <- matrix(as.numeric(cov), ncol = 1L)
  colnames(m) <- "slope1"
  m
}

# Compute back-transformed intercepts on the response scale. The link is read
# per process (`pi$link`, default "logit"): logit-link processes (occupancy /
# detection probabilities) back-transform with plogis(); the log-link
# abundance process (lambda) with exp(), giving the mean per-site abundance.
compute_intercepts <- function(model, means) {
  result <- list()
  offset <- 0
  for (pi in model$process_info) {
    b0 <- means[offset + 1]
    link <- pi$link %||% "logit"
    result[[pi$name]] <- if (identical(link, "log")) exp(b0)
                         else if (identical(link, "identity")) b0
                         else plogis(b0)
    offset <- offset + pi$p
  }
  result
}

# Transform a coefficient covariance from the autoscaled design back to the
# natural scale. The per-process scaling is a block-diagonal linear map T (each
# block's `.scale_transform()`), so vcov_natural = T vcov_scaled T'. Cross-arm
# blocks are transformed exactly (block-diagonal T preserves them).
.unscale_vcov <- function(vcov, scales, process_info) {
  if (is.null(vcov) || is.null(scales) || is.null(process_info)) return(vcov)
  p_tot <- nrow(vcov)
  Tfull <- diag(p_tot)
  off <- 0L
  for (k in seq_along(process_info)) {
    p_k <- as.integer(process_info[[k]]$p)
    if (p_k == 0L) next
    sc <- scales[[k]]
    if (!is.null(sc) && length(sc$cols) > 0L) {
      idx <- off + seq_len(p_k)
      Tfull[idx, idx] <- .scale_transform(sc)
    }
    off <- off + p_k
  }
  out <- Tfull %*% vcov %*% t(Tfull)
  dimnames(out) <- dimnames(vcov)
  out
}

# Build spatial params list for C++ from spatial spec (or NULL)
build_spatial_params <- function(spatial, n_sites) {
  if (is.null(spatial)) return(list(type = "none"))
  if (isTRUE(spatial$is_multifield) || isTRUE(spatial$is_bar) ||
      !is.null(spatial$weight)) {
    stop("A spatially-varying coefficient (a `spatial(~ ... || node)` bar or a ",
         "weighted areal term) is fitted on the nested-Laplace engine ",
         "(method = \"nested_laplace\"), not the NUTS sampler.", call. = FALSE)
  }

  params <- list(type = spatial$type)

  if (spatial$type %in% c("icar", "bym2")) {
    if (spatial$n_units != n_sites) {
      stop(sprintf("spatial has %d units but model has %d sites",
                   spatial$n_units, n_sites))
    }
    params$n_units <- spatial$n_units
    params$adj_row_ptr <- spatial$adj_row_ptr
    params$adj_col_idx <- spatial$adj_col_idx
    params$n_neighbors <- spatial$n_neighbors
    params$spatial_shared_occ <- spatial$shared[1]
    params$spatial_shared_det <- spatial$shared[2]
    # The compiled sampler is tulpa's, and its `scale_factor` multiplies the
    # structured block (`sigma * sqrt(rho) * scale_factor`), where the term
    # carries the Riebler constant itself --.
    if (spatial$type == "bym2") {
      params$scale_factor <- .bym2_engine_scale(spatial$scale_factor)
    }

  } else if (spatial$type == "gp") {
    if (spatial$n_obs != n_sites) {
      stop(sprintf("spatial has %d locations but model has %d sites",
                   spatial$n_obs, n_sites))
    }
    params$n_obs <- spatial$n_obs
    params$nn <- spatial$nn
    params$coords <- spatial$coords
    params$nn_idx <- spatial$nn_idx
    params$nn_dist <- spatial$nn_dist
    params$nn_neighbor_dist <- spatial$nn_neighbor_dist
    params$nn_order <- spatial$nn_order
    params$nn_order_inv <- spatial$nn_order_inv
    params$cov_type <- spatial$cov_type
    params$nu <- spatial$nu
    params$spatial_shared_occ <- spatial$shared[1]
    params$spatial_shared_det <- spatial$shared[2]
    params$sigma2_prior_U <- spatial$sigma2_prior_U
    params$sigma2_prior_alpha <- spatial$sigma2_prior_alpha
    params$phi_prior_U <- spatial$prior_range[1]
    params$phi_prior_alpha <- spatial$prior_range[2]

  } else if (spatial$type == "multiscale_gp") {
    if (spatial$n_obs != n_sites) {
      stop(sprintf("spatial has %d locations but model has %d sites",
                   spatial$n_obs, n_sites))
    }
    params$n_obs <- spatial$n_obs
    params$coords <- spatial$coords
    params$nn_local <- spatial$nn_local
    params$nn_idx_local <- spatial$nn_idx_local
    params$nn_dist_local <- spatial$nn_dist_local
    params$nn_neighbor_dist_local <- spatial$nn_neighbor_dist_local
    params$nn_order_local <- spatial$nn_order_local
    params$nn_order_inv_local <- spatial$nn_order_inv_local
    params$nn_regional <- spatial$nn_regional
    params$nn_idx_regional <- spatial$nn_idx_regional
    params$nn_dist_regional <- spatial$nn_dist_regional
    params$nn_neighbor_dist_regional <- spatial$nn_neighbor_dist_regional
    params$nn_order_regional <- spatial$nn_order_regional
    params$nn_order_inv_regional <- spatial$nn_order_inv_regional
    params$cov_type <- spatial$cov_type
    params$nu <- spatial$nu
    params$spatial_shared_occ <- spatial$shared[1]
    params$spatial_shared_det <- spatial$shared[2]
    params$range_local_lower <- spatial$range_local_lower
    params$range_local_upper <- spatial$range_local_upper
    params$range_regional_lower <- spatial$range_regional_lower
    params$range_regional_upper <- spatial$range_regional_upper
    params$sigma2_local_prior_U <- spatial$sigma2_local_prior_U
    params$sigma2_local_prior_alpha <- spatial$sigma2_local_prior_alpha
    params$sigma2_regional_prior_U <- spatial$sigma2_regional_prior_U
    params$sigma2_regional_prior_alpha <- spatial$sigma2_regional_prior_alpha
  }

  params
}
