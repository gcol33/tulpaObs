# ============================================================================
# Shared helpers
# ============================================================================

# SPDE (continuous mesh field) coverage on the .tobs_laplace path, by model
# type and arm (each cell is "wired" or an honest stop()):
#
#   model_type   state arm (shared[1])   detection arm (shared[2])
#   single       wired                   wired
#   integrated   wired                   wired (per source)
#   jsdm         wired                   n/a (no detection process)
#   community    wired                   stop()
#   dynamic      wired (psi1 only)       stop()
#
# The state field on jsdm / community broadcasts the site-indexed mesh projection A
# onto the N = n_sites * n_species state rows (one shared site-level field across
# the species at a site, via .tobs_spde_broadcast_spec and .tobs_state_block_dims).
# The dynamic state field enters season-1 psi1 only; the colonization / extinction
# transition predictors are separate latent processes whose own mesh fields are not
# wired (a state-arm spde() term maps to psi1). The integrated detection cell is
# reached by writing the spde() term on the `detection` formula (a term's arm
# membership is the formula it is written in). That arm is S source blocks reading
# one shared latent occupancy state, so each source fits its own realization of the
# field at its own sites and `spatial_field_det` is a per-source named list;
# .tobs_validate_integrated_terms() (R/occu.R) gates what the arm accepts, since
# the areal kinds and the temporal / re / svc / latent classes are built into the
# multi-block latent prior the nested-Laplace path attaches to the STATE block. A
# single field shared across both arms at once (shared = c(TRUE, TRUE)) is a stop()
# everywhere: the single-Laplace block fitter fits one field realization per
# submodel block, so a genuinely shared realization needs the copy() path, not two
# independent blocks. The areal path (icar/bym2/car_proper via nested_laplace) is
# wider; this matrix is the continuous-mesh SPDE path only. The continuous
# gp()/spde() fields on the N-mixture arms (abun / em_nested / ms_abun) are tracked
# separately.

# Validate that `spatial` (a `tobs_spatial` or NULL) can be consumed by the
# Laplace path. Wired by arm and model type:
#   state arm (shared[1]): single, integrated, jsdm, community, dynamic (psi1)
#   detection arm (shared[2]): single, integrated (per source)
# jsdm has no detection process; community and dynamic do not yet carry a
# detection-arm field. A single realization shared across both arms at once
# needs the copy() path (the single-Laplace block fitter fits one realization
# per submodel block), so c(TRUE, TRUE) errors here rather than silently fitting
# two independent fields. Other combinations error explicitly.
.validate_spatial_laplace <- function(spatial, model_type) {
  if (is.null(spatial)) return(invisible())
  if (!inherits(spatial, "tobs_spatial")) {
    stop("spatial must be a tobs_spatial term (from a spde() formula term)",
         call. = FALSE)
  }
  if (!identical(spatial$type, "spde")) {
    stop(sprintf(
      ".tobs_laplace currently supports spatial$type == 'spde' only (got '%s'). Use method = 'nuts' for other spatial types.",
      spatial$type), call. = FALSE)
  }
  on_occ <- isTRUE(spatial$shared[1])
  on_det <- length(spatial$shared) >= 2 && isTRUE(spatial$shared[2])
  if (!on_occ && !on_det) {
    stop("SPDE must be attached to the occupancy/state or detection submodel.",
         call. = FALSE)
  }
  if (on_occ && on_det) {
    stop("A single SPDE field shared across the occupancy and detection arms is not plumbed in .tobs_laplace; attach the field to one arm, or use method = 'nuts'.",
         call. = FALSE)
  }
  if (on_det && model_type %in% c("community", "dynamic")) {
    stop(sprintf(
      "SPDE on the detection process is plumbed for single-season and integrated occupancy only in .tobs_laplace (got model_type = '%s'). Attach the field to the state arm, or use method = 'nuts'.",
      model_type), call. = FALSE)
  }
  invisible()
}

# Gate the deterministic random-effect path. The variance-component EM in
# R/em_laplace_re.R fits iid intercept, uncorrelated slopes, and correlated
# slopes (a full RE covariance) on EITHER the occupancy or the detection
# predictor of a single-season model (each arm carries its own RE block). Forms
# it cannot fit -- non-single families, RE + spatial, RE + visit-level
# detection, a single RE shared across both predictors -- error here with a
# pointer to `method = "nuts"` (which fits every RE form) rather than being
# silently dropped. The raw EM variance components (sigma, correlation) carry
# the Laplace small-cluster bias for binary data (the glmer nAGQ=1 regime, not
# Breslow-Clayton PQL); the default re.aghq = TRUE refines them on the
# exact-marginal adaptive Gauss-Hermite likelihood (R/re_aghq.R), removing the
# attenuation, with a default LKJ(re.lkj = 1.5) penalty regularizing a
# weakly-identified RE correlation off the +-1 boundary.
.validate_re_laplace <- function(re, model, spatial) {
  re_list <- if (inherits(re, "tobs_re")) list(re) else re

  if (!identical(model$model_type, "single")) {
    stop(sprintf(
      "Random effects under method = 'laplace' are wired for single-season occupancy only (got model_type = '%s'). Use method = 'nuts' for random effects on this family.",
      model$model_type), call. = FALSE)
  }
  if (!is.null(spatial)) {
    stop("A random effect combined with a spatial term is not supported on the Laplace path. Use method = 'nuts'.",
         call. = FALSE)
  }
  if (!is.null(model$X_det_visit)) {
    stop("Random effects with visit-level detection covariates are not supported on the Laplace path. Use method = 'nuts'.",
         call. = FALSE)
  }
  for (r in re_list) {
    if (length(r$shared) >= 2L && isTRUE(r$shared[1]) && isTRUE(r$shared[2])) {
      stop("A single random effect shared across occupancy and detection is not supported on the Laplace path (each arm fits its own RE block). Use method = 'nuts'.",
           call. = FALSE)
    }
  }
  invisible()
}

# Linear-predictor offset induced by the SPDE mesh field at the current fit.
# After tulpa_laplace returns mode = c(beta, u_mesh), the spatial contribution
# to eta at the observed locations is A %*% u_mesh.
.spatial_eta_offset <- function(spatial, fits_sub, p_fixed) {
  if (is.null(spatial) || is.null(fits_sub) || is.null(fits_sub$mode)) {
    return(rep(0, 0))
  }
  if (!identical(spatial$type, "spde")) return(rep(0, 0))
  mode_vec <- fits_sub$mode
  if (length(mode_vec) <= p_fixed) return(rep(0, 0))
  u <- mode_vec[(p_fixed + 1L):length(mode_vec)]
  as.numeric(spatial$tulpa_spec$A %*% u)
}

# Stop a detection-arm structured term from being fit against the state arm.
#
# The two Laplace routes carry a field differently. The single-Laplace route
# hands each M-step block its own `spatial` spec, so `shared = c(occ, det)`
# decides which block gets the field. The nested-Laplace route builds every
# latent block upstream (`.tobs_to_multi_block_prior`) and attaches the whole
# multi-block prior to the STATE block, with no arm channel -- a detection-arm
# term arriving there is fit against the occupancy predictor, which is not the
# arm it was written on. `latent_prior` is what distinguishes the routes at the
# callbacks builder, so the reject sits there.
.tobs_reject_nested_det_term <- function(model, latent_prior, family) {
  if (is.null(latent_prior)) return(invisible())
  terms <- model$structured_terms
  if (is.null(terms) || !any(vapply(terms, function(t) 2L %in% t$processes,
                                    logical(1)))) {
    return(invisible())
  }
  stop(sprintf(paste0(
    "%s(): a structured term on the `detection` formula is fit by the ",
    "single-Laplace EM (method = \"laplace\"); method = \"nested_laplace\" ",
    "attaches its latent blocks to the state arm only, so the term would be ",
    "fit against occupancy. Re-fit with method = \"laplace\", or move the term ",
    "to the occupancy formula."), family), call. = FALSE)
}

# Latent tail of a fitted M-step block's mode. tulpa_laplace returns
# mode = c(beta (p_fixed), latent units...) for a block carrying an SPDE field
# or a multi-block latent prior. NULL before the first M-step has produced a
# mode, and for a block with no latent tail.
.tobs_block_latent_tail <- function(fit_sub, p_fixed) {
  mode_vec <- fit_sub$mode
  if (is.null(mode_vec) || length(mode_vec) <= p_fixed) return(NULL)
  mode_vec[(p_fixed + 1L):length(mode_vec)]
}

# Attach the tulpa-side spatial spec to an M-step block. The block's `spatial`
# field is forwarded as-is by tulpa_em_laplace -> tulpa_laplace.
.attach_spatial_spde <- function(block, spatial) {
  if (is.null(spatial) || !identical(spatial$type, "spde")) return(block)
  block$spatial <- spatial$tulpa_spec
  block
}

# Encode one aggregated binomial detection block for an EM M-step: `nd`
# detections out of `nv` valid visits per row, each row weighted by the E-step
# occupancy probability `w`. The weight is what makes the M-step maximise the
# expected complete-data log-likelihood -- a row the E-step calls almost
# certainly empty must not feed its all-zero history in as evidence about
# (1 - p)^J.
#
# Row set: without a field, a row carrying no visits or negligible weight
# contributes nothing and is dropped. A mesh field pins the row set instead --
# the block's projection A carries exactly one row per input row (broadcast onto
# the source's sites for an integrated block), so dropping a row would leave the
# remaining rows reading the wrong mesh basis. Those rows enter at weight ~0,
# which is what dropping them expresses.
.tobs_encode_det_block <- function(nd, nv, X, w, keep, spatial) {
  dk <- if (is.null(spatial)) keep & (w > 1e-6) else rep(TRUE, length(w))
  block <- list(y = nd[dk], n_trials = nv[dk], X = X[dk, , drop = FALSE],
                weights = w[dk], family = "binomial")
  .attach_spatial_spde(block, spatial)
}

# Broadcast a site-indexed SPDE field onto a state block whose rows are
# (site, species) -- community and jsdm carry N = n_sites * n_species rows
# ordered site-major, so a single site-level field is shared across the species
# at a site. The mesh / FEM matrices (C, G, n_mesh, nu, priors) describe the
# field on the sites and are unchanged; only the projection A (and its
# pre-extracted CSC slots A_x / A_i / A_p, which `laplace_spde_at()` hands to
# the SPDE solver) is re-rowed so row r of the block projects through the mesh
# basis at `site_of_row[r]`. The eta offset reader `.spatial_eta_offset()` uses
# the same broadcast A, so the field contribution lands on every species row of
# a site identically. `site_of_row` is 1-based into the n_sites mesh rows.
.tobs_spde_broadcast_spec <- function(spatial, site_of_row) {
  if (is.null(spatial) || !identical(spatial$type, "spde")) return(spatial)
  sp <- spatial$tulpa_spec
  A_b <- sp$A[site_of_row, , drop = FALSE]
  A_csc <- methods::as(A_b, "CsparseMatrix")
  sp$A   <- A_b
  sp$A_x <- A_csc@x
  sp$A_i <- A_csc@i
  sp$A_p <- A_csc@p
  spatial$tulpa_spec <- sp
  spatial
}

# Select the SPDE spec for a given arm (1 = occupancy/state, 2 = detection)
# from a `tobs_spatial` term carrying a `$shared = c(occ, det)` membership.
# Returns the spec when the arm carries the field, NULL otherwise. One mesh
# term is shared across the arms it enters, so each arm references the same
# `tulpa_spec` (the field realization is fit independently per arm here -- two
# separate mesh blocks, not a copied realization, which the single-Laplace path
# does not support across submodels).
.spatial_for_arm <- function(spatial, arm) {
  if (is.null(spatial) || !identical(spatial$type, "spde")) return(NULL)
  sh <- spatial$shared
  if (length(sh) >= arm && isTRUE(sh[arm])) spatial else NULL
}

# NUTS sampler-health diagnostics for a non-sampled (Laplace / nested-Laplace)
# fit. No HMC trajectory exists, so acceptance, divergence, tree depth and the
# integrator step size are unavailable; they are NA rather than 0/1 so a user
# inspecting sampler health does not read "no sampler ran" as "sampler ran
# cleanly" (NA-on-unavailable, the same rule as .se_from_laplace_fit). Splice
# into a tobs_fit build with !!! / do.call so the named fields land directly.
.tobs_na_nuts_diagnostics <- function(n_draws) {
  list(
    accept_prob = rep(NA_real_, n_draws),
    divergent   = rep(NA_real_, n_draws),
    treedepth   = rep(NA_integer_, n_draws),
    epsilon     = NA_real_
  )
}

# Maximise a closed-form marginal by BFGS and pack the result as a `tobs_fit`.
#
# The families whose latent state marginalises analytically -- royle_nichols,
# occu_ttd, double_observer, gdistremoval, occu_multi, distsamp_open -- all end
# the same way: minimise the negative log-likelihood, invert the observed
# information at the mode, draw a pseudo-posterior from N(mode, V), and hand the
# S3 layer one slot roster. Only the objective, the parameter names, the sample
# size and a handful of family-specific slots differ.
#
# `gr` is the gradient of `nll`. Supplying it both drives BFGS and changes where
# the observed information comes from: the Jacobian of the negative-log-
# likelihood gradient IS that information, and computing it by finite difference
# of an analytic gradient costs 2p cheap gradient calls, far less than optim's
# numeric Hessian of the value. Without `gr`, optim's own Hessian is used. A
# singular information matrix yields an all-NA covariance rather than an error,
# so the fit still returns with its SEs reading as unavailable.
#
# `extra` is a function of the fitted `means` returning the family-specific
# slots to splice in (an intercept list, a mixture label, a dispersion, ...).
.tobs_bfgs_marginal_fit <- function(nll, init, par_names, model, N,
                                    gr = NULL,
                                    extra = function(means) list(),
                                    control = list(maxit = 500L),
                                    n_draws = 1000L) {
  opt <- stats::optim(init, nll, gr = gr, method = "BFGS",
                      hessian = is.null(gr), control = control)

  means <- opt$par; names(means) <- par_names
  info <- if (is.null(gr)) opt$hessian
          else tryCatch(.tobs_fd_jacobian(gr, opt$par), error = function(e) NULL)
  V <- tryCatch(solve(info), error = function(e) diag(NA_real_, length(means)))
  V <- (V + t(V)) / 2
  dimnames(V) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  draws <- .rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(-opt$value, n_draws),
    log_lik      = -opt$value,
    N            = N),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    n_fixed      = length(means),
    fixed_names  = par_names,
    process_info = model$process_info,
    model        = model,
    spatial      = NULL,
    method       = "laplace",
    convergence  = list(converged = opt$convergence == 0L,
                        n_iter = opt$counts[[1L]])),
    extra(means)),
    class = c("tobs_fit", "tulpa_fit"))
}

extract_beta <- function(sub, p) {
  if (is.null(sub)) return(rep(0, p))
  if (!is.null(sub$beta)) return(sub$beta)
  if (!is.null(sub$mean)) return(sub$mean)
  if (!is.null(sub$mode)) return(sub$mode[seq_len(p)])
  rep(0, p)
}

# SE for the fixed-effect block of a tulpa_laplace() fit. Reads the
# negative-log-posterior Hessian (`H_beta`, the precision matrix), inverts
# it, and returns sqrt(diag(.)) restricted to the first `p` fixed effects.
# When the inner fit had a spatial mesh field attached (`spde` / `gp`),
# tulpa_laplace skips H_beta — return NA so callers can flag the
# uncertainty as unavailable instead of carrying a placeholder.
.se_from_laplace_fit <- function(fi, p) {
  if (!is.null(fi$se)) {
    se <- as.numeric(fi$se)
    if (length(se) >= p) return(se[seq_len(p)])
  }
  H <- fi$H_beta
  if (is.null(H)) return(rep(NA_real_, p))
  cov <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(cov)) return(rep(NA_real_, p))
  d <- sqrt(pmax(diag(cov), 0))
  if (length(d) >= p) return(d[seq_len(p)])
  c(d, rep(NA_real_, p - length(d)))
}

# Louis-corrected observed Fisher info for the occupancy fixed-effect block
# of a single-season occu fit.
#
# Why this is needed. The inner M-step encodes the soft-imputed P(z_i = 1 | y_i)
# as a pseudo-binomial likelihood with n_trials = M (M = 1000 non-spatial,
# M = 4 spatial). The resulting inner Hessian is
#
#   H_inner = M * X' diag(psi (1 - psi)) X + P_prior
#
# i.e. the complete-data Fisher info inflated by the M trick, plus the prior
# precision. This is the wrong object for SE reporting on two counts: the M
# factor is an artefact of the M-step encoding (not data information), and
# the complete-data info ignores the missing-z variance.
#
# Louis identity for the occupancy score s_i = x_i (z_i - psi_i) gives the
# observed Fisher info at the EM stationary point:
#
#   I_obs(beta_psi) = E[-d2 log f / dbeta2 | y] - Var(s_complete | y)
#                   = X' diag(psi (1 - psi)) X - X' diag(w (1 - w)) X
#                   = X' diag(psi (1 - psi) - w (1 - w)) X
#
# where w_i = P(z_i = 1 | y_i, theta_hat) is the converged E-step weight. The
# per-site `psi(1-psi) - w(1-w)` term can be negative (the marginal log-lik
# can be locally convex at a single site), but the aggregate X' D X is PSD at
# the MLE because it equals minus the marginal log-lik Hessian at its max.
#
# The occupancy score x_i (z_i - psi_i) is family-generic, so the same identity
# gives the state-block info of the dynamic (season-1 weights `w[, 1]`) and
# integrated fits. `prior_arms` names the submodel key(s) the prior spec is
# looked up under, in order, taking the first that resolves -- the dynamic fit
# keys its initial-occupancy prior on "psi1", the others on "psi".
.louis_info_psi_single <- function(X_occ, beta_psi, weights,
                                   spatial = NULL, spatial_fit = NULL,
                                   prior_spec = NULL,
                                   coef_names = NULL,
                                   prior_arms = "psi") {
  p_psi <- length(beta_psi)
  if (p_psi == 0L) return(NULL)
  if (is.null(X_occ) || nrow(X_occ) == 0L) return(NULL)
  if (is.null(weights) || length(weights) != nrow(X_occ)) return(NULL)

  eta <- as.numeric(X_occ %*% beta_psi)
  sp_off <- .spatial_eta_offset(spatial, spatial_fit, p_psi)
  if (length(sp_off) == nrow(X_occ)) eta <- eta + sp_off
  eta <- .tobs_clamp_eta(eta)
  psi <- plogis(eta)

  d <- psi * (1 - psi) - weights * (1 - weights)
  I_obs <- as.matrix(crossprod(X_occ, d * X_occ))

  if (!is.null(prior_spec)) {
    if (is.null(coef_names)) coef_names <- colnames(X_occ) %||% paste0("x", seq_len(p_psi))
    pr <- NULL
    for (arm in prior_arms) {
      pr <- .prior_for_submodel(prior_spec, arm, coef_names)
      if (!is.null(pr)) break
    }
    if (!is.null(pr)) {
      pen_prec <- ifelse(is.finite(pr$sd), 1 / (pr$sd^2), 0)
      diag(I_obs) <- diag(I_obs) + pen_prec[seq_len(p_psi)]
    }
  }
  I_obs
}

# Marginal Louis observed Fisher info for the SITE-LEVEL detection block of a
# single-season fit (detection arm). The detection M-step fits a weighted
# binomial whose returned H_beta = X_det' diag(w n_valid p(1-p)) X_det is the
# *complete-data* info: it treats the soft occupancy weight w_i = P(z_i = 1 |
# y) as known and so under-states the SE (the M-step Hessian is the wrong
# object for SEs, exactly as on the psi arm). The occupancy and detection
# estimating equations both depend on the latent z, so the detection SE must
# come from the JOINT (psi, det) Louis observed info, marginalized over psi --
# the diagonal det block alone fixes the intercept but leaves the slope
# under-dispersed.
#
# Complete-data scores: s_psi,i = x_psi,i (z_i - psi_i),
# s_det,i = z_i (n_det_i - n_valid_i p_i) x_det,i. With z_i | y ~ Bern(w_i) the
# Louis identity (E[I_complete | y] - Var(s_complete | y)) gives the joint
# observed info in three blocks (the complete-data cross block is 0):
#
#   I_pp = X_psi' diag( psi(1-psi) - w(1-w) ) X_psi                 (+ psi prior)
#   I_dd = X_det' diag( w n_valid p(1-p) - (n_valid p)^2 w(1-w) ) X_det (+ p prior)
#   I_pd = - X_psi' diag( n_valid p w(1-w) ) X_det
#
# (a detected site has w_i = 1 so its w(1-w) terms vanish.) The marginal
# detection info is the Schur complement I_dd - I_pd' I_pp^{-1} I_pd, whose
# inverse is the (beta_det) block of the full joint covariance.
.louis_info_det_single <- function(X_occ, beta_psi, X_det, beta_det,
                                   weights, n_valid, prior_spec = NULL,
                                   occ_coef_names = NULL, det_coef_names = NULL,
                                   spatial = NULL, spatial_fit = NULL) {
  p_det <- length(beta_det)
  p_psi <- length(beta_psi)
  if (p_det == 0L) return(NULL)
  if (is.null(X_det) || nrow(X_det) == 0L) return(NULL)
  if (is.null(weights) || length(weights) != nrow(X_det)) return(NULL)
  if (is.null(n_valid) || length(n_valid) != nrow(X_det)) return(NULL)

  w  <- weights
  nv <- as.numeric(n_valid)
  p  <- plogis(.tobs_clamp_eta(as.numeric(X_det %*% beta_det)))

  add_prior <- function(I, arm, p_k, coef_names, Xcols) {
    if (is.null(prior_spec)) return(I)
    if (is.null(coef_names)) coef_names <- Xcols %||% paste0("x", seq_len(p_k))
    pr <- .prior_for_submodel(prior_spec, arm, coef_names)
    if (!is.null(pr)) {
      pen <- ifelse(is.finite(pr$sd), 1 / (pr$sd^2), 0)
      diag(I) <- diag(I) + pen[seq_len(p_k)]
    }
    I
  }

  # Detection diagonal block.
  d_dd <- w * nv * p * (1 - p) - (nv * p)^2 * w * (1 - w)
  d_dd[nv <= 0] <- 0
  I_dd <- add_prior(as.matrix(crossprod(X_det, d_dd * X_det)),
                    "p", p_det, det_coef_names, colnames(X_det))

  # Couple with the occupancy block via the joint Louis cross term, then
  # marginalize psi out by Schur complement. Skip when the occupancy inputs are
  # unavailable (fall back to the diagonal block, which still fixes the level).
  if (p_psi > 0L && !is.null(X_occ) && nrow(X_occ) == nrow(X_det)) {
    eta_o <- as.numeric(X_occ %*% beta_psi)
    sp_off <- .spatial_eta_offset(spatial, spatial_fit, p_psi)
    if (length(sp_off) == nrow(X_occ)) eta_o <- eta_o + sp_off
    psi <- plogis(.tobs_clamp_eta(eta_o))

    d_pp <- psi * (1 - psi) - w * (1 - w)
    I_pp <- add_prior(as.matrix(crossprod(X_occ, d_pp * X_occ)),
                      "psi", p_psi, occ_coef_names, colnames(X_occ))
    d_pd <- -nv * p * w * (1 - w)
    d_pd[nv <= 0] <- 0
    I_pd <- as.matrix(crossprod(X_occ, d_pd * X_det))    # p_psi x p_det

    schur <- tryCatch(I_dd - crossprod(I_pd, solve(I_pp, I_pd)),
                      error = function(e) NULL)
    if (!is.null(schur)) I_dd <- schur
  }
  I_dd
}

# SE vector from an observed-info matrix; returns NA of length p on failure.
.se_from_info <- function(I, p) {
  if (is.null(I)) return(rep(NA_real_, p))
  cov <- tryCatch(solve(I), error = function(e) NULL)
  if (is.null(cov)) return(rep(NA_real_, p))
  d <- sqrt(pmax(diag(cov), 0))
  if (length(d) >= p) d[seq_len(p)] else c(d, rep(NA_real_, p - length(d)))
}

# Within-arm covariance block carrying the off-diagonal correlation of the
# observed-information inverse, rescaled so its diagonal matches the reported
# marginal SEs exactly (`sds`). `prec` is the precision / observed-information
# matrix the SEs were derived from (cov = solve(prec)); NULL or a non-invertible
# `prec` yields the diagonal block diag(sds^2), i.e. the previous behaviour. This
# keeps the marginal SEs byte-identical while restoring the joint correlation the
# diagonal pseudo-draws used to discard. NA / non-finite SEs map to a
# zero-variance coordinate (drawn as a point mass downstream).
.cor_scaled_cov <- function(prec, sds) {
  p <- length(sds)
  s <- ifelse(is.finite(sds), sds, 0)
  if (is.null(prec)) return(diag(s^2, nrow = p))
  cov <- tryCatch(solve(prec), error = function(e) NULL)
  if (is.null(cov)) return(diag(s^2, nrow = p))
  cov <- as.matrix(cov)[seq_len(p), seq_len(p), drop = FALSE]
  d <- sqrt(pmax(diag(cov), 0))
  ok <- is.finite(d) & d > 0
  R <- diag(p)
  if (sum(ok) > 1L) R[ok, ok] <- cov[ok, ok, drop = FALSE] / tcrossprod(d[ok])
  outer(s, s) * R
}

# Assemble a block-diagonal covariance from per-block matrices in append order,
# flooring zero / NA variances so the matrix is PD and chol-decomposable in
# .rmvn (mirrors the old `max(sd_j, 1e-4)` point-mass floor for NA-SE columns).
.assemble_block_diag <- function(blocks, n_params) {
  V <- matrix(0, n_params, n_params)
  off <- 0L
  for (B in blocks) {
    b <- nrow(B)
    if (b == 0L) next
    idx <- off + seq_len(b)
    V[idx, idx] <- B
    off <- off + b
  }
  dd <- diag(V)
  dd[!is.finite(dd) | dd <= 0] <- 1e-8
  diag(V) <- dd
  V
}

clamp_w <- function(w) pmin(pmax(w, 0.001), 0.999)

occ_weights <- function(psi, p, N, n_valid, n_det, any_det) {
  weights <- numeric(N)
  for (i in seq_len(N)) {
    if (any_det[i]) { weights[i] <- 1 }
    else if (n_valid[i] == 0) { weights[i] <- psi[i] }
    else {
      prod_1mp <- (1 - p[i])^n_valid[i]
      num <- psi[i] * prod_1mp
      weights[i] <- num / (num + (1 - psi[i]))
    }
  }
  weights
}

glm_init <- function(X_occ, X_det, any_det, n_det, n_valid, keep, p_occ, p_det) {
  tryCatch({
    occ_glm <- glm(any_det ~ X_occ[, -1, drop = FALSE] - 1 + X_occ[, 1], family = binomial)
    det_glm <- glm(cbind(n_det[keep], n_valid[keep] - n_det[keep]) ~
                      X_det[keep, -1, drop = FALSE] - 1 + X_det[keep, 1], family = binomial)
    list(occ = list(beta = unname(coef(occ_glm)), se = rep(1, p_occ)),
         det = list(beta = unname(coef(det_glm)), se = rep(1, p_det)))
  }, error = function(e) {
    list(occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)),
         det = list(beta = rep(0, p_det), se = rep(1, p_det)))
  })
}

# Build tobs_fit from EM result
build_laplace_fit <- function(em_result, model, spatial, p_per_submodel,
                              prior_spec = NULL,
                              approx = "gaussian_laplace",
                              re_block = NULL, latent_prior = NULL) {
  pi_list <- model$process_info

  # Per-arm SPDE membership: the field may sit on the state arm, the detection
  # arm, or both. The fixed-effect SE machinery is arm-specific (the Louis
  # observed-info correction assumes a non-spatial M-step Hessian on that arm),
  # so route each arm with its own field.
  spatial_occ <- .spatial_for_arm(spatial, 1L)
  spatial_det <- .spatial_for_arm(spatial, 2L)

  # Collect betas from correction (if available) or EM fits. Each arm also
  # contributes its within-arm covariance block (`prec_k` = the precision the
  # arm's SEs came from), so the pseudo-draws carry the joint correlation rather
  # than a diagonal stand-in. Cross-arm covariance stays zero, matching the EM's
  # separate-arm M-step factorization.
  means <- numeric()
  sds <- numeric()
  nms <- character()
  louis_psi_se <- NULL
  cov_blocks <- list()

  for (k in seq_along(pi_list)) {
    pi <- pi_list[[k]]
    sub_name <- names(p_per_submodel)[k]
    if (is.null(sub_name)) sub_name <- names(p_per_submodel)[min(k, length(p_per_submodel))]
    prec_k <- NULL

    if (is.list(em_result$pooled) && !is.null(em_result$pooled[[sub_name]])) {
      # MI/Gibbs correction pool from rubins_pool().
      cr <- em_result$pooled[[sub_name]]
      means <- c(means, cr$mean)
      sds_k <- cr$se
      cov_k <- if (!is.null(cr$vcov) &&
                   all(dim(as.matrix(cr$vcov)) == length(sds_k)))
                 as.matrix(cr$vcov) else .cor_scaled_cov(NULL, sds_k)
    } else if (!is.null(em_result$fits[[sub_name]])) {
      fi <- em_result$fits[[sub_name]]
      beta <- extract_beta(fi, pi$p)
      means <- c(means, beta)

      # Louis-corrected SE on the psi block of a single-season fit. The inner
      # M-step Hessian is M * I_complete + P_prior (pseudo-binomial trick); I_obs
      # = X' diag(psi(1-psi) - w(1-w)) X + P_prior is the right object for SEs.
      # See `.louis_info_psi_single`. The Louis-corrected occupancy SE applies
      # to the
      # fixed-effect-only fit. When random effects are present the occupancy
      # block's fixed-effect SE comes from the GLMM marginal precision (`H_beta`,
      # Schur over the RE block) that tulpa_laplace returns, so skip Louis on the
      # RE path. The Louis identity for the state arm is not
      # single-season-specific. The arm's complete-data score is x_i (z_i -
      # psi_i) in every one of these families, so I_obs = X' diag(psi(1-psi) -
      # w(1-w)) X with w = E[z_i | y]. For a dynamic fit the state arm is psi1
      # and its latent is z_1, whose smoothed posterior mean is the first
      # season's weight column; integrated shares one psi across sources and
      # carries a single weight per site. A state arm left out of this list falls
      # through to .se_from_laplace_fit(), which finds no H_beta on a
      # nested-Laplace fit and returns NA, so that family reports no state SEs
      # and no intervals.
      use_louis <- model$model_type %in% c("single", "dynamic", "integrated") &&
                   identical(sub_name, "occ") &&
                   !is.null(em_result$weights) &&
                   is.null(re_block)
      if (!is.null(re_block) && identical(sub_name, "occ")) {
        # Occupancy fixed-effect SE on the RE path: natural-scale observed info
        # marginalised over the random-effect block (the M-step H_beta is
        # M-inflated). Computed in .tobs_re_occ_fixed_se().
        sds_k <- re_block$occ_se
      } else if (use_louis) {
        # Dynamic carries an [n_sites x n_seasons] weight matrix; the state arm
        # needs the season-1 column, and the helper length-checks against
        # nrow(X_occ).
        w_occ <- if (identical(model$model_type, "dynamic"))
                   em_result$weights[, 1L] else em_result$weights
        I_obs <- .louis_info_psi_single(
          X_occ       = model$X_processes[[1]],
          beta_psi    = beta,
          weights     = w_occ,
          spatial     = spatial_occ,
          spatial_fit = fi,
          prior_spec  = prior_spec,
          coef_names  = pi$coef_names
        )
        louis_psi_se <- .se_from_info(I_obs, pi$p)
        sds_k <- louis_psi_se
        prec_k <- I_obs
      } else if (identical(model$model_type, "single") &&
                 identical(sub_name, "det") &&
                 is.null(spatial_det) &&
                 is.null(re_block) &&
                 is.null(model$X_det_visit) &&
                 !is.null(em_result$weights) &&
                 nrow(model$X_processes[[2]]) == length(em_result$weights)) {
        # Site-level detection SE via the marginal Louis observed info
        # (detection arm). The detection M-step's H_beta is the
        # complete-data info (soft occupancy weight treated as known), which
        # under-states the SE the same way the psi arm did; recompute the
        # observed info, marginalizing over the coupled occupancy block.
        beta_psi_fit <- extract_beta(em_result$fits[["occ"]],
                                     ncol(model$X_processes[[1]]))
        I_obs <- .louis_info_det_single(
          X_occ          = model$X_processes[[1]],
          beta_psi       = beta_psi_fit,
          X_det          = model$X_processes[[2]],
          beta_det       = beta,
          weights        = em_result$weights,
          n_valid        = rowSums(model$y >= 0),
          prior_spec     = prior_spec,
          occ_coef_names = pi_list[[1]]$coef_names,
          det_coef_names = pi$coef_names,
          spatial        = spatial_occ,
          spatial_fit    = em_result$fits[["occ"]]
        )
        se_det <- .se_from_info(I_obs, pi$p)
        if (any(!is.finite(se_det))) {
          sds_k <- .se_from_laplace_fit(fi, pi$p)
          prec_k <- fi$H_beta
        } else {
          sds_k <- se_det
          prec_k <- I_obs
        }
      } else {
        sds_k <- .se_from_laplace_fit(fi, pi$p)
        prec_k <- fi$H_beta
      }
      cov_k <- .cor_scaled_cov(prec_k, sds_k)
    } else {
      means <- c(means, rep(0, pi$p))
      sds_k <- rep(NA_real_, pi$p)
      cov_k <- .cor_scaled_cov(NULL, sds_k)
    }
    sds <- c(sds, sds_k)
    cov_blocks[[length(cov_blocks) + 1L]] <- cov_k
    nms <- c(nms, paste0(pi$name, "_", pi$coef_names))
  }

  # Append visit-level detection coefficients when X_det_visit is present.
  # The detection M-step block has X of width p_det + p_det_visit; the main
  # loop above extracts only the first p_det elements (the site-level
  # detection coefs). Pull the visit-level tail and label as `p_visit_<name>`
  # so the public output matches the NUTS engine's column layout.
  if (!is.null(model$det_visit_names) && length(model$det_visit_names) > 0L) {
    p_det_visit <- length(model$det_visit_names)
    pi_p <- pi_list[[2]]  # detection process metadata
    p_det <- pi_p$p
    p_det_total <- as.integer(p_per_submodel[["det"]] %||% (p_det + p_det_visit))
    visit_idx <- (p_det + 1L):p_det_total
    visit_nms <- paste0("p_visit_", model$det_visit_names)

    if (is.list(em_result$pooled) && !is.null(em_result$pooled[["det"]])) {
      cr <- em_result$pooled[["det"]]
      visit_means <- cr$mean[visit_idx]
      visit_sds   <- cr$se[visit_idx]
    } else if (!is.null(em_result$fits[["det"]])) {
      fi_det <- em_result$fits[["det"]]
      beta_full <- extract_beta(fi_det, p_det_total)
      se_full <- .se_from_laplace_fit(fi_det, p_det_total)
      visit_means <- beta_full[visit_idx]
      visit_sds   <- se_full[visit_idx]
    } else {
      visit_means <- rep(0, p_det_visit)
      visit_sds   <- rep(NA_real_, p_det_visit)
    }

    means <- c(means, visit_means)
    sds   <- c(sds, visit_sds)
    nms   <- c(nms, visit_nms)
    # Visit-level detection coefs are carried diagonal (their cross-covariance
    # with the site-level det block is not surfaced separately from the M-step).
    cov_blocks[[length(cov_blocks) + 1L]] <- .cor_scaled_cov(NULL, visit_sds)
  }

  # Append the deterministic random-effect block (sigma hyperparameters +
  # per-group BLUPs) so the public output matches the NUTS column layout and
  # ranef() / summary() can name them.
  if (!is.null(re_block)) {
    means <- c(means, re_block$means)
    sds   <- c(sds, re_block$sds)
    nms   <- c(nms, re_block$names)
    cov_blocks[[length(cov_blocks) + 1L]] <- .cor_scaled_cov(NULL, re_block$sds)
  }

  names(means) <- nms
  names(sds)   <- nms
  n_params <- length(means)

  # Pseudo-draws from the block-diagonal joint covariance: full within each
  # fixed-effect arm (so derived quantities like predicted psi = plogis(X beta)
  # propagate the coefficient correlation), zero across arms. Coordinates with
  # an unavailable SE (NA) floor to a near-constant point mass, the same as the
  # previous per-coefficient draw.
  n_pseudo <- 1000L
  V_draw <- .assemble_block_diag(cov_blocks, n_params)
  dimnames(V_draw) <- list(nms, nms)
  draws <- .rmvn(n_pseudo, means, V_draw)
  colnames(draws) <- nms

  # Simplified-Laplace skewness correction
  # Computes gamma_j at the original observation likelihood (NOT the M-step
  # pseudo-binomial encoding — see dev_notes/simplified_laplace_derivation.md
  # §3 and dev_notes/upstream_tulpa_sla_spec.md §3 for why).
  sla_gamma <- NULL
  sla_status <- "off"
  # Simplified-Laplace skewness correction is not wired for the random-effect
  # path (the gamma derivation assumes a fixed-effect-only M-step).
  if (identical(approx, "simplified_laplace") && is.null(re_block)) {
    sla_res <- switch(model$model_type,
      single     = .sla_compute_occu_single(model, em_result,
                                            spatial = spatial,
                                            prior_spec = prior_spec),
      dynamic    = .sla_compute_dyn_occu(model, em_result,
                                         spatial = spatial,
                                         prior_spec = prior_spec),
      integrated = .sla_compute_int_occu(model, em_result,
                                         spatial = spatial,
                                         prior_spec = prior_spec),
      list(gamma = NULL, valid = FALSE,
           reason = sprintf("simplified Laplace not yet supported for model_type '%s'",
                            model$model_type))
    )
    if (isTRUE(sla_res$valid)) {
      sla_gamma  <- sla_res$gamma
      sla_status <- "simplified_laplace"
      # Align gamma names with the joint parameter ordering used in `means`
      sla_gamma <- sla_gamma[intersect(names(means), names(sla_gamma))]
      gamma_full <- setNames(rep(0, n_params), nms)
      gamma_full[names(sla_gamma)] <- sla_gamma
      sla_gamma <- gamma_full
      draws <- .sla_replace_draws(draws, means, sds, sla_gamma)
    } else {
      sla_status <- paste0("fallback_gaussian (", sla_res$reason, ")")
    }
  } else if (identical(approx, "simplified_laplace")) {
    # A random-effect fit asked for the correction and does not get it: the
    # gamma derivation assumes a fixed-effect-only M-step. Record it as a
    # decline, not as "off" -- "off" is the status print() suppresses, and it
    # says the correction was never requested.
    sla_status <- "fallback_gaussian (random-effect M-step)"
  }

  intercepts <- compute_intercepts(model, means)

  # When SPDE is attached to the occ submodel, the M-step mode is
  # c(beta_occ, u_mesh). Extract u_mesh so callers can inspect or project
  # the latent field to observation locations via A %*% u_mesh.
  # Nested-Laplace: the occ block mode is c(beta_occ, latent units...) where the
  # latent tail is the grid-weighted posterior mean of the multi-block field.
  spatial_field <- if (!is.null(spatial_occ) || !is.null(latent_prior))
    .tobs_block_latent_tail(em_result$fits$occ, pi_list[[1]]$p) else NULL

  # When SPDE is attached to the detection submodel, the det M-step mode is
  # c(beta_det, u_mesh_det); extract the detection field tail. The field is
  # identified off its own proper Matern (range, sigma) PC prior the same way
  # the state field is -- no separate sum-to-zero constraint is imposed (the
  # detection intercept absorbs the field level under the mean-zero prior).
  #
  # An integrated model's detection arm is S source blocks reading one shared
  # latent occupancy state, each fitting its own realization of the field at its
  # own sites, so the slot is a per-source named list there (source names as
  # given on `y`); a single-season fit keeps the plain n_mesh vector.
  spatial_field_det <- NULL
  if (!is.null(spatial_det)) {
    is_int <- identical(model$model_type, "integrated")
    det_blocks <- if (is_int) paste0("det", seq_len(model$n_sources)) else "det"
    tails <- lapply(seq_along(det_blocks), function(s)
      .tobs_block_latent_tail(em_result$fits[[det_blocks[s]]],
                              pi_list[[s + 1L]]$p))
    names(tails) <- vapply(seq_along(det_blocks),
                           function(s) pi_list[[s + 1L]]$name, character(1))
    spatial_field_det <- if (!is_int) tails[[1L]]
      else if (all(vapply(tails, is.null, logical(1)))) NULL else tails
  }

  # Marginal log-likelihood at the fixed-effect mode, so logLik() / AIC() /
  # BIC() / glance() surface a finite value. The EM tracks parameter deltas, not
  # the marginal, so evaluate it here through the shared family pointwise kernel
  # -- the same marginal the WAIC / LOO scoring uses. The single-season
  # exact-marginal refine (.tobs_occu_marginal_refine) moves the mode afterwards
  # and refreshes log_lik / log_prob from the refined means.
  marg_ll <- .tobs_laplace_marginal_loglik(model, means)
  # Free-parameter count for the AIC / BIC penalty (logLik()'s `df`): the fixed
  # coefficients (the process betas plus any visit-level detection betas), which
  # is `means` without the trailing random-effect sigma + BLUP block. Set
  # explicitly because an empty `mode` would otherwise resolve df to 0.
  n_re_block <- if (!is.null(re_block)) length(re_block$means) else 0L
  structure(c(list(
    draws = draws, means = means, sds = sds,
    skew = sla_gamma, sla_status = sla_status,
    n_samples = n_pseudo, n_params = n_params,
    n_fixed  = n_params - n_re_block,
    log_prob = rep(marg_ll$loglik, n_pseudo),
    log_lik  = marg_ll$loglik,
    N        = marg_ll$nobs,
    converged = em_result$convergence$converged %||% em_result$converged),
    .tobs_na_nuts_diagnostics(n_pseudo),
    list(
    col_names = nms, param_names = nms,
    intercepts = intercepts,
    model = model, spatial = spatial,
    spatial_field = spatial_field,
    spatial_field_det = spatial_field_det,
    process_info = model$process_info,
    method = "laplace",
    re_effects = re_block$re_effects,
    aghq = em_result$aghq,
    convergence = em_result$convergence,
    correction = em_result$correction
  )), class = c("tobs_fit", "tulpa_fit"))
}
