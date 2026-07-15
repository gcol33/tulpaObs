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
# The state field on jsdm / community broadcasts the site-indexed mesh
# projection A onto the N = n_sites * n_species state rows (one shared
# site-level field across the species at a site, via .tobs_spde_broadcast_spec
# and .tobs_state_block_dims). The dynamic state field enters season-1 psi1
# only; the colonization / extinction transition predictors are separate latent
# processes whose own mesh fields are not wired (a state-arm spde() term maps to
# psi1). A single field shared across both arms at once (shared = c(TRUE, TRUE))
# is a stop() everywhere: the single-Laplace block fitter fits one field
# realization per submodel block, so a genuinely shared realization needs the
# copy() path, not two independent blocks. The areal path (icar/bym2/car_proper
# via nested_laplace) is wider; this matrix is the continuous-mesh SPDE path
# only. The continuous gp()/spde() fields on the N-mixture arms (abun /
# em_nested / ms_abun) are tracked separately (gcol33/tulpaObs#21).

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
# silently dropped (gcol33/tulpaObs#11). The raw EM variance components (sigma,
# correlation) carry the Laplace small-cluster bias for binary data (the glmer
# nAGQ=1 regime, not Breslow-Clayton PQL); the default re.aghq = TRUE refines
# them on the exact-marginal adaptive Gauss-Hermite likelihood (R/re_aghq.R),
# removing the attenuation, with a default LKJ(re.lkj = 1.5) penalty
# regularizing a weakly-identified RE correlation off the +-1 boundary.
.validate_re_laplace <- function(re, model, spatial, approx) {
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

# Attach the tulpa-side spatial spec to an M-step block. The block's `spatial`
# field is forwarded as-is by tulpa_em_laplace -> tulpa_laplace.
.attach_spatial_spde <- function(block, spatial) {
  if (is.null(spatial) || !identical(spatial$type, "spde")) return(block)
  block$spatial <- spatial$tulpa_spec
  block
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

# TRUE when a nested-Laplace multi-block latent prior carries a continuous
# Matern (SPDE) field block. Used to switch the occupancy M-step to the modest
# pseudo-binomial inflation that keeps the mesh-field prior from being swamped.
.tobs_latent_prior_has_spde <- function(latent_prior) {
  if (is.null(latent_prior)) return(FALSE)
  blocks <- if (!is.null(latent_prior$type)) list(latent_prior) else latent_prior
  any(vapply(blocks, function(b) identical(b$type, "spde"), logical(1)))
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
# of a single-season occu fit (tulpaObs#7).
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
.louis_info_psi_single <- function(X_occ, beta_psi, weights,
                                   spatial = NULL, spatial_fit = NULL,
                                   prior_spec = NULL,
                                   coef_names = NULL) {
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
    pr <- .prior_for_submodel(prior_spec, "psi", coef_names)
    if (!is.null(pr)) {
      pen_prec <- ifelse(is.finite(pr$sd), 1 / (pr$sd^2), 0)
      diag(I_obs) <- diag(I_obs) + pen_prec[seq_len(p_psi)]
    }
  }
  I_obs
}

# Marginal Louis observed Fisher info for the SITE-LEVEL detection block of a
# single-season fit (tulpaObs#7, detection arm). The detection M-step fits a
# weighted binomial whose returned H_beta = X_det' diag(w n_valid p(1-p)) X_det
# is the *complete-data* info: it treats the soft occupancy weight
# w_i = P(z_i = 1 | y) as known and so under-states the SE (the M-step Hessian
# is the wrong object for SEs, exactly as on the psi arm). The occupancy and
# detection estimating equations both depend on the latent z, so the detection
# SE must come from the JOINT (psi, det) Louis observed info, marginalized over
# psi -- the diagonal det block alone fixes the intercept but leaves the slope
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
# diagonal pseudo-draws used to discard (gcol33/tulpaObs#44). NA / non-finite SEs
# map to a zero-variance coordinate (drawn as a point mass downstream).
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
  # than a diagonal stand-in (gcol33/tulpaObs#44). Cross-arm covariance stays
  # zero, matching the EM's separate-arm M-step factorization.
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
      # M-step Hessian is M * I_complete + P_prior (pseudo-binomial trick);
      # I_obs = X' diag(psi(1-psi) - w(1-w)) X + P_prior is the right object
      # for SEs. See `.louis_info_psi_single` and tulpaObs#7.
      # Louis-corrected occupancy SE applies to the fixed-effect-only fit. When
      # random effects are present the occupancy block's fixed-effect SE comes
      # from the GLMM marginal precision (`H_beta`, Schur over the RE block)
      # that tulpa_laplace returns, so skip Louis on the RE path.
      use_louis <- identical(model$model_type, "single") &&
                   identical(sub_name, "occ") &&
                   !is.null(em_result$weights) &&
                   is.null(re_block)
      if (!is.null(re_block) && identical(sub_name, "occ")) {
        # Occupancy fixed-effect SE on the RE path: natural-scale observed info
        # marginalised over the random-effect block (the M-step H_beta is
        # M-inflated). Computed in .tobs_re_occ_fixed_se().
        sds_k <- re_block$occ_se
      } else if (use_louis) {
        I_obs <- .louis_info_psi_single(
          X_occ       = model$X_processes[[1]],
          beta_psi    = beta,
          weights     = em_result$weights,
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
        # (tulpaObs#7, detection arm). The detection M-step's H_beta is the
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
  # ranef() / summary() can name them (gcol33/tulpaObs#11).
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
  # previous per-coefficient draw (gcol33/tulpaObs#44).
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
  }

  intercepts <- compute_intercepts(model, means)

  # When SPDE is attached to the occ submodel, the M-step mode is
  # c(beta_occ, u_mesh). Extract u_mesh so callers can inspect or project
  # the latent field to observation locations via A %*% u_mesh.
  spatial_field <- NULL
  if (!is.null(spatial_occ) && !is.null(em_result$fits$occ$mode)) {
    p_occ <- pi_list[[1]]$p
    mode_vec <- em_result$fits$occ$mode
    if (length(mode_vec) > p_occ) {
      spatial_field <- mode_vec[(p_occ + 1L):length(mode_vec)]
    }
  }
  # Nested-Laplace: the occ block mode is c(beta_occ, latent units...) where the
  # latent tail is the grid-weighted posterior mean of the multi-block field.
  if (is.null(spatial_field) && !is.null(latent_prior) &&
      !is.null(em_result$fits$occ$mode)) {
    p_occ <- pi_list[[1]]$p
    mode_vec <- em_result$fits$occ$mode
    if (length(mode_vec) > p_occ) {
      spatial_field <- mode_vec[(p_occ + 1L):length(mode_vec)]
    }
  }

  # When SPDE is attached to the detection submodel, the det M-step mode is
  # c(beta_det, u_mesh_det); extract the detection field tail. The field is
  # identified off its own proper Matern (range, sigma) PC prior the same way
  # the state field is -- no separate sum-to-zero constraint is imposed (the
  # detection intercept absorbs the field level under the mean-zero prior).
  spatial_field_det <- NULL
  if (!is.null(spatial_det) && !is.null(em_result$fits$det$mode)) {
    p_det <- pi_list[[2]]$p
    mode_det <- em_result$fits$det$mode
    if (length(mode_det) > p_det) {
      spatial_field_det <- mode_det[(p_det + 1L):length(mode_det)]
    }
  }

  # Marginal log-likelihood at the fixed-effect mode, so logLik() / AIC() /
  # BIC() / glance() surface a finite value (gcol33/tulpaObs#87). The EM tracks
  # parameter deltas, not the marginal, so evaluate it here through the shared
  # family pointwise kernel -- the same marginal the WAIC / LOO scoring uses.
  # The single-season exact-marginal refine (.tobs_occu_marginal_refine) moves
  # the mode afterwards and refreshes log_lik / log_prob from the refined means.
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
