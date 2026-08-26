# =============================================================================
# field_offset.R - the fitted latent field as a per-arm eta offset
#
# A fit carrying an areal / temporal / continuous field estimates a latent
# surface that is PART OF THE LINEAR PREDICTOR: eta_k = X_k beta_k + f. The
# coefficient means and draws do not carry it (the field is reported separately,
# on `fit$spatial_field` / `fit$temporal_field` / `fit$svc_field`), so every
# post-fit reader that rebuilds eta from `model$X_processes` -- fitted(),
# residuals(), predict(), simulate(), and the pointwise log-likelihood behind
# WAIC / LOO / DIC / CPO -- has to add it back, or it reports and scores a model
# the fit never made.
#
# The convention: `model$field_eta_offset` is a list with one entry per process,
# each NULL or a numeric vector in THAT process's own design row layout
# (`nrow(model$X_processes[[k]])`). A fitter records the field's per-unit offset
# on the fit (`.tobs_set_field_eta_offset()`); `.tobs_attach_model_eta_offset()`
# expands it to the row layout once, at the dispatch tail where
# `.tobs_fit_model()` has swapped the natural-scale model back in. The field is
# latent, so the offset is invariant under the per-process design autoscaling.
# Readers go through `.tobs_eta_offset()` (a vector predictor) or
# `.tobs_add_eta_offset()` (a [n_draws x n_rows] predictor matrix).
#
# `count()` and the community families reached the same conclusion earlier and
# carry their own named slots (`count_field_offset`, `occu_field_offset`, ...),
# written by `.tobs_latent_attach_field()`; the shared field -> per-site
# contribution reader `.tobs_spatial_field_offset()` below is the one both
# conventions use.
# =============================================================================

# Per-site field contribution to eta: sum_k W[i,k] F_k[u(i)]. The intercept field
# carries weight 1; a varying-coefficient (SVC) field carries its covariate
# column. Returns NULL when the fit carries no reconstructed field.
.tobs_spatial_field_offset <- function(fit, spatial, model) {
  # The NNGP GP field is integrated out on the nested path (no reconstructed
  # per-cell realization), so there is no per-site field offset to add; fitted()
  # / predict() run on the field-integrated fixed effects.
  if (identical(spatial$type, "gp")) return(NULL)
  flds <- c(list(fit$spatial_field), as.list(fit$trend_fields %||% list()))
  flds <- Filter(function(z) !is.null(z) && !all(is.na(z)), flds)
  if (!length(flds)) return(NULL)

  # Continuous SPDE field: the realization lives on the mesh nodes, so the
  # per-site contribution is the barycentric projection A %*% mesh_field (A is
  # n_sites x n_mesh), summed over the intercept + any covariate-weighted fields.
  A <- fit$spatial$tulpa_spec$A
  if (!is.null(A)) {
    off <- numeric(nrow(A))
    for (fk in flds) off <- off + as.numeric(A %*% as.numeric(fk))
    return(off)
  }

  sp_fields <- .tobs_resolve_occu_spatial_fields(spatial, model)
  if (length(sp_fields) < length(flds)) return(NULL)
  n   <- length(as.numeric(flds[[1L]]))
  off <- numeric(n)
  for (k in seq_along(flds)) {
    wk <- sp_fields[[k]]$weight
    w  <- if (is.null(wk)) rep(1, n) else as.numeric(wk)
    if (length(w) != n) return(NULL)
    off <- off + w * as.numeric(flds[[k]])
  }
  off
}

# Record the field's eta offset on the fit, tagged with the process (arm) it
# loads on. `off` is in the field's own unit layout -- per site for an areal
# field mapped one node per site, already mapped through `field_map` for a
# temporal or grouped field. A NULL / empty offset leaves the fit untouched.
.tobs_set_field_eta_offset <- function(fit, arm, off) {
  if (is.null(off) || !length(off)) return(fit)
  fit$field_eta_offset <- as.numeric(off)
  fit$field_eta_arm    <- as.integer(arm)
  fit
}

# Expand a per-site offset onto process k's design rows. The state arms and the
# per-site detection arms are one row per site; removal's detection design is one
# row per PASS, reached through the model's site index (the same expansion its
# fitter applies inside the likelihood).
.tobs_expand_site_offset <- function(model, k, off) {
  off    <- as.numeric(off)
  n_rows <- nrow(model$X_processes[[k]])
  if (length(off) == n_rows) return(off)
  si <- model$site_idx
  if (!is.null(si) && length(si) == n_rows &&
      length(off) == (model$n_sites %||% -1L)) return(off[si])
  stop(sprintf(paste0("Field eta offset has length %d, but process %d's design ",
                      "has %d rows and no site index maps between them."),
               length(off), k, n_rows), call. = FALSE)
}

# Move the recorded offset onto the fit's model, in the arm's row layout. Called
# at the dispatch tail, AFTER `fit$model <- model` swaps the natural-scale model
# in (which drops any slot the fitter set on its own autoscaled copy).
.tobs_attach_model_eta_offset <- function(fit) {
  off <- fit$field_eta_offset
  if (is.null(off)) return(fit)
  k     <- fit$field_eta_arm
  slots <- vector("list", length(fit$model$process_info))
  slots[[k]] <- .tobs_expand_site_offset(fit$model, k, off)
  fit$model$field_eta_offset <- slots
  fit
}

# The recorded offset for process k, or NULL when the fit carries no field there.
.tobs_eta_offset <- function(model, k) {
  off <- model$field_eta_offset
  if (is.null(off) || length(off) < k) return(NULL)
  off[[k]]
}

# `eta` is [n_draws x n_rows] for process k; add the field offset broadcast
# across draws (its posterior draws are not sampled on the deterministic paths,
# and the sampled paths report the posterior-mean surface).
.tobs_add_eta_offset <- function(eta, model, k) {
  off <- .tobs_eta_offset(model, k)
  if (is.null(off)) return(eta)
  if (length(off) != ncol(eta)) {
    stop(sprintf(paste0("Field eta offset for process %d has length %d but the ",
                        "predictor has %d columns."), k, length(off), ncol(eta)),
         call. = FALSE)
  }
  eta + matrix(off, nrow(eta), ncol(eta), byrow = TRUE)
}

# The simulation kernels take each arm's design plus one contiguous leading block
# of draw columns, and build eta by row. A fitted field is carried into them as a
# design COLUMN whose coefficient is pinned at 1 -- eta = X beta + f is
# [X | f] %*% c(beta, 1), the device an `offset()` term is in a GLM -- so the
# kernels stay untouched and an arm with no field is handed back exactly what it
# was given (a field-free fit simulates byte-identically). `n_arms` is how many
# leading processes the kernel reads; returns the augmented designs, the
# augmented draw block, and the per-arm column counts to pass alongside them.
.tobs_sim_arm_block <- function(model, draws, n_arms) {
  p <- vapply(model$process_info[seq_len(n_arms)], function(pp) pp$p, integer(1))
  X <- model$X_processes[seq_len(n_arms)]
  blk  <- vector("list", n_arms)
  base <- 0L
  for (k in seq_len(n_arms)) {
    cols     <- base + seq_len(p[k])
    base     <- base + p[k]
    blk[[k]] <- draws[, cols, drop = FALSE]
    off      <- .tobs_eta_offset(model, k)
    if (!is.null(off)) {
      X[[k]]   <- cbind(X[[k]], as.numeric(off))
      blk[[k]] <- cbind(blk[[k]], 1)
      p[k]     <- p[k] + 1L
    }
  }
  list(X = X, draws = do.call(cbind, blk), p = p)
}

# A posterior interval for one arm's response-scale predictor at a NEW design,
# shared by the observation families whose per-cell marginal thins a Poisson /
# NegBin abundance (distance / fp_occu / dyn_abun -- R/distance.R,
# R/fp_occu.R, R/dyn_abun.R). Design-matrix mode only -- a caller passing no
# `X.0` gets fitted()'s in-sample, field-aware value instead (the
# .tobs_predict_nmix() convention: a new row has no field node of its own, so
# interpolating the surface to it is a separate concern, as on count() and on
# the sampled paths). Design validated against the arm's coefficient count
# (the same message .tobs_predict_nmix() raises), draws propagated through
# .tobs_nmix_response_draws() (R/abun.R). Reports mean/sd/quantile columns via
# .tobs_quantile_df() (R/methods.R), so the column names always match
# `quantiles` rather than a hardcoded q2.5/q50/q97.5.
.tobs_count_arm_predict <- function(object, X.0, arm_idx, quantiles) {
  pi_list <- object$model$process_info
  p       <- vapply(pi_list, function(pp) pp$p, integer(1))
  beta_off <- if (arm_idx > 1L) sum(p[seq_len(arm_idx - 1L)]) else 0L
  p_arm <- p[arm_idx]
  if (ncol(X.0) != p_arm) {
    stop(sprintf("X.0 has %d columns but the %s arm has %d coefficients",
                 ncol(X.0), pi_list[[arm_idx]]$name, p_arm), call. = FALSE)
  }
  pred <- .tobs_nmix_response_draws(object$draws, X.0, beta_off, p_arm,
                                    pi_list[[arm_idx]]$link %||% "log")
  .tobs_quantile_df(pred, .tobs_check_quantiles(quantiles, n = 3L))
}

# A SAMPLED (NUTS) field is not among the coefficient columns the fit keeps, so
# the family's own posterior-mean marginal evaluation runs it at offset 0: the
# reported data log-likelihood -- what `logLik()` / `AIC()` / `BIC()` /
# `glance()` read off `log_prob` -- then describes a model without the field
# while WAIC / LOO score one with it. Re-evaluate it through the same pointwise
# kernel the criteria use, so the two cannot disagree. The deterministic
# (grid-integrated) field paths already evaluate their marginal with the field
# in the predictor, so they keep their own grid-weighted value.
.tobs_nuts_field_loglik <- function(fit) {
  if (!identical(fit$method, "nuts") || is.null(fit$field_eta_offset)) return(fit)
  draws <- fit$draws
  if (is.null(draws) || !is.matrix(draws) || !nrow(draws)) return(fit)
  par <- stats::setNames(colMeans(draws), colnames(draws))
  ll  <- .tobs_laplace_marginal_loglik(.tobs_model_with_nb_size(fit), par)$loglik
  if (!is.finite(ll)) return(fit)
  fit$log_lik  <- ll
  fit$log_prob <- rep(ll, max(1L, length(fit$log_prob)))
  fit
}
