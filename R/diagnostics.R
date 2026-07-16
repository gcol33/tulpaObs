# ============================================================================
# Occupancy-specific diagnostics
# Generic diagnostics (moran_i, durbin_watson, variogram, compare_models,
# model_average) are in tulpa — inherited via tulpa_fit class.
# ============================================================================

# Thread count for the parallel occu_cover pointwise-loglik kernel. An explicit
# n.threads wins; NULL falls back to the occu_cover fit's own outer-grid default
# (all but four logical cores), so WAIC / LOO reuse the machine budget the fit
# used. detectCores() can return NA on exotic platforms, so guard it.
.tobs_ploglik_threads <- function(n.threads = NULL) {
  if (!is.null(n.threads)) return(max(1L, as.integer(n.threads)))
  nc <- parallel::detectCores()
  if (is.na(nc)) nc <- 1L
  max(1L, as.integer(nc) - 4L)
}

#' Model criteria for occupancy / cover models
#'
#' `tobs_waic()`, `tobs_dic()`, and `tobs_cpo()` build the family pointwise
#' log-likelihood matrix once and hand it to the engine's single criteria layer
#' [tulpa::tulpa_criteria()], which derives WAIC / DIC / CPO / LPML / PSIS-LOO.
#' DIC additionally evaluates the deviance at the posterior mean of the
#' parameters, supplied by the family-specific `.tobs_loglik_at_mean()`.
#'
#' @param object A `tobs_fit` object.
#' @param n.draws Posterior draws used to build the pointwise log-likelihood
#'   (the cover / occu_cover paths sample this many; the draw-matrix families use
#'   the first `n.draws` stored draws). Default 1000.
#' @param n.threads Threads for the parallel `occu_cover()` pointwise
#'   log-likelihood (the compact / ragged path). The draw loop is embarrassingly
#'   parallel, so this is the WAIC / LOO analogue of the fit's `n.threads.outer`.
#'   `NULL` (default) uses all but four logical cores, matching the occu_cover
#'   fit's own outer-grid default; other families ignore it.
#' @param loo.unit The cross-validation unit for `tobs_waic()` / `tobs_cpo()`.
#'   `"obs"` (default) is the family's pointwise unit -- one column of the
#'   log-likelihood per plot (cover) or site (occu_cover) -- and is byte-identical
#'   to the call without the argument. `"cell"` switches to leave-one-group-out
#'   cross-validation (LOGO-CV): the fit's own per-observation cell map is
#'   supplied to [tulpa::tulpa_criteria()] as `group`, so each spatial cell is one
#'   fold instead of each plot / site. Implemented for `cover()` (the areal field
#'   node, when sites are grouped via `group_var`) and `occu_cover()` (the
#'   `site_cell` map); a non-spatial fit has no cells, so `"cell"` errors there.
#'   Equivalent to passing `group =` the cell map directly, without hand-building
#'   it.
#' @param ... Forwarded to [tulpa::tulpa_criteria()] (e.g. `chunk_size`, or an
#'   explicit `group =` for a custom leave-one-group-out unit).
#' @return A `tulpa_criteria` object. `tobs_waic()` also carries `$elpd`
#'   (an alias for `elpd_waic`) for back-compatibility.
#' @seealso [tulpa::tulpa_criteria()]
#' @export
tobs_waic <- function(object, n.draws = 1000L, loo.unit = c("obs", "cell"),
                      n.threads = NULL, ...) {
  loo.unit <- match.arg(loo.unit)
  ll_mat <- .tobs_pointwise_loglik(object, n.draws = n.draws,
                                   n.threads = .tobs_ploglik_threads(n.threads))
  dots <- .tobs_criteria_group(object, loo.unit, list(...))
  cr <- do.call(tulpa::tulpa_criteria,
                c(list(ll_mat, criteria = "waic"), dots))
  cr$elpd <- cr$elpd_waic
  cr
}

#' @rdname tobs_waic
#' @export
tobs_dic <- function(object, n.draws = 1000L, n.threads = NULL, ...) {
  ll_mat <- .tobs_pointwise_loglik(object, n.draws = n.draws,
                                   n.threads = .tobs_ploglik_threads(n.threads))
  lam <- .tobs_loglik_at_mean(object, n.draws = n.draws)
  tulpa::tulpa_criteria(ll_mat, criteria = "dic", loglik_at_mean = lam, ...)
}

#' @rdname tobs_waic
#' @export
tobs_cpo <- function(object, n.draws = 1000L, loo.unit = c("obs", "cell"),
                     n.threads = NULL, ...) {
  loo.unit <- match.arg(loo.unit)
  ll_mat <- .tobs_pointwise_loglik(object, n.draws = n.draws,
                                   n.threads = .tobs_ploglik_threads(n.threads))
  dots <- .tobs_criteria_group(object, loo.unit, list(...))
  cr <- do.call(tulpa::tulpa_criteria,
                c(list(ll_mat, criteria = c("loo", "cpo", "lpml"),
                       pointwise = TRUE), dots))
  # LOO-PIT (the INLA cpo$pit analogue): tulpa_criteria does not return it, so
  # add the per-observation leave-one-out PIT for the families that expose the
  # per-draw predictive CDF limits. occu_cover builds them from the field-folded
  # detection-summary marginal, so its LOO-PIT is full-model. The PIT is a
  # per-observation diagnostic, unchanged by the cell-level LOO unit.
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    cr$pit <- .tobs_loo_pit_occu_cover(object, n.draws, ll = ll_mat)
  }
  cr
}

# loo.unit = "cell" -> a leave-one-group-out (LOGO-CV) criteria call: aggregate
# each cell's pointwise log-likelihood columns into one fold. The cell map is
# auto-supplied from the fit, so the caller need not hand-build it. The
# loo.unit -> group plumbing and the explicit-`group` conflict check live here so
# tobs_waic() / tobs_cpo() share one path; "obs" leaves `...` untouched (so an
# explicit `group =` still flows through, byte-identical to the ungrouped call).
.tobs_criteria_group <- function(object, loo.unit, dots) {
  if (identical(loo.unit, "obs")) return(dots)
  if ("group" %in% names(dots)) {
    stop("Pass either `loo.unit = \"cell\"` or an explicit `group = `, ",
         "not both.", call. = FALSE)
  }
  grp <- .tobs_loo_cell_map(object)
  if (is.null(grp)) {
    stop("`loo.unit = \"cell\"` needs a per-observation cell map, but this fit ",
         "carries none. Cell-level (leave-one-group-out) LOO is implemented for ",
         "cover() and occu_cover() with a spatial field; a non-spatial fit has ",
         "no cells, so use the default `loo.unit = \"obs\"`.", call. = FALSE)
  }
  dots$group <- grp
  dots
}

# Per-observation -> cell map matching the column order of the family's pointwise
# log-likelihood (.tobs_pointwise_loglik). tobs_cpo() / tobs_waic() pass it as
# tulpa_criteria(group =) for cell-level LOGO-CV. Each .tobs_ploglik_* builder
# fixes its own column order, so the map is family-specific:
#   * cover()      -- columns are the occupancy-arm rows (enc$occ_data$y order);
#                     spi_full is the per-row spatial node (cell) in that order.
#   * occu_cover() -- columns are the sites (1..n_sites); site_cell maps each
#                     site to its field cell (identity when no group_var).
# Returns NULL when the fit carries no cell structure (a non-spatial cover() fit,
# or a family without a cell-level unit), letting the caller raise a clear error.
.tobs_loo_cell_map <- function(object) {
  if (inherits(object, "cover_fit")) {
    sc <- object$spi_full
    return(if (is.null(sc)) NULL else as.integer(sc))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    n_sites <- object$model$n_sites
    return(as.integer(object$model$site_cell %||% seq_len(n_sites)))
  }
  NULL
}

# Pointwise log-likelihood matrix [n_draws x n_obs], marginalized over the
# latent state, for every family. This is the input WAIC, PSIS-LOO, and LOO
# stacking all consume, so it lives in one place (tobs_stack reuses it). The
# per-family marginal likelihoods mirror the C++ kernels in src/*_likelihood.h.
#
# Fidelity note: the cover() / occu_cover() spatial fits score a *full-model*
# pointwise log-likelihood -- the shared spatial / temporal field is folded into
# the occupancy (and, when copied, cover) predictor per cell and the per-visit
# detection / cover covariates are folded in site-major, so WAIC / DIC / LOO /
# CPO match the INLA / spOccupancy full-model criteria. The joint engine
# samples the grid-integrated field jointly with the arm coefficients
# (.tobs_occu_cover_components -> .tobs_joint_draws), giving the exact integrated
# field uncertainty. The v3 nested-Laplace path stores no joint object, so it
# folds the field from the per-cell marginal posterior (field_table z_mean /
# z_sd, .tobs_occu_cover_v3_field); a single site's pointwise term depends only on
# its own cell's field, so the per-observation marginal is exact, but the joint
# field-coefficient covariance is not reconstructed. The draw-matrix families
# routed through .tobs_eta_draws (single, dynamic, integrated, jsdm, ...) evaluate
# the predictor from the process fixed-effect coefficient draws only; for those
# models any structured field is not added and the score is conditional on the
# fixed-effect predictor.
.tobs_pointwise_loglik <- function(object, n.draws = NULL, n.threads = 1L) {
  nd <- n.draws %||% 1000L
  if (inherits(object, "cover_fit"))
    return(.tobs_ploglik_cover(object, nd, n.threads = n.threads))
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_ploglik_occu_cover(object, nd, n.threads = n.threads))
  }
  # Three-level occupancy + cover: the WAIC / LOO unit is the cell (the top-level
  # marginalised observation); scored over the draw matrix (calibrated under the
  # NUTS path, pseudo-draws under Laplace).
  if (identical(object$model$model_type %||% "NULL", "occu_multiscale_cover")) {
    return(.tobs_ploglik_occu_multiscale_cover(object, nd, n.threads = n.threads))
  }
  # Spatial-factor community occu_cover: the per-cell likelihood needs the latent
  # field, so it is scored over the NUTS draws (calibrated WAIC / LOO -- the point
  # of the NUTS path).
  if (identical(object$model$model_type %||% "NULL", "ms_occu_cover_spatial")) {
    return(.tobs_ploglik_ms_occu_cover_spatial(object, nd, n.threads = n.threads))
  }
  # Community N-mixture (ms_abun): the per-(species, site) likelihood needs the
  # per-species deviations, so it is scored over the NUTS draws (the calibrated
  # WAIC / LOO the NUTS path enables; the Laplace community-mean draws omit them).
  if (identical(object$model$model_type %||% "NULL", "ms_nmix")) {
    return(.tobs_ploglik_ms_nmix(object, nd, n.threads = n.threads))
  }
  # Community count / relative-abundance GLMM: per-(species, site) GLMM density
  # over the community-mean draws with per-species BLUP deviations plugged in
  # (R/ms_count.R).
  if (identical(object$model$model_type %||% "NULL", "ms_count")) {
    return(.tobs_ploglik_ms_count(object, nd))
  }
  # Royle-Nichols: the exact per-site Poisson-marginal over the posterior draws
  # (R/royle_nichols.R).
  if (identical(object$model$model_type %||% "NULL", "royle_nichols")) {
    return(.tobs_ploglik_royle_nichols(object, nd, n.threads = n.threads))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_ttd")) {
    return(.tobs_ploglik_occu_ttd(object, nd, n.threads = n.threads))
  }
  # Community occupancy (ms_occu / ms_dyn_occu / ms_int_occu): per-(species,
  # site) marginal scored over the community-mean pseudo-draws with per-species
  # BLUP deviations plugged in (R/community_ploglik.R).
  if ((object$model$model_type %||% "NULL") %in%
      c("ms_occu", "ms_dyn_occu", "ms_int_occu", "ms_occu_cover")) {
    return(.tobs_ploglik_ms_community(object, nd, n.threads = n.threads))
  }

  draws <- object$draws
  if (is.null(draws) || !is.matrix(draws)) {
    stop("Pointwise log-likelihood needs a posterior draw matrix; ",
         "`object$draws` is missing or not a matrix.", call. = FALSE)
  }
  if (!is.null(n.draws) && n.draws < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  .tobs_ploglik_from_draws(object$model, draws, n.threads = n.threads)
}

# Per-family pointwise log-likelihood given an explicit [n_draws x p] draw
# matrix. Split out from the dispatcher so the posterior-mean evaluation
# (.tobs_loglik_at_mean) drives the same per-family kernels with a one-row mean.
.tobs_ploglik_from_draws <- function(model, draws, n.threads = 1L) {
  mt <- model$model_type %||% "NULL"
  switch(
    mt,
    single     = .tobs_ploglik_replicated(model, draws, n.threads),
    dynamic    = .tobs_ploglik_dynamic(model, draws, n.threads),
    integrated = .tobs_ploglik_integrated(model, draws, n.threads),
    count      = .tobs_ploglik_count(model, draws),
    nmix       = .tobs_ploglik_nmix(model, draws, n.threads),
    removal    = .tobs_ploglik_removal(model, draws, n.threads),
    distance   = .tobs_ploglik_distance(model, draws, n.threads),
    fp_occu    = .tobs_ploglik_fp_occu(model, draws, n.threads),
    dyn_abun   = .tobs_ploglik_dyn_abun(model, draws, n.threads),
    occu_multiscale_cover = .occu_ms_cover_ploglik_core(model, draws, n.threads),
    stop("Pointwise log-likelihood is not implemented for model_type = '",
         mt, "'.", call. = FALSE)
  )
}

# Pointwise log-likelihood at the posterior mean of the parameters, the plug-in
# DIC needs (length n_obs). The draw-matrix families evaluate their per-family
# kernel at the column-mean draw; the cover / occu_cover families plug in the
# posterior-mean linear predictors (and mean dispersion) via family helpers.
.tobs_loglik_at_mean <- function(object, n.draws = 1000L) {
  if (inherits(object, "cover_fit")) {
    return(.tobs_cover_loglik_at_mean(object, n.draws))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_occu_cover_loglik_at_mean(object, n.draws))
  }
  if ((object$model$model_type %||% "NULL") %in%
      c("ms_occu", "ms_dyn_occu", "ms_int_occu", "ms_occu_cover")) {
    return(.tobs_community_loglik_at_mean(object))
  }
  if (identical(object$model$model_type %||% "NULL", "royle_nichols")) {
    mean_draw <- matrix(colMeans(object$draws), nrow = 1L,
                        dimnames = list(NULL, colnames(object$draws)))
    obj_mean  <- object; obj_mean$draws <- mean_draw
    return(as.numeric(.tobs_ploglik_royle_nichols(obj_mean, n.draws = 1L)))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_ttd")) {
    mean_draw <- matrix(colMeans(object$draws), nrow = 1L,
                        dimnames = list(NULL, colnames(object$draws)))
    obj_mean  <- object; obj_mean$draws <- mean_draw
    return(as.numeric(.tobs_ploglik_occu_ttd(obj_mean, n.draws = 1L)))
  }
  draws <- object$draws
  if (is.null(draws) || !is.matrix(draws)) {
    stop("DIC needs a posterior draw matrix to evaluate the deviance at the ",
         "posterior mean; `object$draws` is missing or not a matrix.",
         call. = FALSE)
  }
  mean_draw <- matrix(colMeans(draws), nrow = 1L)
  as.numeric(.tobs_ploglik_from_draws(object$model, mean_draw))
}

# Total marginal log-likelihood at a fixed-effect point estimate -- the value
# logLik() / AIC() / BIC() / glance() report for an EM+Laplace fit (single,
# dynamic, integrated, jsdm). Reuses the family pointwise kernel
# (.tobs_ploglik_from_draws) on a one-row draw matrix, so the reported marginal
# is the same likelihood the WAIC / LOO scoring uses -- one source of truth. The
# unit summed over is the model's marginal observation (the site, with its
# visits / seasons pooled into the closed-form or HMM-forward marginal), so
# `nobs` is the number of those units. `par` is the named fixed-effect vector in
# the process-block layout the draw matrix uses (betas concatenated in process
# order; any trailing visit / random-effect columns are ignored by the kernel,
# matching WAIC). Returns NA when the model_type has no pointwise kernel or the
# evaluation errors, so a fit never fails to assemble over a logLik plumbing
# detail (gcol33/tulpaObs#87).
.tobs_laplace_marginal_loglik <- function(model, par) {
  tryCatch({
    draw <- matrix(as.numeric(par), nrow = 1L)
    colnames(draw) <- names(par)
    ll  <- .tobs_ploglik_from_draws(model, draw)
    val <- sum(ll[1L, ])
    list(loglik = if (is.finite(val)) val else NA_real_,
         nobs   = ncol(ll))
  }, error = function(e) list(loglik = NA_real_, nobs = NA_integer_))
}

# --- shared numerics --------------------------------------------------------

# log(inv_logit(eta)) = log(p) and log(1 - inv_logit(eta)) = log(1 - p),
# computed stably via plogis(log.p) -- the R analogues of the C++
# log_inv_logit / log1m_inv_logit.
.tobs_log_p   <- function(eta) stats::plogis(eta,  log.p = TRUE)
.tobs_log_1mp <- function(eta) stats::plogis(-eta, log.p = TRUE)

# Numerically stable elementwise log(exp(a) + exp(b)).
.tobs_logaddexp <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  both_inf <- is.infinite(m) & (a == b)
  out[both_inf] <- m[both_inf]
  out
}

# Cumulative column offset of process k's beta block within the draw matrix
# (processes are laid out in order, fixed-effect betas first).
.tobs_proc_offset <- function(model, k) {
  if (k == 1L) return(0L)
  sum(vapply(model$process_info[seq_len(k - 1L)],
             function(p) p$p, integer(1)))
}

# Linear-predictor draws for process k: [n_draws x nrow(X_k)].
.tobs_eta_draws <- function(model, draws, k) {
  p_k    <- model$process_info[[k]]$p
  off    <- .tobs_proc_offset(model, k)
  beta_k <- draws[, off + seq_len(p_k), drop = FALSE]
  beta_k %*% t(model$X_processes[[k]])
}

# --- per-family marginal likelihoods ---------------------------------------

# Single-season occupancy: per replicate row, marginalized over z. The per-draw
# linear predictors are built by BLAS here; the per-observation marginal (the
# former R loop, now the C++ oracle in test-occu-family-ploglik-cpp.R) runs in
# cpp_occu_single_ploglik, parallel over observations.
.tobs_ploglik_replicated <- function(model, draws, n.threads = 1L) {
  eta_psi <- .tobs_eta_draws(model, draws, 1L)   # [S x n_obs]
  eta_p   <- .tobs_eta_draws(model, draws, 2L)
  y <- model$y                                   # [n_obs x max_visits], <0 = NA
  storage.mode(y) <- "integer"
  cpp_occu_single_ploglik(eta_psi, eta_p, y, max(1L, as.integer(n.threads)))
}

# count(): a plain GLMM on the observed response. The per-site pointwise
# log-density is Poisson / negative-binomial (log link) or Gaussian (identity),
# evaluated at each draw's linear predictor. The dispersion (negbin size /
# Gaussian variance) is the fixed `count_phi` estimated by the outer loop.
.tobs_ploglik_count <- function(model, draws) {
  eta <- .tobs_eta_draws(model, draws, 1L)       # [S x N]
  # Areal fit: the per-site latent field is part of the log-mean (eta = X beta +
  # f_site), so it must enter the pointwise log-likelihood or WAIC / LOO would
  # score the fixed-effect-only model. The field is a fixed per-site offset (its
  # posterior draws are not sampled on this path), broadcast across the draws.
  fld <- model$count_field_offset
  if (!is.null(fld) && length(fld) == ncol(eta)) {
    eta <- eta + matrix(as.numeric(fld), nrow(eta), length(fld), byrow = TRUE)
  }
  y   <- as.numeric(model$y_count)
  Y   <- matrix(y, nrow(eta), length(y), byrow = TRUE)
  response <- model$response %||% "poisson"
  # Cap the mean to a finite range: an extreme Gaussian-Laplace draw can send
  # exp(eta) to Inf, and dpois(y, Inf) / dnbinom(y, mu = Inf) is NaN. Capping
  # keeps the per-draw density finite (a very negative log-density), so that
  # draw is correctly down-weighted in the WAIC log-mean-exp rather than
  # poisoning the whole column.
  link <- model$link %||% "log"
  # Binomial (logit link): the per-site success probability p = plogis(eta) and
  # the density is the binomial pmf at the site's trial count, not a count mean.
  if (identical(response, "binomial")) {
    p  <- stats::plogis(eta)                       # [S x N]
    nt <- as.numeric(model$n_trials %||% rep(1L, length(y)))
    NT <- matrix(nt, nrow(eta), length(nt), byrow = TRUE)
    return(stats::dbinom(Y, size = NT, prob = pmin(pmax(p, 1e-12), 1 - 1e-12),
                         log = TRUE))
  }
  mu  <- if (identical(link, "log")) exp(eta) else eta
  mu  <- pmin(pmax(mu, 1e-300), 1e8)
  phi <- model$count_phi %||% 1
  switch(response,
    poisson  = stats::dpois(Y, pmax(mu, 1e-300), log = TRUE),
    negbin   = stats::dnbinom(Y, size = phi, mu = pmax(mu, 1e-8), log = TRUE),
    gaussian = stats::dnorm(Y, mu, sqrt(max(phi, 1e-8)), log = TRUE),
    stop("Pointwise log-likelihood: unsupported count response '", response,
         "'.", call. = FALSE))
}

# N-mixture: per site, the latent abundance N integrated out in closed form (the
# Royle 2004 marginal). The observation unit is the site (the per-site marginal
# pools that site's visits), so the pointwise log-likelihood is [n_draws x
# n_sites]. NB is detected by the trailing log_r draw column; the per-draw size is
# r = exp(log_r). Reuses the same nmix_site_marginal() kernel the fit used, so the
# WAIC / LOO scoring is on one source of truth. (For NUTS fits the draws are the
# exact posterior; the laplace path's Gaussian draws also score, the N-mixture
# coefficient marginal being well-behaved.)
.tobs_ploglik_nmix <- function(model, draws, n.threads = 1L) {
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  is_nb <- ("log_r" %in% colnames(draws)) || (ncol(draws) > p_lam + p_p)
  K_max <- as.integer(max(model$y_long) + 100L)
  # Per-draw linear predictors by BLAS; the per-site Royle marginal (the former
  # R loop, still via compute_nmix_site) is batched over draws in the kernel.
  eta_lambda <- draws[, seq_len(p_lam), drop = FALSE] %*% t(X_lambda)  # [S x n_sites]
  eta_p      <- draws[, p_lam + seq_len(p_p), drop = FALSE] %*% t(X_p) # [S x n_obs]
  r_vec <- if (is_nb) exp(draws[, p_lam + p_p + 1L]) else rep(Inf, nrow(draws))
  cpp_nmix_ploglik_batch(as.integer(model$y_long), as.integer(model$site_idx),
                         eta_p, eta_lambda, K_max, as.numeric(r_vec),
                         max(1L, as.integer(n.threads)))
}

# Removal sampling: per site, the latent abundance N integrated out in closed
# form over the depleting-binomial removal likelihood (the same marginal the fit
# used). The observation unit is the site (its passes are pooled), so the
# pointwise log-likelihood is [n_draws x n_sites]. NB is detected by the trailing
# log_r draw column.
.tobs_ploglik_removal <- function(model, draws, n.threads = 1L) {
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  is_nb <- ("log_r" %in% colnames(draws)) || (ncol(draws) > p_lam + p_p)
  y <- as.integer(model$y_long); site_idx <- as.integer(model$site_idx)
  n_sites <- model$n_sites
  # K_max from per-site removal totals, matching .tobs_removal_nuts_marginal.
  site_tot <- tapply(y, factor(site_idx, levels = seq_len(n_sites)), sum)
  site_tot[is.na(site_tot)] <- 0L
  K_max <- as.integer(max(as.integer(site_tot)) + 100L)
  eta_lambda <- draws[, seq_len(p_lam), drop = FALSE] %*% t(X_lambda)
  eta_p      <- draws[, p_lam + seq_len(p_p), drop = FALSE] %*% t(X_p)
  r_vec <- if (is_nb) exp(draws[, p_lam + p_p + 1L]) else rep(Inf, nrow(draws))
  cpp_removal_ploglik_batch(y, site_idx, eta_p, eta_lambda, K_max,
                            as.numeric(r_vec), max(1L, as.integer(n.threads)))
}

# Distance sampling: per site, the latent abundance N integrated out in closed
# form over the binned multinomial-over-N detection likelihood (the same marginal
# the fit used). The observation unit is the site (its bins are pooled), so the
# pointwise log-likelihood is [n_draws x n_sites]. The hazard-rate shape is read
# from the model key; NB is detected by the trailing log_r draw column.
.tobs_ploglik_distance <- function(model, draws, n.threads = 1L) {
  p_lam <- model$process_info[[1]]$p; p_sig <- model$process_info[[2]]$p
  hazard <- identical(model$key, "hazard")
  is_nb  <- "log_r" %in% colnames(draws)
  X_lambda <- model$X_processes[[1]]; X_sigma <- model$X_processes[[2]]
  y <- model$y; storage.mode(y) <- "integer"
  K_max <- as.integer(3L * max(rowSums(y)) + 100L)
  off <- p_lam + p_sig
  eta_lambda <- draws[, seq_len(p_lam), drop = FALSE] %*% t(X_lambda)      # [S x n_sites]
  eta_sigma  <- draws[, p_lam + seq_len(p_sig), drop = FALSE] %*% t(X_sigma)
  eta_b <- if (hazard) draws[, off + 1L] else rep(0, nrow(draws))
  r_vec <- if (is_nb) exp(draws[, off + (if (hazard) 2L else 1L)])
           else rep(Inf, nrow(draws))
  cpp_distance_ploglik_batch(
    y, as.numeric(model$cutpoints), .dist_transect_code(model$transect),
    .dist_key_code(model$key), as.integer(model$quad_order), K_max,
    eta_lambda, eta_sigma, as.numeric(eta_b), as.numeric(r_vec),
    max(1L, as.integer(n.threads)))
}

# False-positive occupancy: per site, the latent occupancy z integrated out in
# closed form over the Miller et al. (2011) multistate marginal. The observation
# unit is the site (its visits are pooled), so the pointwise log-likelihood is
# [n_draws x n_sites]. Coefficient layout: (psi, p11, p10, b) site-level arms.
.tobs_ploglik_fp_occu <- function(model, draws, n.threads = 1L) {
  lay <- .tobs_fp_occu_nuts_layout(model$process_info[[1]]$p,
                                   model$process_info[[2]]$p,
                                   model$process_info[[3]]$p,
                                   model$process_info[[4]]$p)
  X_psi <- model$X_processes[[1]]; X_p11 <- model$X_processes[[2]]
  X_p10 <- model$X_processes[[3]]; X_b   <- model$X_processes[[4]]
  eta_psi <- draws[, lay$psi, drop = FALSE] %*% t(X_psi)
  eta_p11 <- draws[, lay$p11, drop = FALSE] %*% t(X_p11)
  eta_p10 <- draws[, lay$p10, drop = FALSE] %*% t(X_p10)
  eta_b   <- draws[, lay$b,   drop = FALSE] %*% t(X_b)
  cpp_fp_occu_ploglik_batch(as.integer(model$y_long), as.integer(model$site_idx),
                            eta_psi, eta_p11, eta_p10, eta_b,
                            max(1L, as.integer(n.threads)))
}

# Open N-mixture (dyn_abun): per site, the latent abundance sequence integrated
# out by the HMM forward recursion (the same marginal the fit used). The
# observation unit is the site (its seasons / visits are pooled), so the pointwise
# log-likelihood is [n_draws x n_sites]. Layout: (lambda, p, omega, gamma) arms.
.tobs_ploglik_dyn_abun <- function(model, draws, n.threads = 1L) {
  lay <- .tobs_dyn_abun_nuts_layout(model$process_info[[1]]$p,
                                    model$process_info[[2]]$p,
                                    model$process_info[[3]]$p,
                                    model$process_info[[4]]$p)
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  X_omega  <- model$X_processes[[3]]; X_gamma <- model$X_processes[[4]]
  eta_lambda <- draws[, lay$lambda, drop = FALSE] %*% t(X_lambda)  # [S x n_sites]
  eta_p      <- draws[, lay$p,      drop = FALSE] %*% t(X_p)
  eta_omega  <- draws[, lay$omega,  drop = FALSE] %*% t(X_omega)
  eta_gamma  <- draws[, lay$gamma,  drop = FALSE] %*% t(X_gamma)
  use_nb <- identical(model$mixture %||% "poisson", "negbin")
  # eta_logr = 0 mirrors the former loop, which called eval_beta without the log
  # r argument (default 0); the per-site HMM marginal is otherwise unchanged.
  cpp_dyn_abun_ploglik_batch(as.integer(model$y_flat), model$n_sites,
                             model$n_seasons, model$max_visits, model$K_max,
                             eta_lambda, eta_p, eta_omega, eta_gamma, use_nb,
                             rep(0, nrow(draws)), max(1L, as.integer(n.threads)))
}

# Community N-mixture (ms_abun): per (species, site), the latent abundance N
# integrated out in closed form, scored over the NUTS draws of the FULL parameter
# vector (community means + per-species deviations + community covariances). The
# pointwise unit is a (species, site): the per-species marginal pools that
# species-site's visits, so the matrix is [n_draws x (n_species * n_sites)] with
# the per-species blocks laid contiguously. NB is read from the layout. Needs a
# NUTS fit (object$nuts): the Laplace community-mean draws omit the per-species
# deviations, so the per-(species, site) likelihood is not identified from them
# (the same constraint as the spatial-factor community occu_cover path).
.tobs_ploglik_ms_nmix <- function(object, n.draws = 1000L, n.threads = 1L) {
  nd <- object$nuts
  if (is.null(nd) || is.null(nd$draws)) {
    stop("WAIC / LOO for the community N-mixture ms_abun() needs a NUTS fit ",
         "(method = \"nuts\"): the Laplace community-mean draws omit the ",
         "per-species deviations, so the per-(species, site) likelihood is not ",
         "identified from them.", call. = FALSE)
  }
  model <- object$model
  lay   <- nd$layout
  is_nb <- isTRUE(lay$is_nb)
  lf    <- .tobs_ms_nmix_longform(model)
  K_max <- nd$K_max %||% as.integer(max(lf$y) + 100L)
  draws <- nd$draws
  M <- nrow(draws)
  if (!is.null(n.draws) && as.integer(n.draws) < M) {
    idx <- unique(round(seq(1, M, length.out = as.integer(n.draws))))
    draws <- draws[idx, , drop = FALSE]; M <- nrow(draws)
  }
  # The non-centered reconstruction (b = C z per species) and the per-(species,
  # site) Royle marginal (still compute_nmix_site) are batched over draws in the
  # kernel; mirrors the former R loop (via .tobs_ms_abun_nuts_b_from_z) exactly.
  clogr <- if (is_nb) as.integer(lay$chol_logr[1L]) - 1L else 0L
  cpp_ms_nmix_ploglik_batch(
    as.integer(lf$y), as.integer(lf$species_idx), as.integer(lf$site_idx),
    lf$X_p, model$X_processes[[1]], draws,
    as.integer(lay$mu[1L]) - 1L, as.integer(lay$b_off),
    as.integer(lay$chol_lam[1L]) - 1L, as.integer(lay$chol_p[1L]) - 1L, clogr,
    as.integer(lay$p_lam), as.integer(lay$p_p), as.integer(lay$n_species),
    as.integer(model$n_sites), is_nb, as.integer(K_max),
    max(1L, as.integer(n.threads)))
}

# Integrated multi-source: per site, shared psi, detection summed over the
# sources that observed it (src/integrated_occ_likelihood.h). The per-(site,
# source) detection counts are draw-invariant, so they are gathered here and the
# per-site z-marginal (the former R loop, now the C++ oracle in the tests) runs
# in cpp_occu_integrated_ploglik, parallel over sites.
.tobs_ploglik_integrated <- function(model, draws, n.threads = 1L) {
  eta_psi <- .tobs_eta_draws(model, draws, 1L)   # [S x n_sites]
  S <- nrow(eta_psi); n_sites <- model$n_sites
  n_src <- model$n_sources
  K1 <- matrix(0L, n_sites, n_src); K0 <- matrix(0L, n_sites, n_src)
  eta_src <- array(0, dim = c(S, n_sites, n_src))
  for (s in seq_len(n_src)) {
    eta_src[, , s] <- .tobs_eta_draws(model, draws, 1L + s)   # [S x n_sites]
    for (i in seq_len(n_sites)) {
      local <- which(model$site_maps[[s]] + 1L == i)
      if (!length(local)) next
      yvec  <- model$y_sources[[s]][local[1L], ]
      valid <- yvec >= 0
      if (!any(valid)) next
      K1[i, s] <- sum(yvec[valid] == 1)
      K0[i, s] <- sum(valid) - K1[i, s]
    }
  }
  cpp_occu_integrated_ploglik(eta_psi, as.numeric(eta_src), K1, K0, n_src,
                              max(1L, as.integer(n.threads)))
}

# Dynamic (multi-season HMM): per-site forward recursion in log space,
# mirroring src/dyn_occ_likelihood.h. Site-level detection only. The recursion
# (the former R loop, now the C++ oracle in the tests) runs in
# cpp_occu_dynamic_ploglik, parallel over sites.
.tobs_ploglik_dynamic <- function(model, draws, n.threads = 1L) {
  eta_psi1 <- .tobs_eta_draws(model, draws, 1L)
  eta_p    <- .tobs_eta_draws(model, draws, 2L)
  eta_gam  <- .tobs_eta_draws(model, draws, 3L)
  eta_eps  <- .tobs_eta_draws(model, draws, 4L)
  y <- model$y; y[is.na(y)] <- -1L                  # 3D [n_sites x mv x Tn]
  cpp_occu_dynamic_ploglik(
    eta_psi1, eta_p, eta_gam, eta_eps,
    as.integer(y), as.integer(model$n_visits), as.integer(model$any_detected),
    model$n_sites, dim(model$y)[2L], model$n_seasons,
    max(1L, as.integer(n.threads)))
}

#' Posterior predictive check
#' @param object A `tobs_fit` object.
#' @param fit.stat `"freeman-tukey"` (default) or `"chi-squared"`.
#' @param n.samples Number of posterior samples (default 500).
#' @return A list with `fit.y`, `fit.y.rep`, and `bayesian.p`.
#' @export
tobs_ppc <- function(object, fit.stat = c("freeman-tukey", "chi-squared"),
                     n.samples = 500) {
  fit.stat <- match.arg(fit.stat)
  if (inherits(object, "cover_fit")) {
    return(.tobs_ppc_cover(object, fit.stat, n.samples))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_ppc_occu_cover(object, fit.stat, n.samples))
  }
  model <- object$model
  if (model$model_type != "single") {
    stop("tobs_ppc supports single-season occupancy, cover(), and occu_cover() ",
         "fits.", call. = FALSE)
  }

  draws <- object$draws; pi_list <- model$process_info
  X_occ <- model$X_processes[[1]]; X_det <- model$X_processes[[2]]
  y <- model$y; n_sites <- model$n_sites; max_visits <- model$max_visits
  p_occ <- pi_list[[1]]$p; p_det <- pi_list[[2]]$p
  n.samples <- min(n.samples, nrow(draws))
  draw_idx <- sample.int(nrow(draws), n.samples)

  # Per-site valid mask, visit count, and whether the species was ever detected.
  # The latent z is sampled from its full conditional given the detection history
  # (the spOccupancy ppcOcc construction), then the detection replicate y_rep ~
  # Bernoulli(z p). The per-draw simulation (the former R loop) runs in
  # cpp_single_ppc, drawing from R's RNG stream in the same order (the draw
  # selection sample.int stays in R), so under a fixed seed it is byte-identical.
  valid_mat <- y >= 0
  n_valid   <- rowSums(valid_mat)
  any_det   <- rowSums(y * valid_mat) > 0
  yint <- y; storage.mode(yint) <- "integer"
  r <- cpp_single_ppc(X_occ, X_det, draws[, seq_len(p_occ + p_det), drop = FALSE],
                      as.integer(draw_idx), yint, as.integer(n_valid),
                      as.integer(any_det), identical(fit.stat, "freeman-tukey"))
  list(fit.y = r$fit.y, fit.y.rep = r$fit.y.rep,
       bayesian.p = mean(r$fit.y.rep > r$fit.y))
}

#' PIT residuals
#' @param object A `tobs_fit` object.
#' @param n.samples Number of posterior samples (default 250).
#' @return Numeric vector of PIT residuals.
#' @export
tobs_pit_residuals <- function(object, n.samples = 250) {
  if (inherits(object, "cover_fit")) {
    return(.tobs_pit_cover(object, n.samples))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_pit_occu_cover(object, n.samples))
  }
  model <- object$model
  if (model$model_type != "single") {
    stop("tobs_pit_residuals supports single-season occupancy, cover(), and ",
         "occu_cover() fits.", call. = FALSE)
  }
  draws <- object$draws; pi_list <- model$process_info
  X_occ <- model$X_processes[[1]]; X_det <- model$X_processes[[2]]
  y <- model$y
  p_occ <- pi_list[[1]]$p; p_det <- pi_list[[2]]$p
  n_draws <- min(n.samples, nrow(draws))
  draw_idx <- sample.int(nrow(draws), n_draws)   # draw selection stays in R
  yint <- y; storage.mode(yint) <- "integer"
  # Per-site posterior-mean predictive CDF plus a uniform jitter (the former R
  # loop) in cpp_single_pit; the jitter draws from R's stream in the same order.
  cpp_single_pit(X_occ, X_det, draws[, seq_len(p_occ + p_det), drop = FALSE],
                 as.integer(draw_idx), yint)
}

#' Goodness-of-fit tests for a tobs fit
#'
#' Posterior-predictive goodness-of-fit checks for an occupancy / N-mixture
#' `tobs_fit`: dispersion, zero-inflation, and outlier counts compared against
#' the fitted model's simulated replicates, plus a Kolmogorov-Smirnov
#' uniformity test on PIT residuals.
#'
#' @param pit Numeric vector of PIT residuals (`tobs_test_uniformity`).
#' @param object A fitted `tobs_fit`.
#' @param n.samples Number of posterior-predictive replicates to simulate.
#' @return A list with the observed statistic, its posterior-predictive
#'   expectation, and a tail p-value; `tobs_test_uniformity` returns the
#'   `htest` object from [stats::ks.test()].
#' @name tobs_gof_tests
NULL

#' @rdname tobs_gof_tests
#' @export
tobs_test_uniformity <- function(pit) {
  # PIT residuals can have ties when the response is discrete (counts /
  # detections). The asymptotic KS p-value remains valid; only the
  # ties-warning is noisy, so suppress just that warning class.
  withCallingHandlers(
    ks.test(pit, "punif"),
    warning = function(w) {
      if (grepl("ties should not be present", conditionMessage(w),
                fixed = TRUE)) invokeRestart("muffleWarning")
    }
  )
}

# The dispersion / zero-inflation / outlier GOF tests use single-season
# detection-history semantics (per-site totals over an N x J 0/1/NA matrix, and
# `== 1` detection counts for the outlier envelope). They mis-compute silently on
# a 3D / long-form / cover response, so require the single-season occupancy
# model_type and point elsewhere for the shapes that have a dedicated path.
.tobs_gof_require_single <- function(object, fn) {
  mt <- object$model$model_type %||% "NULL"
  if (!identical(mt, "single")) {
    hint <- if (inherits(object, "cover_fit") ||
                identical(mt, "occu_cover")) {
      " Use tobs_ppc() for cover() / occu_cover() fits."
    } else {
      paste0(" Multi-season and community goodness-of-fit tests are not ",
             "implemented; use tobs_waic() for those families. (The count ",
             "families abun / removal / distance / dyn_abun have per-site-total ",
             "dispersion / zero-inflation / outlier tests.)")
    }
    stop(sprintf("%s supports single-season occupancy fits only (model_type = %s).%s",
                 fn, mt, hint), call. = FALSE)
  }
  invisible(object)
}

#' @rdname tobs_gof_tests
#' @export
tobs_test_dispersion <- function(object, n.samples = 250) {
  if ((object$model$model_type %||% "NULL") %in% .tobs_count_gof_families) {
    return(.tobs_test_dispersion_count(object, n.samples))
  }
  .tobs_gof_require_single(object, "tobs_test_dispersion")
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  obs_var <- var(rowSums(y_obs * (y_obs >= 0), na.rm = TRUE))
  sim_vars <- vapply(sims, function(ys) var(rowSums(ys * (ys >= 0), na.rm = TRUE)), double(1))
  list(observed = obs_var, expected = mean(sim_vars),
       ratio = obs_var / mean(sim_vars), p.value = mean(sim_vars >= obs_var))
}

#' @rdname tobs_gof_tests
#' @export
tobs_test_zero_inflation <- function(object, n.samples = 250) {
  if ((object$model$model_type %||% "NULL") %in% .tobs_count_gof_families) {
    return(.tobs_test_zero_inflation_count(object, n.samples))
  }
  .tobs_gof_require_single(object, "tobs_test_zero_inflation")
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  count_zeros <- function(y) sum(apply(y, 1, function(r) { v <- r >= 0; all(r[v] == 0) }))
  obs <- count_zeros(y_obs); sim <- vapply(sims, count_zeros, integer(1))
  list(observed = obs, expected = mean(sim), ratio = obs / max(mean(sim), 1),
       p.value = mean(sim >= obs))
}

#' @rdname tobs_gof_tests
#' @export
tobs_test_outliers <- function(object, n.samples = 250) {
  if ((object$model$model_type %||% "NULL") %in% .tobs_count_gof_families) {
    return(.tobs_test_outliers_count(object, n.samples))
  }
  .tobs_gof_require_single(object, "tobs_test_outliers")
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  n_sites <- nrow(y_obs)
  obs_det <- apply(y_obs, 1, function(r) sum(r[r >= 0] == 1))
  sim_det <- vapply(sims, function(ys) apply(ys, 1, function(r) sum(r[r >= 0] == 1)), double(n_sites))
  lower <- apply(sim_det, 1, quantile, 0.025); upper <- apply(sim_det, 1, quantile, 0.975)
  n_out <- sum(obs_det < lower | obs_det > upper)
  sim_out <- vapply(seq_len(n.samples), function(s) sum(sim_det[,s] < lower | sim_det[,s] > upper), integer(1))
  list(n_outliers = n_out, expected = mean(sim_out), p.value = mean(sim_out >= n_out))
}

#' Comprehensive model checking
#' @param object A `tobs_fit` object.
#' @param coords Optional coordinates for spatial diagnostics.
#' @param n.samples Posterior samples for simulation tests.
#' @return Invisibly, diagnostic results.
#' @export
tobs_check <- function(object, coords = NULL, n.samples = 250) {
  cat("=== tobs Model Diagnostics ===\n\n")
  if (identical(object$method, "nuts")) {
    cat(sprintf("Sampler: %d samples, %d divergent, mean accept = %.3f\n",
                object$n_samples, sum(object$divergent), mean(object$accept_prob)))
  } else {
    cat(sprintf("Fit: %s, %d posterior draws (no NUTS sampler diagnostics)\n",
                object$method %||% "laplace", object$n_samples))
  }

  w <- tryCatch(tobs_waic(object), error = function(e) NULL)
  if (!is.null(w)) cat(sprintf("\nWAIC: %.1f (p_waic = %.1f)\n", w$waic, w$p_waic))

  dic <- tryCatch(tobs_dic(object), error = function(e) NULL)
  if (!is.null(dic) && is.finite(dic$dic))
    cat(sprintf("DIC:  %.1f (p_DIC = %.1f)\n", dic$dic, dic$p_dic))

  cpo <- tryCatch(tobs_cpo(object), error = function(e) NULL)
  if (!is.null(cpo)) cat(sprintf("LPML: %.1f (elpd_loo = %.1f)\n",
                                 cpo$lpml, cpo$elpd_loo))

  ppc <- tryCatch(tobs_ppc(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(ppc)) {
    cat(sprintf("\nPPC: Bayesian p = %.3f\n", ppc$bayesian.p))
    if (ppc$bayesian.p < 0.05 || ppc$bayesian.p > 0.95) cat("  WARNING: poor fit\n")
  }

  zi <- tryCatch(tobs_test_zero_inflation(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(zi)) cat(sprintf("\nZero-inflation: obs=%d, exp=%.1f, p=%.3f\n", zi$observed, zi$expected, zi$p.value))

  if (!is.null(coords)) {
    mi <- tryCatch(tulpa::moran_i(residuals(object)$occ, coords), error = function(e) NULL)
    if (!is.null(mi)) {
      cat(sprintf("\nMoran's I: %.3f (p = %.3f)\n", unname(mi$statistic), mi$p.value))
      if (mi$p.value < 0.05) cat("  WARNING: spatial autocorrelation\n")
    }
  }
  cat("\n")
  invisible(list(waic = w, dic = dic, cpo = cpo, ppc = ppc,
                 zero_inflation = zi))
}
