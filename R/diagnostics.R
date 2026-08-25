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
#' `waic()`, `loo()`, `dic()`, and `cpo()` on a `tobs_fit` build the family
#' pointwise log-likelihood matrix once and hand it to the engine's single
#' criteria layer [tulpa::tulpa_criteria()], which derives WAIC / DIC / CPO /
#' LPML / PSIS-LOO. DIC additionally evaluates the deviance at the posterior
#' mean of the parameters, supplied by the family-specific
#' `.tobs_loglik_at_mean()`.
#'
#' `waic()` and `loo()` are the \pkg{loo} package's generics, so
#' [loo::loo_compare()] and the rest of that ecosystem read a `tobs_fit`
#' directly; `loo()` returns a genuine `loo` object built from the same
#' pointwise matrix, via PSIS with relative effective sample sizes. `dic()` and
#' `cpo()` are \pkg{tulpa}'s, and return a `tulpa_criteria` object.
#'
#' @param x,object A `tobs_fit` object.
#' @param n.draws Posterior draws used to build the pointwise log-likelihood
#'   (the cover / occu_cover paths sample this many; the draw-matrix families use
#'   the first `n.draws` stored draws). Default 1000.
#' @param n.threads Threads for the parallel `occu_cover()` pointwise
#'   log-likelihood (the compact / ragged path). The draw loop is embarrassingly
#'   parallel, so this is the WAIC / LOO analogue of the fit's `n.threads.outer`.
#'   `NULL` (default) uses all but four logical cores, matching the occu_cover
#'   fit's own outer-grid default; other families ignore it.
#' @param loo.unit The cross-validation unit for `waic()` / `loo()` / `cpo()`.
#'   `"obs"` (default) is the family's pointwise unit -- one column of the
#'   log-likelihood per plot (cover) or site (occu_cover) -- and is byte-identical
#'   to the call without the argument. `"cell"` switches to leave-one-group-out
#'   cross-validation (LOGO-CV): the fit's own per-observation cell map folds the
#'   log-likelihood columns of a spatial cell into one, so each cell is a fold
#'   instead of each plot / site. `waic()` and `cpo()` hand the map to
#'   [tulpa::tulpa_criteria()] as `group`; `loo()` applies the same fold itself
#'   and runs PSIS on the resulting `[n_draws x n_cells]` matrix, whose column is
#'   the cell's joint conditional log-likelihood per draw, so the importance
#'   ratio inverts the whole cell's likelihood. A whole cell is a larger
#'   perturbation of the posterior than a single row, so the Pareto k values are
#'   correspondingly higher and are the diagnostic to read before trusting the
#'   cell-level number. Implemented for `cover()` (the areal field node, when
#'   sites are grouped via `group_var`) and `occu_cover()` (the `site_cell` map);
#'   a non-spatial fit has no cells, so `"cell"` errors there. Equivalent to
#'   passing `group =` the cell map directly, without hand-building it. `dic()`
#'   has no cross-validation unit -- it is a plug-in deviance over all
#'   observations -- and rejects `loo.unit`.
#' @param ... Forwarded to [tulpa::tulpa_criteria()] (e.g. `chunk_size`, or an
#'   explicit `group =` for a custom leave-one-group-out unit). `loo()` builds
#'   its object through [loo::loo()] rather than the criteria layer, so it reads
#'   `group =` and ignores the rest.
#' @return `dic()` and `cpo()` return a `tulpa_criteria` object; `waic()`
#'   returns one carrying `$elpd` as an alias for `elpd_waic`. `loo()` returns
#'   a `loo` object from the \pkg{loo} package, carrying the cross-validation
#'   unit it scored in its `"loo_unit"` attribute.
#' @seealso [tulpa::tulpa_criteria()], [loo::loo_compare()]
#' @name tobs_criteria
NULL

#' @rdname tobs_criteria
#' @export
waic.tobs_fit <- function(x, n.draws = 1000L, loo.unit = c("obs", "cell"),
                          n.threads = NULL, ...) {
  loo.unit <- match.arg(loo.unit)
  ll_mat <- .tobs_pointwise_loglik(x, n.draws = n.draws,
                                   n.threads = .tobs_ploglik_threads(n.threads))
  dots <- .tobs_criteria_group(x, loo.unit, list(...))
  cr <- do.call(tulpa::tulpa_criteria,
                c(list(ll_mat, criteria = "waic"), dots))
  cr$elpd <- cr$elpd_waic
  cr
}

# The PSIS door. loo::loo() wants the pointwise matrix and the relative
# effective sample sizes off the chain layout, which is what the stacking path
# assembles per member, so both go through `.tobs_loo_one()`. The LOO unit is one
# column of that matrix, so `loo.unit` acts on the matrix rather than on the
# call: the group resolved by `.tobs_criteria_group()` folds the columns before
# PSIS sees them, and the fold count becomes the number of PSIS folds.
#' @rdname tobs_criteria
#' @export
loo.tobs_fit <- function(x, n.draws = 1000L, loo.unit = c("obs", "cell"),
                         n.threads = NULL, ...) {
  loo.unit <- match.arg(loo.unit)
  ll_mat <- .tobs_pointwise_loglik(x, n.draws = n.draws,
                                   n.threads = .tobs_ploglik_threads(n.threads))
  dots <- .tobs_criteria_group(x, loo.unit, list(...))
  out <- .tobs_loo_one(.tobs_loglik_fold_group(ll_mat, dots$group), x$chain_id)
  attr(out, "loo_unit") <- loo.unit
  out
}

#' @rdname tobs_criteria
#' @export
dic.tobs_fit <- function(object, n.draws = 1000L, n.threads = NULL, ...) {
  # DIC is a plug-in deviance over all observations, so it has no fold and
  # tulpa_criteria() takes no `loo.unit`. Reject it at the door rather than
  # letting it reach that call as an unused argument.
  if ("loo.unit" %in% names(list(...))) {
    stop("`dic()` has no cross-validation unit: DIC is a plug-in deviance over ",
         "all observations, unaffected by the LOO fold. `loo.unit` applies to ",
         "waic(), loo() and cpo().", call. = FALSE)
  }
  ll_mat <- .tobs_pointwise_loglik(object, n.draws = n.draws,
                                   n.threads = .tobs_ploglik_threads(n.threads))
  lam <- .tobs_loglik_at_mean(object, n.draws = n.draws)
  tulpa::tulpa_criteria(ll_mat, criteria = "dic", loglik_at_mean = lam, ...)
}

#' @rdname tobs_criteria
#' @export
cpo.tobs_fit <- function(object, n.draws = 1000L, loo.unit = c("obs", "cell"),
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
# waic() / loo() / cpo() share one path; "obs" leaves `...` untouched (so an
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

# The group fold on the PSIS side. tulpa_criteria() reduces a grouped call
# internally, but loo::loo() scores the matrix it is handed, so the same
# reduction is applied here first: sum each group's per-draw pointwise
# log-likelihoods into one column, giving the group's joint conditional
# log-likelihood given each draw. PSIS on that column is leave-one-group-out --
# the importance ratio is 1 / p(y_group | theta_s) instead of
# 1 / p(y_i | theta_s). Column order follows sort(unique(group)), the same order
# factor(group) gives tulpa_criteria(), so the two agree fold for fold. A NULL
# group is the ungrouped matrix, returned untouched.
.tobs_loglik_fold_group <- function(ll, group = NULL) {
  if (is.null(group)) return(ll)
  if (length(group) != ncol(ll)) {
    stop("`group` must have length n_obs = ", ncol(ll), "; got ", length(group),
         ".", call. = FALSE)
  }
  folded <- t(rowsum(t(ll), group = group, reorder = TRUE))
  dimnames(folded) <- NULL
  folded
}

# Per-observation -> cell map matching the column order of the family's pointwise
# log-likelihood (.tobs_pointwise_loglik). cpo() / waic() pass it as
# tulpa_criteria(group =) and loo() folds the matrix by it, both for cell-level
# LOGO-CV. Each .tobs_ploglik_* builder
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
  if (identical(object$model$model_type %||% "NULL", "occu_multi")) {
    return(.tobs_ploglik_occu_multi(object, nd, n.threads = n.threads))
  }
  if (identical(object$model$model_type %||% "NULL", "gdistremoval")) {
    return(.tobs_ploglik_gdistremoval(object, nd, n.threads = n.threads))
  }
  if (identical(object$model$model_type %||% "NULL", "distsamp_open")) {
    return(.tobs_ploglik_distsamp_open(object, nd, n.threads = n.threads))
  }
  if (identical(object$model$model_type %||% "NULL", "double_observer")) {
    return(.tobs_ploglik_double_observer(object, nd, n.threads = n.threads))
  }
  if (identical(object$model$model_type %||% "NULL", "dyn_int_occu")) {
    return(.tobs_ploglik_dyn_int_occu(object, nd, n.threads = n.threads))
  }
  # Community occupancy (ms_occu / ms_dyn_occu / ms_int_occu) and community
  # binned distance sampling (ms_distance): per-(species, site) marginal scored
  # over the community-mean pseudo-draws with per-species BLUP deviations
  # plugged in (R/community_ploglik.R).
  if ((object$model$model_type %||% "NULL") %in%
      c("ms_occu", "ms_dyn_occu", "ms_int_occu", "ms_occu_cover", "ms_distance")) {
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
  .tobs_ploglik_from_draws(.tobs_model_with_nb_size(object), draws,
                           n.threads = n.threads)
}

# An areal count fit integrates the negative-binomial size over the outer
# hyperparameter grid rather than estimating it as a draw coordinate, so the
# per-family kernels cannot read it off the draw matrix and would otherwise score
# the Poisson marginal for a negbin fit. Hand them a model carrying the
# grid-integrated posterior mean, which .tobs_count_nb_size() reads when there is
# no log_r column. `.tobs_fit_model()` re-attaches the outer model to the fit
# after the family builder returns, so this is the boundary where the size and
# the model are both in scope.
.tobs_model_with_nb_size <- function(object) {
  model <- object[["model"]]
  if (is.null(model) || "log_r" %in% colnames(object[["draws"]])) return(model)
  r <- object[["nmix_dispersion"]]$r %||% NA_real_
  if (is.finite(r)) model$nb_r <- as.numeric(r)
  model
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
    occu_multiscale_cover = .occu_mscale_cover_ploglik_core(model, draws, n.threads),
    stop("Pointwise log-likelihood is not implemented for model_type = '",
         mt, "'.", call. = FALSE)
  )
}

# Families whose pointwise log-likelihood kernel is `.tobs_ploglik_<model_type>(
# fit, n.draws)` -- it reads the whole fit rather than a bare draw matrix, so the
# posterior-mean evaluation hands it a fit whose `draws` is the single mean row.
.TOBS_PLOGLIK_FIT_FAMILIES <- c(
  "royle_nichols", "occu_ttd", "occu_multi", "double_observer",
  "gdistremoval", "distsamp_open", "dyn_int_occu")

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
  # Families whose pointwise kernel takes the FIT (not a raw draw matrix):
  # evaluate it on a one-row fit carrying the posterior-mean draw.
  mt <- object$model$model_type %||% "NULL"
  if (mt %in% .TOBS_PLOGLIK_FIT_FAMILIES) {
    obj_mean <- object
    obj_mean$draws <- matrix(colMeans(object$draws), nrow = 1L,
                             dimnames = list(NULL, colnames(object$draws)))
    fn <- get0(paste0(".tobs_ploglik_", mt), envir = asNamespace("tulpaObs"),
               mode = "function", inherits = FALSE)
    return(as.numeric(fn(obj_mean, n.draws = 1L)))
  }
  draws <- object$draws
  if (is.null(draws) || !is.matrix(draws)) {
    stop("DIC needs a posterior draw matrix to evaluate the deviance at the ",
         "posterior mean; `object$draws` is missing or not a matrix.",
         call. = FALSE)
  }
  # Carry the draw names onto the mean row: the count kernels read their trailing
  # coordinates (log_r, logit_omega, zi_logit) by name, so an unnamed row would
  # score the plain Poisson marginal for a negbin / zero-inflated fit.
  mean_draw <- matrix(colMeans(draws), nrow = 1L,
                      dimnames = list(NULL, colnames(draws)))
  as.numeric(.tobs_ploglik_from_draws(.tobs_model_with_nb_size(object),
                                      mean_draw))
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
# detail.
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

# Negative-binomial size per draw for a count family, read BY NAME. The fit
# appends up to three things after the (lambda, p) coefficient block -- the
# structural-zero logit, the random-effect variance / BLUP block, and log_r --
# each conditionally and in that order, so the coordinate one past the
# coefficients is log_r only when the other two are absent. An areal fit
# integrates the size over the outer hyperparameter grid rather than estimating
# it jointly, so it carries no log_r coordinate at all; its grid-integrated
# posterior mean travels on the model as `nb_r`. Poisson -> Inf.
.tobs_count_nb_size <- function(model, draws) {
  n <- nrow(draws)
  if ("log_r" %in% colnames(draws)) return(exp(as.numeric(draws[, "log_r"])))
  r_grid <- model$nb_r %||% NA_real_
  if (is.finite(r_grid)) return(rep(as.numeric(r_grid), n))
  rep(Inf, n)
}

# Structural-zero mixture layer over a per-site pointwise log-likelihood: a site
# with any detection cannot be a structural zero, an all-zero site mixes the
# point mass in. `ll` is [n_draws x n_sites]; the omega draws are one per row, so
# a vector of that length recycles down each column. This is the same two-
# component mixture the ZIP / ZINB fitters build over the same per-site marginal
# (R/nmix_zip.R, .tobs_fit_dyn_abun_zip), applied in R over what the kernel
# returns rather than inside it. A fit whose draws carry no `coord` column is a
# plain Poisson / NB fit and is returned unchanged.
.tobs_count_zi_mix <- function(ll, all_zero, draws, coord) {
  if (!(coord %in% colnames(draws))) return(ll)
  omega <- stats::plogis(as.numeric(draws[, coord]))
  out <- ll + log1p(-omega)
  if (any(all_zero)) {
    a  <- out[, all_zero, drop = FALSE]
    b  <- log(omega)
    mx <- pmax(a, b)
    out[, all_zero] <- mx + log(exp(a - mx) + exp(b - mx))
  }
  out
}

# Per-site all-zero indicator from a long-form (y, site_idx) response: TRUE where
# the site has no detection at any visit, and for a site contributing no row at
# all. Matches the `az` the ZIP fitter builds (R/nmix_zip.R:30-33).
.tobs_count_all_zero <- function(y, site_idx, n_sites) {
  az  <- rep(TRUE, n_sites)
  agg <- tapply(as.integer(y), as.integer(site_idx), function(v) all(v == 0L))
  az[as.integer(names(agg))] <- as.logical(agg)
  az
}

# The same indicator from dyn_abun's flat response, which the kernel is handed
# directly: `y_flat` is aperm(y, c(2, 3, 1)), so each site owns a contiguous
# block of max_visits * n_seasons entries, and a missing visit is stored as -1
# rather than NA. A count is non-negative, so a block with no positive entry is
# the structural-zero candidate -- which covers a site whose visits are all
# missing. Matches the `az` the ZIP fitter builds (R/dyn_abun.R:237).
.tobs_count_all_zero_flat <- function(y_flat, n_per_site) {
  blocks <- matrix(as.integer(y_flat), nrow = as.integer(n_per_site))
  apply(blocks, 2L, function(v) all(v <= 0L))
}

# N-mixture: per site, the latent abundance N integrated out in closed form (the
# Royle 2004 marginal). The observation unit is the site (the per-site marginal
# pools that site's visits), so the pointwise log-likelihood is [n_draws x
# n_sites]. Reuses the same nmix_site_marginal() kernel the fit used, so the
# WAIC / LOO scoring is on one source of truth. (For NUTS fits the draws are the
# exact posterior; the laplace path's Gaussian draws also score, the N-mixture
# coefficient marginal being well-behaved.) A zero-inflated fit adds the
# structural-zero mixture over that marginal, as the ZIP fitter does.
.tobs_ploglik_nmix <- function(model, draws, n.threads = 1L) {
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  K_max <- as.integer(max(model$y_long) + 100L)
  # Per-draw linear predictors by BLAS; the per-site Royle marginal (the former
  # R loop, still via compute_nmix_site) is batched over draws in the kernel.
  eta_lambda <- draws[, seq_len(p_lam), drop = FALSE] %*% t(X_lambda)  # [S x n_sites]
  eta_p      <- draws[, p_lam + seq_len(p_p), drop = FALSE] %*% t(X_p) # [S x n_obs]
  # An areal / temporal / continuous field is part of the arm it loads on, so it
  # enters the pointwise log-likelihood or WAIC / LOO would score the
  # coefficient-only model (R/field_offset.R).
  eta_lambda <- .tobs_add_eta_offset(eta_lambda, model, 1L)
  eta_p      <- .tobs_add_eta_offset(eta_p, model, 2L)
  r_vec <- .tobs_count_nb_size(model, draws)
  ll <- cpp_nmix_ploglik_batch(as.integer(model$y_long), as.integer(model$site_idx),
                               eta_p, eta_lambda, K_max, as.numeric(r_vec),
                               max(1L, as.integer(n.threads)))
  .tobs_count_zi_mix(ll, .tobs_count_all_zero(model$y_long, model$site_idx,
                                              ncol(ll)),
                     draws, "logit_omega")
}

# Removal sampling: per site, the latent abundance N integrated out in closed
# form over the depleting-binomial removal likelihood (the same marginal the fit
# used). The observation unit is the site (its passes are pooled), so the
# pointwise log-likelihood is [n_draws x n_sites]. The NB size is read from the
# draws by name (.tobs_count_nb_size).
.tobs_ploglik_removal <- function(model, draws, n.threads = 1L) {
  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  y <- as.integer(model$y_long); site_idx <- as.integer(model$site_idx)
  n_sites <- model$n_sites
  # K_max from per-site removal totals, matching .tobs_removal_nuts_marginal.
  site_tot <- tapply(y, factor(site_idx, levels = seq_len(n_sites)), sum)
  site_tot[is.na(site_tot)] <- 0L
  K_max <- as.integer(max(as.integer(site_tot)) + 100L)
  eta_lambda <- .tobs_add_eta_offset(
    draws[, seq_len(p_lam), drop = FALSE] %*% t(X_lambda), model, 1L)
  eta_p      <- .tobs_add_eta_offset(
    draws[, p_lam + seq_len(p_p), drop = FALSE] %*% t(X_p), model, 2L)
  r_vec <- .tobs_count_nb_size(model, draws)
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
  eta_lambda <- .tobs_add_eta_offset(
    draws[, seq_len(p_lam), drop = FALSE] %*% t(X_lambda), model, 1L)  # [S x n_sites]
  eta_sigma  <- .tobs_add_eta_offset(
    draws[, p_lam + seq_len(p_sig), drop = FALSE] %*% t(X_sigma), model, 2L)
  eta_b <- if (hazard) draws[, off + 1L] else rep(0, nrow(draws))
  r_vec <- if (is_nb) exp(draws[, off + (if (hazard) 2L else 1L)])
           else rep(Inf, nrow(draws))
  quad_xptr <- cpp_distance_build_quad(as.numeric(model$cutpoints),
                                       .dist_transect_code(model$transect),
                                       as.integer(model$quad_order))
  cpp_distance_ploglik_batch(
    y, quad_xptr, .dist_key_code(model$key), K_max,
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
  eta_psi <- .tobs_add_eta_offset(draws[, lay$psi, drop = FALSE] %*% t(X_psi),
                                  model, 1L)
  eta_p11 <- .tobs_add_eta_offset(draws[, lay$p11, drop = FALSE] %*% t(X_p11),
                                  model, 2L)
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
  eta_lambda <- .tobs_add_eta_offset(
    draws[, lay$lambda, drop = FALSE] %*% t(X_lambda), model, 1L)  # [S x n_sites]
  eta_p      <- .tobs_add_eta_offset(
    draws[, lay$p,      drop = FALSE] %*% t(X_p), model, 2L)
  eta_omega  <- draws[, lay$omega,  drop = FALSE] %*% t(X_omega)
  eta_gamma  <- draws[, lay$gamma,  drop = FALSE] %*% t(X_gamma)
  # The NB size is estimated jointly with the coefficients and stored as the
  # `log_r` draw column (R/dyn_abun.R), so the kernel is handed the per-draw
  # log size rather than the 0 (size 1) it used before log_r was estimated.
  # "zinb" carries the same dispersion as "negbin" under a structural-zero layer.
  use_nb <- (model$mixture %||% "poisson") %in% c("negbin", "zinb")
  eta_logr <- if ("log_r" %in% colnames(draws)) as.numeric(draws[, "log_r"])
              else rep(0, nrow(draws))
  ll <- cpp_dyn_abun_ploglik_batch(as.integer(model$y_flat), model$n_sites,
                                   model$n_seasons, model$max_visits, model$K_max,
                                   eta_lambda, eta_p, eta_omega, eta_gamma, use_nb,
                                   eta_logr, max(1L, as.integer(n.threads)))
  # Zero-inflation (zip / zinb): a site is a structural-zero candidate when it
  # has no detection in any season. The ZI logit is named `zi_logit`; `omega_*`
  # is dyn_abun's SURVIVAL arm. Passed as an argument rather than bound first,
  # so a fit with no `zi_logit` column never builds a mask it cannot use, and
  # read off the same `y_flat` the kernel scored rather than a second copy of
  # the response.
  .tobs_count_zi_mix(ll,
                     .tobs_count_all_zero_flat(model$y_flat,
                                               model$max_visits * model$n_seasons),
                     draws, "zi_logit")
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
  # kernel; mirrors the former R loop (via .ms_ocs_b_from_z) exactly.
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
#'
#' The discrepancy-statistic check (the spOccupancy `ppcOcc` construction):
#' latent occupancy is drawn from its full conditional, a detection replicate
#' from the fitted model, and the two are compared through a discrepancy
#' statistic to give a Bayesian p-value. This is a test, not a plot, so it is
#' its own generic rather than a method on [tulpa::pp_check()], which draws the
#' graphical check.
#'
#' @param object A `tobs_fit` object.
#' @param fit.stat `"freeman-tukey"` (default) or `"chi-squared"`.
#' @param n.samples Number of posterior samples (default 500).
#' @param ... Passed to methods.
#' @return A list with `fit.y`, `fit.y.rep`, and `bayesian.p`.
#' @export
ppc <- function(object, ...) {
  UseMethod("ppc")
}

#' @rdname ppc
#' @export
ppc.tobs_fit <- function(object, fit.stat = c("freeman-tukey", "chi-squared"),
                         n.samples = 500, ...) {
  fit.stat <- match.arg(fit.stat)
  if (inherits(object, "cover_fit")) {
    return(.tobs_ppc_cover(object, fit.stat, n.samples))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_ppc_occu_cover(object, fit.stat, n.samples))
  }
  model <- object$model
  if (model$model_type != "single") {
    stop("ppc() supports single-season occupancy, cover(), and occu_cover() ",
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

#' PIT residuals for a tobs fit
#'
#' The [tulpa::pit_residuals()] method for `tobs_fit`: the family's own
#' predictive CDF, with the randomisation step a discrete response needs.
#'
#' @param object A `tobs_fit` object.
#' @param n.samples Number of posterior samples (default 250).
#' @param nsim Alias for `n.samples`, and the name [tulpa::test_uniformity()]
#'   forwards under. An explicit `nsim` wins, so a uniformity test run at
#'   `nsim = 2000` draws 2000 samples rather than the default 250.
#' @param seed Optional RNG seed for the posterior-draw selection, forwarded by
#'   [tulpa::test_uniformity()]. Set it and repeated calls return the same
#'   residuals; the caller's RNG stream is restored afterwards.
#' @param observed Accepted and ignored: [tulpa::test_uniformity()] forwards it
#'   for models whose response is supplied separately, and a latent-state fit
#'   carries its own.
#' @param ... Unused.
#' @return Numeric vector of PIT residuals.
#' @seealso [tulpa::test_uniformity()] for the uniformity test on the result.
#' @export
pit_residuals.tobs_fit <- function(object, n.samples = 250, nsim = NULL,
                                   seed = NULL, observed = NULL, ...) {
  # tulpa::test_uniformity() calls pit_residuals(object, observed=, nsim=,
  # seed=). None of those is a prefix of `n.samples`, so without these formals
  # all three fell into `...`: `nsim = 2000` ran at 250 draws and `seed` did not
  # reach the unseeded draw selection, so the same call gave a different KS
  # statistic each time.
  if (!is.null(nsim)) n.samples <- as.integer(nsim)
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(as.integer(seed))
  }
  if (inherits(object, "cover_fit")) {
    return(.tobs_pit_cover(object, n.samples))
  }
  if (identical(object$model$model_type %||% "NULL", "occu_cover")) {
    return(.tobs_pit_occu_cover(object, n.samples))
  }
  model <- object$model
  if (model$model_type != "single") {
    stop("pit_residuals() supports single-season occupancy, cover(), and ",
         "occu_cover() fits.", call. = FALSE)
  }
  draws <- object$draws; pi_list <- model$process_info
  X_occ <- model$X_processes[[1]]; X_det <- model$X_processes[[2]]
  y <- model$y
  p_occ <- pi_list[[1]]$p; p_det <- pi_list[[2]]$p
  n_draws <- min(n.samples, nrow(draws))
  draw_idx <- sample.int(nrow(draws), n_draws)   # draw selection stays in R
  yint <- y; storage.mode(yint) <- "integer"
  # Per-site posterior-mean predictive CDF limits for the detected/all-zero
  # event (cpp_single_pit_cdf); the randomized-PIT interpolation + jitter is
  # tulpa::tulpa_pit()'s job, as for cover()/occu_cover().
  lim <- cpp_single_pit_cdf(X_occ, X_det, draws[, seq_len(p_occ + p_det), drop = FALSE],
                            as.integer(draw_idx), yint)
  tulpa::tulpa_pit(lim$cdf_upper, cdf_lower = lim$cdf_lower)
}

#' Goodness-of-fit tests for a tobs fit
#'
#' The \pkg{tulpa} goodness-of-fit generics for `tobs_fit`: dispersion,
#' zero-inflation, and outlier counts compared against the fitted model's
#' simulated replicates. For the Kolmogorov-Smirnov uniformity test, hand
#' [pit_residuals()] to [tulpa::test_uniformity()], which takes a PIT vector
#' directly.
#'
#' @param object A fitted `tobs_fit`.
#' @param n.samples Number of posterior-predictive replicates to simulate.
#' @param ... Unused.
#' @return A list with the observed statistic (`observed`), its
#'   posterior-predictive expectation (`expected`), their ratio (`ratio`)
#'   and a tail p-value (`p.value`). The names are the same on every
#'   family, so a caller can read the statistic without branching on the
#'   model type.
#' @name tobs_gof_tests
NULL

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
      " Use ppc() for cover() / occu_cover() fits."
    } else {
      paste0(" Multi-season and community goodness-of-fit tests are not ",
             "implemented; use waic() for those families. (The count ",
             "families abun / removal / distance / dyn_abun / count have ",
             "per-site-total dispersion / zero-inflation / outlier tests.)")
    }
    stop(sprintf("%s supports single-season occupancy fits only (model_type = %s).%s",
                 fn, mt, hint), call. = FALSE)
  }
  invisible(object)
}

#' @rdname tobs_gof_tests
#' @export
test_dispersion.tobs_fit <- function(object, n.samples = 250, ...) {
  if ((object$model$model_type %||% "NULL") %in% .tobs_count_gof_families) {
    return(.tobs_test_dispersion_count(object, n.samples))
  }
  .tobs_gof_require_single(object, "test_dispersion()")
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  obs_var <- var(rowSums(y_obs * (y_obs >= 0), na.rm = TRUE))
  sim_vars <- vapply(sims, function(ys) var(rowSums(ys * (ys >= 0), na.rm = TRUE)), double(1))
  list(observed = obs_var, expected = mean(sim_vars),
       ratio = obs_var / mean(sim_vars), p.value = mean(sim_vars >= obs_var),
       sim = sim_vars)
}

#' @rdname tobs_gof_tests
#' @export
test_zero_inflation.tobs_fit <- function(object, n.samples = 250, ...) {
  if ((object$model$model_type %||% "NULL") %in% .tobs_count_gof_families) {
    return(.tobs_test_zero_inflation_count(object, n.samples))
  }
  .tobs_gof_require_single(object, "test_zero_inflation()")
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  count_zeros <- function(y) sum(apply(y, 1, function(r) { v <- r >= 0; all(r[v] == 0) }))
  obs <- count_zeros(y_obs); sim <- vapply(sims, count_zeros, integer(1))
  list(observed = obs, expected = mean(sim), ratio = obs / max(mean(sim), 1),
       p.value = mean(sim >= obs))
}

#' @rdname tobs_gof_tests
#' @export
test_outliers.tobs_fit <- function(object, n.samples = 250, ...) {
  if ((object$model$model_type %||% "NULL") %in% .tobs_count_gof_families) {
    return(.tobs_test_outliers_count(object, n.samples))
  }
  .tobs_gof_require_single(object, "test_outliers()")
  sims <- simulate(object, nsim = n.samples); y_obs <- object$model$y
  n_sites <- nrow(y_obs)
  obs_det <- apply(y_obs, 1, function(r) sum(r[r >= 0] == 1))
  sim_det <- vapply(sims, function(ys) apply(ys, 1, function(r) sum(r[r >= 0] == 1)), double(n_sites))
  lower <- apply(sim_det, 1, quantile, 0.025); upper <- apply(sim_det, 1, quantile, 0.975)
  n_out <- sum(obs_det < lower | obs_det > upper)
  sim_out <- vapply(seq_len(n.samples), function(s) sum(sim_det[,s] < lower | sim_det[,s] > upper), integer(1))
  # `observed` / `ratio` are the names the count branch and the two sibling
  # tests already return, so a caller reads the statistic the same way on
  # every family.
  list(observed = n_out, expected = mean(sim_out),
       ratio = n_out / max(mean(sim_out), 1),
       p.value = mean(sim_out >= n_out))
}

#' Comprehensive model checking
#'
#' The [tulpa::check_model()] method for `tobs_fit`: the criteria, the
#' posterior-predictive check and the goodness-of-fit tests in one call,
#' printed as a report and drawn as the diagnostic panel.
#'
#' The engine's default method builds its panel from `fitted()` and
#' `residuals()` as numeric vectors, which a latent-state fit does not have --
#' `fitted()` returns `list(psi, p, z)` and `residuals()` one series per
#' process, so the response cannot be resolved. This method reads the same
#' quantities through the family's own doors instead: [pit_residuals()],
#' [ppc()], [test_dispersion()] and [tulpa::moran_i()].
#'
#' @param object A `tobs_fit` object.
#' @param coords Optional `n_sites x 2` coordinate matrix. Adds Moran's I on
#'   the occupancy residuals and the correlogram panel.
#' @param n.samples Posterior samples for simulation tests.
#' @param plot Draw the panel. `FALSE` prints the report alone, for a log or a
#'   headless run.
#' @param ... Unused.
#' @return Invisibly, the diagnostic results: `waic`, `dic`, `cpo`, `ppc`,
#'   `zero_inflation`, `dispersion`, `pit`, `uniformity`, and `moran` when
#'   `coords` is supplied.
#' @importFrom graphics abline hist par
#' @importFrom grDevices adjustcolor
#' @export
check_model.tobs_fit <- function(object, coords = NULL, n.samples = 250,
                                 plot = TRUE, ...) {
  cat("=== tobs Model Diagnostics ===\n\n")
  if (identical(object$method, "nuts")) {
    cat(sprintf("Sampler: %d samples, %d divergent, mean accept = %.3f\n",
                object$n_samples, sum(object$divergent), mean(object$accept_prob)))
  } else {
    cat(sprintf("Fit: %s, %d posterior draws (no NUTS sampler diagnostics)\n",
                object$method %||% "laplace", object$n_samples))
  }

  w <- tryCatch(waic(object), error = function(e) NULL)
  if (!is.null(w)) cat(sprintf("\nWAIC: %.1f (p_waic = %.1f)\n", w$waic, w$p_waic))

  d <- tryCatch(dic(object), error = function(e) NULL)
  if (!is.null(d) && is.finite(d$dic))
    cat(sprintf("DIC:  %.1f (p_DIC = %.1f)\n", d$dic, d$p_dic))

  cp <- tryCatch(cpo(object), error = function(e) NULL)
  if (!is.null(cp)) cat(sprintf("LPML: %.1f (elpd_loo = %.1f)\n",
                                cp$lpml, cp$elpd_loo))

  pc <- tryCatch(ppc(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(pc)) {
    cat(sprintf("\nPPC: Bayesian p = %.3f\n", pc$bayesian.p))
    if (pc$bayesian.p < 0.05 || pc$bayesian.p > 0.95) cat("  WARNING: poor fit\n")
  }

  zi <- tryCatch(test_zero_inflation(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(zi)) cat(sprintf("\nZero-inflation: obs=%d, exp=%.1f, p=%.3f\n", zi$observed, zi$expected, zi$p.value))

  dp <- tryCatch(test_dispersion(object, n.samples = n.samples), error = function(e) NULL)
  if (!is.null(dp)) cat(sprintf("Dispersion:     ratio=%.2f, p=%.3f\n",
                                dp$ratio, dp$p.value))

  pit <- tryCatch(pit_residuals(object, n.samples = n.samples),
                  error = function(e) NULL)
  ks <- if (is.null(pit)) NULL else
    tryCatch(test_uniformity(pit), error = function(e) NULL)
  if (!is.null(ks)) cat(sprintf("PIT uniformity: KS D=%.3f, p=%.3f\n",
                                unname(ks$statistic), ks$p.value))

  mi <- NULL
  r_occ <- NULL
  if (!is.null(coords)) {
    # A family with no state-level residual (or none registered) has no Moran's
    # I to report; the panel below is keyed on the statistic rather than on the
    # coordinates, so it is dropped rather than drawn from a residual that is
    # not there.
    r_occ <- tryCatch(residuals(object)$occ, error = function(e) NULL)
    mi <- if (is.null(r_occ)) NULL else
      tryCatch(tulpa::moran_i(r_occ, coords), error = function(e) NULL)
    if (!is.null(mi)) {
      cat(sprintf("\nMoran's I: %.3f (p = %.3f)\n", unname(mi$statistic), mi$p.value))
      if (mi$p.value < 0.05) cat("  WARNING: spatial autocorrelation\n")
    }
  }
  cat("\n")
  if (isTRUE(plot)) .tobs_check_panel(object, pit, ks, pc, dp, coords, r_occ, mi)
  invisible(list(waic = w, dic = d, cpo = cp, ppc = pc, zero_inflation = zi,
                 dispersion = dp, pit = pit, uniformity = ks, moran = mi))
}

# The panel behind check_model(). Every quantity arrives already computed by the
# report above, so nothing is simulated twice. A panel whose quantity could not
# be computed for this family is dropped rather than drawn empty, and the layout
# follows how many survive.
.tobs_check_panel <- function(object, pit, ks, pc, dp, coords, r_occ = NULL,
                              mi = NULL) {
  have <- c(pit = !is.null(pit), ppc = !is.null(pc),
            disp = !is.null(dp) && !is.null(dp$sim),
            moran = !is.null(mi))
  n <- sum(have)
  if (n == 0L) return(invisible(NULL))
  old_par <- par(mfrow = if (n >= 4L) c(2, 2) else c(1, n),
                 mar = c(4, 4, 2.5, 1))
  on.exit(par(old_par))
  grey <- adjustcolor("black", 0.6)

  if (have[["pit"]]) {
    m <- length(pit)
    plot(sort(pit), (seq_len(m) - 0.5) / m,
         xlab = "PIT residuals", ylab = "Expected Uniform",
         main = if (is.null(ks)) "QQ Uniform" else
           sprintf("QQ Uniform (KS p = %.3f)", ks$p.value),
         pch = 19, cex = 0.5, col = grey)
    abline(0, 1, col = "red", lty = 2, lwd = 1.5)
  }

  if (have[["ppc"]]) {
    hist(pc$fit.y.rep, breaks = 30, col = "grey85", border = "white",
         main = sprintf("PPC (Bayesian p = %.3f)", pc$bayesian.p),
         xlab = "Replicate fit statistic")
    abline(v = mean(pc$fit.y), col = "red", lwd = 2)
  }

  if (have[["disp"]]) {
    hist(dp$sim, breaks = 20, col = "grey85", border = "white",
         main = sprintf("Dispersion (ratio = %.2f)", dp$ratio),
         xlab = "Simulated variance")
    abline(v = dp$observed, col = "red", lwd = 2)
  }

  if (have[["moran"]]) {
    xy <- as.matrix(coords)
    r <- r_occ
    ks_vals <- c(3, 5, 8, 12, 20)
    ks_vals <- ks_vals[ks_vals < nrow(xy)]
    mi <- lapply(ks_vals, function(k)
      tryCatch(tulpa::moran_i(r, coords = xy, weights = "knn", k = k),
               error = function(e) NULL))
    keep <- !vapply(mi, is.null, logical(1))
    if (any(keep)) {
      I_vals <- vapply(mi[keep], function(m) unname(m$statistic), numeric(1))
      p_vals <- vapply(mi[keep], function(m) m$p.value, numeric(1))
      sig <- p_vals < 0.05
      plot(ks_vals[keep], I_vals, type = "b",
           xlab = "k neighbours", ylab = "Moran's I",
           main = "Spatial Correlogram",
           pch = ifelse(sig, 19, 1), col = ifelse(sig, "red", "black"))
      abline(h = -1 / (nrow(xy) - 1), col = "grey50", lty = 2)
    }
  }
  invisible(NULL)
}
