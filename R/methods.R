# ============================================================================
# tobs_fit-specific S3 methods.
# Generic S3 (coef, confint, vcov, logLik, tidy, glance, ranef, plot) are
# inherited from tulpa::tulpa_fit via class = c("tobs_fit", "tulpa_fit").
# `summary` is overridden below to surface the simplified-Laplace skewness
# coefficients when present.
# ============================================================================

#' Summary for tobs_fit, with skewness column when simplified Laplace is used
#'
#' Extends `summary.tulpa_fit` with an extra `skew` column populated from
#' `object$skew` when the fit was produced with an SLA method
#' (`method = "laplace_sla"` / `"nested_laplace_sla"`).
#' All quantile / mean / sd columns come from `object$draws`, which under
#' simplified Laplace are skew-normal samples — so 2.5%/97.5% quantiles are
#' already SLA-corrected.
#'
#' @param object A `tobs_fit` object.
#' @param ... Forwarded to `summary.tulpa_fit`.
#' @return Data frame as for `summary.tulpa_fit`, with extra `skew` column
#'   when `$skew` is present.
#' @export
summary.tobs_fit <- function(object, ...) {
  s <- NextMethod()
  if (!is.null(object$skew)) {
    sk <- rep(NA_real_, nrow(s))
    nm <- rownames(s)
    matched <- intersect(nm, names(object$skew))
    sk[match(matched, nm)] <- object$skew[matched]
    s$skew <- sk
  }
  # Surface cross-chain convergence diagnostics (NUTS) alongside the
  # posterior summary, matched by parameter name.
  if (!is.null(object$convergence)) {
    cv  <- object$convergence
    idx <- match(rownames(s), cv$parameter)
    s$rhat     <- cv$rhat[idx]
    s$ess_bulk <- cv$ess_bulk[idx]
    s$ess_tail <- cv$ess_tail[idx]
  }
  s
}

# Per-family S3 handlers follow one naming rule: `.tobs_<generic>_<model_type>`.
# This resolves one, returning NULL when a family has no handler and the
# generic's own inline branch (the single-season path) should run instead.
# `alias` redirects a model type onto the handler of another family that already
# implements the identical surface.
.tobs_s3_handler <- function(generic, model_type, alias = character(0)) {
  if (model_type %in% names(alias)) model_type <- alias[[model_type]]
  get0(paste0(".tobs_", generic, "_", model_type),
       envir = asNamespace("tulpaObs"), mode = "function", inherits = FALSE)
}

# `newdata` for a family predictor: the named argument, or the positional `X.0`
# when that was handed a data frame (`predict(fit, some_df)` is the common call).
.tobs_resolve_newdata <- function(newdata, X.0) {
  if (is.null(newdata) && is.data.frame(X.0)) X.0 else newdata
}

# Families whose predictor takes `newdata` and a family-specific response type,
# with the type each reports by default. The three community occupancy families
# share one predictor (`.tobs_predict_ms_community`).
.TOBS_PREDICT_NEWDATA_TYPE <- c(
  royle_nichols   = "abundance",
  occu_ttd        = "state",
  occu_multi      = "state",
  double_observer = "abundance",
  gdistremoval    = "abundance",
  distsamp_open   = "abundance",
  dyn_int_occu    = "state",
  ms_occu         = "occupancy",
  ms_dyn_occu     = "occupancy",
  ms_int_occu     = "occupancy")

# Model types whose predictor has no `terms` argument. `distance` / `fp_occu` /
# `dyn_abun` vary covariates through `X.0`; the `newdata` families above build
# their grid from a data frame. `terms` reaches none of them, so a caller who
# passes it gets told, rather than an in-sample fit that quietly ignored it.
.TOBS_PREDICT_NO_TERMS <- c("distance", "fp_occu", "dyn_abun",
                            names(.TOBS_PREDICT_NEWDATA_TYPE))

#' Number of observations
#' @param object A `tobs_fit` object.
#' @param ... Ignored.
#' @return Integer count of non-NA detection history entries.
#' @export
nobs.tobs_fit <- function(object, ...) {
  model <- object$model
  switch(model$model_type,
    single  = sum(model$y >= 0),
    dynamic = sum(model$y_flat >= 0),
    # Long-form response: one row per (site, pass / visit).
    nmix = , removal = , fp_occu = length(model$y_long),
    count = length(model$y_count),
    # Padded response grid: NA marks an unsampled cell.
    distance = , dyn_abun = , ms_nmix = , ms_distance = , ms_occu_cover = ,
    ms_occu_cover_spatial = , occu_multiscale_cover = sum(!is.na(model$y)),
    # Explicit validity mask alongside the response.
    ms_occu = , ms_dyn_occu = , ms_count = sum(model$valid),
    # One validity mask per detection source.
    ms_int_occu = sum(vapply(model$valid, sum, integer(1))),
    NA_integer_)
}

# Re-export tulpa's ranef() generic so it is reachable as tulpaObs::ranef()
# / after library(tulpaObs) without attaching tulpa (which is Imports-only).
#' @importFrom tulpa ranef
#' @export
tulpa::ranef

# Re-export tulpa's tidy() / glance() generics for the same reason: a tobs_fit
# inherits tulpa_fit, so tulpa's tidy.tulpa_fit / glance.tulpa_fit handle it, but
# the generics must be reachable after library(tulpaObs) for tidy(fit) /
# glance(fit) to resolve (gcol33/tulpaObs#87).
#' @importFrom tulpa tidy
#' @export
tulpa::tidy

#' @importFrom tulpa glance
#' @export
tulpa::glance

#' One-row model summary for a tobs_fit
#'
#' Extends the generic `tulpa_fit` glance with the joint nested-Laplace outer
#' grid placement and the outer Pareto-k diagnostic when present. The
#' joint-coupled families (`occu_cover()`, `cover()`, `occu()` spatial,
#' `occu_multiscale_cover()`) carry both at the fit top level
#' (gcol33/tulpaObs#104, gcol33/tulpaObs#187); every other family glances exactly
#' as before.
#'
#' @param x A fitted `tobs_fit`.
#' @param ... Ignored.
#' @return The base one-row `glance` data frame. A joint-coupled fit adds two
#'   outer-grid placement columns, reported whether or not the grid moved so an
#'   inert auto-recenter is visible in a batch summary:
#'   \describe{
#'     \item{`outer_grid_placement`}{`"fixed"` (the grid the fit was given) or
#'       `"auto_recentered"` (the engine moved it onto the hyperparameter mode).}
#'     \item{`outer_grid_recenter_declined`}{On a `"fixed"` placement, why the
#'       recenter did not apply -- for instance a user-pinned axis, or
#'       `"auto_recenter_disabled"`. `NA` when the grid was recentered.}
#'   }
#'   A fit that also requested the Pareto-k diagnostic (`control$diagnose.k =
#'   TRUE`, off by default per gcol33/tulpaObs#101) adds three more:
#'   \describe{
#'     \item{`pareto_k`}{The outer importance-sampling \eqn{\hat{k}} for the
#'       hyperparameter Gaussian summary; `< 0.7` indicates a reliable summary,
#'       `NA` when the diagnostic did not run or the proposal was degenerate.}
#'     \item{`pareto_k_is_ess`}{The importance-sampling effective sample size on
#'       the PSIS-smoothed weights (numeric); `pareto_k_is_ess / control$k.samples`
#'       is the relative IS efficiency.}
#'     \item{`pareto_k_proposal_source`}{How the importance proposal was built
#'       (gcol33/tulpa#116, #121): `"mode_hessian"` from the Laplace curvature at
#'       the hyperparameter mode -- curvature-backed, so the \eqn{\hat{k}} stays
#'       trustworthy even when a sharp posterior collapses the integration grid to
#'       ~1 cell; `"grid_moment"` from the grid-weighted covariance of the
#'       integration nodes (with `"moment_matched"` its refinement); or
#'       `"grid_mixture"`, the local-bump-per-cell mixture matching what the engine
#'       samples on a spread grid.}
#'   }
#' @export
glance.tobs_fit <- function(x, ...) {
  g <- NextMethod()
  # Prefer the promoted top-level fields; fall back to the nested joint object so
  # a fit saved before the promotion (gcol33/tulpaObs#104) still glances its k-hat.
  g <- .tobs_glance_outer_grid(g, x)
  pk <-.tobs_promote_pareto_k(x) %||% .tobs_promote_pareto_k(.tobs_joint_fit(x))
  if (is.null(pk)) return(g)
  if (!is.null(pk$pareto_k))        g$pareto_k <- pk$pareto_k
  if (!is.null(pk$pareto_k_is_ess)) g$pareto_k_is_ess <- pk$pareto_k_is_ess
  g$pareto_k_proposal_source <- pk$pareto_k_proposal_source %||% NA_character_
  g
}

#' Convergence record for a fitted model
#'
#' The public accessor for whether a fit converged, with one return shape across
#' every `tobs()` family. Each family stores its optimiser / EM / sampler verdict
#' under `fit$convergence`, but historically the cover hurdle (`cover()`) put the
#' flag at `fit$converged` instead, so a consumer that read one location got `NA`
#' for the other family (gcol33/tulpaObs#88). These accessors normalise both
#' layouts: `convergence()` returns the full record (`converged`, `n_iter`, and
#' `sla_status` when the simplified-Laplace marginals were used), and
#' `converged()` returns the single logical.
#'
#' @param object A fitted `tobs_fit` (occupancy / abundance / cover / ...).
#' @param ... Ignored.
#' @return `convergence()`: a list with `converged` (logical), `n_iter`
#'   (integer, `NA` for grid / closed-form fits with no iteration count), and
#'   `sla_status` (character, when present). `converged()`: a single `TRUE` /
#'   `FALSE`.
#' @examples
#' \dontrun{
#' fit <- tobs(y ~ 1, data = d, family = occu(), detection = ~1)
#' converged(fit)        # TRUE / FALSE, same call for every family
#' convergence(fit)$n_iter
#' }
#' @export
convergence <- function(object, ...) UseMethod("convergence")

#' @rdname convergence
#' @export
convergence.tobs_fit <- function(object, ...) {
  rec <- object$convergence
  if (is.null(rec) || !is.list(rec)) rec <- list()
  # Normalise: prefer the unified record, fall back to the legacy top-level
  # `converged` / `n_iter` so old saved fits and any family still on the flat
  # layout answer through the same accessor.
  rec$converged <- rec$converged %||% object$converged %||% NA
  rec$n_iter    <- rec$n_iter    %||% object$n_iter    %||% NA_integer_
  if (is.null(rec$sla_status) && !is.null(object$sla_status)) {
    rec$sla_status <- object$sla_status
  }
  rec
}

#' @rdname convergence
#' @export
converged <- function(object, ...) UseMethod("converged")

#' @rdname convergence
#' @export
converged.tobs_fit <- function(object, ...) isTRUE(convergence(object)$converged)

#' Random-effect estimates (BLUPs) for a tobs_fit
#'
#' Returns the per-group random-effect posterior summaries. Under NUTS the
#' non-centred draws are reconstructed to the natural BLUP scale
#' (`b_{g,c} = sigma_c * (L z_g)_c`) and summarised; the deterministic Laplace
#' path returns the variance-component EM modes and their Schur-complement
#' standard errors. When the fit carries no formula random effects this falls
#' back to the generic flat random-effect table.
#'
#' @param object A `tobs_fit` object.
#' @param ... Ignored.
#' @return A data frame with one row per group level and coefficient
#'   (`group`, `level`, `term`, `estimate`, `std.error`), or the generic
#'   `ranef` table when no `re_effects` are present.
#' @export
ranef.tobs_fit <- function(object, ...) {
  fn <- .tobs_s3_handler("ranef", object$model$model_type)
  if (!is.null(fn)) return(fn(object))
  if (identical(object$model$model_type, "occu_cover") && !is.null(object$re)) {
    # occu_cover() shared-field + per-group RE (gcol33/tulpaObs#56, #102, #103):
    # `fit$re` is a flat list of random-intercept terms, one per arm for a lone
    # term, several for crossed / nested groupings sharing an arm. Stack them into
    # one table with `arm` + `var` (grouping variable) columns; every arm carries
    # its grouping `level` labels.
    rows <- lapply(object$re, function(re) {
      bl <- re$blup; bsd <- re$blup_sd
      if (is.matrix(bl)) {
        # Random slope: one (group, coefficient) row per cell of the
        # [n_groups x n_coefs] BLUP matrix, tagged by the coefficient `term`.
        cn  <- colnames(bl) %||% re$coef_names %||%
               paste0("coef", seq_len(ncol(bl)))
        lev <- re$levels %||% as.character(seq_len(nrow(bl)))
        data.frame(arm     = re$arm,
                   var     = re$var %||% NA_character_,
                   group   = rep(lev, times = ncol(bl)),
                   term    = rep(cn, each = nrow(bl)),
                   level   = rep(seq_len(nrow(bl)), times = ncol(bl)),
                   blup    = as.numeric(bl),
                   blup_sd = as.numeric(bsd),
                   stringsAsFactors = FALSE)
      } else {
        lev <- re$levels %||% as.character(seq_along(bl))
        data.frame(arm     = re$arm,
                   var     = re$var %||% NA_character_,
                   group   = lev,
                   term    = "(Intercept)",
                   level   = seq_along(bl),
                   blup    = bl,
                   blup_sd = bsd,
                   stringsAsFactors = FALSE)
      }
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }
  if (!is.null(object$re) && !is.null(object$re$blup)) {
    # Single-block per-group RE on a non-occu_cover fit (e.g. abun NUTS): one
    # random intercept per group on the named arm, with per-group BLUP and its SD.
    re <- object$re
    return(data.frame(
      arm   = re$arm,
      group = seq_along(re$blup),
      blup  = re$blup,
      blup_sd = re$blup_sd,
      stringsAsFactors = FALSE))
  }
  if (!is.null(object$re_effects) && length(object$re_effects) > 0L) {
    out <- do.call(rbind, object$re_effects)
    rownames(out) <- NULL
    return(out)
  }
  NextMethod()
}

#' Coefficients for a tobs_fit
#'
#' Returns a per-process coefficient list keyed by the model's process names
#' (`psi`, `p`, `gamma`, `lambda`, ...), splitting the generic flat fixed-effect
#' vector on the `<process>_<coef>` name prefix, and appends the visit-level
#' detection coefficients (`p_visit_<cov>`) carried separately from the
#' site-level detection design. Coordinates with no process
#' prefix (e.g. the `log_r` overdispersion nuisance) are not arm coefficients
#' and are omitted from the list (they remain in `vcov()` / `confint()`).
#'
#' @param object A `tobs_fit` object.
#' @param ... Ignored.
#' @return A per-process coefficient list, one named numeric vector per linear
#'   predictor, with visit-level detection coefficients appended to the
#'   detection process when present.
#' @export
coef.tobs_fit <- function(object, ...) {
  cf <- NextMethod()
  if (!is.list(cf)) cf <- .tobs_coef_by_process(cf, object$process_info)

  vn <- object$model$det_visit_names
  if (!is.null(vn) && length(vn) > 0L && is.list(cf)) {
    det_name <- object$process_info[[2]]$name
    pv <- object$means[paste0("p_visit_", vn)]
    names(pv) <- paste0("visit_", vn)
    if (!is.null(cf[[det_name]])) {
      cf[[det_name]] <- c(cf[[det_name]], pv)
    } else {
      cf[["p_visit"]] <- pv
    }
  }
  cf
}

# Split a flat fixed-effect coefficient vector into a per-process list keyed by
# process name. Each name is "<process>_<coef>"; group by prefix in
# process_info order and strip the prefix from the inner names. Returns the
# input unchanged when there is no process_info or no name matches a prefix.
.tobs_coef_by_process <- function(flat, pi_list) {
  if (is.null(pi_list) || is.null(names(flat))) return(flat)
  nm <- names(flat)
  cf <- list()
  for (pp in pi_list) {
    prefix <- paste0(pp$name, "_")
    hit <- startsWith(nm, prefix)
    if (!any(hit)) next
    vals <- flat[hit]
    names(vals) <- sub(paste0("^", prefix), "", nm[hit])
    cf[[pp$name]] <- vals
  }
  if (length(cf) == 0L) flat else cf
}

#' Fitted values (occupancy and detection probabilities)
#' @param object A `tobs_fit` object.
#' @param ... Ignored.
#' @return A list with `psi` (occupancy probabilities), `p` (detection probabilities),
#'   and `z` (posterior state probability at posterior-mean parameters). For
#'   single-season and community models `z` is the per-site Bayes posterior
#'   `P(z=1 | y)`; for dynamic models `z` is the forward-backward (HMM smoothing)
#'   state posterior `P(z_t=1 | y_{1:T})` as an `[n_sites x n_seasons]` matrix.
#' @export
fitted.tobs_fit <- function(object, ...) {
  model <- object$model
  # removal shares the N-mixture fitted surface: same count marginal, one
  # removal unit per site.
  fn <- .tobs_s3_handler("fitted", model$model_type, c(removal = "nmix"))
  if (!is.null(fn)) return(fn(object))
  means <- object$means
  pi_list <- model$process_info

  # Extract occupancy and detection linear predictors at posterior mean
  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  beta_occ <- means[seq_len(pi_list[[1]]$p)]
  beta_det <- means[pi_list[[1]]$p + seq_len(pi_list[[2]]$p)]

  # A fitted latent surface on the state arm (the continuous NNGP svc() field of
  # a Laplace fit, gcol33/tulpaObs#143) is a per-site offset on the occupancy
  # logit, so the in-sample psi / z read it. NULL on every other route, which
  # leaves eta_occ exactly as it was.
  eta_occ <- as.vector(X_occ %*% beta_occ) + (model$occ_eta_offset %||% 0)
  eta_det <- as.vector(X_det %*% beta_det)

  psi <- plogis(eta_occ)
  p <- plogis(eta_det)

  # Compute z posterior: P(z=1 | y)
  if (identical(model$model_type, "single")) {
    y <- model$y
    n_obs <- nrow(y)
    max_visits <- ncol(y)
    z <- numeric(n_obs)
    for (i in seq_len(n_obs)) {
      yi <- y[i, ]
      valid <- yi >= 0
      if (any(yi[valid] == 1)) {
        z[i] <- 1  # Detected -> occupied
      } else {
        # P(z=1|all zeros) = psi * prod(1-p) / [psi*prod(1-p) + (1-psi)]
        n_valid <- sum(valid)
        prod_1mp <- (1 - p[i])^n_valid
        z[i] <- psi[i] * prod_1mp / (psi[i] * prod_1mp + (1 - psi[i]))
      }
    }
  } else if (identical(model$model_type, "dynamic")) {
    z <- .tobs_dynamic_smoothed_z(model, means, pi_list)
  } else {
    z <- psi
  }

  list(psi = psi, p = p, z = z)
}

# Forward-backward (HMM smoothing) state posterior P(z_{i,t} = 1 | y_{i,1:T})
# for a dynamic (MacKenzie et al. 2003) occupancy fit. The forward filter is the
# same recursion the likelihood evaluates (src/dyn_occ_likelihood.h); the
# backward pass folds in the future detection history so each season's state is
# conditioned on the whole series, not just the marginal occupancy probability.
# Returns an [n_sites x n_seasons] matrix. Detection is indexed per (site, season)
# and colonization / extinction per (site, interval), so a season-varying
# detection covariate ([n_sites x T]) or an interval-varying transition covariate
# ([n_sites x (T-1)]) is honoured (gcol33/tulpaObs#124); constant arms broadcast a
# single value across the periods, matching the engine.
.tobs_dynamic_smoothed_z <- function(model, means, pi_list) {
  n_sites   <- model$n_sites
  T_seasons <- model$n_seasons
  n_int     <- T_seasons - 1L

  X <- model$X_processes
  off <- cumsum(c(0L, vapply(pi_list, function(pp) pp$p, integer(1))))
  beta_psi1 <- means[off[1] + seq_len(pi_list[[1]]$p)]
  beta_p    <- means[off[2] + seq_len(pi_list[[2]]$p)]
  beta_gam  <- means[off[3] + seq_len(pi_list[[3]]$p)]
  beta_eps  <- means[off[4] + seq_len(pi_list[[4]]$p)]

  psi1  <- plogis(as.vector(X[[1]] %*% beta_psi1))
  # Per-(site, period) rate matrices: a long-form design (season / interval
  # varying) flattens site-major period-minor and is reshaped byrow; a per-site
  # design broadcasts the single value across the periods.
  ep <- as.vector(X[[2]] %*% beta_p)
  p_mat <- if (isTRUE(model$det_season_varying))
    matrix(plogis(ep), n_sites, T_seasons, byrow = TRUE)
  else matrix(plogis(ep), n_sites, T_seasons)
  eg <- as.vector(X[[3]] %*% beta_gam)
  gam_mat <- if (isTRUE(model$col_season_varying))
    matrix(plogis(eg), n_sites, n_int, byrow = TRUE)
  else matrix(plogis(eg), n_sites, n_int)
  ee <- as.vector(X[[4]] %*% beta_eps)
  eps_mat <- if (isTRUE(model$ext_season_varying))
    matrix(plogis(ee), n_sites, n_int, byrow = TRUE)
  else matrix(plogis(ee), n_sites, n_int)

  y <- model$y  # [n_sites x max_visits x n_seasons]
  z <- matrix(NA_real_, n_sites, T_seasons)

  for (i in seq_len(n_sites)) {
    # Per-season emission likelihood under each state, and a hard-detection mask.
    em <- matrix(1, T_seasons, 2L)  # columns: state 0 (unocc), state 1 (occ)
    for (t in seq_len(T_seasons)) {
      p_it <- p_mat[i, t]
      raw <- y[i, , t]
      raw <- raw[!is.na(raw) & raw >= 0]
      if (length(raw) == 0L) {
        em[t, ] <- c(1, 1)  # no visits: uninformative
      } else if (any(raw == 1)) {
        # A detection rules out the unoccupied state.
        em[t, 1L] <- 0
        em[t, 2L] <- prod(p_it^raw * (1 - p_it)^(1 - raw))
      } else {
        em[t, 1L] <- 1                       # unoccupied -> all non-detections
        em[t, 2L] <- prod(1 - p_it)^length(raw)
      }
    }

    # Per-interval transition matrices Tr_t[a, b] = P(z_{t+1}=b-1 | z_t=a-1).
    Tr_of <- function(iv) matrix(c(1 - gam_mat[i, iv], eps_mat[i, iv],
                                   gam_mat[i, iv],     1 - eps_mat[i, iv]), 2L, 2L)

    # Forward filtering (scaled). Step t-1 -> t uses interval (t - 1)'s rates.
    fwd <- matrix(0, T_seasons, 2L)
    prior <- c(1 - psi1[i], psi1[i])
    a <- prior * em[1L, ]; a <- a / sum(a); fwd[1L, ] <- a
    if (T_seasons > 1L) {
      for (t in 2L:T_seasons) {
        pred <- as.vector(t(Tr_of(t - 1L)) %*% fwd[t - 1L, ])
        a <- pred * em[t, ]; a <- a / sum(a); fwd[t, ] <- a
      }
    }

    # Backward smoothing (Rauch-Tung-Striebel style for discrete HMM). The
    # transition at backward step t is interval t (seasons t -> t+1).
    sm <- matrix(0, T_seasons, 2L)
    sm[T_seasons, ] <- fwd[T_seasons, ]
    if (T_seasons > 1L) {
      for (t in (T_seasons - 1L):1L) {
        Tr <- Tr_of(t)
        pred <- as.vector(t(Tr) %*% fwd[t, ])  # P(z_{t+1} | y_{1:t})
        ratio <- ifelse(pred > 0, sm[t + 1L, ] / pred, 0)
        sm[t, ] <- fwd[t, ] * as.vector(Tr %*% ratio)
        s <- sum(sm[t, ]); if (s > 0) sm[t, ] <- sm[t, ] / s
      }
    }
    z[i, ] <- sm[, 2L]
  }
  z
}

# Logit-arm probability draws at a design matrix X.0: plogis(X.0 %*% beta) for
# every posterior draw, returned as [n_draws x nrow(X.0)]. `offset` is where the
# arm's `p` coefficients start in the packed draw row -- 0 for the occupancy
# arm, p_occ for detection. Used by predict.tobs_fit (design-matrix mode),
# predict.tobs_stack (stacked predictive) and the spOccupancy-compatible
# `$psi.samples` / `$p.samples` accessors; kept in one place so they share the
# same parameterization.
.tobs_psi_draws <- function(draws, X.0, p_occ, offset = 0L) {
  beta <- draws[, offset + seq_len(p_occ), drop = FALSE]
  plogis(beta %*% t(X.0))
}

#' Residuals from occupancy model
#' @param object A `tobs_fit` object.
#' @param type One of `"deviance"` (default), `"pearson"`, or `"response"`.
#' @param ... Ignored.
#' @return A list with `occ` (site-level) and `det` (visit-level) residuals.
#' @export
residuals.tobs_fit <- function(object, type = c("deviance", "pearson", "response"), ...) {
  type <- match.arg(type)
  # The three community occupancy families share one residual surface.
  # (`ms_count` covers jsdm() too -- one model class, bernoulli response.)
  fn <- .tobs_s3_handler(
    "residuals", object$model$model_type,
    c(ms_occu = "ms_community", ms_dyn_occu = "ms_community",
      ms_int_occu = "ms_community"))
  if (!is.null(fn)) return(fn(object, type))
  fit_vals <- fitted(object)
  model <- object$model

  # Occupancy residuals (site-level; site-by-season for dynamic, matching the
  # smoothed fitted()$z matrix). z_obs is the realized "ever-detected" indicator
  # the smoothed state posterior is compared against (NA where the unit had no
  # visits, so the state is unobserved).
  z_obs <- if (identical(model$model_type, "single")) {
    apply(model$y, 1, function(row) as.integer(any(row[row >= 0] == 1)))
  } else if (identical(model$model_type, "dynamic")) {
    y <- model$y
    n_sites <- dim(y)[1]; T_seasons <- dim(y)[3]
    zo <- matrix(NA_real_, n_sites, T_seasons)
    for (i in seq_len(n_sites)) for (t in seq_len(T_seasons)) {
      raw <- y[i, , t]; raw <- raw[!is.na(raw) & raw >= 0]
      if (length(raw)) zo[i, t] <- as.numeric(any(raw == 1))
    }
    zo
  } else {
    rep(NA_real_, model$n_sites)
  }

  occ_resid <- switch(type,
    response = z_obs - fit_vals$z,
    pearson = (z_obs - fit_vals$z) / sqrt(fit_vals$z * (1 - fit_vals$z) + 1e-10),
    deviance = {
      sign(z_obs - fit_vals$z) * sqrt(2 * abs(
        ifelse(z_obs == 1,
               -log(fit_vals$z + 1e-10),
               -log(1 - fit_vals$z + 1e-10))
      ))
    }
  )

  # Detection residuals (visit-level, single-season)
  det_resid <- NULL
  if (identical(model$model_type, "single")) {
    y <- model$y
    p_hat <- fit_vals$p
    n_obs <- nrow(y)
    max_visits <- ncol(y)
    det_resid <- matrix(NA_real_, n_obs, max_visits)
    for (i in seq_len(n_obs)) {
      for (j in seq_len(max_visits)) {
        if (y[i, j] >= 0) {
          expected <- fit_vals$z[i] * p_hat[i]
          det_resid[i, j] <- switch(type,
            response = y[i, j] - expected,
            pearson = (y[i, j] - expected) / sqrt(expected * (1 - expected) + 1e-10),
            deviance = {
              sign(y[i, j] - expected) * sqrt(2 * abs(
                ifelse(y[i, j] == 1, -log(expected + 1e-10), -log(1 - expected + 1e-10))
              ))
            }
          )
        }
      }
    }
  }

  list(occ = occ_resid, det = det_resid)
}

#' Simulate replicate datasets from posterior
#' @param object A `tobs_fit` object.
#' @param nsim Number of simulated datasets (default 1).
#' @param seed Optional random seed.
#' @param ... Ignored.
#' @return A list of simulated detection history matrices.
#' @export
simulate.tobs_fit <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  model <- object$model
  fn <- .tobs_s3_handler("simulate", model$model_type)
  if (!is.null(fn)) return(fn(object, nsim))
  draws <- object$draws
  n_samples <- nrow(draws)
  pi_list <- model$process_info

  if (model$model_type != "single") {
    stop("simulate() currently only supports single-season models")
  }

  # The posterior-draw selection (R_unif_index, the sample.int primitive) and the
  # z / y_rep Bernoulli draws run in cpp_simulate_single from R's RNG stream in
  # the same order, so the simulation is byte-identical under a fixed seed.
  p_occ <- pi_list[[1]]$p; p_det <- pi_list[[2]]$p
  yint <- model$y; storage.mode(yint) <- "integer"
  res <- cpp_simulate_single(model$X_processes[[1]], model$X_processes[[2]],
                             draws[, seq_len(p_occ + p_det), drop = FALSE],
                             yint, p_occ, p_det, as.integer(nsim))
  if (nsim == 1) res[[1]] else res
}

#' Predict from occupancy model
#'
#' Four modes:
#' - **In-sample**: `predict(fit)` returns fitted values.
#' - **State posterior / NA-response**: `predict(fit, type = "state")` returns
#'   the marginalised per-site occupancy posterior of a `"nested_laplace"` fit:
#'   the per-row eta posterior is a Gaussian mixture over the hyperparameter
#'   grid (per-cell fitted linear predictor and predictive variance), so `psi`
#'   is its Gauss-Hermite mean and `psi_lower` / `psi_upper` are the mixture-CDF
#'   quantiles -- a calibrated 95% credible interval, exact for every latent
#'   prior. Rows with `heldout = TRUE` are the INLA-style NA-response
#'   prediction targets (single-season sites whose detection history was
#'   all-missing), where occupancy is interpolated from the latent field rather
#'   than informed by detections.
#' - **Design-matrix**: `predict(fit, X.0 = ...)` predicts at new covariate values.
#' - **Terms-based**: `predict(fit, terms = "elev")` varies one covariate, others at mean.
#' - **Joint occu_cover**: for a `occu_cover` joint-coupled fit,
#'   `predict(fit, newdata, type = "occurrence" | "cover_cond" | "cover_exp" |
#'   "change")` samples the joint latent from the grid-integrated posterior
#'   (the outer-grid mixture via [tulpa::tulpa_posterior_draws()]) and
#'   marginalises every derived quantity per draw. `type = "change"` with
#'   `times = c(t1, t2)` returns a per-cell change table (`delta_p`,
#'   `delta_cover_cond`, `delta_cover_exp`, the occupancy / abundance
#'   decomposition, and `.lwr` / `.upr` at `level`), plus the start / end
#'   occupancy `p_T1` / `p_T2` with their own `.sd` / `.lwr` / `.upr`, and a
#'   `.prob_pos` column per headline delta giving the directional posterior
#'   probability `P(delta > 0)` per cell. The result is a
#'   `tobs_prediction` table (one row per cell) carrying per-unit `[cell x nsim]`
#'   draw matrices in `attr(, "draws")`; map it yourself, e.g.
#'   `left_join(cents, pr, by = "cell")` then
#'   `geom_tile(aes(x, y, fill = delta_p))` (or `geom_sf()` on polygon cells).
#' - **Spatial-factor community**: for a reduced-rank spatial-factor
#'   `ms_occu_cover()` fit, `predict(fit, type = "occupancy" | "cover_cond" |
#'   "cover_exp")` returns the per-species per-cell posterior as a long table
#'   (one row per cell x species): `"occupancy"` gives `psi`, `"cover_cond"` the
#'   conditional cover mean `E[cover | present]`, `"cover_exp"` the unconditional
#'   expected cover `psi * E[cover | present]`, each with a `_median` and a
#'   `_lower` / `_upper` interval. The maps are marginalised over the loading +
#'   field posterior, so a rare species borrows strength across the shared
#'   factors for a calibrated map. The latent fields are tied to the cell graph,
#'   so `X.0` / `newdata` prediction is not supported.
#'
#' @param object A `tobs_fit` object.
#' @param X.0 Optional design matrix for occupancy prediction.
#' @param type `"occupancy"` (default), `"detection"`, `"both"`, or `"state"`
#'   (nested-Laplace marginalised per-site psi, incl. held-out sites). For an
#'   `occu_cover` fit: `"occurrence"`, `"cover_cond"`, `"cover_exp"`, or
#'   `"change"`.
#' @param quantiles Quantile levels for credible intervals.
#' @param terms Name of a single term to vary: one column of the predicted
#'   process's design matrix, every other column held at its column mean. A
#'   vector of more than one term is an error, since a second term would be
#'   held at its mean rather than grouped over its levels. For a grid over
#'   more than one covariate, build the design matrix and pass it as `X.0`.
#' @param n_points Number of prediction points per continuous term.
#' @param newdata `occu_cover` only: data.frame of prediction units, one row per
#'   field cell (or carrying a `cell` column mapping rows to field cells).
#'   Defaults to the training data.
#' @param times `occu_cover` `type = "change"` only: length-2 numeric
#'   `c(t1, t2)`, the two values of the time covariate to difference.
#' @param level `occu_cover` only: credible level for the interval columns
#'   (default 0.95).
#' @param nsim `occu_cover` only: number of joint posterior draws (default 1000).
#' @param draws `occu_cover` only: if `TRUE` (default), carry the per-unit
#'   `[cell x nsim]` draw matrices in `attr(, "draws")`.
#' @param time_col `occu_cover` only: name of the time covariate weighting the
#'   trend field / driving the change map; auto-resolved from the fit's stored
#'   trend weight when omitted.
#' @param ... Ignored.
#' @return Depends on mode. In-sample: `fitted()` result. `"state"`: a
#'   data.frame with `row`, `psi` (marginalised posterior mean), `psi_lower` /
#'   `psi_upper` (equal-tailed 95% credible interval; `NA` when the engine did
#'   not return per-cell predictive variance), and `heldout`. Design-matrix/
#'   terms: data.frame with estimate and CIs. `occu_cover`: a `tobs_prediction`
#'   table (one row per cell) with per-unit draw matrices in `attr(, "draws")`.
#' @export
predict.tobs_fit <- function(object, X.0 = NULL,
                                 type = c("occupancy", "detection", "both",
                                          "state"),
                                 quantiles = c(0.025, 0.5, 0.975),
                                 terms = NULL, n_points = 50L,
                                 newdata = NULL, times = NULL, level = 0.95,
                                 nsim = 1000L, draws = TRUE, time_col = NULL,
                                 ...) {
  # `terms` names one design column, and every family reading it does so
  # through this generic, so validate here rather than in each fitter.
  if (!is.null(terms)) terms <- .tobs_predict_term(terms)
  # Several families return below without ever reading `terms`; say so instead
  # of handing back an in-sample fit that ignored it.
  if (!is.null(terms)) {
    mt <- object$model$model_type %||% ""
    if (mt %in% .TOBS_PREDICT_NO_TERMS)
      stop(sprintf(paste0("predict(terms = ) is not supported for model type ",
                          "'%s': its predictor varies covariates through %s. ",
                          "Build the grid you want and pass it there."),
                   mt,
                   if (mt %in% names(.TOBS_PREDICT_NEWDATA_TYPE)) "`newdata`"
                   else "`X.0`"),
           call. = FALSE)
  }
  # N-mixture abundance: the response types are "abundance" / "detection", so
  # route before the occupancy-specific match.arg(type) rejects them.
  if (identical(object$model$model_type, "nmix") ||
      identical(object$model$model_type, "removal")) {
    nmix_type <- if (missing(type) || length(type) > 1L) "abundance" else type
    return(.tobs_predict_nmix(object, X.0 = X.0, type = nmix_type,
                              quantiles = quantiles, terms = terms,
                              n_points = n_points))
  }
  # Distance sampling: response types are "lambda" (abundance / density) and
  # "sigma" (detection scale); route before the occupancy match.arg(type).
  if (identical(object$model$model_type, "distance")) {
    dist_type <- if (missing(type) || length(type) > 1L) "lambda" else type
    return(.tobs_predict_distance(object, X.0 = X.0, type = dist_type))
  }
  # False-positive occupancy: response types are "psi" (occupancy) and "p11"
  # (true detection); route before the standard occupancy match.arg(type).
  if (identical(object$model$model_type, "fp_occu")) {
    fp_type <- if (missing(type) || length(type) > 1L) "psi" else type
    if (identical(fp_type, "occupancy")) fp_type <- "psi"
    if (identical(fp_type, "detection")) fp_type <- "p11"
    return(.tobs_predict_fp_occu(object, X.0 = X.0, type = fp_type))
  }
  # Open N-mixture: response types are "lambda" (initial abundance) and "gamma"
  # (recruitment); route before the standard occupancy match.arg(type).
  if (identical(object$model$model_type, "dyn_abun")) {
    da_type <- if (missing(type) || length(type) > 1L) "lambda" else type
    if (identical(da_type, "abundance")) da_type <- "lambda"
    return(.tobs_predict_dyn_abun(object, X.0 = X.0, type = da_type))
  }
  # Families predicting from `newdata` with their own response types, all
  # through `.tobs_predict_<model_type>(object, newdata =, type =)`. Each
  # names the type it reports when the caller did not choose one.
  nd_default <- .TOBS_PREDICT_NEWDATA_TYPE[object$model$model_type %||% ""]
  if (!is.na(nd_default)) {
    fn <- .tobs_s3_handler(
      "predict", object$model$model_type,
      c(ms_occu = "ms_community", ms_dyn_occu = "ms_community",
        ms_int_occu = "ms_community"))
    return(fn(object,
              newdata = .tobs_resolve_newdata(newdata, X.0),
              type = if (missing(type) || length(type) > 1L) unname(nd_default)
                     else type))
  }
  # occu_cover joint fit: the response types are occurrence / cover_cond /
  # cover_exp / change, so route before the occupancy match.arg(type) rejects
  # them. `newdata` (or the positional `X.0` when a data.frame) carries the
  # prediction units; `times` drives the change map. See ?predict.tobs_fit.
  if (identical(object$model$model_type, "occu_cover")) {
    oc_type <- if (missing(type) || length(type) > 1L) "occurrence" else type
    nd <- newdata
    if (is.null(nd) && is.data.frame(X.0)) nd <- X.0
    return(.tobs_predict_joint(object, newdata = nd, type = oc_type,
                               times = times, level = level, nsim = nsim,
                               draws = draws, time_col = time_col))
  }
  if (identical(object$model$model_type, "occu_multiscale_cover")) {
    if (!is.null(X.0) || !is.null(terms) || !is.null(newdata)) {
      stop("predict() for an occu_multiscale_cover() fit returns the in-sample ",
           "per-unit posterior; the areal field is tied to the cell graph, so ",
           "X.0 / newdata / terms prediction is not supported. Call ",
           "predict(fit, type = \"state\" / \"availability\" / \"detection\" / ",
           "\"cover\").", call. = FALSE)
    }
    oms_type <- if (missing(type) || length(type) > 1L) "state" else type
    return(.tobs_predict_occu_multiscale_cover(object, oms_type))
  }
  if (identical(object$model$model_type, "ms_occu_cover_spatial")) {
    # The latent fields are tied to the cell graph, so prediction is the
    # in-sample per-species per-cell posterior (calibrated psi / cover + a
    # central interval, marginalised over the loading + field posterior). New
    # covariate / cell prediction is not supported (no field at an unseen cell).
    if (!is.null(X.0) || !is.null(terms) || !is.null(newdata)) {
      stop("predict() for a spatial-factor ms_occu_cover() fit returns the ",
           "in-sample per-species per-cell posterior; the latent fields are ",
           "tied to the cell graph, so X.0 / newdata / terms prediction is not ",
           "supported. Call predict(fit, type = \"occupancy\" / \"cover_cond\" ",
           "/ \"cover_exp\").", call. = FALSE)
    }
    oc_type <- if (missing(type) || length(type) > 1L) "occupancy" else type
    return(.tobs_ms_ocs_predict_state(object, oc_type))
  }
  # The count GLMMs have a single response, so they take no `type`.
  # (`ms_count` covers jsdm() too -- one model class, bernoulli response.)
  if (object$model$model_type %in% c("ms_count", "count")) {
    fn <- .tobs_s3_handler("predict", object$model$model_type)
    return(fn(object, newdata = .tobs_resolve_newdata(newdata, X.0)))
  }
  # Standalone occu() SVC fit rerouted through the joint direct-grid engine
  # (gcol33/tulpaObs#81): the occupancy psi / detection p / per-cell change carry
  # the shared areal field, so route field-aware prediction through the joint
  # substrate (the occupancy-only twin of the occu_cover joint predict). `newdata`
  # (or the positional `X.0` when a data.frame) carries the prediction units;
  # `times` drives the change map.
  if (isTRUE(object$occu_only_joint)) {
    oc_type <- if (missing(type) || length(type) > 1L) "occupancy" else type
    nd <- newdata
    if (is.null(nd) && is.data.frame(X.0)) nd <- X.0
    return(.tobs_predict_occu_joint(object, newdata = nd, type = oc_type,
                                    times = times, level = level, nsim = nsim,
                                    draws = draws, time_col = time_col))
  }
  type <- match.arg(type)

  # State posterior / NA-response prediction (nested-Laplace only).
  if (type == "state") {
    sp <- object$state_posterior
    if (!is.null(sp)) return(sp)
    if (identical(object$method, "nested_laplace")) {
      stop("predict(type = \"state\") is unavailable for this nested-Laplace ",
           "fit: the engine returned no per-cell fitted predictor and the ",
           "latent field is not reconstructable from the modes alone (a mixed-",
           "scale prior such as bym2 on an older tulpa). Reinstall tulpa so ",
           "`tulpa_nested_laplace()` returns `fitted_eta`.", call. = FALSE)
    }
    stop("predict(type = \"state\") is available for method = ",
         "\"nested_laplace\" fits only (it reads the marginalised per-site psi ",
         "posterior the nested path stores).", call. = FALSE)
  }

  # In-sample mode
  if (is.null(X.0) && is.null(terms)) return(fitted(object))

  model <- object$model
  draws <- object$draws
  pi_list <- model$process_info

  # Terms-based mode
  if (!is.null(terms)) {
    return(predict_terms(object, terms, type, quantiles, n_points))
  }

  # Design-matrix mode
  n_pred <- nrow(X.0)
  n_draws <- nrow(draws)
  p_occ <- pi_list[[1]]$p

  if (ncol(X.0) != p_occ) {
    stop(sprintf("X.0 has %d columns but model has %d occupancy coefficients",
                 ncol(X.0), p_occ))
  }

  # Occupancy probability draws at the new design (shared with the ensemble
  # stacked-predictive path; see .tobs_psi_draws()).
  psi_draws <- .tobs_psi_draws(draws, X.0, p_occ)

  data.frame(
    mean = colMeans(psi_draws),
    sd = apply(psi_draws, 2, sd),
    q2.5 = apply(psi_draws, 2, quantile, quantiles[1]),
    q50 = apply(psi_draws, 2, quantile, quantiles[2]),
    q97.5 = apply(psi_draws, 2, quantile, quantiles[3])
  )
}

# `terms` names ONE column of the fitted design matrix to vary; every other
# column is held at its column mean. A second term could only be honoured the
# same way, at its own mean, which is a grid over one covariate reported as a
# grid over two. Grouping it instead (the ggeffects convention) is not
# expressible here: a factor's levels span several design columns, and this
# path sees column names, not model-frame variables. So a longer vector is
# rejected, and `X.0` carries any grid over more than one covariate.
.tobs_predict_term <- function(terms) {
  if (is.character(terms) && length(terms) == 1L && !is.na(terms)) return(terms)
  got <- if (!is.character(terms))
           paste0("a ", class(terms)[1L], " of length ", length(terms))
         else if (length(terms) == 0L) "an empty character vector"
         else paste0("[", paste(terms, collapse = ", "), "]")
  stop("predict(terms = ) varies one term: the name of a single column of ",
       "the predicted process's design matrix, with every other column held ",
       "at its mean. Got ", got, ". For a grid over more than one covariate, ",
       "build the design matrix yourself and pass it as `X.0` (one row per ",
       "prediction unit).", call. = FALSE)
}

# Terms-based prediction: the response over the range of one design column.
predict_terms <- function(object, terms, type, quantiles, n_points) {
  model <- object$model
  draws <- object$draws
  pi_list <- model$process_info

  proc_idx <- if (type == "detection") 2 else 1
  p_proc <- pi_list[[proc_idx]]$p
  X_orig <- model$X_processes[[proc_idx]]
  coef_names <- pi_list[[proc_idx]]$coef_names
  beta_offset <- if (proc_idx > 1) sum(vapply(pi_list[1:(proc_idx-1)], function(pi) pi$p, integer(1))) else 0

  # Also reached directly from tobs_marginal_effect(), so validate here too.
  term_var <- .tobs_predict_term(terms)
  col_idx <- match(term_var, coef_names)
  if (is.na(col_idx)) {
    stop(sprintf("term '%s' not found in %s coefficients: %s",
                 term_var, pi_list[[proc_idx]]$name,
                 paste(coef_names, collapse = ", ")))
  }

  # Create prediction grid: vary term_var, hold others at mean
  x_range <- range(X_orig[, col_idx])
  x_grid <- seq(x_range[1], x_range[2], length.out = n_points)

  X_pred <- matrix(colMeans(X_orig), nrow = n_points, ncol = p_proc, byrow = TRUE)
  X_pred[, col_idx] <- x_grid

  # Predict from each draw, on the process's response scale (logit -> probability
  # for occupancy / detection, log -> intensity for abundance).
  link <- pi_list[[proc_idx]]$link %||% "logit"
  inv_link <- switch(link, logit = plogis, log = exp, identity = identity,
                     stop("predict_terms: unsupported link '", link, "'.",
                          call. = FALSE))
  n_draws <- nrow(draws)
  pred_draws <- matrix(NA_real_, n_draws, n_points)
  for (s in seq_len(n_draws)) {
    beta <- draws[s, beta_offset + seq_len(p_proc)]
    pred_draws[s, ] <- inv_link(as.vector(X_pred %*% beta))
  }

  result <- data.frame(
    x = x_grid,
    estimate = colMeans(pred_draws),
    lower = apply(pred_draws, 2, quantile, quantiles[1]),
    upper = apply(pred_draws, 2, quantile, quantiles[3])
  )
  attr(result, "term") <- term_var
  attr(result, "process") <- pi_list[[proc_idx]]$name
  class(result) <- c("tobs_prediction", "data.frame")
  result
}

#' @export
plot.tobs_prediction <- function(x, ...) {
  # occu_cover predictions are objects-only (a tidy table + per-unit draw
  # matrices in attr "draws"); the real map is one join away in the user's own
  # ggplot/sf. See ?predict.tobs_fit.
  if (!is.null(attr(x, "quantity"))) {
    q <- attr(x, "quantity")
    message("predict(occu_cover) returns a table (one row per cell) plus ",
            "per-unit draw matrices in attr(x, \"draws\"); join it to your ",
            "spatial geometry and map it yourself (geom_tile / geom_sf). ",
            "Quantity: ", q, ".")
    return(invisible(x))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    plot(x$x, x$estimate, type = "l", ylim = range(c(x$lower, x$upper)),
         xlab = attr(x, "term"), ylab = attr(x, "process"),
         main = sprintf("Effect of %s on %s", attr(x, "term"), attr(x, "process")))
    polygon(c(x$x, rev(x$x)), c(x$lower, rev(x$upper)),
            col = rgb(0, 0, 0, 0.1), border = NA)
    return(invisible(x))
  }
  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data$x, y = .data$estimate)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), alpha = 0.2) +
    ggplot2::geom_line() +
    ggplot2::labs(x = attr(x, "term"), y = attr(x, "process"))
  print(p)
  invisible(x)
}

#' Compute marginal effect of a covariate
#'
#' Returns the fitted response over the range of one covariate, holding the
#' others at their means, on the process's response scale (occupancy /
#' detection probability, or abundance intensity).
#' @param object A `tobs_fit` object.
#' @param covariate Name of covariate.
#' @param process `"occupancy"` (default), `"detection"`, or `"abundance"` (the
#'   state process of a count family, on the intensity scale).
#' @param n_points Number of prediction points (default 100).
#' @return A data.frame with covariate values and the predicted response.
#' @examples
#' \donttest{
#' sim <- simulate_abun(N = 100, J = 4, n_abund_covs = 2, n_det_covs = 1,
#'                      seed = 1)
#' fit <- tobs(~ abund_cov1 + abund_cov2, data = sim$data,
#'             family = abun(K_max = 50), detection = ~ det_cov1, y = sim$y,
#'             control = list(verbose = FALSE))
#' me <- tobs_marginal_effect(fit, "abund_cov1", process = "abundance")
#' head(me)   # estimate is on the abundance (lambda) scale, not a probability
#' }
#' @export
tobs_marginal_effect <- function(object, covariate,
                                 process = c("occupancy", "detection",
                                             "abundance"),
                                 n_points = 100L) {
  process <- match.arg(process)
  type <- switch(process,
                 detection = "detection",
                 abundance = "abundance",   # state process 1 on the log scale
                 "occupancy")
  predict_terms(object, terms = covariate, type = type,
                quantiles = c(0.025, 0.5, 0.975), n_points = n_points)
}

#' Estimate species richness from community model
#'
#' Per-site expected species richness `sum_s psi_{s,i}` with a posterior
#' credible interval propagated from the community-mean occupancy draws. For the
#' dynamic community family this is the first-season richness (from `psi1`).
#' @param object A `tobs_fit` object from a community occupancy model
#'   (`ms_occu()`, `ms_dyn_occu()`, or `ms_int_occu()`).
#' @return A data.frame with site-level richness estimates.
#' @examples
#' \donttest{
#' sim <- simulate_ms_occu(N = 40, J = 3, n_species = 5, seed = 1)
#' fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
#'             y = sim$y, species = paste0("sp", 1:5),
#'             control = list(verbose = FALSE))
#' head(tobs_richness(fit))
#' }
#' @export
tobs_richness <- function(object) {
  if (!(object$model$model_type %in%
        c("ms_occu", "ms_dyn_occu", "ms_int_occu"))) {
    stop("tobs_richness() requires a community occupancy fit (ms_occu(), ",
         "ms_dyn_occu(), or ms_int_occu()).", call. = FALSE)
  }
  .tobs_richness_community(object)
}

# tidy, glance, ranef inherited from tulpa::tulpa_fit

#' Update and refit an occupancy model
#' @param object A `tobs_fit` object.
#' @param ... Named arguments to override in the underlying fit.
#' @param evaluate If TRUE (default), refit the model.
#' @return Updated `tobs_fit` object (or call if `evaluate = FALSE`).
#' @export
update.tobs_fit <- function(object, ..., evaluate = TRUE) {
  # Structured terms (spatial / temporal / re / svc / latent) travel with the
  # model via `model$structured_terms`, so refitting only needs the model plus
  # any overridden fit controls.
  args <- list(model = object$model)
  dots <- list(...)
  for (nm in names(dots)) args[[nm]] <- dots[[nm]]

  if (!evaluate) return(args)
  do.call(.tobs_fit_model, args)
}

#' Check model identifiability
#'
#' Diagnostics for potential identifiability issues in occupancy models.
#' Checks for: confounded covariates, low detection rates, sparse data.
#'
#' @param model A `tobs_model` object (before fitting).
#' @param fit Optional `tobs_fit` object (for post-fit diagnostics).
#' @return A list with diagnostic messages and flags.
#' @export
tobs_check_id <- function(model, fit = NULL) {
  issues <- character()

  if (!inherits(model, "tobs_model")) {
    stop("model must be a tobs_model object")
  }

  # Pre-fit checks
  if (model$model_type %in% c("single", "community")) {
    y <- model$y
    n_sites <- nrow(y)
    max_visits <- ncol(y)

    # Naive occupancy
    n_detected <- sum(apply(y, 1, function(row) any(row[row >= 0] == 1)))
    naive_occ <- n_detected / n_sites
    if (naive_occ < 0.05) {
      issues <- c(issues, sprintf("Very low naive occupancy (%.1f%%). Model may struggle to estimate occupancy coefficients.", 100 * naive_occ))
    }
    if (naive_occ > 0.95) {
      issues <- c(issues, sprintf("Very high naive occupancy (%.1f%%). Little information to estimate occupancy effects.", 100 * naive_occ))
    }

    # Mean visits per site
    n_visits <- apply(y, 1, function(row) sum(row >= 0))
    if (mean(n_visits) < 2) {
      issues <- c(issues, sprintf("Very few visits per site (mean %.1f). Detection and occupancy may be confounded.", mean(n_visits)))
    }

    # Check for collinearity in occupancy covariates
    X_occ <- model$X_processes[[1]]
    if (ncol(X_occ) > 2) {
      cors <- cor(X_occ[, -1, drop = FALSE])
      high_cor <- which(abs(cors) > 0.8 & upper.tri(cors), arr.ind = TRUE)
      if (nrow(high_cor) > 0) {
        names_occ <- colnames(X_occ)[-1]
        for (k in seq_len(nrow(high_cor))) {
          issues <- c(issues, sprintf("High correlation (%.2f) between %s and %s.",
                                      cors[high_cor[k, 1], high_cor[k, 2]],
                                      names_occ[high_cor[k, 1]],
                                      names_occ[high_cor[k, 2]]))
        }
      }
    }
  }

  # Post-fit checks (NUTS only; Laplace fits carry NA sampler diagnostics)
  if (!is.null(fit) && inherits(fit, "tobs_fit") && identical(fit$method, "nuts")) {
    if (sum(fit$divergent) > 0) {
      issues <- c(issues, sprintf("%d divergent transitions. Consider increasing adapt_delta or reparameterizing.", sum(fit$divergent)))
    }
    if (mean(fit$accept_prob) < 0.5) {
      issues <- c(issues, sprintf("Low mean acceptance probability (%.2f). Model may be poorly specified.", mean(fit$accept_prob)))
    }
  }

  result <- list(
    identifiable = length(issues) == 0,
    issues = issues
  )
  if (length(issues) > 0) {
    for (msg in issues) message("- ", msg)
  } else {
    message("No identifiability issues detected.")
  }
  invisible(result)
}

# ============================================================================
# spOccupancy $ compatibility accessor
# ============================================================================

#' Access spOccupancy-compatible fields from tobs fits
#'
#' Allows accessing spOccupancy-style fields (e.g., `$beta.samples`,
#' `$psi.samples`) on tobs_fit objects. Since tobs stores actual posterior
#' draws, this is a thin remapping layer.
#'
#' @param x A `tobs_fit` object.
#' @param name Field name to access.
#' @return The requested field value.
#' @export
`$.tobs_fit` <- function(x, name) {
  # First check native fields
  val <- .subset2(x, name)
  if (!is.null(val)) return(val)

  model <- .subset2(x, "model")
  draws <- .subset2(x, "draws")
  pi_list <- if (!is.null(model)) model$process_info else NULL

  switch(name,
    # Occupancy fixed effect draws
    "beta.samples" = {
      if (is.null(pi_list)) return(NULL)
      p_occ <- pi_list[[1]]$p
      draws[, seq_len(p_occ), drop = FALSE]
    },

    # Detection fixed effect draws
    "alpha.samples" = {
      if (is.null(pi_list) || length(pi_list) < 2) return(NULL)
      p_occ <- pi_list[[1]]$p
      p_det <- pi_list[[2]]$p
      draws[, p_occ + seq_len(p_det), drop = FALSE]
    },

    # Occupancy probabilities (n_draws x n_sites), recomputed from the draws so
    # they carry the posterior uncertainty. Only families whose first process IS
    # occupancy have them; an abundance arm (lambda) yields NULL.
    "psi.samples" = {
      if (is.null(pi_list) || !identical(pi_list[[1]]$name, "psi")) return(NULL)
      .tobs_psi_draws(draws, model$X_processes[[1]], pi_list[[1]]$p)
    },

    # Latent occupancy state
    "z.samples" = {
      psi <- x$psi.samples
      if (!is.null(psi)) {
        matrix(rbinom(length(psi), 1, psi), nrow = nrow(psi))
      }
    },

    # Detection probabilities: the second arm's coefficients start after the
    # first arm's block in each draw row.
    "p.samples" = {
      if (is.null(pi_list) || length(pi_list) < 2) return(NULL)
      .tobs_psi_draws(draws, model$X_processes[[2]], pi_list[[2]]$p,
                      offset = pi_list[[1]]$p)
    },

    # Computation time
    "run.time" = .subset2(x, "elapsed"),

    # Not a compat field
    NULL
  )
}


# ============================================================================
# Spatial prediction at new locations
# ============================================================================

# The fitted spatial field a `tobs_fit` carries, and where its variance comes
# from. Two sources, in order of preference:
#
#   "draws" -- the sampled field columns of `fit$draws`, which the samplers name
#     `spatial_field[i]` on the areal path (gcol33/tulpaObs#142) and `gp_w[i]`
#     on the continuous-GP path. These carry the field's posterior variance, so
#     the predicted spread includes it.
#   "point" -- `fit$spatial_field`, the posterior-mean surface the deterministic
#     backends report. It is a point estimate: it shifts every draw by the same
#     amount, so the predicted spread is the coefficients' alone.
#
# The two prefixes are the ones `src/occu_fit.cpp` actually emits. An SVC term's
# `svc_w[i,j]` columns are deliberately not matched: those are one surface per
# varying coefficient, not a single additive field, so summing them here would
# report a quantity the model does not have.
.tobs_spatial_field_source <- function(object) {
  idx <- grep("^spatial_field\\[|^gp_w\\[", colnames(object$draws))
  if (length(idx))
    return(list(kind = "draws", values = object$draws[, idx, drop = FALSE],
                n_nodes = length(idx)))
  f <- object$spatial_field
  if (!is.null(f) && length(f))
    return(list(kind = "point", values = as.numeric(f), n_nodes = length(f)))
  list(kind = "none", values = NULL, n_nodes = 0L)
}

# Coordinates for the fitted field's nodes, as an n_nodes x 2 matrix, or NULL
# when they cannot be resolved. A caller-supplied `node.coords` wins; otherwise
# the continuous terms (gp / multiscale_gp / svc) carry their own, stored
# flattened row-major by their constructors, so a bare vector is reshaped.
# Areal terms carry none, which is what NULL reports.
.tobs_spatial_node_coords <- function(object, node.coords, n_nodes) {
  co <- node.coords %||% object$spatial$coords
  if (is.null(co) || !length(co)) return(NULL)
  if (is.null(dim(co))) {
    if (length(co) != 2L * n_nodes) return(NULL)
    co <- matrix(as.numeric(co), ncol = 2L, byrow = TRUE)
  }
  co <- as.matrix(co)
  if (ncol(co) < 2L || nrow(co) != n_nodes) return(NULL)
  co
}

#' Predict the state process at new spatial locations
#'
#' Generates state-process predictions at new coordinates, including the fitted
#' spatial field interpolated to those locations by inverse-distance weighting
#' over its five nearest nodes. The returned scale follows the family's
#' state-process link: occupancy probability for occupancy families, abundance
#' intensity (lambda) for the count families.
#'
#' The field is taken from the sampled field columns of `object$draws` when the
#' fit has them (the NUTS paths), so the reported `sd` and quantiles carry the
#' field's own posterior variance. On the deterministic backends
#' (`laplace` / `nested_laplace`) the draws hold only the coefficients, so the
#' posterior-mean surface in `object$spatial_field` is used instead: it enters
#' every draw as the same offset, and the reported spread is then the
#' coefficients' alone.
#'
#' A continuous term (`gp()`, `spde()`, `svc()`) carries the coordinates of its
#' own nodes, so nothing extra is needed. An areal term (`icar()`, `bym2()`,
#' `car()`) has graph nodes and no geometry, so interpolating it to a new point
#' is only defined once the nodes are placed: pass `node.coords`, one row per
#' element of `object$spatial_field`. Without it the call is an error rather
#' than a prediction that quietly drops the field.
#'
#' @param object A `tobs_fit` object fitted with a spatial component.
#' @param newcoords Matrix of new coordinates (n_new x 2).
#' @param newocc.covs Optional data.frame of covariates at new locations.
#' @param quantiles Quantiles for credible intervals (default 0.025, 0.5, 0.975).
#' @param node.coords Optional matrix of coordinates for the fitted field's
#'   nodes (n_nodes x 2), required for an areal field and ignored when the term
#'   already carries its own coordinates.
#' @return A data.frame with `mean`, `sd`, and quantile columns (on the response
#'   scale: occupancy probability or abundance intensity).
#' @export
tobs_predict_spatial <- function(object, newcoords, newocc.covs = NULL,
                                 quantiles = c(0.025, 0.5, 0.975),
                                 node.coords = NULL) {
  if (is.null(object$spatial)) {
    stop("tobs_predict_spatial requires a model fitted with a spatial component", call. = FALSE)
  }

  # Build design matrix at new locations
  if (!is.null(newocc.covs)) {
    # Use model's formula to build X
    model <- object$model
    occ_formula <- model$occ_formula
    if (!is.null(occ_formula)) {
      X.0 <- model.matrix(occ_formula, data = newocc.covs)
    } else {
      X.0 <- as.matrix(cbind(1, newocc.covs))
    }
  } else {
    n_new <- nrow(newcoords)
    p_occ <- object$model$process_info[[1]]$p
    X.0 <- matrix(0, n_new, p_occ)
    X.0[, 1] <- 1  # Intercept only
  }

  # Fixed effect prediction
  draws <- object$draws
  p_occ <- object$model$process_info[[1]]$p
  n_draws <- nrow(draws)
  n_new <- nrow(newcoords)

  eta_draws <- matrix(NA_real_, n_draws, n_new)
  for (s in seq_len(n_draws)) {
    beta <- draws[s, seq_len(p_occ)]
    eta_draws[s, ] <- as.vector(X.0 %*% beta)
  }

  # Interpolate the fitted field to the new locations (IDW over the k = 5
  # nearest nodes). Both the field values and the node coordinates are resolved
  # explicitly, and a field that cannot be placed is an error: silently
  # returning the fixed-effect-only prediction is indistinguishable from a fit
  # whose field is flat (gcol33/tulpaObs#179).
  fld <- .tobs_spatial_field_source(object)
  if (identical(fld$kind, "none"))
    stop("tobs_predict_spatial: the fit declares a spatial component but ",
         "carries no fitted field to interpolate (neither sampled field ",
         "columns in `draws` nor `spatial_field`).", call. = FALSE)

  fit_coords <- .tobs_spatial_node_coords(object, node.coords, fld$n_nodes)
  if (is.null(fit_coords))
    stop(sprintf(paste0("tobs_predict_spatial: the fitted %s field has %d ",
                        "nodes and no coordinates to place them at. An areal ",
                        "field's nodes are graph vertices, so supply ",
                        "`node.coords` (a %d x 2 matrix, one row per element ",
                        "of `fit$spatial_field`)."),
                 object$spatial$type %||% "spatial", fld$n_nodes, fld$n_nodes),
         call. = FALSE)

  # The IDW weights depend only on geometry, so build them once per new
  # location and apply them to every draw at once.
  k <- min(5L, nrow(fit_coords))
  W <- matrix(0, n_new, fld$n_nodes)
  for (i in seq_len(n_new)) {
    dists <- sqrt((fit_coords[, 1] - newcoords[i, 1])^2 +
                  (fit_coords[, 2] - newcoords[i, 2])^2)
    nn <- order(dists)[seq_len(k)]
    w  <- 1 / (dists[nn] + 1e-10)
    W[i, nn] <- w / sum(w)
  }
  eta_draws <- if (identical(fld$kind, "draws")) {
    eta_draws + tcrossprod(fld$values, W)          # per-draw field, real variance
  } else {
    sweep(eta_draws, 2L, as.numeric(W %*% fld$values), "+")  # point surface
  }

  # Apply the state process's inverse link so the returned scale is correct for
  # the family: occupancy (logit -> probability), abundance (log -> intensity).
  # Otherwise a count fit would silently report logit-of-log-lambda.
  link <- object$model$process_info[[1]]$link %||% "logit"
  inv_link <- switch(link,
    logit    = plogis,
    log      = exp,
    identity = identity,
    stop("tobs_predict_spatial: unsupported state-process link '", link, "'.",
         call. = FALSE))
  psi_draws <- inv_link(eta_draws)

  result <- data.frame(
    mean = colMeans(psi_draws),
    sd = apply(psi_draws, 2, sd)
  )
  for (q in quantiles) {
    qname <- paste0("q", gsub("\\.", "", format(q * 100, nsmall = 1)))
    result[[qname]] <- apply(psi_draws, 2, quantile, q)
  }
  result
}
