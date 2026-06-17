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

#' Number of observations
#' @param object A `tobs_fit` object.
#' @param ... Ignored.
#' @return Integer count of non-NA detection history entries.
#' @export
nobs.tobs_fit <- function(object, ...) {
  model <- object$model
  if (model$model_type == "single") {
    y <- model$y
    sum(y >= 0)
  } else if (model$model_type == "dynamic") {
    sum(model$y_flat >= 0)
  } else if (model$model_type == "nmix" || model$model_type == "removal") {
    length(model$y_long)
  } else if (model$model_type == "distance") {
    sum(!is.na(model$y))
  } else if (model$model_type == "fp_occu") {
    length(model$y_long)
  } else if (model$model_type == "dyn_abun") {
    sum(!is.na(model$y))
  } else if (model$model_type == "ms_nmix" ||
             model$model_type == "ms_occu_cover" ||
             model$model_type == "ms_occu_cover_spatial" ||
             model$model_type == "occu_multiscale_cover") {
    sum(!is.na(model$y))
  } else if (model$model_type == "ms_occu" ||
             model$model_type == "ms_dyn_occu") {
    sum(model$valid)
  } else if (model$model_type == "ms_int_occu") {
    sum(vapply(model$valid, sum, integer(1)))
  } else {
    NA_integer_
  }
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
#' Pareto-k diagnostic when present. The joint-coupled families (`occu_cover()`,
#' `occu()` spatial, `occu_multiscale_cover()`) carry the diagnostic at the fit
#' top level (gcol33/tulpaObs#104); every other family glances exactly as before.
#'
#' @param x A fitted `tobs_fit`.
#' @param ... Ignored.
#' @return The base one-row `glance` data frame, with three extra columns on a
#'   joint-coupled fit that requested the diagnostic (`control$diagnose.k = TRUE`,
#'   off by default per gcol33/tulpaObs#101):
#'   \describe{
#'     \item{`pareto_k`}{The outer importance-sampling \eqn{\hat{k}} for the
#'       hyperparameter Gaussian summary; `< 0.7` indicates a reliable summary.
#'       When `pareto_k_is_ess` is `TRUE` this column instead holds the quad-ESS
#'       fallback (the \eqn{\hat{k}} fit declined).}
#'     \item{`pareto_k_is_ess`}{`TRUE` when the `pareto_k` column is the quad-ESS
#'       fallback rather than a fitted \eqn{\hat{k}}.}
#'     \item{`pareto_k_proposal_source`}{How the importance proposal was built
#'       (gcol33/tulpa#116): `"mode_hessian"` from the Laplace curvature at the
#'       hyperparameter mode -- curvature-backed, so the \eqn{\hat{k}} stays
#'       trustworthy even when a sharp posterior collapses the integration grid
#'       to ~1 cell; `"grid_moment"` from the grid-weighted covariance of the
#'       integration nodes -- the regime to watch, since it under-disperses (and
#'       can flag a spurious high \eqn{\hat{k}}) when the grid concentrates.}
#'   }
#' @export
glance.tobs_fit <- function(x, ...) {
  g <- NextMethod()
  # Prefer the promoted top-level fields; fall back to the nested joint object so
  # a fit saved before the promotion (gcol33/tulpaObs#104) still glances its k-hat.
  pk <- .tobs_promote_pareto_k(x) %||% .tobs_promote_pareto_k(.tobs_joint_fit(x))
  if (is.null(pk)) return(g)
  if (!is.null(pk$pareto_k))        g$pareto_k <- pk$pareto_k
  if (!is.null(pk$pareto_k_is_ess)) g$pareto_k_is_ess <- as.logical(pk$pareto_k_is_ess)
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
  if (identical(object$model$model_type, "ms_nmix")) {
    return(.tobs_ranef_ms_nmix(object))
  }
  if (identical(object$model$model_type, "ms_occu")) {
    return(.tobs_ranef_ms_occu(object))
  }
  if (identical(object$model$model_type, "ms_occu_cover")) {
    return(.tobs_ranef_ms_occu_cover(object))
  }
  if (identical(object$model$model_type, "ms_dyn_occu")) {
    return(.tobs_ranef_ms_dyn_occu(object))
  }
  if (identical(object$model$model_type, "ms_int_occu")) {
    return(.tobs_ranef_ms_int_occu(object))
  }
  if (identical(object$model$model_type, "occu_cover") && !is.null(object$re)) {
    # occu_cover() shared-field + per-group RE (gcol33/tulpaObs#56, #102, #103):
    # `fit$re` is a flat list of random-intercept terms, one per arm for a lone
    # term, several for crossed / nested groupings sharing an arm. Stack them into
    # one table with `arm` + `var` (grouping variable) columns; the observation
    # arms carry their grouping `level` labels (the occupancy arm reports
    # 1..n_groups).
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
  if (identical(model$model_type, "nmix") ||
      identical(model$model_type, "removal")) return(.tobs_fitted_nmix(object))
  if (identical(model$model_type, "distance")) return(.tobs_fitted_distance(object))
  if (identical(model$model_type, "fp_occu")) return(.tobs_fitted_fp_occu(object))
  if (identical(model$model_type, "dyn_abun")) return(.tobs_fitted_dyn_abun(object))
  if (identical(model$model_type, "ms_nmix")) return(.tobs_fitted_ms_nmix(object))
  if (identical(model$model_type, "ms_occu")) {
    return(.tobs_fitted_ms_occu(object))
  }
  if (identical(model$model_type, "ms_occu_cover")) {
    return(.tobs_fitted_ms_occu_cover(object))
  }
  if (identical(model$model_type, "ms_occu_cover_spatial")) {
    return(.tobs_fitted_ms_occu_cover_spatial(object))
  }
  if (identical(model$model_type, "ms_dyn_occu")) {
    return(.tobs_fitted_ms_dyn_occu(object))
  }
  if (identical(model$model_type, "ms_int_occu")) {
    return(.tobs_fitted_ms_int_occu(object))
  }
  if (identical(model$model_type, "occu_multiscale_cover")) {
    return(.tobs_fitted_occu_multiscale_cover(object))
  }
  means <- object$means
  pi_list <- model$process_info

  # Extract occupancy and detection linear predictors at posterior mean
  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  beta_occ <- means[seq_len(pi_list[[1]]$p)]
  beta_det <- means[pi_list[[1]]$p + seq_len(pi_list[[2]]$p)]

  eta_occ <- as.vector(X_occ %*% beta_occ)
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
# Returns an [n_sites x n_seasons] matrix. Detection / colonization / extinction
# are site-level (constant across seasons), matching the engine.
.tobs_dynamic_smoothed_z <- function(model, means, pi_list) {
  n_sites   <- model$n_sites
  T_seasons <- model$n_seasons

  X <- model$X_processes
  off <- cumsum(c(0L, vapply(pi_list, function(pp) pp$p, integer(1))))
  beta_psi1 <- means[off[1] + seq_len(pi_list[[1]]$p)]
  beta_p    <- means[off[2] + seq_len(pi_list[[2]]$p)]
  beta_gam  <- means[off[3] + seq_len(pi_list[[3]]$p)]
  beta_eps  <- means[off[4] + seq_len(pi_list[[4]]$p)]

  psi1  <- plogis(as.vector(X[[1]] %*% beta_psi1))
  p     <- plogis(as.vector(X[[2]] %*% beta_p))
  gamma <- plogis(as.vector(X[[3]] %*% beta_gam))
  eps   <- plogis(as.vector(X[[4]] %*% beta_eps))

  y <- model$y  # [n_sites x max_visits x n_seasons]
  z <- matrix(NA_real_, n_sites, T_seasons)

  for (i in seq_len(n_sites)) {
    pi_i <- p[i]; gam_i <- gamma[i]; eps_i <- eps[i]

    # Per-season emission likelihood under each state, and a hard-detection mask.
    em <- matrix(1, T_seasons, 2L)  # columns: state 0 (unocc), state 1 (occ)
    for (t in seq_len(T_seasons)) {
      raw <- y[i, , t]
      raw <- raw[!is.na(raw) & raw >= 0]
      if (length(raw) == 0L) {
        em[t, ] <- c(1, 1)  # no visits: uninformative
      } else if (any(raw == 1)) {
        # A detection rules out the unoccupied state.
        em[t, 1L] <- 0
        em[t, 2L] <- prod(p[i]^raw * (1 - p[i])^(1 - raw))
      } else {
        em[t, 1L] <- 1                       # unoccupied -> all non-detections
        em[t, 2L] <- prod(1 - p[i])^length(raw)
      }
    }

    # Transition matrix Tr[a, b] = P(z_{t+1}=b-1 | z_t=a-1).
    Tr <- matrix(c(1 - gam_i, eps_i,
                   gam_i,     1 - eps_i), 2L, 2L)

    # Forward filtering (scaled).
    fwd <- matrix(0, T_seasons, 2L)
    prior <- c(1 - psi1[i], psi1[i])
    a <- prior * em[1L, ]
    a <- a / sum(a)
    fwd[1L, ] <- a
    if (T_seasons > 1L) {
      for (t in 2L:T_seasons) {
        pred <- as.vector(t(Tr) %*% fwd[t - 1L, ])
        a <- pred * em[t, ]
        a <- a / sum(a)
        fwd[t, ] <- a
      }
    }

    # Backward smoothing (Rauch-Tung-Striebel style for discrete HMM).
    sm <- matrix(0, T_seasons, 2L)
    sm[T_seasons, ] <- fwd[T_seasons, ]
    if (T_seasons > 1L) {
      for (t in (T_seasons - 1L):1L) {
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

# Occupancy-probability draws at a design matrix X.0: plogis(X.0 %*% beta_occ)
# for every posterior draw, returned as [n_draws x nrow(X.0)]. Used by
# predict.tobs_fit (design-matrix mode) and predict.tobs_stack (stacked
# predictive); kept in one place so the two share the same parameterization.
.tobs_psi_draws <- function(draws, X.0, p_occ) {
  beta <- draws[, seq_len(p_occ), drop = FALSE]
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
  if (identical(object$model$model_type, "nmix")) {
    return(.tobs_residuals_nmix(object, type))
  }
  if (identical(object$model$model_type, "removal")) {
    return(.tobs_residuals_removal(object, type))
  }
  if (identical(object$model$model_type, "distance")) {
    return(.tobs_residuals_distance(object, type))
  }
  if (identical(object$model$model_type, "fp_occu")) {
    return(.tobs_residuals_fp_occu(object, type))
  }
  if (identical(object$model$model_type, "dyn_abun")) {
    return(.tobs_residuals_dyn_abun(object, type))
  }
  if (object$model$model_type %in% c("ms_occu", "ms_dyn_occu", "ms_int_occu")) {
    stop(sprintf(paste0("residuals() is not defined for %s() community fits; ",
         "use fitted() / ranef() / coef()."), object$model$model_type),
         call. = FALSE)
  }
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
  if (identical(model$model_type, "nmix")) {
    return(.tobs_simulate_nmix(object, nsim))
  }
  if (identical(model$model_type, "removal")) {
    return(.tobs_simulate_removal(object, nsim))
  }
  if (identical(model$model_type, "distance")) {
    return(.tobs_simulate_distance(object, nsim))
  }
  if (identical(model$model_type, "fp_occu")) {
    return(.tobs_simulate_fp_occu(object, nsim))
  }
  if (identical(model$model_type, "dyn_abun")) {
    return(.tobs_simulate_dyn_abun(object, nsim))
  }
  if (identical(model$model_type, "ms_nmix")) {
    return(.tobs_simulate_ms_nmix(object, nsim))
  }
  if (identical(model$model_type, "ms_occu")) {
    return(.tobs_simulate_ms_occu(object, nsim))
  }
  if (identical(model$model_type, "ms_occu_cover")) {
    return(.tobs_simulate_ms_occu_cover(object, nsim))
  }
  if (identical(model$model_type, "ms_dyn_occu")) {
    return(.tobs_simulate_ms_dyn_occu(object, nsim))
  }
  if (identical(model$model_type, "ms_int_occu")) {
    return(.tobs_simulate_ms_int_occu(object, nsim))
  }
  if (identical(model$model_type, "ms_occu_cover_spatial")) {
    return(.tobs_simulate_ms_occu_cover_spatial(object, nsim))
  }
  draws <- object$draws
  n_samples <- nrow(draws)
  pi_list <- model$process_info

  if (model$model_type != "single") {
    stop("simulate() currently only supports single-season models")
  }

  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  n_sites <- model$n_sites
  max_visits <- model$max_visits

  result <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    # Sample a posterior draw
    draw_idx <- sample.int(n_samples, 1)
    beta_occ <- draws[draw_idx, seq_len(pi_list[[1]]$p)]
    beta_det <- draws[draw_idx, pi_list[[1]]$p + seq_len(pi_list[[2]]$p)]

    psi <- plogis(as.vector(X_occ %*% beta_occ))
    p <- plogis(as.vector(X_det %*% beta_det))
    z <- rbinom(n_sites, 1, psi)

    y_sim <- matrix(NA_integer_, n_sites, max_visits)
    for (i in seq_len(n_sites)) {
      # Respect original visit structure (NA pattern)
      for (j in seq_len(max_visits)) {
        if (model$y[i, j] >= 0) {
          y_sim[i, j] <- rbinom(1, 1, z[i] * p[i])
        }
      }
    }
    result[[s]] <- y_sim
  }
  if (nsim == 1) result[[1]] else result
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
#'   decomposition, and `.lwr` / `.upr` at `level`). The result is a
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
#' @param terms Character vector of terms to vary (ggpredict-style).
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
  if (object$model$model_type %in% c("ms_occu", "ms_dyn_occu", "ms_int_occu")) {
    stop(sprintf(paste0("predict() is not yet implemented for %s(). Use ",
         "fitted() for per-species occupancy / detection probabilities, ",
         "coef() / summary() for the community means, and ranef() for the ",
         "per-species deviations."), object$model$model_type), call. = FALSE)
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

# Terms-based prediction (ggpredict-style)
predict_terms <- function(object, terms, type, quantiles, n_points) {
  model <- object$model
  draws <- object$draws
  pi_list <- model$process_info

  proc_idx <- if (type == "detection") 2 else 1
  p_proc <- pi_list[[proc_idx]]$p
  X_orig <- model$X_processes[[proc_idx]]
  coef_names <- pi_list[[proc_idx]]$coef_names
  beta_offset <- if (proc_idx > 1) sum(vapply(pi_list[1:(proc_idx-1)], function(pi) pi$p, integer(1))) else 0

  # Parse the first term (simple: just a variable name for now)
  term_var <- terms[1]
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

  # Predict from each draw
  n_draws <- nrow(draws)
  pred_draws <- matrix(NA_real_, n_draws, n_points)
  for (s in seq_len(n_draws)) {
    beta <- draws[s, beta_offset + seq_len(p_proc)]
    pred_draws[s, ] <- plogis(as.vector(X_pred %*% beta))
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
#' @param object A `tobs_fit` object.
#' @param covariate Name of covariate.
#' @param process `"occupancy"` (default) or `"detection"`.
#' @param n_points Number of prediction points (default 100).
#' @return A data.frame with covariate values and predicted probabilities.
#' @export
tobs_marginal_effect <- function(object, covariate,
                                 process = c("occupancy", "detection"),
                                 n_points = 100L) {
  process <- match.arg(process)
  type <- if (process == "detection") "detection" else "occupancy"
  predict_terms(object, terms = covariate, type = type,
                quantiles = c(0.025, 0.5, 0.975), n_points = n_points)
}

#' Estimate species richness from community model
#' @param object A `tobs_fit` object from a community model.
#' @return A data.frame with site-level richness estimates.
#' @export
tobs_richness <- function(object) {
  if (!identical(object$model$model_type, "ms_occu")) {
    stop("tobs_richness() requires an ms_occu() community fit.", call. = FALSE)
  }
  .tobs_richness_ms_occu(object)
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

    # Occupancy probabilities (n_draws x n_sites)
    "psi.samples" = {
      fit_vals <- fitted(x)
      if (!is.null(fit_vals$psi)) {
        # Recompute from draws for proper uncertainty
        if (is.null(pi_list)) return(NULL)
        p_occ <- pi_list[[1]]$p
        X_occ <- model$X_processes[[1]]
        n_draws <- nrow(draws)
        psi_mat <- matrix(NA_real_, n_draws, nrow(X_occ))
        for (s in seq_len(n_draws)) {
          beta <- draws[s, seq_len(p_occ)]
          psi_mat[s, ] <- plogis(as.vector(X_occ %*% beta))
        }
        psi_mat
      }
    },

    # Latent occupancy state
    "z.samples" = {
      psi <- x$psi.samples
      if (!is.null(psi)) {
        matrix(rbinom(length(psi), 1, psi), nrow = nrow(psi))
      }
    },

    # Detection probabilities
    "p.samples" = {
      if (is.null(pi_list) || length(pi_list) < 2) return(NULL)
      p_occ <- pi_list[[1]]$p
      p_det <- pi_list[[2]]$p
      X_det <- model$X_processes[[2]]
      n_draws <- nrow(draws)
      p_mat <- matrix(NA_real_, n_draws, nrow(X_det))
      for (s in seq_len(n_draws)) {
        alpha <- draws[s, p_occ + seq_len(p_det)]
        p_mat[s, ] <- plogis(as.vector(X_det %*% alpha))
      }
      p_mat
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

#' Predict occupancy at new spatial locations
#'
#' Generates occupancy predictions at new coordinates, including the
#' spatial random effect interpolated from the fitted field.
#'
#' @param object A `tobs_fit` object fitted with a spatial component.
#' @param newcoords Matrix of new coordinates (n_new x 2).
#' @param newocc.covs Optional data.frame of covariates at new locations.
#' @param quantiles Quantiles for credible intervals (default 0.025, 0.5, 0.975).
#' @return A data.frame with `mean`, `sd`, and quantile columns.
#' @export
tobs_predict_spatial <- function(object, newcoords, newocc.covs = NULL,
                                 quantiles = c(0.025, 0.5, 0.975)) {
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

  # Interpolate spatial field to new locations using nearest-neighbor
  sp_type <- object$spatial$type
  cn <- colnames(draws)
  sp_cols <- grep("^phi_spatial\\[|^w_gp\\[|^gp_w\\[", cn)

  if (length(sp_cols) > 0 && !is.null(object$spatial$coords)) {
    fit_coords <- object$spatial$coords
    n_fit <- nrow(fit_coords)

    # Compute distances from new points to fitted points
    for (s in seq_len(n_draws)) {
      sp_effects <- draws[s, sp_cols]
      # Nearest-neighbor interpolation
      for (i in seq_len(n_new)) {
        dists <- sqrt((fit_coords[, 1] - newcoords[i, 1])^2 +
                       (fit_coords[, 2] - newcoords[i, 2])^2)
        # IDW with k=5 nearest neighbors
        k <- min(5, n_fit)
        nn <- order(dists)[seq_len(k)]
        w <- 1 / (dists[nn] + 1e-10)
        w <- w / sum(w)
        eta_draws[s, i] <- eta_draws[s, i] + sum(w * sp_effects[nn])
      }
    }
  }

  psi_draws <- plogis(eta_draws)

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
