# =============================================================================
# stacking.R -- LOO-stacked ensembles of occupancy fits.
#
# tobs_stack() combines several fitted models into a single predictive
# distribution using leave-one-out predictive stacking (Yao, Vehtari, Simpson &
# Gelman 2018), with the leave-one-out densities approximated by Pareto-smoothed
# importance sampling (Vehtari, Gelman & Gabry 2017). PSIS and the stacking
# weight optimization are the canonical `loo` package routines; tulpaObs only
# supplies the per-model pointwise log-likelihood (.tobs_pointwise_loglik) and
# the prediction-combining methods.
#
# Members must be fit to the same observations (the weights are chosen to
# maximize the stacked leave-one-out predictive density over those points).
# Stacking pays off across *different* models / methods; stacking seed-variants
# of one model yields near-uniform weights (the members are statistically
# identical), so the `n.seeds` sugar in tobs() is mainly a Monte-Carlo
# stability / robustness device, while a genuine model average comes from
# passing distinct specifications here.
# =============================================================================

#' Stack fitted models into a LOO-weighted ensemble
#'
#' Combines two or more fitted models into one predictive distribution with
#' leave-one-out predictive stacking. Each member contributes its
#' Pareto-smoothed importance-sampling LOO predictive density; the weights
#' maximize the stacked leave-one-out predictive density and so favour members
#' that predict held-out observations best. Works across every family
#' (single-season / community / dynamic / integrated occupancy, JSDM, and the
#' cover hurdle), each scored by the marginal likelihood in
#' `.tobs_pointwise_loglik()`.
#'
#' The members must be fit to the same observation set. They may otherwise
#' differ freely -- different covariates, priors, or inference `method`s. (An
#' ensemble of pure seed-variants, e.g. from `tobs(..., control = list(n.seeds
#' = K))`, is a valid but degenerate case: the members are statistically
#' identical, so the weights come out roughly uniform.)
#'
#' @param ... two or more `tobs_fit` objects, optionally named, or a single
#'   list of them. Names label the members in the weight table.
#' @param method weighting scheme passed to [loo::loo_model_weights()]:
#'   `"stacking"` (default) or `"pseudobma"`.
#'
#' @return An object of class `tobs_stack`: a list with `weights` (named
#'   numeric, summing to 1), `fits` (the members), `loo` (per-member
#'   [loo::loo()] objects), `comparison` (a data frame of `elpd_loo`, `weight`,
#'   and worst Pareto-k per member), and `method`. [predict()] / [fitted()] on
#'   the result return the weight-combined predictive.
#'
#' @examples
#' \dontrun{
#' f1 <- tobs(~ elev,          data = sites, family = occu(),
#'            detection = ~ effort, y = y, method = "nuts")
#' f2 <- tobs(~ elev + forest, data = sites, family = occu(),
#'            detection = ~ effort, y = y, method = "nuts")
#' ens <- tobs_stack(simple = f1, full = f2)
#' ens$weights
#' predict(ens)            # weight-combined in-sample psi / p / z
#' }
#' @export
tobs_stack <- function(..., method = c("stacking", "pseudobma")) {
  method <- match.arg(method)

  fits <- list(...)
  # Accept a single list argument: tobs_stack(list(f1, f2)).
  if (length(fits) == 1L && is.list(fits[[1]]) &&
      !inherits(fits[[1]], "tobs_fit")) {
    fits <- fits[[1]]
  }
  if (length(fits) < 2L) {
    stop("tobs_stack() needs at least two fitted models to combine.",
         call. = FALSE)
  }
  is_fit <- vapply(fits, inherits, logical(1), what = "tobs_fit")
  if (!all(is_fit)) {
    stop("tobs_stack() inputs must all be `tobs_fit` objects; element(s) ",
         paste(which(!is_fit), collapse = ", "), " are not.", call. = FALSE)
  }
  if (!requireNamespace("loo", quietly = TRUE)) {
    stop("Package 'loo' is required for LOO stacking. Install it with ",
         "install.packages(\"loo\").", call. = FALSE)
  }

  nms <- names(fits)
  if (is.null(nms) || any(!nzchar(nms))) {
    nms <- paste0("model", seq_along(fits))
  }
  names(fits) <- nms

  # Per-member pointwise log-lik (single-season occupancy). The weights are
  # only comparable if every member scores the same observation set.
  ll_list <- lapply(fits, .tobs_pointwise_loglik)
  n_obs <- vapply(ll_list, ncol, integer(1))
  if (length(unique(n_obs)) != 1L) {
    stop("tobs_stack() members must be fit to the same observations; they ",
         "score differing numbers of observations (",
         paste(n_obs, collapse = ", "), ").", call. = FALSE)
  }

  loo_list <- Map(function(ll, f) .tobs_loo_one(ll, f$chain_id), ll_list, fits)
  names(loo_list) <- nms

  w <- loo::loo_model_weights(loo_list, method = method)
  weights <- stats::setNames(as.numeric(w), nms)

  elpd_loo <- vapply(loo_list,
                     function(l) l$estimates["elpd_loo", "Estimate"], numeric(1))
  pareto_k_max <- vapply(loo_list, function(l) {
    k <- tryCatch(loo::pareto_k_values(l), error = function(e) NA_real_)
    if (length(k) == 0L || all(is.na(k))) NA_real_ else max(k, na.rm = TRUE)
  }, numeric(1))

  comparison <- data.frame(
    model        = nms,
    elpd_loo     = elpd_loo,
    weight       = weights,
    pareto_k_max = pareto_k_max,
    row.names    = NULL,
    stringsAsFactors = FALSE
  )
  comparison <- comparison[order(-comparison$weight), , drop = FALSE]

  structure(
    list(weights = weights, fits = fits, loo = loo_list,
         comparison = comparison, method = method),
    class = "tobs_stack"
  )
}

# One loo::loo() per member. Uses relative effective sample sizes from the
# chain layout when available (NUTS carries fit$chain_id); single-chain /
# Laplace draws fall back to a single pseudo-chain.
.tobs_loo_one <- function(ll, chain_id = NULL) {
  n_draws <- nrow(ll)
  cid <- if (!is.null(chain_id) && length(chain_id) == n_draws) {
    as.integer(chain_id)
  } else {
    rep(1L, n_draws)
  }
  r_eff <- tryCatch(loo::relative_eff(exp(ll), chain_id = cid),
                    error = function(e) NULL)
  loo::loo(ll, r_eff = r_eff)
}


# ---------------------------------------------------------------------------
# Methods on tobs_stack
# ---------------------------------------------------------------------------

#' @param x a `tobs_stack` object.
#' @param ... ignored.
#' @rdname tobs_stack
#' @export
print.tobs_stack <- function(x, ...) {
  cat(sprintf("<tobs_stack: %d members, %s weights>\n",
              length(x$fits), x$method))
  cmp <- x$comparison
  out <- data.frame(
    model        = cmp$model,
    elpd_loo     = sprintf("%.1f", cmp$elpd_loo),
    weight       = sprintf("%.3f", cmp$weight),
    pareto_k_max = ifelse(is.na(cmp$pareto_k_max), "NA",
                          sprintf("%.2f", cmp$pareto_k_max)),
    stringsAsFactors = FALSE
  )
  print(out, row.names = FALSE)
  if (any(cmp$pareto_k_max > 0.7, na.rm = TRUE)) {
    cat("  note: Pareto k > 0.7 for some member(s); their LOO weight may be ",
        "unreliable.\n", sep = "")
  }
  invisible(x)
}

#' @param object a `tobs_stack` object.
#' @rdname tobs_stack
#' @export
fitted.tobs_stack <- function(object, ...) {
  w <- object$weights
  parts <- lapply(object$fits, function(f) {
    fv <- tryCatch(stats::fitted(f), error = function(e) NULL)
    if (!is.list(fv) || is.null(names(fv))) {
      stop("fitted()/predict() on a tobs_stack combines occupancy-family ",
           "members whose fitted() returns named components (psi / p / z). ",
           "For cover or mixed-family ensembles, combine $fits using $weights ",
           "directly.", call. = FALSE)
    }
    fv
  })
  keys <- Reduce(intersect, lapply(parts, names))
  if (!length(keys)) {
    stop("tobs_stack members expose no common fitted() components to combine.",
         call. = FALSE)
  }
  out <- lapply(keys, function(k) {
    Reduce(`+`, Map(function(p, wk) wk * p[[k]], parts, w))
  })
  names(out) <- keys
  out
}

#' @param X.0 optional occupancy design matrix for out-of-sample prediction.
#'   When supplied, all members must share the same occupancy design.
#' @param quantiles credible-interval quantile levels.
#' @param n.draws size of the pooled stacked-predictive sample (design-matrix
#'   mode); each member contributes a share proportional to its weight.
#' @rdname tobs_stack
#' @export
predict.tobs_stack <- function(object, X.0 = NULL,
                               quantiles = c(0.025, 0.5, 0.975),
                               n.draws = 4000L, ...) {
  # In-sample: weight-combine the member fitted() values (psi / p / z).
  if (is.null(X.0)) return(fitted(object))

  fits <- object$fits
  w    <- object$weights

  p_occ <- vapply(fits, function(f) f$model$process_info[[1]]$p, integer(1))
  if (length(unique(p_occ)) != 1L) {
    stop("predict(<tobs_stack>, X.0 = ) needs every member to share the same ",
         "occupancy design; members differ in occupancy coefficient count (",
         paste(p_occ, collapse = ", "), "). Predict per member instead.",
         call. = FALSE)
  }
  p_occ <- p_occ[1L]
  if (ncol(X.0) != p_occ) {
    stop(sprintf("X.0 has %d columns but the occupancy design has %d.",
                 ncol(X.0), p_occ), call. = FALSE)
  }

  # Stacked predictive: pool each member's psi draws in proportion to its
  # weight, then summarize the mixture.
  target <- max(length(fits), as.integer(n.draws))
  pooled <- lapply(seq_along(fits), function(k) {
    pd  <- .tobs_psi_draws(fits[[k]]$draws, X.0, p_occ)   # [S_k x n_pred]
    n_k <- max(1L, as.integer(round(w[k] * target)))
    pd[sample.int(nrow(pd), n_k, replace = TRUE), , drop = FALSE]
  })
  M <- do.call(rbind, pooled)

  data.frame(
    mean  = colMeans(M),
    sd    = apply(M, 2, stats::sd),
    q2.5  = apply(M, 2, stats::quantile, quantiles[1]),
    q50   = apply(M, 2, stats::quantile, quantiles[2]),
    q97.5 = apply(M, 2, stats::quantile, quantiles[3])
  )
}
