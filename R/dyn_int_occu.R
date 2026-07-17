# dyn_int_occu.R - Multi-season integrated occupancy (spOccupancy tIntPGOcc;
# gcol33/tulpaObs#122). The product of the two shipped families: a dynamic
# (multi-season HMM) occupancy state whose per-season emission pools SEVERAL
# detection sources (integrated occupancy). Integrated models exist because a
# single source is too sparse to identify psi; pooling sources across seasons is
# how colonization / extinction estimates come out of opportunistic data.
#
#   z_i1              ~ Bernoulli(psi1_i)
#   z_it | z_i,t-1    : colonization gamma_i (0 -> 1), survival 1 - eps_i (1 -> 1)
#   y_isjt | z_it = 1 ~ Bernoulli(p_s)             detection, source s, visit j, season t
#
# with a shared occupancy process across sources and a per-source detection
# probability. The latent occupancy sequence integrates out by the 2-state HMM
# forward recursion (the same vectorised forward as the dynamic-occupancy
# exact-marginal refine, R/dyn_occu_marginal.R), whose per-season emission is the
# PRODUCT over sources of each source's per-visit detection likelihood. The exact
# marginal is maximised (optim BFGS) with an observed-information vcov -- pure R,
# no new C++.
#
# v1 scope: every source covers all sites and the same T-season grid; constant
# (non-season-varying) colonization / extinction; a shared detection covariate
# design with per-source coefficients. Partial site / season overlap across
# sources (the general tIntPGOcc), season-varying transitions (the #124 recipe),
# an areal psi1 field (stIntPGOcc), and NUTS are documented follow-ups.
#
#   .tobs_build_dyn_int_occu()   data binder -> model_type = "dyn_int_occu"
#   .tobs_fit_dyn_int_occu()     optim over the HMM-forward multi-source marginal
#   .dispatch_dyn_int_occu()     tobs() entry (bind + fit + assemble)

# The exact HMM-forward multi-source marginal is inlined in the fitter and the
# pointwise-log-likelihood (both need the per-site, per-source detection vectors
# `.dio_unpack` returns); the forward recursion is the two-state colext forward
# with the per-season emission pooled over sources.

# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# `y` is a length-S list of [n_sites x max_visits_s x T] detection arrays (0/1/NA;
# NA = visit not conducted). `state_formula` models logit psi1; `col_formula` /
# `ext_formula` the site-level colonization / extinction; `det_formula` the shared
# per-source detection design (each source carries its own coefficients).
.tobs_build_dyn_int_occu <- function(state_formula, col_formula, ext_formula,
                                     det_formula, data, y, sources = NULL) {
  if (!is.list(y) || length(y) < 2L)
    stop("dyn_int_occu() y must be a list of >= 2 detection-source arrays ",
         "[sites x visits x seasons].", call. = FALSE)
  S <- length(y)
  if (is.null(sources)) sources <- names(y) %||% paste0("src", seq_len(S))
  dims <- lapply(y, dim)
  if (any(vapply(dims, length, 0L) != 3L))
    stop("each dyn_int_occu() source must be a 3D [sites x visits x seasons] ",
         "array.", call. = FALSE)
  n_sites <- dims[[1L]][1L]; T_s <- dims[[1L]][3L]
  if (any(vapply(dims, function(d) d[1L] != n_sites || d[3L] != T_s, TRUE)))
    stop("all dyn_int_occu() sources must share the site count and season grid ",
         "(v1: full overlap).", call. = FALSE)
  if (T_s < 2L) stop("dyn_int_occu() needs >= 2 seasons.", call. = FALSE)
  .tobs_check_site_count(n_sites, nrow(data), "sites")

  # Per-source [n_sites x T] detection sufficient statistics.
  nvalid <- ndet <- vector("list", S)
  for (s in seq_len(S)) {
    ys <- y[[s]]
    nv <- nd <- matrix(0L, n_sites, T_s)
    for (t in seq_len(T_s)) {
      yt <- ys[, , t]; v <- !is.na(yt)
      if (any(yt[v] != 0L & yt[v] != 1L))
        stop(sprintf("dyn_int_occu() source '%s' y must be 0/1/NA.", sources[s]),
             call. = FALSE)
      nv[, t] <- rowSums(v); yt0 <- yt; yt0[!v] <- 0L; nd[, t] <- rowSums(yt0)
    }
    nvalid[[s]] <- nv; ndet[[s]] <- nd
  }

  bind    <- .tobs_bind_formulas(list(psi = state_formula, gamma = col_formula,
                                      eps = ext_formula, p = det_formula), data)
  X_psi <- stats::model.matrix(bind$fe$psi,   data)
  X_gam <- stats::model.matrix(bind$fe$gamma, data)
  X_eps <- stats::model.matrix(bind$fe$eps,   data)
  X_det <- stats::model.matrix(bind$fe$p,     data)

  pi_list <- list(
    list(name = "psi1",  p = ncol(X_psi), coef_names = colnames(X_psi), link = "logit"),
    list(name = "gamma", p = ncol(X_gam), coef_names = colnames(X_gam), link = "logit"),
    list(name = "eps",   p = ncol(X_eps), coef_names = colnames(X_eps), link = "logit"))
  for (s in seq_len(S))
    pi_list[[length(pi_list) + 1L]] <- list(
      name = paste0("p_", sources[s]), p = ncol(X_det),
      coef_names = colnames(X_det), link = "logit")

  structure(list(
    model_type  = "dyn_int_occu",
    y           = y,
    sources     = sources,
    S           = S,
    nvalid      = nvalid,
    ndet        = ndet,
    X_psi       = X_psi, X_gam = X_gam, X_eps = X_eps, X_det = X_det,
    formulas    = list(psi = bind$fe$psi, gamma = bind$fe$gamma,
                       eps = bind$fe$eps, p = bind$fe$p),
    structured_terms = bind$terms,
    data        = data,
    n_sites     = n_sites,
    n_seasons   = T_s,
    process_info = pi_list
  ), class = "tobs_model")
}

.dio_unpack <- function(theta, model) {
  p_psi <- ncol(model$X_psi); p_gam <- ncol(model$X_gam)
  p_eps <- ncol(model$X_eps); p_det <- ncol(model$X_det); S <- model$S
  o <- 0L
  psi1  <- stats::plogis(as.vector(model$X_psi %*% theta[o + seq_len(p_psi)])); o <- o + p_psi
  gamma <- stats::plogis(as.vector(model$X_gam %*% theta[o + seq_len(p_gam)])); o <- o + p_gam
  eps   <- stats::plogis(as.vector(model$X_eps %*% theta[o + seq_len(p_eps)])); o <- o + p_eps
  # Per-source, per-site detection probability (covariate-aware).
  p_site <- vector("list", S)
  for (s in seq_len(S)) {
    p_site[[s]] <- stats::plogis(as.vector(model$X_det %*% theta[o + seq_len(p_det)]))
    o <- o + p_det
  }
  list(psi1 = psi1, gamma = gamma, eps = eps, p_site = p_site)
}

# ---------------------------------------------------------------------------
# Fitter
# ---------------------------------------------------------------------------

.tobs_fit_dyn_int_occu <- function(model, verbose = TRUE, ...) {
  S <- model$S; T_s <- model$n_seasons; n_sites <- model$n_sites
  p_psi <- ncol(model$X_psi); p_gam <- ncol(model$X_gam)
  p_eps <- ncol(model$X_eps); p_det <- ncol(model$X_det)
  n_theta <- p_psi + p_gam + p_eps + S * p_det
  nvalid <- model$nvalid; ndet <- model$ndet

  # Per-site HMM-forward marginal with per-source (per-site) detection.
  site_ll <- function(up) {
    lg  <- log(pmax(up$psi1, 1e-12));  l1g <- log(pmax(1 - up$psi1, 1e-12))
    lgam <- log(pmax(up$gamma, 1e-12)); l1gam <- log(pmax(1 - up$gamma, 1e-12))
    leps <- log(pmax(up$eps, 1e-12));   l1eps <- log(pmax(1 - up$eps, 1e-12))
    emit1 <- matrix(0, n_sites, T_s); det_any <- matrix(0L, n_sites, T_s)
    for (s in seq_len(S)) {
      lp  <- log(pmax(up$p_site[[s]], 1e-12))
      l1p <- log(pmax(1 - up$p_site[[s]], 1e-12))
      emit1 <- emit1 + ndet[[s]] * lp + (nvalid[[s]] - ndet[[s]]) * l1p
      det_any <- det_any + ndet[[s]]
    }
    emit0 <- ifelse(det_any > 0L, -Inf, 0)
    lse2 <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }
    la1 <- lg + emit1[, 1L]; la0 <- l1g + emit0[, 1L]
    for (t in 2:T_s) {
      n1 <- lse2(la1 + l1eps, la0 + lgam) + emit1[, t]
      n0 <- lse2(la1 + leps,  la0 + l1gam) + emit0[, t]
      la1 <- n1; la0 <- n0
    }
    lse2(la1, la0)
  }

  nll <- function(theta) {
    up  <- .dio_unpack(theta, model)
    ll  <- site_ll(up)
    val <- -sum(ll[is.finite(ll)])
    if (is.finite(val)) val else 1e10
  }

  # Init: psi1 / detection from the pooled first-season detection; gamma / eps
  # from a modest transition guess.
  any_det1 <- Reduce(`+`, lapply(ndet, function(m) m[, 1L])) > 0
  init <- numeric(n_theta)
  init[1L] <- stats::qlogis(min(max(mean(any_det1), 0.1), 0.9))
  init[p_psi + 1L] <- stats::qlogis(0.3)                 # gamma
  init[p_psi + p_gam + 1L] <- stats::qlogis(0.3)         # eps
  o <- p_psi + p_gam + p_eps
  for (s in seq_len(S)) {
    ph <- sum(ndet[[s]]) / max(sum(nvalid[[s]]), 1)
    init[o + 1L] <- stats::qlogis(min(max(ph, 0.05), 0.95)); o <- o + p_det
  }

  opt <- stats::optim(init, nll, method = "BFGS", hessian = TRUE,
                      control = list(maxit = 800L))
  converged <- opt$convergence == 0L

  par_names <- unlist(lapply(model$process_info, function(pp)
    paste0(pp$name, "_", pp$coef_names)))
  means <- opt$par; names(means) <- par_names
  V <- tryCatch(solve(opt$hessian),
                error = function(e) diag(NA_real_, length(means)))
  V <- (V + t(V)) / 2
  dimnames(V) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  structure(c(list(
    draws = draws, means = means, sds = sds, vcov = V,
    n_samples = n_draws, n_params = length(means),
    log_prob = rep(-opt$value, n_draws), log_lik = -opt$value, N = n_sites),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names = par_names, param_names = par_names,
    n_fixed = length(means), fixed_names = par_names,
    process_info = model$process_info,
    model = model, spatial = NULL, method = "laplace",
    convergence = list(converged = converged, n_iter = opt$counts[[1L]])
  )), class = c("tobs_fit", "tulpa_fit"))
}

# ---------------------------------------------------------------------------
# tobs() dispatcher
# ---------------------------------------------------------------------------

.dispatch_dyn_int_occu <- function(formula, data, family, detection, y, visits,
                                   engine, priors, control,
                                   approx = "gaussian_laplace",
                                   correction = "none", colonization = NULL,
                                   extinction = NULL, sources = NULL, ...) {
  if (is.null(detection))
    stop("dyn_int_occu() requires a `detection` formula (shared per-source ",
         "site-level detection).", call. = FALSE)
  if (is.null(y))
    stop("dyn_int_occu() requires `y` (a list of detection-source arrays ",
         "[sites x visits x seasons]).", call. = FALSE)
  if (is.null(colonization) || is.null(extinction))
    stop("dyn_int_occu() requires `colonization = ~ ...` and `extinction = ~ ...`",
         " transition formulas (as dyn_occu does).", call. = FALSE)
  if (!is.null(visits))
    stop("dyn_int_occu() detection is site-level; visit-level detection ",
         "covariates (`visits`) are not yet supported.", call. = FALSE)
  if (!identical(.map_engine(engine, family = "dyn_int_occu"), "laplace"))
    stop("dyn_int_occu() supports method = \"laplace\" only.", call. = FALSE)
  model <- .tobs_build_dyn_int_occu(
    state_formula = formula, col_formula = colonization,
    ext_formula = extinction, det_formula = detection, data = data, y = y,
    sources = sources)
  .tobs_fit_dyn_int_occu(model, verbose = isTRUE(control$verbose))
}

# ---------------------------------------------------------------------------
# fitted / predict / residuals / pointwise log-likelihood / simulate
# ---------------------------------------------------------------------------

.tobs_fitted_dyn_int_occu <- function(object) {
  up <- .dio_unpack(object$means, object$model)
  p_mat <- do.call(cbind, up$p_site); colnames(p_mat) <- paste0("p_", object$model$sources)
  list(psi1 = up$psi1, gamma = up$gamma, eps = up$eps, p = p_mat)
}

.tobs_predict_dyn_int_occu <- function(object, newdata = NULL,
                                       type = c("state", "colonization",
                                                "extinction")) {
  type  <- match.arg(type)
  model <- object$model
  key   <- switch(type, state = "psi", colonization = "gamma", extinction = "eps")
  fk    <- switch(type, state = "psi1", colonization = "gamma", extinction = "eps")
  X <- if (is.null(newdata)) switch(fk, psi1 = model$X_psi, gamma = model$X_gam,
                                    eps = model$X_eps)
       else stats::model.matrix(model$formulas[[key]], newdata)
  # offset into means for this arm
  off <- switch(fk, psi1 = 0L, gamma = ncol(model$X_psi),
                eps = ncol(model$X_psi) + ncol(model$X_gam))
  pn <- switch(fk, psi1 = ncol(model$X_psi), gamma = ncol(model$X_gam),
               eps = ncol(model$X_eps))
  stats::plogis(as.vector(X %*% object$means[off + seq_len(pn)]))
}

.tobs_residuals_dyn_int_occu <- function(object, type) {
  model <- object$model
  fv    <- .tobs_fitted_dyn_int_occu(object)
  # Observed vs expected any-detection over all sources / seasons per site.
  det_any <- Reduce(`+`, lapply(model$ndet, rowSums)) > 0
  # Marginal P(detected at least once) ~ psi1 * (1 - prod over seasons/sources of
  # (1 - p)); a coarse site-level check.
  pmax_det <- 1 - Reduce(`*`, lapply(seq_len(model$S), function(s)
    (1 - fv$p[, s])^rowSums(model$nvalid[[s]])))
  pdet <- fv$psi1 * pmax_det
  eps <- 1e-10; pc <- pmin(pmax(pdet, eps), 1 - eps)
  obs <- as.numeric(det_any)
  res <- switch(type,
    response = obs - pdet,
    pearson  = (obs - pdet) / sqrt(pc * (1 - pc) + eps),
    deviance = sign(obs - pdet) * sqrt(2 * abs(
      ifelse(obs > 0, obs * log(obs / pc), 0) +
      ifelse(obs < 1, (1 - obs) * log((1 - obs) / (1 - pc)), 0))))
  list(occ = res, det = NULL)
}

.tobs_ploglik_dyn_int_occu <- function(object, n.draws = 1000L, n.threads = 1L) {
  model <- object$model
  draws <- object$draws
  if (!is.null(n.draws) && as.integer(n.draws) < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  t(vapply(seq_len(nrow(draws)), function(d) {
    up <- .dio_unpack(draws[d, ], model)
    lg  <- log(pmax(up$psi1, 1e-12));  l1g <- log(pmax(1 - up$psi1, 1e-12))
    lgam <- log(pmax(up$gamma, 1e-12)); l1gam <- log(pmax(1 - up$gamma, 1e-12))
    leps <- log(pmax(up$eps, 1e-12));   l1eps <- log(pmax(1 - up$eps, 1e-12))
    T_s <- model$n_seasons; n_sites <- model$n_sites
    emit1 <- matrix(0, n_sites, T_s); det_any <- matrix(0L, n_sites, T_s)
    for (s in seq_len(model$S)) {
      lp <- log(pmax(up$p_site[[s]], 1e-12)); l1p <- log(pmax(1 - up$p_site[[s]], 1e-12))
      emit1 <- emit1 + model$ndet[[s]] * lp + (model$nvalid[[s]] - model$ndet[[s]]) * l1p
      det_any <- det_any + model$ndet[[s]]
    }
    emit0 <- ifelse(det_any > 0L, -Inf, 0)
    lse2 <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }
    la1 <- lg + emit1[, 1L]; la0 <- l1g + emit0[, 1L]
    for (t in 2:T_s) {
      n1 <- lse2(la1 + l1eps, la0 + lgam) + emit1[, t]
      n0 <- lse2(la1 + leps,  la0 + l1gam) + emit0[, t]
      la1 <- n1; la0 <- n0
    }
    lse2(la1, la0)
  }, numeric(model$n_sites)))
}

.tobs_simulate_dyn_int_occu <- function(object, nsim = 1) {
  model <- object$model; S <- model$S; T_s <- model$n_seasons
  valid <- lapply(model$y, function(a) !is.na(a))
  draw_one <- function() {
    idx <- sample.int(nrow(object$draws), 1L)
    up  <- .dio_unpack(object$draws[idx, ], model)
    z <- matrix(0L, model$n_sites, T_s)
    z[, 1L] <- stats::rbinom(model$n_sites, 1L, up$psi1)
    for (t in 2:T_s) {
      surv <- stats::rbinom(model$n_sites, 1L, 1 - up$eps)
      col  <- stats::rbinom(model$n_sites, 1L, up$gamma)
      z[, t] <- ifelse(z[, t - 1L] == 1L, surv, col)
    }
    out <- lapply(seq_len(S), function(s) array(NA_integer_, dim(model$y[[s]])))
    for (s in seq_len(S)) {
      for (t in seq_len(T_s)) {
        for (j in seq_len(dim(model$y[[s]])[2L])) {
          vv <- valid[[s]][, j, t]
          det <- (z[, t] == 1L) & vv
          out[[s]][, j, t][vv] <- 0L
          out[[s]][, j, t][det] <- stats::rbinom(sum(det), 1L, up$p_site[[s]][det])
        }
      }
    }
    names(out) <- model$sources
    out
  }
  if (nsim == 1L) return(draw_one())
  lapply(seq_len(nsim), function(s) draw_one())
}

# ---------------------------------------------------------------------------
# Simulator for recovery tests
# ---------------------------------------------------------------------------

#' Simulate a multi-season integrated occupancy data set
#'
#' Draws from the [dyn_int_occu()] model: a dynamic (multi-season) occupancy
#' process (`psi1`, colonization `gamma`, extinction `eps`) observed by `S`
#' detection sources with per-source detection probability.
#'
#' @param N Number of sites (default 200).
#' @param T_seasons Number of seasons (default 4).
#' @param S Number of detection sources (default 2).
#' @param J Visits per source (scalar or length-`S`; default 3).
#' @param psi1,gamma,eps Season-1 occupancy, colonization, extinction (defaults
#'   0.5 / 0.3 / 0.2).
#' @param p Per-source detection probabilities (length-`S`; default 0.4 / 0.6).
#' @param seed Optional random seed.
#' @return A list with `y` (a length-`S` list of `[N x J x T]` arrays), `data`,
#'   `sources`, and `truth`.
#' @export
simulate_dyn_int_occu <- function(N = 200, T_seasons = 4, S = 2, J = 3,
                                  psi1 = 0.5, gamma = 0.3, eps = 0.2,
                                  p = c(0.4, 0.6), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(J) != S) J <- rep(J[1L], S)
  if (length(p) != S) p <- rep(p[1L], S)
  data <- data.frame(row.names = seq_len(N))
  z <- matrix(0L, N, T_seasons)
  z[, 1L] <- stats::rbinom(N, 1L, psi1)
  for (t in 2:T_seasons) {
    surv <- stats::rbinom(N, 1L, 1 - eps); col <- stats::rbinom(N, 1L, gamma)
    z[, t] <- ifelse(z[, t - 1L] == 1L, surv, col)
  }
  y <- lapply(seq_len(S), function(s) {
    arr <- array(0L, c(N, J[s], T_seasons))
    for (t in seq_len(T_seasons)) for (j in seq_len(J[s]))
      arr[, j, t] <- ifelse(z[, t] == 1L, stats::rbinom(N, 1L, p[s]), 0L)
    arr
  })
  names(y) <- paste0("src", seq_len(S))
  list(y = y, data = data, sources = names(y),
       truth = list(psi1 = psi1, gamma = gamma, eps = eps, p = p, z = z))
}
