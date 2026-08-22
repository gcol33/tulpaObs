# dyn_int_occu.R
# - Multi-season integrated occupancy (spOccupancy tIntPGOcc). The product of
# the two shipped families: a dynamic (multi-season HMM) occupancy state whose
# per-season emission pools SEVERAL detection sources (integrated occupancy).
# Integrated models exist because a single source is too sparse to identify psi;
# pooling sources across seasons is how colonization / extinction estimates come
# out of opportunistic data.
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
# NA = visit not conducted). Sources share the site count and season grid, but a
# source that does not observe a (site, season) marks it NA -- so PARTIAL season
# overlap (a staggered survey where sources cover different seasons) is expressed
# by NA-padding each source to the common grid: a source absent at season t
# contributes nothing to that season's emission (nvalid = 0), and a (site, season)
# unobserved by every source is marginalised (e0 = e1 = 1) by the forward.
# `state_formula` models logit psi1; `col_formula` / `ext_formula` the site-level
# colonization / extinction; `det_formula` the shared per-source detection design
# (each source carries its own coefficients).
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
         "([n_sites x max_visits_s x T]); a source that does not cover a ",
         "(site, season) marks it NA (partial season overlap is NA-padded to the ",
         "common grid).", call. = FALSE)
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

# Vectorised scaled forward-backward for the multi-source colext HMM. Returns the
# per-site log marginal AND the exact (Fisher-identity) per-site eta gradients on
# every arm -- the smoothed season-1 occupancy w1[,1] gives the psi1 (and areal
# field) gradient w1 - psi1, the pairwise transition joints give the colonization
# / extinction score, and the per-season smoothed occupancy weights the per-source
# detection binomial score. Emission stays on the probability scale (e1 <= 1,
# e0 in {0, 1}, so at least one state is O(1) per season -- the per-season
# normaliser cs keeps it stable without log-space). psi1 / gamma / eps / detection
# are all constant across a site's seasons here, so the transition is a single
# per-site 2x2.
.dio_fb <- function(up, model, offset_psi1 = NULL) {
  S <- model$S; T_s <- model$n_seasons; n <- model$n_sites
  nvalid <- model$nvalid; ndet <- model$ndet
  psi1 <- up$psi1
  if (!is.null(offset_psi1))
    psi1 <- stats::plogis(stats::qlogis(pmin(pmax(psi1, 1e-12), 1 - 1e-12)) + offset_psi1)
  gamma <- up$gamma; eps <- up$eps

  emit1_log <- matrix(0, n, T_s); det_any <- matrix(0L, n, T_s)
  for (s in seq_len(S)) {
    lp  <- log(pmax(up$p_site[[s]], 1e-12))
    l1p <- log(pmax(1 - up$p_site[[s]], 1e-12))
    emit1_log <- emit1_log + ndet[[s]] * lp + (nvalid[[s]] - ndet[[s]]) * l1p
    det_any <- det_any + ndet[[s]]
  }
  e1 <- exp(emit1_log)                        # [n x T], <= 1
  e0 <- ifelse(det_any > 0L, 0, 1)            # [n x T]

  # Scaled forward.
  a0 <- matrix(0, n, T_s); a1 <- matrix(0, n, T_s); cs <- matrix(0, n, T_s)
  u0 <- (1 - psi1) * e0[, 1]; u1 <- psi1 * e1[, 1]; c1 <- u0 + u1
  a0[, 1] <- u0 / c1; a1[, 1] <- u1 / c1; cs[, 1] <- c1
  if (T_s > 1L) for (t in 2:T_s) {
    pr0 <- a0[, t - 1] * (1 - gamma) + a1[, t - 1] * eps
    pr1 <- a0[, t - 1] * gamma       + a1[, t - 1] * (1 - eps)
    v0 <- e0[, t] * pr0; v1 <- e1[, t] * pr1; ct <- v0 + v1
    a0[, t] <- v0 / ct; a1[, t] <- v1 / ct; cs[, t] <- ct
  }
  loglik <- rowSums(log(pmax(cs, 1e-300)))

  # Scaled backward + smoothed marginals / pairwise joints.
  b0 <- matrix(0, n, T_s); b1 <- matrix(0, n, T_s)
  b0[, T_s] <- 1; b1[, T_s] <- 1
  w0 <- matrix(0, n, T_s); w1 <- matrix(0, n, T_s)
  w0[, T_s] <- a0[, T_s]; w1[, T_s] <- a1[, T_s]
  col_y <- numeric(n); ext_y <- numeric(n)
  if (T_s > 1L) for (t in (T_s - 1):1) {
    bb0 <- e0[, t + 1] * b0[, t + 1]; bb1 <- e1[, t + 1] * b1[, t + 1]
    inv_c <- 1 / cs[, t + 1]
    xi01 <- a0[, t] * gamma * bb1 * inv_c        # colonization event
    xi10 <- a1[, t] * eps   * bb0 * inv_c        # extinction event
    col_y <- col_y + xi01; ext_y <- ext_y + xi10
    b0[, t] <- ((1 - gamma) * bb0 + gamma * bb1) * inv_c
    b1[, t] <- (eps * bb0 + (1 - eps) * bb1) * inv_c
    w0[, t] <- a0[, t] * b0[, t]; w1[, t] <- a1[, t] * b1[, t]
  }

  # Per-site eta gradients (Fisher identity). Transition origin sums run over the
  # T-1 intervals (seasons 1..T-1); detection over occupied season weights.
  g_eta_psi1 <- w1[, 1] - psi1
  keep <- if (T_s > 1L) seq_len(T_s - 1L) else integer(0)
  g_eta_gam <- col_y - gamma * rowSums(w0[, keep, drop = FALSE])
  g_eta_eps <- ext_y - eps   * rowSums(w1[, keep, drop = FALSE])
  g_eta_p <- lapply(seq_len(S), function(s)
    rowSums(w1 * (ndet[[s]] - nvalid[[s]] * matrix(up$p_site[[s]], n, T_s))))

  list(loglik = loglik, g_eta_psi1 = g_eta_psi1, g_eta_gam = g_eta_gam,
       g_eta_eps = g_eta_eps, g_eta_p = g_eta_p, w1 = w1)
}

.tobs_fit_dyn_int_occu <- function(model, verbose = TRUE, ...) {
  S <- model$S; T_s <- model$n_seasons; n_sites <- model$n_sites
  p_psi <- ncol(model$X_psi); p_gam <- ncol(model$X_gam)
  p_eps <- ncol(model$X_eps); p_det <- ncol(model$X_det)
  n_theta <- p_psi + p_gam + p_eps + S * p_det
  nvalid <- model$nvalid; ndet <- model$ndet
  X_psi <- model$X_psi; X_gam <- model$X_gam; X_eps <- model$X_eps
  X_det <- model$X_det

  nll <- function(theta) {
    up  <- .dio_unpack(theta, model)
    ll  <- .dio_fb(up, model)$loglik
    val <- -sum(ll[is.finite(ll)])
    if (is.finite(val)) val else 1e10
  }
  # Analytic gradient over the exact forward-backward smoothing (no finite diff).
  ngr <- function(theta) {
    up <- .dio_unpack(theta, model)
    fb <- .dio_fb(up, model)
    g <- numeric(n_theta); o <- 0L
    g[o + seq_len(p_psi)] <- crossprod(X_psi, fb$g_eta_psi1); o <- o + p_psi
    g[o + seq_len(p_gam)] <- crossprod(X_gam, fb$g_eta_gam); o <- o + p_gam
    g[o + seq_len(p_eps)] <- crossprod(X_eps, fb$g_eta_eps); o <- o + p_eps
    for (s in seq_len(S)) {
      g[o + seq_len(p_det)] <- crossprod(X_det, fb$g_eta_p[[s]]); o <- o + p_det
    }
    -g
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

  opt <- stats::optim(init, nll, ngr, method = "BFGS",
                      control = list(maxit = 800L))
  converged <- opt$convergence == 0L
  # Observed-information vcov from the FD-Jacobian of the analytic gradient
  # (O(p) marginal evals; no numeric Hessian over the forward-backward). `ngr` is
  # the negative-log-likelihood gradient, so its Jacobian at the minimum IS the
  # observed information (positive definite); solve() gives the vcov directly.
  opt$hessian <- .tobs_fd_jacobian(ngr, opt$par)

  par_names <- unlist(lapply(model$process_info, function(pp)
    paste0(pp$name, "_", pp$coef_names)))
  means <- opt$par; names(means) <- par_names
  V <- tryCatch(solve(opt$hessian),
                error = function(e) diag(NA_real_, length(means)))
  V <- (V + t(V)) / 2
  dimnames(V) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V)
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
# Areal field on the first-season occupancy arm (stIntPGOcc)
# ---------------------------------------------------------------------------

# An ICAR field on the initial-occupancy (psi1) arm of the multi-season
# integrated model. psi1 sets ONLY the initial mixing weight of each site's HMM,
# so the per-site marginal is linear in psi1 and the exact per-site field gradient
# is the Fisher-identity score w1[,1] - psi1 (the smoothed season-1 occupancy),
# which .dio_fb already returns. The shared areal-BFGS driver (R/areal_bfgs.R)
# runs BFGS over (all fixed coefficients, field) + the CAR prior and forms the
# Laplace marginal from an FD-Hessian at the mode -- the same recipe as fp_occu /
# dyn_abun, differing only in the family's eval. One field unit per site; the
# colonization / extinction / per-source detection arms carry fixed effects only.
.tobs_fit_dyn_int_occu_spatial <- function(model, spatial, max_iter = 200L,
                                           tol = 1e-8, verbose = TRUE,
                                           integration = "grid") {
  if (!identical(spatial$type, "icar"))
    stop("dyn_int_occu() + a spatial field supports icar() only in v1 ",
         "(bym2 / car_proper are follow-ups).", call. = FALSE)
  S <- model$S; n_sites <- model$n_sites
  X_psi <- model$X_psi; X_gam <- model$X_gam; X_eps <- model$X_eps
  X_det <- model$X_det
  p_psi <- ncol(X_psi); p_gam <- ncol(X_gam); p_eps <- ncol(X_eps); p_det <- ncol(X_det)
  off <- cumsum(c(0L, p_psi, p_gam, p_eps, rep(p_det, S)))
  i_psi <- off[1] + seq_len(p_psi); i_gam <- off[2] + seq_len(p_gam)
  i_eps <- off[3] + seq_len(p_eps)
  i_p   <- lapply(seq_len(S), function(s) off[3 + s] + seq_len(p_det))
  n_fixed <- off[length(off)]
  fb <- .tobs_areal_field_blocks(spatial, n_sites, "dyn_int_occu", model$data)
  field <- if (length(fb$blocks) == 1L) fb$blocks[[1L]] else fb$blocks

  unpack_fix <- function(theta_fix) list(
    psi1  = stats::plogis(as.numeric(X_psi %*% theta_fix[i_psi])),
    gamma = stats::plogis(as.numeric(X_gam %*% theta_fix[i_gam])),
    eps   = stats::plogis(as.numeric(X_eps %*% theta_fix[i_eps])),
    p_site = lapply(seq_len(S), function(s)
      stats::plogis(as.numeric(X_det %*% theta_fix[i_p[[s]]]))))

  eval <- function(theta_fix, offset) {
    fb <- .dio_fb(unpack_fix(theta_fix), model, offset_psi1 = offset)
    g <- numeric(n_fixed)
    g[i_psi] <- crossprod(X_psi, fb$g_eta_psi1)
    g[i_gam] <- crossprod(X_gam, fb$g_eta_gam)
    g[i_eps] <- crossprod(X_eps, fb$g_eta_eps)
    for (s in seq_len(S)) g[i_p[[s]]] <- crossprod(X_det, fb$g_eta_p[[s]])
    ll <- fb$loglik
    list(log_lik = sum(ll[is.finite(ll)]), grad_fixed = g,
         grad_eta = fb$g_eta_psi1)
  }

  warm <- tryCatch(.tobs_fit_dyn_int_occu(model, verbose = FALSE),
                   error = function(e) NULL)
  theta0_fix <- if (!is.null(warm)) as.numeric(warm$means) else numeric(n_fixed)

  res <- .tobs_areal_bfgs_fit(eval, n_fixed, field, theta0_fix,
                              max_iter = max_iter, tol = tol,
                              label = "dyn-int-occu-spatial", integration = integration)

  nm <- unlist(lapply(model$process_info, function(pp)
    paste0(pp$name, "_", pp$coef_names)))
  means <- res$beta_mean; names(means) <- nm
  V <- res$vcov; dimnames(V) <- list(nm, nm)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- nm
  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V); colnames(draws) <- nm

  # Intercept field on the legacy scalar slots; any weighted (SVC) blocks become
  # the trend field(s) -- svcTIntPGOcc.
  fmeans <- res$field_means %||% list(res$field_mean)
  trend_labels <- fb$labels[-1L]
  trend_fields <- if (length(fmeans) > 1L) {
    tf <- fmeans[-1L]; names(tf) <- trend_labels; tf
  } else NULL

  structure(c(list(
    draws = draws, means = means, sds = sds, vcov = V,
    n_samples = n_draws, n_params = length(means),
    log_prob = rep(res$log_lik, n_draws), log_lik = res$log_lik, N = n_sites),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names = nm, param_names = nm, n_fixed = length(means), fixed_names = nm,
    process_info = model$process_info, model = model,
    spatial = spatial, spatial_field = fmeans[[1L]], spatial_hyper = res$hyper,
    trend_field = if (!is.null(trend_fields)) trend_fields[[1L]] else NULL,
    trend_fields = trend_fields,
    spatial_integration = res$integration, spatial_pareto_k = res$pareto_k,
    method = "nested_laplace",
    convergence = list(converged = TRUE, n_iter = NA_integer_)
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
  model <- .tobs_build_dyn_int_occu(
    state_formula = formula, col_formula = colonization,
    ext_formula = extinction, det_formula = detection, data = data, y = y,
    sources = sources)

  # A shared areal field on the first-season occupancy formula routes to the
  # stIntPGOcc fitter under nested_laplace; otherwise the non-spatial Laplace
  # fit. Only a psi1-arm icar() field is supported.
  structs <- .tobs_structures_from_model(model)
  if (!is.null(structs$temporal) || !is.null(structs$re) ||
      !is.null(structs$svc) || !is.null(structs$latent))
    stop("dyn_int_occu(): temporal / re / svc / latent terms are not wired; a ",
         "shared areal field icar() on the first-season occupancy formula is the ",
         "structured term supported.", call. = FALSE)
  if (!is.null(structs$spatial)) {
    if (!isTRUE(structs$spatial$shared[1L]))
      stop("dyn_int_occu() areal field sits on the first-season occupancy arm ",
           "only (a field on colonization / extinction / detection is not ",
           "supported).", call. = FALSE)
    if (!identical(engine, "nested_laplace"))
      stop("a shared areal field on the dyn_int_occu() occupancy formula needs ",
           "method = \"nested_laplace\" (drop the icar() term for the ",
           "non-spatial fit).", call. = FALSE)
    return(.tobs_fit_dyn_int_occu_spatial(
      model, spatial = structs$spatial,
      max_iter = control[["max.iter"]] %||% 200L,
      tol = control[["tol"]] %||% 1e-8,
      verbose = isTRUE(control$verbose),
      integration = control[["integration"]] %||% "grid"))
  }
  if (!identical(.map_engine(engine, family = "dyn_int_occu"), "laplace"))
    stop("dyn_int_occu() supports method = \"laplace\" (non-spatial) or ",
         "\"nested_laplace\" (with an icar() field on the occupancy formula).",
         call. = FALSE)
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
#' @param field Optional per-site shared areal field (length `N`) added to the
#'   first-season occupancy logit -- the shared field of the multi-season
#'   integrated spatial model (stIntPGOcc). Default `NULL` (no field).
#' @param trend Optional per-site varying-coefficient (SVC) areal field (length
#'   `N`); with a covariate `w ~ N(0,1)` stored in `data`, the first-season logit
#'   gains `w * trend` on top of `field` -- the svcTIntPGOcc surface. Default
#'   `NULL`. When set, `data` carries a `cell` node index (`1..N`) and the
#'   covariate `w` for the bar `spatial(~ 1 + w || cell, graph)`.
#' @param source_seasons Optional length-`S` list; `source_seasons[[s]]` is the
#'   integer vector of seasons source `s` observes (partial season overlap). The
#'   seasons a source does not cover are set to `NA` in its array -- the staggered
#'   survey where sources rarely share the full season grid. Default `NULL` (every
#'   source observes every season).
#' @param seed Optional random seed.
#' @return A list with `y` (a length-`S` list of `[N x J x T]` arrays), `data`,
#'   `sources`, and `truth`.
#' @export
simulate_dyn_int_occu <- function(N = 200, T_seasons = 4, S = 2, J = 3,
                                  psi1 = 0.5, gamma = 0.3, eps = 0.2,
                                  p = c(0.4, 0.6), field = NULL, trend = NULL,
                                  source_seasons = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(J) != S) J <- rep(J[1L], S)
  if (length(p) != S) p <- rep(p[1L], S)
  if (!is.null(field) && length(field) != N)
    stop("simulate_dyn_int_occu(): `field` must have length N.", call. = FALSE)
  if (!is.null(trend) && length(trend) != N)
    stop("simulate_dyn_int_occu(): `trend` must have length N.", call. = FALSE)
  if (!is.null(source_seasons)) {
    if (!is.list(source_seasons) || length(source_seasons) != S)
      stop("simulate_dyn_int_occu(): `source_seasons` must be a length-S list.",
           call. = FALSE)
    if (any(unlist(source_seasons) < 1L | unlist(source_seasons) > T_seasons))
      stop("simulate_dyn_int_occu(): `source_seasons` indices must be in 1..T.",
           call. = FALSE)
  }
  data <- data.frame(row.names = seq_len(N))
  z <- matrix(0L, N, T_seasons)
  # Season-1 occupancy carries the optional shared field + a varying-coefficient
  # field (weighted by a covariate w) on the logit scale. No field and no trend
  # keeps the exact constant-psi1 path (byte-identical to the pre-SVC simulator).
  w <- NULL
  if (is.null(field) && is.null(trend)) {
    psi1_i <- rep(psi1, N)
  } else {
    eta1 <- stats::qlogis(psi1)
    if (!is.null(field)) eta1 <- eta1 + field
    if (!is.null(trend)) {
      w <- stats::rnorm(N)
      data$cell <- seq_len(N); data$w <- w
      eta1 <- eta1 + w * trend
    }
    psi1_i <- stats::plogis(eta1)
  }
  z[, 1L] <- stats::rbinom(N, 1L, psi1_i)
  for (t in 2:T_seasons) {
    surv <- stats::rbinom(N, 1L, 1 - eps); col <- stats::rbinom(N, 1L, gamma)
    z[, t] <- ifelse(z[, t - 1L] == 1L, surv, col)
  }
  y <- lapply(seq_len(S), function(s) {
    arr <- array(0L, c(N, J[s], T_seasons))
    for (t in seq_len(T_seasons)) for (j in seq_len(J[s]))
      arr[, j, t] <- ifelse(z[, t] == 1L, stats::rbinom(N, 1L, p[s]), 0L)
    # Partial season overlap: source s observes only source_seasons[[s]]; the
    # seasons it does not cover are NA (absent), the per-source-season-map form
    # of a staggered survey where sources rarely share the full season grid.
    if (!is.null(source_seasons)) {
      miss <- setdiff(seq_len(T_seasons), source_seasons[[s]])
      if (length(miss)) arr[, , miss] <- NA_integer_
    }
    arr
  })
  names(y) <- paste0("src", seq_len(S))
  list(y = y, data = data, sources = names(y),
       truth = list(psi1 = psi1, gamma = gamma, eps = eps, p = p, z = z,
                    field = field, trend = trend, w = w,
                    source_seasons = source_seasons))
}
