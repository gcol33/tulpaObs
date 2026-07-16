# =============================================================================
# ms_distance.R - community / multispecies binned distance sampling
# (the spAbundance msDS analogue), with optional latent factors (lfMsDS) and a
# shared field (the spatial-factor case, sfMsDS). gcol33/tulpaObs#117. Poisson.
#
# Per-species binned distance sampling with Gaussian community hyperpriors on the
# per-species abundance and detection-scale coefficients:
#
#   N_{s,i}          ~ Poisson(lambda_{s,i})
#   y_{s,i,.} | N    ~ Multinomial(N_{s,i}; pi_{s,i,1..B}, 1 - p_det)
#   log lambda_{s,i} = X_i     . (mu_lambda + b_lambda_s) [+ f_i]
#                                              [+ sum_q lambda_{s,q} zeta_{q,i}]
#   log sigma_{s,i}  = X_sig_i . (mu_sigma  + b_sigma_s)
#   b_lambda_s ~ N(0, Sigma_lambda),  b_sigma_s ~ N(0, Sigma_sigma)
#
# with pi_b = int_bin g(x; sigma) f(x) dx the per-bin detection probability
# (half-normal / hazard key, line / point transect density f). The latent
# N_{s,i} integrates out in closed form per species-site, exactly as for the
# single-species distance() family, so this file needs no new C++: every fit is
# driven by cpp_distance_site_sweep, which already returns the per-site marginal
# pieces
#
#   log_lik, grad_lam, info_lam, var_N, swl (= score_wt_lambda)
#
# from which the community EM reads its per-species score and the latent driver
# reads its working oracle
#
#   score = grad_lam,  curv = info_lam - var_N * swl^2
#
# -- the (1,1) entry of the per-site marginal observed information
# B_i = diag(info_lam, info_sig) - Var(N_i|y) v v' (Louis 1982), i.e. the
# abundance curvature with the detection arm profiled out. That is the SAME
# formula the community N-mixture uses (R/ms_abun_latent.R): both are count
# marginals with the same Louis structure, so the two families differ only in
# which kernel supplies the pieces.
#
# Under the hazard-rate key the scalar log-shape eta_b is shared across species,
# so it rides the community EM's `global` slot rather than the per-species arms.
#
#   .tobs_build_ms_distance()  data binder -> model_type = "ms_distance"
#   .tobs_ms_distance_engine() per-species sweep closures (single source of truth)
#   .tobs_ms_distance_oracle() the working oracle for the latent driver
#   .tobs_fit_ms_distance()    the community fit (plain / factors / field)
#   build_ms_distance_fit()    pack into a tobs_fit
#   simulate_ms_distance()     community distance simulator
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# `y` is a 3D array [n_sites x n_bins x n_species] or a (named) list of
# n_sites x n_bins per-bin count matrices, one per species. The abundance design
# X_lambda and the detection-scale design X_sigma are shared across species
# (community covariates); the binned counts differ per species. The cut points,
# key and transect geometry are shared.
.tobs_build_ms_distance <- function(abund_formula, det_formula, data, y, species,
                                    cutpoints, key = "halfnorm",
                                    transect = "line", mixture = "poisson",
                                    quad_order = 64L) {
  if (is.list(y) && !is.array(y)) {
    n_species <- length(y)
    species_names <- if (is.character(species)) species
                     else if (!is.null(names(y))) names(y)
                     else paste0("sp", seq_len(n_species))
    n_sites <- nrow(y[[1]]); n_bins <- ncol(y[[1]])
    y_arr <- array(NA_integer_, dim = c(n_sites, n_bins, n_species))
    for (s in seq_len(n_species)) y_arr[, , s] <- as.integer(round(y[[s]]))
    y <- y_arr
  } else {
    if (length(dim(y)) != 3L) {
      stop("y must be a 3D array [n_sites x n_bins x n_species] or a list of ",
           "per-bin count matrices.", call. = FALSE)
    }
    species_names <- if (is.character(species)) species
                     else paste0("sp", seq_len(dim(y)[3]))
    storage.mode(y) <- "integer"
  }
  n_sites <- dim(y)[1]; n_bins <- dim(y)[2]; n_species <- dim(y)[3]
  .tobs_check_site_count(n_sites, nrow(data), "sites")
  if (length(species_names) != n_species) {
    stop(sprintf("species has %d names but y has %d species",
                 length(species_names), n_species), call. = FALSE)
  }
  if (anyNA(y)) {
    stop("ms_distance() needs complete bin counts: y must not contain NA.",
         call. = FALSE)
  }
  if (any(y < 0L)) {
    stop("y must contain nonnegative integer counts.", call. = FALSE)
  }
  if (is.null(cutpoints) || length(cutpoints) != n_bins + 1L) {
    stop(sprintf("cutpoints must have length dim(y)[2] + 1 = %d (the bin edges, ",
                 "0 = c_0 < c_1 < ... < c_B).", n_bins + 1L), call. = FALSE)
  }
  cutpoints <- as.numeric(cutpoints)
  if (any(diff(cutpoints) <= 0) || cutpoints[1] < 0) {
    stop("cutpoints must be strictly increasing and start at >= 0.", call. = FALSE)
  }

  bind     <- .tobs_bind_formulas(list(lambda = abund_formula,
                                       sigma = det_formula), data)
  X_lambda <- model.matrix(bind$fe$lambda, data)
  X_sigma  <- model.matrix(bind$fe$sigma, data)

  structure(list(
    model_type  = "ms_distance",
    y           = y,
    X_processes = list(X_lambda, X_sigma),
    formulas    = list(lambda = bind$fe$lambda, sigma = bind$fe$sigma),
    structured_terms = bind$terms,
    data        = data,
    n_sites     = n_sites,
    n_bins      = n_bins,
    n_species   = n_species,
    species_names = species_names,
    cutpoints   = cutpoints,
    key         = key,
    transect    = transect,
    mixture     = mixture,
    quad_order  = as.integer(quad_order),
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda), coef_names = colnames(X_lambda),
           link = "log"),
      list(name = "sigma",  p = ncol(X_sigma),  coef_names = colnames(X_sigma),
           link = "log")
    )
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Per-species marginal engine
# ---------------------------------------------------------------------------

# One closure per species over that species' bin counts, sharing the cut points /
# key / transect. `sweep(s, eta_lam, eta_sig, eta_b)` returns the raw
# cpp_distance_site_sweep list: the per-site marginal log-likelihood, the
# eta-level scores, and the (info_lam, var_N, swl) pieces the marginal curvature
# is built from. Single source of truth for the community EM and the latent
# driver's oracle.
.tobs_ms_distance_engine <- function(model, K_max = NULL) {
  S   <- model$n_species
  cut <- as.numeric(model$cutpoints)
  tc  <- .dist_transect_code(model$transect)
  kc  <- .dist_key_code(model$key)
  qo  <- as.integer(model$quad_order %||% 64L)
  ys  <- lapply(seq_len(S), function(s)
    matrix(as.integer(model$y[, , s]), model$n_sites, model$n_bins))
  R_max <- max(vapply(ys, function(m) if (length(m)) max(rowSums(m)) else 0, 0))
  K_max <- if (is.null(K_max)) as.integer(3L * R_max + 100L) else as.integer(K_max)
  list(
    y_s = ys, K_max = K_max, hazard = identical(model$key, "hazard"),
    sweep = function(s, eta_lam, eta_sig, eta_b = 0) {
      cpp_distance_site_sweep(ys[[s]], eta_lam, eta_sig, cut, tc, qo, K_max,
                              nb = FALSE, r = Inf, key = kc,
                              eta_b = as.numeric(eta_b))
    })
}


# Working oracle at fixed detection-scale predictors: the score and curvature of
# each species' binned-distance marginal with respect to an additive offset on
# its abundance log-predictor.
.tobs_ms_distance_oracle <- function(eng, eta_sig_list, eta_b, Ns, S) {
  eval_all <- function(eta) lapply(seq_len(S), function(s)
    eng$sweep(s, eta[, s], eta_sig_list[[s]], eta_b))
  list(
    n_sites = Ns, n_species = S,
    working = function(eta) {
      sws <- eval_all(eta)
      list(
        score = vapply(sws, function(sw) as.numeric(sw$grad_lam), numeric(Ns)),
        curv  = vapply(sws, function(sw)
          pmax(as.numeric(sw$info_lam - sw$var_N * sw$swl^2), 1e-8), numeric(Ns)))
    },
    data_ll = function(eta) sum(vapply(eval_all(eta),
                                       function(sw) sum(sw$log_lik), 0)))
}


# ---------------------------------------------------------------------------
# Fitter
# ---------------------------------------------------------------------------

# The community distance fit. `latent` / `spatial` are optional: with neither the
# community EM runs once (msDS); with either, the shared block-coordinate driver
# alternates the EM (latent as a per-species offset on log lambda) with the field
# / factor updates (lfMsDS / sfMsDS). One fitter, so the plain community model is
# the same code path as the latent ones.
.tobs_fit_ms_distance <- function(model, spatial = NULL, latent = NULL,
                                  mixture = "poisson", K_max = NULL,
                                  max.iter = 100L, tol = 1e-4,
                                  sigma.beta = 5, priors = NULL,
                                  max.outer = 25L, verbose = FALSE, ...) {
  if (!identical(mixture, "poisson")) {
    stop("ms_distance() is Poisson-only: the negative-binomial size is a ",
         "per-site dispersion that the community distance fit does not yet ",
         "carry as a per-species random effect (gcol33/tulpaObs#117). Use ",
         "mixture = \"poisson\".", call. = FALSE)
  }
  S  <- model$n_species
  Ns <- model$n_sites
  eng <- .tobs_ms_distance_engine(model, K_max = K_max)
  X_lam <- model$X_processes[[1L]]; X_sig <- model$X_processes[[2L]]
  p_lam <- ncol(X_lam); p_sig <- ncol(X_sig)
  P <- p_lam + p_sig
  lam_idx <- seq_len(p_lam); sig_idx <- p_lam + seq_len(p_sig)
  arm_idx <- list(lambda = lam_idx, sigma = sig_idx)
  # Under the hazard-rate key the scalar log-shape is shared across species, so
  # it is a community `global` rather than a per-species arm coordinate.
  hazard <- eng$hazard
  G <- if (hazard) 1L else 0L

  # Warm start: fit each species on its own (the single-species distance Laplace)
  # and pool. The distance marginal trades lambda against sigma from a naive
  # start, and this is the same kernel the community EM then refines.
  warm <- lapply(seq_len(S), function(s) tryCatch(
    distance_laplace(y = eng$y_s[[s]], X_lambda = X_lam, X_sigma = X_sig,
                     cutpoints = model$cutpoints, key = model$key,
                     transect = model$transect, mixture = "P",
                     K_max = eng$K_max, quad_order = model$quad_order,
                     max_iter = as.integer(max.iter), tol = 1e-6, verbose = FALSE),
    error = function(e) NULL))
  ok <- !vapply(warm, is.null, logical(1))
  mu0 <- if (any(ok)) {
    c(colMeans(do.call(rbind, lapply(warm[ok], function(w) w$beta_lambda))),
      colMeans(do.call(rbind, lapply(warm[ok], function(w) w$beta_sigma))))
  } else {
    c(log(max(mean(vapply(eng$y_s, function(m) mean(rowSums(m)), 0)), 0.5) / 0.5),
      rep(0, p_lam - 1L),
      log(max(model$cutpoints) / 2), rep(0, p_sig - 1L))
  }
  g0 <- if (hazard) {
    gb <- vapply(warm[ok], function(w) w$eta_b %||% NA_real_, 0)
    if (any(is.finite(gb))) mean(gb[is.finite(gb)]) else 0
  } else numeric(0)

  em_fit <- function(site_off, fac_off, em_prev) {
    eta_of <- function(s, theta) list(
      lam = as.numeric(X_lam %*% theta[lam_idx]) + site_off + fac_off[, s],
      sig = as.numeric(X_sig %*% theta[sig_idx]))
    sp_ll <- function(s, theta, global) {
      e <- eta_of(s, theta)
      sum(eng$sweep(s, e$lam, e$sig, if (hazard) global[1L] else 0)$log_lik)
    }
    sp_grad <- function(s, theta, global) {
      e  <- eta_of(s, theta)
      sw <- eng$sweep(s, e$lam, e$sig, if (hazard) global[1L] else 0)
      c(as.numeric(crossprod(X_lam, sw$grad_lam)),
        as.numeric(crossprod(X_sig, sw$grad_sig)),
        if (hazard) sw$grad_b else numeric(0))
    }
    # Each sweep sums over the latent N, so a cold restart every outer pass
    # dominates the fit: resume the whole state from the previous pass.
    .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      init_mu     = if (is.null(em_prev)) mu0 else em_prev$mu,
      init_global = if (is.null(em_prev)) g0  else em_prev$global,
      init_b      = em_prev$b_list,
      init_Sigma  = em_prev$Sigma,
      penalize_global = FALSE,
      sigma_beta = sigma.beta, priors = priors, sigma_init = 0.3,
      max_iter = min(as.integer(max.iter), 40L), tol = as.numeric(tol),
      newton_max = 20L, verbose = FALSE)
  }
  offset_of <- function(em) {
    vapply(seq_len(S), function(s)
      as.numeric(X_lam %*% (em$mu + em$b_list[[s]])[lam_idx]), numeric(Ns))
  }
  make_oracle <- function(em) {
    eta_sig_list <- lapply(seq_len(S), function(s)
      as.numeric(X_sig %*% (em$mu + em$b_list[[s]])[sig_idx]))
    .tobs_ms_distance_oracle(eng, eta_sig_list,
                             if (hazard) em$global[1L] else 0, Ns, S)
  }

  # With no latent structure the block-coordinate ascent degenerates to its EM
  # block, so run that directly rather than looping it against nothing.
  res <- if (is.null(spatial) && is.null(latent)) {
    list(em = em_fit(numeric(Ns), matrix(0, Ns, S), NULL),
         field = NULL, factor = NULL, geom = NULL)
  } else {
    .tobs_community_latent_ascent(
      spatial = spatial, latent = latent, model = model, what = "ms_distance()",
      make_oracle = make_oracle, em_fit = em_fit, offset_of = offset_of,
      allow = c("icar", "car_proper", "bym2", "spde"),
      tol = tol, max.outer = max.outer, verbose = verbose)
  }

  fit <- build_ms_distance_fit(res$em, model, lam_idx, sig_idx, hazard)
  fit <- .tobs_latent_attach_field(fit, res, spatial, "distance_field_offset")
  fit <- .tobs_latent_attach_factor(fit, res, latent, model,
                                    "distance_factor_offset")
  fit
}


# ---------------------------------------------------------------------------
# Wrap into a tobs_fit
# ---------------------------------------------------------------------------

build_ms_distance_fit <- function(em, model, lam_idx, sig_idx, hazard = FALSE) {
  pi_list <- model$process_info
  lam_nm  <- pi_list[[1]]$coef_names
  sig_nm  <- pi_list[[2]]$coef_names
  nms <- c(paste0("lambda_", lam_nm), paste0("sigma_", sig_nm))
  means <- as.numeric(em$mu)
  if (hazard) { means <- c(means, as.numeric(em$global[1L])); nms <- c(nms, "log_b") }
  names(means) <- nms

  vcov <- as.matrix(em$Vf)
  if (nrow(vcov) != length(nms)) {
    stop(sprintf("build_ms_distance_fit(): vcov is %dx%d but there are %d ",
                 nrow(vcov), ncol(vcov), length(nms)),
         "coefficient names.", call. = FALSE)
  }
  rownames(vcov) <- colnames(vcov) <- nms
  sds <- sqrt(pmax(diag(vcov), 0)); names(sds) <- nms

  Sigma_lambda <- as.matrix(em$Sigma$lambda)
  Sigma_sigma  <- as.matrix(em$Sigma$sigma)
  dimnames(Sigma_lambda) <- list(lam_nm, lam_nm)
  dimnames(Sigma_sigma)  <- list(sig_nm, sig_nm)

  blup_lambda <- do.call(rbind, lapply(em$b_list, function(b) b[lam_idx]))
  blup_sigma  <- do.call(rbind, lapply(em$b_list, function(b) b[sig_idx]))
  coef_lambda <- sweep(blup_lambda, 2, em$mu[lam_idx], "+")
  coef_sigma  <- sweep(blup_sigma,  2, em$mu[sig_idx], "+")
  rownames(coef_lambda) <- rownames(coef_sigma) <- model$species_names
  rownames(blup_lambda) <- rownames(blup_sigma) <- model$species_names
  colnames(coef_lambda) <- colnames(blup_lambda) <- lam_nm
  colnames(coef_sigma)  <- colnames(blup_sigma)  <- sig_nm

  n_pseudo <- 1000L
  draws <- .rmvn(n_pseudo, means, vcov)
  colnames(draws) <- nms

  structure(list(
    draws = draws, means = means, sds = sds, vcov = vcov,
    n_samples = n_pseudo, n_params = length(means),
    col_names = nms, param_names = nms,
    n_fixed = length(means), fixed_names = nms,
    process_info = pi_list,
    model = model, spatial = NULL, spatial_field = NULL,
    method = "laplace", mixture = model$mixture %||% "poisson",
    log_lik = em$logML %||% NA_real_,
    ms_community = list(
      Sigma_lambda = Sigma_lambda, Sigma_sigma = Sigma_sigma,
      sd_lambda = sqrt(pmax(diag(Sigma_lambda), 0)),
      sd_sigma  = sqrt(pmax(diag(Sigma_sigma), 0)),
      coef_lambda = coef_lambda, coef_sigma = coef_sigma,
      blup_lambda = blup_lambda, blup_sigma = blup_sigma),
    convergence = list(converged = isTRUE(em$converged),
                       n_iter = em$n_iter %||% NA_integer_)
  ), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed to from methods.R by model_type == "ms_distance")
# ---------------------------------------------------------------------------

.tobs_ranef_ms_distance <- function(object) {
  .tobs_ranef_ms_long(object$ms_community,
                      c(lambda = "blup_lambda", sigma = "blup_sigma"))
}

# Per-species site-level expected abundance lambda_{s,i} and detection scale
# sigma_{s,i} at the posterior-mean (mu + BLUP) coefficients, each an
# [n_sites x n_species] matrix. The field / factor offsets enter log lambda.
.tobs_fitted_ms_distance <- function(object) {
  model <- object$model
  cm    <- object$ms_community
  eta_lambda <- model$X_processes[[1L]] %*% t(cm$coef_lambda)
  fld <- model$distance_field_offset
  if (!is.null(fld)) {
    eta_lambda <- sweep(eta_lambda, 1, as.numeric(fld), "+")
  }
  if (!is.null(model$distance_factor_offset)) {
    eta_lambda <- eta_lambda + model$distance_factor_offset
  }
  lambda <- exp(eta_lambda)
  sigma  <- exp(model$X_processes[[2L]] %*% t(cm$coef_sigma))
  dimnames(lambda) <- dimnames(sigma) <- list(NULL, model$species_names)
  list(lambda = lambda, sigma = sigma)
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate community (multispecies) binned distance-sampling data
#'
#' Per-species binned distance sampling with Gaussian community hyperpriors:
#' `beta_lambda_s ~ N(mu_lambda, diag(sd_lambda^2))`,
#' `beta_sigma_s ~ N(mu_sigma, diag(sd_sigma^2))`, then
#' `N_{s,i} ~ Poisson(lambda_{s,i})` and the detected individuals allocated to
#' distance bins by the half-normal detection function. The returned `y` is a 3D
#' array `[n_sites x n_bins x n_species]` suitable for [tobs()] with
#' [ms_distance()].
#'
#' @param n_species Number of species (default 10).
#' @param N Number of sites (default 100).
#' @param cutpoints Distance-bin edges (default `c(0, 25, 50, 75, 100)`).
#' @param transect Transect geometry: `"line"` (default) or `"point"`.
#' @param key Detection function: `"halfnorm"` (default) or `"hazard"`.
#' @param shape Hazard-rate log-shape, shared across species (ignored under the
#'   half-normal key). Default 0.
#' @param n_abund_covs,n_det_covs Number of abundance / detection-scale
#'   covariates (default 1 and 0).
#' @param mu_lambda Community-mean abundance coefficients on the log scale.
#'   Default `c(log(30), rep(0.4, n_abund_covs))`.
#' @param mu_sigma Community-mean detection-scale coefficients on the log scale.
#'   Default `c(log(40), rep(0, n_det_covs))`.
#' @param sd_lambda,sd_sigma Per-coefficient community SDs. Default 0.4 and 0.2.
#' @param n_factors If `> 0`, add `n_factors` per-site latent factors with
#'   per-species loadings to `log lambda` (the lfMsDS truth). Default 0.
#' @param load_sd SD of the factor loadings (default 0.5).
#' @param field Optional length-`N` shared spatial field added to every species'
#'   `log lambda`. When given alongside `n_factors > 0` the loadings are centred
#'   across species, so the field owns the shared spatial mean.
#' @param seed Optional random seed.
#' @return A list with `y`, `data`, `species`, `cutpoints`, and `truth`
#'   (community means / SDs, per-species coefficients, `lambda`, `sigma`, and --
#'   with factors -- the `loadings`, `factors` and the implied residual
#'   correlation `cor_res`). The latent `N` is drawn inside the shared C++
#'   simulator and is not returned.
#' @export
simulate_ms_distance <- function(n_species = 10, N = 100,
                                 cutpoints = c(0, 25, 50, 75, 100),
                                 transect = c("line", "point"),
                                 key = c("halfnorm", "hazard"), shape = 0,
                                 n_abund_covs = 1, n_det_covs = 0,
                                 mu_lambda = NULL, mu_sigma = NULL,
                                 sd_lambda = 0.4, sd_sigma = 0.2,
                                 n_factors = 0, load_sd = 0.5,
                                 field = NULL, seed = NULL) {
  transect <- match.arg(transect)
  key      <- match.arg(key)
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(field)) N <- length(field)
  if (is.null(mu_lambda)) mu_lambda <- c(log(30), rep(0.4, n_abund_covs))
  if (is.null(mu_sigma))  mu_sigma  <- c(log(40), rep(0, n_det_covs))
  p_lam <- length(mu_lambda); p_sig <- length(mu_sigma)
  sd_lambda <- if (length(sd_lambda) == 1L) rep(sd_lambda, p_lam) else sd_lambda
  sd_sigma  <- if (length(sd_sigma)  == 1L) rep(sd_sigma,  p_sig) else sd_sigma

  make_covs <- function(n_covs, prefix) {
    if (n_covs <= 0L) return(data.frame(row.names = seq_len(N)))
    m <- matrix(stats::rnorm(N * n_covs), N, n_covs)
    df <- as.data.frame(m); names(df) <- paste0(prefix, seq_len(n_covs)); df
  }
  abund_covs <- make_covs(n_abund_covs, "abund_cov")
  det_covs   <- make_covs(n_det_covs,   "det_cov")
  data <- data.frame(row.names = seq_len(N))
  if (ncol(abund_covs)) data <- cbind(data, abund_covs)
  if (ncol(det_covs))   data <- cbind(data, det_covs)
  design_of <- function(df) {
    if (ncol(df)) stats::model.matrix(~ ., df)
    else stats::model.matrix(~ 1, data.frame(row.names = seq_len(N)))
  }
  X_lambda <- design_of(abund_covs)
  X_sigma  <- design_of(det_covs)

  beta_lambda <- matrix(stats::rnorm(n_species * p_lam, 0,
                                     rep(sd_lambda, each = n_species)),
                        n_species, p_lam) +
                 matrix(mu_lambda, n_species, p_lam, byrow = TRUE)
  beta_sigma <- matrix(stats::rnorm(n_species * p_sig, 0,
                                    rep(sd_sigma, each = n_species)),
                       n_species, p_sig) +
                matrix(mu_sigma, n_species, p_sig, byrow = TRUE)

  f <- if (is.null(field)) numeric(N) else as.numeric(field)
  Q <- as.integer(n_factors)
  loadings <- NULL; factors <- NULL; cor_res <- NULL
  fac_off <- matrix(0, N, n_species)
  if (Q > 0L) {
    loadings <- matrix(stats::rnorm(n_species * Q, 0, load_sd), n_species, Q)
    if (!is.null(field)) loadings <- scale(loadings, scale = FALSE)
    factors  <- matrix(stats::rnorm(N * Q), N, Q)
    fac_off  <- factors %*% t(loadings)
    cor_res  <- stats::cov2cor(tcrossprod(loadings) + diag(1e-8, n_species))
  }

  B <- length(cutpoints) - 1L
  species_names <- paste0("sp", seq_len(n_species))
  y <- array(0L, dim = c(N, B, n_species),
             dimnames = list(NULL, NULL, species_names))
  lambda <- matrix(NA_real_, N, n_species)
  sigma  <- matrix(NA_real_, N, n_species)
  # Draw through cpp_simulate_distance, the same kernel the likelihood
  # integrates against (src/distance_quad.h): a separate R-side quadrature for
  # the per-bin probabilities would simulate from a pi the model is not fit
  # against, and bias every recovery check. The shared field / factor offset is
  # folded into the abundance design as a column with coefficient 1, so
  # log lambda = X beta_s + f + sum_q l_{s,q} z_{q,i} without touching the
  # kernel. A one-row `draws` matrix pins the coefficients (no draw selection).
  tc <- if (identical(transect, "point")) 1L else 0L
  kc <- .dist_key_code(key)
  for (s in seq_len(n_species)) {
    off_s  <- f + fac_off[, s]
    has_off <- any(off_s != 0)
    Xl_s   <- if (has_off) cbind(X_lambda, .offset = off_s) else X_lambda
    beta_s <- if (has_off) c(beta_lambda[s, ], 1) else beta_lambda[s, ]
    res <- cpp_simulate_distance(
      Xl_s, X_sigma,
      matrix(c(beta_s, beta_sigma[s, ]), nrow = 1L),
      as.numeric(cutpoints), kc, tc, as.numeric(shape),
      N, B, ncol(Xl_s), p_sig, FALSE, NA_real_, 1L)
    y[, , s] <- res[[1L]]
    lambda[, s] <- exp(as.numeric(X_lambda %*% beta_lambda[s, ]) + off_s)
    sigma[, s]  <- exp(as.numeric(X_sigma %*% beta_sigma[s, ]))
  }

  list(
    y = y, data = data, species = species_names, cutpoints = cutpoints,
    truth = list(
      mu_lambda = mu_lambda, mu_sigma = mu_sigma,
      sd_lambda = sd_lambda, sd_sigma = sd_sigma,
      beta_lambda = beta_lambda, beta_sigma = beta_sigma,
      lambda = lambda, sigma = sigma,
      field = if (!is.null(field)) f else NULL,
      loadings = loadings, factors = factors, cor_res = cor_res,
      transect = transect, key = key, shape = shape)
  )
}
