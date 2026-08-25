# =============================================================================
# ms_distance.R
# - community / multispecies binned distance sampling (the spAbundance
# msDS analogue), with optional latent factors (lfMsDS) and a shared
# field (the spatial-factor case, sfMsDS).. Poisson.
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
    stop(sprintf(paste0("cutpoints must have length dim(y)[2] + 1 = %d (the bin edges, ",
                        "0 = c_0 < c_1 < ... < c_B)."), n_bins + 1L), call. = FALSE)
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
.tobs_ms_distance_engine <- function(model, K_max = NULL, headroom = NULL) {
  S   <- model$n_species
  cut <- as.numeric(model$cutpoints)
  tc  <- .dist_transect_code(model$transect)
  kc  <- .dist_key_code(model$key)
  qo  <- as.integer(model$quad_order %||% 64L)
  # Built ONCE for the whole fit (every species, every EM/block-coordinate
  # iteration) rather than Newton-Raphson root-finding the Gauss-Legendre nodes
  # fresh on every `sweep()` call -- `sweep()` is the community EM's inner-loop
  # oracle, called far more than once per species.
  qptr <- cpp_distance_build_quad(cut, tc, qo)
  ys  <- lapply(seq_len(S), function(s)
    matrix(as.integer(model$y[, , s]), model$n_sites, model$n_bins))
  site_tot <- unlist(lapply(ys, rowSums))
  trunc    <- .dist_truncation(K_max, site_tot)
  K_max    <- trunc$K_max
  headroom <- if (is.null(headroom)) trunc$headroom else as.integer(headroom)
  list(
    y_s = ys, K_max = K_max, headroom = headroom, key_code = kc,
    hazard = identical(model$key, "hazard"),
    quad_xptr = qptr,
    # `idx` (a subset of site rows) lets a caller sweep only those sites --
    # `eta_lam` / `eta_sig` must already be subsetted to the same rows by the
    # caller; `ys[[s]]` is subsetted here since it is this closure's own fixed
    # data. `value_only = TRUE` skips every gradient/Fisher computation in the
    # compiled kernel and returns just `log_lik` -- for a caller that only reads
    # `sw$log_lik` (the oracle's `ll_cell` / `data_ll`), not `working()`, which
    # needs the full sweep.
    sweep = function(s, eta_lam, eta_sig, eta_b = 0, idx = NULL, value_only = FALSE) {
      y_use <- if (is.null(idx)) ys[[s]] else ys[[s]][idx, , drop = FALSE]
      cpp_distance_site_sweep(y_use, eta_lam, eta_sig, qptr, K_max,
                              nb = FALSE, r = Inf, key = kc,
                              eta_b = as.numeric(eta_b), value_only = value_only,
                              headroom = headroom)
    })
}


# Working oracle at fixed detection-scale predictors: the score and curvature of
# each species' binned-distance marginal with respect to an additive offset on
# its abundance log-predictor.
.tobs_ms_distance_oracle <- function(eng, eta_sig_list, eta_b, Ns, S) {
  # `idx` (a subset of site rows) restricts the sweep to those sites --
  # `eta_sig_list[[s]]` and `eta_b` (fixed at this oracle's construction) are
  # subsetted here to match, and `eng$sweep`'s own `ys[[s]]` subsets itself.
  # `eta_b` is a scalar (0, or the hazard shape parameter) in every caller, so
  # it needs no subsetting.
  eval_all <- function(eta, idx = NULL, value_only = FALSE) {
    ii <- idx %||% seq_len(Ns)
    lapply(seq_len(S), function(s)
      eng$sweep(s, eta[ii, s], eta_sig_list[[s]][ii], eta_b, idx = ii,
               value_only = value_only))
  }
  list(
    n_sites = Ns, n_species = S,
    working = function(eta) {
      sws <- eval_all(eta)
      list(
        score = vapply(sws, function(sw) as.numeric(sw$grad_lam), numeric(Ns)),
        curv  = vapply(sws, function(sw)
          pmax(as.numeric(sw$info_lam - sw$var_N * sw$swl^2), 1e-8), numeric(Ns)))
    },
    # value_only = TRUE: ll_cell only ever reads log_lik.
    ll_cell = function(eta, idx = NULL) {
      nn <- length(idx %||% seq_len(Ns))
      # vapply degenerates to a plain vector (not an nn x S matrix) when
      # nn == 1, since a length-1 FUN.VALUE never triggers its matrix path --
      # only pending on one site is common in the mode-adaptation backtracking
      # tail, so this must be forced back into a matrix rather than assumed.
      matrix(vapply(eval_all(eta, idx, value_only = TRUE),
                    function(sw) as.numeric(sw$log_lik),
                    numeric(nn)), nrow = nn)
    },
    data_ll = function(eta) sum(vapply(eval_all(eta, value_only = TRUE),
                                       function(sw) sum(sw$log_lik), 0)))
}


# Per-species marginal observed information from one sweep, assembled rather than
# finite-differenced. The per-site Louis (1982) block is
#
#   B_i = diag(info_lam_i, info_sig_i) - Var(N_i|y) v_i v_i',
#         v_i = (-swl_i, -vN_sig_i)
#
# the same block nmix_site_marginal()'s obs_info_block() forms for the community
# N-mixture (R/ms_abun_latent.R); only the kernel supplying the pieces differs.
# The second component takes a further sign flip because the distance kernel
# stores `vN_d` already negated (-p_k/(1-p), src/distance_kernel.h), where the
# N-mixture's v carries +p. That is invisible on the diagonal and wrong by a
# factor of ~2 on the cross block, so it is asserted against the finite
# difference rather than reasoned from the sibling family
# (test-ms-distance-info.R).
# The complete-data log-likelihood splits as log P(N_i | lambda) +
# log P(y_i | N_i, sigma), so its expected information carries no lambda-sigma
# cross term and the whole off-diagonal comes from the rank-1 score-variance
# subtraction. Distance carries ONE unit per site (the bins sit inside it), so
# D_i is the 2 x P row [X_lam[i, ] | 0 ; 0 | X_sig[i, ]] and sum_i D_i' B_i D_i
# is three crossprods against diagonal weights, with no per-site loop -- the
# N-mixture needs one only because its visit count varies by site.
#
# Half-normal key only: under the hazard key the shared log-shape is a community
# global, and the sweep returns its score and information already summed over
# sites, with no per-site cross terms to sandwich.
.tobs_ms_distance_info_block <- function(sw, X_lam, X_sig, lam_idx, sig_idx, P) {
  # diag(info) - Var(N|y) v v', expanded. The two diagonal entries cannot see the
  # relative sign inside v; the off-diagonal is where it lands.
  b11 <- sw$info_lam     - sw$var_N * sw$swl^2
  b22 <- sw$info_sig_obs - sw$var_N * sw$vN_sig^2
  b12 <-                 - sw$var_N * sw$swl * sw$vN_sig
  I <- matrix(0, P, P)
  I[lam_idx, lam_idx] <- crossprod(X_lam, X_lam * b11)
  I[sig_idx, sig_idx] <- crossprod(X_sig, X_sig * b22)
  cross <- crossprod(X_lam, X_sig * b12)
  I[lam_idx, sig_idx] <- cross
  I[sig_idx, lam_idx] <- t(cross)
  (I + t(I)) / 2
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
                                  headroom = NULL,
                                  max.iter = 100L, tol = 1e-4,
                                  sigma.beta = 5, priors = NULL,
                                  max.outer = NULL, factor.starts = NULL,
                                  n.quad = NULL, verbose = FALSE, ...) {
  if (!identical(mixture, "poisson")) {
    stop("ms_distance() is Poisson-only: the negative-binomial size is a ",
         "per-site dispersion that the community distance fit does not yet ",
         "carry as a per-species random effect. Use ",
         "mixture = \"poisson\".", call. = FALSE)
  }
  S  <- model$n_species
  Ns <- model$n_sites
  eng <- .tobs_ms_distance_engine(model, K_max = K_max, headroom = headroom)
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
                     K_max = eng$K_max, headroom = eng$headroom,
                     quad_order = model$quad_order,
                     max_iter = as.integer(max.iter), tol = 1e-6, verbose = FALSE,
                     quad_xptr = eng$quad_xptr),
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
    # The community EM otherwise central-differences this, at 2U sweeps per
    # species per Newton step, and every sweep sums over the latent N.
    sp_info <- if (hazard) NULL else function(s, theta, global) {
      e <- eta_of(s, theta)
      .tobs_ms_distance_info_block(eng$sweep(s, e$lam, e$sig, 0),
                                   X_lam, X_sig, lam_idx, sig_idx, P)
    }
    # Each sweep sums over the latent N, so a cold restart every outer pass
    # dominates the fit: resume the whole state from the previous pass.
    .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      sp_info = sp_info,
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
      tol = tol, max.outer = max.outer, factor.starts = factor.starts,
      n.quad = n.quad, verbose = verbose)
  }

  # Guard the per-site truncation at the converged community coefficients: the
  # field / factor offsets the block-coordinate ascent converged to
  # (`res$field$site_off`, `res$factor$offset`) are the SAME ones `eta_of()`
  # inside `em_fit()` folds into eta_lambda, so the check runs at the
  # predictor the fit actually used. Escalates to the shared ceiling on
  # disagreement and re-fits the whole community model, exactly as the
  # single-species and areal-spatial guards do.
  if (eng$headroom >= 0L) {
    site_off_mode <- if (!is.null(res$field)) res$field$site_off else NULL
    fac_off_mode  <- if (!is.null(res$factor)) res$factor$offset else NULL
    gap <- tryCatch(
      .dist_community_score_gap(eng, res$em, X_lam, X_sig, lam_idx, sig_idx,
                                hazard, S, site_off = site_off_mode,
                                fac_off = fac_off_mode),
      error = function(e) NA_real_)
    if (is.finite(gap) && gap > .NMIX_SCORE_TOL) {
      h_next <- .nmix_widen_headroom(eng$headroom, eng$K_max)
      if (!is.null(h_next)) {
        cl <- match.call()
        cl$headroom <- h_next
        return(eval(cl, parent.frame()))
      }
    }
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
      blup_lambda = blup_lambda, blup_sigma = blup_sigma,
      # Per-species posterior covariance Cov(b_s|y) (Louis 1982, from the
      # community EM's own Newton solve, conditional on the converged
      # community mean) -- what a per-species-coefficient consumer (SBC's
      # "rank a fixed species set" design, a calibrated per-species CI) needs
      # beyond the point BLUP; not previously exposed on the fit object. Bf =
      # the (mu,global)-b_s cross-Hessian block from the same Newton solve:
      # mu/global and b_s are NOT independent in the posterior, and Bf is
      # what lets a consumer draw them jointly instead -- see
      # .tobs_sbc_community_b_draws (R/sbc.R).
      Cinv = em$Cinv, Bf = em$Bf),
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

# Pointwise log-likelihood [n_draws x (n_species * n_sites)] for WAIC / LOO /
# DIC / CPO: per species, the exact binned-distance marginal
# (cpp_distance_site_sweep, via the SAME .tobs_ms_distance_engine() the fit
# used) scored over the community-mean pseudo-draws with the per-species BLUP
# plugged into both arms, the same shape .tobs_ploglik_community_occu_cover()
# uses. K_max / headroom are rebuilt from the data at their fit-time defaults,
# the convention .tobs_ploglik_distance() already uses for the single-species
# family (the fit's own converged truncation is not carried on the model).
.tobs_ploglik_ms_distance <- function(object, n.draws = 1000L, n.threads = 1L) {
  model <- object$model
  cm    <- object$ms_community
  draws <- object$draws
  if (!is.null(n.draws) && n.draws < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  eng    <- .tobs_ms_distance_engine(model)
  hazard <- eng$hazard
  p_lam  <- model$process_info[[1L]]$p; p_sig <- model$process_info[[2L]]$p
  lam_cols <- seq_len(p_lam); sig_cols <- p_lam + seq_len(p_sig)
  b_col  <- if (hazard) p_lam + p_sig + 1L else NA_integer_
  X_lam  <- model$X_processes[[1L]]; X_sig <- model$X_processes[[2L]]
  field_off  <- model$distance_field_offset
  factor_off <- model$distance_factor_offset
  M <- nrow(draws); n_sites <- model$n_sites
  cols <- vector("list", model$n_species)
  for (s in seq_len(model$n_species)) {
    ll_s <- matrix(0, M, n_sites)
    fo_s <- if (is.null(factor_off)) NULL else factor_off[, s]
    for (m in seq_len(M)) {
      eta_lam <- as.numeric(X_lam %*% (draws[m, lam_cols] + cm$blup_lambda[s, ]))
      if (!is.null(field_off)) eta_lam <- eta_lam + as.numeric(field_off)
      if (!is.null(fo_s))      eta_lam <- eta_lam + fo_s
      eta_sig <- as.numeric(X_sig %*% (draws[m, sig_cols] + cm$blup_sigma[s, ]))
      eta_b   <- if (hazard) draws[m, b_col] else 0
      ll_s[m, ] <- eng$sweep(s, eta_lam, eta_sig, eta_b,
                             value_only = TRUE)$log_lik
    }
    cols[[s]] <- ll_s
  }
  do.call(cbind, cols)
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

# residuals() for the community distance sampling fit. Per species the bin-b
# count is marginally y_ibs ~ Poisson(lambda_is pi_ibs) under Poisson abundance
# (NB(r_s, lambda pi) under negbin), the same marginal the single-species
# .tobs_residuals_distance() scores, with the detection scale sigma_is driving
# each site's own bin probabilities. A count family has no state-level residual,
# so `occ` is NULL.
.tobs_residuals_ms_distance <- function(object, type) {
  model  <- object$model
  fitv   <- .tobs_fitted_ms_distance(object)
  shape  <- object$distance_shape$shape
  y      <- model$y
  size_s <- .tobs_ms_count_size(object, model$n_species)
  out <- array(NA_real_, dim = dim(y),
               dimnames = list(NULL, NULL, model$species_names))
  for (s in seq_len(model$n_species)) {
    pi_mat <- t(vapply(fitv$sigma[, s], function(sg)
      .distance_pi(sg, model$cutpoints, model$key, model$transect, shape),
      numeric(model$n_bins)))
    mu <- pmax(fitv$lambda[, s] * pi_mat, 1e-10)
    out[, , s] <- .tobs_count_residual(y[, , s], mu, type, size = size_s[[s]],
                                       eps = 1e-10)
  }
  list(occ = NULL, det = out)
}

# Posterior replicate binned-count arrays: per species, draw the latent N and
# allocate detections to distance bins through `cpp_simulate_distance` -- the
# SAME kernel the likelihood integrates against (as `simulate_ms_distance()`'s
# own docstring notes: a separate R-side quadrature would simulate from a pi
# the model is not fit against). Reads `ms_community$coef_lambda`/`coef_sigma`
# (deterministic per-species matrices), not `object$draws`, matching
# ms_occu()/ms_count()'s own `simulate()` handlers. `key = "halfnorm"` only --
# the hazard key's log-shape is a community `global` scalar not carried by
# `ms_community`, a documented follow-up.
.tobs_simulate_ms_distance <- function(object, nsim = 1) {
  model <- object$model
  if (!identical(model$key, "halfnorm")) {
    stop("simulate() for ms_distance() fits is only wired for key = ",
         "\"halfnorm\" (the hazard key's log-shape is a community global not ",
         "carried by ms_community).", call. = FALSE)
  }
  cm <- object$ms_community
  n_sites <- model$n_sites; n_bins <- model$n_bins; n_species <- model$n_species
  X_lambda <- model$X_processes[[1L]]; X_sigma <- model$X_processes[[2L]]
  p_lam <- ncol(X_lambda); p_sig <- ncol(X_sigma)
  kc <- .dist_key_code(model$key)
  tc <- .dist_transect_code(model$transect)
  one <- function() {
    y <- array(0L, dim = c(n_sites, n_bins, n_species),
              dimnames = list(NULL, NULL, model$species_names))
    for (s in seq_len(n_species)) {
      res <- cpp_simulate_distance(
        X_lambda, X_sigma,
        matrix(c(cm$coef_lambda[s, ], cm$coef_sigma[s, ]), nrow = 1L),
        as.numeric(model$cutpoints), kc, tc,
        as.integer(model$quad_order %||% 64L), 0,
        n_sites, n_bins, p_lam, p_sig,
        FALSE, NA_real_, 1L)
      y[, , s] <- res[[1L]]
    }
    y
  }
  if (nsim == 1L) one() else lapply(seq_len(nsim), function(i) one())
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
#' @param quad_order Gauss-Legendre nodes per bin used to integrate the per-bin
#'   detection probabilities (default 64, matching [ms_distance()]). Set it to
#'   the `quad_order` the model will be fit at: the rule is
#'   `(cutpoints, transect, quad_order)`, so a different order integrates a
#'   different pi and the data would come from a model the fit does not use.
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
                                 field = NULL, quad_order = 64L, seed = NULL) {
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
      as.numeric(cutpoints), kc, tc, as.integer(quad_order),
      as.numeric(shape),
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
