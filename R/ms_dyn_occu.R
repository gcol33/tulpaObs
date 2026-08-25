# =============================================================================
# ms_dyn_occu.R - community / multispecies DYNAMIC (multi-season) occupancy
#
# The community version of dyn_occu(): a per-species dynamic (HMM) occupancy
# model with Gaussian community hyperpriors on the per-species first-season
# occupancy and detection coefficients, and SHARED (community-wide) season
# transition coefficients.
#
#   z_{s,i,1}        ~ Bernoulli(psi1_{s,i})                    (season-1 state)
#   z_{s,i,t}|z,..   ~ transition(gamma_i, eps_i)               (t = 2..T)
#   y_{s,i,t,j}|z=1  ~ Bernoulli(p_{s,i})                       (detection)
#   logit psi1_{s,i} = X_psi1_i . (mu_psi1 + b_psi1_s)
#   logit p_{s,i}    = X_p_i    . (mu_p    + b_p_s)
#   logit gamma_i    = X_gamma_i . beta_gamma                   (colonization)
#   logit eps_i      = X_eps_i   . beta_eps                     (extinction)
#   b_psi1_s ~ N(0, Sigma_psi1), b_p_s ~ N(0, Sigma_p)          (community RE)
#
# The first-season state and the season-to-season transitions integrate out the
# latent occupancy path z exactly by an HMM forward filter (the same recursion
# as the single-species dynamic model in build_dynamic_callbacks); the
# per-species first-season occupancy / detection deviations b_s = (b_psi1_s,
# b_p_s) are the random effects. The colonization (gamma) and extinction (eps)
# coefficients carry no per-species random effect: they are shared community
# globals, weakly identified at small per-season transition counts, so they get
# a weak Gaussian prior to stabilize them.
#
# Detection, gamma, and eps are SITE-LEVEL (constant across visits and seasons),
# matching the single-species dynamic occupancy model.
#
# Fit: the shared community Laplace-EM engine (.tobs_community_em), with the RE
# arms (psi1, p) in `theta` and the transition arms (gamma, eps) in `global`.
# The per-species marginal log-likelihood is the HMM forward filter summed over
# sites; the engine drives the community-covariance M-step and the marginal
# fixed-effect information. Non-spatial Laplace only.
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind a community dynamic occupancy model. `y` is a 4D array
# [n_sites x max_visits x n_seasons x n_species] or a named list of n_species
# 3D arrays [n_sites x max_visits x n_seasons]. The first-season occupancy
# design X_psi1, the detection design X_p, the colonization design X_gamma, and
# the extinction design X_eps are all site-level (one row per site). Detection
# is 0/1/NA; values that are NA or negative are missing visits.
.tobs_build_ms_dyn_occu <- function(occ_formula, det_formula,
                                    col_formula = NULL, ext_formula = NULL,
                                    data, y, species = NULL,
                                    structured_terms = list()) {
  to_array <- function(z) {
    if (is.list(z) && !is.array(z)) {
      n_sp <- length(z)
      d <- dim(z[[1L]])
      if (length(d) != 3L) {
        stop("each element of the y list must be a 3D array ",
             "[n_sites x max_visits x n_seasons].", call. = FALSE)
      }
      arr <- array(NA_real_, dim = c(d, n_sp))
      for (s in seq_len(n_sp)) arr[, , , s] <- as.array(z[[s]])
      attr(arr, "names_from") <- names(z)
      return(arr)
    }
    if (length(dim(z)) != 4L) {
      stop("y must be a 4D array ",
           "[n_sites x max_visits x n_seasons x n_species] ",
           "or a list of 3D arrays.", call. = FALSE)
    }
    z
  }
  y <- to_array(y)

  n_sites    <- dim(y)[1L]
  max_visits <- dim(y)[2L]
  n_seasons  <- dim(y)[3L]
  n_species  <- dim(y)[4L]

  species_names <- if (is.character(species)) species
                   else if (!is.null(attr(y, "names_from"))) attr(y, "names_from")
                   else paste0("sp", seq_len(n_species))
  if (length(species_names) != n_species) {
    stop(sprintf("species has %d names but y has %d species",
                 length(species_names), n_species), call. = FALSE)
  }
  .tobs_check_site_count(n_sites, nrow(data), "sites")

  if (is.null(col_formula)) col_formula <- ~ 1
  if (is.null(ext_formula)) ext_formula <- ~ 1

  X_psi1  <- stats::model.matrix(occ_formula, data)
  X_p     <- stats::model.matrix(det_formula, data)
  X_gamma <- stats::model.matrix(col_formula, data)
  X_eps   <- stats::model.matrix(ext_formula, data)

  # Per-species integer detection array and per-cell validity mask (a visit is
  # missing when it is NA or negative). Detections must be 0/1.
  y_int <- array(0L,    dim = dim(y))
  valid <- array(FALSE, dim = dim(y))
  for (s in seq_len(n_species)) {
    ys <- array(as.integer(round(y[, , , s])),
                dim = c(n_sites, max_visits, n_seasons))
    vs <- !is.na(ys) & ys >= 0L
    if (any(ys[vs] != 0L & ys[vs] != 1L)) {
      stop(sprintf("species '%s': y must contain only 0, 1, NA, or negative ",
                   species_names[s]), "(missing) values.", call. = FALSE)
    }
    ys[!vs] <- 0L
    y_int[, , , s] <- ys
    valid[, , , s] <- vs
  }

  structure(list(
    model_type    = "ms_dyn_occu",
    y             = y_int,
    valid         = valid,
    n_sites       = n_sites,
    max_visits    = max_visits,
    n_seasons     = n_seasons,
    n_species     = n_species,
    species_names = species_names,
    X_psi1        = X_psi1,
    X_p           = X_p,
    X_gamma       = X_gamma,
    X_eps         = X_eps,
    X_occ         = X_psi1,      # alias for the shared field-setup helpers
    formulas      = list(occ = occ_formula, det = det_formula,
                         col = col_formula, ext = ext_formula),
    structured_terms = structured_terms,
    data          = data,
    process_info  = list(
      list(name = "psi1",  p = ncol(X_psi1),
           coef_names = colnames(X_psi1),  link = "logit"),
      list(name = "p",     p = ncol(X_p),
           coef_names = colnames(X_p),     link = "logit"),
      list(name = "gamma", p = ncol(X_gamma),
           coef_names = colnames(X_gamma), link = "logit"),
      list(name = "eps",   p = ncol(X_eps),
           coef_names = colnames(X_eps),   link = "logit")
    )
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Laplace-EM fitter (shared community engine)
# ---------------------------------------------------------------------------

# Fit the community dynamic occupancy model via the shared community Laplace-EM
# engine. The RE-bearing coefficients (psi1, p) live in `theta`; the shared
# transition coefficients (gamma, eps) live in `global`. Returns a `tobs_fit`
# (via build_ms_dyn_occu_fit).
.tobs_fit_ms_dyn_occu <- function(model,
                                  priors     = NULL,
                                  max.iter   = 200L,
                                  tol        = 1e-4,
                                  sigma.beta = 5,
                                  verbose    = TRUE,
                                  ...) {
  dots       <- list(...)
  newton.max <- as.integer(dots$newton.max %||% 30L)

  pi_list <- model$process_info
  P_psi1 <- pi_list[[1L]]$p
  P_p    <- pi_list[[2L]]$p
  P_gam  <- pi_list[[3L]]$p
  P_eps  <- pi_list[[4L]]$p
  P      <- P_psi1 + P_p
  G      <- P_gam + P_eps
  S      <- model$n_species

  psi1_idx <- seq_len(P_psi1)
  p_idx    <- P_psi1 + seq_len(P_p)
  arm_idx  <- list(psi1 = psi1_idx, p = p_idx)

  gam_idx <- seq_len(P_gam)
  eps_idx <- P_gam + seq_len(P_eps)

  X_psi1  <- model$X_psi1
  X_p     <- model$X_p
  X_gamma <- model$X_gamma
  X_eps   <- model$X_eps
  n_sites   <- model$n_sites
  n_seasons <- model$n_seasons

  # Per-species detection / validity arrays.
  ys_list <- lapply(seq_len(S), function(s) model$y[, , , s])
  vs_list <- lapply(seq_len(S), function(s) model$valid[, , , s])

  # Per-(site, season) detection sufficient statistics, once per species: they
  # carry no parameter, so the forward marginal's only per-step input is the
  # emission built from them at the current p.
  em_stats <- lapply(seq_len(S), function(s)
    .ms_dyn_occu_emit_stats(ys_list[[s]], vs_list[[s]], n_sites, n_seasons))

  sp_ll <- function(s, theta, global) {
    beta_psi1 <- theta[psi1_idx]
    beta_p    <- theta[p_idx]
    beta_gam  <- global[gam_idx]
    beta_eps  <- global[eps_idx]
    psi1  <- stats::plogis(as.numeric(X_psi1  %*% beta_psi1))
    p     <- stats::plogis(as.numeric(X_p     %*% beta_p))
    gamma <- stats::plogis(as.numeric(X_gamma %*% beta_gam))
    eps   <- stats::plogis(as.numeric(X_eps   %*% beta_eps))
    em    <- .ms_dyn_occu_emissions(p, em_stats[[s]]$nvalid, em_stats[[s]]$ndet)
    .ms_dyn_occu_fwd_ll_vec(psi1, gamma, eps, em, n_sites, n_seasons)
  }

  # ---- warm start ----
  clamp01 <- function(q) min(max(q, 1e-3), 1 - 1e-3)
  # Naive season-1 occupancy proportion: fraction of (species, site) with any
  # detection in season 1; naive detection rate among detected visits.
  occ_props <- numeric(S); det_rates <- numeric(S)
  for (s in seq_len(S)) {
    v <- vs_list[[s]]; yy <- ys_list[[s]]
    v1 <- v[, , 1L, drop = FALSE]; y1 <- yy[, , 1L, drop = FALSE]
    site_det1 <- vapply(seq_len(n_sites), function(i) {
      any(y1[i, , 1L][v1[i, , 1L]] == 1L)
    }, logical(1))
    occ_props[s] <- mean(site_det1)
    detected <- yy[v]
    det_rates[s] <- if (length(detected)) mean(detected == 1L) else NA_real_
  }
  init_mu <- numeric(P)
  init_mu[psi1_idx][1L] <- stats::qlogis(clamp01(mean(occ_props)))
  dr <- mean(det_rates[is.finite(det_rates)])
  if (!is.finite(dr)) dr <- 0.3
  init_mu[p_idx][1L] <- stats::qlogis(clamp01(dr))

  init_global <- numeric(G)
  init_global[gam_idx][1L] <- stats::qlogis(0.15)
  init_global[eps_idx][1L] <- stats::qlogis(0.10)

  res <- .tobs_community_em(
    S = S, P = P, arm_idx = arm_idx,
    sp_ll = sp_ll, sp_grad = NULL,
    init_mu = init_mu, init_global = init_global,
    penalize_global = TRUE, sigma_beta = sigma.beta, priors = priors,
    sigma_init = 0.3, max_iter = as.integer(max.iter), tol = as.numeric(tol),
    newton_max = newton.max, verbose = isTRUE(verbose)
  )

  build_ms_dyn_occu_fit(model, res, arm_idx, gam_idx, eps_idx)
}


# ---------------------------------------------------------------------------
# Wrap the EM output into a tobs_fit
# ---------------------------------------------------------------------------

build_ms_dyn_occu_fit <- function(model, res, arm_idx, gam_idx, eps_idx) {
  pi_list <- model$process_info
  P_psi1 <- pi_list[[1L]]$p
  P_p    <- pi_list[[2L]]$p
  P_gam  <- pi_list[[3L]]$p
  P_eps  <- pi_list[[4L]]$p
  P      <- P_psi1 + P_p

  mu     <- res$mu
  global <- res$global

  beta_names <- c(
    paste0("psi1_",  pi_list[[1L]]$coef_names),
    paste0("p_",     pi_list[[2L]]$coef_names),
    paste0("gamma_", pi_list[[3L]]$coef_names),
    paste0("eps_",   pi_list[[4L]]$coef_names)
  )
  par_names <- beta_names

  means <- c(mu, global); names(means) <- par_names
  V <- res$Vf; dimnames(V) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  # Per-species community structure (mu + BLUP deviations) for the RE arms.
  B <- do.call(rbind, res$b_list)              # S x P
  arm_block <- function(arm) {
    idx  <- arm_idx[[arm]]
    blup <- B[, idx, drop = FALSE]
    coef <- sweep(blup, 2L, mu[idx], "+")
    rownames(blup) <- rownames(coef) <- model$species_names
    list(blup = blup, coef = coef)
  }
  psi1_b <- arm_block("psi1"); p_b <- arm_block("p")
  colnames(psi1_b$blup) <- colnames(psi1_b$coef) <- pi_list[[1L]]$coef_names
  colnames(p_b$blup)    <- colnames(p_b$coef)    <- pi_list[[2L]]$coef_names

  Sigma_psi1 <- res$Sigma$psi1; Sigma_p <- res$Sigma$p
  dimnames(Sigma_psi1) <- list(pi_list[[1L]]$coef_names, pi_list[[1L]]$coef_names)
  dimnames(Sigma_p)    <- list(pi_list[[2L]]$coef_names, pi_list[[2L]]$coef_names)

  F_val <- res$logML

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(F_val, n_draws),
    log_lik      = F_val,
    N            = sum(model$valid)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    n_fixed      = length(means),
    fixed_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = NULL,
    method       = "laplace",
    ms_community = list(
      Sigma_psi1 = Sigma_psi1, Sigma_p = Sigma_p,
      sd_psi1 = sqrt(pmax(diag(Sigma_psi1), 0)),
      sd_p    = sqrt(pmax(diag(Sigma_p),    0)),
      coef_psi1 = psi1_b$coef, coef_p = p_b$coef,
      blup_psi1 = psi1_b$blup, blup_p = p_b$blup,
      # Per-species posterior covariance Cov(b_s|y) (Louis 1982, from the
      # community EM's own Newton solve, conditional on the converged
      # community mean) -- what a per-species-coefficient consumer (SBC's
      # "rank a fixed species set" design, a calibrated per-species CI) needs
      # beyond the point BLUP; not previously exposed on the fit object. Bf =
      # the (mu,global)-b_s cross-Hessian block from the same Newton solve:
      # mu/global and b_s are NOT independent in the posterior, and Bf is
      # what lets a consumer draw them jointly instead -- see
      # .tobs_sbc_community_b_draws (R/sbc.R).
      Cinv = res$Cinv, Bf = res$Bf
    ),
    convergence  = list(converged = isTRUE(res$converged), n_iter = res$n_iter)
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "ms_dyn_occu")
# ---------------------------------------------------------------------------

# Per-species BLUP deviations, long form: one row per (species, arm, term) over
# the RE arms (psi1, p). The shared gamma / eps coefficients carry no per-
# species random effect and so do not appear here.
.tobs_ranef_ms_dyn_occu <- function(object) {
  .tobs_ranef_ms_long(object$ms_community,
                      c(psi1 = "blup_psi1", p = "blup_p"))
}

# Per-species posterior-mean linear predictors: site-level first-season
# occupancy psi1 [n_sites x n_species] and site-level detection p
# [n_sites x n_species]; plus the community colonization gamma and extinction
# eps as length-n_sites vectors (no species dimension, since they are shared).
.tobs_fitted_ms_dyn_occu <- function(object) {
  model <- object$model
  cm    <- object$ms_community
  X_psi1  <- model$X_psi1
  X_p     <- model$X_p
  X_gamma <- model$X_gamma
  X_eps   <- model$X_eps
  pi_list <- model$process_info
  P_gam   <- pi_list[[3L]]$p
  P_eps   <- pi_list[[4L]]$p

  beta_gam <- object$means[paste0("gamma_", pi_list[[3L]]$coef_names)]
  beta_eps <- object$means[paste0("eps_",   pi_list[[4L]]$coef_names)]

  psi1 <- stats::plogis(X_psi1 %*% t(cm$coef_psi1))
  p    <- stats::plogis(X_p    %*% t(cm$coef_p))
  gamma <- as.numeric(stats::plogis(X_gamma %*% beta_gam))
  eps   <- as.numeric(stats::plogis(X_eps   %*% beta_eps))
  dimnames(psi1) <- dimnames(p) <- list(NULL, model$species_names)
  list(psi1 = psi1, p = p, gamma = gamma, eps = eps)
}

# Draw community dynamic occupancy data under the fitted per-species first-season
# / detection coefficients and the shared transition coefficients, at the
# observed season / visit pattern. Returns a 4D array matching the input y.
.tobs_simulate_ms_dyn_occu <- function(object, nsim = 1) {
  model <- object$model
  cm    <- object$ms_community
  n_sites    <- model$n_sites
  max_visits <- model$max_visits
  n_seasons  <- model$n_seasons
  n_species  <- model$n_species
  pi_list <- model$process_info

  beta_gam <- object$means[paste0("gamma_", pi_list[[3L]]$coef_names)]
  beta_eps <- object$means[paste0("eps_",   pi_list[[4L]]$coef_names)]
  gamma <- as.numeric(stats::plogis(model$X_gamma %*% beta_gam))
  eps   <- as.numeric(stats::plogis(model$X_eps   %*% beta_eps))

  # Per-species psi1 / p (community means) + per-site gamma / eps (above); the
  # season-1 state, transitions, and detections run in cpp_simulate_ms_dyn_occu
  # from R's RNG stream in the former order (byte-identical).
  psi1 <- vapply(seq_len(n_species),
                 function(s) as.numeric(stats::plogis(model$X_psi1 %*% cm$coef_psi1[s, ])),
                 numeric(n_sites))
  p <- vapply(seq_len(n_species),
              function(s) as.numeric(stats::plogis(model$X_p %*% cm$coef_p[s, ])),
              numeric(n_sites))
  res <- cpp_simulate_ms_dyn_occu(psi1, p, gamma, eps, as.integer(model$valid),
    n_sites, max_visits, n_seasons, n_species, as.integer(nsim))
  dn <- list(NULL, NULL, NULL, model$species_names)
  if (nsim == 1L) { a <- array(res[, , , , 1], dim = dim(res)[1:4]); dimnames(a) <- dn; return(a) }
  lapply(seq_len(nsim), function(s) { a <- array(res[, , , , s], dim = dim(res)[1:4]); dimnames(a) <- dn; a })
}


# ---------------------------------------------------------------------------
# Family constructor
# ---------------------------------------------------------------------------

#' Community (multispecies) dynamic occupancy family
#'
#' Per-species dynamic (HMM) occupancy with Gaussian community hyperpriors on
#' the per-species first-season occupancy and detection coefficients, and shared
#' community-wide colonisation / extinction transition coefficients.
#'
#' @section Scope:
#' The Laplace engine is the supported route: the shared colonisation /
#' extinction dynamics and the per-species first-season occupancy / detection
#' components recover across seeds (see `tests/testthat/test-ms-dyn-occu.R`,
#' community-mean 95% CI coverage measured ~0.98). A shared areal field on the
#' first-season occupancy formula (`~ 1 + icar(graph = adj)`) fits the
#' `spOccupancy` `stMsPGOcc` model under `method = "nested_laplace"`: the field
#' is shared across species, and because the first-season occupancy `psi1` only
#' sets the initial mixing weight of each species' HMM, the block-coordinate
#' driver alternates the community EM (field as a `psi1` offset) with an areal
#' field Newton, the field recovering cleanly (`cor` ~0.94). A
#' spatially-varying-coefficient bar (`~ spatial(~ 1 + w || cell, graph = adj)`)
#' adds a shared covariate-weighted field alongside the intercept field, fitting
#' the `svcTMsPGOcc` model through the same K-field weighted-ICAR solve as the
#' community count SVC (both fields recover, `cor` ~0.90 / ~0.89). `icar()` only;
#' a NUTS sampler and `bym2()` / `car_proper()` fields are a deliberate
#' follow-up; `method = "nuts"` errors from the dispatcher with a pointer rather
#' than silently downgrading.
#'
#' @return A `tobs_family` object.
#' @seealso [dyn_occu()], [ms_occu()]
#' @export
ms_dyn_occu <- function() {
  obs_family(
    name           = "ms_dyn_occu",
    class_long     = "community dynamic occupancy",
    latent         = "bernoulli_hmm",
    observation    = "binomial_detection",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    # The shared psi1 field is fit by the block-coordinate driver, but the call
    # passes `latent = NULL`: a field block and no factors, so no candidate
    # starting directions and no `factor.starts`.
    control_groups = "block_coordinate"
  )
}
