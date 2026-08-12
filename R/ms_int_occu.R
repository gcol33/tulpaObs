# =============================================================================
# ms_int_occu.R - community / multispecies integrated occupancy
#
# The community version of int_occu(): multiple detection data sources share one
# latent occupancy state per species, with Gaussian community hyperpriors on the
# per-species coefficients of the occupancy arm and of every per-source
# detection arm.
#
#   z_{s,i}          ~ Bernoulli(psi_{s,i})                    (latent presence)
#   y_{s,i,d,j}|z=1  ~ Bernoulli(p_{s,i,d})                    (source-d detection)
#   logit psi_{s,i}      = X_psi_i  . (mu_psi  + b_psi_s)
#   logit p_{s,i,d}      = X_p_d_i  . (mu_p_d  + b_p_d_s)
#   b_psi_s ~ N(0, Sigma_psi), b_p_d_s ~ N(0, Sigma_p_d)       (community RE)
#
# The latent presence z marginalises out per species-site in closed form (a
# two-state mixture, the multi-source generalisation of single-season
# occupancy); the per-species coefficient deviations b_s = (b_psi_s, b_p1_s,
# ..., b_pD_s) are the random effects, fit by the shared community Laplace-EM
# engine (.tobs_community_em in R/community_em.R). There are no global (shared,
# RE-free) coefficients here, so G = 0.
#
# Detection designs are site-level per source (one detection probability per
# source-site, constant across that source's visits), matching the single-
# species integrated model in R/laplace.R::build_integrated_callbacks.
#
# Site coverage: each source covers a (possibly partial, possibly overlapping)
# subset of the n_sites in `data`, declared by a per-source `site_map` (the
# global site index of each of that source's rows; full overlap if omitted, the
# back-compatible default). Sources are scattered into a full n_sites detection
# array padded with NA (no visit) at the sites a source does not cover; an
# uncovered (site, source) cell has n_valid = 0 and so drops out of the marginal
# (zero contribution), which is exactly the partial/overlapping coverage model --
# no special-casing in the per-species likelihood. Non-spatial Laplace only.
# =============================================================================


# ---------------------------------------------------------------------------
# Per-source detection summaries
# ---------------------------------------------------------------------------

# Precompute, per source and site, the number of non-missing visits and the
# number of detections for one species, plus the per-site any-detection flag
# across all sources. `y_list` is a list of D integer matrices [n_sites x J_d]
# with 0/1 entries and a parallel `valid_list` of logical matrices.
.ms_int_occu_sp_summary <- function(y_list, valid_list) {
  D <- length(y_list)
  n_sites <- nrow(y_list[[1L]])
  n_valid <- matrix(0L, n_sites, D)
  n_det   <- matrix(0L, n_sites, D)
  for (d in seq_len(D)) {
    vs <- valid_list[[d]]
    yd <- y_list[[d]]
    n_valid[, d] <- rowSums(vs)
    n_det[, d]   <- rowSums(yd * vs)
  }
  any_det <- rowSums(n_det) > 0L
  list(n_valid = n_valid, n_det = n_det, any_det = any_det)
}


# ---------------------------------------------------------------------------
# Per-species marginal log-likelihood and gradient
# ---------------------------------------------------------------------------

# Per-species marginal log-likelihood. `eta_psi` is the per-site occupancy
# linear predictor (length n_sites); `eta_p` is a list of D per-site detection
# linear predictors (each length n_sites). `summ` is the .ms_int_occu_sp_summary
# for this species. The latent z is integrated out exactly (det / no-det branch
# per site), accumulated in log space with a log-sum-exp on the no-detection
# branch for stability. `per_site = TRUE` returns the per-site vector instead of
# its sum, which is what the community latent driver's joint site marginal reads.
.ms_int_occu_sp_ll <- function(eta_psi, eta_p, summ, per_site = FALSE) {
  clp <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)
  psi <- clp(stats::plogis(eta_psi))
  D <- length(eta_p)
  n_valid <- summ$n_valid; n_det <- summ$n_det; any_det <- summ$any_det

  # sum_d [ n_det * log p_d + (n_valid - n_det) * log(1 - p_d) ]
  log_det_term  <- numeric(length(psi))
  # sum_d n_valid * log(1 - p_d)  (the occupied-undetected log mass)
  log_undet_mass <- numeric(length(psi))
  for (d in seq_len(D)) {
    pd  <- clp(stats::plogis(eta_p[[d]]))
    lpd <- log(pd); l1m <- log1p(-pd)
    log_det_term   <- log_det_term + n_det[, d] * lpd +
                      (n_valid[, d] - n_det[, d]) * l1m
    log_undet_mass <- log_undet_mass + n_valid[, d] * l1m
  }

  ll <- numeric(length(psi))
  # Detection sites: z = 1 forced, LL = log psi + detection log-likelihood.
  if (any(any_det)) {
    ll[any_det] <- log(psi[any_det]) + log_det_term[any_det]
  }
  # No-detection sites: marginal = psi * prod(1 - p)^n_valid + (1 - psi),
  # via log-sum-exp of the two log masses.
  nd <- !any_det
  if (any(nd)) {
    a <- log(psi[nd]) + log_undet_mass[nd]    # occupied-undetected log mass
    b <- log1p(-psi[nd])                      # unoccupied log mass
    m <- pmax(a, b)
    ll[nd] <- m + log(exp(a - m) + exp(b - m))
  }
  if (per_site) ll else sum(ll)
}

# Per-species packed coefficient gradient over theta = (beta_psi, beta_p1, ...,
# beta_pD). `X_psi` [n_sites x P_psi]; `X_p` a list of D site-level designs.
# The multi-source generalisation of the single-season occupancy marginal
# gradient, chained from the eta-gradient to the coefficients.
.ms_int_occu_sp_grad <- function(eta_psi, eta_p, summ, X_psi, X_p) {
  clp <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)
  psi <- clp(stats::plogis(eta_psi))
  D <- length(eta_p)
  n_valid <- summ$n_valid; n_det <- summ$n_det; any_det <- summ$any_det
  n_sites <- length(psi)

  p_list <- vector("list", D)
  log_prod <- numeric(n_sites)                # sum_d n_valid * log(1 - p_d)
  for (d in seq_len(D)) {
    pd <- clp(stats::plogis(eta_p[[d]]))
    p_list[[d]] <- pd
    log_prod <- log_prod + n_valid[, d] * log1p(-pd)
  }
  prodterm <- exp(log_prod)                    # prod_d (1 - p_d)^n_valid
  A <- psi * prodterm                          # occupied-undetected mass
  B <- 1 - psi                                 # unoccupied mass
  L <- A + B

  # ---- occupancy predictor ----
  d_eta_psi <- numeric(n_sites)
  d_eta_psi[any_det] <- 1 - psi[any_det]
  nd <- !any_det
  d_eta_psi[nd] <- psi[nd] * (1 - psi[nd]) * (prodterm[nd] - 1) / L[nd]

  # ---- per-source detection predictors ----
  g_psi <- as.numeric(crossprod(X_psi, d_eta_psi))
  g_p   <- vector("list", D)
  for (d in seq_len(D)) {
    pd <- p_list[[d]]
    d_eta_pd <- numeric(n_sites)
    d_eta_pd[any_det] <- n_det[any_det, d] - n_valid[any_det, d] * pd[any_det]
    d_eta_pd[nd] <- -(A[nd] / L[nd]) * n_valid[nd, d] * pd[nd]
    g_p[[d]] <- as.numeric(crossprod(X_p[[d]], d_eta_pd))
  }
  c(g_psi, unlist(g_p))
}


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind a community integrated occupancy model. `y` is a list of length D
# (n_sources); element d is a 3D array [n_sites x J_d x n_species] or a named
# list of n_species matrices [n_sites x J_d]. The occupancy design X_psi is
# site-level; each source's detection design X_p_d is site-level (one detection
# probability per source-site). `det_formula` is one formula applied to every
# source (separate coefficient vector per source) or a list of D formulas.
.tobs_build_ms_int_occu <- function(occ_formula, det_formula, data, y, species,
                                    site_map = NULL) {
  if (!is.list(y) || is.data.frame(y)) {
    stop("y must be a list of length n_sources, each a 3D array ",
         "[n_sites x J_d x n_species] or a list of per-species matrices.",
         call. = FALSE)
  }
  D <- length(y)
  if (D < 1L) stop("y must contain at least one data source.", call. = FALSE)

  to_array <- function(z, label) {
    if (is.list(z) && !is.array(z)) {
      n_sp <- length(z)
      arr <- array(NA_real_, dim = c(nrow(z[[1L]]), ncol(z[[1L]]), n_sp))
      for (s in seq_len(n_sp)) arr[, , s] <- as.matrix(z[[s]])
      attr(arr, "names_from") <- names(z)
      return(arr)
    }
    if (length(dim(z)) != 3L) {
      stop(sprintf("%s must be a 3D array [n_sites x J_d x n_species] ", label),
           "or a list of matrices.", call. = FALSE)
    }
    z
  }
  y_arr <- lapply(seq_len(D), function(d) to_array(y[[d]], sprintf("y[[%d]]", d)))

  n_sites   <- nrow(data)
  n_species <- dim(y_arr[[1L]])[3L]
  J_d       <- vapply(y_arr, function(a) dim(a)[2L], integer(1))

  # Resolve each source's site coverage. `site_map` (a list of D integer vectors,
  # one global site index per source row) supports partial / overlapping coverage;
  # omitted, each source must span all n_sites in declaration order (full overlap,
  # the historical default).
  if (!is.null(site_map)) {
    if (!is.list(site_map) || length(site_map) != D) {
      stop(sprintf("site_map must be a list of %d integer vectors (one per source).",
                   D), call. = FALSE)
    }
  }
  maps <- vector("list", D)
  for (d in seq_len(D)) {
    n_d <- dim(y_arr[[d]])[1L]
    if (dim(y_arr[[d]])[3L] != n_species) {
      stop(sprintf("source %d has %d species but source 1 has %d.",
                   d, dim(y_arr[[d]])[3L], n_species), call. = FALSE)
    }
    if (is.null(site_map)) {
      if (n_d != n_sites) {
        stop(sprintf("source %d has %d sites but data has %d rows. Supply ",
                     d, n_d, n_sites),
             "`site_map` (one global site index per source row) for partial / ",
             "overlapping coverage.", call. = FALSE)
      }
      maps[[d]] <- seq_len(n_sites)
    } else {
      m <- as.integer(site_map[[d]])
      if (length(m) != n_d) {
        stop(sprintf("site_map[[%d]] has %d entries but source %d has %d rows.",
                     d, length(m), d, n_d), call. = FALSE)
      }
      if (anyNA(m) || any(m < 1L) || any(m > n_sites)) {
        stop(sprintf("site_map[[%d]] must index sites in 1..%d.", d, n_sites),
             call. = FALSE)
      }
      if (anyDuplicated(m)) {
        stop(sprintf("site_map[[%d]] has duplicate site indices (a source maps ",
                     d), "each of its rows to a distinct site).", call. = FALSE)
      }
      maps[[d]] <- m
    }
  }

  src_names <- if (!is.null(names(y)) && all(nzchar(names(y)))) names(y)
               else paste0("src", seq_len(D))
  proc_names <- paste0("p", seq_len(D))

  species_names <-
    if (is.character(species)) species
    else if (!is.null(attr(y_arr[[1L]], "names_from"))) attr(y_arr[[1L]], "names_from")
    else paste0("sp", seq_len(n_species))
  if (length(species_names) != n_species) {
    stop(sprintf("species has %d names but y has %d species",
                 length(species_names), n_species), call. = FALSE)
  }

  det_formulas <- if (inherits(det_formula, "formula")) {
    rep(list(det_formula), D)
  } else if (is.list(det_formula)) {
    if (length(det_formula) != D) {
      stop(sprintf("det_formula list has %d entries but y has %d sources.",
                   length(det_formula), D), call. = FALSE)
    }
    det_formula
  } else {
    stop("det_formula must be a formula or a list of D formulas.", call. = FALSE)
  }

  X_psi <- stats::model.matrix(occ_formula, data)
  X_p   <- lapply(det_formulas, function(f) stats::model.matrix(f, data))

  # Clean each source / species detection matrix: binary, NA -> 0 with a
  # parallel validity mask. Each source's rows are scattered into the full
  # n_sites grid at its `maps[[d]]` global sites; sites the source does not cover
  # stay all-invalid (no visit), so they contribute nothing to that source's
  # marginal terms.
  y_int <- vector("list", D)
  valid <- vector("list", D)
  for (d in seq_len(D)) {
    md  <- maps[[d]]
    n_d <- length(md)
    yi <- array(0L,    dim = c(n_sites, J_d[d], n_species))
    vi <- array(FALSE, dim = c(n_sites, J_d[d], n_species))
    for (s in seq_len(n_species)) {
      ys <- matrix(as.integer(round(y_arr[[d]][, , s])), n_d, J_d[d])
      vs <- !is.na(ys)
      if (any(ys[vs] != 0L & ys[vs] != 1L)) {
        stop(sprintf("source %d, species '%s': y must contain only 0, 1, or NA.",
                     d, species_names[s]), call. = FALSE)
      }
      ys[!vs] <- 0L
      yi[md, , s] <- ys
      vi[md, , s] <- vs
    }
    y_int[[d]] <- yi
    valid[[d]] <- vi
  }

  # Per-species precomputed detection summaries (shared by ll / grad / fit).
  summaries <- lapply(seq_len(n_species), function(s) {
    y_s <- lapply(seq_len(D), function(d) y_int[[d]][, , s])
    v_s <- lapply(seq_len(D), function(d) valid[[d]][, , s])
    .ms_int_occu_sp_summary(y_s, v_s)
  })

  process_info <- c(
    list(list(name = "psi", p = ncol(X_psi),
              coef_names = colnames(X_psi), link = "logit")),
    lapply(seq_len(D), function(d) {
      list(name = proc_names[d], p = ncol(X_p[[d]]),
           coef_names = colnames(X_p[[d]]), link = "logit")
    })
  )

  structure(list(
    model_type    = "ms_int_occu",
    y             = y_int,
    valid         = valid,
    n_sites       = n_sites,
    n_sources     = D,
    n_visits      = J_d,
    site_maps     = maps,
    n_species     = n_species,
    species_names = species_names,
    source_names  = src_names,
    process_names = proc_names,
    X_psi         = X_psi,
    X_p           = X_p,
    summaries     = summaries,
    formulas      = list(occ = occ_formula, det = det_formulas),
    data          = data,
    process_info  = process_info
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Laplace-EM fitter
# ---------------------------------------------------------------------------

# Fit the community integrated occupancy model via the shared community
# Laplace-EM engine. `model` is the bound ms_int_occu model. Returns a
# `tobs_fit` (via build_ms_int_occu_fit).
.tobs_fit_ms_int_occu <- function(model,
                                  priors     = NULL,
                                  max.iter   = 200L,
                                  tol        = 1e-4,
                                  sigma.beta = 5,
                                  verbose    = TRUE,
                                  ...) {
  dots <- list(...)
  newton.max <- as.integer(dots$newton.max %||% 30L)

  pi_list <- model$process_info
  D       <- model$n_sources
  P_psi   <- pi_list[[1L]]$p
  P_p     <- vapply(seq_len(D), function(d) pi_list[[d + 1L]]$p, integer(1))
  P       <- P_psi + sum(P_p)
  S       <- model$n_species

  # arm_idx over theta: psi coefs, then each source's detection coefs.
  psi_idx <- seq_len(P_psi)
  p_idx   <- vector("list", D)
  off <- P_psi
  for (d in seq_len(D)) {
    p_idx[[d]] <- off + seq_len(P_p[d])
    off <- off + P_p[d]
  }
  arm_idx <- c(list(psi = psi_idx),
               stats::setNames(p_idx, model$process_names))

  X_psi <- model$X_psi
  X_p   <- model$X_p
  summaries <- model$summaries

  eta_from_theta <- function(theta) {
    eta_psi <- as.numeric(X_psi %*% theta[psi_idx])
    eta_p   <- lapply(seq_len(D), function(d) as.numeric(X_p[[d]] %*% theta[p_idx[[d]]]))
    list(psi = eta_psi, p = eta_p)
  }
  sp_ll <- function(s, theta, global) {
    e <- eta_from_theta(theta)
    .ms_int_occu_sp_ll(e$psi, e$p, summaries[[s]])
  }
  sp_grad <- function(s, theta, global) {
    e <- eta_from_theta(theta)
    .ms_int_occu_sp_grad(e$psi, e$p, summaries[[s]], X_psi, X_p)
  }

  # ---- warm start ----
  clp <- function(x) min(max(x, 1e-3), 1 - 1e-3)
  any_det_prop <- mean(vapply(summaries, function(z) mean(z$any_det), numeric(1)))
  mu <- numeric(P)
  mu[psi_idx][1L] <- stats::qlogis(clp(any_det_prop))
  for (d in seq_len(D)) {
    # Naive per-source detection rate among detected sites.
    det_sites <- vapply(summaries, function(z) sum(z$n_det[, d]), numeric(1))
    val_sites <- vapply(summaries, function(z) sum(z$n_valid[z$any_det, d]),
                        numeric(1))
    rate <- if (sum(val_sites) > 0) sum(det_sites) / sum(val_sites) else 0.5
    mu[p_idx[[d]]][1L] <- stats::qlogis(clp(rate))
  }

  fit <- .tobs_community_em(
    S = S, P = P, arm_idx = arm_idx,
    sp_ll = sp_ll, sp_grad = sp_grad,
    init_mu = mu, init_global = numeric(0),
    penalize_global = FALSE, sigma_beta = sigma.beta, priors = priors,
    sigma_init = 0.3, max_iter = as.integer(max.iter), tol = as.numeric(tol),
    newton_max = newton.max, verbose = isTRUE(verbose)
  )

  build_ms_int_occu_fit(model, fit, arm_idx)
}


# ---------------------------------------------------------------------------
# Wrap the EM output into a tobs_fit
# ---------------------------------------------------------------------------

build_ms_int_occu_fit <- function(model, fit, arm_idx) {
  pi_list   <- model$process_info
  D         <- model$n_sources
  arm_names <- names(arm_idx)              # psi, p1, ..., pD

  beta_names <- unlist(lapply(seq_along(arm_names), function(k) {
    paste0(arm_names[k], "_", pi_list[[k]]$coef_names)
  }))
  par_names <- beta_names

  means <- fit$mu; names(means) <- par_names
  V <- fit$Vf; dimnames(V) <- list(par_names, par_names)
  V <- (V + t(V)) / 2
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  # Per-species community structure (mu + BLUP deviations) per arm.
  B <- do.call(rbind, fit$b_list)          # S x P
  Sigma_list <- list(); sd_list <- list()
  coef_list  <- list(); blup_list <- list()
  for (k in seq_along(arm_names)) {
    arm <- arm_names[k]
    idx <- arm_idx[[arm]]
    cn  <- pi_list[[k]]$coef_names
    blup <- B[, idx, drop = FALSE]
    coef <- sweep(blup, 2L, means[idx], "+")
    rownames(blup) <- rownames(coef) <- model$species_names
    colnames(blup) <- colnames(coef) <- cn
    Sig <- fit$Sigma[[arm]]
    dimnames(Sig) <- list(cn, cn)
    Sigma_list[[paste0("Sigma_", arm)]] <- Sig
    sd_list[[paste0("sd_", arm)]]       <- sqrt(pmax(diag(Sig), 0))
    coef_list[[paste0("coef_", arm)]]   <- coef
    blup_list[[paste0("blup_", arm)]]   <- blup
  }
  # Per-species posterior covariance Cov(b_s|y) (Louis 1982, from the
  # community EM's own Newton solve, conditional on the converged community
  # mean) -- what a per-species-coefficient consumer (SBC's "rank a fixed
  # species set" design, a calibrated per-species CI) needs beyond the point
  # BLUP; not previously exposed on the fit object. Covers the FULL b_s
  # vector across every arm (psi + all D detection sources), matching
  # `B <- do.call(rbind, fit$b_list)` above. Bf = the mu-b_s cross-Hessian
  # block from the same Newton solve (gcol33/tulpaObs#226): mu and b_s are
  # NOT independent in the posterior, and Bf is what lets a consumer draw
  # them jointly instead.
  ms_community <- c(Sigma_list, sd_list, coef_list, blup_list,
                    list(Cinv = fit$Cinv, Bf = fit$Bf))

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(fit$logML, n_draws),
    log_lik      = fit$logML,
    N            = sum(vapply(model$valid, sum, integer(1)))),
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
    ms_community = ms_community,
    convergence  = list(converged = isTRUE(fit$converged), n_iter = fit$n_iter)
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "ms_int_occu")
# ---------------------------------------------------------------------------

# Per-species BLUP deviations, long form: one row per (species, arm, term),
# over the occupancy arm and every per-source detection arm.
.tobs_ranef_ms_int_occu <- function(object) {
  cm <- object$ms_community
  arms <- c("psi", object$model$process_names)
  to_long <- function(B, arm) {
    sp <- rownames(B); tm <- colnames(B)
    data.frame(species = rep(sp, times = ncol(B)), arm = arm,
               term = rep(tm, each = nrow(B)),
               estimate = as.numeric(B), stringsAsFactors = FALSE)
  }
  parts <- lapply(arms, function(arm) {
    to_long(cm[[paste0("blup_", arm)]], arm)
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

# Per-species posterior-mean linear predictors: site-level occupancy psi
# [n_sites x n_species] and a per-source list of site-level detection
# probabilities p[[d]] [n_sites x n_species].
.tobs_fitted_ms_int_occu <- function(object) {
  model <- object$model
  cm    <- object$ms_community
  D     <- model$n_sources

  psi <- stats::plogis(model$X_psi %*% t(cm$coef_psi))
  dimnames(psi) <- list(NULL, model$species_names)

  p_list <- vector("list", D)
  for (d in seq_len(D)) {
    arm <- model$process_names[d]
    pd  <- stats::plogis(model$X_p[[d]] %*% t(cm[[paste0("coef_", arm)]]))
    dimnames(pd) <- list(NULL, model$species_names)
    p_list[[d]] <- pd
  }
  names(p_list) <- model$process_names
  list(psi = psi, p = p_list)
}

# Draw community integrated data under the fitted per-species coefficients, at
# the observed per-source visit pattern (one z per species-site shared across
# sources, then per-source detection). Returns a list of per-source 3D arrays
# matching the input y.
.tobs_simulate_ms_int_occu <- function(object, nsim = 1) {
  model <- object$model
  cm    <- object$ms_community
  D     <- model$n_sources
  n_sites <- model$n_sites
  n_species <- model$n_species
  J_d <- model$n_visits

  # Per-species psi + per-source detection (community means, deterministic); the
  # z + per-source detections run in cpp_simulate_ms_int_occu from R's RNG stream
  # in the former order (byte-identical).
  psi <- vapply(seq_len(n_species),
                function(s) stats::plogis(as.numeric(model$X_psi %*% cm$coef_psi[s, ])),
                numeric(n_sites))
  pd_list <- lapply(seq_len(D), function(d) {
    arm <- model$process_names[d]
    vapply(seq_len(n_species),
           function(s) stats::plogis(as.numeric(model$X_p[[d]] %*% cm[[paste0("coef_", arm)]][s, ])),
           numeric(n_sites))
  })
  valid_list <- lapply(seq_len(D), function(d) as.integer(model$valid[[d]]))
  res <- cpp_simulate_ms_int_occu(psi, pd_list, valid_list, as.integer(J_d),
                                  n_sites, n_species, D, as.integer(nsim))
  add_names <- function(srcs) {
    names(srcs) <- model$process_names
    lapply(srcs, function(a) { dimnames(a) <- list(NULL, NULL, model$species_names); a })
  }
  res <- lapply(res, add_names)
  if (nsim == 1L) res[[1]] else res
}


# ---------------------------------------------------------------------------
# Family constructor
# ---------------------------------------------------------------------------

#' Community (multispecies) integrated occupancy family
#'
#' Multiple detection data sources share one latent occupancy state per species,
#' with Gaussian community hyperpriors on the per-species coefficients of the
#' occupancy arm and of every per-source detection arm. Fit by a shared
#' community Laplace-EM.
#'
#' @section Scope:
#' The Laplace engine is the supported route: the shared occupancy mean and the
#' per-source detection community components recover near nominal across seeds
#' and more than one source (see `tests/testthat/test-ms-int-occu.R`). A NUTS
#' sampler and an areal-field path are a deliberate follow-up, not part of this
#' family's working surface; `method = "nuts"` / `"nested_laplace"` error from the
#' dispatcher with a pointer rather than silently downgrading. The binary
#' community-mean intervals carry the mild Laplace under-dispersion typical of
#' occupancy data (measured 95% CI coverage ~0.89 at small per-species n).
#'
#' @return A `tobs_family` object.
#' @seealso [int_occu()], [ms_occu()]
#' @export
ms_int_occu <- function() {
  obs_family(
    name           = "ms_int_occu",
    class_long     = "community integrated occupancy",
    latent         = "bernoulli",
    observation    = "multisource_detection",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working"
  )
}
