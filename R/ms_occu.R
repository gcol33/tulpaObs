# =============================================================================
# ms_occu.R - community / multispecies single-season occupancy
#
# Per-species occupancy and detection coefficient RE with per-arm Gaussian
# community covariances (the spOccupancy msPGOcc model):
#
#   z_{s,i}        ~ Bernoulli(psi_{s,i})                      (latent presence)
#   y_{s,i,j}|z=1  ~ Bernoulli(p_{s,i})                        (detection)
#   logit psi_{s,i} = X_occ_i . (mu_occ + b_occ_s)
#   logit p_{s,i}   = X_det_i . (mu_p   + b_p_s)
#   b_occ_s ~ N(0, Sigma_occ), b_p_s ~ N(0, Sigma_p)          (community RE)
#
# The latent presence z marginalises out per species-site in closed form (the
# two-state mixture); the per-species coefficient deviations b_s = (b_occ_s,
# b_p_s) are the random effects, fit by the shared community Laplace-EM
# (.tobs_community_em in R/community_em.R). Occupancy and detection deviations
# are INDEPENDENT per arm, each with its own community covariance -- matching
# simulate_ms_occu() and spOccupancy::msPGOcc.
#
# Detection is site-level (one detection probability per species-site), as in
# simulate_ms_occu(). This is the single-source case of the integrated community
# marginal, so the per-species log-likelihood / gradient kernel is reused from
# R/ms_int_occu.R (single source of truth). Non-spatial Laplace only.
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind a community single-season occupancy model. `y` is a 3D array
# [n_sites x max_visits x n_species] or a named list of n_species matrices
# [n_sites x max_visits]. The occupancy design X_occ and detection design X_det
# are site-level. Detection is 0/1/NA.
.tobs_build_ms_occu <- function(occ_formula, det_formula, data, y, species,
                                structured_terms = list()) {
  to_array <- function(z) {
    if (is.list(z) && !is.array(z)) {
      n_sp <- length(z)
      arr <- array(NA_real_, dim = c(nrow(z[[1L]]), ncol(z[[1L]]), n_sp))
      for (s in seq_len(n_sp)) arr[, , s] <- as.matrix(z[[s]])
      attr(arr, "names_from") <- names(z)
      return(arr)
    }
    if (length(dim(z)) != 3L) {
      stop("y must be a 3D array [n_sites x max_visits x n_species] ",
           "or a list of matrices.", call. = FALSE)
    }
    z
  }
  y <- to_array(y)

  n_sites    <- dim(y)[1L]
  max_visits <- dim(y)[2L]
  n_species  <- dim(y)[3L]

  species_names <- if (is.character(species)) species
                   else if (!is.null(attr(y, "names_from"))) attr(y, "names_from")
                   else paste0("sp", seq_len(n_species))
  if (length(species_names) != n_species) {
    stop(sprintf("species has %d names but y has %d species",
                 length(species_names), n_species), call. = FALSE)
  }
  .tobs_check_site_count(n_sites, nrow(data), "sites")

  X_occ <- stats::model.matrix(occ_formula, data)
  X_det <- stats::model.matrix(det_formula, data)

  # Per-species integer detection / validity arrays.
  y_int <- array(0L,    dim = dim(y))
  valid <- array(FALSE, dim = dim(y))
  for (s in seq_len(n_species)) {
    ys <- matrix(as.integer(round(y[, , s])), n_sites, max_visits)
    vs <- !is.na(ys)
    if (any(ys[vs] != 0L & ys[vs] != 1L)) {
      stop(sprintf("species '%s': y must contain only 0, 1, or NA.",
                   species_names[s]), call. = FALSE)
    }
    ys[!vs] <- 0L
    y_int[, , s] <- ys
    valid[, , s] <- vs
  }

  # Per-species detection summaries (reuse the integrated single-source kernel).
  summaries <- lapply(seq_len(n_species), function(s) {
    .ms_int_occu_sp_summary(list(y_int[, , s]), list(valid[, , s]))
  })

  structure(list(
    model_type    = "ms_occu",
    y             = y_int,
    valid         = valid,
    n_sites       = n_sites,
    max_visits    = max_visits,
    n_species     = n_species,
    species_names = species_names,
    X_occ         = X_occ,
    X_det         = X_det,
    summaries     = summaries,
    structured_terms = structured_terms,
    formulas      = list(occ = occ_formula, det = det_formula),
    data          = data,
    process_info  = list(
      list(name = "psi", p = ncol(X_occ),
           coef_names = colnames(X_occ), link = "logit"),
      list(name = "p",   p = ncol(X_det),
           coef_names = colnames(X_det), link = "logit")
    )
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Laplace-EM fitter (shared community engine)
# ---------------------------------------------------------------------------

# Fit the community single-season occupancy model via the shared community
# Laplace-EM engine. The RE arms (psi, p) live in `theta`; there are no shared
# globals. Returns a `tobs_fit` (via build_ms_occu_fit).
.tobs_fit_ms_occu <- function(model,
                              priors     = NULL,
                              max.iter   = 200L,
                              tol        = 1e-4,
                              sigma.beta = 5,
                              verbose    = TRUE,
                              ...) {
  dots       <- list(...)
  newton.max <- as.integer(dots$newton.max %||% 30L)

  pi_list <- model$process_info
  P_occ <- pi_list[[1L]]$p
  P_p   <- pi_list[[2L]]$p
  P     <- P_occ + P_p
  S     <- model$n_species

  occ_idx <- seq_len(P_occ)
  p_idx   <- P_occ + seq_len(P_p)
  arm_idx <- list(psi = occ_idx, p = p_idx)

  X_occ <- model$X_occ
  X_det <- model$X_det
  summaries <- model$summaries

  # Single detection source: wrap the per-site detection eta in a length-1 list
  # so the integrated-community kernel (R/ms_int_occu.R) applies directly.
  sp_ll <- function(s, theta, global) {
    eta_psi <- as.numeric(X_occ %*% theta[occ_idx])
    eta_p   <- as.numeric(X_det %*% theta[p_idx])
    .ms_int_occu_sp_ll(eta_psi, list(eta_p), summaries[[s]])
  }
  sp_grad <- function(s, theta, global) {
    eta_psi <- as.numeric(X_occ %*% theta[occ_idx])
    eta_p   <- as.numeric(X_det %*% theta[p_idx])
    .ms_int_occu_sp_grad(eta_psi, list(eta_p), summaries[[s]],
                         X_occ, list(X_det))
  }

  # ---- warm start ----
  clp <- function(q) min(max(q, 1e-3), 1 - 1e-3)
  any_det_prop <- mean(vapply(summaries, function(z) mean(z$any_det), numeric(1)))
  det_n <- sum(vapply(summaries, function(z) sum(z$n_det[, 1L]), numeric(1)))
  val_n <- sum(vapply(summaries, function(z) sum(z$n_valid[z$any_det, 1L]),
                      numeric(1)))
  rate  <- if (val_n > 0) det_n / val_n else 0.5
  mu <- numeric(P)
  mu[occ_idx][1L] <- stats::qlogis(clp(any_det_prop))
  mu[p_idx][1L]   <- stats::qlogis(clp(rate))

  fit <- .tobs_community_em(
    S = S, P = P, arm_idx = arm_idx,
    sp_ll = sp_ll, sp_grad = sp_grad,
    init_mu = mu, init_global = numeric(0),
    penalize_global = FALSE, sigma_beta = sigma.beta, priors = priors,
    sigma_init = 0.3, max_iter = as.integer(max.iter), tol = as.numeric(tol),
    newton_max = newton.max, verbose = isTRUE(verbose)
  )

  build_ms_occu_fit(model, fit, arm_idx)
}


# ---------------------------------------------------------------------------
# Wrap the EM output into a tobs_fit
# ---------------------------------------------------------------------------

build_ms_occu_fit <- function(model, fit, arm_idx) {
  pi_list <- model$process_info

  beta_names <- c(
    paste0("psi_", pi_list[[1L]]$coef_names),
    paste0("p_",   pi_list[[2L]]$coef_names)
  )
  par_names <- beta_names

  means <- fit$mu; names(means) <- par_names
  V <- fit$Vf; dimnames(V) <- list(par_names, par_names)
  V <- (V + t(V)) / 2
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  B <- do.call(rbind, fit$b_list)              # S x P
  arm_block <- function(arm) {
    idx  <- arm_idx[[arm]]
    blup <- B[, idx, drop = FALSE]
    coef <- sweep(blup, 2L, means[idx], "+")
    rownames(blup) <- rownames(coef) <- model$species_names
    list(blup = blup, coef = coef)
  }
  occ_b <- arm_block("psi"); p_b <- arm_block("p")
  colnames(occ_b$blup) <- colnames(occ_b$coef) <- pi_list[[1L]]$coef_names
  colnames(p_b$blup)   <- colnames(p_b$coef)   <- pi_list[[2L]]$coef_names

  Sigma_occ <- fit$Sigma$psi; Sigma_p <- fit$Sigma$p
  dimnames(Sigma_occ) <- list(pi_list[[1L]]$coef_names, pi_list[[1L]]$coef_names)
  dimnames(Sigma_p)   <- list(pi_list[[2L]]$coef_names, pi_list[[2L]]$coef_names)

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(fit$logML, n_draws),
    log_lik      = fit$logML,
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
      Sigma_psi = Sigma_occ, Sigma_p = Sigma_p,
      sd_psi = sqrt(pmax(diag(Sigma_occ), 0)),
      sd_p   = sqrt(pmax(diag(Sigma_p),   0)),
      coef_psi = occ_b$coef, coef_p = p_b$coef,
      blup_psi = occ_b$blup, blup_p = p_b$blup
    ),
    convergence  = list(converged = isTRUE(fit$converged), n_iter = fit$n_iter)
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "ms_occu")
# ---------------------------------------------------------------------------

# Per-species BLUP deviations, long form: one row per (species, arm, term).
.tobs_ranef_ms_occu <- function(object) {
  .tobs_ranef_ms_long(object$ms_community,
                      c(psi = "blup_psi", p = "blup_p"))
}

# Per-species posterior-mean linear predictors: site-level occupancy psi and
# detection p [n_sites x n_species], and the per-species occupancy posterior z
# = P(z = 1 | y) [n_sites x n_species].
.tobs_fitted_ms_occu <- function(object) {
  model <- object$model
  cm    <- object$ms_community
  eta_psi <- model$X_occ %*% t(cm$coef_psi)
  # Shared areal field (nested_laplace path): add the per-site occupancy offset
  # to every species' occupancy predictor.
  if (!is.null(object$spatial_field)) {
    eta_psi <- sweep(eta_psi, 1L, as.numeric(object$spatial_field), "+")
  }
  psi <- stats::plogis(eta_psi)
  p   <- stats::plogis(model$X_det %*% t(cm$coef_p))

  z <- matrix(0, model$n_sites, model$n_species)
  for (s in seq_len(model$n_species)) {
    summ <- model$summaries[[s]]
    nv <- summ$n_valid[, 1L]; ad <- summ$any_det
    prod1mp <- (1 - p[, s])^nv
    num <- psi[, s] * prod1mp
    z[, s] <- ifelse(ad, 1, num / (num + (1 - psi[, s])))
  }
  dimnames(psi) <- dimnames(p) <- dimnames(z) <-
    list(NULL, model$species_names)
  list(psi = psi, p = p, z = z)
}

# Draw community single-season data under the fitted per-species coefficients,
# at the observed visit pattern. Returns a 3D array matching the input y.
.tobs_simulate_ms_occu <- function(object, nsim = 1) {
  model <- object$model
  cm    <- object$ms_community
  n_sites <- model$n_sites; max_visits <- model$max_visits
  n_species <- model$n_species

  one <- function() {
    y_sim <- array(NA_integer_, dim = c(n_sites, max_visits, n_species),
                   dimnames = list(NULL, NULL, model$species_names))
    for (s in seq_len(n_species)) {
      psi <- stats::plogis(as.numeric(model$X_occ %*% cm$coef_psi[s, ]))
      p   <- stats::plogis(as.numeric(model$X_det %*% cm$coef_p[s, ]))
      z   <- stats::rbinom(n_sites, 1L, psi)
      vs  <- model$valid[, , s]
      for (i in seq_len(n_sites)) {
        vis <- which(vs[i, ])
        if (!length(vis)) next
        y_sim[i, vis, s] <- stats::rbinom(length(vis), 1L, z[i] * p[i])
      }
    }
    y_sim
  }
  if (nsim == 1L) one() else lapply(seq_len(nsim), function(i) one())
}


# ---------------------------------------------------------------------------
# Species richness
# ---------------------------------------------------------------------------

# Per-site expected species richness sum_s psi_{s,i} under the fitted per-species
# occupancy coefficients, with a posterior credible interval propagated from the
# community-mean coefficient draws (per-species deviations held at their BLUPs).
.tobs_richness_ms_occu <- function(object) {
  model <- object$model
  cm    <- object$ms_community
  draws <- object$draws
  X_occ <- model$X_occ
  P_occ <- model$process_info[[1L]]$p
  n_sites   <- model$n_sites
  n_species <- model$n_species
  n_draws   <- nrow(draws)

  # Per-draw community-mean occupancy coefficient; add each species' BLUP
  # deviation to get the per-species occupancy linear predictor.
  beta_occ_draws <- draws[, seq_len(P_occ), drop = FALSE]
  blup <- cm$blup_psi                                # S x P_occ
  richness <- matrix(0, n_draws, n_sites)
  for (d in seq_len(n_draws)) {
    mu_occ <- beta_occ_draws[d, ]
    psi_si <- vapply(seq_len(n_species), function(s) {
      stats::plogis(as.numeric(X_occ %*% (mu_occ + blup[s, ])))
    }, numeric(n_sites))                              # n_sites x n_species
    richness[d, ] <- rowSums(psi_si)
  }

  data.frame(
    site  = seq_len(n_sites),
    mean  = colMeans(richness),
    sd    = apply(richness, 2, stats::sd),
    q2.5  = apply(richness, 2, stats::quantile, 0.025),
    q97.5 = apply(richness, 2, stats::quantile, 0.975)
  )
}


# ---------------------------------------------------------------------------
# Family constructor
# ---------------------------------------------------------------------------

#' Multispecies (community) single-season occupancy family
#'
#' Per-species occupancy and detection with Gaussian community hyperpriors on the
#' per-species coefficients of each arm (the spOccupancy `msPGOcc` model). The
#' occupancy and detection per-species deviations are independent, each with its
#' own community covariance, fit by a community Laplace-EM.
#'
#' @return A `tobs_family` object.
#' @seealso [occu()], [ms_dyn_occu()], [ms_int_occu()]
#' @examples
#' \donttest{
#' sim <- simulate_ms_occu(N = 80, J = 4, n_species = 8, seed = 1)
#' fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
#'             y = sim$y, species = paste0("sp", seq_len(8)), method = "laplace")
#' summary(fit)
#' }
#' @export
ms_occu <- function() {
  obs_family(
    name           = "ms_occu",
    class_long     = "multispecies occupancy",
    latent         = "bernoulli",
    observation    = "binomial_detection",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working"
  )
}
