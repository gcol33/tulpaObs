# =============================================================================
# ms_occu_field.R - community occupancy with a shared latent structure on the
# occupancy arm: a shared areal field / varying-coefficient field (the
# spOccupancy sfMsPGOcc / svcMsPGOcc analogues), latent factors (lfMsPGOcc), or
# BOTH (the spatial-factor case). gcol33/tulpaObs#117, #118, #119.
#
#   logit psi_{s,i} = X_i . (mu_psi + b_psi_s) + sum_k W[i,k] F[u(i),k]
#                                              + sum_q lambda_{s,q} zeta_{q,i}
#   logit p_{s,i}   = Xdet_i . (mu_p + b_p_s)
#
# The block coordinate ascent, the areal Newton, the factor update, and the field
# hyperparameter grids live in R/community_latent.R and are shared with every
# other community family. This file supplies only the occupancy two-state
# working oracle and the model wiring. One field node per site (or a group_var
# site -> cell map).
#
# Information limit when a field and factors are combined (sfMsPGOcc). The
# centred loadings (sum_s lambda_sq = 0) make the shared field (loading == 1) and
# the factors orthogonal in species space, so both are identified, but they
# compete for the same per-site signal. A measured sweep
# (dev_notes/probe_sf_occu_diagnose.R) shows the residual-correlation recovery
# tracks FIELD STRENGTH -- median 0.83 / 0.72 / 0.65 at field SD 0.5 / 1.0 / 2.0
# -- because a strong field saturates psi toward 0/1, where the binary curvature
# psi (1 - psi) collapses and a detection history carries almost no information
# about species-specific deviations. It is not a convergence failure (20 vs 60
# outer iterations move it by 0.002) and not an artefact of carrying a field
# (a field absent from the data costs ~nothing). More species / visits restores
# it (J = 10, S = 30 -> median 0.91). The field itself recovers throughout
# (>= 0.95), so a field-only fit is unaffected. The Poisson community analogue
# (ms_count) does not hit this: a count is information-rich per (site, species).
# =============================================================================


# Occupancy two-state working oracle for the occupancy arm. Integrating the
# latent z out per (species, site) leaves the marginal
#
#   D = 1 (any detection): L = psi (1 - q),   q = (1 - p)^n_valid
#   D = 0                : L = psi q + (1 - psi)
#
# whose score / curvature wrt an additive offset on logit psi are the pair below.
# `eta_p` [n_sites x n_species] is the per-species detection linear predictor
# from the current coefficients, held fixed across a field / factor update.
.tobs_ms_occu_oracle <- function(su, eta_p) {
  Ns <- nrow(eta_p); S <- ncol(eta_p)
  p_site <- stats::plogis(eta_p)
  nv <- vapply(su, function(z) z$n_valid[, 1L], numeric(Ns))    # n_sites x n_species
  ad <- vapply(su, function(z) z$any_det,       logical(Ns))    # n_sites x n_species
  qmat <- (1 - p_site)^nv
  # `idx` (a subset of site rows) restricts the per-species marginal to those
  # sites; `su[[s]]` (fixed per-site detection summaries) is subsetted to match
  # (gcol33/tulpaObs#162 lever 2).
  summ_subset <- function(summ, idx)
    list(n_valid = summ$n_valid[idx, , drop = FALSE],
         n_det   = summ$n_det[idx, , drop = FALSE],
         any_det = summ$any_det[idx])
  ll_cell <- function(eta, idx = NULL) {
    ii <- idx %||% seq_len(Ns)
    # vapply degenerates to a plain vector (not a length(ii) x S matrix) when
    # length(ii) == 1, since a length-1 FUN.VALUE never triggers its matrix
    # path -- only pending on one site is common in the mode-adaptation
    # backtracking tail, so this must be forced back into a matrix.
    matrix(vapply(seq_len(S), function(s)
      .ms_int_occu_sp_ll(eta[ii, s], list(eta_p[ii, s]),
                         summ_subset(su[[s]], ii), per_site = TRUE),
      numeric(length(ii))), nrow = length(ii))
  }
  list(
    n_sites   = Ns,
    n_species = S,
    working = function(eta) {
      psi <- stats::plogis(eta); s1 <- psi * (1 - psi)
      L   <- psi * qmat + (1 - psi)
      list(score = ifelse(ad, 1 - psi, s1 * (qmat - 1) / L),
           curv  = ifelse(ad, s1, pmax(s1 * (1 - qmat) / L, 1e-8)))
    },
    # The marginal curvature collapses toward zero at undetected sites (a site
    # with no detection barely constrains psi), which would blow up the Laplace
    # field covariance feeding the tau M-step. The complete-data Fisher curvature
    # psi (1 - psi) keeps that solve conditioned.
    cov_curv = function(eta) {
      psi <- stats::plogis(eta); psi * (1 - psi)
    },
    # Per-(site, species) marginal log-likelihood. The joint site marginal that
    # sets the factor magnitude integrates zeta with the species at a site kept
    # together, so it needs the cells rather than their sum.
    ll_cell = ll_cell,
    data_ll = function(eta) sum(ll_cell(eta)))
}


# Fit the community occupancy model with a shared areal / varying-coefficient
# field (`spatial`), latent factors (`latent`), or both. Single source of truth
# for every latent community-occupancy route.
.tobs_fit_ms_occu_field <- function(model, spatial = NULL, latent = NULL,
                                    max.iter = 200L, tol = 1e-4,
                                    sigma.beta = 5, priors = NULL,
                                    max.outer = NULL, factor.starts = NULL,
                                    n.quad = NULL, verbose = FALSE, ...) {
  pi_list <- model$process_info
  P_occ <- pi_list[[1L]]$p; P_p <- pi_list[[2L]]$p; P <- P_occ + P_p
  S <- model$n_species; Ns <- model$n_sites
  occ_idx <- seq_len(P_occ); p_idx <- P_occ + seq_len(P_p)
  arm_idx <- list(psi = occ_idx, p = p_idx)
  Xocc <- model$X_occ; Xdet <- model$X_det; su <- model$summaries

  clp <- function(z) min(max(z, 1e-3), 1 - 1e-3)
  mu0 <- numeric(P)
  mu0[occ_idx][1L] <- stats::qlogis(clp(mean(vapply(su, function(z)
    mean(z$any_det), numeric(1)))))

  em_fit <- function(site_off, fac_off, em_prev) {
    sp_ll <- function(s, theta, global) {
      ep <- as.numeric(Xocc %*% theta[occ_idx]) + site_off + fac_off[, s]
      .ms_int_occu_sp_ll(ep, list(as.numeric(Xdet %*% theta[p_idx])), su[[s]])
    }
    sp_grad <- function(s, theta, global) {
      ep <- as.numeric(Xocc %*% theta[occ_idx]) + site_off + fac_off[, s]
      .ms_int_occu_sp_grad(ep, list(as.numeric(Xdet %*% theta[p_idx])), su[[s]],
                           Xocc, list(Xdet))
    }
    .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      init_mu = if (is.null(em_prev)) mu0 else em_prev$mu,
      init_b = em_prev$b_list, init_Sigma = em_prev$Sigma,
      init_global = numeric(0), penalize_global = FALSE,
      sigma_beta = sigma.beta, priors = priors, sigma_init = 0.3,
      max_iter = min(as.integer(max.iter), 50L), tol = as.numeric(tol),
      newton_max = 30L, verbose = FALSE)
  }
  offset_of <- function(em) {
    vapply(seq_len(S), function(s)
      as.numeric(Xocc %*% (em$mu[occ_idx] + em$b_list[[s]][occ_idx])),
      numeric(Ns))
  }
  make_oracle <- function(em) {
    eta_p <- vapply(seq_len(S), function(s)
      as.numeric(Xdet %*% (em$mu[p_idx] + em$b_list[[s]][p_idx])), numeric(Ns))
    .tobs_ms_occu_oracle(su, eta_p)
  }

  res <- .tobs_community_latent_ascent(
    spatial = spatial, latent = latent, model = model, what = "ms_occu()",
    make_oracle = make_oracle, em_fit = em_fit, offset_of = offset_of,
    # 150 is the ms_count-measured factor budget (gcol33/tulpaObs#156); this
    # family's 16-seed recovery was measured AT it (magnitude median 1.021,
    # slope z 0.71) rather than inheriting it untested.
    allow = "icar", tol = tol, max.outer = max.outer, factor.outer = 150L,
    # One starting direction, from this family's own measurement
    # (gcol33/tulpaObs#164). Over 10 seeds at N=250, S=16, Q=2, J=5 the
    # multi-start's 7 extra candidates cost 3.5-4.3x and, on 4 of the 8
    # reseeded seeds, landed on the EXACT SAME loadings as the 1-start fit
    # (d_mag = d_res = 0.0000) -- the eight directions are converging to one
    # basin, not escaping it. The one seed flagged by mag_ratio (1.40x truth)
    # reproduced identically at 8 starts too, so it is a hard fixture the
    # multi-start was never going to reach (the same character as the
    # ms_abun measurement's own seed 314). Family default 1;
    # `control$factor.starts` still overrides it.
    factor.starts = if (is.null(factor.starts)) 1L else factor.starts,
    n.quad = n.quad, verbose = verbose)

  fit <- build_ms_occu_fit(model, res$em, arm_idx)
  fit$method <- "laplace"
  fit <- .tobs_latent_attach_field(fit, res, spatial, "occu_field_offset")
  fit <- .tobs_latent_attach_factor(fit, res, latent, model,
                                    "occu_factor_offset")
  fit
}
