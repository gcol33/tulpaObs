# =============================================================================
# ms_occu_field.R - community occupancy with a shared areal field / varying-
# coefficient field on the occupancy arm via block coordinate ascent (the
# spOccupancy sfMsPGOcc / svcMsPGOcc analogue; gcol33/tulpaObs#117, #118).
#
#   logit psi_{s,i} = X_i . (mu_psi + b_psi_s) + sum_k W[i,k] F[u(i), k]
#   logit p_{s,i}   = Xdet_i . (mu_p + b_p_s)
#
# Reuses the pure-R community occupancy Laplace-EM (.tobs_community_em with the
# ms_int_occu two-state marginal) unchanged: the shared occupancy field enters
# each species' psi as a fixed per-site offset in the sp_ll closure, so a
# coefficient update is an ordinary community occupancy EM; the field update
# given the coefficients is a self-contained occupancy-marginal Laplace (a Newton
# over the two-state per-(species, site) marginal + the ICAR prior). One field
# node per site. This is the occupancy counterpart of ms_count_spatial.R.
# =============================================================================


# Occupancy-marginal field update for K covariate-weighted ICAR fields on the
# occupancy arm. `o_occ` [n_sites x n_species] is the per-species psi linear
# predictor from the coefficients (X . (mu + b_s)); `p_site` the per-species
# detection probabilities; `su` the per-species detection summaries (n_valid,
# any_det). Per site the two-state marginal score / curvature aggregate over
# species, weighted by W. Joint Newton over the K*n field vector, then a per-
# field tau M-step. The field is sum-to-zero (the intercept owns the level).
.ms_occu_field_solve <- function(o_occ, p_site, su, F, tau, W, Q,
                                 max_iter = 40L, tol = 1e-7) {
  Ns <- nrow(o_occ); S <- ncol(o_occ); K <- ncol(W)
  build_H <- function(cw) {
    blocks <- vector("list", K * K)
    for (k in seq_len(K)) for (l in seq_len(K)) {
      D <- Matrix::Diagonal(x = W[, k] * W[, l] * cw)
      blocks[[(k - 1L) * K + l]] <- if (k == l) D + tau[k] * Q else D
    }
    do.call(rbind, lapply(seq_len(K), function(k)
      do.call(cbind, blocks[((k - 1L) * K + 1L):(k * K)])))
  }
  for (it in seq_len(max_iter)) {
    field_off <- rowSums(W * F)
    sc <- numeric(Ns); cw <- numeric(Ns)
    for (s in seq_len(S)) {
      psi <- stats::plogis(o_occ[, s] + field_off); s1 <- psi * (1 - psi)
      nv  <- su[[s]]$n_valid[, 1L]; ad <- su[[s]]$any_det
      q   <- (1 - p_site[, s])^nv
      L   <- psi * q + (1 - psi)
      # D=1 (any detection): score 1-psi, neg-curvature psi(1-psi).
      # D=0: L = psi q + (1-psi); a PSD Fisher-like curvature keeps Newton stable.
      sc  <- sc + ifelse(ad, 1 - psi, s1 * (q - 1) / L)
      cw  <- cw + ifelse(ad, s1, pmax(s1 * (1 - q) / L, 1e-8))
    }
    g <- unlist(lapply(seq_len(K), function(k)
      W[, k] * sc - tau[k] * as.numeric(Q %*% F[, k])))
    step <- as.numeric(Matrix::solve(build_H(cw), g))
    F <- F + matrix(step, Ns, K)
    for (k in seq_len(K)) F[, k] <- F[, k] - mean(F[, k])
    if (max(abs(step)) < tol) break
  }
  field_off <- rowSums(W * F)
  cw <- numeric(Ns)
  for (s in seq_len(S)) {
    psi <- stats::plogis(o_occ[, s] + field_off); cw <- cw + psi * (1 - psi)
  }
  Cov <- Matrix::solve(build_H(cw))
  tau_new <- numeric(K)
  for (k in seq_len(K)) {
    idx  <- (k - 1L) * Ns + seq_len(Ns)
    quad <- as.numeric(t(F[, k]) %*% (Q %*% F[, k])) +
            sum(Matrix::diag(Q %*% Cov[idx, idx, drop = FALSE]))
    tau_new[k] <- (Ns - 1) / max(quad, 1e-8)
  }
  list(F = F, tau = tau_new)
}


# Fit the community-spatial / SVC occupancy model. `model` is the ms_occu model;
# `spatial` the resolved icar term or bar (intercept + varying-coefficient
# fields). Block coordinate ascent between the community occupancy EM (field as a
# psi offset) and the occupancy-marginal field update.
.tobs_fit_ms_occu_field <- function(model, spatial, max.iter = 200L, tol = 1e-4,
                                    sigma.beta = 5, priors = NULL,
                                    max.outer = 20L, verbose = FALSE, ...) {
  pi_list <- model$process_info
  P_occ <- pi_list[[1L]]$p; P_p <- pi_list[[2L]]$p; P <- P_occ + P_p
  S <- model$n_species; Ns <- model$n_sites
  occ_idx <- seq_len(P_occ); p_idx <- P_occ + seq_len(P_p)
  arm_idx <- list(psi = occ_idx, p = p_idx)
  Xocc <- model$X_occ; Xdet <- model$X_det; su <- model$summaries

  fields <- .tobs_resolve_occu_spatial_fields(spatial, model)
  A <- fields[[1L]]$graph
  if (is.null(A) || nrow(A) != Ns) {
    stop("ms_occu() shared occupancy field needs an icar() graph with one node ",
         "per site (gcol33/tulpaObs#118).", call. = FALSE)
  }
  if (!all(vapply(fields, function(fd) identical(fd$type %||% "icar", "icar"),
                  logical(1)))) {
    stop("ms_occu() occupancy field / SVC supports icar() (bym2()/car_proper() ",
         "are follow-ups, gcol33/tulpaObs#118).", call. = FALSE)
  }
  K <- length(fields); W <- matrix(1, Ns, K); field_labels <- character(K)
  for (k in seq_len(K)) {
    wk <- fields[[k]]$weight
    if (!is.null(wk)) {
      if (length(wk) != Ns)
        stop("ms_occu() varying-coefficient field weight must be one per site.",
             call. = FALSE)
      W[, k] <- as.numeric(wk)
      field_labels[k] <- fields[[k]]$weight_label %||% paste0("trend", k - 1L)
    } else field_labels[k] <- "intercept"
  }
  Q <- Matrix::Diagonal(x = rowSums(A)) - methods::as(A, "CsparseMatrix")

  Ffield <- matrix(0, Ns, K); tau <- rep(1, K); em <- NULL
  clp <- function(z) min(max(z, 1e-3), 1 - 1e-3)
  mu0 <- numeric(P)
  mu0[occ_idx][1L] <- stats::qlogis(clp(mean(vapply(su, function(z)
    mean(z$any_det), numeric(1)))))

  for (outer in seq_len(max.outer)) {
    field_off <- rowSums(W * Ffield)
    sp_ll <- function(s, theta, global) {
      ep <- as.numeric(Xocc %*% theta[occ_idx]) + field_off
      .ms_int_occu_sp_ll(ep, list(as.numeric(Xdet %*% theta[p_idx])), su[[s]])
    }
    sp_grad <- function(s, theta, global) {
      ep <- as.numeric(Xocc %*% theta[occ_idx]) + field_off
      .ms_int_occu_sp_grad(ep, list(as.numeric(Xdet %*% theta[p_idx])), su[[s]],
                           Xocc, list(Xdet))
    }
    em <- .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      init_mu = if (is.null(em)) mu0 else em$mu, init_global = numeric(0),
      penalize_global = FALSE, sigma_beta = sigma.beta, priors = priors,
      sigma_init = 0.3, max_iter = min(as.integer(max.iter), 50L),
      tol = as.numeric(tol), newton_max = 30L, verbose = FALSE)

    o_occ  <- vapply(seq_len(S), function(s)
      as.numeric(Xocc %*% (em$mu[occ_idx] + em$b_list[[s]][occ_idx])), numeric(Ns))
    p_site <- vapply(seq_len(S), function(s)
      stats::plogis(as.numeric(Xdet %*% (em$mu[p_idx] + em$b_list[[s]][p_idx]))),
      numeric(Ns))
    fu <- .ms_occu_field_solve(o_occ, p_site, su, Ffield, tau, W, Q)
    delta <- max(abs(fu$F - Ffield)); Ffield <- fu$F; tau <- fu$tau
    if (isTRUE(verbose))
      message(sprintf("[ms_occu field %d] delta=%.2e", outer, delta))
    if (outer > 2L && delta < tol) break
  }

  fit <- build_ms_occu_fit(model, em, arm_idx)
  fit$method <- "nested_laplace"
  fit$spatial <- spatial
  fit$spatial_field <- as.numeric(Ffield[, 1L])
  sigma_field <- 1 / sqrt(tau)
  fit$spatial_hyper <- list(type = "icar", tau = tau, sigma = sigma_field,
                            field_labels = field_labels)
  if (K > 1L) {
    fit$trend_fields <- lapply(2:K, function(k) as.numeric(Ffield[, k]))
    names(fit$trend_fields) <- field_labels[-1L]
    fit$trend_field <- fit$trend_fields[[1L]]
  }
  fit
}
