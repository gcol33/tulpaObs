# =============================================================================
# distance_spatial.R - areal-spatial binned distance-sampling abundance (#51)
#
# An ICAR or proper-CAR field on the abundance arm of the distance marginal, fit
# by nested Laplace (outer grid over tau[, rho][, NB size r], inner Newton over
# (beta_lambda, beta_sigma, z)). The distance marginal is NOT the binomial-
# detection N-mixture (the count-spatial driver is built around that), so this is
# a dedicated fitter -- but it reuses the exact per-site moments the distance
# kernel already computes (cpp_distance_site_sweep) and the CAR prior helpers.
#
# Per-site marginal observed-information in eta-space is the form distance_kernel.h
# documents: diag(info_lambda, info_sigma) - var_N v v', v = (score_wt_lambda,
# vN_sigma) with vN_sigma = -p_sigma/(1-p) the detection N-coupling. z loads onto
# eta_lambda exactly like the abundance intercept (one unit per site), so the
# field couples to (beta_lambda, beta_sigma) through the same per-site blocks.
# Half-normal key only (the hazard shape is a global coordinate, deferred).
#
#   .tobs_fit_distance_spatial()   dispatch from .tobs_fit_model (icar / car_proper)
# =============================================================================

# CAR precision Q(rho) = D - rho W from a dense adjacency (D = degree). ICAR is
# rho = 1 (rank-deficient); proper CAR uses rho < 1 (full rank).
.dist_spatial_Q <- function(adj, rho) {
  deg <- rowSums(adj != 0)
  diag(deg) - rho * adj
}

# One inner Newton at fixed (tau, rho) over x = (beta_lambda, beta_sigma, z).
.dist_spatial_inner <- function(sweep_fn, X_lam, X_sig, map, Q, tau, rho,
                                kind, log_det_Q, beta_lam, beta_sig, z,
                                max_iter, tol, verbose) {
  p_lam <- ncol(X_lam); p_sig <- ncol(X_sig); n_sp <- length(z)
  n_sites <- nrow(X_lam); n_x <- p_lam + p_sig + n_sp
  li <- seq_len(p_lam); si <- p_lam + seq_len(p_sig); zi <- p_lam + p_sig + seq_len(n_sp)
  cl <- function(e) pmin(pmax(e, -30), 30)

  car_logprior <- function(zz) {
    quad <- as.numeric(t(zz) %*% Q %*% zz)
    if (kind == "icar") -0.5 * tau * quad + 0.5 * (n_sp - 1) * log(tau)
    else 0.5 * log_det_Q + 0.5 * n_sp * log(tau) - 0.5 * tau * quad
  }
  obj_ll <- function(bl, bs, zz) {
    eta_lam <- cl(as.numeric(X_lam %*% bl) + zz[map])
    eta_sig <- cl(as.numeric(X_sig %*% bs))
    sw <- sweep_fn(eta_lam, eta_sig)
    list(ll = sum(sw$log_lik), sw = sw)
  }

  grad_norm <- Inf; conv <- FALSE; n_iter <- 0L
  for (iter in seq_len(max_iter)) {
    eta_lam <- cl(as.numeric(X_lam %*% beta_lam) + z[map])
    eta_sig <- cl(as.numeric(X_sig %*% beta_sig))
    sw <- sweep_fn(eta_lam, eta_sig)
    h_ll <- sw$info_lam - sw$var_N * sw$swl^2
    h_ss <- sw$info_sig_obs - sw$var_N * sw$vN_sig^2
    h_ls <- -sw$var_N * sw$swl * sw$vN_sig
    # PSD Fisher-scoring blocks (no var_N correction, no cross) for the fallback.
    f_ll <- sw$info_lam; f_ss <- sw$info_sig_fs

    g <- numeric(n_x)
    g[li] <- as.numeric(crossprod(X_lam, sw$grad_lam))
    g[si] <- as.numeric(crossprod(X_sig, sw$grad_sig))
    gz <- numeric(n_sp)
    for (s in seq_len(n_sites)) gz[map[s]] <- gz[map[s]] + sw$grad_lam[s]
    g[zi] <- gz - tau * as.numeric(Q %*% z)

    asm <- function(hll, hss, hls, with_cross) {
      H <- matrix(0, n_x, n_x)
      H[li, li] <- crossprod(X_lam, X_lam * hll)
      H[si, si] <- crossprod(X_sig, X_sig * hss)
      if (with_cross) { H[li, si] <- crossprod(X_lam, X_sig * hls); H[si, li] <- t(H[li, si]) }
      ZL <- matrix(0, p_lam, n_sp); ZS <- matrix(0, p_sig, n_sp); zz_diag <- numeric(n_sp)
      for (s in seq_len(n_sites)) {
        u <- map[s]
        ZL[, u] <- ZL[, u] + X_lam[s, ] * hll[s]
        if (with_cross) ZS[, u] <- ZS[, u] + X_sig[s, ] * hls[s]
        zz_diag[u] <- zz_diag[u] + hll[s]
      }
      H[li, zi] <- ZL; H[zi, li] <- t(ZL)
      if (with_cross) { H[si, zi] <- ZS; H[zi, si] <- t(ZS) }
      H[zi, zi] <- diag(zz_diag, n_sp) + tau * Q
      H
    }
    H <- asm(h_ll, h_ss, h_ls, TRUE)
    ridge <- max(1e-10 * mean(diag(H)), 1e-12)
    diag(H) <- diag(H) + ridge
    ch <- tryCatch(chol(H), error = function(e) NULL)
    if (is.null(ch)) {
      Hf <- asm(f_ll, f_ss, NULL, FALSE)
      diag(Hf) <- diag(Hf) + ridge
      ch <- tryCatch(chol(Hf), error = function(e) NULL)
      if (is.null(ch)) break
    }
    delta <- backsolve(ch, backsolve(ch, g, transpose = TRUE))
    grad_norm <- sqrt(sum(g^2))
    if (verbose) cat(sprintf("  iter %d  ll? grad_norm %.3e\n", iter, grad_norm))
    if (grad_norm < tol) { conv <- TRUE; n_iter <- iter; break }

    cur <- obj_ll(beta_lam, beta_sig, z)
    obj_cur <- cur$ll + car_logprior(z)
    step <- 1; stepped <- FALSE
    for (h in 0:11) {
      bl <- beta_lam + step * delta[li]; bs <- beta_sig + step * delta[si]
      zz <- z + step * delta[zi]
      tr <- obj_ll(bl, bs, zz); obj_try <- tr$ll + car_logprior(zz)
      if (is.finite(obj_try) && obj_try >= obj_cur - 1e-10) {
        beta_lam <- bl; beta_sig <- bs; z <- zz
        if (kind == "icar") z <- z - mean(z)   # sum-to-zero
        stepped <- TRUE; break
      }
      step <- step * 0.5
    }
    if (!stepped) break
    n_iter <- iter
  }

  # log marginal at the mode
  eta_lam <- cl(as.numeric(X_lam %*% beta_lam) + z[map])
  eta_sig <- cl(as.numeric(X_sig %*% beta_sig))
  sw <- sweep_fn(eta_lam, eta_sig)
  h_ll <- sw$info_lam - sw$var_N * sw$swl^2
  h_ss <- sw$info_sig_obs - sw$var_N * sw$vN_sig^2
  h_ls <- -sw$var_N * sw$swl * sw$vN_sig
  H <- {
    Hm <- matrix(0, n_x, n_x)
    Hm[li, li] <- crossprod(X_lam, X_lam * h_ll)
    Hm[si, si] <- crossprod(X_sig, X_sig * h_ss)
    Hm[li, si] <- crossprod(X_lam, X_sig * h_ls); Hm[si, li] <- t(Hm[li, si])
    ZL <- matrix(0, p_lam, n_sp); ZS <- matrix(0, p_sig, n_sp); zzd <- numeric(n_sp)
    for (s in seq_len(n_sites)) { u <- map[s]
      ZL[, u] <- ZL[, u] + X_lam[s, ] * h_ll[s]
      ZS[, u] <- ZS[, u] + X_sig[s, ] * h_ls[s]
      zzd[u] <- zzd[u] + h_ll[s] }
    Hm[li, zi] <- ZL; Hm[zi, li] <- t(ZL)
    Hm[si, zi] <- ZS; Hm[zi, si] <- t(ZS)
    Hm[zi, zi] <- diag(zzd, n_sp) + tau * Q
    Hm
  }
  ridge <- max(1e-10 * mean(diag(H)), 1e-12); diag(H) <- diag(H) + ridge
  p_beta <- p_lam + p_sig
  ch <- tryCatch(chol(H), error = function(e) NULL)
  if (is.null(ch)) {
    return(list(ok = FALSE))
  }
  log_det_H <- 2 * sum(log(diag(ch)))
  # constrained coef cov: penalise sum(z)=0 for the rank-deficient ICAR field.
  Hc <- H
  if (kind == "icar") {
    pen <- 1e6 * mean(diag(H))
    Hc[zi, zi] <- Hc[zi, zi] + pen   # add pen * 1 1' on the field block
  }
  cov_full <- tryCatch(solve(Hc), error = function(e) NULL)
  if (is.null(cov_full)) return(list(ok = FALSE))
  cov_beta <- cov_full[seq_len(p_beta), seq_len(p_beta), drop = FALSE]

  list(ok = TRUE, beta_lam = beta_lam, beta_sig = beta_sig, z = z,
       cov_beta = cov_beta, log_lik = sum(sw$log_lik),
       log_marginal = sum(sw$log_lik) + car_logprior(z) - 0.5 * log_det_H,
       boundary_max = max(sw$boundary, na.rm = TRUE),
       converged = conv, n_iter = n_iter)
}

# Areal-spatial distance fit dispatched from .tobs_fit_model. icar() / car_proper()
# on one spatial unit per site (the field on log lambda); bym2 / spde / hazard key
# / weighted (SVC) fields are deferred. Poisson or NB (the NB size integrated over
# the outer grid). Packs through build_distance_fit with the grid-integrated
# coefficient means / covariance and the posterior-mean field.
.tobs_fit_distance_spatial <- function(model, spatial, mixture = "poisson",
                                       K_max = NULL, max_iter = 100L, tol = 1e-6,
                                       verbose = TRUE) {
  .tobs_reject_weighted_spatial(spatial, "distance abundance spatial")
  if (!identical(model$key, "halfnorm")) {
    stop("distance() areal spatial supports the half-normal key only; the ",
         "hazard-rate shape is a global coordinate not yet wired into the ",
         "spatial path. (tulpaObs#51)", call. = FALSE)
  }
  if (spatial$type %in% c("bym2", "spde", "gp", "multiscale_gp")) {
    stop(sprintf(paste0("distance() areal spatial supports icar() / car_proper() ",
                        "under method = \"nested_laplace\"; the '%s' field is not ",
                        "yet wired for distance. (tulpaObs#51)"), spatial$type),
         call. = FALSE)
  }
  if (!spatial$type %in% c("icar", "car_proper")) {
    stop(sprintf("distance() areal spatial supports icar() / car_proper(); got '%s'.",
                 spatial$type), call. = FALSE)
  }
  n_sites <- model$n_sites
  if (spatial$n_units != n_sites) {
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for distance.",
                 spatial$n_units, n_sites), call. = FALSE)
  }
  adj <- if (!is.null(spatial$graph)) as.matrix(spatial$graph) else
    stop("distance() spatial term must carry an adjacency graph.", call. = FALSE)
  n_sp <- spatial$n_units
  map <- seq_len(n_sites)

  X_lam <- model$X_processes[[1]]; X_sig <- model$X_processes[[2]]
  y <- matrix(as.integer(model$y), nrow(model$y), ncol(model$y))
  R_max <- if (length(y)) max(rowSums(y)) else 0L
  K_max <- if (is.null(K_max)) as.integer(3L * R_max + 100L) else as.integer(K_max)
  transect_code <- .dist_transect_code(model$transect)
  quad_order <- model$quad_order %||% 64L
  is_nb <- mixture %in% c("negbin", "NB")
  r_grid <- if (is_nb) exp(seq(log(0.5), log(40), length.out = 6L)) else Inf
  tau_grid <- exp(seq(log(0.3), log(30), length.out = 9L))
  rho_grid <- if (identical(spatial$type, "car_proper")) seq(0.1, 0.95, length.out = 6L) else 1.0
  kind <- if (identical(spatial$type, "icar")) "icar" else "car_proper"

  beta_lam0 <- c(log(max(mean(rowSums(y)), 0.5) + 0.5), rep(0, ncol(X_lam) - 1L))
  beta_sig0 <- c(log(stats::median(model$cutpoints[-1])), rep(0, ncol(X_sig) - 1L))

  grid <- expand.grid(tau = tau_grid, rho = rho_grid, r = r_grid)
  n_grid <- nrow(grid)
  log_det_Q <- vapply(rho_grid, function(rho) {
    if (kind == "icar") 0 else {
      Qr <- .dist_spatial_Q(adj, rho)
      ch <- tryCatch(chol(Qr), error = function(e) NULL)
      if (is.null(ch)) -Inf else 2 * sum(log(diag(ch)))
    }
  }, numeric(1))
  names(log_det_Q) <- as.character(rho_grid)

  p_lam <- ncol(X_lam); p_sig <- ncol(X_sig); p_beta <- p_lam + p_sig
  modes <- matrix(NA_real_, n_grid, p_beta); fields <- matrix(NA_real_, n_grid, n_sp)
  cov_blocks <- vector("list", n_grid)
  logm <- rep(-Inf, n_grid); bmax <- rep(0, n_grid)
  for (k in seq_len(n_grid)) {
    rho <- grid$rho[k]; tau <- grid$tau[k]; rr <- grid$r[k]
    ldQ <- log_det_Q[[as.character(rho)]]
    if (!is.finite(ldQ)) next
    Q <- .dist_spatial_Q(adj, rho)
    sweep_fn <- function(eta_lam, eta_sig)
      cpp_distance_site_sweep(y, eta_lam, eta_sig, as.numeric(model$cutpoints),
                              transect_code, as.integer(quad_order), K_max,
                              nb = is_nb, r = rr)
    ir <- .dist_spatial_inner(sweep_fn, X_lam, X_sig, map, Q, tau, rho, kind, ldQ,
                              beta_lam0, beta_sig0, rep(0, n_sp), max_iter, tol, verbose)
    if (!isTRUE(ir$ok)) next
    modes[k, ] <- c(ir$beta_lam, ir$beta_sig); fields[k, ] <- ir$z
    cov_blocks[[k]] <- ir$cov_beta; logm[k] <- ir$log_marginal; bmax[k] <- ir$boundary_max
  }
  ok <- is.finite(logm)
  if (!any(ok)) stop("distance() areal spatial fit produced no usable grid point.", call. = FALSE)
  w <- tulpa:::.nl_normalise_weights_safe(logm, "tau_grid / data")
  w[!ok] <- 0; w <- w / sum(w)

  beta_mean <- as.numeric(crossprod(w, modes)); beta_mean[!is.finite(beta_mean)] <- 0
  field_mean <- as.numeric(crossprod(w, ifelse(is.finite(fields), fields, 0)))
  nm <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
          paste0("sigma_",  model$process_info[[2]]$coef_names))
  # grid-integrated coef cov (law of total covariance over the grid).
  mbar <- beta_mean
  V <- matrix(0, p_beta, p_beta)
  for (k in which(ok)) {
    Vk <- cov_blocks[[k]]; if (is.null(Vk)) next
    dk <- modes[k, ] - mbar
    V <- V + w[k] * (Vk + outer(dk, dk))
  }
  dimnames(V) <- list(nm, nm)
  names(beta_mean) <- nm

  raw <- list(
    mixture = if (is_nb) "negbin" else "poisson",
    beta_lambda = beta_mean[seq_len(p_lam)], beta_sigma = beta_mean[p_lam + seq_len(p_sig)],
    log_r = NA_real_, r = NA_real_, vcov = V,
    log_lik = sum(w * ifelse(ok, logm, 0)), converged = TRUE,
    key = model$key, transect = model$transect, hazard = FALSE, K_max = K_max)
  fit <- build_distance_fit(raw, model)
  fit$method <- "nested_laplace"
  fit$spatial_field <- field_mean
  if (any(bmax > 1e-4, na.rm = TRUE))
    warning(sprintf("Max posterior weight on N = K_max is %.2e at one or more grid points; raise K_max.",
                    max(bmax, na.rm = TRUE)), call. = FALSE)
  fit
}
