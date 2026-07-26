# =============================================================================
# occu_cover_nested.R - v3 nested-Laplace spatial path for occu_cover().
#
# Profile the latent field z[1..n_cells] OUT of the outer optimisation. For
# each candidate (betas, alpha, log_sigma), find z* by inner Newton on the
# conditional log-posterior, then evaluate the Laplace approximation of the
# marginal log-likelihood
#
#   log p(y | betas, alpha, sigma) approx
#     log p(y | betas, alpha, z*) + log p(z* | sigma)
#       + (n_eff / 2) * log(2 pi) - 0.5 * log |H_z(z*)|
#
# Outer optim over (betas, alpha, log_sigma) is then ~10-dimensional, so
# BFGS handles it without the (alpha, sigma) ridge that v2's joint Laplace
# slid down.
#
# v3 vs v2:
#   v2: joint Laplace MAP on c(betas, z, alpha, log_sigma) -- ridge issue.
#   v3: nested Laplace -- inner Newton on z given outer; outer BFGS on small
#       hyperparameter+coef block. Same family API (`method = "nested_laplace"`).
#
# Limitations still pending v4:
#   - Lognormal positive arm only (beta arm pending; needs the cleaner
#     d log f_pos / d eta closed form). occu_cover("beta") under
#     nested_laplace currently errors with a pointer.
#   - ICAR (besag) field only; BYM2 with free rho mixing pending.
#   - Outer-grid integration of (alpha, sigma) for proper posterior summaries
#     pending; for now the outer Laplace gives the point + observed-Fisher SE.
# =============================================================================


# ---------------------------------------------------------------------------
# Time-throttled Stan-style progress reporter. Returns a closure called once
# per outer BFGS call; emits a line when more than `throttle` seconds have
# passed since the last emit. Text-only (no cursor control), so it works in
# batch logs, file redirects, and CI as well as in interactive R.
# Motivation: collaborators cancel a long fit when nothing prints for
# minutes -- a 5 s heartbeat eliminates that failure mode.
# ---------------------------------------------------------------------------
.tobs_progress_reporter <- function(label, throttle = 5, max_calls = NULL) {
  start <- Sys.time()
  last  <- start
  iter  <- 0L
  best  <- Inf
  fmt_secs <- function(s) {
    if (!is.finite(s) || s < 0) return("?")
    if (s < 90)        sprintf("%.0fs",   s)
    else if (s < 5400) sprintf("%.1fmin", s / 60)
    else                sprintf("%.1fh",  s / 3600)
  }
  function(value = NA_real_, extra = "") {
    iter <<- iter + 1L
    if (is.finite(value) && value < best) best <<- value
    now <- Sys.time()
    if (iter == 1L ||
        as.numeric(difftime(now, last, units = "secs")) >= throttle) {
      elapsed <- as.numeric(difftime(now, start, units = "secs"))
      rate    <- if (iter > 1L) elapsed / iter else NA_real_
      rate_str <- if (is.finite(rate)) sprintf("  %.2fs/call", rate) else ""
      eta_str <- ""
      if (!is.null(max_calls) && is.finite(rate)) {
        # Upper bound: assume the optim runs the full max-call ceiling
        # (calls/iter * max.iter). Ceiling shrinks each line; the fit
        # almost always converges before it.
        remaining <- max(0, max_calls - iter)
        eta_str <- sprintf("  max ETA %s", fmt_secs(remaining * rate))
      }
      msg <- sprintf("[%s] call %d  nlp = %.4f  best = %.4f  elapsed %s%s%s",
                      label, iter, value, best, fmt_secs(elapsed),
                      rate_str, eta_str)
      if (nzchar(extra)) msg <- paste0(msg, "  ", extra)
      cat(msg, "\n", sep = "")
      utils::flush.console()
      last <<- now
    }
    invisible(NULL)
  }
}


# ---------------------------------------------------------------------------
# Conditional log-density of z given (betas, alpha, sigma), unnormalised.
# Returns -log p_unnorm(z | other), the negative log NLP for inner Newton.
# Vectorised across cells; per-cell data contributions are independent.
# ---------------------------------------------------------------------------
.occu_cover_z_nlp <- function(z, model, Q, scale_q,
                               beta_psi, beta_p_vec, beta_pos_vec,
                               log_disp, alpha, log_sigma) {
  cl <- .tobs_clamp_eta
  n_cells    <- model$n_sites
  max_visits <- model$max_visits

  # Psi: cell-level + field.
  eta_psi <- as.numeric(model$X_occ %*% beta_psi) + z
  psi     <- stats::plogis(cl(eta_psi))

  # Detection: site-level + visit-level (no field).
  bp_site  <- beta_p_vec[seq_len(ncol(model$X_det_site))]
  bp_visit <- if (!is.null(model$X_det_visit))
                 beta_p_vec[ncol(model$X_det_site) + seq_len(ncol(model$X_det_visit))]
              else numeric(0)
  eta_p_site <- as.numeric(model$X_det_site %*% bp_site)
  p_mat <- matrix(eta_p_site, n_cells, max_visits)
  if (length(bp_visit)) {
    p_mat <- p_mat + matrix(as.numeric(model$X_det_visit %*% bp_visit),
                             n_cells, max_visits, byrow = TRUE)
  }
  p_mat <- stats::plogis(cl(p_mat))

  # Cover-arm linear predictor (site + visit + alpha * field).
  bpos_site  <- beta_pos_vec[seq_len(ncol(model$X_pos_site))]
  bpos_visit <- if (!is.null(model$X_pos_visit))
                   beta_pos_vec[ncol(model$X_pos_site) + seq_len(ncol(model$X_pos_visit))]
                else numeric(0)
  eta_pos_site <- as.numeric(model$X_pos_site %*% bpos_site)
  ep_mat <- matrix(eta_pos_site, n_cells, max_visits)
  if (length(bpos_visit)) {
    ep_mat <- ep_mat + matrix(as.numeric(model$X_pos_visit %*% bpos_visit),
                               n_cells, max_visits, byrow = TRUE)
  }
  ep_mat <- ep_mat + matrix(alpha * z, n_cells, max_visits)

  valid <- model$valid
  y     <- model$y
  y_pos <- model$y_pos

  log_p   <- ifelse(valid, log(p_mat),     0)
  log_1mp <- ifelse(valid, log(1 - p_mat), 0)

  # Positive-arm log-density (lognormal).
  sigma_pos <- exp(log_disp)
  pos_mask <- valid & (y == 1L)
  log_f_pos <- matrix(0, n_cells, max_visits)
  log_f_pos[pos_mask] <- (-log(y_pos[pos_mask]) - log(sigma_pos)
                          - 0.5 * log(2 * pi)
                          - 0.5 * ((log(y_pos[pos_mask]) - ep_mat[pos_mask]) /
                                     sigma_pos)^2)

  log_h <- ifelse(valid,
                  ifelse(y == 1L, log_p + log_f_pos, log_1mp),
                  0)

  any_det <- rowSums(y * valid, na.rm = FALSE) > 0
  log_psi   <- log(pmax(psi, 1e-300))
  log_1mpsi <- log(pmax(1 - psi, 1e-300))

  det_ll <- log_psi + rowSums(log_h)
  ln_a <- log_psi   + rowSums(log_1mp)
  ln_b <- log_1mpsi
  nodet_ll <- .tobs_logsumexp2(ln_a, ln_b)

  ll <- sum(ifelse(any_det, det_ll, nodet_ll))

  inv_sig2 <- exp(-2 * log_sigma)
  z_prior  <- 0.5 * inv_sig2 * scale_q * as.numeric(crossprod(z, Q %*% z))

  -ll + z_prior
}


# ---------------------------------------------------------------------------
# Analytical gradient of -log p(z | other) wrt z.
# Used by the inner BFGS so we do not pay finite-difference cost per cell.
# ---------------------------------------------------------------------------
.occu_cover_z_grad <- function(z, model, Q, scale_q,
                                beta_psi, beta_p_vec, beta_pos_vec,
                                log_disp, alpha, log_sigma) {
  cl <- .tobs_clamp_eta
  n_cells    <- model$n_sites
  max_visits <- model$max_visits

  eta_psi <- as.numeric(model$X_occ %*% beta_psi) + z
  psi     <- stats::plogis(cl(eta_psi))

  bp_site  <- beta_p_vec[seq_len(ncol(model$X_det_site))]
  bp_visit <- if (!is.null(model$X_det_visit))
                 beta_p_vec[ncol(model$X_det_site) + seq_len(ncol(model$X_det_visit))]
              else numeric(0)
  eta_p_site <- as.numeric(model$X_det_site %*% bp_site)
  p_mat <- matrix(eta_p_site, n_cells, max_visits)
  if (length(bp_visit)) {
    p_mat <- p_mat + matrix(as.numeric(model$X_det_visit %*% bp_visit),
                             n_cells, max_visits, byrow = TRUE)
  }
  p_mat <- stats::plogis(cl(p_mat))

  bpos_site  <- beta_pos_vec[seq_len(ncol(model$X_pos_site))]
  bpos_visit <- if (!is.null(model$X_pos_visit))
                   beta_pos_vec[ncol(model$X_pos_site) + seq_len(ncol(model$X_pos_visit))]
                else numeric(0)
  eta_pos_site <- as.numeric(model$X_pos_site %*% bpos_site)
  ep_mat <- matrix(eta_pos_site, n_cells, max_visits)
  if (length(bpos_visit)) {
    ep_mat <- ep_mat + matrix(as.numeric(model$X_pos_visit %*% bpos_visit),
                               n_cells, max_visits, byrow = TRUE)
  }
  ep_mat <- ep_mat + matrix(alpha * z, n_cells, max_visits)

  valid <- model$valid
  y     <- model$y
  y_pos <- model$y_pos
  sigma_pos <- exp(log_disp)

  any_det <- rowSums(y * valid, na.rm = FALSE) > 0

  # Cover-arm: d log f_pos / d eta_pos = (log y - eta_pos) / sigma^2 at
  # detected visits (lognormal). Multiply by alpha for d/dz_i then sum across
  # visits j (only detected ones contribute).
  pos_grad_mat <- matrix(0, n_cells, max_visits)
  pos_mask <- valid & (y == 1L)
  pos_grad_mat[pos_mask] <- (log(y_pos[pos_mask]) - ep_mat[pos_mask]) /
                              (sigma_pos^2)
  cover_grad <- alpha * rowSums(pos_grad_mat)

  # Psi-arm gradient varies by case.
  psi_grad <- numeric(n_cells)

  # any_det case: d log L / dz = (1 - psi) + cover_grad
  psi_grad[any_det] <- (1 - psi[any_det])

  # no_det case: d log L / dz = -psi (1 - psi) (1 - P_0) / L
  if (any(!any_det)) {
    log_1mp <- ifelse(valid, log(1 - p_mat), 0)
    log_P0  <- rowSums(log_1mp)
    P0      <- exp(log_P0)
    L_nodet <- psi * P0 + (1 - psi)
    psi_grad[!any_det] <- -psi[!any_det] * (1 - psi[!any_det]) *
                            (1 - P0[!any_det]) / pmax(L_nodet[!any_det], 1e-300)
  }

  data_grad <- psi_grad + cover_grad

  inv_sig2  <- exp(-2 * log_sigma)
  prior_grad <- inv_sig2 * scale_q * as.numeric(Q %*% z)

  -data_grad + prior_grad
}


# ---------------------------------------------------------------------------
# Diagonal data-Hessian wrt z, evaluated analytically for the any_det case
# (lognormal) and via 3-point FD for the no_det case (cleaner than chain
# rule on `log[psi P_0 + 1 - psi]`).
# Returns the diagonal of `-d^2 log L_data / dz^2` per cell.
# ---------------------------------------------------------------------------
.occu_cover_data_hess_diag <- function(z, model, Q, scale_q,
                                        beta_psi, beta_p_vec, beta_pos_vec,
                                        log_disp, alpha, log_sigma) {
  cl <- .tobs_clamp_eta
  n_cells    <- model$n_sites
  max_visits <- model$max_visits

  eta_psi <- as.numeric(model$X_occ %*% beta_psi) + z
  psi     <- stats::plogis(cl(eta_psi))
  any_det <- rowSums(model$y * model$valid, na.rm = FALSE) > 0

  H_diag <- numeric(n_cells)

  # any_det case: d^2/dz^2 log L = d/dz (1-psi) + alpha^2 * sum d^2 log f / d eta^2
  #   = -psi*(1-psi)  +  alpha^2 * n_pos_i * (-1/sigma_pos^2)
  # n_pos_i = number of detected visits in cell i.
  sigma_pos <- exp(log_disp)
  n_pos <- rowSums(model$valid & (model$y == 1L))
  H_diag[any_det] <- psi[any_det] * (1 - psi[any_det]) +
                       (alpha^2) * n_pos[any_det] / (sigma_pos^2)

  # no_det case: FD on the NLP wrt z_i. O(n_nodet) cell-level evaluations.
  # Each evaluation re-walks the whole likelihood; this is the v3 bottleneck.
  # For small/medium n_cells (<= ~1000) it is still much cheaper than the
  # joint optim's BFGS line search.
  if (any(!any_det)) {
    delta <- 1e-4
    base_nlp <- .occu_cover_z_nlp(z, model, Q, scale_q, beta_psi,
                                   beta_p_vec, beta_pos_vec,
                                   log_disp, alpha, log_sigma)
    # Strip out the prior term (it has its own analytic Hessian -- scale_q
    # * Q -- added at the caller).
    inv_sig2 <- exp(-2 * log_sigma)
    base_prior <- 0.5 * inv_sig2 * scale_q * as.numeric(crossprod(z, Q %*% z))
    base_loss  <- base_nlp - base_prior

    nodet_idx <- which(!any_det)
    for (k in nodet_idx) {
      z_p <- z; z_p[k] <- z[k] + delta
      z_m <- z; z_m[k] <- z[k] - delta
      nlp_p <- .occu_cover_z_nlp(z_p, model, Q, scale_q, beta_psi,
                                   beta_p_vec, beta_pos_vec,
                                   log_disp, alpha, log_sigma)
      nlp_m <- .occu_cover_z_nlp(z_m, model, Q, scale_q, beta_psi,
                                   beta_p_vec, beta_pos_vec,
                                   log_disp, alpha, log_sigma)
      pr_p <- 0.5 * inv_sig2 * scale_q * as.numeric(crossprod(z_p, Q %*% z_p))
      pr_m <- 0.5 * inv_sig2 * scale_q * as.numeric(crossprod(z_m, Q %*% z_m))
      loss_p <- nlp_p - pr_p
      loss_m <- nlp_m - pr_m
      # FD second derivative on -log L_data only (positive at the mode).
      H_diag[k] <- max(0, (loss_p + loss_m - 2 * base_loss) / delta^2)
    }
  }

  H_diag
}


# ---------------------------------------------------------------------------
# Inner Newton: find z* given outer params. Returns z mode, value at mode,
# and the diagonal data-Hessian (so the outer step can add the prior).
# ---------------------------------------------------------------------------
.occu_cover_inner_newton <- function(outer_par, model, Q, scale_q, z_warm,
                                      max_inner = 30L, tol = 1e-6) {
  p_psi <- model$process_info[[1L]]$p
  p_p   <- model$process_info[[2L]]$p
  p_pos <- model$process_info[[3L]]$p

  off <- 0L
  beta_psi    <- outer_par[off + seq_len(p_psi)]; off <- off + p_psi
  beta_p_vec  <- outer_par[off + seq_len(p_p)];   off <- off + p_p
  beta_pos_vec<- outer_par[off + seq_len(p_pos)]; off <- off + p_pos
  log_disp    <- outer_par[off + 1L];             off <- off + 1L
  alpha       <- outer_par[off + 1L];             off <- off + 1L
  log_sigma   <- outer_par[off + 1L]

  z <- if (is.null(z_warm)) numeric(model$n_sites) else z_warm
  inv_sig2 <- exp(-2 * log_sigma)
  # Prior Hessian piece (constant for the inner Newton at fixed sigma).
  H_prior <- inv_sig2 * scale_q * Q
  kappa_sum <- 1e4   # soft sum-to-zero penalty (matches v2)

  for (iter in seq_len(max_inner)) {
    g <- .occu_cover_z_grad(z, model, Q, scale_q, beta_psi, beta_p_vec,
                             beta_pos_vec, log_disp, alpha, log_sigma)
    # Add gradient of soft sum-to-zero.
    g <- g + kappa_sum * sum(z)

    H_diag_data <- .occu_cover_data_hess_diag(z, model, Q, scale_q,
                                                beta_psi, beta_p_vec,
                                                beta_pos_vec, log_disp,
                                                alpha, log_sigma)
    H <- diag(H_diag_data) + H_prior + kappa_sum   # sum-to-zero Hessian = kappa*J
    delta <- tryCatch(solve(H, g), error = function(e) NULL)
    if (is.null(delta)) break

    # Damped step if it overshoots.
    step <- 1.0
    z_new <- z - step * delta
    if (max(abs(delta)) < tol) {
      z <- z_new
      break
    }
    z <- z_new
  }

  H_diag_data <- .occu_cover_data_hess_diag(z, model, Q, scale_q,
                                              beta_psi, beta_p_vec,
                                              beta_pos_vec, log_disp,
                                              alpha, log_sigma)

  list(z = z, H_data_diag = H_diag_data,
       beta_psi = beta_psi, beta_p_vec = beta_p_vec,
       beta_pos_vec = beta_pos_vec, log_disp = log_disp,
       alpha = alpha, log_sigma = log_sigma)
}


# ---------------------------------------------------------------------------
# Outer objective: negative marginal log-likelihood at given outer params.
# Side effect: write the z mode to `z_warm_env$z` so the next outer call
# starts warm from the previous mode.
# ---------------------------------------------------------------------------
.occu_cover_outer_nlp <- function(outer_par, model, Q, scale_q, z_warm_env) {
  inner <- .occu_cover_inner_newton(outer_par, model, Q, scale_q,
                                      z_warm = z_warm_env$z)
  z_warm_env$z <- inner$z

  # Evaluate NLP at inner mode (data NLP + prior).
  nlp_at_mode <- .occu_cover_z_nlp(inner$z, model, Q, scale_q,
                                     inner$beta_psi, inner$beta_p_vec,
                                     inner$beta_pos_vec, inner$log_disp,
                                     inner$alpha, inner$log_sigma)

  # Laplace correction: -log marginal approx
  #   = NLP_at_mode + 0.5 * log|H_z| - (n_eff/2) * log(2 pi)
  # H_z = diag(H_data_diag) + (1/sigma^2) * scale_q * Q + kappa_sum * J(soft)
  inv_sig2  <- exp(-2 * inner$log_sigma)
  n_cells   <- model$n_sites
  n_eff     <- n_cells - 1L
  kappa_sum <- 1e4
  H_z <- diag(inner$H_data_diag) + inv_sig2 * scale_q * Q + kappa_sum

  log_det_H <- tryCatch(as.numeric(determinant(H_z, logarithm = TRUE)$modulus),
                         error = function(e) NA_real_)
  if (!is.finite(log_det_H)) return(.Machine$double.xmax / 2)

  # Sigma log-determinant correction from the field prior.
  #   log p(z|sigma) at the mode includes -0.5 * log|Sigma_prior|
  #   = -n_eff * log_sigma + 0.5 * log|Q_pseudo_inv|  (const in sigma)
  # The "-n_eff * log_sigma" is the only sigma-dependent normaliser; bake
  # it into the NLP for the outer-step sigma score to be correct.
  sigma_logdet_term <- n_eff * inner$log_sigma

  # Outer NLP. Drop -0.5 * (n_eff) * log(2*pi) (constant).
  nlp_at_mode + sigma_logdet_term + 0.5 * log_det_H
}


# ---------------------------------------------------------------------------
# Public fitter -- the v3 replacement for v2's joint Laplace.
# ---------------------------------------------------------------------------
.tobs_fit_occu_cover_nested <- function(model, adj,
                                          priors    = NULL,
                                          max.iter  = 200L,
                                          tol       = 1e-6,
                                          verbose   = TRUE,
                                          sigma.beta = 5,
                                          ...) {
  if (!identical(model$positive, "lognormal")) {
    stop("occu_cover() v3 nested-Laplace currently supports positive = ",
         "\"lognormal\" only. Beta positive arm coming in v4 (needs the ",
         "cleaner closed-form d log f / d eta for the inner Newton).",
         call. = FALSE)
  }

  pi_list <- model$process_info
  p_psi   <- pi_list[[1L]]$p
  p_p     <- pi_list[[2L]]$p
  p_pos   <- pi_list[[3L]]$p
  n_cells <- model$n_sites

  if (nrow(adj) != n_cells) {
    stop(sprintf("Spatial graph has %d nodes but data has %d cells.",
                 nrow(adj), n_cells), call. = FALSE)
  }

  Q       <- .occu_cover_icar_Q(adj)
  scale_q <- .occu_cover_icar_scale(adj)

  n_outer <- p_psi + p_p + p_pos + 1L + 2L   # +1 log_disp, +2 (alpha, log_sigma)
  par_names <- c(
    paste0("psi_", pi_list[[1L]]$coef_names),
    paste0("p_",   pi_list[[2L]]$coef_names),
    paste0("pos_", pi_list[[3L]]$coef_names),
    "log_sigma_pos", "alpha", "log_sigma"
  )

  start <- numeric(n_outer)
  any_det <- rowSums(model$y * model$valid) > 0
  det_rate <- max(mean(any_det), 1e-3)
  start[1L] <- stats::qlogis(min(max(det_rate, 1e-3), 1 - 1e-3))

  pos_vals <- model$y_pos[model$valid & model$y == 1L]
  pos_vals <- pos_vals[is.finite(pos_vals)]
  disp_idx <- p_psi + p_p + p_pos + 1L
  pos_int_idx <- p_psi + p_p + 1L
  if (length(pos_vals) > 0L) {
    start[pos_int_idx] <- mean(log(pos_vals))
    start[disp_idx]    <- log(stats::sd(log(pos_vals)) + 0.1)
  } else {
    start[disp_idx] <- log(0.4)
  }
  start[n_outer - 1L] <- 0    # alpha
  start[n_outer]      <- 0    # log_sigma -> sigma = 1

  # Gaussian-prior precision on fixed effects only.
  pmean <- numeric(n_outer)
  pprec <- numeric(n_outer)
  if (isTRUE(is.null(priors)) || !isFALSE(priors)) {
    beta_idx <- seq_len(p_psi + p_p + p_pos)
    pprec[beta_idx] <- 1 / (sigma.beta^2)
  }
  # Weakly-informative anchors on hypers (much wider than v2's: v3 doesn't
  # need them tight, the (alpha, sigma) ridge is broken by profiling).
  pprec[n_outer - 1L] <- 1 / (5^2)   # alpha
  pprec[n_outer]      <- 1 / (2^2)   # log_sigma

  # Warm-storage env for the inner z mode across outer calls.
  warm_env <- new.env(parent = emptyenv())
  warm_env$z <- numeric(n_cells)

  # Progress heartbeat (default on; silenced by verbose = FALSE).
  # max_calls ceiling = max.iter * (FD-gradient calls + ~2 line-search calls).
  # Gradient takes n_outer + 1 calls under stats::optim's default FD;
  # line search is typically 1-4 calls per iter. Use n_outer + 3 as a tight
  # upper bound on calls / iter.
  calls_per_iter_est <- n_outer + 3L
  max_calls_est <- as.integer(max.iter) * calls_per_iter_est
  report <- if (isTRUE(verbose) || is.null(verbose) || verbose != FALSE)
              .tobs_progress_reporter("occu_cover v3", throttle = 5,
                                       max_calls = max_calls_est)
            else
              function(...) invisible(NULL)

  outer_nlp_with_prior <- function(par) {
    val <- .occu_cover_outer_nlp(par, model, Q, scale_q, warm_env) +
             0.5 * sum(pprec * (par - pmean)^2)
    report(val)
    val
  }

  opt <- stats::optim(start, outer_nlp_with_prior,
                       method = "BFGS", hessian = TRUE,
                       control = list(maxit = max.iter, reltol = tol,
                                      trace = if (isTRUE(verbose)) 1L else 0L))

  V <- tryCatch(solve(opt$hessian), error = function(e) NULL)
  if (is.null(V)) {
    warning("occu_cover nested: outer Hessian not invertible; SEs unreliable.",
            call. = FALSE)
    V <- matrix(NA_real_, n_outer, n_outer)
  }
  se <- sqrt(pmax(diag(V), 0))

  means <- opt$par
  names(means) <- par_names
  names(se)    <- par_names
  dimnames(V)  <- list(par_names, par_names)

  # Run one more inner pass at the converged outer params to get the final z
  # and the inner Hessian for field SE.
  inner_final <- .occu_cover_inner_newton(opt$par, model, Q, scale_q,
                                            z_warm = warm_env$z,
                                            max_inner = 50L)
  z_final  <- inner_final$z - mean(inner_final$z)  # enforce sum-to-zero exactly
  inv_sig2 <- exp(-2 * inner_final$log_sigma)
  H_z_final <- diag(inner_final$H_data_diag) + inv_sig2 * scale_q * Q + 1e4
  V_z       <- tryCatch(solve(H_z_final), error = function(e) NULL)
  z_sd      <- if (is.null(V_z)) rep(NA_real_, n_cells) else sqrt(pmax(diag(V_z), 0))

  field_table <- data.frame(
    cell    = seq_len(n_cells),
    z_mean  = z_final,
    z_sd    = z_sd,
    z_lower = z_final - 1.96 * z_sd,
    z_upper = z_final + 1.96 * z_sd
  )

  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = se,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = n_outer,
    log_prob     = rep(-opt$value, n_draws),
    log_lik      = -opt$value,
    N            = sum(model$valid)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = list(type = "icar", graph = adj,
                        sigma_mean = exp(means["log_sigma"]),
                        alpha_mean = means["alpha"]),
    spatial_field = z_final,
    field_table  = field_table,
    method       = "nested_laplace",
    positive     = model$positive,
    convergence  = list(converged = opt$convergence == 0L,
                        n_iter    = opt$counts[1L])
  )), class = c("tobs_fit", "tulpa_fit"))
}
