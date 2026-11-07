# R oracles for the compiled targets and kernels: each is the byte-exact (or
# FD-validated) R mirror of a C++ log-posterior / likelihood that the tests
# compare against. Test-only, so they live beside their only consumers
# rather than in the package namespace.

#' Third derivative of Gaussian log-likelihood (identity link)
#'
#' Always zero -- the SLA correction reduces to the Laplace approximation
#' for Gaussian likelihoods. Required regression test.
#'
#' @keywords internal
.l3_gaussian_identity <- function(eta) {
  rep(0, length(eta))
}

#' Third derivative of Poisson log-likelihood wrt log(lambda)
#'
#' For y_i ~ Poisson(exp(eta_i)),
#'   d^3 log p / d eta^3 = -lambda_i = -exp(eta_i)
#'
#' @keywords internal
.l3_poisson_log <- function(eta) {
  -exp(eta)
}

#' R-side single-season occupancy log-likelihood
#'
#' Computes the marginal observation log-likelihood
#'   `log P(y | beta) = Sum_i log[ psi_i * P(y_i | z=1) + (1 - psi_i) I(no detections) ]`
#' with psi_i = plogis(X_occ_i beta_occ), p_i = plogis(X_det_i beta_det).
#'
#' Used to cross-check `.sla_gamma_fd()` against the analytical
#' `.sla_gamma_diag()` on single-season fits.
#'
#' @keywords internal
.loglik_occu_single <- function(beta, model) {
  pi_list <- model$process_info
  p_occ <- pi_list[[1]]$p
  p_det <- pi_list[[2]]$p
  beta_occ <- beta[seq_len(p_occ)]
  beta_det <- beta[p_occ + seq_len(p_det)]

  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  y     <- model$y
  n_sites <- model$n_sites

  psi <- plogis(as.numeric(X_occ %*% beta_occ))
  p   <- plogis(as.numeric(X_det %*% beta_det))

  ll <- 0
  for (i in seq_len(n_sites)) {
    valid <- y[i, ] >= 0
    n_v <- sum(valid)
    if (n_v == 0L) next
    n_d <- sum(y[i, valid] == 1L)
    if (n_d > 0L) {
      # any-detection: z=1 forced
      ll <- ll + log(psi[i]) + n_d * log(p[i]) + (n_v - n_d) * log1p(-p[i])
    } else {
      # all zeros: marginalise over z
      q_i <- (1 - p[i])^n_v
      ll <- ll + log(psi[i] * q_i + (1 - psi[i]))
    }
  }
  ll
}

# Full-vector unconstrained joint log-posterior and its gradient for the
# spatial-factor community occu_cover model. `theta` packs
#   c(par_inner, chol_occ, chol_p, chol_pos, log_tau_w[, logit_h_w])
# with par_inner = c(mu, {b_s}, L-block, [Lpos], W, log_disp) (constrained L-block
# when `constrain`). The trailing logit_h_w block is present iff the field carries
# a hyperparameter (proper-CAR rho / BYM2 phi); h_k = plogis(logit_h_k) in (0, 1).
# Returns list(lp, grad) over the whole vector.
.ms_ocs_joint_logpost <- function(model, theta, priors = .ms_ocs_nuts_priors(),
                                  constrain = FALSE, sigma.beta = 5, sd_L = 1.0,
                                  grad = TRUE) {
  spec <- model$field_spec %||%
    .ms_ocs_field_spec(model$adj, model$field_type %||% "icar")
  has_hyper <- isTRUE(spec$has_hyper)
  d   <- .ms_ocs_dims(model)
  lay <- .ms_ocs_nuts_layout(d, constrain, has_hyper)
  S   <- d$S; K <- d$K
  rank_k <- spec$rank                                  # N (proper) / N - 1 (icar)

  par_inner <- theta[lay$inner]
  C <- list(occ = .ms_ocs_chol_unpack(theta[lay$chol_occ], d$P_occ),
            p   = .ms_ocs_chol_unpack(theta[lay$chol_p],   d$P_p),
            pos = .ms_ocs_chol_unpack(theta[lay$chol_pos], d$P_pos))
  log_tau <- theta[lay$log_tau]
  tau_w   <- exp(log_tau)
  h_w <- if (has_hyper) stats::plogis(theta[lay$logit_h]) else NULL

  # Per-factor field structure R_k(h_k) (the fixed ICAR Q for every factor when
  # the field has no hyperparameter); the inner objective reads it from field_R.
  Rk <- .ms_ocs_build_field_R(model, h_w)
  model_eff <- model; model_eff$field_R <- Rk

  Si <- lapply(C, function(Cm) chol2inv(t(Cm)))        # Sigma_arm^{-1}
  Sinv <- matrix(0, d$P, d$P)
  Sinv[d$occ_idx, d$occ_idx] <- Si$occ
  Sinv[d$p_idx,   d$p_idx]   <- Si$p
  Sinv[d$pos_idx, d$pos_idx] <- Si$pos
  Pmu      <- diag(1 / sigma.beta^2, d$P)
  inv_sdL2 <- 1 / sd_L^2

  objective <- if (constrain) .ms_ocs_penll_grad_c else .ms_ocs_penll_grad
  res <- objective(model_eff, par_inner, Sinv, Pmu, inv_sdL2, tau_w, grad = grad)

  # ---- (B) MVN(b) log-determinant normalisers + (C) GMRF tau normalisers ----
  logdet <- c(occ = 2 * sum(log(diag(C$occ))),
              p   = 2 * sum(log(diag(C$p))),
              pos = 2 * sum(log(diag(C$pos))))
  lp <- res$ll - 0.5 * S * sum(logdet) + 0.5 * rank_k * sum(log_tau)
  # The field-structure log|R(h_k)| enters only for a hyperparameterised field;
  # for icar it is the constant pseudo-determinant and drops.
  if (has_hyper) {
    lp <- lp + 0.5 * sum(vapply(h_w, function(h)
      .ms_ocs_field_logdetR(spec, h), 0))
  }

  # ---- (D) hyperpriors on the sampled coordinates ----
  pr_occ <- .ms_ocs_chol_logprior(theta[lay$chol_occ], d$P_occ, priors)
  pr_p   <- .ms_ocs_chol_logprior(theta[lay$chol_p],   d$P_p,   priors)
  pr_pos <- .ms_ocs_chol_logprior(theta[lay$chol_pos], d$P_pos, priors)
  pr_tau <- -0.5 * sum(((log_tau - priors$log_tau_mean) / priors$log_tau_sd)^2)
  ld     <- par_inner[lay$n_inner]                     # log_disp = last inner coord
  pr_ld  <- -0.5 * ((ld - priors$log_disp_mean) / priors$log_disp_sd)^2
  lp <- lp + pr_occ$lp + pr_p$lp + pr_pos$lp + pr_tau + pr_ld
  if (has_hyper) {
    lp <- lp - 0.5 * sum(((theta[lay$logit_h] - priors$logit_h_mean) /
                            priors$logit_h_sd)^2)
  }

  if (!grad) return(list(lp = lp))

  g <- numeric(lay$total)
  g_inner <- res$grad
  g_inner[lay$n_inner] <- g_inner[lay$n_inner] -
    (ld - priors$log_disp_mean) / priors$log_disp_sd^2  # log_disp prior
  g[lay$inner] <- g_inner

  # ---- Cholesky-block gradient, per arm ----
  # b-block second moment M_arm = sum_s b_{s,arm} b_{s,arm}'. The packed inner b
  # for species s sits at offset P + (s-1)P; extract each arm's sub-vector.
  up_b <- .ms_ocs_inner_b(par_inner, d, constrain)     # list length S of length-P
  arm_idx <- list(occ = d$occ_idx, p = d$p_idx, pos = d$pos_idx)
  chol_slot <- list(occ = lay$chol_occ, p = lay$chol_p, pos = lay$chol_pos)
  Pa <- list(occ = d$P_occ, p = d$P_p, pos = d$P_pos)
  for (arm in c("occ", "p", "pos")) {
    ai <- arm_idx[[arm]]
    M  <- matrix(0, Pa[[arm]], Pa[[arm]])
    for (s in seq_len(S)) { bb <- up_b[[s]][ai]; M <- M + outer(bb, bb) }
    g[chol_slot[[arm]]] <- .ms_ocs_chol_block_grad(
      C[[arm]], M, S, theta[chol_slot[[arm]]], Pa[[arm]], priors)
  }

  # ---- log tau_w (and, for a hyper field, logit_h_w) gradient ----
  # d/d log_tau_k of [ -0.5 tau_k W_k' R_k W_k + 0.5 rank log tau_k - prior ];
  # d/d logit_h_k chains the field hyperparameter through h_k = plogis().
  W  <- .ms_ocs_inner_W(par_inner, d, constrain)       # N x K
  for (k in seq_len(K)) {
    quad <- as.numeric(W[, k] %*% (Rk[[k]] %*% W[, k]))
    g[lay$log_tau[k]] <- -0.5 * tau_w[k] * quad + 0.5 * rank_k -
      (log_tau[k] - priors$log_tau_mean) / priors$log_tau_sd^2
    if (has_hyper) {
      ht <- .ms_ocs_field_hyper_terms(spec, h_w[k], W[, k])
      dlp_dh <- -0.5 * tau_w[k] * ht$dquad + 0.5 * ht$dlogdet
      dh     <- h_w[k] * (1 - h_w[k])                  # plogis derivative
      g[lay$logit_h[k]] <- dlp_dh * dh -
        (theta[lay$logit_h[k]] - priors$logit_h_mean) / priors$logit_h_sd^2
    }
  }

  list(lp = lp, grad = g)
}

# Per-species community N-mixture oracle for the shared RE integrator. Each
# species' marginal (nmix_site_marginal) is built once; the oracle
# evaluates it at the full coefficients coef = (mu + b_s), mapping the RE
# deviation through the abundance / detection designs. The b-space curvature is
# the design-sandwiched per-site observed-information block (with the latent-N
# coupling) -- the covariate generalization of the intercept-only assembly.
.nmix_re_oracle <- function(y, site_idx, species_idx, X_lambda, X_p,
                            n_sites, n_species, p_lam, p_p, K_max,
                            headroom = NULL) {
  cl <- .tobs_clamp_eta
  d  <- p_lam + p_p
  rows_by_sp <- split(seq_len(length(y)), species_idx)
  marg <- lapply(seq_len(n_species), function(s) {
    sel <- rows_by_sp[[as.character(s)]]
    if (is.null(sel)) sel <- integer(0)
    nmix_site_marginal(
      y = y[sel], site_idx = site_idx[sel],
      X_lambda = X_lambda, X_p = X_p[sel, , drop = FALSE],
      mixture = "P", K_max = K_max, headroom = headroom)
  })

  eta_of <- function(m, coef) {
    list(lambda = cl(as.numeric(m$X_lambda %*% coef[seq_len(p_lam)])),
         p      = cl(as.numeric(m$X_p %*% coef[p_lam + seq_len(p_p)])))
  }

  list(
    grad_hess = function(s, coef) {
      m  <- marg[[s]]
      e  <- eta_of(m, coef)
      ev <- m$eval(e$lambda, e$p)
      grad <- c(as.numeric(crossprod(m$X_lambda, ev$grad_eta_lambda)),
                as.numeric(crossprod(m$X_p, ev$grad_eta_p)))
      # negH: marginal observed info (with the Var[N|y] abundance/detection
      # coupling), used by the Laplace marginal. fisher: complete-data Fisher
      # (block-diagonal, PSD), supplied for the mode-find Newton -- the observed
      # info can be indefinite away from the mode (latent-N coupling).
      negH <- matrix(0, d, d); fisher <- matrix(0, d, d)
      for (i in seq_len(m$n_sites)) {
        obs <- m$obs_by_site[[i]]; Ji <- length(obs)
        Zi  <- matrix(0, 1L + Ji, d)
        Zi[1L, seq_len(p_lam)] <- m$X_lambda[i, ]
        if (Ji > 0L) Zi[-1L, p_lam + seq_len(p_p)] <- m$X_p[obs, , drop = FALSE]
        negH   <- negH + crossprod(Zi, m$obs_info_block(i, ev) %*% Zi)
        Fdiag  <- c(ev$info_eta_lambda[i], if (Ji > 0L) ev$info_eta_p[obs] else numeric(0))
        fisher <- fisher + crossprod(Zi, Fdiag * Zi)
      }
      list(logL = ev$log_lik, grad = grad, negH = negH, fisher = fisher)
    },
    node_ll = function(s, COEF) {
      m <- marg[[s]]
      vapply(seq_len(nrow(COEF)), function(k) {
        e <- eta_of(m, COEF[k, ])
        m$eval(e$lambda, e$p)$log_lik
      }, numeric(1))
    })
}

# Compact (ragged) counterpart of .occu_cover_site_ll: the per-site occu-cover
# marginal from length-V per-visit predictors grouped by `model$site_of_visit`,
# returned as a length-n_sites vector. Algebraically identical to the dense
# version (same MacKenzie mixture, same per-visit Beta/lognormal cover term);
# rowsum() accumulates each site's visits instead of a matrix rowSums. Scoped to
# cover_aggregate = "none" (the only mode the ragged path builds).
.occu_cover_site_ll_ragged <- function(model, psi, p_vec, ep_vec, log_disp) {
  site    <- model$site_of_visit
  y       <- model$y_det_visit
  n_sites <- model$n_sites
  g       <- factor(site, levels = seq_len(n_sites))

  log_p   <- log(p_vec)
  log_1mp <- log(1 - p_vec)
  log_h_det <- ifelse(y == 1L, log_p, log_1mp)

  positive <- model$positive %||% "lognormal"
  pos_ld   <- numeric(length(y))
  # Cover term at detected visits with an observed cover; a missing (NA) cover
  # at a detected visit drops out (missing-at-random cover), the detection term
  # above still counts.
  cobs    <- (y == 1L) & is.finite(model$y_pos_visit)
  if (any(cobs)) {
    pos_ld[cobs] <- .occu_cover_pos_logdens(model$y_pos_visit[cobs], ep_vec[cobs],
                                            exp(log_disp), positive)
  }

  sum_log_h  <- as.numeric(rowsum(log_h_det, g))
  sum_log1mp <- as.numeric(rowsum(log_1mp,   g))
  cover_term <- as.numeric(rowsum(pos_ld,    g))
  n_det      <- as.numeric(rowsum(as.numeric(y), g))

  log_psi   <- log(pmax(psi, 1e-300))
  log_1mpsi <- log(pmax(1 - psi, 1e-300))
  det_ll    <- log_psi + sum_log_h + cover_term
  nodet_ll  <- .tobs_logsumexp2(log_psi + sum_log1mp, log_1mpsi)
  ifelse(n_det > 0, det_ll, nodet_ll)
}

# Inverse of `.unscale_beta_vec`: map natural-scale betas back into the
# centered+scaled parameterization. Used for round-trip tests and for
# evaluating arm log-likelihoods (`.loglik_cover_*`) at user-facing
# natural-scale coefficients when the encoded design X is scaled.
.scale_beta_vec <- function(beta_nat, scale) {
  if (length(scale$cols) == 0L) return(beta_nat)
  beta_sc <- beta_nat
  i0 <- scale$intercept_col
  for (k in seq_along(scale$cols)) {
    j <- scale$cols[k]
    beta_sc[j] <- scale$sds[k] * beta_nat[j]
  }
  if (!is.na(i0)) {
    shift <- 0
    for (k in seq_along(scale$cols)) {
      j <- scale$cols[k]
      shift <- shift + scale$means[k] * beta_nat[j]
    }
    beta_sc[i0] <- beta_nat[i0] + shift
  }
  beta_sc
}

# Joint log-posterior + gradient of the N-mixture coefficient vector. `marg` is a
# `nmix_marginal` object (nmix_site_marginal()); `lay` the layout above. Weak
# Gaussian priors N(0, sigma.beta^2) on every coefficient and N(0, sigma.logr^2)
# on log_r keep the ridge off the boundary without materially shifting the
# data-dominated optimum, matching the other Laplace paths. Returns
# list(lp, grad) over the packed coordinates.
.tobs_abun_nuts_logpost <- function(theta, marg, lay,
                                    sigma.beta = 10, sigma.logr = 1.5) {
  beta_lambda <- theta[lay$lambda]
  beta_p      <- theta[lay$p]
  r <- if (lay$is_nb) exp(min(theta[lay$log_r], .TOBS_LOG_R_MAX)) else Inf
  ev <- marg$eval_beta(beta_lambda, beta_p, r = r)

  lp   <- ev$log_lik
  grad <- c(as.numeric(crossprod(marg$X_lambda, ev$grad_eta_lambda)),
            as.numeric(crossprod(marg$X_p,      ev$grad_eta_p)))
  if (lay$is_nb) grad <- c(grad, sum(ev$grad_theta))

  # Weak Gaussian coefficient priors (data-dominated ridge).
  ib2 <- 1 / sigma.beta^2
  lp  <- lp - 0.5 * ib2 * (sum(beta_lambda^2) + sum(beta_p^2))
  grad[lay$lambda] <- grad[lay$lambda] - ib2 * beta_lambda
  grad[lay$p]      <- grad[lay$p]      - ib2 * beta_p
  if (lay$is_nb) {
    lr  <- theta[lay$log_r]; ilr2 <- 1 / sigma.logr^2
    lp  <- lp - 0.5 * ilr2 * lr^2
    grad[lay$log_r] <- grad[lay$log_r] - ilr2 * lr
  }
  list(lp = lp, grad = grad)
}

# Joint log-posterior + gradient of the non-spatial cover-hurdle coefficient
# vector theta = c(beta_presence, beta_pos, log_disp). Weak Gaussian priors
# N(0, sigma.beta^2) on every coefficient and a broad N(0, sigma.logdisp^2) on
# log_disp keep the optimum proper without materially shifting the data-dominated
# mode. `enc` is a natural-scale cover encoding (encode_cover_hurdle with
# autoscale = FALSE): enc$occ_data carries the presence Bernoulli arm
# (`y` = 1{cover > 0}, `X`), enc$pos_data the positive arm (`y` = the arm
# response the policy consumes, `X`). The positive-arm response is
# exp(enc$pos_data$y) for lognormal (the encoder stores log(cover); the policy
# takes the log internally) and enc$pos_data$y for beta. Returns list(lp, grad)
# over the packed coordinates. This is the oracle the C++ FullGradFn mirrors.
.tobs_cover_nuts_logpost <- function(theta, enc, sigma.beta = 5,
                                     sigma.logdisp = 5) {
  X_pres <- enc$occ_data$X
  X_pos  <- enc$pos_data$X
  p_pres <- ncol(X_pres)
  p_pos  <- ncol(X_pos)
  total  <- p_pres + p_pos + 1L
  pos_code <- .occu_cover_pos_code(enc$positive)

  b_pres   <- theta[seq_len(p_pres)]
  b_pos    <- theta[p_pres + seq_len(p_pos)]
  log_disp <- theta[total]
  disp     <- exp(log_disp)

  present <- enc$occ_data$y
  y_pos   <- .tobs_cover_pos_response(enc)
  sgm     <- function(e) 1 / (1 + exp(-e))

  # Presence arm: plain Bernoulli on 1{cover > 0}.
  eta_pres <- as.numeric(X_pres %*% b_pres)
  pr       <- sgm(eta_pres)
  lp       <- sum(ifelse(present == 1L, log(pr), log(1 - pr)))
  g_eta_pres <- ifelse(present == 1L, 1 - pr, -pr)
  grad_pres  <- as.numeric(crossprod(X_pres, g_eta_pres))

  # Positive arm: beta / lognormal density at the present rows.
  eta_pos    <- as.numeric(X_pos %*% b_pos)
  g_eta_pos  <- numeric(length(eta_pos))
  g_logdisp  <- 0
  if (pos_code == 3L) {              # beta
    mu   <- sgm(eta_pos)
    a    <- mu * disp; b <- (1 - mu) * disp
    ly   <- log(y_pos); l1my <- log(1 - y_pos)
    lp   <- lp + sum(lgamma(disp) - lgamma(a) - lgamma(b) +
                     (a - 1) * ly + (b - 1) * l1my)
    g_eta_pos <- disp * mu * (1 - mu) * (-digamma(a) + digamma(b) + ly - l1my)
    g_logdisp <- sum(disp * (digamma(disp) - mu * digamma(a) -
                             (1 - mu) * digamma(b) + mu * ly + (1 - mu) * l1my))
  } else if (pos_code == 4L) {       # identity-Gaussian (#112): raw response
    sig <- disp; r <- (y_pos - eta_pos) / sig
    lp  <- lp + sum(-log(sig) - 0.5 * log(2 * pi) - 0.5 * r * r)
    g_eta_pos <- r / sig
    g_logdisp <- sum(r * r - 1)
  } else {                           # lognormal: Gaussian on log(cover)
    sig <- disp; r <- (log(y_pos) - eta_pos) / sig
    lp  <- lp + sum(-log(y_pos) - log(sig) - 0.5 * log(2 * pi) - 0.5 * r * r)
    g_eta_pos <- r / sig
    g_logdisp <- sum(r * r - 1)
  }
  grad_pos <- as.numeric(crossprod(X_pos, g_eta_pos))

  grad <- c(grad_pres, grad_pos, g_logdisp)

  # Weak Gaussian priors.
  ib2 <- 1 / sigma.beta^2
  nb  <- total - 1L
  bv  <- theta[seq_len(nb)]
  lp  <- lp - 0.5 * ib2 * sum(bv^2)
  grad[seq_len(nb)] <- grad[seq_len(nb)] - ib2 * bv
  ild2 <- 1 / sigma.logdisp^2
  lp   <- lp - 0.5 * ild2 * log_disp^2
  grad[total] <- grad[total] - ild2 * log_disp
  list(lp = lp, grad = grad)
}

# Joint log-posterior + gradient of the distance coefficient vector (R oracle).
.tobs_distance_nuts_logpost <- function(theta, marg, lay,
                                        sigma.beta = 10, sigma.shape = 1.5,
                                        sigma.logr = 1.5) {
  beta_lambda <- theta[lay$lambda]
  beta_sigma  <- theta[lay$sigma]
  eta_b <- if (lay$hazard) theta[lay$log_shape] else 0
  r     <- if (lay$is_nb)  exp(min(theta[lay$log_r], .TOBS_LOG_R_MAX)) else Inf
  ev <- marg$eval_beta(beta_lambda, beta_sigma, eta_b = eta_b, r = r)

  lp   <- ev$log_lik
  grad <- c(as.numeric(crossprod(marg$X_lambda, ev$grad_eta_lambda)),
            as.numeric(crossprod(marg$X_sigma,  ev$grad_eta_sigma)))
  if (lay$hazard) grad <- c(grad, ev$grad_eta_b)
  if (lay$is_nb)  grad <- c(grad, sum(ev$grad_theta))

  ib2 <- 1 / sigma.beta^2
  lp  <- lp - 0.5 * ib2 * (sum(beta_lambda^2) + sum(beta_sigma^2))
  grad[lay$lambda] <- grad[lay$lambda] - ib2 * beta_lambda
  grad[lay$sigma]  <- grad[lay$sigma]  - ib2 * beta_sigma
  if (lay$hazard) {
    is2 <- 1 / sigma.shape^2
    lp  <- lp - 0.5 * is2 * eta_b^2
    grad[lay$log_shape] <- grad[lay$log_shape] - is2 * eta_b
  }
  if (lay$is_nb) {
    lr  <- theta[lay$log_r]; ilr2 <- 1 / sigma.logr^2
    lp  <- lp - 0.5 * ilr2 * lr^2
    grad[lay$log_r] <- grad[lay$log_r] - ilr2 * lr
  }
  list(lp = lp, grad = grad)
}

.tobs_dyn_abun_nuts_logpost <- function(theta, marg, lay, sigma.beta = 10,
                                        sigma.logr = 1.5) {
  use_nb <- isTRUE(lay$use_nb)
  eta_logr <- if (use_nb) min(theta[lay$logr], .TOBS_LOG_R_MAX) else 0
  ev <- marg$eval_beta(theta[lay$lambda], theta[lay$p], theta[lay$omega],
                       theta[lay$gamma], eta_logr)
  arms <- list(
    list(idx = lay$lambda, X = marg$X_lambda, grad = "grad_eta_lambda"),
    list(idx = lay$p,      X = marg$X_p,      grad = "grad_eta_p"),
    list(idx = lay$omega,  X = marg$X_omega,  grad = "grad_eta_omega"),
    list(idx = lay$gamma,  X = marg$X_gamma,  grad = "grad_eta_gamma"))
  # The dispersion log r has no design: its gradient is the scalar
  # grad_eta_logr (already summed over sites), folded in as a 1x1 arm.
  if (use_nb) {
    ev$grad_eta_logr <- as.numeric(ev$grad_eta_logr)
    arms <- c(arms, list(list(idx = lay$logr, X = matrix(1, 1, 1),
                              grad = "grad_eta_logr")))
  }
  # Per-coordinate prior SD: sigma.beta on the coefficients, sigma.logr on log_r.
  prior_sd <- rep(sigma.beta, lay$total)
  if (use_nb) prior_sd[lay$logr] <- sigma.logr
  .tobs_nuts_logpost_k(theta, ev, arms, lay$total, prior_sd)
}

.tobs_fp_occu_nuts_logpost <- function(theta, marg, lay, sigma.beta = 10) {
  ev <- marg$eval_beta(theta[lay$psi], theta[lay$p11], theta[lay$p10], theta[lay$b])
  arms <- list(
    list(idx = lay$psi, X = marg$X_psi, grad = "grad_eta_psi"),
    list(idx = lay$p11, X = marg$X_p11, grad = "grad_eta_p11"),
    list(idx = lay$p10, X = marg$X_p10, grad = "grad_eta_p10"),
    list(idx = lay$b,   X = marg$X_b,   grad = "grad_eta_b"))
  .tobs_nuts_logpost_k(theta, ev, arms, lay$total, sigma.beta)
}

# Full-vector joint log-posterior and its gradient for the community N-mixture.
# `margs` is the per-species list from .tobs_ms_abun_nuts_marginals(); `lay` the
# layout. Returns list(lp, grad) over the packed coordinates. Mirrors the C++
# ms_abun_nuts_eval (src/ms_abun_nuts.cpp).
#
# NON-CENTERED parameterisation: the per-species block holds standard-normal
# z_s, and the deviation is reconstructed per arm as b_{s,arm} = C_arm z_{s,arm}
# (C_arm the log-Cholesky factor of Sigma_arm). The community covariance then
# leaves the b-prior entirely -- z ~ N(0, I), no b' Sigma^{-1} b quadratic and no
# -0.5 S log|Sigma| normaliser -- and enters ONLY the data term through b = C z.
# This breaks the centered b<->Sigma funnel that otherwise saturates the NUTS
# treedepth (mean ~7, eps ~0.025 measured for the centered map); the chol
# gradient now flows from the data term via b = C z plus the hyperprior.
.tobs_ms_abun_nuts_logpost <- function(theta, margs, lay, priors,
                                       sigma.beta = 10, sigma.logr = 1.5,
                                       grad = TRUE) {
  P <- lay$P; S <- lay$n_species; is_nb <- lay$is_nb
  mu <- theta[lay$mu]
  g    <- numeric(lay$total)
  g_mu <- numeric(P)
  lp   <- 0

  # Cholesky factors per arm (Sigma_arm = C_arm C_arm'); logr arm is the 1x1
  # scalar SD = exp(chol_logr).
  C_lam <- .ms_ocs_chol_unpack(theta[lay$chol_lam], lay$p_lam)
  C_p   <- .ms_ocs_chol_unpack(theta[lay$chol_p],   lay$p_p)
  C_lr  <- if (is_nb) exp(theta[lay$chol_logr]) else NULL

  # chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
  A_lam <- matrix(0, lay$p_lam, lay$p_lam)
  A_p   <- matrix(0, lay$p_p,   lay$p_p)
  A_lr  <- 0

  # ---- data log-lik + inner gradient (per species), non-centered b = C z ----
  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]
    zl <- z_s[lay$lambda]; zp <- z_s[lay$p]
    bl   <- mu[lay$lambda] + as.numeric(C_lam %*% zl)
    bp   <- mu[lay$p]      + as.numeric(C_p   %*% zp)
    if (is_nb) {
      zr <- z_s[lay$logr]
      r  <- exp(min(mu[lay$logr] + C_lr * zr, .TOBS_LOG_R_MAX))
    } else {
      r <- Inf
    }
    ev   <- margs[[s]]$eval_beta(bl, bp, r = r)
    lp   <- lp + ev$log_lik
    if (grad) {
      gl <- as.numeric(crossprod(margs[[s]]$X_lambda, ev$grad_eta_lambda)) # grad_b lambda
      gp <- as.numeric(crossprod(margs[[s]]$X_p,      ev$grad_eta_p))      # grad_b p
      g_mu[lay$lambda]    <- g_mu[lay$lambda] + gl
      g_mu[lay$p]         <- g_mu[lay$p]      + gp
      # z gradient (data part) = C' grad_b; the -z prior is added below.
      g[bidx[lay$lambda]] <- g[bidx[lay$lambda]] + as.numeric(crossprod(C_lam, gl))
      g[bidx[lay$p]]      <- g[bidx[lay$p]]      + as.numeric(crossprod(C_p,   gp))
      A_lam <- A_lam + outer(gl, zl)
      A_p   <- A_p   + outer(gp, zp)
      if (is_nb) {
        gt <- sum(ev$grad_theta)                                          # grad_b logr
        g_mu[lay$logr]    <- g_mu[lay$logr] + gt
        g[bidx[lay$logr]] <- g[bidx[lay$logr]] + C_lr * gt
        A_lr <- A_lr + gt * zr
      }
    }
  }

  # ---- z prior: standard normal over the entire per-species block ----
  z_idx <- lay$b_off + seq_len(S * P)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior ----
  arms <- list(list(chol = lay$chol_lam, A = A_lam, C = C_lam, Pa = lay$p_lam),
               list(chol = lay$chol_p,   A = A_p,   C = C_p,   Pa = lay$p_p))
  if (is_nb) arms <- c(arms, list(list(chol = lay$chol_logr,
                                       A = matrix(A_lr, 1, 1),
                                       C = matrix(C_lr, 1, 1), Pa = 1L)))
  for (arm in arms) {
    Pa <- arm$Pa; if (Pa == 0L) next
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_abun_nuts_chol_data_grad(arm$A, arm$C, Pa) +
        pr$grad
  }

  # ---- community-mean priors ----
  ib2 <- 1 / sigma.beta^2
  coef_idx <- c(lay$lambda, lay$p)
  lp <- lp - 0.5 * ib2 * sum(mu[coef_idx]^2)
  g_mu[coef_idx] <- g_mu[coef_idx] - ib2 * mu[coef_idx]
  if (is_nb) {
    il2 <- 1 / sigma.logr^2
    lp  <- lp - 0.5 * il2 * mu[lay$logr]^2
    g_mu[lay$logr] <- g_mu[lay$logr] - il2 * mu[lay$logr]
  }
  g[lay$mu] <- g[lay$mu] + g_mu

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Full-vector joint log-posterior + gradient (the NUTS target / oracle). `Y` is
# the n_sites x n_species response matrix; `X` the site-level design. An NA entry
# Y[i, s] is a missing site x species observation: that (species, site) is dropped
# species keeps only its observed sites. Mirrors the C++ ms_count_nuts_eval
# (src/ms_count_nuts.cpp), including the NA skip.
.tobs_ms_count_nuts_logpost <- function(theta, X, Y, lay, priors,
                                        sigma.beta = 10, sigma.logr = 1.5,
                                        grad = TRUE) {
  pb <- lay$p_beta; S <- lay$n_species; P <- lay$P
  is_nb <- lay$is_nb; is_gauss <- lay$is_gauss
  is_bern <- identical(lay$family, "bernoulli")
  mu     <- theta[lay$mu]
  C_beta <- .ms_ocs_chol_unpack(theta[lay$chol_beta], pb)
  C_lr   <- if (is_nb) exp(theta[lay$chol_logr]) else NULL
  logphi <- if (is_gauss) theta[lay$logphi] else NULL

  g <- numeric(lay$total); g_mu <- numeric(P)
  A_beta <- matrix(0, pb, pb); A_lr <- 0; lp <- 0

  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]; zb <- z_s[lay$beta]
    b_beta <- mu[lay$beta] + as.numeric(C_beta %*% zb)
    # Drop missing (NA) site x species observations for this species.
    ys_all <- Y[, s]
    obs    <- !is.na(ys_all)
    Xs     <- X[obs, , drop = FALSE]
    ys     <- ys_all[obs]
    eta    <- as.numeric(Xs %*% b_beta)
    if (is_nb) {
      zr  <- z_s[lay$logr]
      r   <- exp(min(mu[lay$logr] + C_lr * zr, .TOBS_LOG_R_MAX))
      muv <- pmax(exp(pmin(eta, 700)), 1e-10)
      lp  <- lp + sum(stats::dnbinom(ys, size = r, mu = muv, log = TRUE))
      if (grad) {
        gl     <- as.numeric(crossprod(Xs, r * (ys - muv) / (r + muv)))
        dLL_dr <- digamma(ys + r) - digamma(r) + log(r / (r + muv)) + 1 -
                  (ys + r) / (r + muv)
        g_logr <- r * sum(dLL_dr)
        g_mu[lay$beta] <- g_mu[lay$beta] + gl
        g_mu[lay$logr] <- g_mu[lay$logr] + g_logr
        g[bidx[lay$beta]] <- g[bidx[lay$beta]] + as.numeric(crossprod(C_beta, gl))
        g[bidx[lay$logr]] <- g[bidx[lay$logr]] + C_lr * g_logr
        A_beta <- A_beta + outer(gl, zb)
        A_lr   <- A_lr + g_logr * zr
      }
    } else if (is_gauss) {
      phi_s <- exp(min(logphi[s], 30))
      lp    <- lp + sum(stats::dnorm(ys, mean = eta, sd = sqrt(phi_s), log = TRUE))
      if (grad) {
        gl <- as.numeric(crossprod(Xs, (ys - eta) / phi_s))
        g_mu[lay$beta] <- g_mu[lay$beta] + gl
        g[bidx[lay$beta]] <- g[bidx[lay$beta]] + as.numeric(crossprod(C_beta, gl))
        g[lay$logphi[s]]  <- g[lay$logphi[s]] + sum(-0.5 + (ys - eta)^2 / (2 * phi_s))
        A_beta <- A_beta + outer(gl, zb)
      }
    } else if (is_bern) {
      # jsdm(): log p = y eta - log(1 + exp(eta)); d log p / d eta = y - plogis(eta).
      # Written via plogis(log.p) to match the C++ stable log-sum-exp branch.
      lp <- lp + sum(ifelse(ys > 0, stats::plogis(eta,  log.p = TRUE),
                                    stats::plogis(-eta, log.p = TRUE)))
      if (grad) {
        gl <- as.numeric(crossprod(Xs, ys - stats::plogis(eta)))
        g_mu[lay$beta] <- g_mu[lay$beta] + gl
        g[bidx[lay$beta]] <- g[bidx[lay$beta]] + as.numeric(crossprod(C_beta, gl))
        A_beta <- A_beta + outer(gl, zb)
      }
    } else {
      lam <- exp(pmin(eta, 700))
      lp  <- lp + sum(ys * eta - lam - lgamma(ys + 1))
      if (grad) {
        gl <- as.numeric(crossprod(Xs, ys - lam))
        g_mu[lay$beta] <- g_mu[lay$beta] + gl
        g[bidx[lay$beta]] <- g[bidx[lay$beta]] + as.numeric(crossprod(C_beta, gl))
        A_beta <- A_beta + outer(gl, zb)
      }
    }
  }
  # z prior N(0, I)
  z_idx <- lay$b_off + seq_len(S * P); z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all
  # chol_beta coords: data gradient (b = C z) + log-Cholesky hyperprior
  pr <- .ms_ocs_chol_logprior(theta[lay$chol_beta], pb, priors)
  lp <- lp + pr$lp
  if (grad) g[lay$chol_beta] <- .ms_abun_nuts_chol_data_grad(A_beta, C_beta, pb) +
      pr$grad
  if (is_nb) {                                # scalar log-dispersion covariance
    pr_r <- .ms_ocs_chol_logprior(theta[lay$chol_logr], 1L, priors)
    lp   <- lp + pr_r$lp
    if (grad) g[lay$chol_logr] <-
      .ms_abun_nuts_chol_data_grad(matrix(A_lr, 1, 1), matrix(C_lr, 1, 1), 1L) +
      pr_r$grad
  }
  if (is_gauss) {                             # free log_phi prior N(mean, sd^2)
    lpm <- priors$logphi_mean; lps <- priors$logphi_sd
    lp  <- lp - 0.5 * sum(((logphi - lpm) / lps)^2)
    if (grad) g[lay$logphi] <- g[lay$logphi] - (logphi - lpm) / lps^2
  }
  # community-mean priors
  ib2 <- 1 / sigma.beta^2
  lp  <- lp - 0.5 * ib2 * sum(mu[lay$beta]^2)
  g_mu[lay$beta] <- g_mu[lay$beta] - ib2 * mu[lay$beta]
  if (is_nb) {
    il2 <- 1 / sigma.logr^2
    lp  <- lp - 0.5 * il2 * mu[lay$logr]^2
    g_mu[lay$logr] <- g_mu[lay$logr] - il2 * mu[lay$logr]
  }
  g[lay$mu] <- g[lay$mu] + g_mu
  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Full-vector joint log-posterior and its gradient for the community dynamic
# occupancy model. `X_psi1` / `X_p` / `X_gamma` / `X_eps` are the site-level
# designs; `em_stats` a per-species list of precomputed (nvalid, ndet)
# [n_sites x n_seasons] detection sufficient statistics (from
# .ms_dyn_occu_emit_stats). `lay` the layout. Returns list(lp, grad) over the
# packed coordinates. Mirrors the C++ ms_dyn_occu_nuts_eval.
#
# NON-CENTERED: the per-species block holds standard-normal z_s, the deviation is
# b_{s,arm} = C_arm z_{s,arm}, so the community covariance enters ONLY the data
# term through b = C z. The gamma / eps globals are shared across species and so
# accumulate their gradient over the full species loop.
.tobs_ms_dyn_occu_nuts_logpost <- function(theta, X_psi1, X_p, X_gamma, X_eps,
                                           em_stats, n_sites, n_seasons, lay,
                                           priors, sigma.beta = 5, grad = TRUE) {
  P <- lay$P; S <- lay$n_species
  mu     <- theta[lay$mu]
  global <- theta[lay$global]
  beta_gam <- global[lay$gam]
  beta_eps <- global[lay$eps]
  eta_gam0 <- as.numeric(X_gamma %*% beta_gam)
  eta_eps0 <- as.numeric(X_eps   %*% beta_eps)
  gamma <- stats::plogis(eta_gam0)
  eps   <- stats::plogis(eta_eps0)

  g      <- numeric(lay$total)
  g_mu   <- numeric(P)
  g_gam  <- numeric(n_sites)         # accumulated d ll / d eta_gam over species
  g_eps  <- numeric(n_sites)
  lp     <- 0

  C_psi1 <- .ms_ocs_chol_unpack(theta[lay$chol_psi1], lay$p_psi1)
  C_p    <- .ms_ocs_chol_unpack(theta[lay$chol_p],    lay$p_p)

  # chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
  A_psi1 <- matrix(0, lay$p_psi1, lay$p_psi1)
  A_p    <- matrix(0, lay$p_p,    lay$p_p)

  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]
    zpsi1 <- z_s[lay$psi1]; zp <- z_s[lay$p]
    bpsi1 <- mu[lay$psi1] + as.numeric(C_psi1 %*% zpsi1)
    bp    <- mu[lay$p]    + as.numeric(C_p    %*% zp)
    eta_psi1 <- as.numeric(X_psi1 %*% bpsi1)
    eta_p    <- as.numeric(X_p    %*% bp)
    psi1 <- stats::plogis(eta_psi1)
    p    <- stats::plogis(eta_p)

    em <- .ms_dyn_occu_emissions(p, em_stats[[s]]$nvalid, em_stats[[s]]$ndet)
    if (grad) {
      fb <- .ms_dyn_occu_fb_vec(psi1, p, gamma, eps, em,
                                em_stats[[s]]$nvalid, em_stats[[s]]$ndet,
                                n_sites, n_seasons)
      lp <- lp + fb$ll
      # per-arm eta scores (Fisher identity).
      s_psi1 <- fb$w1 - psi1
      s_p    <- fb$p_score
      # design-sandwiched b gradients (psi1, p arms).
      gpsi1 <- as.numeric(crossprod(X_psi1, s_psi1))
      gp    <- as.numeric(crossprod(X_p,    s_p))
      g_mu[lay$psi1] <- g_mu[lay$psi1] + gpsi1
      g_mu[lay$p]    <- g_mu[lay$p]    + gp
      # z gradient (data part) = C' grad_b.
      g[bidx[lay$psi1]] <- g[bidx[lay$psi1]] + as.numeric(crossprod(C_psi1, gpsi1))
      g[bidx[lay$p]]    <- g[bidx[lay$p]]    + as.numeric(crossprod(C_p,    gp))
      A_psi1 <- A_psi1 + outer(gpsi1, zpsi1)
      A_p    <- A_p    + outer(gp,    zp)
      # shared transition eta scores, accumulated across species.
      g_gam <- g_gam + (fb$col_y - gamma * fb$col_n)
      g_eps <- g_eps + (fb$ext_y - eps   * fb$ext_n)
    } else {
      lp <- lp + .ms_dyn_occu_fwd_ll_vec(psi1, gamma, eps, em,
                                         n_sites, n_seasons)
    }
  }

  # ---- z prior: standard normal over the entire per-species block ----
  z_idx <- lay$b_off + seq_len(S * P)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior ----
  arms <- list(
    list(chol = lay$chol_psi1, A = A_psi1, C = C_psi1, Pa = lay$p_psi1),
    list(chol = lay$chol_p,    A = A_p,    C = C_p,    Pa = lay$p_p))
  for (arm in arms) {
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], arm$Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_abun_nuts_chol_data_grad(arm$A, arm$C, arm$Pa) +
        pr$grad
  }

  # ---- community-mean priors ----
  ib2 <- 1 / sigma.beta^2
  lp <- lp - 0.5 * ib2 * sum(mu^2)
  g_mu <- g_mu - ib2 * mu
  if (grad) g[lay$mu] <- g[lay$mu] + g_mu

  # ---- shared-transition (global) priors + design-sandwiched gradient ----
  lp <- lp - 0.5 * ib2 * sum(global^2)
  if (grad) {
    g_glob <- numeric(lay$G)
    g_glob[lay$gam] <- as.numeric(crossprod(X_gamma, g_gam)) - ib2 * beta_gam
    g_glob[lay$eps] <- as.numeric(crossprod(X_eps,   g_eps)) - ib2 * beta_eps
    g[lay$global] <- g_glob
  }

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Full-vector joint log-posterior and its gradient for the community integrated
# occupancy model. `X_psi` [n_sites x P_psi]; `X_p` a list of D site-level
# detection designs; `summaries` the per-species .ms_int_occu_sp_summary list;
# `lay` the layout. Returns list(lp, grad) over the packed coordinates.
#
# NON-CENTERED: the per-species block holds standard-normal z_s, the deviation is
# b_{s,arm} = C_arm z_{s,arm}, so each community covariance leaves the b-prior
# (z ~ N(0, I)) and enters ONLY the data term through b = C z.
.tobs_ms_int_occu_nuts_logpost <- function(theta, X_psi, X_p, summaries, lay,
                                           priors, sigma.beta = 5, grad = TRUE) {
  P <- lay$P; S <- lay$n_species; D <- lay$D
  mu <- theta[lay$mu]
  g    <- numeric(lay$total)
  g_mu <- numeric(P)
  lp   <- 0

  # Unpack the per-arm Cholesky factors (psi + D detection).
  C_psi <- .ms_ocs_chol_unpack(theta[lay$chol_psi], lay$P_psi)
  C_p   <- lapply(seq_len(D),
                  function(d) .ms_ocs_chol_unpack(theta[lay$chol_p[[d]]], lay$P_p[d]))

  # chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
  A_psi <- matrix(0, lay$P_psi, lay$P_psi)
  A_p   <- lapply(seq_len(D), function(d) matrix(0, lay$P_p[d], lay$P_p[d]))

  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]
    zpsi <- z_s[lay$psi]
    bpsi <- mu[lay$psi] + as.numeric(C_psi %*% zpsi)
    eta_psi <- as.numeric(X_psi %*% bpsi)
    zp   <- vector("list", D); eta_p <- vector("list", D)
    for (d in seq_len(D)) {
      zp[[d]]    <- z_s[lay$p[[d]]]
      bpd        <- mu[lay$p[[d]]] + as.numeric(C_p[[d]] %*% zp[[d]])
      eta_p[[d]] <- as.numeric(X_p[[d]] %*% bpd)
    }
    lp <- lp + .ms_int_occu_sp_ll(eta_psi, eta_p, summaries[[s]])
    if (grad) {
      gvec <- .ms_int_occu_sp_grad(eta_psi, eta_p, summaries[[s]], X_psi, X_p)
      gpsi <- gvec[lay$psi]                       # grad_b psi (coefficient space)
      g_mu[lay$psi] <- g_mu[lay$psi] + gpsi
      g[bidx[lay$psi]] <- g[bidx[lay$psi]] + as.numeric(crossprod(C_psi, gpsi))
      A_psi <- A_psi + outer(gpsi, zpsi)
      for (d in seq_len(D)) {
        gpd <- gvec[lay$p[[d]]]
        g_mu[lay$p[[d]]] <- g_mu[lay$p[[d]]] + gpd
        g[bidx[lay$p[[d]]]] <- g[bidx[lay$p[[d]]]] +
          as.numeric(crossprod(C_p[[d]], gpd))
        A_p[[d]] <- A_p[[d]] + outer(gpd, zp[[d]])
      }
    }
  }

  # ---- z prior: standard normal over the entire per-species block ----
  z_idx <- lay$b_off + seq_len(S * P)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior, per arm ----
  arms <- c(list(list(chol = lay$chol_psi, A = A_psi, C = C_psi, Pa = lay$P_psi)),
            lapply(seq_len(D), function(d)
              list(chol = lay$chol_p[[d]], A = A_p[[d]], C = C_p[[d]], Pa = lay$P_p[d])))
  for (arm in arms) {
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], arm$Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_abun_nuts_chol_data_grad(arm$A, arm$C, arm$Pa) +
        pr$grad
  }

  # ---- community-mean priors ----
  ib2 <- 1 / sigma.beta^2
  lp <- lp - 0.5 * ib2 * sum(mu^2)
  g_mu <- g_mu - ib2 * mu
  g[lay$mu] <- g[lay$mu] + g_mu

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Full-vector joint log-posterior and its gradient. `views` is the per-species
# .ms_occu_cover_species_view list (the single-species occu_cover model views the
# laplace path builds); `lay` the layout. Returns list(lp, grad) over the packed
# coordinates.
#
# NON-CENTERED: the per-species block holds standard-normal z_s, the deviation is
# b_{s,arm} = C_arm z_{s,arm}, so each community covariance leaves the b-prior and
# enters ONLY the data term. The shared log_disp carries a weakly-informative
# Normal prior on the log scale and accumulates its gradient over all species.
.tobs_ms_occu_cover_nuts_logpost <- function(theta, views, lay, priors,
                                             sigma.beta = 5, grad = TRUE) {
  P <- lay$P; S <- lay$n_species
  mu <- theta[lay$mu]
  log_disp <- theta[lay$log_disp]
  g    <- numeric(lay$total)
  g_mu <- numeric(P)
  g_ld <- 0
  lp   <- 0

  C_occ <- .ms_ocs_chol_unpack(theta[lay$chol_occ], lay$P_occ)
  C_p   <- .ms_ocs_chol_unpack(theta[lay$chol_p],   lay$P_p)
  C_pos <- .ms_ocs_chol_unpack(theta[lay$chol_pos], lay$P_pos)

  A_occ <- matrix(0, lay$P_occ, lay$P_occ)
  A_p   <- matrix(0, lay$P_p,   lay$P_p)
  A_pos <- matrix(0, lay$P_pos, lay$P_pos)

  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]
    zocc <- z_s[lay$occ]; zp <- z_s[lay$p]; zpos <- z_s[lay$pos]
    bocc <- mu[lay$occ] + as.numeric(C_occ %*% zocc)
    bp   <- mu[lay$p]   + as.numeric(C_p   %*% zp)
    bpos <- mu[lay$pos] + as.numeric(C_pos %*% zpos)

    lp <- lp + .occu_cover_sp_ll(views[[s]], bocc, bp, bpos, log_disp)
    if (grad) {
      gvec <- .occu_cover_sp_grad(views[[s]], bocc, bp, bpos, log_disp)
      gocc <- gvec[lay$occ]; gp <- gvec[lay$p]; gpos <- gvec[lay$pos]
      g_ld <- g_ld + gvec[P + 1L]
      g_mu[lay$occ] <- g_mu[lay$occ] + gocc
      g_mu[lay$p]   <- g_mu[lay$p]   + gp
      g_mu[lay$pos] <- g_mu[lay$pos] + gpos
      g[bidx[lay$occ]] <- g[bidx[lay$occ]] + as.numeric(crossprod(C_occ, gocc))
      g[bidx[lay$p]]   <- g[bidx[lay$p]]   + as.numeric(crossprod(C_p,   gp))
      g[bidx[lay$pos]] <- g[bidx[lay$pos]] + as.numeric(crossprod(C_pos, gpos))
      A_occ <- A_occ + outer(gocc, zocc)
      A_p   <- A_p   + outer(gp,   zp)
      A_pos <- A_pos + outer(gpos, zpos)
    }
  }

  # ---- z prior: standard normal over the entire per-species block ----
  z_idx <- lay$b_off + seq_len(S * P)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior, per arm ----
  arms <- list(list(chol = lay$chol_occ, A = A_occ, C = C_occ, Pa = lay$P_occ),
               list(chol = lay$chol_p,   A = A_p,   C = C_p,   Pa = lay$P_p),
               list(chol = lay$chol_pos, A = A_pos, C = C_pos, Pa = lay$P_pos))
  for (arm in arms) {
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], arm$Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_abun_nuts_chol_data_grad(arm$A, arm$C, arm$Pa) +
        pr$grad
  }

  # ---- community-mean priors ----
  ib2 <- 1 / sigma.beta^2
  lp <- lp - 0.5 * ib2 * sum(mu^2)
  g_mu <- g_mu - ib2 * mu
  g[lay$mu] <- g[lay$mu] + g_mu

  # ---- shared log-dispersion prior (weakly-informative Normal on log scale) ----
  ld_mean <- priors$log_disp_mean %||% log(0.5)
  ld_sd   <- priors$log_disp_sd   %||% 2.0
  lp <- lp - 0.5 * ((log_disp - ld_mean) / ld_sd)^2
  g_ld <- g_ld - (log_disp - ld_mean) / ld_sd^2
  if (grad) g[lay$log_disp] <- g_ld

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Joint log-posterior + gradient for the dispersion-RE variant. `views` the
# per-species occu_cover model views; `lay` the RE layout. Non-centered on all
# four arms (occ/p/pos coefficients + the 1-D log-dispersion).
.tobs_ms_occu_cover_re_disp_logpost <- function(theta, views, lay, priors,
                                                sigma.beta = 5, grad = TRUE) {
  S <- lay$n_species
  mu_coef <- theta[lay$mu_coef]
  mu_ld   <- theta[lay$mu_ld]
  g    <- numeric(lay$total)
  g_mc <- numeric(lay$P_coef)
  g_mld <- 0
  lp   <- 0

  C_occ <- .ms_ocs_chol_unpack(theta[lay$chol_occ], lay$P_occ)
  C_p   <- .ms_ocs_chol_unpack(theta[lay$chol_p],   lay$P_p)
  C_pos <- .ms_ocs_chol_unpack(theta[lay$chol_pos], lay$P_pos)
  sigma_ld <- exp(theta[lay$chol_ld])            # 1x1 Cholesky factor

  A_occ <- matrix(0, lay$P_occ, lay$P_occ)
  A_p   <- matrix(0, lay$P_p,   lay$P_p)
  A_pos <- matrix(0, lay$P_pos, lay$P_pos)
  A_ld  <- 0

  for (s in seq_len(S)) {
    bidx <- .tobs_ms_occu_cover_re_disp_b_idx(lay, s)
    z_s  <- theta[bidx]
    zocc <- z_s[lay$occ]; zp <- z_s[lay$p]; zpos <- z_s[lay$pos]; zld <- z_s[lay$ld]
    bocc <- mu_coef[lay$occ] + as.numeric(C_occ %*% zocc)
    bp   <- mu_coef[lay$p]   + as.numeric(C_p   %*% zp)
    bpos <- mu_coef[lay$pos] + as.numeric(C_pos %*% zpos)
    log_disp_s <- mu_ld + sigma_ld * zld

    lp <- lp + .occu_cover_sp_ll(views[[s]], bocc, bp, bpos, log_disp_s)
    if (grad) {
      gvec <- .occu_cover_sp_grad(views[[s]], bocc, bp, bpos, log_disp_s)
      gocc <- gvec[lay$occ]; gp <- gvec[lay$p]; gpos <- gvec[lay$pos]
      g_ld_s <- gvec[lay$P_coef + 1L]
      g_mc[lay$occ] <- g_mc[lay$occ] + gocc
      g_mc[lay$p]   <- g_mc[lay$p]   + gp
      g_mc[lay$pos] <- g_mc[lay$pos] + gpos
      g_mld <- g_mld + g_ld_s
      g[bidx[lay$occ]] <- g[bidx[lay$occ]] + as.numeric(crossprod(C_occ, gocc))
      g[bidx[lay$p]]   <- g[bidx[lay$p]]   + as.numeric(crossprod(C_p,   gp))
      g[bidx[lay$pos]] <- g[bidx[lay$pos]] + as.numeric(crossprod(C_pos, gpos))
      g[bidx[lay$ld]]  <- g[bidx[lay$ld]]  + sigma_ld * g_ld_s
      A_occ <- A_occ + outer(gocc, zocc)
      A_p   <- A_p   + outer(gp,   zp)
      A_pos <- A_pos + outer(gpos, zpos)
      A_ld  <- A_ld  + g_ld_s * zld
    }
  }

  # ---- z prior: standard normal over the whole per-species block ----
  z_idx <- lay$b_off + seq_len(S * lay$Pz)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior, per arm ----
  arms <- list(list(chol = lay$chol_occ, A = A_occ, C = C_occ, Pa = lay$P_occ),
               list(chol = lay$chol_p,   A = A_p,   C = C_p,   Pa = lay$P_p),
               list(chol = lay$chol_pos, A = A_pos, C = C_pos, Pa = lay$P_pos),
               list(chol = lay$chol_ld,  A = matrix(A_ld, 1, 1),
                    C = matrix(sigma_ld, 1, 1), Pa = 1L))
  for (arm in arms) {
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], arm$Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_abun_nuts_chol_data_grad(arm$A, arm$C, arm$Pa) +
        pr$grad
  }

  # ---- community-mean priors: coefficients ridge + log-dispersion prior ----
  ib2 <- 1 / sigma.beta^2
  lp <- lp - 0.5 * ib2 * sum(mu_coef^2)
  g_mc <- g_mc - ib2 * mu_coef
  ld_mean <- priors$log_disp_mean %||% log(0.5)
  ld_sd   <- priors$log_disp_sd   %||% 2.0
  lp <- lp - 0.5 * ((mu_ld - ld_mean) / ld_sd)^2
  g_mld <- g_mld - (mu_ld - ld_mean) / ld_sd^2
  if (grad) {
    g[lay$mu_coef] <- g[lay$mu_coef] + g_mc
    g[lay$mu_ld]   <- g_mld
  }

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Full-vector joint log-posterior and its gradient for the community occupancy
# model. `X_psi` / `X_p` are the site-level designs; `summaries` the per-species
# .ms_int_occu_sp_summary list (single-source). `lay` the layout. Returns
# list(lp, grad) over the packed coordinates. Mirrors the C++ ms_occu_nuts_eval
# (src/ms_occu_nuts.cpp).
#
# NON-CENTERED: the per-species block holds standard-normal z_s, the deviation is
# b_{s,arm} = C_arm z_{s,arm} (C_arm the log-Cholesky factor of Sigma_arm), so the
# community covariance leaves the b-prior (z ~ N(0, I)) and enters ONLY the data
# term through b = C z.
.tobs_ms_occu_nuts_logpost <- function(theta, X_psi, X_p, summaries, lay, priors,
                                       sigma.beta = 5, grad = TRUE) {
  P <- lay$P; S <- lay$n_species
  mu <- theta[lay$mu]
  g    <- numeric(lay$total)
  g_mu <- numeric(P)
  lp   <- 0

  C_psi <- .ms_ocs_chol_unpack(theta[lay$chol_psi], lay$p_psi)
  C_p   <- .ms_ocs_chol_unpack(theta[lay$chol_p],   lay$p_p)

  # chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
  A_psi <- matrix(0, lay$p_psi, lay$p_psi)
  A_p   <- matrix(0, lay$p_p,   lay$p_p)

  for (s in seq_len(S)) {
    bidx <- .ms_ocs_b_idx(lay, s)
    z_s  <- theta[bidx]
    zpsi <- z_s[lay$psi]; zp <- z_s[lay$p]
    bpsi <- mu[lay$psi] + as.numeric(C_psi %*% zpsi)
    bp   <- mu[lay$p]   + as.numeric(C_p   %*% zp)
    eta_psi <- as.numeric(X_psi %*% bpsi)
    eta_p   <- as.numeric(X_p   %*% bp)
    lp <- lp + .ms_int_occu_sp_ll(eta_psi, list(eta_p), summaries[[s]])
    if (grad) {
      gvec <- .ms_int_occu_sp_grad(eta_psi, list(eta_p), summaries[[s]],
                                   X_psi, list(X_p))
      gpsi <- gvec[lay$psi]              # grad_b psi
      gp   <- gvec[lay$p]               # grad_b p
      g_mu[lay$psi] <- g_mu[lay$psi] + gpsi
      g_mu[lay$p]   <- g_mu[lay$p]   + gp
      # z gradient (data part) = C' grad_b.
      g[bidx[lay$psi]] <- g[bidx[lay$psi]] + as.numeric(crossprod(C_psi, gpsi))
      g[bidx[lay$p]]   <- g[bidx[lay$p]]   + as.numeric(crossprod(C_p,   gp))
      A_psi <- A_psi + outer(gpsi, zpsi)
      A_p   <- A_p   + outer(gp,   zp)
    }
  }

  # ---- z prior: standard normal over the entire per-species block ----
  z_idx <- lay$b_off + seq_len(S * P)
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all

  # ---- chol coords: data gradient (via b = C z) + hyperprior ----
  arms <- list(list(chol = lay$chol_psi, A = A_psi, C = C_psi, Pa = lay$p_psi),
               list(chol = lay$chol_p,   A = A_p,   C = C_p,   Pa = lay$p_p))
  for (arm in arms) {
    Pa <- arm$Pa
    pr <- .ms_ocs_chol_logprior(theta[arm$chol], Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[arm$chol] <- .ms_abun_nuts_chol_data_grad(arm$A, arm$C, Pa) +
        pr$grad
  }

  # ---- community-mean priors ----
  ib2 <- 1 / sigma.beta^2
  coef_idx <- c(lay$psi, lay$p)
  lp <- lp - 0.5 * ib2 * sum(mu[coef_idx]^2)
  g_mu[coef_idx] <- g_mu[coef_idx] - ib2 * mu[coef_idx]
  g[lay$mu] <- g[lay$mu] + g_mu

  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

# Joint log-posterior + gradient of the non-spatial occu_cover coefficient vector
# theta = c(beta_psi, beta_p, beta_pos, log_disp) (beta_p / beta_pos each packed as
# the site-level block then the optional visit-level block, exactly as the fitter
# stacks them). Weak Gaussian priors N(0, sigma.beta^2) on every coefficient and a
# broad N(0, sigma.logdisp^2) on log_disp keep the ridge / dispersion proper
# without materially shifting the data-dominated optimum. Returns list(lp, grad)
# over the packed coordinates. This is the oracle the C++ FullGradFn mirrors.
.tobs_occu_cover_nuts_logpost <- function(theta, model, sigma.beta = 5,
                                          sigma.logdisp = 5, field = NULL,
                                          re = NULL) {
  pin         <- model$process_info
  p_occ       <- pin[[1L]]$p
  p_det_site  <- ncol(model$X_det_site)
  p_det_visit <- if (!is.null(model$X_det_visit)) ncol(model$X_det_visit) else 0L
  p_pos_site  <- ncol(model$X_pos_site)
  p_pos_visit <- if (!is.null(model$X_pos_visit)) ncol(model$X_pos_visit) else 0L
  p_p   <- p_det_site + p_det_visit
  p_pos <- p_pos_site + p_pos_visit
  n_coef <- p_occ + p_p + p_pos
  total  <- n_coef + 1L

  # Optional coupled field(s): each block's z over its own trailing slice of
  # theta, entering psi weighted by its per-site design column and the cover arm
  # scaled by its own copy amplitude alpha. sigma / rho / alpha are sampled
  # coordinates or pinned, as each block declares. The R oracle mirrors the C++
  # FullGradFn byte-for-byte; absent, the non-spatial target.
  fvs       <- .ochf_views(field, total)
  has_field <- length(fvs) > 0L
  # Per-block state, the per-site loading on psi, and the coordinates the field
  # blocks occupy in total (so any trailing block -- the RE blocks -- starts
  # where the C++ layout starts it).
  fss     <- lapply(fvs, .ochf_forward, theta = theta)
  f_site  <- lapply(seq_along(fvs), function(b) .ochf_site_value(fvs[[b]], fss[[b]]))
  n_field <- sum(vapply(fvs, `[[`, numeric(1), "n_coord"))
  psi_off <- if (has_field) Reduce(`+`, f_site) else numeric(model$n_sites)
  pos_off <- if (has_field)
    Reduce(`+`, lapply(seq_along(fvs),
                       function(b) fss[[b]]$alpha$value * f_site[[b]]))
    else numeric(model$n_sites)

  bo         <- theta[seq_len(p_occ)]
  bp_site    <- theta[p_occ + seq_len(p_det_site)]
  bp_visit   <- if (p_det_visit > 0L) theta[p_occ + p_det_site + seq_len(p_det_visit)] else numeric(0)
  bpos_site  <- theta[p_occ + p_p + seq_len(p_pos_site)]
  bpos_visit <- if (p_pos_visit > 0L) theta[p_occ + p_p + p_pos_site + seq_len(p_pos_visit)] else numeric(0)
  log_disp   <- theta[total]
  disp       <- exp(log_disp)
  pos_code   <- .occu_cover_pos_code(model$positive)

  N <- model$n_sites; J <- model$max_visits
  sgm <- function(e) 1 / (1 + exp(-e))

  eta_psi <- as.numeric(model$X_occ %*% bo) + psi_off
  psi     <- sgm(eta_psi)

  eta_p <- matrix(as.numeric(model$X_det_site %*% bp_site), N, J)   # site block broadcast
  if (p_det_visit > 0L) {
    eta_p <- eta_p + matrix(as.numeric(model$X_det_visit %*% bp_visit), N, J, byrow = TRUE)
  }
  eta_pos <- matrix(as.numeric(model$X_pos_site %*% bpos_site) + pos_off, N, J)
  if (p_pos_visit > 0L) {
    eta_pos <- eta_pos + matrix(as.numeric(model$X_pos_visit %*% bpos_visit), N, J, byrow = TRUE)
  }

  # Observation-arm random intercepts. Each block owns n_groups whitened z
  # coordinates and its own SAMPLED log_sigma_re, laid out back to back AFTER
  # the field block, and adds the per-row offset sigma_re * z[code] to its arm.
  # `codes` is site-major over the padded grid (row (i - 1) * J + v), so
  # `matrix(..., byrow = TRUE)` puts each code on its own (site, visit); a 0
  # code carries no effect.
  re_view <- .ocre_view(re, total + n_field)
  for (rb in re_view) {
    off <- matrix(.ocre_offset(rb, theta), N, J, byrow = TRUE)
    if (rb$arm == 1L) eta_p   <- eta_p   + off
    else              eta_pos <- eta_pos + off
  }

  valid <- model$valid; y <- model$y; y_pos <- model$y_pos

  g_eta_psi <- numeric(N)
  g_eta_p   <- matrix(0, N, J)
  g_eta_pos <- matrix(0, N, J)
  g_logdisp <- 0
  lp_data   <- 0

  for (i in seq_len(N)) {
    vv <- which(valid[i, ])
    if (length(vv) == 0L) next
    any_det <- any(y[i, vv] == 1L)
    if (any_det) {
      lp_data <- lp_data + log(psi[i])
      g_eta_psi[i] <- 1 - psi[i]
      for (v in vv) {
        pv <- sgm(eta_p[i, v])
        if (y[i, v] == 1L) { lp_data <- lp_data + log(pv);     g_eta_p[i, v] <- 1 - pv }
        else               { lp_data <- lp_data + log(1 - pv); g_eta_p[i, v] <- -pv }
      }
      # Cover factor at detected visits with an observed cover; a missing (NA)
      # cover drops out (missing-at-random cover), the detection loop above still
      # counts the visit.
      for (v in vv[y[i, vv] == 1L & is.finite(y_pos[i, vv])]) {
        ev <- eta_pos[i, v]; yy <- y_pos[i, v]
        if (pos_code == 3L) {            # beta
          mu <- sgm(ev); a <- mu * disp; b <- (1 - mu) * disp
          # `.tobs_log_safe`, matching src/tobs_math.h::log_safe, which the
          # C++ target this oracle is asserted byte-for-byte against uses. A
          # raw log() made the two differ by -Inf against a finite value at a
          # cover of exactly 0 or 1, so an FD gradient check on boundary data
          # compared a different function than the sampler ran.
          ly <- .tobs_log_safe(yy); l1my <- .tobs_log_safe(1 - yy)
          lp_data <- lp_data + lgamma(disp) - lgamma(a) - lgamma(b) +
                     (a - 1) * ly + (b - 1) * l1my
          g_eta_pos[i, v] <- disp * mu * (1 - mu) *
                             (-digamma(a) + digamma(b) + ly - l1my)
          g_logdisp <- g_logdisp + disp * (digamma(disp) - mu * digamma(a) -
                       (1 - mu) * digamma(b) + mu * ly + (1 - mu) * l1my)
        } else if (pos_code == 4L) {     # identity-Gaussian (#112): raw response
          sig <- disp; r <- (yy - ev) / sig
          lp_data <- lp_data - .tobs_log_safe(sig) -
                     0.5 * log(2 * pi) - 0.5 * r * r
          g_eta_pos[i, v] <- r / sig
          g_logdisp <- g_logdisp + (r * r - 1)
        } else {                         # lognormal: Gaussian on log(cover)
          sig <- disp; lyy <- .tobs_log_safe(yy); r <- (lyy - ev) / sig
          lp_data <- lp_data - lyy - .tobs_log_safe(sig) -
                     0.5 * log(2 * pi) - 0.5 * r * r
          g_eta_pos[i, v] <- r / sig
          g_logdisp <- g_logdisp + (r * r - 1)
        }
      }
    } else {
      P0 <- exp(sum(log(1 - sgm(eta_p[i, vv]))))
      L  <- psi[i] * P0 + (1 - psi[i])
      lp_data <- lp_data + log(L)
      g_eta_psi[i] <- -psi[i] * (1 - psi[i]) * (1 - P0) / L
      for (v in vv) g_eta_p[i, v] <- -psi[i] * P0 * sgm(eta_p[i, v]) / L
    }
  }

  # Design-sandwich the eta-gradients onto the coefficient blocks. The visit-level
  # long vector follows site-major order (row i visit v -> (i-1)*J + v), matching
  # the byrow = TRUE design broadcast and the C++ row map i * J + v.
  grad_bo  <- as.numeric(crossprod(model$X_occ, g_eta_psi))
  grad_bp  <- as.numeric(crossprod(model$X_det_site, rowSums(g_eta_p)))
  if (p_det_visit > 0L)
    grad_bp <- c(grad_bp, as.numeric(crossprod(model$X_det_visit, as.numeric(t(g_eta_p)))))
  grad_bpos <- as.numeric(crossprod(model$X_pos_site, rowSums(g_eta_pos)))
  if (p_pos_visit > 0L)
    grad_bpos <- c(grad_bpos, as.numeric(crossprod(model$X_pos_visit, as.numeric(t(g_eta_pos)))))
  grad <- c(grad_bo, grad_bp, grad_bpos, g_logdisp)

  # Field blocks: per block accumulate the per-cell field score -- the psi-arm
  # score plus alpha times the cover-arm score, both weighted by that block's own
  # per-site design column -- and its copy amplitude's own data-score, then hand
  # both to the field backward pass. Each block's gradient (whitened field, then
  # any sampled hyper) appends in block order, as the C++ layout lays them out.
  lp_field <- 0
  if (has_field) {
    g_pos_site <- rowSums(g_eta_pos)
    for (b in seq_along(fvs)) {
      fv <- fvs[[b]]
      w  <- fv$weight %||% rep(1, N)
      g_f_site <- w * (g_eta_psi + fss[[b]]$alpha$value * g_pos_site)
      g_field  <- numeric(fv$n_units)
      for (s in seq_len(N)) {
        cell <- fv$field_map[s]
        g_field[cell] <- g_field[cell] + g_f_site[s]
      }
      bk       <- .ochf_backward(fv, fss[[b]], g_field,
                                 sum(f_site[[b]] * g_pos_site))
      grad     <- c(grad, bk$grad)
      lp_field <- lp_field + bk$lp
    }
  }

  # RE blocks trail the field block, each scoring against its own arm's per-row
  # eta-gradient (site-major, matching `codes`).
  lp_re <- 0
  for (rb in re_view) {
    g_row <- if (rb$arm == 1L) as.numeric(t(g_eta_p)) else as.numeric(t(g_eta_pos))
    bk    <- .ocre_backward(rb, theta, g_row)
    grad  <- c(grad, bk$grad)
    lp_re <- lp_re + bk$lp
  }

  lp  <- lp_data + lp_field + lp_re
  ib2 <- 1 / sigma.beta^2
  nb  <- n_coef
  bv  <- theta[seq_len(nb)]
  lp  <- lp - 0.5 * ib2 * sum(bv^2)
  grad[seq_len(nb)] <- grad[seq_len(nb)] - ib2 * bv
  ild2 <- 1 / sigma.logdisp^2
  lp   <- lp - 0.5 * ild2 * log_disp^2
  grad[total] <- grad[total] - ild2 * log_disp
  list(lp = lp, grad = grad)
}

# Byte-exact R oracle for the C++ FullGradFn: the exact three-level marginal
# log-likelihood plus weak Gaussian coefficient priors, with the same value the
# kernel returns. Reuses .occu_mscale_cover_nonspatial_ll (the Laplace path's LL) for
# the log-density and a finite-difference-free analytic gradient is supplied by
# the C++ side; here we numerically check the value and gradient via the C++
# joint-logpost cross-check, so the oracle only needs the (penalised) log-density.
#
# `sigma.logdisp` keeps log_disp proper, matching the two sibling NUTS targets
# (occu_cover_nuts.cpp, cover_nuts.cpp), which both sample it under this prior.
.tobs_occu_mscale_cover_nuts_logpost <- function(theta, model, idx, sigma.beta = 5,
                                                 sigma.logdisp = 5) {
  ll  <- .occu_mscale_cover_nonspatial_ll(theta, model, idx)
  ib2 <- 1 / sigma.beta^2
  beta_idx <- setdiff(seq_len(idx$total), idx$disp)
  ild2 <- 1 / sigma.logdisp^2
  ll - 0.5 * ib2 * sum(theta[beta_idx]^2) - 0.5 * ild2 * theta[idx$disp]^2
}

# Joint log-posterior + gradient of the removal coefficient vector. `marg` is a
# `.tobs_removal_nuts_marginal()` object; `lay` the shared abun layout. Weak
# Gaussian priors keep the ridge off the boundary. Returns list(lp, grad). This
# R version is the oracle the C++ FullGradFn (cpp_removal_nuts_joint_logpost) is
# cross-checked against.
.tobs_removal_nuts_logpost <- function(theta, marg, lay,
                                       sigma.beta = 10, sigma.logr = 1.5) {
  beta_lambda <- theta[lay$lambda]
  beta_p      <- theta[lay$p]
  r <- if (lay$is_nb) exp(min(theta[lay$log_r], .TOBS_LOG_R_MAX)) else Inf
  ev <- marg$eval_beta(beta_lambda, beta_p, r = r)

  lp   <- ev$log_lik
  grad <- c(as.numeric(crossprod(marg$X_lambda, ev$grad_eta_lambda)),
            as.numeric(crossprod(marg$X_p,      ev$grad_eta_p)))
  if (lay$is_nb) grad <- c(grad, sum(ev$grad_theta))

  ib2 <- 1 / sigma.beta^2
  lp  <- lp - 0.5 * ib2 * (sum(beta_lambda^2) + sum(beta_p^2))
  grad[lay$lambda] <- grad[lay$lambda] - ib2 * beta_lambda
  grad[lay$p]      <- grad[lay$p]      - ib2 * beta_p
  if (lay$is_nb) {
    lr  <- theta[lay$log_r]; ilr2 <- 1 / sigma.logr^2
    lp  <- lp - 0.5 * ilr2 * lr^2
    grad[lay$log_r] <- grad[lay$log_r] - ilr2 * lr
  }
  list(lp = lp, grad = grad)
}

