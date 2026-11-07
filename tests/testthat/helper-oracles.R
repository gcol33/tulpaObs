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

# Joint log-posterior + gradient of a single-species count-marginal coefficient
# vector (abun, removal, distance): the abundance arm, the detection arm `det`
# (`"p"`, or `"sigma"` for distance), and, where the layout carries them, the
# hazard log-shape and the NB log r, each scalar entering as a one-by-one arm.
# Weak Gaussian priors: sigma.beta on the coefficients, sigma.shape on the
# log-shape, sigma.logr on log r. `marg` is the family's marginal object
# (`eval_beta()` plus the arm designs); `lay` its layout. Returns list(lp, grad),
# the value the family's C++ FullGradFn is cross-checked against.
.tobs_count_nuts_logpost <- function(theta, marg, lay, det = "p",
                                     sigma.beta = 10, sigma.shape = 1.5,
                                     sigma.logr = 1.5) {
  hazard <- isTRUE(lay$hazard)
  is_nb  <- isTRUE(lay$is_nb)
  extra  <- list(r = if (is_nb) exp(min(theta[lay$log_r], .TOBS_LOG_R_MAX)) else Inf)
  if (hazard) extra$eta_b <- theta[lay$log_shape]
  ev <- do.call(marg$eval_beta, c(list(theta[lay$lambda], theta[lay[[det]]]), extra))
  arms <- list(
    list(idx = lay$lambda, X = marg$X_lambda, grad = "grad_eta_lambda"),
    list(idx = lay[[det]], X = marg[[paste0("X_", det)]],
         grad = paste0("grad_eta_", det)))
  prior_sd <- rep(sigma.beta, length(theta))
  one <- matrix(1, 1L, 1L)
  if (hazard) {
    arms <- c(arms, list(list(idx = lay$log_shape, X = one, grad = "grad_eta_b")))
    prior_sd[lay$log_shape] <- sigma.shape
  }
  if (is_nb) {
    ev$grad_log_r <- sum(ev$grad_theta)
    arms <- c(arms, list(list(idx = lay$log_r, X = one, grad = "grad_log_r")))
    prior_sd[lay$log_r] <- sigma.logr
  }
  .tobs_nuts_logpost_k(theta, ev, arms, length(theta), prior_sd)
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

# Community NUTS oracles. Each mirrors its family's C++ FullGradFn
# (src/ms_*_nuts.cpp) and is the value the C++ target is cross-checked against.
# NON-CENTERED parameterisation: the per-species block holds standard-normal
# z_s and the deviation is b_{s,arm} = C_arm z_{s,arm}, so the community
# covariance leaves the b-prior (z ~ N(0, I), no b' Sigma^{-1} b and no
# -0.5 S log|Sigma|) and enters only the data term through b = C z. This breaks
# the centered b <-> Sigma funnel that otherwise saturates the NUTS treedepth.
# The shared skeleton is .community_nc_logpost(); each family supplies its
# per-species likelihood and gradient.

# Joint log-posterior + gradient of a non-centred community NUTS target: the
# shared skeleton every community oracle below fills. Each species s carries a
# standard-normal block z_s (theta[b_idx(s)]); its coefficient vector is
# b_s = mu + C z_s arm by arm (C_arm the log-Cholesky factor of Sigma_arm,
# exp() of the coordinate for a scalar arm), so each community covariance
# enters only the data term. `mu_idx` gives the theta positions of the
# community means, in species-block order; `arms` lists each arm's positions
# `idx` within the block, its log-Cholesky coordinates `chol` and dimension
# `Pa`; `z_idx` all standard-normal coordinates. `sp_llgrad(s, b, grad)`
# returns the species' data log-likelihood `ll` at b and, when `grad`, its
# gradient `grad_b` in b and `grad_other` (a full-length vector, or NULL) for
# any coordinates outside the species blocks (shared globals, per-species
# nuisances). `mu_prior_sd` is the community-mean prior SD per mean (Inf for
# none); `extra_prior(theta)` returns list(lp, grad) for the coordinates the
# skeleton does not price. Returns list(lp, grad).
.community_nc_logpost <- function(theta, mu_idx, b_idx, arms, z_idx, sp_llgrad,
                                  n_species, priors, mu_prior_sd,
                                  extra_prior = NULL, grad = TRUE) {
  mu <- theta[mu_idx]
  P  <- length(mu)
  g    <- numeric(length(theta))
  g_mu <- numeric(P)
  lp   <- 0
  Cs <- lapply(arms, function(a) .ms_ocs_chol_unpack(theta[a$chol], a$Pa))
  As <- lapply(arms, function(a) matrix(0, a$Pa, a$Pa))
  for (s in seq_len(n_species)) {
    bidx <- b_idx(s)
    z_s  <- theta[bidx]
    b    <- mu
    for (k in seq_along(arms)) {
      ix <- arms[[k]]$idx
      b[ix] <- mu[ix] + as.numeric(Cs[[k]] %*% z_s[ix])
    }
    r  <- sp_llgrad(s, b, grad)
    lp <- lp + r$ll
    if (grad) {
      g_mu <- g_mu + r$grad_b
      for (k in seq_along(arms)) {
        ix <- arms[[k]]$idx
        gk <- r$grad_b[ix]
        g[bidx[ix]] <- g[bidx[ix]] + as.numeric(crossprod(Cs[[k]], gk))
        As[[k]] <- As[[k]] + outer(gk, z_s[ix])
      }
      if (!is.null(r$grad_other)) g <- g + r$grad_other
    }
  }
  z_all <- theta[z_idx]
  lp <- lp - 0.5 * sum(z_all^2)
  if (grad) g[z_idx] <- g[z_idx] - z_all
  for (k in seq_along(arms)) {
    a  <- arms[[k]]
    pr <- .ms_ocs_chol_logprior(theta[a$chol], a$Pa, priors)
    lp <- lp + pr$lp
    if (grad) g[a$chol] <- .ms_abun_nuts_chol_data_grad(As[[k]], Cs[[k]], a$Pa) + pr$grad
  }
  ip <- ifelse(is.finite(mu_prior_sd), 1 / mu_prior_sd^2, 0)
  lp <- lp - 0.5 * sum(ip * mu^2)
  if (grad) g[mu_idx] <- g[mu_idx] + g_mu - ip * mu
  if (!is.null(extra_prior)) {
    ep <- extra_prior(theta)
    lp <- lp + ep$lp
    if (grad) g <- g + ep$grad
  }
  if (!grad) return(list(lp = lp))
  list(lp = lp, grad = g)
}

.tobs_ms_abun_nuts_logpost <- function(theta, margs, lay, priors,
                                       sigma.beta = 10, sigma.logr = 1.5,
                                       grad = TRUE) {
  is_nb <- lay$is_nb
  arms <- c(list(list(idx = lay$lambda, chol = lay$chol_lam, Pa = lay$p_lam),
                 list(idx = lay$p,      chol = lay$chol_p,   Pa = lay$p_p)),
            if (is_nb) list(list(idx = lay$logr, chol = lay$chol_logr, Pa = 1L)))
  sd <- rep(Inf, lay$P)
  sd[c(lay$lambda, lay$p)] <- sigma.beta
  if (is_nb) sd[lay$logr] <- sigma.logr
  sp <- function(s, b, grad) {
    r  <- if (is_nb) exp(min(b[lay$logr], .TOBS_LOG_R_MAX)) else Inf
    ev <- margs[[s]]$eval_beta(b[lay$lambda], b[lay$p], r = r)
    if (!grad) return(list(ll = ev$log_lik))
    gb <- numeric(lay$P)
    gb[lay$lambda] <- as.numeric(crossprod(margs[[s]]$X_lambda, ev$grad_eta_lambda))
    gb[lay$p]      <- as.numeric(crossprod(margs[[s]]$X_p,      ev$grad_eta_p))
    if (is_nb) gb[lay$logr] <- sum(ev$grad_theta)
    list(ll = ev$log_lik, grad_b = gb)
  }
  .community_nc_logpost(theta, lay$mu, function(s) .ms_ocs_b_idx(lay, s), arms,
                        lay$b_off + seq_len(lay$n_species * lay$P), sp,
                        lay$n_species, priors, sd, grad = grad)
}

.tobs_ms_count_nuts_logpost <- function(theta, X, Y, lay, priors,
                                        sigma.beta = 10, sigma.logr = 1.5,
                                        grad = TRUE) {
  is_nb <- lay$is_nb; is_gauss <- lay$is_gauss
  is_bern <- identical(lay$family, "bernoulli")
  arms <- c(list(list(idx = lay$beta, chol = lay$chol_beta, Pa = lay$p_beta)),
            if (is_nb) list(list(idx = lay$logr, chol = lay$chol_logr, Pa = 1L)))
  sd <- rep(Inf, lay$P)
  sd[lay$beta] <- sigma.beta
  if (is_nb) sd[lay$logr] <- sigma.logr
  logphi <- if (is_gauss) theta[lay$logphi] else NULL
  sp <- function(s, b, grad) {
    # A missing (NA) site x species observation drops out of that species.
    obs <- !is.na(Y[, s])
    Xs  <- X[obs, , drop = FALSE]
    ys  <- Y[obs, s]
    eta <- as.numeric(Xs %*% b[lay$beta])
    gb  <- numeric(lay$P); go <- NULL
    if (is_nb) {
      r   <- exp(min(b[lay$logr], .TOBS_LOG_R_MAX))
      muv <- pmax(exp(pmin(eta, 700)), 1e-10)
      ll  <- sum(stats::dnbinom(ys, size = r, mu = muv, log = TRUE))
      if (grad) {
        gb[lay$beta] <- as.numeric(crossprod(Xs, r * (ys - muv) / (r + muv)))
        gb[lay$logr] <- r * sum(digamma(ys + r) - digamma(r) + log(r / (r + muv)) +
                                  1 - (ys + r) / (r + muv))
      }
    } else if (is_gauss) {
      phi_s <- exp(min(logphi[s], 30))
      ll <- sum(stats::dnorm(ys, mean = eta, sd = sqrt(phi_s), log = TRUE))
      if (grad) {
        gb[lay$beta] <- as.numeric(crossprod(Xs, (ys - eta) / phi_s))
        go <- numeric(length(theta))
        go[lay$logphi[s]] <- sum(-0.5 + (ys - eta)^2 / (2 * phi_s))
      }
    } else if (is_bern) {
      # jsdm(): log p = y eta - log(1 + exp(eta)), through plogis(log.p) as the
      # C++ stable log-sum-exp branch does.
      ll <- sum(ifelse(ys > 0, stats::plogis(eta,  log.p = TRUE),
                               stats::plogis(-eta, log.p = TRUE)))
      if (grad) gb[lay$beta] <- as.numeric(crossprod(Xs, ys - stats::plogis(eta)))
    } else {
      lam <- exp(pmin(eta, 700))
      ll  <- sum(ys * eta - lam - lgamma(ys + 1))
      if (grad) gb[lay$beta] <- as.numeric(crossprod(Xs, ys - lam))
    }
    list(ll = ll, grad_b = gb, grad_other = go)
  }
  extra <- if (is_gauss) function(th) {
    lpm <- priors$logphi_mean; lps <- priors$logphi_sd
    gr <- numeric(length(th)); gr[lay$logphi] <- -(logphi - lpm) / lps^2
    list(lp = -0.5 * sum(((logphi - lpm) / lps)^2), grad = gr)
  }
  .community_nc_logpost(theta, lay$mu, function(s) .ms_ocs_b_idx(lay, s), arms,
                        lay$b_off + seq_len(lay$n_species * lay$P), sp,
                        lay$n_species, priors, sd, extra, grad = grad)
}

.tobs_ms_dyn_occu_nuts_logpost <- function(theta, X_psi1, X_p, X_gamma, X_eps,
                                           em_stats, n_sites, n_seasons, lay,
                                           priors, sigma.beta = 5, grad = TRUE) {
  global <- theta[lay$global]
  beta_gam <- global[lay$gam]; beta_eps <- global[lay$eps]
  gamma <- stats::plogis(as.numeric(X_gamma %*% beta_gam))
  eps   <- stats::plogis(as.numeric(X_eps   %*% beta_eps))
  arms <- list(list(idx = lay$psi1, chol = lay$chol_psi1, Pa = lay$p_psi1),
               list(idx = lay$p,    chol = lay$chol_p,    Pa = lay$p_p))
  sp <- function(s, b, grad) {
    psi1 <- stats::plogis(as.numeric(X_psi1 %*% b[lay$psi1]))
    p    <- stats::plogis(as.numeric(X_p    %*% b[lay$p]))
    em <- .ms_dyn_occu_emissions(p, em_stats[[s]]$nvalid, em_stats[[s]]$ndet)
    if (!grad) {
      return(list(ll = .ms_dyn_occu_fwd_ll_vec(psi1, gamma, eps, em, n_sites, n_seasons)))
    }
    fb <- .ms_dyn_occu_fb_vec(psi1, p, gamma, eps, em, em_stats[[s]]$nvalid,
                              em_stats[[s]]$ndet, n_sites, n_seasons)
    gb <- numeric(lay$P)
    gb[lay$psi1] <- as.numeric(crossprod(X_psi1, fb$w1 - psi1))
    gb[lay$p]    <- as.numeric(crossprod(X_p,    fb$p_score))
    # The shared transition arms take every species' score.
    go <- numeric(length(theta)); gg <- numeric(lay$G)
    gg[lay$gam] <- as.numeric(crossprod(X_gamma, fb$col_y - gamma * fb$col_n))
    gg[lay$eps] <- as.numeric(crossprod(X_eps,   fb$ext_y - eps   * fb$ext_n))
    go[lay$global] <- gg
    list(ll = fb$ll, grad_b = gb, grad_other = go)
  }
  ib2 <- 1 / sigma.beta^2
  extra <- function(th) {
    gr <- numeric(length(th)); gr[lay$global] <- -ib2 * global
    list(lp = -0.5 * ib2 * sum(global^2), grad = gr)
  }
  .community_nc_logpost(theta, lay$mu, function(s) .ms_ocs_b_idx(lay, s), arms,
                        lay$b_off + seq_len(lay$n_species * lay$P), sp,
                        lay$n_species, priors, rep(sigma.beta, lay$P), extra,
                        grad = grad)
}

.tobs_ms_int_occu_nuts_logpost <- function(theta, X_psi, X_p, summaries, lay,
                                           priors, sigma.beta = 5, grad = TRUE) {
  D <- lay$D
  arms <- c(list(list(idx = lay$psi, chol = lay$chol_psi, Pa = lay$P_psi)),
            lapply(seq_len(D), function(d)
              list(idx = lay$p[[d]], chol = lay$chol_p[[d]], Pa = lay$P_p[d])))
  sp <- function(s, b, grad) {
    eta_psi <- as.numeric(X_psi %*% b[lay$psi])
    eta_p   <- lapply(seq_len(D), function(d) as.numeric(X_p[[d]] %*% b[lay$p[[d]]]))
    ll <- .ms_int_occu_sp_ll(eta_psi, eta_p, summaries[[s]])
    if (!grad) return(list(ll = ll))
    list(ll = ll, grad_b = .ms_int_occu_sp_grad(eta_psi, eta_p, summaries[[s]],
                                                X_psi, X_p)[seq_len(lay$P)])
  }
  .community_nc_logpost(theta, lay$mu, function(s) .ms_ocs_b_idx(lay, s), arms,
                        lay$b_off + seq_len(lay$n_species * lay$P), sp,
                        lay$n_species, priors, rep(sigma.beta, lay$P), grad = grad)
}

.tobs_ms_occu_cover_nuts_logpost <- function(theta, views, lay, priors,
                                             sigma.beta = 5, grad = TRUE) {
  P <- lay$P
  log_disp <- theta[lay$log_disp]
  arms <- list(list(idx = lay$occ, chol = lay$chol_occ, Pa = lay$P_occ),
               list(idx = lay$p,   chol = lay$chol_p,   Pa = lay$P_p),
               list(idx = lay$pos, chol = lay$chol_pos, Pa = lay$P_pos))
  sp <- function(s, b, grad) {
    bo <- b[lay$occ]; bp <- b[lay$p]; bpos <- b[lay$pos]
    ll <- .occu_cover_sp_ll(views[[s]], bo, bp, bpos, log_disp)
    if (!grad) return(list(ll = ll))
    gv <- .occu_cover_sp_grad(views[[s]], bo, bp, bpos, log_disp)
    go <- numeric(length(theta)); go[lay$log_disp] <- gv[P + 1L]
    list(ll = ll, grad_b = gv[seq_len(P)], grad_other = go)
  }
  # The shared log-dispersion: a weakly-informative Normal on the log scale.
  ld_mean <- priors$log_disp_mean %||% log(0.5)
  ld_sd   <- priors$log_disp_sd   %||% 2.0
  extra <- function(th) {
    gr <- numeric(length(th)); gr[lay$log_disp] <- -(log_disp - ld_mean) / ld_sd^2
    list(lp = -0.5 * ((log_disp - ld_mean) / ld_sd)^2, grad = gr)
  }
  .community_nc_logpost(theta, lay$mu, function(s) .ms_ocs_b_idx(lay, s), arms,
                        lay$b_off + seq_len(lay$n_species * P), sp,
                        lay$n_species, priors, rep(sigma.beta, P), extra,
                        grad = grad)
}

.tobs_ms_occu_cover_re_disp_logpost <- function(theta, views, lay, priors,
                                                sigma.beta = 5, grad = TRUE) {
  Pc <- lay$P_coef
  arms <- list(list(idx = lay$occ, chol = lay$chol_occ, Pa = lay$P_occ),
               list(idx = lay$p,   chol = lay$chol_p,   Pa = lay$P_p),
               list(idx = lay$pos, chol = lay$chol_pos, Pa = lay$P_pos),
               list(idx = lay$ld,  chol = lay$chol_ld,  Pa = 1L))
  sp <- function(s, b, grad) {
    bo <- b[lay$occ]; bp <- b[lay$p]; bpos <- b[lay$pos]; ld <- b[lay$ld]
    ll <- .occu_cover_sp_ll(views[[s]], bo, bp, bpos, ld)
    if (!grad) return(list(ll = ll))
    gv <- .occu_cover_sp_grad(views[[s]], bo, bp, bpos, ld)
    gb <- numeric(lay$Pz)
    gb[seq_len(Pc)] <- gv[seq_len(Pc)]
    gb[lay$ld] <- gv[Pc + 1L]
    list(ll = ll, grad_b = gb)
  }
  # Coefficient means: the sigma.beta ridge; the log-dispersion mean: a
  # weakly-informative Normal on the log scale.
  ld_mean <- priors$log_disp_mean %||% log(0.5)
  ld_sd   <- priors$log_disp_sd   %||% 2.0
  mu_ld <- theta[lay$mu_ld]
  extra <- function(th) {
    gr <- numeric(length(th)); gr[lay$mu_ld] <- -(mu_ld - ld_mean) / ld_sd^2
    list(lp = -0.5 * ((mu_ld - ld_mean) / ld_sd)^2, grad = gr)
  }
  .community_nc_logpost(theta, c(lay$mu_coef, lay$mu_ld),
                        function(s) .tobs_ms_occu_cover_re_disp_b_idx(lay, s), arms,
                        lay$b_off + seq_len(lay$n_species * lay$Pz), sp,
                        lay$n_species, priors, c(rep(sigma.beta, Pc), Inf), extra,
                        grad = grad)
}

.tobs_ms_occu_nuts_logpost <- function(theta, X_psi, X_p, summaries, lay, priors,
                                       sigma.beta = 5, grad = TRUE) {
  arms <- list(list(idx = lay$psi, chol = lay$chol_psi, Pa = lay$p_psi),
               list(idx = lay$p,   chol = lay$chol_p,   Pa = lay$p_p))
  sp <- function(s, b, grad) {
    eta_psi <- as.numeric(X_psi %*% b[lay$psi])
    eta_p   <- as.numeric(X_p   %*% b[lay$p])
    ll <- .ms_int_occu_sp_ll(eta_psi, list(eta_p), summaries[[s]])
    if (!grad) return(list(ll = ll))
    gv <- .ms_int_occu_sp_grad(eta_psi, list(eta_p), summaries[[s]], X_psi, list(X_p))
    gb <- numeric(lay$P); gb[c(lay$psi, lay$p)] <- gv[c(lay$psi, lay$p)]
    list(ll = ll, grad_b = gb)
  }
  sd <- rep(Inf, lay$P); sd[c(lay$psi, lay$p)] <- sigma.beta
  .community_nc_logpost(theta, lay$mu, function(s) .ms_ocs_b_idx(lay, s), arms,
                        lay$b_off + seq_len(lay$n_species * lay$P), sp,
                        lay$n_species, priors, sd, grad = grad)
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


