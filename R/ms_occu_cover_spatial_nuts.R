# ms_occu_cover_spatial_nuts.R - full-vector joint log-posterior for the
# reduced-rank spatial-factor community occu_cover (gcol33/tulpa#67), the NUTS
# target density.
#
# The Laplace-EM fitter (ms_occu_cover_spatial.R) profiles the community
# covariances Sigma, the field precisions tau_w, and the field hyperparameters
# out by closed-form M-steps, then reports a Gaussian Laplace posterior over the
# inner latent c(mu, {b_s}, L, W, log_disp). NUTS instead samples EVERYTHING
# jointly -- the inner latent AND the hyperparameters -- from the exact joint
# posterior, which removes the Gaussian-approximation over-dispersion that makes
# the Laplace draws unusable for WAIC / LOO and gives calibrated hyperparameter
# intervals.
#
# The target factorises so the inner objective is reused verbatim. Write the
# joint log-posterior as
#
#   log p = [ data log-lik + log N(mu|0,sigma.beta) + log N(L|0,sd_L)
#             + log N(b|0,Sigma) quadratic + GMRF(W|tau) quadratic ]      (A)
#         + [ -0.5 S log|Sigma_arm|  (per arm) ]                          (B)
#         + [  0.5 rank_k log tau_w[k]  (per factor) ]                    (C)
#         + log p(Sigma coords) + log p(log tau) + log p(log_disp)        (D)
#
# Block (A) is EXACTLY .ms_ocs_penll_grad evaluated with the block-diagonal
# Sinv(Sigma) and tau_w derived from the sampled hyperparameters: the b- and
# field-quadratics, the mu / loading Gaussian priors, and the data likelihood all
# live there already, and its gradient over the inner coordinates is complete and
# FD-verified. The inner solve treats Sigma / tau as constants, so it omits the
# normalisers (B) and (C); those depend only on the hyperparameter coordinates and
# are added here with their analytic gradients. (D) are weakly-informative
# hyperpriors placed DIRECTLY on the sampled unconstrained coordinates (the
# log-Cholesky entries of Sigma, log tau_w, log_disp), so the target is the
# density in those coordinates with no further change of variables.
#
# The community covariance of an arm is parameterised by the lower-triangular
# Cholesky factor C (Sigma = C C') with a log-diagonal, the standard unconstrained
# map for a PD matrix; the field precisions by log tau_w. The icar field has no
# extra hyperparameter; the proper-CAR rho / BYM2 phi axes are a later increment.

# ---------------------------------------------------------------------------
# Log-Cholesky packing for a PD community covariance
# ---------------------------------------------------------------------------
#
# A P x P community covariance Sigma = C C' is carried by its lower-triangular
# Cholesky factor C, packed column-major over the lower triangle (column j
# contributes rows j..P) with the diagonal stored on the log scale (so the map
# is onto all of R^{P(P+1)/2} and the diagonal stays positive). The packed length
# is P(P+1)/2.

.ms_ocs_chol_dim <- function(P) as.integer(P * (P + 1L) / 2L)

# Packed positions (1-based) of the diagonal log-entries within the vector, one
# per column: column j's block starts after sum_{j'<j}(P-j'+1) entries and its
# first entry is the diagonal.
.ms_ocs_chol_diag_pos <- function(P) {
  pos <- integer(P); off <- 0L
  for (j in seq_len(P)) { pos[j] <- off + 1L; off <- off + (P - j + 1L) }
  pos
}

# Lower-triangular Cholesky factor C (diagonal positive) -> packed vector.
.ms_ocs_chol_pack <- function(C) {
  P <- nrow(C); out <- numeric(.ms_ocs_chol_dim(P)); pos <- 0L
  for (j in seq_len(P)) {
    out[pos + 1L] <- log(C[j, j]); pos <- pos + 1L
    if (j < P) {
      ni <- P - j
      out[pos + seq_len(ni)] <- C[(j + 1L):P, j]; pos <- pos + ni
    }
  }
  out
}

# Packed vector -> lower-triangular Cholesky factor C (diagonal = exp()).
.ms_ocs_chol_unpack <- function(vec, P) {
  C <- matrix(0, P, P); pos <- 0L
  for (j in seq_len(P)) {
    C[j, j] <- exp(vec[pos + 1L]); pos <- pos + 1L
    if (j < P) {
      ni <- P - j
      C[(j + 1L):P, j] <- vec[pos + seq_len(ni)]; pos <- pos + ni
    }
  }
  C
}


# ---------------------------------------------------------------------------
# Hyperprior specification
# ---------------------------------------------------------------------------
#
# Weakly-informative, proper priors placed directly on the sampled unconstrained
# coordinates. The Cholesky log-diagonal carries a Normal centred at log(0.5)
# (community SDs of order 0.5 on the link scale, the recovery-harness regime) with
# a wide SD; the Cholesky off-diagonals a mean-zero Normal (shrinking the implied
# community correlations toward independence); log tau_w and log_disp wide Normals.
# These only regularise the boundaries (near-singular Sigma, runaway field
# precision) and are swamped by the data at the sample sizes the family targets.
.ms_ocs_nuts_priors <- function() {
  list(chol_logdiag_mean = log(0.5), chol_logdiag_sd = 1.5,
       chol_offdiag_sd   = 1.0,
       log_tau_mean = 0,  log_tau_sd  = 2.0,
       log_disp_mean = log(0.5), log_disp_sd = 2.0,
       # Field hyperparameter (car_proper rho / bym2 phi) on the logit scale: a
       # weakly-informative Normal centred at logit(0.5) = 0 (no a-priori pull
       # toward weak or strong spatial structure), wide enough to span the unit
       # interval. The field shape, well-identified by the shared loadings,
       # dominates this at the family's sample sizes.
       logit_h_mean = 0, logit_h_sd = 1.5)
}

# log density (up to an additive constant) of one arm's Cholesky coordinates and
# its gradient w.r.t. the packed vector.
.ms_ocs_chol_logprior <- function(vec, P, priors) {
  dpos <- .ms_ocs_chol_diag_pos(P)
  is_d <- logical(length(vec)); is_d[dpos] <- TRUE
  sd_v <- ifelse(is_d, priors$chol_logdiag_sd, priors$chol_offdiag_sd)
  mu_v <- ifelse(is_d, priors$chol_logdiag_mean, 0)
  lp   <- -0.5 * sum(((vec - mu_v) / sd_v)^2)
  list(lp = lp, grad = -(vec - mu_v) / sd_v^2)
}


# ---------------------------------------------------------------------------
# Inner-parameter width and hyperparameter layout
# ---------------------------------------------------------------------------

# Length of the inner-latent prefix c(mu, {b_s}, L-block, [Lpos], W, log_disp),
# i.e. length(fit$par): the loading block is S*K unconstrained or the triangular
# free count when constrained.
.ms_ocs_npar_inner <- function(d, constrain) {
  L_width <- if (constrain) .ms_ocs_lfree_dim(d$S, d$K) else d$S * d$K
  d$P + d$S * d$P + L_width + d$Lpos_w + d$N * d$K + 1L
}

# Full NUTS coordinate layout: the inner prefix, then the three arm Cholesky
# blocks (occ, p, pos), then log tau_w (length K), then -- for a hyperparameterised
# field (proper-CAR rho / BYM2 phi) -- the logit field-hyperparameter block
# logit_h_w (length K). Returns the block offsets and the total length.
.ms_ocs_nuts_layout <- function(d, constrain, has_hyper = FALSE) {
  n_inner <- .ms_ocs_npar_inner(d, constrain)
  q_occ <- .ms_ocs_chol_dim(d$P_occ)
  q_p   <- .ms_ocs_chol_dim(d$P_p)
  q_pos <- .ms_ocs_chol_dim(d$P_pos)
  off <- n_inner
  occ <- off + seq_len(q_occ); off <- off + q_occ
  p   <- off + seq_len(q_p);   off <- off + q_p
  pos <- off + seq_len(q_pos); off <- off + q_pos
  tau <- off + seq_len(d$K);   off <- off + d$K
  logit_h <- integer(0)
  if (has_hyper) { logit_h <- off + seq_len(d$K); off <- off + d$K }
  list(n_inner = n_inner, inner = seq_len(n_inner),
       chol_occ = occ, chol_p = p, chol_pos = pos, log_tau = tau,
       logit_h = logit_h, has_hyper = has_hyper,
       total = off)
}


# ---------------------------------------------------------------------------
# Joint log-posterior + gradient (the NUTS target density)
# ---------------------------------------------------------------------------

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

# Field-hyperparameter derivatives at one factor: the quadratic derivative
# dquad = W_k'(dR/dh)W_k and the log-determinant derivative dlogdet = d log|R(h)|/dh,
# both closed form. car_proper: R = D - h A so dR/dh = -A, and log|R| = log|D| +
# sum log(1 - h gamma) so dlogdet = -sum gamma/(1 - h gamma). bym2: R = V diag(1/m)
# V' with m = (1-h) + h s, so W'R W = sum a_j/m_j (a = (V'W)^2) gives dquad =
# -sum a (s-1)/m^2, and log|R| = -sum log m gives dlogdet = -sum (s-1)/m.
.ms_ocs_field_hyper_terms <- function(spec, h, Wk) {
  if (identical(spec$type, "car_proper")) {
    AW    <- as.numeric(spec$A %*% Wk)
    dquad <- -as.numeric(Wk %*% AW)
    dlogdet <- -sum(spec$gamma / (1 - h * spec$gamma))
    return(list(dquad = dquad, dlogdet = dlogdet))
  }
  if (identical(spec$type, "bym2")) {
    Vw <- as.numeric(crossprod(spec$V, Wk))
    m  <- (1 - h) + h * spec$s
    a  <- Vw^2
    dquad   <- -sum(a * (spec$s - 1) / m^2)
    dlogdet <- -sum((spec$s - 1) / m)
    return(list(dquad = dquad, dlogdet = dlogdet))
  }
  list(dquad = 0, dlogdet = 0)
}

# d log p / d (packed Cholesky coords) for one arm, from the term
#   T = -0.5 tr(Sigma^{-1} M) - 0.5 S log|Sigma|
# plus the coordinate hyperprior. With Sigma = C C' and G = dT/dSigma (symmetric),
# dT/dC = 2 G C; the log-diagonal coords carry the extra chain factor C_jj.
.ms_ocs_chol_block_grad <- function(C, M, S, vec, P, priors) {
  Si <- chol2inv(t(C))
  G  <- 0.5 * (Si %*% M %*% Si) - 0.5 * S * Si          # dT/dSigma (symmetric)
  dC <- 2 * (G %*% C)                                   # dT/dC
  pr <- .ms_ocs_chol_logprior(vec, P, priors)
  out <- numeric(length(vec)); pos <- 0L
  for (j in seq_len(P)) {
    out[pos + 1L] <- dC[j, j] * C[j, j]; pos <- pos + 1L  # diag: log-link chain
    if (j < P) {
      ni <- P - j
      out[pos + seq_len(ni)] <- dC[(j + 1L):P, j]; pos <- pos + ni
    }
  }
  out + pr$grad
}

# Extract the per-species b list (length S, each length P) from the inner par,
# honouring the constrained vs unconstrained loading-block width (b sits in the
# shared mu + b prefix, ahead of the loadings, so the width does not matter, but
# the unpackers branch on it for clarity).
.ms_ocs_inner_b <- function(par_inner, d, constrain) {
  up <- if (constrain) .ms_ocs_unpack_c(par_inner, d) else .ms_ocs_unpack(par_inner, d)
  up$b
}

# Extract the N x K field matrix W from the inner par.
.ms_ocs_inner_W <- function(par_inner, d, constrain) {
  up <- if (constrain) .ms_ocs_unpack_c(par_inner, d) else .ms_ocs_unpack(par_inner, d)
  up$W
}


# ---------------------------------------------------------------------------
# C++ marginal-likelihood marshalling
# ---------------------------------------------------------------------------

# Flatten the bound model into the plain list the C++ marginal-likelihood entry
# (cpp_ms_ocs_marginal_ll, and the NUTS gradient_fn to come) reads: the cell-level
# designs, the per-species detection / cover / valid matrices, and the model
# dimensions. The shared field W is supplied per evaluation inside the packed
# parameter vector, not here.
.ms_ocs_nuts_spec <- function(model) {
  d <- .ms_ocs_dims(model)
  N <- model$n_sites; J <- model$max_visits
  per_sp <- function(arr) lapply(seq_len(d$S), function(s)
    matrix(as.double(arr[, , s]), N, J))
  fspec <- model$field_spec %||%
    .ms_ocs_field_spec(model$adj, model$field_type %||% "icar")
  out <- list(n_sites = N, max_visits = J, S = d$S, K = d$K,
       P_occ = d$P_occ, P_p = d$P_p, P_pos = d$P_pos,
       cover_factor = isTRUE(d$cover_factor),
       is_beta = identical(model$positive, "beta"),
       X_occ = model$X_occ, X_p = model$X_det_site, X_pos = model$X_pos_site,
       y = per_sp(model$y), y_pos = per_sp(model$y_pos),
       valid = per_sp(model$valid),
       Q = model$icar_Q, field_rank = fspec$rank,
       field_type = fspec$type, has_hyper = isTRUE(fspec$has_hyper))
  # Proper-field primitives for the C++ R(h) construction (car / bym2). The
  # icar path uses Q only; these are the minimal pieces the per-factor R(h),
  # log|R(h)|, and dR/dh evaluations need.
  if (identical(fspec$type, "car_proper")) {
    out$A <- fspec$A; out$deg <- fspec$deg
    out$gamma <- fspec$gamma; out$logdetD <- fspec$logdetD
  } else if (identical(fspec$type, "bym2")) {
    out$V <- fspec$V; out$s <- fspec$s
  }
  out
}


# ---------------------------------------------------------------------------
# NUTS fit (front-door method = "nuts")
# ---------------------------------------------------------------------------

# Inverse-mass diagonal for the NUTS warm-start: the posterior variance per
# coordinate, read off the finite-difference diagonal of the joint
# log-posterior Hessian at the mode (the Laplace metric). Uses the fast C++
# gradient so the K*(few) gradient evaluations are cheap.
.ms_ocs_nuts_metric <- function(spec, theta, pri, sigma.beta, sd_L, constrain,
                                h = 1e-4) {
  np <- length(theta); md <- numeric(np)
  g <- function(th) cpp_ms_ocs_joint_logpost(spec, th, pri, sigma.beta, sd_L,
                                             constrain)$grad
  for (j in seq_len(np)) {
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    md[j] <- -(g(tp)[j] - g(tm)[j]) / (2 * h)
  }
  1 / pmax(md, 1e-3)
}

# Reconstruct a Laplace-EM-shaped `fit` object from the NUTS draws so the shared
# builder (build_ms_occu_cover_spatial_fit) packages the tobs_fit unchanged. The
# inner-latent posterior is summarised by its mean (`par`) and covariance (`cov`)
# over the draws; the hyperparameters (community covariances, field precisions,
# field hyperparameter) by their posterior means. The shared builder's
# association / map posteriors then draw from N(par, cov) -- a moment-matched
# Gaussian to the NUTS inner-latent posterior, marginalising the loadings / field
# the same way the Laplace path does.
.ms_ocs_nuts_fit_from_draws <- function(model, res, em, constrain) {
  d   <- .ms_ocs_dims(model)
  spec_f <- model$field_spec %||%
    .ms_ocs_field_spec(model$adj, model$field_type %||% "icar")
  has_hyper <- isTRUE(spec_f$has_hyper)
  lay <- .ms_ocs_nuts_layout(d, constrain, has_hyper)
  draws <- res$draws
  inner <- draws[, lay$inner, drop = FALSE]

  par <- colMeans(inner)
  cov <- stats::cov(inner)
  up  <- if (constrain) .ms_ocs_unpack_c(par, d) else .ms_ocs_unpack(par, d)

  # Community covariances: posterior mean of Sigma_arm = C C' over the draws.
  arm_chol <- list(occ = lay$chol_occ, p = lay$chol_p, pos = lay$chol_pos)
  arm_P    <- list(occ = d$P_occ, p = d$P_p, pos = d$P_pos)
  Sigma <- lapply(c("occ", "p", "pos"), function(a) {
    cols <- arm_chol[[a]]; P <- arm_P[[a]]
    acc <- matrix(0, P, P)
    for (i in seq_len(nrow(draws))) {
      C <- .ms_ocs_chol_unpack(draws[i, cols], P); acc <- acc + tcrossprod(C)
    }
    acc <- acc / nrow(draws); (acc + t(acc)) / 2
  })
  names(Sigma) <- c("occ", "p", "pos")

  tau_w <- exp(colMeans(draws[, lay$log_tau, drop = FALSE]))
  field_hyper <- if (has_hyper)
    stats::plogis(colMeans(draws[, lay$logit_h, drop = FALSE])) else NULL

  list(par = par, cov = cov, mu = up$mu, b = up$b, ld = up$ld,
       L = up$L, w = up$W, Lpos = up$Lpos, constrained = constrain,
       Sigma = Sigma, tau_w = tau_w, field_hyper = field_hyper,
       field_type = spec_f$type, hyper_name = spec_f$hyper_name,
       sd_L = em$sd_L, em_logpen = em$em_logpen,
       convergence = 0L, d = d,
       nuts = list(draws = draws, accept_prob = res$accept_prob,
                   divergent = res$divergent, treedepth = res$treedepth,
                   epsilon = res$epsilon, layout = lay))
}

# Fit the spatial-factor community occu_cover by NUTS: a short Laplace-EM warm
# start sets the initial position and the inverse-mass metric, then tulpa's NUTS
# (driven by the C++ FullGradFn) samples the exact joint posterior. Returns the
# reconstructed fit list consumed by build_ms_occu_cover_spatial_fit.
.tobs_fit_ms_occu_cover_spatial_nuts <- function(model, sd_L = 1.0,
                                                 sigma.beta = 5, constrain = FALSE,
                                                 control = list()) {
  em <- .tobs_fit_ms_occu_cover_spatial(
    model, sd_L = sd_L, sigma.beta = sigma.beta, constrain = constrain,
    max.em = control[["max.iter"]] %||% 30L,
    tol = control[["tol"]] %||% 1e-3,
    verbose = isTRUE(control[["verbose"]]))

  theta0 <- .ms_ocs_nuts_pack_init(em)
  spec   <- .ms_ocs_nuts_spec(model)
  pri    <- .ms_ocs_nuts_priors()
  inv_metric <- .ms_ocs_nuts_metric(spec, theta0, pri, sigma.beta, sd_L,
                                    constrain)

  n_warmup <- control[["n.warmup"]] %||% 500L
  n_sample <- control[["n.iter"]]   %||% 1000L
  res <- cpp_ms_ocs_nuts(
    spec, theta0, pri, sigma.beta, sd_L, inv_metric,
    n_iter = as.integer(n_warmup + n_sample), n_warmup = as.integer(n_warmup),
    max_treedepth = as.integer(control[["max.treedepth"]] %||% 10L),
    adapt_delta = control[["adapt.delta"]] %||% 0.95,
    seed = as.integer(control[["seed"]] %||% 1L),
    verbose = isTRUE(control[["verbose"]]), constrain = constrain)

  .ms_ocs_nuts_fit_from_draws(model, res, em, constrain)
}


# ---------------------------------------------------------------------------
# Warm-start packing from an EM fit
# ---------------------------------------------------------------------------

# Pack the Laplace-EM solution into the full NUTS coordinate vector: the inner
# mode (fit$par, already in the constrained vs unconstrained layout the fit used),
# the arm covariances as log-Cholesky coordinates, log tau_w, and -- for a
# hyperparameterised field -- the logit of the EM field hyperparameter (rho/phi).
# This is the NUTS initial position (the chain starts at the EM mode) and the
# reference point for the FD gradient check.
.ms_ocs_nuts_pack_init <- function(fit) {
  chol_coords <- function(Sig) .ms_ocs_chol_pack(t(chol(Sig)))  # Sigma = C C'
  hyper <- if (!is.null(fit$field_hyper)) stats::qlogis(fit$field_hyper) else numeric(0)
  c(fit$par,
    chol_coords(fit$Sigma$occ),
    chol_coords(fit$Sigma$p),
    chol_coords(fit$Sigma$pos),
    log(fit$tau_w),
    hyper)
}
