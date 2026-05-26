# =============================================================================
# Adaptive Gauss-Hermite quadrature (AGHQ) debias for occupancy random-effect
# variance components (gcol33/tulpaObs#11 follow-up).
#
# The variance-component EM in R/em_laplace_re.R integrates the random-effect
# block b by Laplace (Gaussian curvature at the mode). For binary occupancy that
# integral attenuates the variance components (sigma, and the correlation of a
# full block) toward zero at small per-group sample size -- the lme4 glmer
# nAGQ=1 regime. This pass refines the fit by replacing the per-group Laplace
# integral with Q-point ADAPTIVE Gauss-Hermite quadrature: as Q grows the
# per-group marginal becomes exact, so the attenuation is removed. It is the
# variance-component analogue of .tobs_occu_marginal_refine() (which debiases the
# fixed effects on the no-RE single-season path) and reuses the same closed-form
# occupancy site marginal (z integrated exactly).
#
#   M_g(beta, Sigma) = int [ prod_{i in g} L_i(eta_i + Z_i b_g) ] N(b_g; 0, Sigma) db_g
#
# where L_i is the exact single-season site marginal
#   site with a detection : L_i = psi_i * p_i^{n_det} (1-p_i)^{J-n_det}
#   site with no detection: L_i = psi_i (1-p_i)^J + (1 - psi_i)
# and only psi_i = sigmoid(eta_i) depends on b_g. The marginal log-likelihood
# sum_g log M_g is maximized over (beta_occ, beta_det, chol(Sigma)); standard
# errors for the fixed effects are read from the exact-marginal Hessian.
#
# Scope. The per-group integral factorizes only under a SINGLE grouping factor,
# so this pass applies when every RE term shares one group index (one or more
# terms combined into a per-group block of dimension <= 3): the recovery-tested
# forms (1 | g), (x || g), (0 + x | g), (1 + x || g), (1 + x | g). Crossed /
# nested groupings (distinct group indices) do not factorize -- the pass returns
# NULL and the caller keeps the EM result. n_quad = 1 is the plain Laplace
# (nAGQ = 1) marginal; higher n_quad refines it toward the exact marginal. The
# fit is the marginal ML, so it does not reproduce the EM's REML-style variance
# update exactly (they agree to Laplace-level, ~a few percent). Any failure
# falls back to the EM result via tryCatch at the call site.
#
# Measured recovery (dev_notes/probe_re_aghq*.R, 10-15 seeds): at per-group
# n = 8 the EM attenuates sigma by ~18% (bias -0.16 at truth 0.9); AGHQ cuts
# that to ~4% (bias -0.04), matching NUTS. Correlated-slope sigmas recover to
# within ~1% on average.
#
# Correlation regularization. Pure ML has no prior on Sigma, so a weakly
# identified RE correlation can run to the +-1 boundary in a single fit. An
# LKJ(eta) penalty on each correlated block's correlation matrix R -- log-density
# (eta - 1) log det R, maximized at independence -- pulls the correlation off the
# boundary toward 0 without touching the marginal SDs (the Sigma = diag(sd) R
# diag(sd) decomposition). It is O(1) against the O(n_groups) likelihood, so it
# regularizes only weakly-identified correlations and is dominated by informative
# data; eta = 1 disables it; uncorrelated / intercept blocks have R = I and are
# unpenalized. The fit is then a MAP, and the reported Hessian is the posterior
# precision. Default eta = 1.5: on the recovery sim (dev_notes/_eta_sweep.R,
# truth rho = 0.61) it removes every +-1 boundary hit while keeping rho
# near-unbiased (bias -0.00 at per-group n = 12, +0.02 at n = 25), whereas
# eta = 2 over-shrinks (bias -0.07) and eta = 1 leaves boundary hits.
# =============================================================================


# Gauss-Hermite nodes/weights for the physicists' weight exp(-x^2) on the real
# line: int f(x) exp(-x^2) dx ~= sum_q w_q f(x_q). Golub-Welsch -- eigenvalues of
# the symmetric tridiagonal Jacobi matrix (zero diagonal, off-diagonal
# sqrt(k/2)) are the nodes; weights are sqrt(pi) times the squared first
# component of each normalized eigenvector. Native, no dependency.
.gauss_hermite <- function(n) {
  n <- as.integer(n)
  if (n <= 1L) return(list(x = 0, w = sqrt(pi)))
  k <- seq_len(n - 1L)
  off <- sqrt(k / 2)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1L)] <- off
  J[cbind(k + 1L, k)] <- off
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(x = e$values[ord], w = (sqrt(pi) * e$vectors[1, ]^2)[ord])
}


# AGHQ refinement of an RE occupancy fit. `design` is the per-term RE design
# from .tobs_re_design() (each element tagged with its `arm`: "occ" if the
# random effect enters the occupancy predictor, "det" if it enters detection);
# `beta_occ` / `beta_det` / `Sigma_list` / `b` the converged EM estimates (on
# the same scaled-X / natural-Z coordinates the EM used). The random effect b
# enters psi (occ arm) or p (det arm); the per-group marginal and its
# eta-derivatives branch on that, the rest of the quadrature is shared. Returns
# a list of refined estimates (point estimates, BLUPs, per-group BLUP variances,
# marginal fixed-effect SEs, refreshed occupancy weights, and the `arm`) or NULL
# when the pass does not apply (caller keeps the EM result) -- including when the
# RE is split across both arms (no single-factor factorization here).
.tobs_re_aghq <- function(model, design, beta_occ, beta_det,
                          Sigma_list, b, n_quad = 9L, lkj_eta = 1.5) {
  # ---- applicability: one arm, one shared grouping factor, low RE dim ----
  arm <- unique(vapply(design, function(d) d$arm %||% "occ", character(1)))
  if (length(arm) != 1L) return(NULL)
  idx1 <- as.integer(design[[1]]$idx)
  ng   <- as.integer(design[[1]]$n_groups)
  one_group <- all(vapply(design, function(d)
    identical(as.integer(d$idx), idx1) &&
      identical(as.integer(d$n_groups), ng), logical(1)))
  if (!one_group) return(NULL)
  if (!is.null(model$X_det_visit)) return(NULL)
  nc_terms <- vapply(design, function(d) as.integer(d$n_coefs), integer(1))
  dtot <- sum(nc_terms)
  if (dtot > 3L) return(NULL)

  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  y <- model$y
  p_occ <- ncol(X_occ); p_det <- ncol(X_det)

  # ---- per-site detection summary (fixed across the optimization) ----
  vmask   <- y >= 0
  n_valid <- rowSums(vmask)
  Y <- y; Y[!vmask] <- 0
  n_det   <- rowSums(Y)
  any_det <- n_det > 0
  keep    <- n_valid > 0

  # Combined per-group design Z (N x dtot), term-major columns, and the offset
  # of each term within a group's coefficient vector.
  Zc <- do.call(cbind, lapply(design, function(d) d$Z))
  coef_off <- cumsum(c(0L, nc_terms))
  rows_by_g <- lapply(seq_len(ng), function(g) which(idx1 == g & keep))

  # ---- RE covariance parametrization: per-term Cholesky ----
  # Unconstrained params per term: log of the Cholesky diagonal (nc), plus the
  # free lower off-diagonal (nc*(nc-1)/2) for a correlated term. An uncorrelated
  # term keeps Sigma diagonal (off-diagonal forced to 0), preserving the `||`
  # semantics.
  term_meta <- lapply(design, function(d) {
    nc <- as.integer(d$n_coefs)
    corr <- isTRUE(d$correlated) && nc > 1L
    list(nc = nc, corr = corr, noff = if (corr) nc * (nc - 1L) / 2L else 0L)
  })
  re_par0 <- numeric(0)
  for (k in seq_along(design)) {
    nc <- term_meta[[k]]$nc
    Sk <- Sigma_list[[k]]
    if (term_meta[[k]]$corr) {
      L <- t(chol(Sk + diag(1e-8, nc)))
      re_par0 <- c(re_par0, log(pmax(diag(L), 1e-3)), L[lower.tri(L)])
    } else {
      re_par0 <- c(re_par0, 0.5 * log(pmax(diag(Sk), 1e-3)))
    }
  }

  rebuild_sigma <- function(re_par) {
    out <- vector("list", length(design)); pos <- 0L
    for (k in seq_along(design)) {
      nc <- term_meta[[k]]$nc
      dge <- exp(pmin(re_par[pos + seq_len(nc)], 5)); pos <- pos + nc
      if (term_meta[[k]]$corr) {
        noff <- term_meta[[k]]$noff
        L <- diag(dge, nc); L[lower.tri(L)] <- re_par[pos + seq_len(noff)]
        pos <- pos + noff
        out[[k]] <- tcrossprod(L)
      } else {
        out[[k]] <- diag(dge^2, nc)
      }
    }
    out
  }

  block_diag <- function(mats) {
    if (length(mats) == 1L) return(mats[[1]])
    out <- matrix(0, dtot, dtot); pos <- 0L
    for (m in mats) {
      nc <- nrow(m); out[pos + seq_len(nc), pos + seq_len(nc)] <- m; pos <- pos + nc
    }
    out
  }

  # ---- quadrature tensor grid for dimension dtot ----
  gh <- .gauss_hermite(n_quad)
  Q <- length(gh$x)
  grid <- as.matrix(expand.grid(rep(list(seq_len(Q)), dtot)))
  Tnodes <- matrix(gh$x[grid], nrow = nrow(grid), ncol = dtot)
  logw_q <- rowSums(matrix(log(gh$w)[grid], nrow = nrow(grid)))
  t2_q   <- rowSums(Tnodes^2)

  cl <- function(e) pmin(pmax(e, -30), 30)

  # ---- arm-specific per-site marginal --------------------------------------
  # The random effect b enters the RE-arm linear predictor; the OTHER arm's
  # contribution is fixed per site (`fx`). `re_eta_base(bo, bd)` is the RE arm's
  # fixed predictor (the part without b): X_occ %*% beta_occ when b moves psi
  # (occ arm), X_det %*% beta_det when b moves p (det arm).
  #
  #   occ arm: psi = sigmoid(eta), p fixed.
  #     detected      L_i = psi p^d (1-p)^{n-d}  -> log L = log psi + dterm
  #     non-detected  L_i = 1 - a psi,  a = 1 - (1-p)^n
  #   det arm: p = sigmoid(eta), psi fixed.
  #     detected      log L = log psi + d log p + (n-d) log(1-p)
  #     non-detected  L_i = psi (1-p)^n + (1 - psi)
  if (arm == "occ") {
    re_eta_base <- function(bo, bd) cl(as.numeric(X_occ %*% bo))
    make_fixed  <- function(bo, bd) {              # detection per-site, no b
      p <- plogis(cl(as.numeric(X_det %*% bd)))
      list(dterm = ifelse(any_det, n_det * log(p) + (n_valid - n_det) * log1p(-p), 0),
           a_nodet = ifelse(any_det, 0, 1 - exp(n_valid * log1p(-p))))
    }
    # Per-site logL + eta-derivatives at the mode (eta is the occ predictor).
    site_deriv <- function(rows, eta, fx) {
      s <- plogis(eta); ad <- any_det[rows]
      logL <- d1 <- d2 <- numeric(length(rows))
      logL[ad] <- log(s[ad]) + fx$dterm[rows][ad]
      d1[ad]   <- 1 - s[ad]
      d2[ad]   <- -s[ad] * (1 - s[ad])
      nd <- !ad
      aa <- fx$a_nodet[rows][nd]; sn <- s[nd]; sp <- sn * (1 - sn)
      g_ <- 1 - aa * sn; Nn <- -aa * sp
      logL[nd] <- log(g_)
      d1[nd] <- Nn / g_
      d2[nd] <- (-aa * sp * (1 - 2 * sn) * g_ - Nn * (-aa * sp)) / g_^2
      list(logL = logL, d1 = d1, d2 = d2)
    }
    # Vectorized logL over the quadrature node matrix ETA (rows x nodes).
    site_lmat <- function(rows, ETA, fx) {
      S <- plogis(ETA); ad <- any_det[rows]
      out <- matrix(0, length(rows), ncol(ETA))
      if (any(ad))  out[ad, ]  <- log(S[ad, , drop = FALSE]) + fx$dterm[rows][ad]
      if (any(!ad)) out[!ad, ] <- log(1 - fx$a_nodet[rows][!ad] * S[!ad, , drop = FALSE])
      out
    }
    # Refreshed P(z = 1 | y): psi carries b, (1-p)^n is fixed.
    w_refresh <- function(b_out, bo, bd, fx) {
      psi <- plogis(cl(re_eta_base(bo, bd) + .tobs_re_offset(design, b_out)))
      prod0 <- ifelse(any_det, 0, 1 - fx$a_nodet)   # (1-p)^n
      ifelse(any_det, 1, psi * prod0 / (psi * prod0 + (1 - psi)))
    }
  } else {  # arm == "det"
    re_eta_base <- function(bo, bd) cl(as.numeric(X_det %*% bd))
    make_fixed  <- function(bo, bd) list(psi = plogis(cl(as.numeric(X_occ %*% bo))))
    site_deriv <- function(rows, eta, fx) {
      p <- plogis(eta); ad <- any_det[rows]
      nv <- n_valid[rows]; nd <- n_det[rows]; psi <- fx$psi[rows]
      logL <- d1 <- d2 <- numeric(length(rows))
      # detected: binomial-in-p; psi enters only as an additive log constant.
      logL[ad] <- log(psi[ad]) + nd[ad] * log(p[ad]) + (nv[ad] - nd[ad]) * log1p(-p[ad])
      d1[ad]   <- nd[ad] - nv[ad] * p[ad]
      d2[ad]   <- -nv[ad] * p[ad] * (1 - p[ad])
      # non-detected: L = psi q + (1 - psi), q = (1-p)^n, dq/deta = -n p q.
      o <- !ad
      pn <- p[o]; nvn <- nv[o]; psn <- psi[o]
      q <- exp(nvn * log1p(-pn)); L <- psn * q + (1 - psn)
      A  <- -psn * nvn * pn * q                      # = psi * dq/deta
      dA <- -psn * nvn * pn * q * ((1 - pn) - nvn * pn)
      logL[o] <- log(L)
      d1[o] <- A / L
      d2[o] <- (dA * L - A^2) / L^2
      list(logL = logL, d1 = d1, d2 = d2)
    }
    site_lmat <- function(rows, ETA, fx) {
      P <- plogis(ETA); ad <- any_det[rows]
      nv <- n_valid[rows]; nd <- n_det[rows]; psi <- fx$psi[rows]
      out <- matrix(0, length(rows), ncol(ETA))
      if (any(ad)) {
        out[ad, ] <- log(psi[ad]) + nd[ad] * log(P[ad, , drop = FALSE]) +
          (nv[ad] - nd[ad]) * log1p(-P[ad, , drop = FALSE])
      }
      if (any(!ad)) {
        Q <- exp(nv[!ad] * log1p(-P[!ad, , drop = FALSE]))   # (1-p)^n
        out[!ad, ] <- log(psi[!ad] * Q + (1 - psi[!ad]))
      }
      out
    }
    # Refreshed P(z = 1 | y): psi is fixed, (1-p)^n carries b.
    w_refresh <- function(b_out, bo, bd, fx) {
      psi <- fx$psi
      p   <- plogis(cl(re_eta_base(bo, bd) + .tobs_re_offset(design, b_out)))
      prod0 <- ifelse(any_det, 0, exp(n_valid * log1p(-p)))  # (1-p)^n
      ifelse(any_det, 1, psi * prod0 / (psi * prod0 + (1 - psi)))
    }
  }

  # Per-group posterior mode of b (damped Newton on the penalized integrand) and
  # the curvature -H there (positive definite). Arm-agnostic via site_deriv.
  grp_mode <- function(rows, a_g, Zg, P, fx) {
    bb <- numeric(dtot)
    for (it in seq_len(50L)) {
      gv <- site_deriv(rows, cl(a_g + as.numeric(Zg %*% bb)), fx)
      grad <- as.numeric(crossprod(Zg, gv$d1)) - as.numeric(P %*% bb)
      negH <- P - crossprod(Zg, gv$d2 * Zg)
      step <- tryCatch(solve(negH, grad), error = function(e) NULL)
      if (is.null(step)) break
      bb <- bb + step
      if (max(abs(step)) < 1e-9) break
    }
    gv <- site_deriv(rows, cl(a_g + as.numeric(Zg %*% bb)), fx)
    list(b = bb, negH = P - crossprod(Zg, gv$d2 * Zg))
  }

  # log M_g via adaptive GHQ centred at the mode, scaled by the mode curvature.
  grp_logM <- function(rows, a_g, Zg, P, logdetS, fx, want_post = FALSE) {
    m <- grp_mode(rows, a_g, Zg, P, fx)
    Lc <- tryCatch(t(chol(solve(m$negH))), error = function(e) NULL)  # C = Lc Lc'
    if (is.null(Lc)) return(NULL)
    B <- matrix(m$b, nrow(Tnodes), dtot, byrow = TRUE) + sqrt(2) * Tnodes %*% t(Lc)
    ETA <- cl(matrix(a_g, length(rows), nrow(B)) + Zg %*% t(B))
    logLmat <- site_lmat(rows, ETA, fx)
    hvals <- colSums(logLmat) - 0.5 * rowSums((B %*% P) * B)
    terms <- logw_q + hvals + t2_q
    mx <- max(terms)
    lse <- mx + log(sum(exp(terms - mx)))
    logM <- -0.5 * dtot * log(2 * pi) - 0.5 * logdetS +
            0.5 * dtot * log(2) + sum(log(diag(Lc))) + lse
    if (want_post) list(logM = logM, b = m$b, var = diag(solve(m$negH)))
    else logM
  }

  # LKJ(eta) log-prior on the correlated blocks: sum_k (eta - 1) log det R_k,
  # R_k the correlation matrix of Sigma_k. log det R_k = log det Sigma_k -
  # sum_c log Sigma_k[c, c]; -> 0 (no penalty) for an uncorrelated / intercept
  # block (R_k = I) and -> -Inf as a correlation approaches +-1.
  lkj_logprior <- function(Slist) {
    if (lkj_eta == 1) return(0)
    s <- 0
    for (k in seq_along(design)) {
      if (!term_meta[[k]]$corr) next
      Sk <- Slist[[k]]
      ld <- as.numeric(determinant(Sk, logarithm = TRUE)$modulus)
      s <- s + (lkj_eta - 1) * (ld - sum(log(pmax(diag(Sk), 1e-12))))
    }
    s
  }

  nll <- function(theta) {
    bo <- theta[seq_len(p_occ)]
    bd <- theta[p_occ + seq_len(p_det)]
    Slist <- rebuild_sigma(theta[-seq_len(p_occ + p_det)])
    S <- block_diag(Slist) + diag(1e-10, dtot)
    P <- tryCatch(solve(S), error = function(e) NULL)
    if (is.null(P)) return(1e10)
    logdetS <- as.numeric(determinant(S, logarithm = TRUE)$modulus)
    fx <- make_fixed(bo, bd)
    aRE <- re_eta_base(bo, bd)
    total <- 0
    for (g in seq_len(ng)) {
      rows <- rows_by_g[[g]]
      if (!length(rows)) next
      lm <- grp_logM(rows, aRE[rows], Zc[rows, , drop = FALSE], P, logdetS, fx)
      if (is.null(lm) || !is.finite(lm)) return(1e10)
      total <- total + lm
    }
    -(total + lkj_logprior(Slist))
  }

  start <- c(beta_occ, beta_det, re_par0)
  opt <- stats::optim(start, nll, method = "BFGS", hessian = TRUE,
                      control = list(maxit = 200L, reltol = 1e-9))

  V <- tryCatch(solve(opt$hessian), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(opt$par))) return(NULL)
  se_all <- sqrt(pmax(diag(V), 0))

  beta_occ_ref <- opt$par[seq_len(p_occ)]
  beta_det_ref <- opt$par[p_occ + seq_len(p_det)]
  Sigma_ref <- rebuild_sigma(opt$par[-seq_len(p_occ + p_det)])

  # ---- per-group BLUPs + variances at the optimum ----
  S <- block_diag(Sigma_ref) + diag(1e-10, dtot)
  P <- solve(S)
  logdetS <- as.numeric(determinant(S, logarithm = TRUE)$modulus)
  fx <- make_fixed(beta_occ_ref, beta_det_ref)
  aRE <- re_eta_base(beta_occ_ref, beta_det_ref)
  BHAT <- matrix(0, ng, dtot)
  BVAR <- matrix(rep(diag(S), each = ng), ng, dtot)  # empty groups -> prior var
  for (g in seq_len(ng)) {
    rows <- rows_by_g[[g]]
    if (!length(rows)) next
    post <- grp_logM(rows, aRE[rows], Zc[rows, , drop = FALSE], P, logdetS, fx,
                     want_post = TRUE)
    if (is.null(post)) return(NULL)
    BHAT[g, ] <- post$b
    BVAR[g, ] <- post$var
  }

  # Re-pack BLUPs / variances into the term-major, group-major (byrow) layout
  # .tobs_re_param_block() consumes (matching .tobs_re_design term order).
  b_out <- numeric(0); bvar_out <- numeric(0)
  for (k in seq_along(design)) {
    cols <- coef_off[k] + seq_len(nc_terms[k])
    b_out    <- c(b_out,    as.numeric(t(BHAT[, cols, drop = FALSE])))
    bvar_out <- c(bvar_out, as.numeric(t(BVAR[, cols, drop = FALSE])))
  }

  # Refreshed occupancy weights P(z = 1 | y) at the refined estimate, for
  # fitted() / residuals() consistency.
  w_ref <- w_refresh(b_out, beta_occ_ref, beta_det_ref, fx)

  list(
    ok           = TRUE,
    arm          = arm,
    beta_occ     = beta_occ_ref,
    beta_det     = beta_det_ref,
    Sigma_list   = Sigma_ref,
    b            = b_out,
    b_var        = bvar_out,
    beta_occ_se  = se_all[seq_len(p_occ)],
    det_se       = se_all[p_occ + seq_len(p_det)],
    weights      = w_ref,
    n_quad       = as.integer(n_quad),
    lkj_eta      = lkj_eta,
    converged    = isTRUE(opt$convergence == 0L)
  )
}
