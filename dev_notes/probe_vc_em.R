# Validate a deterministic variance-component EM for occupancy random effects
# (gcol33/tulpaObs#11) END-TO-END on simulated data, before wiring into the
# package. Method:
#   E-step:   psi_i = logit^{-1}(X_i b_occ + RE_offset_i); w_i = P(z_i=1|y).
#   M-step:   pseudo-binomial M-trick (y=round(wM), n=M) recovers the penalized
#             MAP of (beta, b) IF the RE prior is rescaled sigma_eff=sigma/sqrt(M)
#             so penalty and data term share the M scale (argmax invariant).
#   VC step:  recompute RE posterior covariance at NATURAL scale (n=1 working
#             weights) via Schur, update sigma^2_c = mean_g(b_gc^2 + Var(b_gc|y)).
suppressMessages(library(tulpa))

# --- RE posterior covariance (b-block of joint Hessian inverse), natural scale.
# X: N x p, Zlist: list of (idx, nc, Zmat[N x nc]); sig: list per term (len nc);
# eta: linear predictor at mode. Returns list of per-term diag(Cov_bb).
re_post_var <- function(X, eta, terms) {
  W <- as.numeric(plogis(eta) * (1 - plogis(eta)))   # n=1 Bernoulli info
  XtWX <- crossprod(X, W * X)
  # Build combined Z and D^{-1}.
  Zc <- list(); Dinv <- numeric(0); span <- list(); off <- 0L
  for (tm in terms) {
    nc <- tm$nc; ng <- tm$n_groups; idx <- tm$idx
    Zfull <- if (is.null(tm$Z)) matrix(1, nrow(X), 1) else tm$Z
    ii <- rep(seq_len(nrow(X)), each = nc)
    jj <- rep((idx - 1L) * nc, each = nc) + rep(seq_len(nc), nrow(X))
    Zc[[length(Zc) + 1L]] <- Matrix::sparseMatrix(i = ii, j = jj,
      x = as.numeric(t(Zfull)), dims = c(nrow(X), ng * nc))
    di <- numeric(ng * nc)
    for (c in seq_len(nc)) di[seq(c, ng * nc, by = nc)] <- 1 / (tm$sig[c]^2 + 1e-10)
    Dinv <- c(Dinv, di)
    span[[length(span) + 1L]] <- (off + 1L):(off + ng * nc); off <- off + ng * nc
  }
  Z <- do.call(cbind, Zc)
  ZtWZ <- as.matrix(Matrix::crossprod(Z, W * Z))
  XtWZ <- as.matrix(crossprod(X, W * Z))
  P_bb <- ZtWZ + diag(Dinv, length(Dinv))
  Schur <- P_bb - t(XtWZ) %*% solve(XtWX, XtWZ)
  Cov <- solve(Schur)
  dvar <- diag(Cov)
  lapply(span, function(s) dvar[s])
}

fit_vc_em <- function(y, X_occ, X_det, g, nc = 1L, Zmat = NULL,
                      M = 1000L, max_iter = 60L, tol = 1e-5) {
  N <- nrow(X_occ); ng <- max(g)
  n_valid <- rowSums(y >= 0); n_det <- rowSums(y == 1); any_det <- n_det > 0
  p_occ <- ncol(X_occ); p_det <- ncol(X_det)
  Zfull <- if (nc == 1L) matrix(1, N, 1) else Zmat
  beta_occ <- rep(0, p_occ); beta_det <- rep(0, p_det)
  b <- numeric(ng * nc); sig <- rep(0.7, nc)
  for (it in seq_len(max_iter)) {
    re_off <- vapply(seq_len(N), function(i)
      sum(Zfull[i, ] * b[(g[i] - 1L) * nc + seq_len(nc)]), numeric(1))
    psi <- plogis(as.numeric(X_occ %*% beta_occ) + re_off)
    p   <- plogis(as.numeric(X_det %*% beta_det))
    w <- numeric(N)
    for (i in seq_len(N)) {
      if (any_det[i]) w[i] <- 1
      else if (n_valid[i] == 0) w[i] <- psi[i]
      else { num <- psi[i] * (1 - p[i])^n_valid[i]; w[i] <- num / (num + 1 - psi[i]) }
    }
    # M-step occ via M-trick + rescaled sigma.
    y_occ <- pmin(pmax(ifelse(any_det, M, as.integer(round(w * M))), 0L), M)
    sig_eff <- sig / sqrt(M)
    re_list <- list(list(idx = g, n_groups = ng, sigma = sig_eff, n_coefs = nc,
                         Z = if (nc > 1L) Zfull else NULL))
    fo <- tulpa::tulpa_laplace(y = y_occ, n_trials = rep(M, N), X = X_occ,
                               re_list = re_list, family = "binomial")
    beta_new <- fo$mode[seq_len(p_occ)]; b_new <- fo$mode[-seq_len(p_occ)]
    # VC update at natural scale.
    eta_mode <- as.numeric(X_occ %*% beta_new) +
      vapply(seq_len(N), function(i)
        sum(Zfull[i, ] * b_new[(g[i] - 1L) * nc + seq_len(nc)]), numeric(1))
    dv <- re_post_var(X_occ, eta_mode,
                      list(list(idx = g, n_groups = ng, nc = nc,
                                Z = if (nc > 1L) Zfull else NULL, sig = sig)))[[1]]
    sig_new <- numeric(nc)
    for (c in seq_len(nc)) {
      bc <- b_new[seq(c, ng * nc, by = nc)]; vc <- dv[seq(c, ng * nc, by = nc)]
      sig_new[c] <- sqrt(mean(bc^2 + vc))
    }
    # Detection M-step (weighted binomial).
    w_det <- w; w_det[any_det] <- 1; keep <- n_valid > 0 & w_det > 1e-6
    fd <- tulpa::tulpa_laplace(y = n_det[keep], n_trials = n_valid[keep],
                               X = X_det[keep, , drop = FALSE],
                               weights = w_det[keep], family = "binomial")
    beta_det_new <- fd$mode[seq_len(p_det)]
    delta <- max(abs(c(beta_new - beta_occ, beta_det_new - beta_det,
                       sig_new - sig)))
    beta_occ <- beta_new; b <- b_new; sig <- sig_new; beta_det <- beta_det_new
    if (delta < tol) break
  }
  list(beta_occ = beta_occ, beta_det = beta_det, sigma = sig, b = b,
       n_iter = it, converged = delta < tol)
}

## ---- Test 1: iid intercept RE ----
set.seed(101)
ng <- 30L; per <- 25L; N <- ng * per; J <- 6L
g <- rep(seq_len(ng), each = per)
x <- rnorm(N)
sigma_true <- 0.9
b_true <- rnorm(ng, 0, sigma_true)
psi <- plogis(0.3 - 0.6 * x + b_true[g]); p_true <- 0.45
z <- rbinom(N, 1, psi)
y <- matrix(0L, N, J); for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p_true)
X_occ <- cbind(1, x); X_det <- matrix(1, N, 1)
f1 <- fit_vc_em(y, X_occ, X_det, g, nc = 1L)
cat("== iid intercept RE ==\n")
cat(sprintf("  beta_occ: %.3f %.3f (truth 0.30 -0.60)\n", f1$beta_occ[1], f1$beta_occ[2]))
cat(sprintf("  p: %.3f (truth %.2f)\n", plogis(f1$beta_det), p_true))
cat(sprintf("  sigma: %.3f (truth %.2f)\n", f1$sigma, sigma_true))
cat(sprintf("  cor(b_hat, b_true): %.3f   iters=%d conv=%s\n",
            cor(f1$b, b_true), f1$n_iter, f1$converged))

## ---- Test 2: uncorrelated random slope (1 + x || g) ----
set.seed(202)
ng <- 30L; per <- 30L; N <- ng * per; J <- 6L
g <- rep(seq_len(ng), each = per)
x <- rnorm(N)
sig_true <- c(0.7, 0.5)
b0 <- rnorm(ng, 0, sig_true[1]); b1 <- rnorm(ng, 0, sig_true[2])
psi <- plogis(0.2 + b0[g] + (-0.4 + b1[g]) * x); p_true <- 0.5
z <- rbinom(N, 1, psi)
y <- matrix(0L, N, J); for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p_true)
X_occ <- cbind(1, x); X_det <- matrix(1, N, 1); Z <- cbind(1, x)
f2 <- fit_vc_em(y, X_occ, X_det, g, nc = 2L, Zmat = Z)
b0h <- f2$b[seq(1, ng * 2, by = 2)]; b1h <- f2$b[seq(2, ng * 2, by = 2)]
cat("\n== uncorrelated slope RE ==\n")
cat(sprintf("  beta_occ: %.3f %.3f (truth 0.20 -0.40)\n", f2$beta_occ[1], f2$beta_occ[2]))
cat(sprintf("  p: %.3f (truth %.2f)\n", plogis(f2$beta_det), p_true))
cat(sprintf("  sigma: %.3f %.3f (truth 0.70 0.50)\n", f2$sigma[1], f2$sigma[2]))
cat(sprintf("  cor(b0): %.3f  cor(b1): %.3f  iters=%d conv=%s\n",
            cor(b0h, b0), cor(b1h, b1), f2$n_iter, f2$converged))

cat("\n=== done ===\n")
