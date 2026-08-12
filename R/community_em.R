# Generic family-agnostic community Laplace-EM engine.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Log-determinant of a symmetric matrix via Cholesky, with an eigenvalue
# fallback flooring eigenvalues at 1e-12 for near-singular input.
.tobs_cem_logdet <- function(M) {
  ch <- tryCatch(chol((M + t(M)) / 2), error = function(e) NULL)
  if (is.null(ch)) {
    ev <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    return(sum(log(pmax(ev, 1e-12))))
  }
  2 * sum(log(diag(ch)))
}

# Moore-Penrose pseudo-inverse fallback (avoids a hard MASS dependency for the
# rare singular block; symmetric PSD input).
.tobs_cem_ginv <- function(M) {
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  pos <- e$values > max(e$values) * 1e-10
  e$vectors[, pos, drop = FALSE] %*%
    (t(e$vectors[, pos, drop = FALSE]) / e$values[pos])
}

# Per-species BLUP deviations in long form: one row per (species, arm, term).
# `arms` is a named character vector mapping the arm label (name) to the
# `ms_community` blup field (value); fields absent from `cm` are skipped (e.g.
# the optional per-species log_r RE of an NB N-mixture). Shared ranef() body for
# every community family (ms_occu / ms_dyn_occu / ms_occu_cover / ms_nmix).
.tobs_ranef_ms_long <- function(cm, arms) {
  to_long <- function(B, arm) {
    sp <- rownames(B); tm <- colnames(B)
    data.frame(species = rep(sp, times = ncol(B)), arm = arm,
               term = rep(tm, each = nrow(B)),
               estimate = as.numeric(B), stringsAsFactors = FALSE)
  }
  blocks <- Map(function(field, arm) {
    B <- cm[[field]]
    if (is.null(B)) NULL else to_long(B, arm)
  }, arms, names(arms))
  out <- do.call(rbind, blocks)
  rownames(out) <- NULL
  out
}


# AGHQ variance-component debias, shared by every .tobs_community_em()
# consumer (generalized from R/ms_occu_cover.R's own bespoke copy,
# gcol33/tulpaObs#226 part 2). The EM's Sigma/Cinv carry the documented
# Laplace small-cluster attenuation; this integrates the EXACT per-species RE
# posterior by adaptive Gauss-Hermite quadrature at the EM's own per-species
# mode/curvature (b_list, Cinv_list) and recovers the raw second moment
# E[b b'] the EM's own M-step targets -- the quantity Sigma should equal at
# the true (non-attenuated) posterior. `arm_idx` is an arbitrary NAMED list
# (not hardcoded arm names): the debiased Sigma comes back keyed the same way,
# so this serves any family's arm layout. `sp_ll(s, theta, global)` is the
# SAME per-species callback .tobs_community_em() already threads through.
.tobs_cem_aghq_sigma <- function(sp_ll, mu, global, b_list, Cinv_list,
                                 Sinv, arm_idx, P, n_quad = 5L) {
  S  <- length(b_list)
  gh <- .ms_gh_quad(as.integer(n_quad))
  # Tensor grid of z in R^P (P small): rows = nodes, plus the product GH weight.
  grid <- as.matrix(do.call(expand.grid, rep(list(gh$x), P)))   # (n_quad^P) x P
  lw   <- as.matrix(do.call(expand.grid, rep(list(log(gh$w)), P)))
  log_w_gh <- rowSums(lw)                                        # product weight
  sqrt2 <- sqrt(2)

  Sigma_acc <- matrix(0, P, P)
  for (s in seq_len(S)) {
    bhat <- b_list[[s]]
    C_s  <- (Cinv_list[[s]] + t(Cinv_list[[s]])) / 2
    L_s  <- tryCatch(t(chol(C_s)),
                     error = function(e) t(chol(C_s + diag(1e-8, P))))
    # b_q = bhat + sqrt(2) L z_q ; AGHQ exp(+sum z^2) undoes the e^{-z^2} weight.
    B_q  <- sweep(grid %*% (sqrt2 * t(L_s)), 2L, bhat, "+")      # n_node x P
    logp <- numeric(nrow(B_q))
    for (q in seq_len(nrow(B_q))) {
      bq <- B_q[q, ]
      logp[q] <- sp_ll(s, mu + bq, global) -
                 0.5 * as.numeric(crossprod(bq, Sinv %*% bq)) +
                 sum(grid[q, ]^2) + log_w_gh[q]
    }
    logp <- logp - max(logp)
    w <- exp(logp); w <- w / sum(w)
    # Raw second moment about zero E[b b'] = sum_q w_q b_q b_q' (the EM M-step
    # quantity, since b ~ N(0, Sigma)).
    M <- matrix(0, P, P)
    for (q in seq_len(nrow(B_q))) M <- M + w[q] * tcrossprod(B_q[q, ])
    Sigma_acc <- Sigma_acc + M
  }
  Sigma_acc <- Sigma_acc / S

  # Keep the per-arm block-diagonal structure of the community covariance.
  lapply(arm_idx, function(idx) {
    blk <- Sigma_acc[idx, idx, drop = FALSE]
    blk <- (blk + t(blk)) / 2
    ev  <- eigen(blk, symmetric = TRUE)
    ev$values <- pmax(ev$values, 1e-4)
    ev$vectors %*% diag(ev$values, length(idx)) %*% t(ev$vectors)
  })
}

# Reproject Cov(b_s|y) at a NEW Sinv without re-deriving the Newton solve
# (gcol33/tulpaObs#226 part 2). Cinv_s = solve(Htt_s + Sinv); Htt_s is pure
# likelihood curvature (independent of Sigma) and is recoverable from the
# EM's own output as solve(Cinv_old_s) - Sinv_old. Exact, no new
# approximation -- validated to machine precision against a direct
# construction of the joint (mu, b_s) arrowhead precision (dev_notes probe,
# 2026-08-12; first applied in R/ms_occu_cover.R, commit `03b87ad`).
.tobs_cem_reproject_cinv <- function(Cinv_list, Sinv_old, Sinv_new) {
  lapply(Cinv_list, function(C_s) {
    Htt_s <- solve(C_s) - Sinv_old
    tryCatch(solve(Htt_s + Sinv_new),
             error = function(e) .tobs_cem_ginv(Htt_s + Sinv_new))
  })
}


# ---------------------------------------------------------------------------
# Community Laplace-EM engine
# ---------------------------------------------------------------------------

#' Generic community Laplace-EM engine
#'
#' @keywords internal
#' @param sp_info Optional `function(s, theta, global)` returning species `s`'s
#'   `(P + G) x (P + G)` marginal observed information at `c(theta, global)`.
#'   Defaults to `NULL`, which finite-differences `sp_grad` -- passing neither
#'   leaves a fit byte-identical. Supply it when the family's kernel already
#'   exposes the per-site marginal observed information (the N-mixture / distance
#'   Louis block): the FD path spends `2 (P + G)` full marginal sweeps per species
#'   per Newton step to rediscover it.
#' @param init_b,init_Sigma Optional warm starts for the per-species deviations
#'   and the per-arm community covariances. Both default to `NULL`, which is the
#'   cold start (`b_s = 0`, `Sigma = sigma_init^2 I`) -- passing neither leaves a
#'   fit byte-identical. A block-coordinate caller that re-enters this EM once per
#'   outer pass (R/community_latent.R) passes the previous pass's state, so the
#'   EM resumes from a near-converged point instead of rediscovering it; it is the
#'   dominant cost when the per-species likelihood is expensive (an N-mixture /
#'   distance marginal sums over the latent count).
#' @param re_aghq Debias `Sigma`/`Cinv` by adaptive Gauss-Hermite quadrature of
#'   the exact per-species RE posterior (gcol33/tulpaObs#226 part 2), gated to
#'   `P <= re_aghq_maxdim` (tensor AGHQ over the joint b couples the arms).
#'   Defaults `FALSE` -- every existing caller is byte-identical unless it
#'   opts in. `n_quad` sets the per-dimension node count.
.tobs_community_em <- function(S, P, arm_idx, sp_ll, sp_grad = NULL,
                               init_mu, init_global = numeric(0),
                               penalize_global = FALSE, sigma_beta = 5,
                               priors = NULL,
                               sigma_init = 0.3, max_iter = 200L, tol = 1e-4,
                               newton_max = 30L, verbose = TRUE,
                               sp_info = NULL,
                               init_b = NULL, init_Sigma = NULL,
                               re_aghq = FALSE, n_quad = 5L,
                               re_aghq_maxdim = 4L) {
  S          <- as.integer(S)
  P          <- as.integer(P)
  G          <- length(init_global)
  U          <- P + G
  newton_max <- as.integer(newton_max)
  max_iter   <- as.integer(max_iter)
  tol        <- as.numeric(tol)
  arm_names  <- names(arm_idx)

  # Prior active when priors is NULL or any not-FALSE value.
  priors_active <- isTRUE(is.null(priors)) || !isFALSE(priors)

  # U x U diagonal prior precision. The mu-block is penalized when priors are
  # active; the global block only where penalize_global is TRUE (recycled).
  P_u <- matrix(0, U, U)
  if (priors_active) {
    diag(P_u)[seq_len(P)] <- 1 / (sigma_beta^2)
    if (G > 0L) {
      pen_g <- rep(as.logical(penalize_global), length.out = G)
      diag(P_u)[P + seq_len(G)][pen_g] <- 1 / (sigma_beta^2)
    }
  }

  glob_seq <- seq_len(G)            # empty integer vector when G == 0

  # Block-diagonal P x P inverse community precision from per-arm Sigma.
  blockdiag_inv <- function(Sig) {
    out <- matrix(0, P, P)
    for (arm in arm_names) {
      idx <- arm_idx[[arm]]
      out[idx, idx] <- solve(Sig[[arm]])
    }
    out
  }

  # Per-species gradient of the marginal log-likelihood wrt c(theta, global)
  # (length P + G). Supplied analytically, or central finite difference of
  # sp_ll with step h = 1e-5.
  sp_grad_fn <- if (!is.null(sp_grad)) {
    function(s, theta, global) sp_grad(s, theta, global)
  } else {
    function(s, theta, global) {
      u <- c(theta, global); h <- 1e-5
      g <- numeric(U)
      for (k in seq_len(U)) {
        up <- u; up[k] <- up[k] + h
        dn <- u; dn[k] <- dn[k] - h
        fp <- sp_ll(s, up[seq_len(P)], up[P + glob_seq])
        fm <- sp_ll(s, dn[seq_len(P)], dn[P + glob_seq])
        g[k] <- (fp - fm) / (2 * h)
      }
      g
    }
  }

  # Per-species observed information (U x U). Supplied analytically, or by
  # central finite difference of the gradient (step h = 1e-4), symmetrized:
  # -0.5 (J + J'). The FD path costs 2U gradient evaluations per species per
  # Newton step, which dominates when the per-species marginal is expensive (an
  # N-mixture / distance marginal sums over the latent count); a family whose
  # kernel already exposes the per-site marginal observed information should pass
  # `sp_info` and skip it.
  sp_info_fn <- if (!is.null(sp_info)) {
    function(s, theta, global) sp_info(s, theta, global)
  } else {
    function(s, theta, global) {
      u <- c(theta, global); h <- 1e-4
      J <- matrix(0, U, U)
      for (k in seq_len(U)) {
        up <- u; up[k] <- up[k] + h
        dn <- u; dn[k] <- dn[k] - h
        gp <- sp_grad_fn(s, up[seq_len(P)], up[P + glob_seq])
        gm <- sp_grad_fn(s, dn[seq_len(P)], dn[P + glob_seq])
        J[, k] <- (gp - gm) / (2 * h)
      }
      -0.5 * (J + t(J))
    }
  }

  # Penalized objective: sum of per-species log-likelihoods minus the RE and
  # fixed-effect Gaussian penalties.
  total_F <- function(mu, global, b_list, Sinv) {
    ll <- 0
    for (s in seq_len(S)) ll <- ll + sp_ll(s, mu + b_list[[s]], global)
    pen_b <- 0
    for (s in seq_len(S)) {
      bs <- b_list[[s]]
      pen_b <- pen_b + as.numeric(crossprod(bs, Sinv %*% bs))
    }
    u <- c(mu, global)
    ll - 0.5 * pen_b - 0.5 * as.numeric(crossprod(u, P_u %*% u))
  }

  # Laplace approximation to the integrated (marginal) log-likelihood at fixed
  # Sigma -- the quantity the EM increases monotonically. Per species the RE
  # integral is exp(ll_s) N(b; 0, Sigma) integrated by Laplace at the mode:
  #   logML_s = ll_s - 0.5 b' Sinv b + 0.5 log|Cinv_s|,
  # with Cinv_s = (H_tt_s + Sinv)^-1 the per-species posterior covariance.
  compute_logML <- function(mu, global, b_list, Cinv_list, Sigma, Sinv) {
    logdet_Sigma <- 0
    for (arm in arm_names) logdet_Sigma <- logdet_Sigma + .tobs_cem_logdet(Sigma[[arm]])
    acc <- 0
    for (s in seq_len(S)) {
      bs <- b_list[[s]]
      acc <- acc + sp_ll(s, mu + bs, global) -
        0.5 * as.numeric(crossprod(bs, Sinv %*% bs)) +
        0.5 * .tobs_cem_logdet(Cinv_list[[s]])
    }
    u <- c(mu, global)
    acc - 0.5 * S * logdet_Sigma -
      0.5 * as.numeric(crossprod(u, P_u %*% u))
  }

  # One joint-Newton mode-find of (mu, global, {b_s}) at fixed Sigma. Returns the
  # updated mode, per-species posterior covariances Cinv (Cov(b_s|y)), the
  # u-b_s cross-Hessian blocks Bf (u=(mu,global); gcol33/tulpaObs#226 -- needed
  # to draw (u, b_s) jointly rather than independently: Cov(u,b_s) =
  # -Vf %*% Bf_s %*% Cinv_s, and conditional on a draw of u, b_s's mean shifts
  # by -Cinv_s %*% t(Bf_s) %*% (u_draw - u_hat) while its covariance stays
  # exactly Cinv_s), and the marginal fixed-effect information Sf (Schur
  # complement of the b-block, Vf = solve(Sf)).
  solve_mode <- function(mu, global, b_list, Sinv) {
    F_cur <- total_F(mu, global, b_list, Sinv)
    Cinv_list <- vector("list", S)
    Sf <- NULL
    for (it in seq_len(newton_max)) {
      sumg_theta <- numeric(P); sumg_glob <- numeric(G)
      A11  <- matrix(0, P, P)
      A_tg <- matrix(0, P, G)
      A_gg <- matrix(0, G, G)
      Bf_list <- vector("list", S); gb_list <- vector("list", S)
      for (s in seq_len(S)) {
        theta <- mu + b_list[[s]]
        g    <- sp_grad_fn(s, theta, global)
        info <- sp_info_fn(s, theta, global)
        Htt <- info[seq_len(P), seq_len(P), drop = FALSE]
        Htg <- info[seq_len(P), P + glob_seq, drop = FALSE]      # P x G
        Hgg <- info[P + glob_seq, P + glob_seq, drop = FALSE]    # G x G
        g_theta <- g[seq_len(P)]
        g_glob  <- g[P + glob_seq]                               # length G
        C_s  <- Htt + Sinv
        Cinv <- tryCatch(solve(C_s), error = function(e) .tobs_cem_ginv(C_s))
        Cinv_list[[s]] <- Cinv
        Bf_list[[s]]   <- rbind(Htt, t(Htg))                     # (P+G) x P
        gb_list[[s]]   <- g_theta - Sinv %*% b_list[[s]]
        sumg_theta <- sumg_theta + g_theta
        sumg_glob  <- sumg_glob + g_glob
        A11  <- A11 + Htt
        A_tg <- A_tg + Htg
        A_gg <- A_gg + Hgg
      }
      u  <- c(mu, global)
      gf <- c(sumg_theta, sumg_glob) - as.numeric(P_u %*% u)     # length U
      if (G == 0L) {
        A_uu <- A11 + P_u
      } else {
        A_uu <- rbind(cbind(A11, A_tg), cbind(t(A_tg), A_gg)) + P_u
      }
      Sf <- A_uu; rhs <- gf
      for (s in seq_len(S)) {
        M   <- Bf_list[[s]] %*% Cinv_list[[s]]
        Sf  <- Sf  - M %*% t(Bf_list[[s]])
        rhs <- rhs - as.numeric(M %*% gb_list[[s]])
      }
      du <- tryCatch(solve(Sf, rhs), error = function(e) {
        solve(Sf + diag(1e-6, U), rhs)
      })
      db <- lapply(seq_len(S), function(s) {
        as.numeric(Cinv_list[[s]] %*% (gb_list[[s]] -
                     t(Bf_list[[s]]) %*% du))
      })
      # Backtracking line search on the penalized objective.
      step <- 1; ok <- FALSE
      for (ls in 1:25) {
        u_new    <- u + step * du
        mu_n     <- u_new[seq_len(P)]
        global_n <- u_new[P + glob_seq]
        b_n      <- lapply(seq_len(S), function(s) b_list[[s]] + step * db[[s]])
        F_n      <- total_F(mu_n, global_n, b_n, Sinv)
        if (is.finite(F_n) && F_n >= F_cur - 1e-8) { ok <- TRUE; break }
        step <- step / 2
      }
      if (!ok) break
      delta <- max(abs(c(step * du, unlist(db) * step)))
      mu <- mu_n; global <- global_n; b_list <- b_n; F_cur <- F_n
      if (delta < 1e-7) break
    }
    list(mu = mu, global = global, b_list = b_list, Cinv = Cinv_list,
         Bf = Bf_list, Sf = Sf, F = F_cur)
  }

  # ---- initialization ----
  # A warm start (init_b / init_Sigma) resumes a block-coordinate caller's
  # previous outer pass; absent, the cold start below is unchanged.
  mu     <- as.numeric(init_mu)
  global <- as.numeric(init_global)
  b_list <- if (is.null(init_b)) replicate(S, numeric(P), simplify = FALSE)
            else init_b
  if (is.null(init_Sigma)) {
    Sigma <- vector("list", length(arm_idx))
    names(Sigma) <- arm_names
    for (arm in arm_names) {
      k <- length(arm_idx[[arm]])
      Sigma[[arm]] <- diag(sigma_init^2, k)
    }
  } else {
    Sigma <- init_Sigma[arm_names]
  }

  # ---- EM loop ----
  converged <- FALSE; n_iter <- 0L; logML_prev <- -Inf
  # Progress + ETA for the community EM iterations (gcol33/tulpaObs#43); ON by
  # default, reusing tulpa's shared reporter so the heartbeat file matches every
  # other fitting loop. ETA is the upper bound to max_iter, finalised on convergence.
  .prog <- tulpa:::.tulpa_iter_progress("community-em", max_iter, unit = "iter")
  for (em in seq_len(max_iter)) {
    n_iter <- em
    Sinv <- blockdiag_inv(Sigma)
    res  <- solve_mode(mu, global, b_list, Sinv)
    mu <- res$mu; global <- res$global; b_list <- res$b_list

    logML <- compute_logML(mu, global, b_list, res$Cinv, Sigma, Sinv)

    # M-step: closed-form community covariance per arm.
    for (arm in arm_names) {
      idx <- arm_idx[[arm]]
      k   <- length(idx)
      acc <- matrix(0, k, k)
      for (s in seq_len(S)) {
        bs  <- b_list[[s]][idx]
        cov <- res$Cinv[[s]][idx, idx, drop = FALSE]
        acc <- acc + tcrossprod(bs) + cov
      }
      acc <- acc / S
      acc <- (acc + t(acc)) / 2
      ev  <- eigen(acc, symmetric = TRUE)
      ev$values <- pmax(ev$values, 1e-4)          # floor off the singular boundary
      Sigma[[arm]] <- ev$vectors %*% diag(ev$values, k) %*% t(ev$vectors)
    }

    rel <- abs(logML - logML_prev) / (abs(logML_prev) + 1)
    .prog$tick()
    if (isTRUE(verbose)) {
      message(sprintf("[community EM %d] logML=%.4f  rel_change=%.2e",
                      em, logML, rel))
    }
    if (em > 1L && rel < tol) { converged <- TRUE; break }
    logML_prev <- logML
  }
  .prog$finish()

  # ---- final marginal fixed-effect information -> community-mean covariance ----
  Sinv <- blockdiag_inv(Sigma)
  res  <- solve_mode(mu, global, b_list, Sinv)
  mu <- res$mu; global <- res$global; b_list <- res$b_list
  Vf <- tryCatch(solve(res$Sf), error = function(e) .tobs_cem_ginv(res$Sf))
  Vf <- (Vf + t(Vf)) / 2
  logML <- compute_logML(mu, global, b_list, res$Cinv, Sigma, Sinv)

  # AGHQ variance-component debias (gcol33/tulpaObs#226 part 2), opt-in
  # (re_aghq defaults FALSE -- every existing caller is byte-identical unless
  # it explicitly asks for this). Sigma carries the documented Laplace
  # small-cluster attenuation; Cinv must be reprojected at the debiased Sinv
  # or it silently keeps the pre-debias value even though Sigma changed (the
  # exact bug found and fixed in R/ms_occu_cover.R, commit `03b87ad`).
  Cinv_out <- res$Cinv
  debias_method <- "none"
  if (isTRUE(re_aghq) && P <= as.integer(re_aghq_maxdim)) {
    aghq_out <- tryCatch({
      Sinv_old <- Sinv
      Sd <- .tobs_cem_aghq_sigma(sp_ll, mu, global, b_list, res$Cinv, Sinv,
                                 arm_idx, P, n_quad = n_quad)
      Sinv_new <- blockdiag_inv(Sd)
      list(Sigma = Sd,
           Cinv = .tobs_cem_reproject_cinv(res$Cinv, Sinv_old, Sinv_new),
           ok = TRUE)
    }, error = function(e) {
      warning("community EM AGHQ variance debias failed (",
              conditionMessage(e), "); reporting the EM (Laplace) variance ",
              "components.", call. = FALSE)
      list(Sigma = Sigma, Cinv = res$Cinv, ok = FALSE)
    })
    Sigma <- aghq_out$Sigma; Cinv_out <- aghq_out$Cinv
    debias_method <- if (aghq_out$ok) "aghq" else "none"
  }

  list(mu = mu, global = global, b_list = b_list, Sigma = Sigma,
       Cinv = Cinv_out, Bf = res$Bf, Vf = Vf, logML = logML,
       converged = converged, n_iter = n_iter, debias_method = debias_method)
}
