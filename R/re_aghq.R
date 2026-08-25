# =============================================================================
# Adaptive Gauss-Hermite quadrature (AGHQ) debias for occupancy random-effect
# variance components ( follow-up).
#
# The variance-component EM in R/em_laplace_re.R integrates the random-effect
# block b by Laplace (Gaussian curvature at the mode). For binary occupancy that
# integral attenuates the variance components (sigma, and the correlation of a
# full block) toward zero at small per-group sample size -- the lme4 glmer
# nAGQ=1 regime. This pass refines the fit by replacing the per-group Laplace
# integral with adaptive Gauss-Hermite quadrature.
#
# DIVISION OF LABOUR. The generic AGHQ machinery -- the quadrature grid,
# per-group mode-finding, the log-Cholesky covariance parametrization, the LKJ
# correlation penalty, the joint (beta, Sigma) optimization, and the marginal
# Hessian -- lives in tulpa (`tulpa::tulpa_re_aghq()`), reused across families.
# tulpaObs owns only the FAMILY-SPECIFIC part: the closed-form occupancy site
# marginal (latent z integrated out) and its eta-derivatives, supplied to the
# engine as a `make_site` callback. This mirrors the EM-Laplace boundary
# (tulpa owns the engine, tulpaObs the likelihood callbacks).
#
#   M_g(beta, Sigma) = int [ prod_{i in g} L_i(eta_i + Z_i b_g) ] N(b_g; 0, Sigma) db_g
#
# The random effect b enters EITHER the occupancy or the detection predictor:
#
#   occ arm  (b moves psi, p fixed):
#     site with a detection : L_i = psi_i * p_i^{n_det} (1-p_i)^{J-n_det}
#     site with no detection: L_i = psi_i (1-p_i)^J + (1 - psi_i)
#   det arm  (b moves p, psi fixed):
#     site with a detection : log L_i = log psi_i + n_det log p_i + (J-n_det) log(1-p_i)
#     site with no detection: L_i = psi_i (1-p_i)^J + (1 - psi_i)
#
# `make_site()` returns, for the current (beta_occ, beta_det): the RE-arm fixed
# predictor `eta_re`, the per-site `deriv(rows, eta)` (logL + first/second
# eta-derivatives, for the per-group mode) and `lmat(rows, ETA)` (logL over the
# quadrature nodes). The detection-arm derivatives are the binomial-in-p form
# (FD-verified in dev_notes/probe_re_det_aghq_deriv.R).
#
# Scope. One shared grouping factor on one arm, RE dim <= 3 (the engine's
# quadrature grid is n_quad^dim). RE split across both arms or crossed / nested
# groupings -> the pass returns NULL and the caller keeps the EM result.
# n_quad = 1 is the plain Laplace (nAGQ = 1) marginal; higher n_quad refines it.
#
# Measured recovery (dev_notes/probe_re_*.R): occupancy arm, per-group n = 8,
# EM attenuates sigma ~18% -> AGHQ ~4% (matching NUTS), correlated sigmas to
# ~1%. Detection arm, 24 seeds, per-group n ~ 10: EM attenuates ~70% (only
# occupied sites inform p) -> AGHQ within ~1% of truth, 88-96% CI coverage.
#
# Correlation regularization. A default LKJ(re.lkj = 1.5) penalty on each
# correlated block (engine `lkj_eta`) pulls a weakly-identified RE correlation
# off the +-1 boundary toward 0 without touching the marginal SDs; eta = 1
# disables it.
# =============================================================================


# AGHQ refinement of an RE occupancy fit. `design` is the per-term RE design
# from .tobs_re_design() (each element tagged with its `arm`: "occ" if the
# random effect enters the occupancy predictor, "det" if it enters detection);
# `beta_occ` / `beta_det` / `Sigma_list` / `b` the converged EM estimates (on
# the same scaled-X / natural-Z coordinates the EM used). Builds the occupancy /
# detection site-marginal callback and delegates the quadrature to
# `tulpa::tulpa_re_aghq()`. Returns a list of refined estimates (point
# estimates, BLUPs, per-group BLUP variances, marginal fixed-effect SEs,
# refreshed occupancy weights, and the `arm`) or NULL when the pass does not
# apply (caller keeps the EM result) -- including when the RE is split across
# both arms or shares no single grouping factor.
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
  dtot <- sum(vapply(design, function(d) as.integer(d$n_coefs), integer(1)))
  if (dtot > 3L) return(NULL)

  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  y <- model$y
  p_occ <- ncol(X_occ); p_det <- ncol(X_det)
  N <- nrow(X_occ)

  # ---- per-site detection summary (fixed across the optimization) ----
  vmask   <- y >= 0
  n_valid <- rowSums(vmask)
  Y <- y; Y[!vmask] <- 0
  n_det   <- rowSums(Y)
  any_det <- n_det > 0
  keep    <- which(n_valid > 0)
  cl <- .tobs_clamp_eta

  # ---- family-specific site marginal (the only occupancy-aware part) ----
  # make_site(theta) closes over the current (beta_occ, beta_det); the engine
  # supplies the RE offset Z b through the eta passed to deriv / lmat.
  if (arm == "occ") {
    make_site <- function(theta) {
      bo <- theta[seq_len(p_occ)]; bd <- theta[p_occ + seq_len(p_det)]
      p  <- plogis(cl(as.numeric(X_det %*% bd)))      # detection, no b
      dterm   <- ifelse(any_det, n_det * log(p) + (n_valid - n_det) * log1p(-p), 0)
      a_nodet <- ifelse(any_det, 0, 1 - exp(n_valid * log1p(-p)))   # 1 - (1-p)^n
      list(
        eta_re = as.numeric(X_occ %*% bo),
        deriv = function(rows, eta) {
          s <- plogis(eta); ad <- any_det[rows]
          logL <- d1 <- d2 <- numeric(length(rows))
          logL[ad] <- log(s[ad]) + dterm[rows][ad]
          d1[ad]   <- 1 - s[ad]
          d2[ad]   <- -s[ad] * (1 - s[ad])
          nd <- !ad
          aa <- a_nodet[rows][nd]; sn <- s[nd]; sp <- sn * (1 - sn)
          g_ <- 1 - aa * sn; Nn <- -aa * sp
          logL[nd] <- log(g_)
          d1[nd] <- Nn / g_
          d2[nd] <- (-aa * sp * (1 - 2 * sn) * g_ - Nn * (-aa * sp)) / g_^2
          list(logL = logL, d1 = d1, d2 = d2)
        },
        lmat = function(rows, ETA) {
          S <- plogis(ETA); ad <- any_det[rows]
          out <- matrix(0, length(rows), ncol(ETA))
          if (any(ad))  out[ad, ]  <- log(S[ad, , drop = FALSE]) + dterm[rows][ad]
          if (any(!ad)) out[!ad, ] <- log(1 - a_nodet[rows][!ad] * S[!ad, , drop = FALSE])
          out
        })
    }
  } else {  # arm == "det": b moves p, psi fixed
    make_site <- function(theta) {
      bo <- theta[seq_len(p_occ)]; bd <- theta[p_occ + seq_len(p_det)]
      psi <- plogis(cl(as.numeric(X_occ %*% bo)))     # occupancy, no b
      list(
        eta_re = as.numeric(X_det %*% bd),
        deriv = function(rows, eta) {
          p <- plogis(eta); ad <- any_det[rows]
          nv <- n_valid[rows]; nd <- n_det[rows]; ps <- psi[rows]
          logL <- d1 <- d2 <- numeric(length(rows))
          logL[ad] <- log(ps[ad]) + nd[ad] * log(p[ad]) + (nv[ad] - nd[ad]) * log1p(-p[ad])
          d1[ad]   <- nd[ad] - nv[ad] * p[ad]
          d2[ad]   <- -nv[ad] * p[ad] * (1 - p[ad])
          o <- !ad
          pn <- p[o]; nvn <- nv[o]; psn <- ps[o]
          q <- exp(nvn * log1p(-pn)); L <- psn * q + (1 - psn)
          A  <- -psn * nvn * pn * q
          dA <- -psn * nvn * pn * q * ((1 - pn) - nvn * pn)
          logL[o] <- log(L)
          d1[o] <- A / L
          d2[o] <- (dA * L - A^2) / L^2
          list(logL = logL, d1 = d1, d2 = d2)
        },
        lmat = function(rows, ETA) {
          P <- plogis(ETA); ad <- any_det[rows]
          nv <- n_valid[rows]; nd <- n_det[rows]; ps <- psi[rows]
          out <- matrix(0, length(rows), ncol(ETA))
          if (any(ad)) {
            out[ad, ] <- log(ps[ad]) + nd[ad] * log(P[ad, , drop = FALSE]) +
              (nv[ad] - nd[ad]) * log1p(-P[ad, , drop = FALSE])
          }
          if (any(!ad)) {
            Q <- exp(nv[!ad] * log1p(-P[!ad, , drop = FALSE]))
            out[!ad, ] <- log(ps[!ad] * Q + (1 - ps[!ad]))
          }
          out
        })
    }
  }

  # Per-term RE spec for the engine (correlated blocks keep their Z).
  re_terms <- lapply(design, function(d) list(
    idx = as.integer(d$idx), n_groups = as.integer(d$n_groups),
    n_coefs = as.integer(d$n_coefs),
    Z = if (d$n_coefs > 1L) d$Z else NULL,
    correlated = isTRUE(d$correlated)))

  ref <- tulpa::tulpa_re_aghq(
    theta0 = c(beta_occ, beta_det), re_terms = re_terms,
    Sigma0 = Sigma_list, make_site = make_site, n_obs = N,
    keep = keep, n_quad = n_quad, lkj_eta = lkj_eta)
  if (is.null(ref)) return(NULL)

  beta_occ_ref <- ref$theta[seq_len(p_occ)]
  beta_det_ref <- ref$theta[p_occ + seq_len(p_det)]

  # Re-pack the per-term BLUP matrices (n_groups x n_coefs) into the term-major,
  # group-major (byrow) layout .tobs_re_param_block() / .tobs_re_offset() use.
  b_out <- unlist(lapply(ref$blup,     function(M) as.numeric(t(M))), use.names = FALSE)
  bvar_out <- unlist(lapply(ref$blup_var, function(M) as.numeric(t(M))), use.names = FALSE)

  # Refreshed P(z = 1 | y) at the refined estimate, for fitted() / residuals().
  if (arm == "occ") {
    psi <- plogis(cl(as.numeric(X_occ %*% beta_occ_ref) + .tobs_re_offset(design, b_out)))
    p   <- plogis(cl(as.numeric(X_det %*% beta_det_ref)))
  } else {
    psi <- plogis(cl(as.numeric(X_occ %*% beta_occ_ref)))
    p   <- plogis(cl(as.numeric(X_det %*% beta_det_ref) + .tobs_re_offset(design, b_out)))
  }
  prod0 <- ifelse(any_det, 0, exp(n_valid * log1p(-p)))   # (1-p)^n
  w_ref <- ifelse(any_det, 1, psi * prod0 / (psi * prod0 + (1 - psi)))

  # Per-group solve status: a group the engine could not solve has NA BLUP and
  # BLUP-variance rows, which would otherwise replace the EM's finite values for
  # that group. The caller declines the whole refinement on it
  # (gcol33/tulpaObs#281).
  gstat <- .tobs_aghq_group_status(ref)

  list(
    ok           = TRUE,
    arm          = arm,
    beta_occ     = beta_occ_ref,
    beta_det     = beta_det_ref,
    Sigma_list   = ref$Sigma_list,
    b            = b_out,
    b_var        = bvar_out,
    beta_occ_se  = ref$theta_se[seq_len(p_occ)],
    det_se       = ref$theta_se[p_occ + seq_len(p_det)],
    weights      = w_ref,
    n_quad       = ref$n_quad,
    lkj_eta      = ref$lkj_eta,
    converged    = .tobs_aghq_converged(ref, gstat),
    group_ok     = gstat$group_ok,
    groups_failed = gstat$failed,
    n_iter       = .tobs_aghq_n_iter(ref)
  )
}


# =============================================================================
# Per-group solve status on an AGHQ fit (gcol33/tulpaObs#281).
#
# tulpa::tulpa_re_aghq() reports, per group, whether that group's posterior mode
# search and precision factorization succeeded (`group_ok`, gcol33/tulpa#605).
# A group it marks FALSE comes back with NA `blup` / `blup_var` / `blup_cov_g`
# rows, so every per-group quantity read off it is NA and every community-level
# quantity summing over groups is not a function of that group's data. The
# engine also raises a warning naming the failures, but a warning is the one
# channel a caller cannot read without parsing indices out of message text, so
# every AGHQ consumer here reads the status through these two functions and
# carries it onto the fit it returns.
# =============================================================================

# Normalize the engine's per-group status into the record a fit carries.
# `group_names` labels the groups where the caller has names for them (species,
# sites, observers); it is dropped when it does not match the group count rather
# than guessing an alignment.
.tobs_aghq_group_status <- function(ref, group_names = NULL) {
  ok <- as.logical(ref$group_ok)
  failed <- which(!ok)
  named <- !is.null(group_names) && length(group_names) == length(ok)
  list(group_ok     = ok,
       n_groups     = length(ok),
       failed       = failed,
       failed_names = if (named) as.character(group_names)[failed] else NULL,
       all_ok       = length(failed) == 0L)
}

# The verdict a fit reports. A fit carrying a failed group is not an estimate of
# what that group contributes to, so it is not reported as converged even where
# the optimizer itself stopped cleanly.
.tobs_aghq_converged <- function(ref, status = .tobs_aghq_group_status(ref)) {
  isTRUE(ref$converged) && isTRUE(status$all_ok)
}

# The effort behind an AGHQ fit, in the units the package's other optim-driven
# fitters already report as `n_iter`: stats::optim's function-evaluation count
# (the joint AGHQ driver is one BFGS call, which counts evaluations rather than
# iterations).
.tobs_aghq_n_iter <- function(ref) {
  cnt <- ref$counts
  if (is.null(cnt) || !length(cnt)) return(NA_integer_)
  # optim names the entries; an unnamed vector still has the function count
  # first. (`cnt[["function"]]` on a vector without that name is an error, not
  # NULL, so the name is checked rather than defaulted.)
  as.integer(if ("function" %in% names(cnt)) cnt[["function"]] else cnt[[1L]])
}

# The `fit$convergence` record for a family whose estimates came from the AGHQ
# engine. One constructor, so a caller filters the same way whatever it fitted:
# `converged`, `n_iter`, and the per-group solve status. `raw` carrying no
# status (every non-AGHQ path) yields the two fields it always had.
# `group_names` labels the failures where the family has names for its groups.
.tobs_aghq_convergence_record <- function(raw, group_names = NULL,
                                          converged = isTRUE(raw$converged)) {
  failed <- raw$groups_failed
  named  <- length(failed) && !is.null(group_names) &&
    length(group_names) >= max(failed)
  list(converged = converged,
       n_iter    = raw$n_iter %||% NA_integer_,
       group_ok  = raw$group_ok,
       groups_failed = failed,
       groups_failed_names = if (named) as.character(group_names)[failed] else NULL)
}
