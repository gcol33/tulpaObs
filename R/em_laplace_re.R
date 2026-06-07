# =============================================================================
# Deterministic variance-component EM for occupancy random effects
# (gcol33/tulpaObs#11).
#
# tulpa's Laplace engine (tulpa_laplace) finds the latent mode at a FIXED RE
# covariance. This driver wraps it in the two loops the engine does not carry
# itself:
#
#   * the occupancy missing-data EM (latent occupancy z), feeding the RE
#     posterior mode back into psi at every E-step (so the group effects enter
#     the occupancy weights, not just the M-step fit);
#   * a variance-component update for the per-term RE covariance Sigma, from the
#     RE posterior mode b-hat and its per-group posterior covariance (the
#     standard EM/REML update Sigma_k <- mean_g [ b_g b_g' + Cov(b_g | y) ]).
#
# The M-step reuses the package's pseudo-binomial encoding (y = round(w*M),
# n = M) because tulpa_laplace rejects a fractional binomial response and the
# direct two-row weighted encoding does not converge in the RE-aware Newton
# solve. The M-inflation scales the data term by M; rescaling the RE prior to
# Sigma/M scales the penalty by the same M, so the penalized MAP of (beta, b)
# is unchanged (the argmax is scale-invariant). The joint Hessian then scales
# by M exactly -- data information via n = M, prior precision via (Sigma/M)^{-1}
# = M Sigma^{-1} -- so the per-group posterior covariance at NATURAL scale is
# M times the block tulpa returns from the inflated fit (return_re_cov = TRUE).
#
# Scope: iid intercept RE, uncorrelated random slopes (diagonal Sigma), and
# CORRELATED random slopes (a full Sigma, lme4 `(1 + x | g)`) on the occupancy
# predictor of a single-season occupancy model. The per-term Sigma is full for
# a correlated block and projected to its diagonal each M-step otherwise, so the
# uncorrelated path is the special case Sigma = diag(sigma^2).
#
# Bias note. The random-effect block b is integrated by Laplace (mode + Gaussian
# curvature of the TRUE binomial conditional). The occupancy state z is
# integrated exactly by the outer EM, and the M-step linearizes nothing -- it
# fits the real likelihood -- so this is NOT Breslow-Clayton PQL (no
# working-response linearization). It is the lme4 glmer nAGQ=1 regime: a
# Laplace-approximate marginal likelihood. For binary data that Laplace integral
# attenuates the VARIANCE COMPONENTS (sigma, and the correlation of a full block)
# toward zero at small per-group sample size. The fixed-effect estimate is the
# conditional mode and its SE is read at natural scale (.tobs_re_occ_fixed_se),
# not off the M-inflated M-step Hessian, so the attenuation is confined to the RE
# covariance.
#
# By DEFAULT (aghq = TRUE) that attenuation is then removed: after the EM
# converges, .tobs_re_aghq() (R/re_aghq.R) refines the variance components on the
# exact-marginal adaptive Gauss-Hermite likelihood (the nAGQ > 1 fix), cutting
# the per-group-n = 8 sigma bias from ~18% to ~4%. A weakly-identified RE
# correlation is regularized off the +-1 boundary by a default LKJ(eta = 1.5)
# penalty (control re.lkj; see R/re_aghq.R). Set re.aghq = FALSE for the raw
# nAGQ = 1 EM. NUTS integrates b by MCMC and is available for a full posterior
# treatment of the RE correlation. See NEWS and ?tobs.
# =============================================================================


# Build the per-term RE design the deterministic driver and the variance-
# component update share. Each element carries the 1-based group index, the
# n_obs x n_coefs design Z (intercept column first when the block keeps the
# group intercept, then the slope columns), and the coefficient labels.
.tobs_re_design <- function(re_list, model) {
  N <- model$N %||% model$n_sites
  lapply(seq_along(re_list), function(t) {
    re <- re_list[[t]]
    idx <- as.integer(re$group_idx)
    has_int <- isTRUE(re$intercept)
    if (identical(re$type, "slope") && !is.null(re$covariate)) {
      Xs <- .tobs_re_slope_matrix(re$covariate, model$data)
      Z <- if (has_int) cbind(1, Xs) else Xs
      coef_names <- if (has_int) c("(Intercept)", colnames(Xs)) else colnames(Xs)
    } else {
      Z <- matrix(1, N, 1L)
      coef_names <- "(Intercept)"
      has_int <- TRUE
    }
    colnames(Z) <- coef_names
    # Group label for naming sigma / BLUP rows: the term's group is an integer
    # code (the original factor levels are not retained), so label terms g1,
    # g2, ... in formula order.
    list(idx = idx, n_groups = as.integer(re$n_groups %||% max(idx)),
         n_coefs = ncol(Z), Z = Z, coef_names = coef_names,
         has_intercept = has_int,
         correlated = isTRUE(re$correlated) && ncol(Z) > 1L,
         group_label = sprintf("g%d", t))
  })
}


# Per-observation random-effect offset sum_t sum_c Z_t[i, c] * b_t[group, c].
# `b` is the concatenated latent vector (term-major, then group-major within a
# term: [t1_g1_c1, t1_g1_c2, ..., t1_g2_c1, ...]).
.tobs_re_offset <- function(design, b) {
  if (!length(design)) return(0)  # no RE on this arm -> scalar 0 broadcasts
  N <- nrow(design[[1]]$Z)
  off <- numeric(N)
  pos <- 0L
  for (d in design) {
    nc <- d$n_coefs
    bt <- b[pos + seq_len(d$n_groups * nc)]
    # b_block[g, c] with group-major layout -> matrix [n_groups x nc].
    Bm <- matrix(bt, nrow = d$n_groups, ncol = nc, byrow = TRUE)
    off <- off + rowSums(d$Z * Bm[d$idx, , drop = FALSE])
    pos <- pos + d$n_groups * nc
  }
  off
}


# RE precision Sigma^{-1} for one term, with a tiny ridge so a degenerate
# (e.g. collapsed-variance) block stays invertible. For a diagonal Sigma this
# is diag(1 / (sigma_c^2 + 1e-10)) -- the uncorrelated path -- so both the
# correlated and uncorrelated blocks come from one expression.
.tobs_re_precision <- function(Sigma) {
  nc <- nrow(Sigma)
  chol2inv(chol(Sigma + diag(1e-10, nc)))
}


# Sparse RE design Z (n_obs x sum(n_groups*n_coefs)) and the block-diagonal RE
# prior precision Q (Sigma^{-1} per group, full for a correlated term, diagonal
# otherwise) in the term-/group-major latent order. Single source of truth for
# the marginal-SE Schur complement below.
.tobs_re_ZD <- function(X, design, Sigma_list) {
  Z_parts <- list()
  Q_blocks <- list()
  for (k in seq_along(design)) {
    d <- design[[k]]; nc <- d$n_coefs; ng <- d$n_groups
    ii <- rep(seq_len(nrow(X)), each = nc)
    jj <- rep((d$idx - 1L) * nc, each = nc) + rep(seq_len(nc), nrow(X))
    Z_parts[[k]] <- Matrix::sparseMatrix(
      i = ii, j = jj, x = as.numeric(t(d$Z)), dims = c(nrow(X), ng * nc))
    Qk <- .tobs_re_precision(Sigma_list[[k]])
    Q_blocks[[k]] <- Matrix::bdiag(
      rep(list(Matrix::Matrix(Qk, sparse = TRUE)), ng))
  }
  list(Z = do.call(cbind, Z_parts), Q = Matrix::bdiag(Q_blocks))
}


# Per-group posterior covariance blocks Cov(b_{k,g} | y, Sigma) at NATURAL
# scale, from tulpa's M-inflated fit (return_re_cov = TRUE). tulpa returns the
# diagonal blocks of the inflated inverse Hessian in term-major then group-major
# order; the inflated Hessian is M times the natural one (n = M data weights,
# Sigma/M prior), so multiplying each block by M recovers the natural-scale
# covariance. Returns a list aligned with `design`: per term a list of n_groups
# nc x nc matrices.
.tobs_re_cov_natural <- function(cov_blocks, design, M) {
  out <- vector("list", length(design))
  idx <- 0L
  for (k in seq_along(design)) {
    ng <- design[[k]]$n_groups
    blk <- vector("list", ng)
    for (g in seq_len(ng)) {
      idx <- idx + 1L
      blk[[g]] <- M * as.matrix(cov_blocks[[idx]])
    }
    out[[k]] <- blk
  }
  out
}


# EM/REML variance-component update: Sigma_k <- mean_g [ b_g b_g' + Cov(b_g) ].
# `cov_nat` is the natural-scale per-group posterior covariance from
# .tobs_re_cov_natural(). A correlated term keeps the full (off-diagonal)
# Sigma; an uncorrelated term is projected to its diagonal. Variances are
# floored so a collapsed block stays positive-definite.
.tobs_re_sigma_update <- function(design, b, cov_nat) {
  Sigma_new <- vector("list", length(design))
  pos <- 0L
  for (k in seq_along(design)) {
    d <- design[[k]]; nc <- d$n_coefs; ng <- d$n_groups
    Bm <- matrix(b[pos + seq_len(ng * nc)], nrow = ng, ncol = nc, byrow = TRUE)
    S <- matrix(0, nc, nc)
    for (g in seq_len(ng)) S <- S + tcrossprod(Bm[g, ]) + cov_nat[[k]][[g]]
    S <- S / ng
    if (!isTRUE(d$correlated)) S <- diag(diag(S), nc)  # project to diagonal
    diag(S) <- pmax(diag(S), 1e-6)                      # floor variances
    Sigma_new[[k]] <- S
    pos <- pos + ng * nc
  }
  Sigma_new
}


# Diagonal of the per-group posterior covariance, in the term-/group-major
# latent order of `b` (the BLUP standard errors). Derived from the same
# natural-scale blocks as the variance update -- single source of truth.
.tobs_re_bvar_from_cov <- function(design, cov_nat) {
  out <- numeric(0)
  for (k in seq_along(design)) {
    d <- design[[k]]; nc <- d$n_coefs; ng <- d$n_groups
    v <- numeric(ng * nc)
    for (g in seq_len(ng)) v[(g - 1L) * nc + seq_len(nc)] <- diag(cov_nat[[k]][[g]])
    out <- c(out, v)
  }
  out
}


# Standard errors of the occupancy fixed effects on the RE path. The M-step
# fits a pseudo-binomial (n = M) design, so its returned H_beta is M-inflated
# and unusable for SEs. Compute the observed-data information at NATURAL scale
# instead: the Louis observed weights D = psi(1-psi) - w(1-w) (complete-data
# info minus the missing-occupancy score variance, mirroring
# .louis_info_psi_single) marginalised over the random-effect block via a Schur
# complement with the full RE precision Q. `w` is the converged E-step weight.
.tobs_re_occ_fixed_se <- function(X, eta, w, design, Sigma_list) {
  psi <- plogis(eta)
  D <- as.numeric(psi * (1 - psi) - w * (1 - w))
  ZD <- .tobs_re_ZD(X, design, Sigma_list)
  Z <- ZD$Z
  XtDX <- crossprod(X, D * X)
  ZtDZ <- as.matrix(Matrix::crossprod(Z, D * Z))
  XtDZ <- as.matrix(crossprod(X, D * Z))
  C <- ZtDZ + as.matrix(ZD$Q)
  H_beta <- tryCatch(XtDX - XtDZ %*% solve(C, t(XtDZ)),
                     error = function(e) XtDX)
  .se_from_info(as.matrix(H_beta), ncol(X))
}


# Build the per-term `re_list` tulpa_laplace() consumes from a design and the
# current per-term covariance. `rows` optionally restricts the design to a row
# subset (the detection arm drops near-empty sites before fitting); the group
# count is left intact so tulpa still returns one latent block per group.
# `inflate` is the M-step pseudo-binomial inflation factor: the occupancy arm
# encodes the soft weight as a binomial with n = M trials and so passes the RE
# prior at covariance Sigma / M (the penalty scales with the data term, leaving
# the penalised MAP unchanged); the detection arm is a genuine weighted binomial
# (inflate = 1) and passes Sigma at natural scale. A correlated term passes the
# full covariance, an uncorrelated term its per-coefficient marginal SD.
.re_list_for_tulpa <- function(design, Sigma_list, rows = NULL, inflate = 1) {
  lapply(seq_along(design), function(k) {
    d   <- design[[k]]
    idx <- if (is.null(rows)) d$idx else d$idx[rows]
    Z   <- if (d$n_coefs > 1L) {
      if (is.null(rows)) d$Z else d$Z[rows, , drop = FALSE]
    } else NULL
    el <- list(idx = idx, n_groups = d$n_groups, n_coefs = d$n_coefs, Z = Z)
    if (isTRUE(d$correlated)) {
      el$cov <- Sigma_list[[k]] / inflate
    } else {
      el$sigma <- sqrt(diag(Sigma_list[[k]])) / sqrt(inflate)
    }
    el
  })
}


# Partition a list of `tobs_re` specs into two predictor arms by their
# `$shared = c(arm1_on, arm2_on)` membership, tagging each design element with
# its arm name and giving the second arm's terms a `p<t>` group label (matching
# the `p_` detection fixed-effect prefix) so the two arms never collide in the
# parameter block. A term flagged on BOTH arms is rejected with `both_msg` --
# the deterministic paths fit a separate RE block per arm, not one realization
# shared across them.
.tobs_re_split_two_arms <- function(re_list, model, arm1, arm2, both_msg) {
  arm_of <- function(r) {
    sh <- r$shared
    on1 <- length(sh) >= 1L && isTRUE(sh[1])
    on2 <- length(sh) >= 2L && isTRUE(sh[2])
    if (on1 && on2) stop(both_msg, call. = FALSE)
    if (on2) arm2 else arm1
  }
  arms <- vapply(re_list, arm_of, character(1))
  tag <- function(sub, arm) {
    if (!length(sub)) return(list())
    design <- .tobs_re_design(sub, model)
    lapply(seq_along(design), function(i) {
      d <- design[[i]]
      d$arm <- arm
      if (arm == arm2) d$group_label <- sprintf("p%d", i)
      d
    })
  }
  stats::setNames(list(tag(re_list[arms == arm1], arm1),
                       tag(re_list[arms == arm2], arm2)),
                  c(arm1, arm2))
}

# Occupancy/detection arm split (Laplace path). A shared term routes to NUTS.
.tobs_re_split_arms <- function(re_list, model) {
  .tobs_re_split_two_arms(
    re_list, model, "occ", "det",
    paste0("A random effect shared across occupancy and detection is not ",
           "supported on the Laplace path. Use method = 'nuts'."))
}


#' Fit single-season occupancy with formula random effects via Laplace + a
#' variance-component EM (internal).
#'
#' @param model A single-season `tobs_model`.
#' @param re A list of `tobs_re` specs. Each enters either the occupancy or the
#'   detection predictor (its `$shared = c(occ, det)` membership); the iid
#'   intercept, uncorrelated-slope, and correlated-slope forms are supported on
#'   both arms. A term shared across both predictors is rejected (use NUTS).
#' @param priors Prior spec; not applied on this path (warns when active).
#' @param max_iter,tol,damping EM controls.
#' @param aghq Logical; run the adaptive Gauss-Hermite debias pass on the
#'   variance components after the EM converges (default `TRUE`). See
#'   `R/re_aghq.R`.
#' @param n_quad Quadrature points per random-effect dimension for the AGHQ
#'   debias (default 9).
#' @param lkj_eta LKJ shape for the RE-correlation regularization in the AGHQ
#'   debias (default 1.5; `1` disables it). See `R/re_aghq.R`.
#' @param verbose Print per-iteration progress.
#' @keywords internal
.tobs_em_laplace_re <- function(model, re, priors = NULL,
                                max_iter = 100L, tol = 1e-5, damping = 0.3,
                                aghq = TRUE, n_quad = 9L, lkj_eta = 1.5,
                                verbose = TRUE) {
  if (inherits(re, "tobs_re")) re <- list(re)
  arms <- .tobs_re_split_arms(re, model)
  design_occ <- arms$occ        # RE on the occupancy predictor (may be empty)
  design_det <- arms$det        # RE on the detection predictor (may be empty)

  y <- model$y
  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  N <- nrow(X_occ)
  p_occ <- ncol(X_occ)
  p_det <- ncol(X_det)
  M <- 1000L

  n_valid <- integer(N); n_det <- integer(N); any_det <- logical(N)
  for (i in seq_len(N)) {
    v <- y[i, ] >= 0
    n_valid[i] <- sum(v)
    n_det[i]   <- sum(y[i, v] == 1)
    any_det[i] <- n_det[i] > 0
  }
  keep <- n_valid > 0

  # State. Each arm carries its own latent block b and per-term covariance;
  # a starting diagonal sigma = 0.5 mirrors the historical single-arm path.
  init <- glm_init(X_occ, X_det, any_det, n_det, n_valid, keep, p_occ, p_det)
  beta_occ <- init$occ$beta
  beta_det <- init$det$beta
  n_lat <- function(d) sum(vapply(d, function(x) x$n_groups * x$n_coefs, integer(1)))
  b_occ <- numeric(n_lat(design_occ))
  b_det <- numeric(n_lat(design_det))
  Sigma_occ <- lapply(design_occ, function(d) diag(0.25, d$n_coefs))
  Sigma_det <- lapply(design_det, function(d) diag(0.25, d$n_coefs))

  weights <- NULL
  converged <- FALSE
  occ_fit <- NULL; det_fit <- NULL
  # Progress + ETA for the variance-component RE-EM iterations
  # (gcol33/tulpaObs#43); ON by default, finalised on convergence.
  .prog <- tulpa:::.tulpa_iter_progress("re-em", max_iter, unit = "iter")
  for (it in seq_len(max_iter)) {
    # ---- E-step: psi and p both carry their arm's RE posterior mode. ----
    eta_occ <- as.numeric(X_occ %*% beta_occ) + .tobs_re_offset(design_occ, b_occ)
    eta_det <- as.numeric(X_det %*% beta_det) + .tobs_re_offset(design_det, b_det)
    psi <- plogis(eta_occ)
    p_site <- plogis(eta_det)
    w <- occ_weights(psi, p_site, N, n_valid, n_det, any_det)
    weights <- w

    # ---- M-step (occupancy): pseudo-binomial, RE prior rescaled by M. ----
    y_occ <- pmin(pmax(ifelse(any_det, M, as.integer(round(w * M))), 0L), M)
    fo <- tulpa::tulpa_laplace(
      y = y_occ, n_trials = rep(M, N), X = X_occ,
      re_list = .re_list_for_tulpa(design_occ, Sigma_occ, inflate = M),
      family = "binomial", return_hessian = TRUE,
      return_re_cov = length(design_occ) > 0L)
    beta_occ_new <- fo$mode[seq_len(p_occ)]
    b_occ_new <- if (length(design_occ)) fo$mode[-seq_len(p_occ)] else numeric(0)
    Sigma_occ_new <- if (length(design_occ)) {
      .tobs_re_sigma_update(design_occ, b_occ_new,
                            .tobs_re_cov_natural(fo$cov_blocks, design_occ, M))
    } else list()

    # ---- M-step (detection): weighted binomial, RE prior at natural scale. ----
    # Near-empty sites (w ~ 0) drop out; the group count stays intact so tulpa
    # still returns one latent block per detection group. The fit is a genuine
    # binomial (no M-inflation), so its posterior cov is natural-scale already.
    w_det <- w; w_det[any_det] <- 1
    keep_det <- keep & (w_det > 1e-6)
    fd <- tulpa::tulpa_laplace(
      y = n_det[keep_det], n_trials = n_valid[keep_det],
      X = X_det[keep_det, , drop = FALSE], weights = w_det[keep_det],
      re_list = .re_list_for_tulpa(design_det, Sigma_det, rows = keep_det,
                                   inflate = 1),
      family = "binomial", return_hessian = TRUE,
      return_re_cov = length(design_det) > 0L)
    beta_det_new <- fd$mode[seq_len(p_det)]
    b_det_new <- if (length(design_det)) fd$mode[-seq_len(p_det)] else numeric(0)
    Sigma_det_new <- if (length(design_det)) {
      .tobs_re_sigma_update(design_det, b_det_new,
                            .tobs_re_cov_natural(fd$cov_blocks, design_det, 1))
    } else list()

    delta <- max(abs(c(beta_occ_new - beta_occ,
                       beta_det_new - beta_det,
                       unlist(Sigma_occ_new) - unlist(Sigma_occ),
                       unlist(Sigma_det_new) - unlist(Sigma_det))))
    beta_occ <- beta_occ_new; b_occ <- b_occ_new; Sigma_occ <- Sigma_occ_new
    beta_det <- beta_det_new; b_det <- b_det_new; Sigma_det <- Sigma_det_new
    occ_fit <- fo; det_fit <- fd

    .prog$tick()
    if (verbose) cat(sprintf("  RE-EM iter %d: delta = %.6g\n", it, delta))
    if (is.finite(delta) && delta < tol) { converged <- TRUE; break }
  }
  .prog$finish()

  # ---- Combine the two arms into one design / latent layout (occ then det). ----
  # .tobs_re_param_block(), ranef(), and the AGHQ pass all consume one ordered
  # set of (design, b, b_var, Sigma) parallel lists; occupancy terms come first.
  design <- c(design_occ, design_det)

  # Final per-group posterior covariance (BLUP SEs + reporting): each arm reads
  # its own M-step fit (occ inflated by M, det natural scale).
  bvar_occ <- if (length(design_occ))
    .tobs_re_bvar_from_cov(design_occ,
                           .tobs_re_cov_natural(occ_fit$cov_blocks, design_occ, M))
  else numeric(0)
  bvar_det <- if (length(design_det))
    .tobs_re_bvar_from_cov(design_det,
                           .tobs_re_cov_natural(det_fit$cov_blocks, design_det, 1))
  else numeric(0)
  b <- c(b_occ, b_det); b_var <- c(bvar_occ, bvar_det)
  Sigma_list <- c(Sigma_occ, Sigma_det)

  # Occupancy fixed-effect SE at natural scale (the M-step H_beta is inflated);
  # the detection SE is tulpa's RE-marginalised H_beta on the det arm (read by
  # .se_from_laplace_fit from det_fit$se). AGHQ recalibrates both when it runs.
  eta_mode <- as.numeric(X_occ %*% beta_occ) + .tobs_re_offset(design_occ, b_occ)
  beta_occ_se <- if (length(design_occ))
    .tobs_re_occ_fixed_se(X_occ, eta_mode, weights, design_occ, Sigma_occ)
  else .se_from_laplace_fit(occ_fit, p_occ)
  aghq_status <- list(applied = FALSE)

  # ---- AGHQ debias of the variance components (R/re_aghq.R). ----
  # The EM integrates b by Laplace, which attenuates sigma/correlation for
  # binary data at small per-group n. Refine on the exact-marginal (adaptive
  # Gauss-Hermite) likelihood: removes the attenuation, recalibrates the
  # fixed-effect SEs off the marginal Hessian, and refreshes the BLUPs. Applies
  # to a single grouping factor on one arm (RE dim <= 3); falls back to the EM
  # result on any failure, on RE split across both arms, or on crossed / nested
  # groupings. The arm is read from the combined design.
  if (isTRUE(aghq)) {
    ref <- tryCatch(
      .tobs_re_aghq(model, design, beta_occ, beta_det, Sigma_list, b,
                    n_quad = n_quad, lkj_eta = lkj_eta),
      error = function(e) {
        warning("AGHQ variance-component refine failed; keeping the EM result. ",
                "Cause: ", conditionMessage(e), call. = FALSE)
        NULL
      })
    if (!is.null(ref) && isTRUE(ref$ok)) {
      beta_occ <- ref$beta_occ; beta_det <- ref$beta_det
      Sigma_list <- ref$Sigma_list; b <- ref$b; b_var <- ref$b_var
      weights <- ref$weights
      occ_fit$beta <- beta_occ            # build_laplace_fit reads $beta first
      det_fit$beta <- beta_det
      # Both fixed-effect SEs come from the one joint marginal Hessian (psi and p
      # are coupled through the occupancy weight), independent of which arm
      # carries the RE.
      beta_occ_se <- ref$beta_occ_se
      det_fit$se  <- ref$det_se           # .se_from_laplace_fit reads $se first
      aghq_status <- list(applied = TRUE, arm = ref$arm, n_quad = ref$n_quad,
                          lkj_eta = ref$lkj_eta, converged = ref$converged)
    }
  }

  list(
    fits = list(occ = occ_fit, det = det_fit),
    weights = weights,
    convergence = list(converged = converged, n_iter = it),
    correction = "none",
    aghq = aghq_status,
    re_post = list(design = design, b = b, b_var = b_var, Sigma = Sigma_list,
                   beta_occ_se = beta_occ_se)
  )
}


# Assemble the named RE parameter block (sigma + correlation hyperparameters and
# per-group BLUPs) and the structured `re_effects` summary the deterministic RE
# fit appends to a tobs_fit. Single source of truth for the deterministic-path
# RE layout and labels.
#
# Returns a list with `means`, `sds`, `names` (to append to the fixed-effect
# parameter vector) and `re_effects` (per-term BLUP tables consumed by
# ranef.tobs_fit()). For a correlated term the off-diagonal of the estimated RE
# covariance is reported as `cor_<g>_<ci>_<cj>` correlations.
.tobs_re_param_block <- function(re_post) {
  design <- re_post$design
  b      <- re_post$b
  b_var  <- re_post$b_var
  Sigma  <- re_post$Sigma

  means <- numeric(0); sds <- numeric(0); nms <- character(0)
  re_effects <- list()
  pos <- 0L
  for (k in seq_along(design)) {
    d <- design[[k]]; nc <- d$n_coefs; ng <- d$n_groups
    g <- d$group_label
    Sk <- Sigma[[k]]
    sig_k <- sqrt(pmax(diag(Sk), 0))
    # sigma hyperparameters (marginal SD, one per coefficient).
    means <- c(means, sig_k)
    sds   <- c(sds, rep(NA_real_, nc))
    nms   <- c(nms, sprintf("sigma_%s_%s", g, d$coef_names))
    # Correlation hyperparameters for a correlated block (upper off-diagonal of
    # the estimated RE covariance, on the correlation scale).
    if (isTRUE(d$correlated) && nc > 1L) {
      for (ci in seq_len(nc - 1L)) for (cj in seq(ci + 1L, nc)) {
        denom <- sig_k[ci] * sig_k[cj]
        means <- c(means, if (denom > 0) Sk[ci, cj] / denom else 0)
        sds   <- c(sds, NA_real_)
        nms   <- c(nms, sprintf("cor_%s_%s_%s", g,
                                d$coef_names[ci], d$coef_names[cj]))
      }
    }
    pos_b <- pos
    # group BLUPs, group-major.
    blk_b <- b[pos_b + seq_len(ng * nc)]
    blk_v <- b_var[pos_b + seq_len(ng * nc)]
    Bm <- matrix(blk_b, nrow = ng, ncol = nc, byrow = TRUE)
    Vm <- matrix(blk_v, nrow = ng, ncol = nc, byrow = TRUE)
    for (c in seq_len(nc)) {
      cn <- d$coef_names[c]
      means <- c(means, Bm[, c])
      sds   <- c(sds, sqrt(pmax(Vm[, c], 0)))
      nms   <- c(nms, sprintf("re_%s_%s[%d]", g, cn, seq_len(ng)))
    }
    re_effects[[g]] <- data.frame(
      group = g,
      level = rep(seq_len(ng), times = nc),
      term  = rep(d$coef_names, each = ng),
      estimate = as.numeric(Bm),
      std.error = as.numeric(sqrt(pmax(Vm, 0))),
      stringsAsFactors = FALSE)
    pos <- pos + ng * nc
  }
  list(means = means, sds = sds, names = nms, re_effects = re_effects,
       occ_se = re_post$beta_occ_se)
}


# NA-aware Gaussian pseudo-draws for an appended random-effect parameter block:
# each column draws rnorm(mean, sd) when its SD is finite (the per-group BLUPs,
# which carry an AGHQ marginal posterior SD), else stays NA (the variance-
# component / correlation columns have no surfaced joint posterior SD, so a
# fabricated near-degenerate column that would read as "known almost exactly" is
# left NA -- the same NA-on-unavailable rule as the NUTS-only sampler
# diagnostics). Shared by the count-family fit packers (build_nmix_fit /
# build_distance_fit) so the appended-RE draw layout has one definition.
.tobs_re_pseudo_draws <- function(re_means, re_sds, re_names, n_pseudo) {
  n_re <- length(re_means)
  re_draws <- matrix(NA_real_, n_pseudo, n_re)
  for (j in seq_len(n_re)) {
    if (is.finite(re_sds[j]))
      re_draws[, j] <- stats::rnorm(n_pseudo, re_means[j], re_sds[j])
  }
  colnames(re_draws) <- re_names
  re_draws
}
