# =============================================================================
# Deterministic variance-component EM for occupancy random effects
# (gcol33/tulpaObs#11).
#
# tulpa's Laplace engine (tulpa_laplace) finds the latent mode at a FIXED RE
# standard deviation and uses a diagonal RE precision. This driver wraps it in
# the two loops the engine does not carry itself:
#
#   * the occupancy missing-data EM (latent occupancy z), feeding the RE
#     posterior mode back into psi at every E-step (so the group effects enter
#     the occupancy weights, not just the M-step fit);
#   * a variance-component update for sigma, per term and per coefficient, from
#     the RE posterior mode b-hat and its Schur-complement posterior variance
#     (the standard EM/REML update Sigma <- mean_g [ b_g b_g' + Var(b_g | y) ]).
#
# The M-step reuses the package's pseudo-binomial encoding (y = round(w*M),
# n = M) because tulpa_laplace rejects a fractional binomial response and the
# direct two-row weighted encoding does not converge in the RE-aware Newton
# solve. The M-inflation scales the data term by M; rescaling the RE prior to
# sigma/sqrt(M) scales the penalty by the same M, so the penalized MAP of
# (beta, b) is unchanged (the argmax is scale-invariant). The variance
# components are then estimated from curvature recomputed at the natural
# (n = 1) scale.
#
# Scope: iid intercept RE and UNCORRELATED random slopes (diagonal Sigma) on
# the occupancy/state predictor of a single-season occupancy model. Correlated
# slopes (a Cholesky-factored covariance) exist only in the NUTS sampler;
# .tobs_laplace() rejects them before reaching this driver. Deterministic
# Laplace variance estimates for binary occupancy carry the usual small-cluster
# (PQL) bias; see NEWS and ?tobs.
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


# Sparse RE design Z (n_obs x sum(n_groups*n_coefs)) and the diagonal prior
# precision D^{-1} for the current per-term, per-coef sigma. Single source of
# truth for the Schur computations below.
.tobs_re_ZD <- function(X, design, sig_list) {
  Z_parts <- list()
  Dinv <- numeric(0)
  for (k in seq_along(design)) {
    d <- design[[k]]; nc <- d$n_coefs; ng <- d$n_groups
    ii <- rep(seq_len(nrow(X)), each = nc)
    jj <- rep((d$idx - 1L) * nc, each = nc) + rep(seq_len(nc), nrow(X))
    Z_parts[[k]] <- Matrix::sparseMatrix(
      i = ii, j = jj, x = as.numeric(t(d$Z)), dims = c(nrow(X), ng * nc))
    di <- numeric(ng * nc)
    for (c in seq_len(nc)) {
      di[seq(c, ng * nc, by = nc)] <- 1 / (sig_list[[k]][c]^2 + 1e-10)
    }
    Dinv <- c(Dinv, di)
  }
  list(Z = do.call(cbind, Z_parts), Dinv = Dinv)
}

# Diagonal of the RE posterior covariance (the b-block of the joint Hessian
# inverse), in the term-/group-major order of the latent vector. Uses the
# complete-data Bernoulli weights W = psi(1-psi): this is the PQL-style EM
# variance-component update, which carries the usual small-cluster downward
# bias for binary responses (documented; NUTS is the calibrated route).
.tobs_re_posterior_var <- function(X, eta, design, sig_list) {
  W <- as.numeric(plogis(eta) * (1 - plogis(eta)))
  ZD <- .tobs_re_ZD(X, design, sig_list)
  Z <- ZD$Z
  XtWX <- crossprod(X, W * X)
  ZtWZ <- as.matrix(Matrix::crossprod(Z, W * Z))
  XtWZ <- as.matrix(crossprod(X, W * Z))
  Schur <- ZtWZ + diag(ZD$Dinv, length(ZD$Dinv)) - t(XtWZ) %*% solve(XtWX, XtWZ)
  cov_bb <- tryCatch(solve(Schur), error = function(e) {
    diag(1 / (diag(Schur) + 1e-8), length(ZD$Dinv))
  })
  diag(cov_bb)
}

# Standard errors of the occupancy fixed effects on the RE path. The M-step
# fits a pseudo-binomial (n = M) design, so its returned H_beta is M-inflated
# and unusable for SEs. Compute the observed-data information at NATURAL scale
# instead: the Louis observed weights D = psi(1-psi) - w(1-w) (complete-data
# info minus the missing-occupancy score variance, mirroring
# .louis_info_psi_single) marginalised over the random-effect block via a Schur
# complement. `w` is the converged E-step occupancy weight.
.tobs_re_occ_fixed_se <- function(X, eta, w, design, sig_list) {
  psi <- plogis(eta)
  D <- as.numeric(psi * (1 - psi) - w * (1 - w))
  ZD <- .tobs_re_ZD(X, design, sig_list)
  Z <- ZD$Z
  XtDX <- crossprod(X, D * X)
  ZtDZ <- as.matrix(Matrix::crossprod(Z, D * Z))
  XtDZ <- as.matrix(crossprod(X, D * Z))
  C <- ZtDZ + diag(ZD$Dinv, length(ZD$Dinv))
  H_beta <- tryCatch(XtDX - XtDZ %*% solve(C, t(XtDZ)),
                     error = function(e) XtDX)
  .se_from_info(as.matrix(H_beta), ncol(X))
}


#' Fit single-season occupancy with formula random effects via Laplace + a
#' variance-component EM (internal).
#'
#' @param model A single-season `tobs_model`.
#' @param re A list of `tobs_re` specs on the occupancy predictor (iid
#'   intercept or uncorrelated slopes).
#' @param priors Prior spec; not applied on this path (warns when active).
#' @param max_iter,tol,damping EM controls.
#' @param verbose Print per-iteration progress.
#' @keywords internal
.tobs_em_laplace_re <- function(model, re, priors = NULL,
                                max_iter = 100L, tol = 1e-5, damping = 0.3,
                                verbose = TRUE) {
  if (inherits(re, "tobs_re")) re <- list(re)
  design <- .tobs_re_design(re, model)

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

  # State.
  init <- glm_init(X_occ, X_det, any_det, n_det, n_valid, keep, p_occ, p_det)
  beta_occ <- init$occ$beta
  beta_det <- init$det$beta
  n_latent <- sum(vapply(design, function(d) d$n_groups * d$n_coefs, integer(1)))
  b <- numeric(n_latent)
  sig_list <- lapply(design, function(d) rep(0.5, d$n_coefs))

  weights <- NULL
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    # ---- E-step: psi includes the RE posterior mode. ----
    eta_occ <- as.numeric(X_occ %*% beta_occ) + .tobs_re_offset(design, b)
    psi <- plogis(eta_occ)
    p_site <- plogis(as.numeric(X_det %*% beta_det))
    w <- occ_weights(psi, p_site, N, n_valid, n_det, any_det)
    weights <- w

    # ---- M-step (occupancy): pseudo-binomial + rescaled RE prior. ----
    y_occ <- pmin(pmax(ifelse(any_det, M, as.integer(round(w * M))), 0L), M)
    re_list_tulpa <- lapply(seq_along(design), function(k) {
      d <- design[[k]]
      list(idx = d$idx, n_groups = d$n_groups,
           sigma = sig_list[[k]] / sqrt(M), n_coefs = d$n_coefs,
           Z = if (d$n_coefs > 1L) d$Z else NULL)
    })
    fo <- tulpa::tulpa_laplace(
      y = y_occ, n_trials = rep(M, N), X = X_occ,
      re_list = re_list_tulpa, family = "binomial", return_hessian = TRUE)
    beta_new <- fo$mode[seq_len(p_occ)]
    b_new <- fo$mode[-seq_len(p_occ)]

    # ---- Variance-component update at the natural scale. ----
    eta_mode <- as.numeric(X_occ %*% beta_new) + .tobs_re_offset(design, b_new)
    dvar <- .tobs_re_posterior_var(X_occ, eta_mode, design, sig_list)
    sig_new <- sig_list
    pos <- 0L
    for (k in seq_along(design)) {
      d <- design[[k]]; nc <- d$n_coefs; ng <- d$n_groups
      blk_b <- b_new[pos + seq_len(ng * nc)]
      blk_v <- dvar[pos + seq_len(ng * nc)]
      for (c in seq_len(nc)) {
        bc <- blk_b[seq(c, ng * nc, by = nc)]
        vc <- blk_v[seq(c, ng * nc, by = nc)]
        sig_new[[k]][c] <- sqrt(max(mean(bc^2 + vc), 1e-6))
      }
      pos <- pos + ng * nc
    }

    # ---- M-step (detection): weighted binomial, site level. ----
    w_det <- w; w_det[any_det] <- 1
    keep_det <- keep & (w_det > 1e-6)
    fd <- tulpa::tulpa_laplace(
      y = n_det[keep_det], n_trials = n_valid[keep_det],
      X = X_det[keep_det, , drop = FALSE], weights = w_det[keep_det],
      family = "binomial", return_hessian = TRUE)
    beta_det_new <- fd$mode[seq_len(p_det)]

    delta <- max(abs(c(beta_new - beta_occ,
                       beta_det_new - beta_det,
                       unlist(sig_new) - unlist(sig_list))))
    beta_occ <- beta_new; b <- b_new
    beta_det <- beta_det_new; sig_list <- sig_new
    occ_fit <- fo; det_fit <- fd

    if (verbose) cat(sprintf("  RE-EM iter %d: delta = %.6g\n", it, delta))
    if (is.finite(delta) && delta < tol) { converged <- TRUE; break }
  }

  # Final posterior variance of the RE modes (for BLUP SEs and naming) and the
  # natural-scale occupancy fixed-effect SE (the M-step H_beta is M-inflated).
  eta_mode <- as.numeric(X_occ %*% beta_occ) + .tobs_re_offset(design, b)
  dvar <- .tobs_re_posterior_var(X_occ, eta_mode, design, sig_list)
  beta_occ_se <- .tobs_re_occ_fixed_se(X_occ, eta_mode, weights, design, sig_list)

  list(
    fits = list(occ = occ_fit, det = det_fit),
    weights = weights,
    convergence = list(converged = converged, n_iter = it),
    correction = "none",
    re_post = list(design = design, b = b, b_var = dvar, sigma = sig_list,
                   beta_occ_se = beta_occ_se)
  )
}


# Assemble the named RE parameter block (sigma hyperparameters + per-group
# BLUPs) and the structured `re_effects` summary the deterministic RE fit
# appends to a tobs_fit. Single source of truth for the deterministic-path RE
# layout and labels.
#
# Returns a list with `means`, `sds`, `names` (to append to the fixed-effect
# parameter vector) and `re_effects` (per-term BLUP tables consumed by
# ranef.tobs_fit()).
.tobs_re_param_block <- function(re_post) {
  design <- re_post$design
  b      <- re_post$b
  b_var  <- re_post$b_var
  sigma  <- re_post$sigma

  means <- numeric(0); sds <- numeric(0); nms <- character(0)
  re_effects <- list()
  pos <- 0L
  for (k in seq_along(design)) {
    d <- design[[k]]; nc <- d$n_coefs; ng <- d$n_groups
    g <- d$group_label
    # sigma hyperparameters (one per coefficient).
    means <- c(means, sigma[[k]])
    sds   <- c(sds, rep(NA_real_, nc))
    nms   <- c(nms, sprintf("sigma_%s_%s", g, d$coef_names))
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
