# =============================================================================
# fp_occu.R — multistate false-positive occupancy family (Miller et al. 2011)
#
# Detections are classified into three states y in {0, 1, 2}: no detection,
# ambiguous detection (a true detection OR a false positive), and certain /
# confirmed detection (only possible when the site is occupied). The latent
# occupancy z marginalises in closed form (two states), so there is no EM: the
# package-internal fit maximises the exact marginal likelihood directly with an
# analytic gradient (BFGS), and the observed-information covariance is the inverse
# of the negative finite-difference Jacobian of that analytic gradient at the mode
# (exact to finite-difference precision). A NUTS path samples the same marginal.
#
# Four site-level logit arms: occupancy psi (formula), true detection p11
# (detection), false-positive rate p10 (fp_formula), and the probability a true
# detection is certain b (b_formula). Per-site math is src/fp_occu_kernel.h.
#
#   .tobs_build_fp_occu()  data binder -> model_type = "fp_occu"
#   .tobs_fit_fp_occu()    dispatch to the false-positive occupancy fit
#   fp_occu_laplace()      analytic-gradient BFGS fit over the exact marginal
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind a false-positive occupancy model. `y` is an n_sites x J matrix of
# detection states in {0, 1, 2}; NA visits are dropped (like occu()). All four
# arms are site-level: occupancy (occ_formula), true detection (det_formula),
# false-positive rate (fp_formula), and certain-classification (b_formula).
.tobs_build_fp_occu <- function(occ_formula, det_formula, data, y,
                                fp_formula = ~1, b_formula = ~1) {
  if (!is.matrix(y)) {
    stop("y must be a matrix (n_sites x J) of detection states in {0, 1, 2}.",
         call. = FALSE)
  }
  if (nrow(y) != nrow(data)) {
    stop(sprintf("y has %d rows but data has %d rows", nrow(y), nrow(data)),
         call. = FALSE)
  }
  y_int <- matrix(as.integer(round(y)), nrow(y), ncol(y))
  ok <- is.na(y_int) | (y_int %in% 0:2)
  if (!all(ok)) {
    stop("y must contain only the detection states 0 (none), 1 (ambiguous), ",
         "2 (certain), or NA.", call. = FALSE)
  }
  n_sites <- nrow(y_int)

  bind <- .tobs_bind_formulas(list(psi = occ_formula, p11 = det_formula,
                                   p10 = fp_formula, b = b_formula), data)
  X_psi <- model.matrix(bind$fe$psi, data)
  X_p11 <- model.matrix(bind$fe$p11, data)
  X_p10 <- model.matrix(bind$fe$p10, data)
  X_b   <- model.matrix(bind$fe$b,   data)

  # Long form (valid visits only), site-major.
  site_mat <- matrix(seq_len(n_sites), n_sites, ncol(y_int))
  keep <- !is.na(as.vector(t(y_int)))
  y_long   <- as.vector(t(y_int))[keep]
  site_idx <- as.vector(t(site_mat))[keep]

  structure(list(
    model_type = "fp_occu",
    y          = y_int,
    y_long     = as.integer(y_long),
    site_idx   = as.integer(site_idx),
    X_processes = list(X_psi, X_p11, X_p10, X_b),
    formulas   = list(psi = bind$fe$psi, p11 = bind$fe$p11,
                      p10 = bind$fe$p10, b = bind$fe$b),
    structured_terms = bind$terms,
    data       = data,
    n_sites    = n_sites,
    max_visits = ncol(y_int),
    process_info = list(
      list(name = "psi", p = ncol(X_psi), coef_names = colnames(X_psi), link = "logit"),
      list(name = "p11", p = ncol(X_p11), coef_names = colnames(X_p11), link = "logit"),
      list(name = "p10", p = ncol(X_p10), coef_names = colnames(X_p10), link = "logit"),
      list(name = "b",   p = ncol(X_b),   coef_names = colnames(X_b),   link = "logit")
    )
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Fitter (called from .tobs_fit_model for model_type == "fp_occu")
# ---------------------------------------------------------------------------

.tobs_fit_fp_occu <- function(model, max_iter = 200L, tol = 1e-8,
                              sigma.beta = NULL, verbose = TRUE) {
  raw <- fp_occu_laplace(
    y        = model$y_long,
    site_idx = model$site_idx,
    X_psi    = model$X_processes[[1]],
    X_p11    = model$X_processes[[2]],
    X_p10    = model$X_processes[[3]],
    X_b      = model$X_processes[[4]],
    max_iter = as.integer(max_iter), tol = as.numeric(tol),
    sigma_beta = sigma.beta, verbose = isTRUE(verbose))
  build_fp_occu_fit(raw, model)
}


# ---------------------------------------------------------------------------
# Laplace fit (analytic-gradient BFGS over the exact marginal)
# ---------------------------------------------------------------------------

# Central-difference Jacobian of a vector-valued function (the analytic gradient),
# used to assemble the observed information at the mode. Exact to FD precision.
.fp_fd_jacobian <- function(fn, x, h = 1e-5) {
  p <- length(x)
  J <- matrix(0, p, p)
  for (i in seq_len(p)) {
    xp <- x; xm <- x; xp[i] <- xp[i] + h; xm[i] <- xm[i] - h
    J[, i] <- (fn(xp) - fn(xm)) / (2 * h)
  }
  0.5 * (J + t(J))      # symmetrise (the marginal Hessian is symmetric)
}

#' Maximum-likelihood fit of the multistate false-positive occupancy model
#'
#' @description
#' Fits the Miller et al. (2011) false-positive occupancy model with confirmed
#' detections (`y` in `{0, 1, 2}`) by maximising the exact two-state marginal
#' likelihood with an analytic gradient (BFGS). The observed-information
#' covariance is the inverse of the negative finite-difference Jacobian of the
#' analytic gradient at the mode. All four arms (occupancy `psi`, true detection
#' `p11`, false-positive `p10`, certain-classification `b`) are site-level logit
#' predictors.
#'
#' @param y Integer vector of detection states (`0`/`1`/`2`), long form (valid
#'   visits only).
#' @param site_idx Integer vector, same length as `y`, 1-based site index.
#' @param X_psi,X_p11,X_p10,X_b Numeric `[n_sites x p_arm]` design matrices for
#'   the four arms.
#' @param sigma_beta Optional Gaussian prior SD on the coefficients (a mild ridge
#'   for stability); `NULL` (default) is the unpenalised MLE.
#' @param max_iter BFGS iteration budget (default 200).
#' @param tol Convergence tolerance (`optim` `reltol`, default 1e-8).
#' @param verbose Print convergence status.
#'
#' @return A list of class `fp_occu_fit` with `beta_psi`, `beta_p11`,
#'   `beta_p10`, `beta_b`, `log_lik`, `vcov`, `H_obs`, per-site posterior
#'   occupancy `w1`, and `converged`.
#'
#' @references
#' Miller, D. A. W., et al. (2011). Improving occupancy estimation when two types
#'   of observational error occur. *Ecology* 92, 1422-1428.
fp_occu_laplace <- function(y, site_idx, X_psi, X_p11, X_p10, X_b,
                            sigma_beta = NULL, max_iter = 200L, tol = 1e-8,
                            verbose = FALSE) {
  y <- as.integer(y); site_idx <- as.integer(site_idx)
  for (nm in c("X_psi", "X_p11", "X_p10", "X_b")) {
    if (!is.matrix(get(nm))) stop(sprintf("`%s` must be a numeric matrix.", nm), call. = FALSE)
  }
  n_sites <- nrow(X_psi)
  p <- c(ncol(X_psi), ncol(X_p11), ncol(X_p10), ncol(X_b))
  off <- cumsum(c(0L, p))
  idx <- list(psi = off[1] + seq_len(p[1]), p11 = off[2] + seq_len(p[2]),
              p10 = off[3] + seq_len(p[3]), b = off[4] + seq_len(p[4]))

  ridge <- if (is.null(sigma_beta)) 0 else 1 / sigma_beta^2
  eval_cpp <- function(theta) {
    eta_psi <- as.numeric(X_psi %*% theta[idx$psi])
    eta_p11 <- as.numeric(X_p11 %*% theta[idx$p11])
    eta_p10 <- as.numeric(X_p10 %*% theta[idx$p10])
    eta_b   <- as.numeric(X_b   %*% theta[idx$b])
    cpp_fp_occu_total_log_lik(y, site_idx, eta_psi, eta_p11, eta_p10, eta_b)
  }
  grad_design <- function(out, theta) {
    g <- numeric(length(theta))
    g[idx$psi] <- as.numeric(crossprod(X_psi, out$grad_eta_psi))
    g[idx$p11] <- as.numeric(crossprod(X_p11, out$grad_eta_p11))
    g[idx$p10] <- as.numeric(crossprod(X_p10, out$grad_eta_p10))
    g[idx$b]   <- as.numeric(crossprod(X_b,   out$grad_eta_b))
    if (ridge > 0) g <- g - ridge * theta
    g
  }
  neg_ll  <- function(theta) {
    val <- -eval_cpp(theta)$log_lik
    if (ridge > 0) val <- val + 0.5 * ridge * sum(theta^2)
    val
  }
  neg_grad <- function(theta) -grad_design(eval_cpp(theta), theta)

  # Warm start: naive occupancy, p11 ~ 0.5, p10 small, b moderate.
  ny <- tapply(y, factor(site_idx, levels = seq_len(n_sites)), function(v) any(v >= 1))
  naive <- mean(ny[!is.na(ny)])
  theta0 <- numeric(sum(p))
  theta0[idx$psi[1]] <- stats::qlogis(min(max(naive, 0.05), 0.95))
  theta0[idx$p11[1]] <- 0                          # p11 ~ 0.5
  theta0[idx$p10[1]] <- stats::qlogis(0.05)        # small false-positive rate
  theta0[idx$b[1]]   <- 0                          # b ~ 0.5

  # Progress + ETA (gcol33/tulpaObs#43); ON by default. BFGS calls the gradient
  # ~once per quasi-Newton iteration, so ticking there approximates iteration
  # progress (maxit is the ETA denominator); finalised after optim returns.
  .prog <- tulpa:::.tulpa_iter_progress("fp-occu-laplace", as.integer(max_iter), unit = "iter")
  neg_grad_p <- function(theta) { .prog$tick(); neg_grad(theta) }
  opt <- stats::optim(theta0, neg_ll, neg_grad_p, method = "BFGS",
                      control = list(maxit = as.integer(max_iter), reltol = tol))
  .prog$finish()
  theta <- opt$par
  converged <- opt$convergence == 0L
  if (!converged && verbose) {
    warning(sprintf("fp_occu_laplace BFGS did not converge (code %d).", opt$convergence),
            call. = FALSE)
  }

  out  <- eval_cpp(theta)
  Hobs <- -.fp_fd_jacobian(function(th) grad_design(eval_cpp(th), th), theta)
  vcov <- tryCatch(solve(Hobs), error = function(e) matrix(NA_real_, length(theta), length(theta)))

  nm <- c(paste0("psi_", colnames(X_psi)), paste0("p11_", colnames(X_p11)),
          paste0("p10_", colnames(X_p10)), paste0("b_", colnames(X_b)))
  dimnames(vcov) <- list(nm, nm); dimnames(Hobs) <- list(nm, nm)

  structure(list(
    beta_psi = theta[idx$psi], beta_p11 = theta[idx$p11],
    beta_p10 = theta[idx$p10], beta_b = theta[idx$b],
    means = theta, vcov = vcov, H_obs = Hobs,
    log_lik = out$log_lik, log_lik_site = out$log_lik_site, w1 = out$w1,
    converged = converged, n_iter = opt$counts[[1]],
    grad_norm = sqrt(sum(neg_grad(theta)^2)), n_sites = n_sites,
    coef_names = nm), class = c("fp_occu_fit", "list"))
}


# ---------------------------------------------------------------------------
# Grouped random effect on the occupancy (psi) arm (gcol33/tulpaObs#51)
# ---------------------------------------------------------------------------

# AGHQ refinement of an fp_occu fit with a site-level grouped RE on the
# occupancy (psi) arm. The latent z marginalises in closed form to
# L_i = psi_i A_i + (1 - psi_i) B_i, with A_i the occupied-state emission product
# (over the site's visits) and B_i the unoccupied-state product (B_i = 0 when a
# certain detection y = 2 is present). With the false-positive arms (p11, p10, b)
# held at their current values, this is exactly the occupancy occ-arm marginal
# (a Bernoulli-in-psi mixture of two fixed weights), so the psi-arm RE goes
# through the same pure-R `make_site` AGHQ path as occu() (R/re_aghq.R) -- no
# native oracle. tulpa owns the quadrature / mode-finding / log-Cholesky / LKJ /
# marginal Hessian; tulpaObs supplies only the per-site emission products and the
# psi-arm eta-derivatives. `design` is the psi-arm RE design; `beta_*` the warm
# starts. Returns refined estimates or NULL when the pass does not apply.
.tobs_fp_occu_re_aghq <- function(model, design, beta_psi, beta_p11, beta_p10,
                                  beta_b, Sigma_list, b, n_quad = 9L,
                                  lkj_eta = 1.5) {
  idx1 <- as.integer(design[[1]]$idx)
  ng   <- as.integer(design[[1]]$n_groups)
  one_group <- all(vapply(design, function(d)
    identical(as.integer(d$idx), idx1) &&
      identical(as.integer(d$n_groups), ng), logical(1)))
  if (!one_group) return(NULL)
  dtot <- sum(vapply(design, function(d) as.integer(d$n_coefs), integer(1)))
  if (dtot > 3L) return(NULL)

  X_psi <- model$X_processes[[1]]; X_p11 <- model$X_processes[[2]]
  X_p10 <- model$X_processes[[3]]; X_b   <- model$X_processes[[4]]
  p_psi <- ncol(X_psi); p_p11 <- ncol(X_p11)
  p_p10 <- ncol(X_p10); p_b   <- ncol(X_b)
  N <- nrow(X_psi)
  if (any(vapply(design, function(d) length(d$idx) != N, logical(1))) ||
      any(vapply(design, function(d) nrow(d$Z) != N, logical(1)))) {
    return(NULL)
  }
  off <- cumsum(c(0L, p_psi, p_p11, p_p10, p_b))
  i_psi <- off[1] + seq_len(p_psi); i_p11 <- off[2] + seq_len(p_p11)
  i_p10 <- off[3] + seq_len(p_p10); i_b   <- off[4] + seq_len(p_b)

  # Per-site detection-state counts (fixed across the optimization).
  si <- as.integer(model$site_idx); yl <- as.integer(model$y_long)
  n0 <- tabulate(si[yl == 0L], nbins = N)
  n1 <- tabulate(si[yl == 1L], nbins = N)
  n2 <- tabulate(si[yl == 2L], nbins = N)
  n_valid <- n0 + n1 + n2
  bzero   <- n2 > 0L
  keep    <- which(n_valid > 0L)
  cl <- function(e) pmin(pmax(e, -30), 30)

  # make_site(theta) closes over the current betas; the engine supplies the RE
  # offset Z b through the eta passed to deriv / lmat (psi predictor only).
  make_site <- function(theta) {
    bp11 <- theta[i_p11]; bp10 <- theta[i_p10]; bb <- theta[i_b]
    e11 <- cl(as.numeric(X_p11 %*% bp11))
    e10 <- cl(as.numeric(X_p10 %*% bp10))
    eb  <- cl(as.numeric(X_b   %*% bb))
    lp11 <- stats::plogis(e11, log.p = TRUE); l1mp11 <- stats::plogis(-e11, log.p = TRUE)
    lp10 <- stats::plogis(e10, log.p = TRUE); l1mp10 <- stats::plogis(-e10, log.p = TRUE)
    lb   <- stats::plogis(eb,  log.p = TRUE); l1mb   <- stats::plogis(-eb,  log.p = TRUE)
    logA <- (n1 + n2) * lp11 + n0 * l1mp11 + n1 * l1mb + n2 * lb
    logB <- n0 * l1mp10 + n1 * lp10               # unoccupied; B = 0 where bzero
    list(
      eta_re = as.numeric(X_psi %*% theta[i_psi]),
      deriv = function(rows, eta) {
        s <- stats::plogis(eta)
        lA <- logA[rows]; lB <- logB[rows]; bz <- bzero[rows]
        logL <- d1 <- d2 <- numeric(length(rows))
        # certain detection present -> L = psi * A  (occupancy-style branch)
        logL[bz] <- log(s[bz]) + lA[bz]
        d1[bz]   <- 1 - s[bz]
        d2[bz]   <- -s[bz] * (1 - s[bz])
        io <- !bz
        cc <- pmax(lA[io], lB[io])
        A_ <- exp(lA[io] - cc); B_ <- exp(lB[io] - cc)
        sn <- s[io]; sp <- sn * (1 - sn); spp <- sp * (1 - 2 * sn)
        g <- sn * A_ + (1 - sn) * B_; u <- A_ - B_
        dd1 <- u * sp / g
        logL[io] <- cc + log(g)
        d1[io] <- dd1
        d2[io] <- u * spp / g - dd1^2
        list(logL = logL, d1 = d1, d2 = d2)
      },
      lmat = function(rows, ETA) {
        S <- stats::plogis(ETA); lA <- logA[rows]; lB <- logB[rows]; bz <- bzero[rows]
        out <- matrix(0, length(rows), ncol(ETA))
        if (any(bz)) out[bz, ] <- log(S[bz, , drop = FALSE]) + lA[bz]
        if (any(!bz)) {
          lS   <- log(S[!bz, , drop = FALSE]); l1mS <- log1p(-S[!bz, , drop = FALSE])
          t1 <- lS + lA[!bz]; t0 <- l1mS + lB[!bz]
          mx <- pmax(t1, t0)
          out[!bz, ] <- mx + log(exp(t1 - mx) + exp(t0 - mx))
        }
        out
      })
  }

  re_terms <- lapply(design, function(d) list(
    idx = as.integer(d$idx), n_groups = as.integer(d$n_groups),
    n_coefs = as.integer(d$n_coefs),
    Z = if (d$n_coefs > 1L) d$Z else NULL,
    correlated = isTRUE(d$correlated)))

  ref <- tulpa::tulpa_re_aghq(
    theta0 = c(beta_psi, beta_p11, beta_p10, beta_b), re_terms = re_terms,
    Sigma0 = Sigma_list, make_site = make_site, n_obs = N,
    keep = keep, n_quad = n_quad, lkj_eta = lkj_eta)
  if (is.null(ref)) return(NULL)

  bpsi <- ref$theta[i_psi]; bp11 <- ref$theta[i_p11]
  bp10 <- ref$theta[i_p10]; bb <- ref$theta[i_b]
  b_out    <- unlist(lapply(ref$blup,     function(M) as.numeric(t(M))), use.names = FALSE)
  bvar_out <- unlist(lapply(ref$blup_var, function(M) as.numeric(t(M))), use.names = FALSE)

  # Refreshed posterior occupancy w1 at the refined estimate (for fitted()).
  eta_psi <- cl(as.numeric(X_psi %*% bpsi) + .tobs_re_offset(design, b_out))
  e11 <- cl(as.numeric(X_p11 %*% bp11)); e10 <- cl(as.numeric(X_p10 %*% bp10))
  eb  <- cl(as.numeric(X_b %*% bb))
  lp11 <- stats::plogis(e11, log.p = TRUE); l1mp11 <- stats::plogis(-e11, log.p = TRUE)
  lp10 <- stats::plogis(e10, log.p = TRUE); l1mp10 <- stats::plogis(-e10, log.p = TRUE)
  lb   <- stats::plogis(eb,  log.p = TRUE); l1mb   <- stats::plogis(-eb,  log.p = TRUE)
  logA <- (n1 + n2) * lp11 + n0 * l1mp11 + n1 * l1mb + n2 * lb
  logB <- n0 * l1mp10 + n1 * lp10
  s <- stats::plogis(eta_psi)
  w1 <- numeric(N)
  w1[bzero] <- 1
  io <- !bzero
  t1 <- log(s[io]) + logA[io]; t0 <- log1p(-s[io]) + logB[io]
  mx <- pmax(t1, t0)
  w1[io] <- exp(t1 - (mx + log(exp(t1 - mx) + exp(t0 - mx))))

  # Fixed-effect covariance: the full marginal cov when the engine surfaces it,
  # else the diagonal of the per-coefficient marginal SEs (the make_site AGHQ
  # path reports SEs, not the cross-coefficient covariance -- the same form the
  # occupancy RE path uses; no fabricated off-diagonal correlations).
  p_tot <- p_psi + p_p11 + p_p10 + p_b
  vcov <- ref$theta_cov
  if (is.null(vcov) || any(dim(vcov) != p_tot)) {
    se <- ref$theta_se; if (length(se) != p_tot) se <- rep(NA_real_, p_tot)
    vcov <- diag(pmax(se, 0)^2, p_tot)
  }

  list(
    ok = TRUE, arm = "psi",
    beta_psi = bpsi, beta_p11 = bp11, beta_p10 = bp10, beta_b = bb,
    Sigma_list = ref$Sigma_list, b = b_out, b_var = bvar_out,
    theta_se = ref$theta_se, vcov = vcov, w1 = w1,
    log_marginal = ref$log_marginal %||% NA_real_,
    n_quad = ref$n_quad, lkj_eta = ref$lkj_eta, converged = ref$converged)
}


# Fit an fp_occu model with a site-level grouped RE on the occupancy (psi) arm
# under the Laplace / AGHQ path (one grouping factor, RE dim <= 3; tulpaObs#51).
# The false-positive (p10) and certain (b) arms never carry structured terms
# (rejected upstream); a detection (p11) RE is not yet wired here, so this fits
# only a psi-arm RE.
.tobs_fit_fp_occu_re <- function(model, re, max_iter = 200L, tol = 1e-8,
                                 verbose = TRUE, n_quad = 9L, lkj_eta = 1.5,
                                 sigma.beta = NULL) {
  if (inherits(re, "tobs_re")) re <- list(re)
  arms <- .tobs_re_split_two_arms(
    re, model, "psi", "p11",
    "An fp_occu random effect shared across the occupancy and detection arms is not supported.")
  if (length(arms$p11)) {
    stop("fp_occu() random effects are supported on the occupancy (psi) arm ",
         "only; a detection-arm random effect is not yet wired on either ",
         "engine. (tulpaObs#51)", call. = FALSE)
  }
  design <- arms$psi
  if (!length(design)) {
    stop("fp_occu() found no occupancy-arm random effect to fit.", call. = FALSE)
  }

  warm <- tryCatch(
    fp_occu_laplace(y = model$y_long, site_idx = model$site_idx,
                    X_psi = model$X_processes[[1]], X_p11 = model$X_processes[[2]],
                    X_p10 = model$X_processes[[3]], X_b = model$X_processes[[4]],
                    sigma_beta = sigma.beta, max_iter = as.integer(max_iter),
                    tol = as.numeric(tol), verbose = FALSE),
    error = function(e) NULL)
  beta_psi_init <- if (!is.null(warm)) warm$beta_psi else c(0, rep(0, ncol(model$X_processes[[1]]) - 1L))
  beta_p11_init <- if (!is.null(warm)) warm$beta_p11 else rep(0, ncol(model$X_processes[[2]]))
  beta_p10_init <- if (!is.null(warm)) warm$beta_p10 else c(stats::qlogis(0.05), rep(0, ncol(model$X_processes[[3]]) - 1L))
  beta_b_init   <- if (!is.null(warm)) warm$beta_b   else rep(0, ncol(model$X_processes[[4]]))

  Sigma_init <- lapply(design, function(d) diag(0.25, d$n_coefs))
  b_init <- numeric(sum(vapply(design,
                               function(d) as.integer(d$n_groups * d$n_coefs),
                               integer(1))))

  ref <- .tobs_fp_occu_re_aghq(model, design,
                               beta_psi = beta_psi_init, beta_p11 = beta_p11_init,
                               beta_p10 = beta_p10_init, beta_b = beta_b_init,
                               Sigma_list = Sigma_init, b = b_init,
                               n_quad = as.integer(n_quad), lkj_eta = lkj_eta)
  if (is.null(ref) || !isTRUE(ref$ok)) {
    stop("fp_occu() AGHQ random-effect refinement did not produce a usable fit ",
         "(singular marginal Hessian or non-finite optimum). Simplify the RE ",
         "structure.", call. = FALSE)
  }

  raw <- list(
    beta_psi = ref$beta_psi, beta_p11 = ref$beta_p11,
    beta_p10 = ref$beta_p10, beta_b = ref$beta_b,
    means = c(ref$beta_psi, ref$beta_p11, ref$beta_p10, ref$beta_b),
    vcov = ref$vcov, theta_se = ref$theta_se,
    log_lik = ref$log_marginal, w1 = ref$w1,
    converged = ref$converged, n_iter = NA_integer_,
    coef_names = c(paste0("psi_", model$process_info[[1]]$coef_names),
                   paste0("p11_", model$process_info[[2]]$coef_names),
                   paste0("p10_", model$process_info[[3]]$coef_names),
                   paste0("b_",   model$process_info[[4]]$coef_names)))
  re_post <- list(arm = ref$arm, design = design, Sigma_list = ref$Sigma_list,
                  b = ref$b, b_var = ref$b_var,
                  n_quad = ref$n_quad, lkj_eta = ref$lkj_eta)
  build_fp_occu_fit(raw, model, re_post = re_post)
}


# ---------------------------------------------------------------------------
# Fit packer
# ---------------------------------------------------------------------------

build_fp_occu_fit <- function(raw, model, re_post = NULL) {
  pi_list <- model$process_info
  nms <- raw$coef_names
  means <- raw$means; names(means) <- nms
  vcov <- as.matrix(raw$vcov); dimnames(vcov) <- list(nms, nms)
  sds <- sqrt(pmax(diag(vcov), 0)); names(sds) <- nms
  n_fixed <- length(nms); fixed_names <- nms

  n_pseudo <- 1000L
  draws <- .rmvn(n_pseudo, means, vcov); colnames(draws) <- nms
  ll <- raw$log_lik %||% NA_real_

  # Grouped random effect on the occupancy (psi) arm (gcol33/tulpaObs#51):
  # append the variance components (sigma_g_*, cor_g_*_* for a correlated block)
  # and per-group BLUPs after the fixed block, exactly as the count families do.
  # The fixed block (n_fixed leading coords) still governs coef() / vcov() /
  # confint(); ranef() / summary() read the trailing RE columns by name.
  re_block <- NULL
  if (!is.null(re_post) && length(re_post$design)) {
    re_block <- .tobs_re_param_block(list(design = re_post$design,
                                          b      = re_post$b,
                                          b_var  = re_post$b_var,
                                          Sigma  = re_post$Sigma_list))
    means <- c(means, re_block$means); sds <- c(sds, re_block$sds)
    nms   <- c(nms, re_block$names)
    names(means) <- nms; names(sds) <- nms
    draws <- cbind(draws, .tobs_re_pseudo_draws(re_block$means, re_block$sds,
                                                re_block$names, n_pseudo))
  }

  structure(c(list(
    draws = draws, means = means, sds = sds, vcov = vcov,
    n_samples = n_pseudo, n_params = length(means),
    log_prob = rep(ll, n_pseudo),
    N = length(model$y_long)),
    .tobs_na_nuts_diagnostics(n_pseudo),
    list(
    col_names = nms, param_names = nms,
    n_fixed = n_fixed, fixed_names = fixed_names,
    process_info = pi_list,
    model = model, spatial = NULL, method = "laplace",
    log_lik = ll, w1 = raw$w1,
    re_effects = re_block$re_effects,
    fp_re = if (!is.null(re_post))
      list(arm = re_post$arm, n_quad = re_post$n_quad,
           lkj_eta = re_post$lkj_eta, Sigma_list = re_post$Sigma_list)
      else NULL,
    convergence = list(converged = raw$converged %||% TRUE,
                       n_iter = raw$n_iter %||% NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "fp_occu")
# ---------------------------------------------------------------------------

# Per-site fitted occupancy psi, true detection p11, false-positive p10, certain
# prob b, and posterior occupancy z = P(z = 1 | y).
.tobs_fitted_fp_occu <- function(object) {
  model <- object$model; means <- object$means
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  lp <- function(k) stats::plogis(as.vector(model$X_processes[[k]] %*%
                                            means[off[k] + seq_len(p[k])]))
  list(psi = lp(1), p11 = lp(2), p10 = lp(3), b = lp(4), z = object$w1)
}

# simulate() for fp_occu: draw z, then per visit emit a multistate detection.
.tobs_simulate_fp_occu <- function(object, nsim = 1) {
  model <- object$model; draws <- object$draws; n_draws <- nrow(draws)
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  n_sites <- model$n_sites; J <- model$max_visits
  result <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    di <- sample.int(n_draws, 1L)
    th <- draws[di, ]
    psi <- stats::plogis(as.vector(model$X_processes[[1]] %*% th[off[1] + seq_len(p[1])]))
    p11 <- stats::plogis(as.vector(model$X_processes[[2]] %*% th[off[2] + seq_len(p[2])]))
    p10 <- stats::plogis(as.vector(model$X_processes[[3]] %*% th[off[3] + seq_len(p[3])]))
    b   <- stats::plogis(as.vector(model$X_processes[[4]] %*% th[off[4] + seq_len(p[4])]))
    z <- stats::rbinom(n_sites, 1L, psi)
    y_sim <- matrix(0L, n_sites, J)
    for (i in seq_len(n_sites)) {
      if (z[i] == 1L) {
        det <- stats::rbinom(J, 1L, p11[i])
        cert <- stats::rbinom(J, 1L, b[i])
        y_sim[i, ] <- ifelse(det == 1L, ifelse(cert == 1L, 2L, 1L), 0L)
      } else {
        y_sim[i, ] <- stats::rbinom(J, 1L, p10[i])    # false positives -> state 1
      }
    }
    result[[s]] <- y_sim
  }
  if (nsim == 1L) result[[1]] else result
}

# residuals() for fp_occu, on the marginal probability of any detection
# (y >= 1): P(y_ij >= 1) = psi * p11 + (1 - psi) * p10.
.tobs_residuals_fp_occu <- function(object, type = c("deviance", "pearson",
                                                   "response")) {
  type  <- match.arg(type)
  fitv  <- .tobs_fitted_fp_occu(object)
  model <- object$model
  mu <- fitv$psi * fitv$p11 + (1 - fitv$psi) * fitv$p10
  y_any <- (model$y >= 1) * 1
  mu_mat <- matrix(mu, model$n_sites, model$max_visits)
  mu_mat <- pmin(pmax(mu_mat, 1e-10), 1 - 1e-10)
  switch(type,
    response = y_any - mu_mat,
    pearson  = (y_any - mu_mat) / sqrt(mu_mat * (1 - mu_mat)),
    deviance = {
      d <- 2 * (ifelse(y_any == 1, log(1 / mu_mat), log(1 / (1 - mu_mat))))
      sign(y_any - mu_mat) * sqrt(pmax(d, 0))
    })
}

# predict() for fp_occu: occupancy psi (default) or true detection p11 at new X.
.tobs_predict_fp_occu <- function(object, X.0 = NULL, type = c("psi", "p11")) {
  type  <- match.arg(type)
  model <- object$model
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  k <- if (identical(type, "psi")) 1L else 2L
  X <- X.0 %||% model$X_processes[[k]]
  stats::plogis(as.vector(X %*% object$means[off[k] + seq_len(p[k])]))
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate multistate false-positive occupancy data
#'
#' Latent occupancy `z_i ~ Bernoulli(psi_i)` with `logit psi = X_psi beta_psi`,
#' observed through `J` replicate visits emitting a detection state `y in {0,1,2}`:
#' at an occupied site a visit detects with probability `p11` and the detection is
#' certain (state 2) with probability `b` else ambiguous (state 1); at an
#' unoccupied site a visit yields a false-positive ambiguous detection (state 1)
#' with probability `p10`. Returns an `N x J` integer matrix in `{0, 1, 2}`
#' suitable for [tobs()] with [fp_occu()].
#'
#' @param N Number of sites (default 300).
#' @param J Number of visits (default 5).
#' @param n_occ_covs Number of occupancy covariates (default 1).
#' @param beta_psi Occupancy coefficients (logit). Default
#'   `c(qlogis(0.5), runif(n_occ_covs, -0.6, 0.6))`.
#' @param p11,p10,b True detection, false-positive, and certain-classification
#'   probabilities (scalars; defaults 0.6, 0.05, 0.5).
#' @param seed Optional random seed.
#' @return A list with `y` (N x J state matrix), `data` (occupancy covariates),
#'   and `truth` (coefficients, per-site `psi`, scalar `p11`/`p10`/`b`, latent `z`).
#' @export
simulate_fp_occu <- function(N = 300, J = 5, n_occ_covs = 1, beta_psi = NULL,
                             p11 = 0.6, p10 = 0.05, b = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_psi)) beta_psi <- c(stats::qlogis(0.5), stats::runif(n_occ_covs, -0.6, 0.6))
  occ_covs <- data.frame(matrix(stats::rnorm(N * n_occ_covs), N, n_occ_covs))
  names(occ_covs) <- paste0("occ_cov", seq_len(n_occ_covs))
  X_psi <- stats::model.matrix(~ ., occ_covs)
  psi <- stats::plogis(as.vector(X_psi %*% beta_psi))
  z <- stats::rbinom(N, 1L, psi)
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) {
    if (z[i] == 1L) {
      det <- stats::rbinom(J, 1L, p11); cert <- stats::rbinom(J, 1L, b)
      y[i, ] <- ifelse(det == 1L, ifelse(cert == 1L, 2L, 1L), 0L)
    } else {
      y[i, ] <- stats::rbinom(J, 1L, p10)
    }
  }
  list(y = y, data = occ_covs,
       truth = list(beta_psi = beta_psi, psi = psi, p11 = p11, p10 = p10,
                    b = b, z = z))
}
