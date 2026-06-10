# =============================================================================
# dyn_abun.R — Dail-Madsen open-population N-mixture family
#
# Latent abundance evolves across primary seasons: N_1 ~ Poisson(lambda); for
# t >= 2, N_t = survivors S_t ~ Binomial(N_{t-1}, omega) plus recruits
# G_t ~ Poisson(gamma); observed via Binomial(N_t, p) over secondary visits. The
# latent abundance sequence is NOT closed-form (unlike the static N-mixture) — it
# is summed out by an exact HMM forward recursion over the abundance states
# (src/dyn_abun_kernel.h), with analytic gradients by forward-mode differentiation
# of the scaled forward algorithm. The fit maximises the exact marginal with that
# analytic gradient (BFGS) and an observed-information covariance; a NUTS path
# samples the same marginal.
#
# Four site-level arms: initial abundance lambda (log; the formula), detection p
# (logit; detection), apparent survival omega (logit; omega_formula), and
# recruitment gamma (log; gamma_formula).
#
#   .tobs_build_dyn_abun()  data binder -> model_type = "dyn_abun"
#   .tobs_fit_dyn_abun()    dispatch to the open N-mixture fit
#   dyn_abun_laplace()      analytic-gradient BFGS over the forward marginal
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Build the survival / recruitment arm design. The arm spans the T-1 transition
# intervals; a covariate that varies by interval is supplied as a matrix column
# of `data` with one column per interval ([n_sites x (T-1)]). When the arm
# references such a column the design is long-form [(site x interval) x p],
# site-major interval-minor (row (i-1)*(T-1) + iv); otherwise it collapses to the
# per-site design [n_sites x p] (the constant-rate path). `arm` names the arm for
# error messages. Returns the model matrix and a season_varying flag.
.tobs_dyn_abun_arm_design <- function(fe_formula, data, n_sites, n_intervals,
                                      arm) {
  vars <- all.vars(fe_formula)
  vars <- vars[vars %in% names(data)]
  is_sv <- vapply(vars, function(v) {
    col <- data[[v]]
    is.matrix(col) && ncol(col) == n_intervals
  }, logical(1))
  bad <- vapply(vars, function(v) {
    col <- data[[v]]
    is.matrix(col) && ncol(col) != n_intervals
  }, logical(1))
  if (any(bad)) {
    stop(sprintf(paste0("dyn_abun() %s covariate '%s' is a matrix with %d ",
                        "columns; a season-varying %s covariate must have one ",
                        "column per transition interval (T-1 = %d)."),
                 arm, vars[bad][1], ncol(data[[vars[bad][1]]]), arm, n_intervals),
         call. = FALSE)
  }
  if (!any(is_sv)) {
    return(list(X = model.matrix(fe_formula, data), season_varying = FALSE))
  }
  # Unroll to long form: season-varying matrix columns flatten site-major
  # interval-minor; per-site columns repeat across the site's intervals.
  rep_site <- rep(seq_len(n_sites), each = n_intervals)
  iv_idx   <- rep(seq_len(n_intervals), times = n_sites)
  long <- list()
  for (v in names(data)) {
    col <- data[[v]]
    if (is.matrix(col) && ncol(col) == n_intervals) {
      long[[v]] <- col[cbind(rep_site, iv_idx)]
    } else if (!is.matrix(col)) {
      long[[v]] <- col[rep_site]
    }
  }
  data_long <- as.data.frame(long, stringsAsFactors = FALSE)
  list(X = model.matrix(fe_formula, data_long), season_varying = TRUE)
}

# Bind a Dail-Madsen open N-mixture. `y` is a 3D array [n_sites x max_visits x
# n_seasons] (or a list of n_sites x max_visits matrices, one per season); a
# missing visit is NA (dropped from its season). Initial abundance (occ_formula =
# the tobs formula) and detection (det_formula) are site-level; survival
# (omega_formula) and recruitment (gamma_formula) span the T-1 transition
# intervals and may carry a season-varying covariate (a [n_sites x (T-1)] matrix
# column of `data`), in which case their design is long-form over (site, interval).
.tobs_build_dyn_abun <- function(occ_formula, det_formula, data, y,
                                 omega_formula = ~1, gamma_formula = ~1,
                                 mixture = "poisson", K_max = NULL) {
  if (is.list(y) && !is.array(y)) {
    n_seasons <- length(y)
    n_sites <- nrow(y[[1]]); max_visits <- ncol(y[[1]])
    ya <- array(NA_integer_, dim = c(n_sites, max_visits, n_seasons))
    for (t in seq_len(n_seasons)) ya[, , t] <- as.matrix(y[[t]])
    y <- ya
  }
  if (length(dim(y)) != 3) {
    stop("y must be a 3D array [n_sites x max_visits x n_seasons] or a list of ",
         "per-season matrices.", call. = FALSE)
  }
  n_sites <- dim(y)[1]; max_visits <- dim(y)[2]; n_seasons <- dim(y)[3]
  if (n_sites != nrow(data)) {
    stop(sprintf("y has %d sites but data has %d rows", n_sites, nrow(data)),
         call. = FALSE)
  }
  if (n_seasons < 2L) {
    stop("dyn_abun() needs >= 2 primary seasons; for a single season use abun().",
         call. = FALSE)
  }
  if (!mixture %in% c("poisson", "negbin")) {
    stop("dyn_abun(mixture = '", mixture, "') is not supported. ",
         "Use 'poisson' or 'negbin'.", call. = FALSE)
  }

  bind <- .tobs_bind_formulas(list(lambda = occ_formula, p = det_formula,
                                   omega = omega_formula, gamma = gamma_formula), data)
  X_lambda <- model.matrix(bind$fe$lambda, data)
  X_p      <- model.matrix(bind$fe$p, data)
  n_intervals <- n_seasons - 1L
  om <- .tobs_dyn_abun_arm_design(bind$fe$omega, data, n_sites, n_intervals, "omega")
  gm <- .tobs_dyn_abun_arm_design(bind$fe$gamma, data, n_sites, n_intervals, "gamma")
  X_omega <- om$X; X_gamma <- gm$X

  # Flat layout site-major then season then visit: j + J*t + J*T*i (0-based),
  # matching compute_dyn_abun_site (dyn_occu's aperm(y, c(2,3,1)) convention).
  y_flat <- as.integer(aperm(y, c(2, 3, 1)))
  y_flat[is.na(y_flat)] <- -1L

  obs <- y[!is.na(y)]
  if (any(obs < 0)) stop("y must contain nonnegative integer counts.", call. = FALSE)
  K <- if (is.null(K_max)) as.integer(max(obs) + 40L) else as.integer(K_max)

  structure(list(
    model_type = "dyn_abun",
    y          = y,
    y_flat     = y_flat,
    X_processes = list(X_lambda, X_p, X_omega, X_gamma),
    formulas   = list(lambda = bind$fe$lambda, p = bind$fe$p,
                      omega = bind$fe$omega, gamma = bind$fe$gamma),
    structured_terms = bind$terms,
    data       = data,
    n_sites    = n_sites, n_seasons = n_seasons, max_visits = max_visits,
    n_intervals = n_intervals,
    omega_season_varying = om$season_varying,
    gamma_season_varying = gm$season_varying,
    K_max      = K, mixture = mixture,
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda), coef_names = colnames(X_lambda), link = "log"),
      list(name = "p",      p = ncol(X_p),      coef_names = colnames(X_p),      link = "logit"),
      list(name = "omega",  p = ncol(X_omega),  coef_names = colnames(X_omega),  link = "logit"),
      list(name = "gamma",  p = ncol(X_gamma),  coef_names = colnames(X_gamma),  link = "log")
    ),
    mean_count = mean(obs), max_count = max(obs)
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Fitter (called from .tobs_fit_model for model_type == "dyn_abun")
# ---------------------------------------------------------------------------

.tobs_fit_dyn_abun <- function(model, max_iter = 300L, tol = 1e-8, verbose = TRUE) {
  raw <- dyn_abun_laplace(
    y_flat = model$y_flat, n_sites = model$n_sites, T = model$n_seasons,
    J = model$max_visits, K_max = model$K_max,
    X_lambda = model$X_processes[[1]], X_p = model$X_processes[[2]],
    X_omega = model$X_processes[[3]], X_gamma = model$X_processes[[4]],
    mixture = model$mixture %||% "poisson",
    max_iter = as.integer(max_iter), tol = as.numeric(tol), verbose = isTRUE(verbose))
  build_dyn_abun_fit(raw, model)
}


# ---------------------------------------------------------------------------
# Grouped random effect on the initial-abundance arm (gcol33/tulpaObs#51)
# ---------------------------------------------------------------------------

# AGHQ refinement of a dyn_abun fit with a site-level grouped RE on the
# initial-abundance (lambda) arm OR the detection (p) arm. The Dail-Madsen
# marginal is NOT a closed-form per-site mixture (unlike the static N-mixture /
# fp_occu families), so it does not factorise into the closed-form A / B weights
# the count oracle and the fp_occu make_site path use. But a single grouping
# factor's RE shifts a per-site SCALAR offset on one arm, and the per-site log
# marginal (and its first and second derivatives in that offset) propagate through
# the exact HMM forward recursion (src/dyn_abun_kernel.h) -- exactly the per-row
# separability the AGHQ engine's make_site contract needs. tulpa owns the
# quadrature / per-group mode / log-Cholesky / LKJ / marginal Hessian; tulpaObs
# supplies the per-site marginal and its eta-derivatives.
#
# The two arms differ only in WHERE the offset enters and HOW cheaply the marginal
# re-evaluates per quadrature node:
#   - lambda (tulpaObs#51): eta_lambda enters ONLY the season-1 initial
#     distribution, so the data-conditional weights c(n1) = P(all data | N_1 = n1)
#     are eta-independent and precomputed ONCE (the O(K^2 T) HMM backward pass,
#     cpp_dyn_abun_init_weights_mat); each node is an O(K) dot (cpp_dyn_abun_init_-
#     loglik).
#   - p (tulpaObs#82): eta_p enters every season's observation pmf, so c cannot be
#     precomputed -- each node re-evaluates the full O(K^2 T) forward marginal with
#     a closed-form second-order eta_p forward-mode pass (cpp_dyn_abun_p_loglik).
# The non-offset arms and the dispersion are held fixed during the integration
# (closed over per make_site call). `arm` is "lambda" or "p"; `design` the
# corresponding RE design; `beta_*` / `log_r` the warm starts. Returns refined
# estimates or NULL when the pass does not apply (caller keeps the no-RE fit).
.tobs_dyn_abun_re_aghq <- function(model, design, beta_lambda, beta_p, beta_omega,
                                   beta_gamma, log_r = NA_real_, Sigma_list, b,
                                   arm = "lambda", n_quad = 9L, lkj_eta = 1.5,
                                   theta_prior_sd = 100) {
  idx1 <- as.integer(design[[1]]$idx)
  ng   <- as.integer(design[[1]]$n_groups)
  one_group <- all(vapply(design, function(d)
    identical(as.integer(d$idx), idx1) &&
      identical(as.integer(d$n_groups), ng), logical(1)))
  if (!one_group) return(NULL)
  dtot <- sum(vapply(design, function(d) as.integer(d$n_coefs), integer(1)))
  if (dtot > 3L) return(NULL)

  X_lambda <- model$X_processes[[1]]; X_p <- model$X_processes[[2]]
  X_omega  <- model$X_processes[[3]]; X_gamma <- model$X_processes[[4]]
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  p_om  <- ncol(X_omega);  p_gm <- ncol(X_gamma)
  N <- model$n_sites
  if (any(vapply(design, function(d) length(d$idx) != N, logical(1))) ||
      any(vapply(design, function(d) nrow(d$Z) != N, logical(1)))) {
    return(NULL)
  }

  T <- model$n_seasons; J <- model$max_visits; K <- model$K_max
  y_flat <- as.integer(model$y_flat)
  use_nb <- identical(model$mixture %||% "poisson", "negbin")
  off <- cumsum(c(0L, p_lam, p_p, p_om, p_gm))
  i_lam <- off[1] + seq_len(p_lam); i_p  <- off[2] + seq_len(p_p)
  i_om  <- off[3] + seq_len(p_om);  i_gm <- off[4] + seq_len(p_gm)
  i_logr <- if (use_nb) off[5] + 1L else NA_integer_
  nIv <- T - 1L
  cl <- function(e) pmin(pmax(e, -30), 30)
  keep <- seq_len(N)

  # Subset a (possibly season-varying) survival / recruitment arm eta to a row
  # subset, in the site-major interval-minor order cpp_dyn_abun_p_loglik expects.
  sub_arm <- function(e, rows) {
    if (length(e) == N) return(e[rows])
    as.numeric(t(matrix(e, nrow = N, ncol = nIv, byrow = TRUE)[rows, , drop = FALSE]))
  }

  if (identical(arm, "lambda")) {
    # eta_lambda enters only the season-1 initial distribution; c(n1) precomputed
    # ONCE, each node an O(K) dot.
    make_site <- function(theta) {
      e_p  <- cl(as.numeric(X_p %*% theta[i_p]))
      e_om <- cl(as.numeric(X_omega %*% theta[i_om]))
      e_gm <- cl(as.numeric(X_gamma %*% theta[i_gm]))
      elr  <- if (use_nb) as.numeric(theta[i_logr]) else 0
      Cmat <- cpp_dyn_abun_init_weights_mat(
        y_flat, N, T, J, K, site = as.integer(seq_len(N) - 1L),
        eta_p = e_p, eta_omega = e_om, eta_gamma = e_gm)
      list(
        eta_re = as.numeric(X_lambda %*% theta[i_lam]),
        deriv = function(rows, eta) {
          r <- cpp_dyn_abun_init_loglik(
            Cmat[rows, , drop = FALSE], as.numeric(eta),
            use_nb = use_nb, eta_logr = elr, deriv = TRUE)
          list(logL = r$logL, d1 = r$d1, d2 = r$d2)
        },
        lmat = function(rows, ETA) {
          Csub <- Cmat[rows, , drop = FALSE]
          out <- matrix(0, length(rows), ncol(ETA))
          for (cc in seq_len(ncol(ETA))) {
            out[, cc] <- cpp_dyn_abun_init_loglik(
              Csub, as.numeric(ETA[, cc]), use_nb = use_nb, eta_logr = elr,
              deriv = FALSE)$logL
          }
          out
        })
    }
  } else {
    # eta_p enters every season's observation pmf; the full forward marginal is
    # re-evaluated per node (closed-form second-order eta_p forward-mode pass). The
    # initial-abundance / survival / recruitment predictors are held fixed.
    make_site <- function(theta) {
      e_lam <- cl(as.numeric(X_lambda %*% theta[i_lam]))
      e_om  <- cl(as.numeric(X_omega %*% theta[i_om]))
      e_gm  <- cl(as.numeric(X_gamma %*% theta[i_gm]))
      elr   <- if (use_nb) as.numeric(theta[i_logr]) else 0
      list(
        eta_re = as.numeric(X_p %*% theta[i_p]),
        deriv = function(rows, eta) {
          r <- cpp_dyn_abun_p_loglik(
            y_flat, N, T, J, K, site = as.integer(rows - 1L),
            eta_lambda = e_lam[rows], eta_p = as.numeric(eta),
            eta_omega = sub_arm(e_om, rows), eta_gamma = sub_arm(e_gm, rows),
            use_nb = use_nb, eta_logr = elr, deriv = TRUE)
          list(logL = r$logL, d1 = r$d1, d2 = r$d2)
        },
        lmat = function(rows, ETA) {
          el <- e_lam[rows]; eo <- sub_arm(e_om, rows); eg <- sub_arm(e_gm, rows)
          si <- as.integer(rows - 1L)
          out <- matrix(0, length(rows), ncol(ETA))
          for (cc in seq_len(ncol(ETA))) {
            out[, cc] <- cpp_dyn_abun_p_loglik(
              y_flat, N, T, J, K, site = si, eta_lambda = el,
              eta_p = as.numeric(ETA[, cc]), eta_omega = eo, eta_gamma = eg,
              use_nb = use_nb, eta_logr = elr, deriv = FALSE)$logL
          }
          out
        })
    }
  }

  re_terms <- lapply(design, function(d) list(
    idx = as.integer(d$idx), n_groups = as.integer(d$n_groups),
    n_coefs = as.integer(d$n_coefs),
    Z = if (d$n_coefs > 1L) d$Z else NULL,
    correlated = isTRUE(d$correlated)))

  theta0 <- c(beta_lambda, beta_p, beta_omega, beta_gamma)
  if (use_nb) theta0 <- c(theta0, if (is.finite(log_r)) log_r else log(2))

  ref <- tulpa::tulpa_re_aghq(
    theta0 = theta0, re_terms = re_terms, Sigma0 = Sigma_list,
    make_site = make_site, n_obs = N, keep = keep,
    n_quad = as.integer(n_quad), lkj_eta = lkj_eta,
    theta_prior_sd = theta_prior_sd)
  if (is.null(ref)) return(NULL)

  bl <- ref$theta[i_lam]; bp <- ref$theta[i_p]
  bo <- ref$theta[i_om];  bg <- ref$theta[i_gm]
  log_r_ref <- if (use_nb) ref$theta[i_logr] else NA_real_
  b_out    <- unlist(lapply(ref$blup,     function(M) as.numeric(t(M))), use.names = FALSE)
  bvar_out <- unlist(lapply(ref$blup_var, function(M) as.numeric(t(M))), use.names = FALSE)

  # Marginal fixed-effect covariance when the engine surfaces it, else the
  # diagonal of the per-coefficient marginal SEs (the make_site AGHQ path reports
  # SEs, not the cross-coefficient covariance -- as the occupancy / fp_occu RE
  # paths do; no fabricated off-diagonal correlations).
  p_tot <- p_lam + p_p + p_om + p_gm + if (use_nb) 1L else 0L
  vcov <- ref$theta_cov
  if (is.null(vcov) || any(dim(vcov) != p_tot)) {
    se <- ref$theta_se; if (length(se) != p_tot) se <- rep(NA_real_, p_tot)
    vcov <- diag(pmax(se, 0)^2, p_tot)
  }

  # mean_N1 / log_lik at the refined estimate + posterior-mode RE offset (for
  # fitted()); the offset rides whichever arm carries the term. The reported
  # marginal log-likelihood is the engine's integrated value.
  off_re <- .tobs_re_offset(design, b_out)
  eta_lambda <- cl(as.numeric(X_lambda %*% bl) +
                     if (identical(arm, "lambda")) off_re else 0)
  eta_p_vec  <- as.numeric(X_p %*% bp) + if (identical(arm, "p")) off_re else 0
  ev <- cpp_dyn_abun_total_log_lik(
    y_flat, N, T, J, K, eta_lambda, eta_p_vec,
    as.numeric(X_omega %*% bo), as.numeric(X_gamma %*% bg),
    use_nb = use_nb, eta_logr = if (use_nb) log_r_ref else 0)

  list(
    ok = TRUE, arm = arm,
    beta_lambda = bl, beta_p = bp, beta_omega = bo, beta_gamma = bg,
    log_r = log_r_ref, r = if (use_nb) exp(log_r_ref) else NA_real_,
    mixture = model$mixture %||% "poisson",
    Sigma_list = ref$Sigma_list, b = b_out, b_var = bvar_out,
    theta_se = ref$theta_se, vcov = vcov,
    mean_N1 = ev$mean_N1, log_marginal = ref$log_marginal %||% NA_real_,
    n_quad = ref$n_quad, lkj_eta = ref$lkj_eta, converged = ref$converged)
}


# Fit a Dail-Madsen open N-mixture with a site-level grouped RE on the
# initial-abundance (lambda, tulpaObs#51) OR the detection (p, tulpaObs#82) arm
# under the Laplace / AGHQ path (one grouping factor, RE dim <= 3). The survival
# (omega) and recruitment (gamma) arms never carry structured terms (rejected
# upstream -- they are processes > 2). RE on BOTH arms at once is rejected: the
# AGHQ path integrates one arm at a time.
.tobs_fit_dyn_abun_re <- function(model, re, max_iter = 300L, tol = 1e-8,
                                  verbose = TRUE, n_quad = 9L, lkj_eta = 1.5,
                                  theta_prior_sd = 100) {
  if (inherits(re, "tobs_re")) re <- list(re)
  arms <- .tobs_re_split_two_arms(
    re, model, "lambda", "p",
    paste0("A dyn_abun random effect shared across the initial-abundance and ",
           "detection arms is not supported; fit a separate RE per arm."))
  if (length(arms$lambda) && length(arms$p)) {
    stop("Random effects on BOTH the initial-abundance (lambda) and detection ",
         "(p) arms in one dyn_abun fit are not supported; the AGHQ path ",
         "integrates one arm at a time. Put the RE on lambda OR p. (tulpaObs#82)",
         call. = FALSE)
  }
  arm    <- if (length(arms$lambda)) "lambda" else "p"
  design <- if (length(arms$lambda)) arms$lambda else arms$p
  if (!length(design)) {
    stop("dyn_abun() found no initial-abundance or detection random effect to fit.",
         call. = FALSE)
  }

  warm <- tryCatch(
    dyn_abun_laplace(
      y_flat = model$y_flat, n_sites = model$n_sites, T = model$n_seasons,
      J = model$max_visits, K_max = model$K_max,
      X_lambda = model$X_processes[[1]], X_p = model$X_processes[[2]],
      X_omega = model$X_processes[[3]], X_gamma = model$X_processes[[4]],
      mixture = model$mixture %||% "poisson",
      max_iter = as.integer(max_iter), tol = as.numeric(tol), verbose = FALSE),
    error = function(e) NULL)
  beta_lambda_init <- if (!is.null(warm)) warm$beta_lambda
                      else c(log(max(mean(model$mean_count %||% 1), 0.5) + 0.5),
                             rep(0, ncol(model$X_processes[[1]]) - 1L))
  beta_p_init     <- if (!is.null(warm)) warm$beta_p     else rep(0, ncol(model$X_processes[[2]]))
  beta_omega_init <- if (!is.null(warm)) warm$beta_omega else c(stats::qlogis(0.6), rep(0, ncol(model$X_processes[[3]]) - 1L))
  beta_gamma_init <- if (!is.null(warm)) warm$beta_gamma else c(log(0.5), rep(0, ncol(model$X_processes[[4]]) - 1L))
  log_r_init <- if (!is.null(warm) && is.finite(warm$log_r %||% NA_real_)) warm$log_r else NA_real_

  Sigma_init <- lapply(design, function(d) diag(0.25, d$n_coefs))
  b_init <- numeric(sum(vapply(design,
                               function(d) as.integer(d$n_groups * d$n_coefs),
                               integer(1))))

  ref <- .tobs_dyn_abun_re_aghq(
    model, design, beta_lambda = beta_lambda_init, beta_p = beta_p_init,
    beta_omega = beta_omega_init, beta_gamma = beta_gamma_init, log_r = log_r_init,
    Sigma_list = Sigma_init, b = b_init, arm = arm,
    n_quad = as.integer(n_quad), lkj_eta = lkj_eta, theta_prior_sd = theta_prior_sd)
  if (is.null(ref) || !isTRUE(ref$ok)) {
    stop("dyn_abun() AGHQ random-effect refinement did not produce a usable fit ",
         "(singular marginal Hessian or non-finite optimum). Try a different ",
         "K_max or simplify the RE structure.", call. = FALSE)
  }

  is_nb <- identical(ref$mixture, "negbin")
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("p_",      model$process_info[[2]]$coef_names),
           paste0("omega_",  model$process_info[[3]]$coef_names),
           paste0("gamma_",  model$process_info[[4]]$coef_names))
  means <- c(ref$beta_lambda, ref$beta_p, ref$beta_omega, ref$beta_gamma)
  if (is_nb) { nms <- c(nms, "log_r"); means <- c(means, ref$log_r) }
  raw <- list(
    beta_lambda = ref$beta_lambda, beta_p = ref$beta_p,
    beta_omega = ref$beta_omega, beta_gamma = ref$beta_gamma,
    log_r = ref$log_r, r = ref$r, mixture = ref$mixture,
    means = means, vcov = ref$vcov, log_lik = ref$log_marginal,
    mean_N1 = ref$mean_N1, K_max = model$K_max,
    converged = ref$converged, n_iter = NA_integer_, coef_names = nms)
  re_post <- list(arm = ref$arm, design = design, Sigma_list = ref$Sigma_list,
                  b = ref$b, b_var = ref$b_var,
                  n_quad = ref$n_quad, lkj_eta = ref$lkj_eta)
  build_dyn_abun_fit(raw, model, re_post = re_post)
}


# ---------------------------------------------------------------------------
# Laplace fit (analytic-gradient BFGS over the forward marginal)
# ---------------------------------------------------------------------------

#' Maximum-likelihood fit of the Dail-Madsen open N-mixture
#'
#' @description
#' Fits the Dail-Madsen (2011) open-population N-mixture (Poisson initial
#' abundance, binomial survival, Poisson recruitment, binomial detection) by
#' maximising the exact HMM forward marginal likelihood with an analytic gradient
#' (BFGS). The observed-information covariance is the inverse of the negative
#' finite-difference Jacobian of the analytic gradient at the mode. All four arms
#' (initial abundance `lambda`, detection `p`, survival `omega`, recruitment
#' `gamma`) are site-level.
#'
#' @param y_flat Integer vector, flattened detection counts (site-major, then
#'   season, then visit; `-1` marks a missing visit).
#' @param n_sites,T,J Sites, primary seasons, secondary visits per season.
#' @param K_max Abundance-state truncation (states `0..K_max`).
#' @param X_lambda,X_p,X_omega,X_gamma Numeric `[n_sites x p_arm]` design matrices.
#' @param mixture Initial-abundance distribution: `"poisson"` (default) or
#'   `"negbin"` (negative-binomial).
#' @param log_r_init Optional starting value (log scale) for the
#'   negative-binomial size parameter `r`; used only when `mixture = "negbin"`.
#'   Default `NULL` starts from `log(2)`.
#' @param max_iter BFGS iteration budget (default 300).
#' @param tol Convergence tolerance (`optim` `reltol`, default 1e-8).
#' @param verbose Print convergence status.
#'
#' @return A list of class `dyn_abun_fit` with `beta_lambda`, `beta_p`,
#'   `beta_omega`, `beta_gamma`, `log_lik`, `vcov`, `H_obs`, per-site `mean_N1`,
#'   and `converged`.
#'
#' @references Dail, D., Madsen, L. (2011). Models for estimating abundance from
#'   repeated counts of an open metapopulation. *Biometrics* 67, 577-587.
dyn_abun_laplace <- function(y_flat, n_sites, T, J, K_max,
                             X_lambda, X_p, X_omega, X_gamma,
                             mixture = "poisson", log_r_init = NULL,
                             max_iter = 300L, tol = 1e-8, verbose = FALSE) {
  y_flat <- as.integer(y_flat)
  K_max <- as.integer(K_max)
  is_nb <- identical(mixture, "negbin")
  p <- c(ncol(X_lambda), ncol(X_p), ncol(X_omega), ncol(X_gamma))
  off <- cumsum(c(0L, p))
  idx <- list(lambda = off[1] + seq_len(p[1]), p = off[2] + seq_len(p[2]),
              omega = off[3] + seq_len(p[3]), gamma = off[4] + seq_len(p[4]))
  # Under NB the dispersion log r is a single trailing coordinate (shared across
  # sites), jointly estimated with the betas like abun()'s nmix dispersion.
  n_par  <- sum(p) + if (is_nb) 1L else 0L
  ir     <- if (is_nb) n_par else NA_integer_

  eval_cpp <- function(theta) {
    cpp_dyn_abun_total_log_lik(
      y_flat, n_sites, T, J, K_max,
      as.numeric(X_lambda %*% theta[idx$lambda]),
      as.numeric(X_p      %*% theta[idx$p]),
      as.numeric(X_omega  %*% theta[idx$omega]),
      as.numeric(X_gamma  %*% theta[idx$gamma]),
      use_nb = is_nb, eta_logr = if (is_nb) theta[ir] else 0.0)
  }
  grad_design <- function(out) {
    g <- numeric(n_par)
    g[idx$lambda] <- as.numeric(crossprod(X_lambda, out$grad_eta_lambda))
    g[idx$p]      <- as.numeric(crossprod(X_p,      out$grad_eta_p))
    g[idx$omega]  <- as.numeric(crossprod(X_omega,  out$grad_eta_omega))
    g[idx$gamma]  <- as.numeric(crossprod(X_gamma,  out$grad_eta_gamma))
    if (is_nb) g[ir] <- as.numeric(out$grad_eta_logr)
    g
  }
  neg_ll   <- function(theta) -eval_cpp(theta)$log_lik
  neg_grad <- function(theta) -grad_design(eval_cpp(theta))

  # Warm start: naive initial abundance from the first-season max count, p ~ 0.5,
  # omega ~ 0.6, gamma ~ 0.5. NB dispersion log r warm-started at log(2) (mild
  # overdispersion) unless overridden.
  J0 <- J; first_counts <- numeric(n_sites)
  for (i in seq_len(n_sites)) {
    seg <- y_flat[((i - 1L) * T * J0) + seq_len(J0)]      # season 1 visits
    seg <- seg[seg >= 0]
    first_counts[i] <- if (length(seg)) max(seg) else 0
  }
  theta0 <- numeric(n_par)
  theta0[idx$lambda[1]] <- log(max(mean(first_counts) / 0.5, 0.5) + 0.5)
  theta0[idx$p[1]]      <- 0
  theta0[idx$omega[1]]  <- stats::qlogis(0.6)
  theta0[idx$gamma[1]]  <- log(0.5)
  if (is_nb) theta0[ir] <- if (is.null(log_r_init)) log(2) else as.numeric(log_r_init)

  # Progress + ETA (gcol33/tulpaObs#43); ON by default. BFGS calls the gradient
  # ~once per quasi-Newton iteration, so ticking there approximates iteration
  # progress (maxit is the ETA denominator); finalised after optim returns.
  .prog <- tulpa:::.tulpa_iter_progress("dyn-abun-laplace", as.integer(max_iter), unit = "iter")
  neg_grad_p <- function(theta) { .prog$tick(); neg_grad(theta) }
  opt <- stats::optim(theta0, neg_ll, neg_grad_p, method = "BFGS",
                      control = list(maxit = as.integer(max_iter), reltol = tol))
  .prog$finish()
  theta <- opt$par
  converged <- opt$convergence == 0L
  out <- eval_cpp(theta)
  Hobs <- -.fp_fd_jacobian(function(th) grad_design(eval_cpp(th)), theta)
  vcov <- tryCatch(solve(Hobs), error = function(e)
    matrix(NA_real_, length(theta), length(theta)))

  nm <- c(paste0("lambda_", colnames(X_lambda)), paste0("p_", colnames(X_p)),
          paste0("omega_", colnames(X_omega)), paste0("gamma_", colnames(X_gamma)))
  if (is_nb) nm <- c(nm, "log_r")
  dimnames(vcov) <- list(nm, nm); dimnames(Hobs) <- list(nm, nm)
  if (!converged) {
    warning(sprintf("dyn_abun_laplace BFGS did not converge (code %d).", opt$convergence),
            call. = FALSE)
  }

  log_r <- if (is_nb) as.numeric(theta[ir]) else NA_real_
  structure(list(
    beta_lambda = theta[idx$lambda], beta_p = theta[idx$p],
    beta_omega = theta[idx$omega], beta_gamma = theta[idx$gamma],
    log_r = log_r, r = if (is_nb) exp(log_r) else NA_real_, mixture = mixture,
    means = theta, vcov = vcov, H_obs = Hobs,
    log_lik = out$log_lik, log_lik_site = out$log_lik_site, mean_N1 = out$mean_N1,
    converged = converged, n_iter = opt$counts[[1]], K_max = K_max,
    n_inadmissible = out$n_inadmissible, coef_names = nm),
    class = c("dyn_abun_fit", "list"))
}


# ---------------------------------------------------------------------------
# Fit packer
# ---------------------------------------------------------------------------

build_dyn_abun_fit <- function(raw, model, re_post = NULL) {
  pi_list <- model$process_info
  mixture <- raw$mixture %||% model$mixture %||% "poisson"
  is_nb   <- identical(mixture, "negbin")
  nms <- raw$coef_names
  means <- raw$means; names(means) <- nms
  vcov <- as.matrix(raw$vcov); dimnames(vcov) <- list(nms, nms)
  sds <- sqrt(pmax(diag(vcov), 0)); names(sds) <- nms
  n_fixed <- length(nms); fixed_names <- nms
  n_pseudo <- 1000L
  draws <- .rmvn(n_pseudo, means, vcov); colnames(draws) <- nms
  ll <- raw$log_lik %||% NA_real_

  # Grouped random effect on the initial-abundance (gcol33/tulpaObs#51) or the
  # detection arm (gcol33/tulpaObs#82): append the variance components (sigma_g_*,
  # cor_g_*_* for a correlated block) and the per-group BLUPs after the fixed
  # block, exactly as the other count families do. The fixed block (n_fixed
  # leading coords) still governs coef() / vcov() / confint(); ranef() / summary()
  # read the trailing RE columns by name.
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

  # NB dispersion summary on the natural (r) scale: r = exp(log_r), SE by the
  # delta method r * SE(log_r). Mirrors abun()'s nmix dispersion slot.
  dispersion <- NULL
  if (is_nb && is.finite(raw$log_r %||% NA_real_)) {
    se_logr <- if ("log_r" %in% nms) sqrt(pmax(vcov["log_r", "log_r"], 0)) else NA_real_
    dispersion <- list(r = as.numeric(raw$r %||% exp(raw$log_r)),
                       log_r = as.numeric(raw$log_r),
                       r_sd = as.numeric(exp(raw$log_r) * se_logr))
  }

  # Fixed-effect block excludes the dispersion coordinate from being treated as a
  # process coefficient; coef()/confint() still surface log_r as the trailing
  # fixed coordinate (it has a name and an SE), matching abun().
  structure(c(list(
    draws = draws, means = means, sds = sds, vcov = vcov,
    n_samples = n_pseudo, n_params = length(means),
    log_prob = rep(ll, n_pseudo), N = sum(model$y >= 0, na.rm = TRUE)),
    .tobs_na_nuts_diagnostics(n_pseudo),
    list(
    col_names = nms, param_names = nms,
    n_fixed = n_fixed, fixed_names = fixed_names,
    process_info = pi_list, model = model, spatial = NULL, method = "laplace",
    log_lik = ll, mean_N1 = raw$mean_N1, K_max = raw$K_max,
    mixture = mixture, dispersion = dispersion,
    re_effects = re_block$re_effects,
    dyn_abun_re = if (!is.null(re_post))
      list(arm = re_post$arm, n_quad = re_post$n_quad,
           lkj_eta = re_post$lkj_eta, Sigma_list = re_post$Sigma_list)
      else NULL,
    convergence = list(converged = raw$converged %||% TRUE,
                       n_iter = raw$n_iter %||% NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "dyn_abun")
# ---------------------------------------------------------------------------

# Reshape a survival / recruitment arm eta vector into a per-site [n_sites x
# (T-1)] interval matrix: a long-form vector (season-varying) is row-major (site,
# interval); a per-site vector is recycled across the site's intervals.
.tobs_dyn_abun_arm_matrix <- function(eta, n_sites, n_intervals) {
  if (length(eta) == n_sites * n_intervals) {
    matrix(eta, n_sites, n_intervals, byrow = TRUE)
  } else {
    matrix(eta, n_sites, n_intervals)
  }
}

# Per-site arm values and the marginal expected-abundance trajectory
# E[N_t] = omega_{t-1} E[N_{t-1}] + gamma_{t-1}, E[N_1] = lambda. omega / gamma
# are returned as per-site [n_sites x (T-1)] interval matrices (a single column
# under constant rates).
.tobs_fitted_dyn_abun <- function(object) {
  model <- object$model; means <- object$means
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  T <- model$n_seasons; n_sites <- model$n_sites; nIv <- T - 1L
  lambda <- exp(as.vector(model$X_processes[[1]] %*% means[off[1] + seq_len(p[1])]))
  pdet   <- stats::plogis(as.vector(model$X_processes[[2]] %*% means[off[2] + seq_len(p[2])]))
  omega  <- .tobs_dyn_abun_arm_matrix(
    stats::plogis(as.vector(model$X_processes[[3]] %*% means[off[3] + seq_len(p[3])])),
    n_sites, nIv)
  gamma  <- .tobs_dyn_abun_arm_matrix(
    exp(as.vector(model$X_processes[[4]] %*% means[off[4] + seq_len(p[4])])),
    n_sites, nIv)
  EN <- matrix(0, n_sites, T)
  EN[, 1] <- lambda
  for (t in seq_len(nIv)) EN[, t + 1L] <- omega[, t] * EN[, t] + gamma[, t]
  list(lambda = lambda, p = pdet, omega = omega, gamma = gamma, EN = EN)
}

# simulate() for dyn_abun: draw the abundance trajectory and emit counts.
.tobs_simulate_dyn_abun <- function(object, nsim = 1) {
  model <- object$model; draws <- object$draws; n_draws <- nrow(draws)
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  n_sites <- model$n_sites; T <- model$n_seasons; J <- model$max_visits
  is_nb <- identical(object$mixture %||% "poisson", "negbin")
  r_disp <- object$dispersion$r %||% NA_real_
  nIv <- T - 1L
  result <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    di <- sample.int(n_draws, 1L); th <- draws[di, ]
    lambda <- exp(as.vector(model$X_processes[[1]] %*% th[off[1] + seq_len(p[1])]))
    pdet   <- stats::plogis(as.vector(model$X_processes[[2]] %*% th[off[2] + seq_len(p[2])]))
    omega  <- .tobs_dyn_abun_arm_matrix(
      stats::plogis(as.vector(model$X_processes[[3]] %*% th[off[3] + seq_len(p[3])])),
      n_sites, nIv)
    gamma  <- .tobs_dyn_abun_arm_matrix(
      exp(as.vector(model$X_processes[[4]] %*% th[off[4] + seq_len(p[4])])),
      n_sites, nIv)
    ya <- array(0L, dim = c(n_sites, J, T))
    for (i in seq_len(n_sites)) {
      N <- if (is_nb && is.finite(r_disp)) stats::rnbinom(1L, mu = lambda[i], size = r_disp)
           else stats::rpois(1L, lambda[i])
      for (t in seq_len(T)) {
        if (t > 1L) N <- stats::rbinom(1L, N, omega[i, t - 1L]) + stats::rpois(1L, gamma[i, t - 1L])
        ya[i, , t] <- stats::rbinom(J, N, pdet[i])
      }
    }
    result[[s]] <- ya
  }
  if (nsim == 1L) result[[1]] else result
}

# residuals() for dyn_abun, on the marginal expected count mu_itj = E[N_t] * p.
.tobs_residuals_dyn_abun <- function(object, type = c("deviance", "pearson",
                                                    "response")) {
  type  <- match.arg(type)
  fitv  <- .tobs_fitted_dyn_abun(object)
  model <- object$model
  T <- model$n_seasons; J <- model$max_visits
  mu_t <- fitv$EN * fitv$p                       # [n_sites x T]
  y <- model$y                                   # [n_sites x J x T]
  out <- array(NA_real_, dim = dim(y))
  for (t in seq_len(T)) {
    mu <- pmax(matrix(mu_t[, t], model$n_sites, J), 1e-10)
    yt <- y[, , t]
    out[, , t] <- switch(type,
      response = yt - mu,
      pearson  = (yt - mu) / sqrt(mu),
      deviance = {
        d <- 2 * (ifelse(yt > 0, yt * log(yt / mu), 0) - (yt - mu))
        sign(yt - mu) * sqrt(pmax(d, 0))
      })
  }
  out
}

# predict() for dyn_abun: initial abundance lambda at new X (default).
.tobs_predict_dyn_abun <- function(object, X.0 = NULL, type = c("lambda", "gamma")) {
  type  <- match.arg(type)
  model <- object$model
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  k <- if (identical(type, "lambda")) 1L else 4L
  X <- X.0 %||% model$X_processes[[k]]
  exp(as.vector(X %*% object$means[off[k] + seq_len(p[k])]))
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate Dail-Madsen open N-mixture data
#'
#' Latent `N_1 ~ Poisson(lambda)`; for `t >= 2`, `N_t = Binomial(N_(t-1),
#' omega_{t-1}) + Poisson(gamma_{t-1})`; observed via `Binomial(N_t, p)` over `J`
#' secondary visits in each of `T` primary seasons. Returns a 3D array `N x J x T`
#' suitable for [tobs()] with [dyn_abun()].
#'
#' Survival and recruitment are constant across seasons by default. Supplying
#' `beta_omega` (logit-scale) or `beta_gamma` (log-scale) instead makes the rate
#' depend on a per-(site, interval) season covariate `z`: `omega_{i,t} =
#' plogis(beta_omega[1] + beta_omega[2] * z_{i,t})`, `gamma_{i,t} =
#' exp(beta_gamma[1] + beta_gamma[2] * z_{i,t})`. The covariate is returned as a
#' `[N x (T-1)]` matrix column `season_cov` of `data`, ready for
#' `omega_formula = ~ season_cov`.
#'
#' @param N Number of sites (default 150).
#' @param T Number of primary seasons (default 4).
#' @param J Number of secondary visits per season (default 3).
#' @param n_abund_covs Number of initial-abundance covariates (default 1).
#' @param beta_lambda Initial-abundance coefficients (log). Default
#'   `c(log(5), runif(n_abund_covs, -0.4, 0.4))`.
#' @param p,omega,gamma Detection, apparent-survival, and recruitment parameters
#'   (scalars; defaults 0.5, 0.6, 1.0). `omega` / `gamma` set the constant rate
#'   used when `beta_omega` / `beta_gamma` are `NULL`.
#' @param beta_omega,beta_gamma Optional length-2 coefficients
#'   `c(intercept, slope)` for a season-varying survival (logit link) or
#'   recruitment (log link) rate driven by a per-(site, interval) season
#'   covariate. `NULL` (default) keeps the constant `omega` / `gamma`.
#' @param mixture Initial-abundance distribution: `"poisson"` (default) or
#'   `"negbin"` (negative-binomial `N_1 ~ NB(mean = lambda, size = r)`).
#' @param r Negative-binomial size for `mixture = "negbin"` (default 2).
#' @param seed Optional random seed.
#' @return A list with `y` (N x J x T count array), `data` (covariates, including
#'   the `[N x (T-1)]` `season_cov` matrix column when a season-varying rate is
#'   used), and `truth` (coefficients, per-site `lambda`, the realised
#'   `omega` / `gamma` rates, and `r` under `"negbin"`).
#' @export
simulate_dyn_abun <- function(N = 150, T = 4, J = 3, n_abund_covs = 1,
                              beta_lambda = NULL, p = 0.5, omega = 0.6,
                              gamma = 1.0, beta_omega = NULL, beta_gamma = NULL,
                              mixture = c("poisson", "negbin"),
                              r = 2, seed = NULL) {
  mixture <- match.arg(mixture)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_lambda)) beta_lambda <- c(log(5), stats::runif(n_abund_covs, -0.4, 0.4))
  abund_covs <- data.frame(matrix(stats::rnorm(N * n_abund_covs), N, n_abund_covs))
  names(abund_covs) <- paste0("abund_cov", seq_len(n_abund_covs))
  X_lambda <- stats::model.matrix(~ ., abund_covs)
  lambda <- exp(as.vector(X_lambda %*% beta_lambda))
  is_nb <- identical(mixture, "negbin")
  nIv <- T - 1L

  # Per-(site, interval) survival / recruitment. Constant rates broadcast the
  # scalar; a supplied beta_omega / beta_gamma drives the rate off a shared
  # season covariate carried as a [N x (T-1)] matrix column of `data`.
  season_cov <- NULL
  if (!is.null(beta_omega) || !is.null(beta_gamma)) {
    season_cov <- matrix(stats::rnorm(N * nIv), N, nIv)
  }
  omega_mat <- if (!is.null(beta_omega))
                 stats::plogis(beta_omega[1] + beta_omega[2] * season_cov)
               else matrix(omega, N, nIv)
  gamma_mat <- if (!is.null(beta_gamma))
                 exp(beta_gamma[1] + beta_gamma[2] * season_cov)
               else matrix(gamma, N, nIv)

  y <- array(0L, dim = c(N, J, T))
  Nmat <- matrix(0L, N, T)
  for (i in seq_len(N)) {
    Ni <- if (is_nb) stats::rnbinom(1L, mu = lambda[i], size = r)
          else stats::rpois(1L, lambda[i])
    for (t in seq_len(T)) {
      if (t > 1L) Ni <- stats::rbinom(1L, Ni, omega_mat[i, t - 1L]) +
                        stats::rpois(1L, gamma_mat[i, t - 1L])
      Nmat[i, t] <- Ni
      y[i, , t] <- stats::rbinom(J, Ni, p)
    }
  }
  data <- abund_covs
  if (!is.null(season_cov)) data$season_cov <- season_cov
  list(y = y, data = data,
       truth = list(beta_lambda = beta_lambda, lambda = lambda, p = p,
                    omega = omega, gamma = gamma,
                    beta_omega = beta_omega, beta_gamma = beta_gamma,
                    omega_mat = omega_mat, gamma_mat = gamma_mat,
                    season_cov = season_cov, mixture = mixture,
                    r = if (is_nb) r else NA_real_, N = Nmat))
}
