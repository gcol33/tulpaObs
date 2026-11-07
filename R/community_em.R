# Generic family-agnostic community Laplace-EM engine.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# The `mu_` prefix on a community-mean coefficient name (fit-level `means` /
# `vcov` / `sds` / draws), distinguishing it from the bare design-column name
# used on the per-species blocks (`coef_mu` / `blup_mu` / `Sigma_mu`). One
# definition for every community output that reports a mean-coefficient arm
# (EM, NUTS, PG-Gibbs) so the three cannot drift apart.
.tobs_community_mean_names <- function(cn) paste0("mu_", cn)
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

# Per-species BLUP deviations in the ranef() layout: one row per (species, arm,
# term), `group = "species"`, `level` the species name. `arms` is a named
# character vector mapping the arm label (name) to the `ms_community` blup
# field (value); fields absent from `cm` are skipped (e.g. the optional
# per-species log_r RE of an NB N-mixture). `std.error` is the per-species
# posterior SD from `.tobs_cem_blup_sd()`, NA for a block the fit carries no
# covariance for. Shared ranef() body for every community family (ms_occu /
# ms_dyn_occu / ms_int_occu / ms_occu_cover / ms_nmix / ...).
.tobs_ranef_ms_long <- function(cm, arms) {
  blocks <- Map(function(field, arm) {
    B <- cm[[field]]
    if (is.null(B)) return(NULL)
    sp <- rownames(B) %||% as.character(seq_len(nrow(B)))
    tm <- colnames(B) %||% paste0("coef", seq_len(ncol(B)))
    SD <- .tobs_cem_blup_sd(cm, field)
    if (!is.null(SD) && !identical(dim(SD), dim(as.matrix(B)))) {
      stop(sprintf(paste0(
        "ms_community$%s is %d x %d but its posterior SD block is %d x %d; ",
        "the index map `blup_idx` does not match the BLUP layout."),
        field, nrow(B), ncol(B), nrow(SD), ncol(SD)), call. = FALSE)
    }
    .tobs_ranef_table(estimate = as.numeric(B),
                      level = rep(sp, times = ncol(B)),
                      term  = rep(tm, each = nrow(B)),
                      arm = arm, group = "species",
                      std.error = if (is.null(SD)) NA_real_ else as.numeric(SD))
  }, arms, names(arms))
  .tobs_ranef_stack(blocks)
}

# Per-species posterior SD of one community BLUP block. `cm$Cinv[[s]]` is
# Cov(b_s | y), the posterior covariance of species s's stacked RE vector, and
# `cm$blup_idx[[field]]` the positions of the `field` block's columns
# (`blup_psi`, `blup_p`, ...) in that vector, recorded by the builder from the
# same index map it slices `b_s` with. The SD of row s, column j of
# `cm[[field]]` is sqrt(Cinv[[s]][idx[j], idx[j]]). The covariance is
# conditional on the community means on the Laplace routes and marginal over
# them on the sampler routes. Returns an S x length(idx) matrix, or NULL when
# the fit records no covariance; a covariance without an index map for the
# block is an error, so a route that has the SD cannot report NA. A species
# with no covariance (a failed per-species solve) gets NA.
.tobs_cem_blup_sd <- function(cm, field) {
  idx <- cm$blup_idx[[field]]
  C   <- cm$Cinv
  if (!length(C)) return(NULL)
  if (is.null(idx)) {
    stop(sprintf(paste0(
      "ms_community carries a per-species covariance but no index map for ",
      "`%s` (ms_community$blup_idx$%s)."), field, field), call. = FALSE)
  }
  rows <- lapply(C, function(Cs) {
    if (is.null(Cs)) return(rep(NA_real_, length(idx)))
    Cs <- as.matrix(Cs)
    if (max(idx) > nrow(Cs)) {
      stop(sprintf(paste0(
        "ms_community$blup_idx$%s addresses position %d of a %d x %d ",
        "per-species covariance."), field, max(idx), nrow(Cs), ncol(Cs)),
        call. = FALSE)
    }
    sqrt(pmax(diag(Cs)[idx], 0))
  })
  matrix(unlist(rows, use.names = FALSE), nrow = length(C), byrow = TRUE)
}

# Running first and second moments of the per-species deviation matrix over
# sampler draws, for the sampler routes' per-species posterior covariance.
# `add` folds one draw's S x P deviation matrix in, `merge` pools the
# accumulators of several chains, `mean` is the posterior-mean deviation and
# `cov` the per-species sample covariance of the stacked deviation vector
# (n - 1 denominator, as stats::cov), NULL below two draws.
.tobs_cem_moments_new <- function(S, P) {
  list(n = 0L, sum = matrix(0, S, P), sq = rep(list(matrix(0, P, P)), S))
}

.tobs_cem_moments_add <- function(acc, B) {
  acc$n   <- acc$n + 1L
  acc$sum <- acc$sum + B
  for (s in seq_len(nrow(B))) acc$sq[[s]] <- acc$sq[[s]] + tcrossprod(B[s, ])
  acc
}

.tobs_cem_moments_merge <- function(accs) {
  Reduce(function(a, b) {
    a$n <- a$n + b$n; a$sum <- a$sum + b$sum; a$sq <- Map(`+`, a$sq, b$sq); a
  }, accs)
}

.tobs_cem_moments_mean <- function(acc) acc$sum / acc$n

.tobs_cem_moments_cov <- function(acc) {
  n <- acc$n
  if (n < 2L) return(NULL)
  M <- acc$sum / n
  lapply(seq_len(nrow(M)), function(s) {
    C <- (acc$sq[[s]] - n * tcrossprod(M[s, ])) / (n - 1)
    (C + t(C)) / 2
  })
}


# ---------------------------------------------------------------------------
# Community Laplace-EM engine
# ---------------------------------------------------------------------------

#' Generic community Laplace-EM engine
#'
#' A community model's per-species coefficients are `theta_s = mu + b_s`, with
#' `b_s ~ N(0, Sigma)` block-diagonal over the arms in `arm_idx` and the
#' community means `mu` (and any `global` parameters) shared. The fit is
#' tulpa's Laplace-EM on the per-species oracle (`tulpa::tulpa_re_em()`,
#' `theta_update = "latent"`): `(mu, global)` and every `b_s` solved jointly by
#' a damped Newton on the penalized joint log-density, the closed-form
#' covariance update per arm, and the community-mean covariance from the Schur
#' complement of the species blocks. This function supplies the oracle: species
#' `s`'s log-likelihood `sp_ll(s, theta, global)` with its derivatives.
#'
#' @keywords internal
#' @param sp_grad Optional `function(s, theta, global)` returning the gradient
#'   of `sp_ll` in `c(theta, global)`; `NULL` central-differences `sp_ll`.
#' @param sp_info Optional `function(s, theta, global)` returning species `s`'s
#'   `(P + G) x (P + G)` marginal observed information at `c(theta, global)`.
#'   `NULL` central-differences the gradient, which costs `2 (P + G)` gradient
#'   evaluations per species per Newton step; supply it when the family's
#'   kernel already exposes the per-site marginal observed information (the
#'   N-mixture / distance Louis block).
#' @param init_b,init_Sigma Optional warm starts for the per-species deviations
#'   and the per-arm community covariances; `NULL` is the cold start
#'   (`b_s = 0`, `Sigma = sigma_init^2 I`). A block-coordinate caller that
#'   re-enters this EM once per outer pass (R/community_latent.R) passes the
#'   previous pass's state, so the EM resumes from a near-converged point.
#' @param re_aghq Debias the community covariances by one adaptive
#'   Gauss-Hermite moment step at the converged mode (tulpa_re_em()'s
#'   `aghq_debias`), gated to `P <= re_aghq_maxdim` (the tensor grid spans the
#'   species' whole coefficient vector, coupling the arms). `n_quad` sets the
#'   per-axis node count. The per-species covariances and the community-mean
#'   covariance are recomputed at the debiased covariance, so `Sigma`, `Cinv`,
#'   `Bf` and `Vf` describe one fit (the joint (u, b_s) draw of
#'   .tobs_sbc_community_b_draws() assumes it).
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
  arm_names  <- names(arm_idx)
  glob_seq   <- seq_len(G)
  # The oracle's RE vector stacks the arms in `arm_idx` order; `perm` maps its
  # coordinates back to the species coefficient vector.
  perm <- as.integer(unlist(arm_idx, use.names = FALSE))
  if (!identical(sort(perm), seq_len(P))) {
    stop(".tobs_community_em: `arm_idx` must partition 1:P.", call. = FALSE)
  }
  to_coef <- function(b) { o <- numeric(P); o[perm] <- b; o }

  # Prior precision on (mu, global): the mu-block when priors are active, the
  # global block only where penalize_global is TRUE (recycled).
  priors_active <- isTRUE(is.null(priors)) || !isFALSE(priors)
  prec <- numeric(U)
  if (priors_active) {
    prec[seq_len(P)] <- 1 / sigma_beta^2
    if (G > 0L) {
      pen_g <- rep(as.logical(penalize_global), length.out = G)
      prec[P + glob_seq][pen_g] <- 1 / sigma_beta^2
    }
  }

  # Per-species gradient of the marginal log-likelihood wrt c(theta, global)
  # (length P + G): supplied, or central difference of sp_ll (h = 1e-5).
  sp_grad_fn <- if (!is.null(sp_grad)) sp_grad else {
    function(s, theta, global) {
      u <- c(theta, global); h <- 1e-5
      g <- numeric(U)
      for (k in seq_len(U)) {
        up <- u; up[k] <- up[k] + h
        dn <- u; dn[k] <- dn[k] - h
        g[k] <- (sp_ll(s, up[seq_len(P)], up[P + glob_seq]) -
                 sp_ll(s, dn[seq_len(P)], dn[P + glob_seq])) / (2 * h)
      }
      g
    }
  }
  # Per-species observed information (U x U): supplied, or central difference
  # of the gradient (h = 1e-4), symmetrized.
  sp_info_fn <- if (!is.null(sp_info)) sp_info else {
    function(s, theta, global) {
      u <- c(theta, global); h <- 1e-4
      J <- matrix(0, U, U)
      for (k in seq_len(U)) {
        up <- u; up[k] <- up[k] + h
        dn <- u; dn[k] <- dn[k] - h
        J[, k] <- (sp_grad_fn(s, up[seq_len(P)], up[P + glob_seq]) -
                   sp_grad_fn(s, dn[seq_len(P)], dn[P + glob_seq])) / (2 * h)
      }
      -0.5 * (J + t(J))
    }
  }

  # The per-species oracle at theta = c(mu, global): ell_s(b) = sp_ll(s, mu + b,
  # global). theta enters as a location shift of b, so the theta score is the
  # species gradient itself and the theta-b cross block its information read at
  # the b columns.
  make_group <- function(theta) {
    mu <- theta[seq_len(P)]; gl <- theta[P + glob_seq]
    joint <- function(s, b) {
      th  <- mu + to_coef(b)
      gr  <- sp_grad_fn(s, th, gl)
      inf <- sp_info_fn(s, th, gl)
      list(logL = sp_ll(s, th, gl), grad = gr[perm],
           negH = inf[perm, perm, drop = FALSE],
           theta_grad = gr, theta_info = inf,
           cross = inf[, perm, drop = FALSE])
    }
    list(grad_hess = function(s, b) joint(s, b)[c("logL", "grad", "negH")],
         node_ll = function(s, B) apply(B, 1L, function(b)
           sp_ll(s, mu + to_coef(b), gl)),
         joint = joint)
  }

  Sigma0 <- if (is.null(init_Sigma)) {
    lapply(arm_idx, function(idx) diag(sigma_init^2, length(idx)))
  } else init_Sigma[arm_names]
  b0 <- if (is.null(init_b)) NULL else
    do.call(rbind, lapply(init_b, function(b) as.numeric(b)[perm]))
  re_terms <- lapply(arm_idx, function(idx)
    list(n_groups = S, n_coefs = length(idx), correlated = TRUE))

  debias <- isTRUE(re_aghq) && P <= as.integer(re_aghq_maxdim)
  em <- tulpa::tulpa_re_em(
    theta0 = c(as.numeric(init_mu), as.numeric(init_global)),
    re_terms = unname(re_terms), Sigma0 = unname(lapply(Sigma0, as.matrix)),
    make_group = make_group,
    theta_prior_sd = ifelse(prec > 0, 1 / sqrt(prec), Inf),
    theta_update = "latent", b0 = b0,
    # Eigenvalue floor 1e-4 on each community covariance.
    sigma_floor = 1e-2, tol = as.numeric(tol),
    max_iter = as.integer(max_iter), newton_max = as.integer(newton_max),
    progress_label = "community-em",
    aghq_debias = if (debias) as.integer(n_quad) else 0L)
  if (is.null(em)) {
    stop("community EM: the joint log-density is not finite at the starting ",
         "values.", call. = FALSE)
  }
  converged <- isTRUE(em$converged)
  if (isTRUE(em$stuck)) {
    converged <- FALSE
    warning(paste0(
      "community EM: the Newton mode-find never accepted a step (undamped or ",
      "Levenberg-damped); the coefficients did not move from their starting ",
      "values. Reporting converged = FALSE rather than certifying the fit off ",
      "the M-step's objective drift alone."), call. = FALSE)
  }

  theta  <- em$theta
  mu     <- theta[seq_len(P)]
  global <- theta[P + glob_seq]
  B      <- do.call(cbind, em$blup)
  b_list <- lapply(seq_len(S), function(s) to_coef(B[s, ]))
  Cinv   <- lapply(em$blup_cov_g, function(C) {
    o <- matrix(0, P, P); o[perm, perm] <- C; o
  })
  # Mode/theta cross blocks -d2 ell_s / d(mu, global) db_s, (P + G) x P.
  Bf <- lapply(em$blup_cross_g, function(X) {
    o <- matrix(0, U, P); o[, perm] <- X; o
  })
  Sigma <- stats::setNames(lapply(em$Sigma_list, as.matrix), arm_names)
  Sf <- em$theta_hessian
  Vf <- tryCatch(solve(Sf), error = function(e) .tobs_cem_ginv(Sf))
  Vf <- (Vf + t(Vf)) / 2
  logML <- em$log_marginal - 0.5 * sum(prec * theta^2)
  if (debias && !isTRUE(em$aghq_debiased)) {
    warning("community EM AGHQ variance debias failed (a per-species ",
            "factorization or quadrature sum); reporting the EM (Laplace) ",
            "variance components.", call. = FALSE)
  }
  debias_method <- if (isTRUE(em$aghq_debiased)) "aghq" else "none"

  list(mu = mu, global = global, b_list = b_list, Sigma = Sigma,
       Cinv = Cinv, Bf = Bf, Vf = Vf, logML = logML,
       converged = converged, n_iter = em$n_iter, debias_method = debias_method)
}


# ---------------------------------------------------------------------------
# Shared tobs_fit skeleton
# ---------------------------------------------------------------------------

# Assemble the common `tobs_fit` skeleton for a Laplace-EM community fit
# (`.tobs_community_em()` or a block-coordinate spatial-factor variant built on
# top of it) -- the PG-Gibbs analogue of `.tobs_pg_finalize_fit()`
# (R/pg_gibbs_shared.R). `means`/`V` are the raw (unnamed) community-mean point
# estimate and covariance; naming, symmetrizing, `sds`, and the posterior draws
# are done here so every caller does it once. `extra` carries the
# family-specific tail (`ms_community`, `ms_dispersion`, `positive`, ...);
# `spatial`/`method` override the single-species-Laplace defaults.
# `grid_mixture` (`.tobs_grid_mixture()`) is the outer grid of a grid-integrated
# fit, which the draws then come from instead of the Gaussian at (means, V).
.tobs_cem_finalize_fit <- function(means, V, par_names, model, process_info, N,
                                   log_prob_val, converged, n_iter,
                                   spatial = NULL, method = "laplace",
                                   grid_mixture = NULL, extra = list()) {
  means <- as.numeric(means); names(means) <- par_names
  V <- (V + t(V)) / 2
  dimnames(V) <- list(par_names, par_names)
  sds <- .tobs_sds_from_vcov(V, par_names)

  n_draws <- 1000L
  draws <- .tobs_grid_mixture_draws(n_draws, grid_mixture$weights,
                                    grid_mixture$modes, grid_mixture$covs,
                                    means, V)
  colnames(draws) <- par_names

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(log_prob_val, n_draws),
    log_lik      = log_prob_val,
    N            = N),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    n_fixed      = length(means),
    fixed_names  = par_names,
    process_info = process_info,
    model        = model,
    spatial      = spatial,
    method       = method),
    extra,
    list(convergence = list(converged = isTRUE(converged), n_iter = n_iter))
  ), class = c("tobs_fit", "tulpa_fit"))
}
