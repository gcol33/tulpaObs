# =============================================================================
# Penalized binomial IRLS for occupancy / detection M-step fits.
#
# tulpa_laplace() does not currently support a coefficient prior, and the
# psi-p ridge in single-season occupancy at small J needs a weakly-informative
# prior to break the identifiability disaster (D1 sweep, 2026-05-15: Laplace
# bias on `psi_(Intercept)` was ~ -1.5 across 30 seeds at N=600, J=6 even
# after the SE fix in tulpaObs#1 / commit 2699daa).
#
# Math. For each M-step submodel block we minimise the penalised binomial
# negative log-posterior
#
#       Q(beta) = - sum_i w_i [ y_i log mu_i + (n_i - y_i) log(1 - mu_i) ]
#                 + sum_j (beta_j - mu_prior_j)^2 / (2 * sd_prior_j^2)
#
# where mu_i = plogis(X_i beta), w_i is an observation weight, n_i a trial
# count. The gradient is
#
#       grad Q = - X^T (w * (y - n * mu))  +  P * (beta - mu_prior)
#
# and the Hessian
#
#       H Q = X^T diag(w * n * mu * (1 - mu)) X  +  P,
#
# with P = diag(1 / sd_prior_j^2) a positive-definite penalty matrix. The
# penalty contribution to H is positive-definite, so it both regularises
# the MAP toward the prior mode and tightens / well-conditions the Laplace
# precision used for SEs.
#
# Setting sd_prior_j = Inf makes 1 / sd_prior_j^2 = 0 and recovers the
# unpenalised MAP (no fallback path: this is the same single objective).
#
# This module covers the non-spatial, no-extra-RE binomial blocks that
# every occupancy M-step uses; spatial blocks still route through
# tulpa_laplace() (priors-on-fixed-effects through the SPDE / NNGP path
# are not yet plumbed, tracked in tulpaObs#5).
# =============================================================================


# Normalise a prior block to a (mean vector, sd vector) of length p.
# `prior` is a list with $mean (scalar or length-p) and $sd (scalar or
# length-p), where sd may include +Inf entries (no penalty on that coef).
.expand_prior <- function(prior, p) {
  if (is.null(prior)) {
    return(list(mean = rep(0, p), sd = rep(Inf, p)))
  }
  m <- prior$mean
  s <- prior$sd
  if (is.null(m)) m <- 0
  if (is.null(s)) s <- Inf
  if (length(m) == 1L) m <- rep(m, p)
  if (length(s) == 1L) s <- rep(s, p)
  if (length(m) != p) {
    stop(sprintf(
      ".expand_prior(): prior$mean has length %d, expected %d", length(m), p),
      call. = FALSE)
  }
  if (length(s) != p) {
    stop(sprintf(
      ".expand_prior(): prior$sd has length %d, expected %d", length(s), p),
      call. = FALSE)
  }
  if (any(!is.finite(m))) {
    stop(".expand_prior(): prior$mean must be finite (got NA / +-Inf)",
         call. = FALSE)
  }
  if (any(s <= 0)) {
    stop(".expand_prior(): prior$sd must be positive (Inf is allowed for no penalty)",
         call. = FALSE)
  }
  list(mean = as.numeric(m), sd = as.numeric(s))
}


# Build a prior block for the occupancy / detection submodel given the
# user-facing tobs prior list and the coefficient names. The prior spec
# accepts four named buckets:
#   $p_intercept       — list(mean, sd) applied to the detection intercept
#   $p_slope           — list(mean, sd) applied to every non-intercept detection coef
#   $beta_occ_intercept — list(mean, sd) applied to the occupancy / psi intercept
#   $beta_occ_slope    — list(mean, sd) applied to every non-intercept occupancy coef
#
# The "intercept" column is detected by name `(Intercept)`, which is what
# model.matrix() uses; falls back to the first column otherwise.
.prior_for_submodel <- function(prior_spec, sub_name, coef_names) {
  p <- length(coef_names)
  if (p == 0L) return(.expand_prior(NULL, 0L))
  if (is.null(prior_spec)) return(.expand_prior(NULL, p))

  is_intercept <- coef_names == "(Intercept)"
  if (!any(is_intercept)) {
    is_intercept <- c(TRUE, rep(FALSE, p - 1L))[seq_len(p)]
  }

  if (sub_name %in% c("p")) {
    int_b <- prior_spec$p_intercept
    slp_b <- prior_spec$p_slope
  } else if (sub_name %in% c("psi", "psi1")) {
    int_b <- prior_spec$beta_occ_intercept
    slp_b <- prior_spec$beta_occ_slope
  } else if (sub_name %in% c("gamma")) {
    int_b <- prior_spec$gamma_intercept %||% prior_spec$beta_occ_intercept
    slp_b <- prior_spec$gamma_slope %||% prior_spec$beta_occ_slope
  } else if (sub_name %in% c("epsilon")) {
    int_b <- prior_spec$epsilon_intercept %||% prior_spec$beta_occ_intercept
    slp_b <- prior_spec$epsilon_slope %||% prior_spec$beta_occ_slope
  } else if (grepl("^det[0-9]+$", sub_name) || sub_name == "det") {
    int_b <- prior_spec$p_intercept
    slp_b <- prior_spec$p_slope
  } else {
    int_b <- prior_spec$p_intercept
    slp_b <- prior_spec$p_slope
  }

  mean_vec <- numeric(p)
  sd_vec   <- numeric(p)
  for (j in seq_len(p)) {
    if (is_intercept[j]) {
      m <- if (is.null(int_b)) 0     else int_b$mean %||% 0
      s <- if (is.null(int_b)) Inf   else int_b$sd   %||% Inf
    } else {
      m <- if (is.null(slp_b)) 0     else slp_b$mean %||% 0
      s <- if (is.null(slp_b)) Inf   else slp_b$sd   %||% Inf
    }
    mean_vec[j] <- m
    sd_vec[j]   <- s
  }
  list(mean = mean_vec, sd = sd_vec)
}


# `%||%` is defined in R/family_cover_hurdle.R; reused here.


# Penalised binomial IRLS for a single M-step block.
#
# Returns the same shape as tulpa_laplace() so build_laplace_fit() can read
# `mode`, `H_beta`, `se` without branching on the fitter.
#
# Per-row observation weight is forwarded from the block (`weights`). For
# binomial-with-trials this enters as the multiplicative scalar on the
# log-likelihood, *not* as a replicate count — matching the existing
# tulpa_laplace() semantics so M-step encoded blocks behave identically to
# the unpenalised path when prior_sd = Inf.
.penalized_irls_binomial <- function(y, n_trials, X, weights = NULL,
                                     prior_mean, prior_sd,
                                     beta_init = NULL,
                                     max_iter = 100L, tol = 1e-8) {
  n_obs   <- length(y)
  n_fixed <- ncol(X)
  if (is.null(weights)) weights <- rep(1.0, n_obs)
  if (is.null(n_trials)) n_trials <- rep(1L, n_obs)

  # Penalty precision diag(1/sd^2); +Inf sd -> 0 precision (no penalty).
  pen_prec <- ifelse(is.finite(prior_sd), 1 / (prior_sd^2), 0)

  if (is.null(beta_init)) {
    beta <- rep(0, n_fixed)
  } else {
    beta <- as.numeric(beta_init)
    if (length(beta) != n_fixed) {
      beta <- rep(0, n_fixed)
    }
  }

  converged <- FALSE
  n_iter <- 0L
  prev_q <- Inf

  for (iter in seq_len(max_iter)) {
    eta <- as.numeric(X %*% beta)
    eta <- pmin(pmax(eta, -30), 30)        # logit clamp
    mu  <- 1 / (1 + exp(-eta))
    mu  <- pmin(pmax(mu, 1e-10), 1 - 1e-10)

    # Working weights and residual for IRLS.
    w_diag <- weights * n_trials * mu * (1 - mu)
    grad_ll <- as.numeric(crossprod(X, weights * (y - n_trials * mu)))
    grad <- -grad_ll + pen_prec * (beta - prior_mean)

    # Hessian: X^T W X + diag(pen_prec).
    XtWX <- crossprod(X, w_diag * X)
    H <- as.matrix(XtWX)
    diag(H) <- diag(H) + pen_prec

    # Newton step. tryCatch protects against singular H in pathological
    # blocks (e.g. an all-zero design column with sd = Inf on that coef).
    step <- tryCatch(solve(H, grad), error = function(e) NULL)
    if (is.null(step)) {
      # Add tiny jitter and retry once; if still singular, return current
      # beta with NA-filled H so build_laplace_fit can NA-fill the SE.
      H_j <- H; diag(H_j) <- diag(H_j) + 1e-8
      step <- tryCatch(solve(H_j, grad), error = function(e) NULL)
      if (is.null(step)) {
        H <- matrix(NA_real_, n_fixed, n_fixed)
        break
      }
    }

    # Damped step for first few iterations to avoid overshoot on flat
    # likelihoods (large logit linear predictors).
    damp <- if (iter <= 3L) 0.5 else 1.0
    beta_new <- beta - damp * as.numeric(step)

    # Convergence: max abs param change relative to current scale.
    max_change <- max(abs(beta_new - beta) / pmax(abs(beta), 1))
    beta <- beta_new
    n_iter <- iter
    if (max_change < tol) {
      converged <- TRUE
      break
    }
  }

  # Final Hessian at the converged beta.
  eta <- pmin(pmax(as.numeric(X %*% beta), -30), 30)
  mu  <- 1 / (1 + exp(-eta))
  mu  <- pmin(pmax(mu, 1e-10), 1 - 1e-10)
  w_diag <- weights * n_trials * mu * (1 - mu)
  H_final <- as.matrix(crossprod(X, w_diag * X))
  diag(H_final) <- diag(H_final) + pen_prec

  se <- tryCatch({
    cov <- solve(H_final)
    sqrt(pmax(diag(cov), 0))
  }, error = function(e) rep(NA_real_, n_fixed))

  list(
    mode      = beta,
    H_beta    = H_final,
    se        = se,
    n_iter    = n_iter,
    converged = converged,
    prior_mean = prior_mean,
    prior_sd   = prior_sd
  )
}


# Fit one M-step block. If the block carries a `prior` element and is not
# spatial / has no extra RE, the penalised IRLS is used. Otherwise the
# block is forwarded to tulpa_laplace() unchanged (spatial / RE blocks
# don't yet have prior support in this package; see tulpaObs#5).
#
# `block$prior` should be a list with $mean, $sd (length = ncol(X)). When
# absent or NULL, this routes to the unpenalised path (same as before).
.fit_block_penalized <- function(block) {
  has_spatial <- !is.null(block$spatial)
  has_re      <- !is.null(block$re_list) && length(block$re_list) > 0
  has_prior   <- !is.null(block$prior)
  is_binomial <- identical(block$family, "binomial")

  if (has_prior && is_binomial && !has_spatial && !has_re) {
    p <- ncol(block$X)
    pr <- .expand_prior(block$prior, p)
    .penalized_irls_binomial(
      y          = as.numeric(block$y),
      n_trials   = block$n_trials,
      X          = block$X,
      weights    = block$weights,
      prior_mean = pr$mean,
      prior_sd   = pr$sd,
      beta_init  = block$beta_init
    )
  } else {
    # No prior, or block kind we don't yet penalise — forward to the
    # generic tulpa fitter. This is the same single objective: the prior
    # term in `Q` simply has all sd = +Inf, giving zero penalty.
    n_trials <- if (is.null(block$n_trials)) rep(1L, length(block$y)) else block$n_trials
    args <- list(
      y         = block$y,
      n_trials  = n_trials,
      X         = block$X,
      re_list   = if (is.null(block$re_list)) list() else block$re_list,
      family    = block$family,
      spatial   = block$spatial,
      weights   = block$weights,
      offset    = block$offset,
      n_threads = 1L,
      return_hessian = TRUE
    )
    if (!is.null(block$phi)) args$phi <- block$phi
    do.call(tulpa::tulpa_laplace, args)
  }
}
