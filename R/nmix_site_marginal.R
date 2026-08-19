# Latent-N truncation policy, shared by every N-mixture marginal.
#
# The per-site sum runs on [max(y_i), K]. The lower end is already the site's own
# maximum (Binom(y, N, p) is zero below it), so only the ceiling is shared. The
# states above max(y_i) carry the individuals never detected there, and that
# posterior decays geometrically, so truncation error is controlled by the
# HEADROOM above a site's own maximum, not by the absolute ceiling.
#
# A shared ceiling of max(y) + h gives the site holding the largest count exactly
# h states of headroom, and forces every other site out to the same absolute
# ceiling: one site with y = 2248 makes a site with y = 5 evaluate 2344 states
# whose combined posterior mass is ~0. Cost is linear in the state count, so the
# whole fit slows by that ratio.
#
# Default (`K_max = NULL`): the same global ceiling as before, plus a per-site
# cap at max(y_i) + h, which gives every site the headroom the binding site
# already had. An explicit `K_max` keeps its documented meaning -- a hard global
# truncation -- and is never capped, since a caller raising it is compensating
# for exactly the case (high abundance, low detection) the headroom cannot see.
.NMIX_HEADROOM <- 100L

.nmix_truncation <- function(K_max, y) {
  if (is.null(K_max)) {
    return(list(K_max = as.integer(max(y) + .NMIX_HEADROOM),
                headroom = .NMIX_HEADROOM))
  }
  list(K_max = as.integer(K_max), headroom = -1L)
}

# A per-site cap is exact wherever the posterior over N decays inside it, and no
# fixed headroom survives p -> 0: the count of never-detected individuals is
# ~Poisson(lambda (1-p)^J), so on the lambda/p ridge an N-mixture fit explores,
# the window a cap leaves can be crossed. Measured on a community fit, the
# max(y_i) + 100 cap starts leaving mass at its truncation once lambda sits ~5x
# above the data with detection near 0.1, while the shared ceiling holds far
# longer -- so the cap cannot simply be assumed and has to be VERIFIED.
#
# The criterion is NOT the boundary mass nmix_laplace() warns on. That bounds the
# truncation error at the point the fit stopped, and it is not enough: on a
# weakly identified fixture a capped fit passed the boundary check at 4.6e-06
# and still sat 0.032 nats below the uncapped optimum, 0.57 away in the
# coefficients -- the optimiser's PATH ran through the truncated region even
# though its endpoint did not.
#
# What has to hold is that the answer is also a stationary point of the UNCAPPED
# likelihood, so that is what gets tested: the score at the fitted coefficients
# under the capped and the shared-ceiling truncation, which must agree. One
# extra marginal evaluation against a fit that took tens of them. On the fixture
# above the two scores differ by 7.4e-02 where the boundary masses differed by
# nothing; on every fixture where the cap is harmless they agree exactly.
#
# Above the tolerance the fit is redone with a wider window, escalating to the
# uncapped shared ceiling, so a guarded fit is never worse than the
# shared-ceiling fit it replaces -- at worst it costs the extra fits.
.NMIX_BOUNDARY_TOL <- 1e-4

# Largest disagreement between the capped and uncapped score, at the answer,
# that still counts as the same optimum. Absolute, on the log-likelihood
# gradient in coefficient units.
.NMIX_SCORE_TOL <- 1e-4

# Next window to try after one was exhausted: widen 4x, and once widening would
# reach the shared ceiling, drop the cap entirely. NULL when nothing is left to
# widen (already uncapped), so a caller knows to stop rather than loop.
.nmix_widen_headroom <- function(headroom, K_max) {
  if (is.null(headroom) || is.na(headroom) || headroom < 0L) return(NULL)
  if (headroom * 4L >= K_max) return(-1L)
  as.integer(headroom * 4L)
}

# Largest disagreement between the capped and the uncapped score at the given
# per-species coefficients: the quantity that decides whether a capped community
# fit is also the uncapped one. Both scores come from the same per-site marginal
# at the same coefficients, differing only in the truncation, so the comparison
# isolates the truncation and nothing else.
#
# `eta_lambda_off` carries any additive abundance offset the fit ran with (a
# latent factor surface, a shared field): the truncation has to hold at the
# predictor the fit actually used, not at the fixed-effect part of it, and a
# latent surface can push a site's abundance well above what its own counts
# suggest -- the direction that exhausts a window.
.nmix_community_score_gap <- function(lf, X_lambda, coef_lambda, coef_p,
                                      K_max, headroom, eta_lambda_off = NULL) {
  worst <- 0
  for (s in seq_len(nrow(coef_lambda))) {
    k <- lf$species_idx == s
    if (!any(k)) next
    Xp <- lf$X_p[k, , drop = FALSE]
    el <- as.numeric(X_lambda %*% coef_lambda[s, ])
    if (!is.null(eta_lambda_off)) el <- el + eta_lambda_off[, s]
    ep <- as.numeric(Xp %*% coef_p[s, ])
    sc <- function(h) {
      ev <- nmix_site_marginal(lf$y[k], lf$site_idx[k], X_lambda, Xp,
                               mixture = "P", K_max = K_max,
                               headroom = h)$eval(el, ep)
      c(as.numeric(crossprod(X_lambda, ev$grad_eta_lambda)),
        as.numeric(crossprod(Xp, ev$grad_eta_p)))
    }
    d <- suppressWarnings(max(abs(sc(headroom) - sc(-1L))))
    if (is.finite(d) && d > worst) worst <- d
  }
  worst
}

#' Per-site N-mixture marginal as a composable random-effect callback
#'
#' @description
#' Exposes the Royle (2004) N-mixture per-site marginal -- the latent abundance
#' \eqn{N_i} summed out in closed form -- as a reusable building block for
#' integrating grouped random effects over the abundance and/or detection
#' linear predictors. It is the bridge between the marginal kernel (which
#' already differentiates through both arms) and a random-effect integrator
#' such as [tulpa_re_aghq()].
#'
#' Unlike [nmix_laplace()], which fixes the coefficients and returns a
#' point fit, this helper returns a *closure object* that evaluates the marginal
#' and its eta-level derivatives at arbitrary linear predictors. A random-effect
#' integrator perturbs the predictors by \eqn{Z b} per group and calls `eval()`
#' at each quadrature point; the per-site value, gradient and observed-
#' information block are everything the per-group integral needs.
#'
#' The abundance arm is per site (\eqn{\log\lambda_i = \eta^\lambda_i}) and the
#' detection arm is per visit (\eqn{\mathrm{logit}\,p_{ij} = \eta^p_{ij}}); the
#' two are coupled through the shared latent \eqn{N_i}. The per-site marginal
#' observed information in the eta coordinates
#' \eqn{(\eta^\lambda_i, \eta^p_{i1}, \dots, \eta^p_{iJ_i})} is
#' \deqn{B_i = \mathrm{diag}(I^\lambda_i, I^p_{ij}) - \mathrm{Var}(N_i\mid y_i)\,
#'       v_i v_i^\top, \qquad v_i = (-w_i, p_{i1}, \dots, p_{iJ_i}),}
#' where \eqn{I^\lambda_i}, \eqn{I^p_{ij}} are the complete-data Fisher diagonal,
#' \eqn{w_i} is the abundance score weight (`1` Poisson, \eqn{1-q_i} NB) and
#' \eqn{p_{ij} = \mathrm{plogis}(\eta^p_{ij})}. The off-diagonal
#' \eqn{\mathrm{Var}(N_i)\,w_i\,p_{ij}} is the abundance/detection coupling an
#' integrator placing random effects on both arms must carry (Louis 1982). This
#' is the eta-level form of the curvature [nmix_laplace()] sandwiches with
#' the design matrices.
#'
#' The single- vs multi-arm `make_site` adapter that wires this into a specific
#' integrator is deliberately not built here -- this object is the arm-agnostic
#' foundation it sits on.
#'
#' @inheritParams nmix_laplace
#' @param headroom Latent-N states summed above each site's own `max(y_i)`.
#'   `NULL` (default) derives it from `K_max`: an unset `K_max` caps each site at
#'   `max(y_i) + 100`, an explicit one truncates globally and uncapped. A
#'   negative value disables the per-site cap. Callers that resolve the ceiling
#'   themselves (the community fitters share one `K_max` across species) pass
#'   both.
#' @param K_site Per-site latent-N ceiling, an integer vector of length
#'   `n_sites` (`NULL` = none). The same ceiling `headroom` derives from
#'   `max(y_i)`, stated directly instead: a caller holding a fitted
#'   \eqn{\lambda_i} keys each site's ceiling to that site's abundance scale,
#'   which `max(y_i)` understates wherever detection is low. Takes precedence
#'   over `headroom`, and must be at least each site's own `max(y_i)`.
#'
#' @return An object of class `nmix_marginal`: a list of the validated
#'   data plus closures
#'   * `eval(eta_lambda, eta_p, r = Inf)` -- evaluate at the abundance log
#'     predictor `eta_lambda` (length `n_sites`) and detection logit predictor
#'     `eta_p` (length `n_obs`, visit order matching `y`). Returns a list with
#'     `log_lik` (total), `log_lik_site`, `grad_eta_lambda`, `grad_eta_p`,
#'     `grad_theta`, the complete-data Fisher `info_eta_lambda` / `info_eta_p`,
#'     the abundance score weight `score_wt_lambda`, `p` (= `plogis(eta_p)`),
#'     `mean_N`, `var_N`, `boundary_weight`, and (NB) the dispersion-coupling
#'     pieces `info_theta` / `info_lambda_theta` / `cov_N_stheta` / `var_stheta`.
#'   * `eval_beta(beta_lambda, beta_p, r = Inf)` -- the same, with the predictors
#'     formed from the stored design matrices.
#'   * `obs_info_block(s, ev)` -- the \eqn{(1+J_s)\times(1+J_s)} per-site
#'     marginal observed-information matrix \eqn{B_s} (above) for site `s`, from
#'     an `eval()` / `eval_beta()` result `ev`. Coordinates are
#'     `(eta_lambda_s, eta_p over site s's visits in input order)`.
#'   Plus `n_sites`, `n_obs`, `p_lambda`, `p_p`, `mixture`, `K_max`,
#'   `obs_by_site` (per-site visit row indices).
#'
#' @references
#' Royle, J. A. (2004). N-mixture models for estimating population size from
#'   spatially replicated counts. *Biometrics* 60, 108-115.
#' Louis, T. A. (1982). Finding the observed information matrix when using the
#'   EM algorithm. *JRSS-B* 44, 226-233.
#'
#' @seealso [nmix_laplace()] for the fixed-effects fit, [tulpa_re_aghq()]
#'   for the grouped random-effect integrator this feeds.
nmix_site_marginal <- function(y,
                                     site_idx,
                                     X_lambda,
                                     X_p,
                                     mixture = c("P", "NB"),
                                     K_max = NULL,
                                     headroom = NULL,
                                     K_site = NULL) {
  mixture <- match.arg(mixture)
  y        <- as.integer(y)
  site_idx <- as.integer(site_idx)
  if (!is.matrix(X_lambda)) stop("`X_lambda` must be a numeric matrix.", call. = FALSE)
  if (!is.matrix(X_p))      stop("`X_p` must be a numeric matrix.", call. = FALSE)
  if (length(y) != length(site_idx)) {
    stop("length(y) must equal length(site_idx).", call. = FALSE)
  }
  if (length(y) != nrow(X_p)) {
    stop("length(y) must equal nrow(X_p).", call. = FALSE)
  }
  if (any(y < 0L) || anyNA(y)) {
    stop("`y` must be nonnegative integers with no NA.", call. = FALSE)
  }
  n_sites <- nrow(X_lambda)
  n_obs   <- length(y)
  if (min(site_idx) < 1L || max(site_idx) > n_sites) {
    stop("`site_idx` values must lie in [1, nrow(X_lambda)].", call. = FALSE)
  }
  p_lambda <- ncol(X_lambda)
  p_p      <- ncol(X_p)

  # A caller that has already resolved the ceiling (the community fitters share
  # one K_max across species) passes `headroom` alongside it, so resolving here
  # does not read the pre-resolved K_max as an explicit user truncation.
  trunc <- .nmix_truncation(K_max, y)
  K_max <- trunc$K_max
  if (K_max < max(y)) stop("`K_max` must be >= max(y).", call. = FALSE)
  headroom <- if (is.null(headroom)) trunc$headroom else as.integer(headroom)

  if (!is.null(K_site)) {
    K_site <- as.integer(K_site)
    if (length(K_site) != n_sites) {
      stop(sprintf("`K_site` must have length n_sites (%d), got %d.",
                   n_sites, length(K_site)), call. = FALSE)
    }
    if (anyNA(K_site)) stop("`K_site` must not contain NA.", call. = FALSE)
    y_max_site <- rep(0L, n_sites)
    ms <- tapply(y, site_idx, max)
    y_max_site[as.integer(names(ms))] <- as.integer(ms)
    if (any(K_site < y_max_site)) {
      bad <- which(K_site < y_max_site)[1L]
      stop(sprintf("`K_site`[%d] = %d is below that site's max(y) = %d.",
                   bad, K_site[bad], y_max_site[bad]), call. = FALSE)
    }
  }

  obs_by_site <- split(seq_len(n_obs), site_idx)
  # Re-key to a dense 1..n_sites list (sites with no visits get integer(0)).
  obs_by_site <- lapply(seq_len(n_sites), function(s) {
    o <- obs_by_site[[as.character(s)]]
    if (is.null(o)) integer(0) else o
  })

  resolve_r <- function(r) {
    if (identical(mixture, "P")) return(Inf)
    if (is.null(r) || !is.finite(r) || r <= 0) {
      stop("NB marginal requires a finite positive `r` (NB size).", call. = FALSE)
    }
    as.numeric(r)
  }

  eval_eta <- function(eta_lambda, eta_p, r = Inf) {
    eta_lambda <- as.numeric(eta_lambda)
    eta_p      <- as.numeric(eta_p)
    if (length(eta_lambda) != n_sites) {
      stop("length(eta_lambda) must equal n_sites.", call. = FALSE)
    }
    if (length(eta_p) != n_obs) {
      stop("length(eta_p) must equal n_obs.", call. = FALSE)
    }
    out <- cpp_nmix_total_log_lik(y, site_idx, eta_p, eta_lambda, K_max,
                                  r = resolve_r(r), headroom = headroom,
                                  K_site = K_site)
    out$p <- plogis(eta_p)
    out
  }

  eval_beta <- function(beta_lambda, beta_p, r = Inf) {
    beta_lambda <- as.numeric(beta_lambda)
    beta_p      <- as.numeric(beta_p)
    if (length(beta_lambda) != p_lambda) {
      stop("length(beta_lambda) must equal ncol(X_lambda).", call. = FALSE)
    }
    if (length(beta_p) != p_p) {
      stop("length(beta_p) must equal ncol(X_p).", call. = FALSE)
    }
    eval_eta(as.numeric(X_lambda %*% beta_lambda),
             as.numeric(X_p %*% beta_p), r = r)
  }

  # Per-site marginal observed-information block B_s from an eval() result.
  # Coordinates: (eta_lambda_s, eta_p over site s's visits in input order).
  # B_s = diag(I_lam, I_p) - var_N * v v', v = (-score_wt_lambda, p_visits).
  obs_info_block <- function(s, ev) {
    obs <- obs_by_site[[s]]
    J <- length(obs)
    d <- 1L + J
    B <- matrix(0, d, d)
    B[1, 1] <- ev$info_eta_lambda[s]
    if (J > 0L) {
      diag(B)[-1] <- ev$info_eta_p[obs]
      vN <- ev$var_N[s]
      v  <- c(-ev$score_wt_lambda[s], ev$p[obs])
      B  <- B - vN * tcrossprod(v)
    }
    B
  }

  structure(
    list(
      y = y, site_idx = site_idx, X_lambda = X_lambda, X_p = X_p,
      mixture = mixture, K_max = K_max, headroom = headroom, K_site = K_site,
      n_sites = n_sites, n_obs = n_obs, p_lambda = p_lambda, p_p = p_p,
      obs_by_site = obs_by_site,
      eval = eval_eta, eval_beta = eval_beta, obs_info_block = obs_info_block
    ),
    class = c("nmix_marginal", "list")
  )
}

#' @exportS3Method print nmix_marginal
print.nmix_marginal <- function(x, ...) {
  cat(sprintf("tulpa N-mixture per-site marginal (mixture = %s)\n", x$mixture))
  cat(sprintf("  n_sites = %d   n_obs = %d   K_max = %d\n",
              x$n_sites, x$n_obs, x$K_max))
  cat(sprintf("  p_lambda = %d   p_p = %d\n", x$p_lambda, x$p_p))
  cat("  closures: eval(eta_lambda, eta_p, r), eval_beta(beta_lambda, beta_p, r), obs_info_block(s, ev)\n")
  invisible(x)
}
