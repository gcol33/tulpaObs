# distance_nuts.R - NUTS target density for the binned distance-sampling family.
#
# The flat coefficient vector is
#   theta = (beta_lambda [p_lam], beta_sigma [p_sig],
#            log_shape [hazard-rate only], log_r [NB only])
# and the joint log-posterior is the distance marginal (cpp_distance_total_log_lik)
# plus weak Gaussian priors. The C++ FullGradFn (src/distance_nuts.cpp) mirrors
# this R target and is cross-checked against it; this R version is the oracle.

# Parameter layout.
# `re_groups` > 0 appends a trailing [z_1..z_G, log_sigma_re] block (a single
# intercept RE on the abundance arm, tulpaObs#51).
.tobs_distance_nuts_layout <- function(p_lam, p_sig, hazard, is_nb, re_groups = 0L) {
  base <- p_lam + p_sig + (if (hazard) 1L else 0L) + (if (is_nb) 1L else 0L)
  off <- p_lam + p_sig
  log_shape <- if (hazard) { off <- off + 1L; off } else integer(0)
  log_r     <- if (is_nb)  off + 1L else integer(0)
  out <- list(p_lam = p_lam, p_sig = p_sig, hazard = isTRUE(hazard),
              is_nb = isTRUE(is_nb), lambda = seq_len(p_lam),
              sigma = p_lam + seq_len(p_sig), log_shape = log_shape, log_r = log_r,
              re_groups = as.integer(re_groups))
  if (re_groups > 0L) {
    out$z <- base + seq_len(re_groups); out$log_sigma <- base + re_groups + 1L
    out$total <- base + re_groups + 1L
  } else { out$z <- integer(0); out$log_sigma <- integer(0); out$total <- base }
  out
}

# Per-site distance marginal closure (the NUTS oracle's data + eval_beta).
# `headroom`: the per-site K_hi cap the caller's warm-start Laplace fit already
# verified against the shared ceiling (gcol33/tulpaObs#168); `eval_beta`'s own
# `headroom` argument lets a caller re-evaluate at a DIFFERENT cap (-1L, the
# uncapped comparison) without rebuilding the marginal, for the post-hoc guard
# in .tobs_fit_distance_nuts().
.tobs_distance_nuts_marginal <- function(model, mixture = "P", K_max = NULL,
                                         headroom = -1L) {
  headroom0 <- as.integer(headroom)
  y         <- model$y
  X_lambda  <- model$X_processes[[1]]
  X_sigma   <- model$X_processes[[2]]
  key_code  <- .dist_key_code(model$key)
  trans_code <- .dist_transect_code(model$transect)
  quad_order <- model$quad_order
  quad_xptr <- cpp_distance_build_quad(as.numeric(model$cutpoints), trans_code,
                                       as.integer(quad_order))
  if (is.null(K_max)) K_max <- 3L * max(rowSums(y)) + 100L
  K_max <- as.integer(K_max)
  resolve_r <- function(r) {
    if (identical(mixture, "P")) return(Inf)
    if (is.null(r) || !is.finite(r) || r <= 0)
      stop("NB distance marginal requires a finite positive `r`.", call. = FALSE)
    as.numeric(r)
  }
  eval_beta <- function(beta_lambda, beta_sigma, eta_b = 0, r = Inf,
                        headroom = headroom0) {
    eta_lambda <- as.numeric(X_lambda %*% beta_lambda)
    eta_sigma  <- as.numeric(X_sigma  %*% beta_sigma)
    cpp_distance_total_log_lik(y, eta_lambda, eta_sigma, as.numeric(eta_b),
                               quad_xptr, key_code, K_max, resolve_r(r),
                               headroom = as.integer(headroom))
  }
  list(X_lambda = X_lambda, X_sigma = X_sigma, K_max = K_max, mixture = mixture,
       headroom = headroom0, eval_beta = eval_beta)
}

# Joint log-posterior + gradient of the distance coefficient vector (R oracle).
.tobs_distance_nuts_logpost <- function(theta, marg, lay,
                                        sigma.beta = 10, sigma.shape = 1.5,
                                        sigma.logr = 1.5) {
  beta_lambda <- theta[lay$lambda]
  beta_sigma  <- theta[lay$sigma]
  eta_b <- if (lay$hazard) theta[lay$log_shape] else 0
  r     <- if (lay$is_nb)  exp(theta[lay$log_r]) else Inf
  ev <- marg$eval_beta(beta_lambda, beta_sigma, eta_b = eta_b, r = r)

  lp   <- ev$log_lik
  grad <- c(as.numeric(crossprod(marg$X_lambda, ev$grad_eta_lambda)),
            as.numeric(crossprod(marg$X_sigma,  ev$grad_eta_sigma)))
  if (lay$hazard) grad <- c(grad, ev$grad_eta_b)
  if (lay$is_nb)  grad <- c(grad, sum(ev$grad_theta))

  ib2 <- 1 / sigma.beta^2
  lp  <- lp - 0.5 * ib2 * (sum(beta_lambda^2) + sum(beta_sigma^2))
  grad[lay$lambda] <- grad[lay$lambda] - ib2 * beta_lambda
  grad[lay$sigma]  <- grad[lay$sigma]  - ib2 * beta_sigma
  if (lay$hazard) {
    is2 <- 1 / sigma.shape^2
    lp  <- lp - 0.5 * is2 * eta_b^2
    grad[lay$log_shape] <- grad[lay$log_shape] - is2 * eta_b
  }
  if (lay$is_nb) {
    lr  <- theta[lay$log_r]; ilr2 <- 1 / sigma.logr^2
    lp  <- lp - 0.5 * ilr2 * lr^2
    grad[lay$log_r] <- grad[lay$log_r] - ilr2 * lr
  }
  list(lp = lp, grad = grad)
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the distance family
# ---------------------------------------------------------------------------

.tobs_fit_distance_nuts <- function(model, mixture = "poisson", K_max = NULL,
                                    headroom = NULL,
                                    sigma.beta = NULL, sigma.shape = 1.5,
                                    sigma.logr = NULL, re = NULL,
                                    n.iter = NULL, n.warmup = NULL, n.chains = NULL,
                                    max.treedepth = NULL, adapt.delta = NULL,
                                    seed = NULL, verbose = FALSE) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188).
  .tobs_fill_sampler(environment(), "nuts", single_species = TRUE)

  is_nb    <- identical(mixture, "negbin")
  mix_code <- if (is_nb) "NB" else "P"
  hazard   <- identical(model$key, "hazard")
  X_lambda <- model$X_processes[[1]]
  X_sigma  <- model$X_processes[[2]]
  y        <- model$y
  p_lam <- ncol(X_lambda); p_sig <- ncol(X_sigma)
  trunc    <- .dist_truncation(K_max, rowSums(y))
  K_max    <- trunc$K_max
  headroom <- if (is.null(headroom)) trunc$headroom else as.integer(headroom)
  # NB's heavier tail is not sized for the per-site window; keep the shared
  # ceiling under NB (mirrors nmix_laplace(), R/nmix_laplace.R).
  if (is_nb) headroom <- -1L

  # Single intercept RE on the abundance arm (tulpaObs#51), via the shared
  # count-NUTS RE helpers. distance's detection arm is the log-sigma scale; an RE
  # there is not wired, so only the abundance (lambda) arm is supported.
  re_info <- .tobs_count_nuts_re_info(re, model)
  if (!is.null(re_info) && re_info$arm != 0L)
    stop("distance() NUTS supports a random effect on the abundance arm only; ",
         "a detection-scale (sigma) RE is not wired. Put the RE on the state ",
         "formula, or use method = \"laplace\".", call. = FALSE)
  n_re_groups <- if (!is.null(re_info)) re_info$n_groups else 0L
  lay <- .tobs_distance_nuts_layout(p_lam, p_sig, hazard, is_nb,
                                    re_groups = n_re_groups)

  warm <- distance_laplace(y = y, X_lambda = X_lambda, X_sigma = X_sigma,
                           cutpoints = model$cutpoints, key = model$key,
                           transect = model$transect, mixture = mix_code,
                           K_max = K_max, headroom = headroom,
                           quad_order = model$quad_order,
                           max_iter = 100L, verbose = FALSE)
  # The warm-start Laplace fit's OWN score-gap guard already verified (and, if
  # needed, widened) this headroom at its mode; every leapfrog gradient below
  # reuses that verified window rather than re-checking it thousands of times
  # per chain (gcol33/tulpaObs#168). A single post-hoc check at the posterior
  # mean, after sampling, is the guard for the region NUTS actually explored.
  headroom <- warm$headroom %||% headroom
  theta0 <- c(as.numeric(warm$beta_lambda), as.numeric(warm$beta_sigma))
  if (hazard) theta0 <- c(theta0, as.numeric(warm$eta_b))
  if (is_nb)  theta0 <- c(theta0, as.numeric(warm$log_r))
  V <- as.matrix(warm$vcov)
  base_metric <- if (!is.null(V) && nrow(V) == length(theta0)) pmax(diag(V), 1e-6)
                 else rep(1, length(theta0))
  init <- .tobs_count_nuts_re_init(list(theta0 = theta0, inv_metric = base_metric),
                                   lay, re_info)
  theta0 <- init$theta0; inv_metric <- init$inv_metric

  quad_xptr <- cpp_distance_build_quad(as.numeric(model$cutpoints),
                                       .dist_transect_code(model$transect),
                                       as.integer(model$quad_order))
  spec <- .tobs_count_nuts_re_spec(
    list(y = y, X_lambda = X_lambda, X_sigma = X_sigma,
         quad_xptr = quad_xptr,
         key = .dist_key_code(model$key), K_max = K_max,
         is_nb = is_nb, headroom = as.integer(headroom)),
    re_info, sigma.logr)

  run_chain <- function(ch) {
    cpp_distance_nuts(spec, theta0 = theta0,
                      sigma_beta = sigma.beta, sigma_shape = sigma.shape,
                      sigma_logr = sigma.logr, inv_metric = inv_metric,
                      n_iter = as.integer(n.iter + n.warmup),
                      n_warmup = as.integer(n.warmup),
                      max_treedepth = as.integer(max.treedepth),
                      adapt_delta = adapt.delta,
                      seed = as.integer(seed + ch - 1L), verbose = isTRUE(verbose))
  }
  nms <- c(paste0("lambda_", model$process_info[[1]]$coef_names),
           paste0("sigma_",  model$process_info[[2]]$coef_names),
           if (hazard) "log_shape", if (is_nb) "log_r",
           .tobs_count_nuts_re_names(re_info))
  run <- .tobs_count_nuts_run(run_chain, n.chains, nms)
  par <- run$par; cov <- run$cov

  marg <- .tobs_distance_nuts_marginal(model, mixture = mix_code, K_max = K_max,
                                       headroom = headroom)
  eta_b_hat <- if (hazard) par[lay$log_shape] else 0
  r_hat <- if (is_nb) exp(par[lay$log_r]) else Inf
  ll_mean <- marg$eval_beta(par[lay$lambda], par[lay$sigma],
                            eta_b = eta_b_hat, r = r_hat)$log_lik

  # Post-hoc guard at the posterior mean (gcol33/tulpaObs#168): the capped
  # headroom was verified at the WARM-START mode, not at the region NUTS
  # actually sampled, so check once more here and re-run the whole chain set
  # at a wider window on disagreement -- the same escalate-and-refit the other
  # distance fitters use, just run once per chain set instead of per leapfrog
  # step (a fixed headroom is cheap to evaluate at every gradient; verifying it
  # at every gradient is not).
  if (headroom >= 0L) {
    gap <- tryCatch({
      ev_h <- marg$eval_beta(par[lay$lambda], par[lay$sigma], eta_b = eta_b_hat,
                             r = r_hat, headroom = headroom)
      ev_u <- marg$eval_beta(par[lay$lambda], par[lay$sigma], eta_b = eta_b_hat,
                             r = r_hat, headroom = -1L)
      .dist_score_gap(ev_h$grad_eta_lambda, ev_u$grad_eta_lambda, X_lambda,
                      ev_h$grad_eta_sigma,  ev_u$grad_eta_sigma,  X_sigma)
    }, error = function(e) NA_real_)
    if (is.finite(gap) && gap > .NMIX_SCORE_TOL) {
      h_next <- .nmix_widen_headroom(headroom, K_max)
      if (!is.null(h_next)) {
        cl <- match.call()
        cl$headroom <- h_next
        return(eval(cl, parent.frame()))
      }
    }
  }

  raw <- list(
    mixture = mix_code, key = model$key, transect = model$transect,
    hazard = hazard, nb = is_nb,
    beta_lambda = unname(par[lay$lambda]),
    beta_sigma  = unname(par[lay$sigma]),
    eta_b  = if (hazard) unname(par[lay$log_shape]) else NA_real_,
    shape  = if (hazard) exp(unname(par[lay$log_shape])) else NA_real_,
    log_r  = if (is_nb) unname(par[lay$log_r]) else NA_real_,
    r      = if (is_nb) exp(unname(par[lay$log_r])) else NA_real_,
    vcov = cov, log_lik = ll_mean, converged = TRUE, K_max = K_max,
    mean_N = warm$mean_N, var_N = warm$var_N, p_det = warm$p_det,
    boundary_weight = warm$boundary_weight)
  fit <- build_distance_fit(raw, model)
  fit$headroom <- headroom

  .tobs_count_nuts_attach(
    fit, run, ll_mean, n.chains, re_info,
    extra = list(is_nb = is_nb, hazard = hazard, K_max = K_max,
                 headroom = headroom,
                 sigma_beta = sigma.beta, sigma_shape = sigma.shape,
                 sigma_logr = sigma.logr))
}
