# royle_nichols.R - Royle-Nichols occupancy (Royle & Nichols 2003; unmarked
# occuRN). Detection heterogeneity is induced by latent abundance:
#
#   N_i           ~ Poisson(lambda_i),   log lambda_i = X_lambda_i . beta_lambda
#   detect at j   ~ Bernoulli(1 - (1 - r_i)^{N_i}),  logit r_i = X_r_i . beta_r
#
# r_i is the per-individual per-visit detection probability. With r_i constant
# across a site's visits the visits are exchangeable, so the site sufficient
# statistics are (k_i detections, n_i valid visits) and the latent N marginalises
# in closed form (a Poisson sum to K_max):
#
#   L_i = sum_{N=0}^{K} dpois(N, lambda_i)
#           * [1 - (1 - r_i)^N]^{k_i} * [(1 - r_i)^N]^{n_i - k_i}
#
# The fit maximises the exact marginal (optim BFGS) and takes the observed-
# information vcov from the Hessian -- the same closed-form-marginal + observed-
# Fisher recipe as abun()/removal()/distance(), here in R because the marginal is
# a light Poisson sum over per-site sufficient statistics. Detection is
# site-level; visit-varying detection is a documented follow-up.
#
#   .tobs_build_royle_nichols()   data binder -> model_type = "royle_nichols"
#   .tobs_fit_royle_nichols()     optim over the closed-form marginal
#   .dispatch_royle_nichols()     tobs() entry (bind + fit + assemble)

# ---------------------------------------------------------------------------
# Marginal log-likelihood (per site), the single source of truth reused by the
# fitter, the pointwise log-likelihood, and simulate().
# ---------------------------------------------------------------------------

# Per-site marginal log-likelihood at given lambda / r vectors and per-site
# sufficient statistics (kk detections, nn valid visits). Returns length-n_sites.
.rn_site_loglik <- function(lambda, r, kk, nn, K) {
  Ns     <- 0:K
  log1mr <- log1p(-r)                                   # log(1 - r), length S
  # (1 - r)^N and log(1 - (1 - r)^N) as [S x (K+1)] matrices.
  N_log1mr <- outer(log1mr, Ns, "*")                    # N * log(1 - r)
  qN       <- exp(N_log1mr)                             # (1 - r)^N in [0, 1]
  a        <- log1p(-qN)                                # log(1 - (1 - r)^N); -Inf at N=0
  ka       <- sweep(a, 1L, kk, "*")                     # k * log(1 - (1-r)^N)
  ka[is.nan(ka)] <- 0                                   # 0 * -Inf (k=0, N=0) -> 0
  logpois  <- outer(lambda, Ns, function(l, N) stats::dpois(N, l, log = TRUE))
  term     <- logpois + ka + sweep(N_log1mr, 1L, (nn - kk), "*")
  m  <- apply(term, 1L, max)                            # log-sum-exp over N
  m + log(rowSums(exp(term - m)))
}

# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# `y` is an n_sites x max_visits detection matrix (0/1/NA). occ_formula models
# log lambda (abundance); det_formula models logit r (site-level detection).
.tobs_build_royle_nichols <- function(occ_formula, det_formula, data, y,
                                      K_max = NULL) {
  y_int <- matrix(as.integer(round(y)), nrow(y), ncol(y))
  vs    <- !is.na(y_int)
  if (any(y_int[vs] != 0L & y_int[vs] != 1L)) {
    stop("royle_nichols() y must contain only 0, 1, or NA (detection history).",
         call. = FALSE)
  }
  n_sites    <- nrow(y_int)
  max_visits <- ncol(y_int)
  .tobs_check_site_count(n_sites, nrow(data), "sites")

  kk <- rowSums(y_int == 1L & vs)                       # detections per site
  nn <- rowSums(vs)                                     # valid visits per site

  bind     <- .tobs_bind_formulas(list(lambda = occ_formula, r = det_formula),
                                  data)
  X_lambda <- stats::model.matrix(bind$fe$lambda, data)
  X_r      <- stats::model.matrix(bind$fe$r, data)

  # K_max must cover the Poisson upper tail of the largest plausible lambda; a
  # naive moment guess (below) informs a floor if the caller gives none.
  if (is.null(K_max)) {
    p_hat  <- mean(kk / pmax(nn, 1))
    r0     <- min(max(p_hat, 0.05), 0.7)
    lam0   <- -log(max(1 - p_hat, 1e-3)) / r0
    K_max  <- max(30L, as.integer(stats::qpois(0.9999, max(lam0, 1)) + 10L))
  }

  structure(list(
    model_type  = "royle_nichols",
    y           = y_int,
    k_site      = as.integer(kk),
    n_site      = as.integer(nn),
    K_max       = as.integer(K_max),
    X_processes = list(X_lambda, X_r),
    formulas    = list(lambda = bind$fe$lambda, r = bind$fe$r),
    structured_terms = bind$terms,
    data        = data,
    n_sites     = n_sites,
    max_visits  = max_visits,
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda),
           coef_names = colnames(X_lambda), link = "log"),
      list(name = "r",      p = ncol(X_r),
           coef_names = colnames(X_r), link = "logit")
    )
  ), class = "tobs_model")
}

# ---------------------------------------------------------------------------
# Fitter: maximise the closed-form marginal, observed-information vcov
# ---------------------------------------------------------------------------

.tobs_fit_royle_nichols <- function(model, verbose = TRUE, ...) {
  X_lambda <- model$X_processes[[1L]]
  X_r      <- model$X_processes[[2L]]
  p_lam    <- ncol(X_lambda); p_r <- ncol(X_r)
  kk <- model$k_site; nn <- model$n_site; K <- model$K_max

  nll <- function(theta) {
    lambda <- exp(as.vector(X_lambda %*% theta[seq_len(p_lam)]))
    r      <- stats::plogis(as.vector(X_r %*% theta[p_lam + seq_len(p_r)]))
    r      <- pmin(pmax(r, 1e-8), 1 - 1e-8)
    ll     <- .rn_site_loglik(lambda, r, kk, nn, K)
    val    <- -sum(ll)
    if (!is.finite(val)) 1e10 else val
  }

  # Moment initialisation: p_detect = 1 - exp(-lambda * r) at the intercept.
  p_hat <- mean(kk / pmax(nn, 1))
  r0    <- min(max(p_hat, 0.05), 0.7)
  lam0  <- -log(max(1 - p_hat, 1e-3)) / r0
  init  <- c(log(max(lam0, 1e-2)), rep(0, p_lam - 1L),
             stats::qlogis(r0),     rep(0, p_r - 1L))

  opt <- stats::optim(init, nll, method = "BFGS", hessian = TRUE,
                      control = list(maxit = 500L))
  converged <- opt$convergence == 0L

  par_names <- c(paste0("lambda_", model$process_info[[1L]]$coef_names),
                 paste0("r_",      model$process_info[[2L]]$coef_names))
  means <- opt$par; names(means) <- par_names
  V <- tryCatch(solve(opt$hessian), error = function(e) {
    diag(NA_real_, length(means))
  })
  V <- (V + t(V)) / 2
  dimnames(V) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  intercepts <- list(
    lambda = stats::setNames(means[1L], par_names[1L]),
    r      = stats::setNames(means[p_lam + 1L], par_names[p_lam + 1L]))

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(-opt$value, n_draws),
    log_lik      = -opt$value,
    N            = sum(nn > 0L)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    n_fixed      = length(means),
    fixed_names  = par_names,
    intercepts   = intercepts,
    process_info = model$process_info,
    model        = model,
    spatial      = NULL,
    method       = "laplace",
    convergence  = list(converged = converged, n_iter = opt$counts[[1L]])
  )), class = c("tobs_fit", "tulpa_fit"))
}

# ---------------------------------------------------------------------------
# tobs() dispatcher
# ---------------------------------------------------------------------------

.dispatch_royle_nichols <- function(formula, data, family, detection, y, visits,
                                    engine, priors, control,
                                    approx = "gaussian_laplace",
                                    correction = "none", ...) {
  if (is.null(detection)) {
    stop("royle_nichols() requires a `detection` formula (site-level ",
         "per-individual detection r).", call. = FALSE)
  }
  if (is.null(y)) {
    stop("royle_nichols() requires `y` (an N x J 0/1 detection-history ",
         "matrix).", call. = FALSE)
  }
  if (!is.null(visits)) {
    stop("royle_nichols() detection is site-level; visit-level detection ",
         "covariates (`visits`) are not yet supported.", call. = FALSE)
  }
  if (!identical(.map_engine(engine, family = "royle_nichols"), "laplace")) {
    stop("royle_nichols() supports method = \"laplace\" only.", call. = FALSE)
  }
  model <- .tobs_build_royle_nichols(
    occ_formula = formula, det_formula = detection, data = data, y = y,
    K_max = family$params$K_max)
  .tobs_fit_royle_nichols(model, verbose = isTRUE(control$verbose))
}

# ---------------------------------------------------------------------------
# fitted / predict / residuals / pointwise log-likelihood
# ---------------------------------------------------------------------------

# fitted(): per-site lambda (abundance), r (per-individual detection), and the
# marginal per-visit detection probability p = 1 - exp(-lambda * r).
.tobs_fitted_royle_nichols <- function(object) {
  model  <- object$model
  p_lam  <- model$process_info[[1L]]$p
  beta   <- object$means
  lambda <- exp(as.vector(model$X_processes[[1L]] %*% beta[seq_len(p_lam)]))
  r      <- stats::plogis(as.vector(model$X_processes[[2L]] %*%
                                      beta[p_lam + seq_len(model$process_info[[2L]]$p)]))
  list(lambda = lambda, r = r, p = 1 - exp(-lambda * r))
}

# predict(): lambda (abundance, default) or r (detection) at the fitted or new
# design, on the response scale.
.tobs_predict_royle_nichols <- function(object, newdata = NULL,
                                        type = c("abundance", "detection")) {
  type  <- match.arg(type)
  model <- object$model
  k     <- if (identical(type, "detection")) 2L else 1L
  X     <- if (is.null(newdata)) model$X_processes[[k]]
           else stats::model.matrix(model$formulas[[if (k == 1L) "lambda" else "r"]],
                                    newdata)
  off   <- if (k == 1L) 0L else model$process_info[[1L]]$p
  eta   <- as.vector(X %*% object$means[off + seq_len(model$process_info[[k]]$p)])
  if (identical(type, "detection")) stats::plogis(eta) else exp(eta)
}

# residuals(): per-site response / pearson / deviance on the observed per-site
# detection frequency (k / n) against the marginal per-visit detection prob p.
.tobs_residuals_royle_nichols <- function(object, type) {
  model <- object$model
  fv    <- .tobs_fitted_royle_nichols(object)
  nn    <- model$n_site; kk <- model$k_site
  obs   <- ifelse(nn > 0L, kk / nn, NA_real_)
  p     <- fv$p; eps <- 1e-10
  res <- switch(type,
    response = obs - p,
    pearson  = (obs - p) / sqrt(p * (1 - p) / pmax(nn, 1) + eps),
    deviance = sign(obs - p) * sqrt(2 * pmax(nn, 1) * abs(
      ifelse(obs > 0, obs * log(pmax(obs, eps) / p), 0) +
      ifelse(obs < 1, (1 - obs) * log((1 - obs + eps) / (1 - p + eps)), 0))))
  list(occ = res, det = NULL)
}

# Posterior replicate detection histories: draw a coefficient vector, then per
# site N ~ Poisson(lambda) and y_ij ~ Bernoulli(1 - (1 - r)^N) at the observed
# visit pattern. Returns an n_sites x max_visits matrix (or a list for nsim > 1).
.tobs_simulate_royle_nichols <- function(object, nsim = 1) {
  model <- object$model
  p_lam <- model$process_info[[1L]]$p; p_r <- model$process_info[[2L]]$p
  Xl <- model$X_processes[[1L]]; Xr <- model$X_processes[[2L]]
  valid <- !is.na(model$y)
  draw_one <- function() {
    idx    <- sample.int(nrow(object$draws), 1L)
    beta   <- object$draws[idx, ]
    lambda <- exp(as.vector(Xl %*% beta[seq_len(p_lam)]))
    r      <- stats::plogis(as.vector(Xr %*% beta[p_lam + seq_len(p_r)]))
    Ni     <- stats::rpois(model$n_sites, lambda)
    p_i    <- 1 - (1 - r)^Ni
    yy     <- matrix(0L, model$n_sites, model$max_visits)
    for (i in seq_len(model$n_sites)) {
      yy[i, ] <- stats::rbinom(model$max_visits, 1L, p_i[i])
    }
    yy[!valid] <- NA_integer_
    yy
  }
  if (nsim == 1L) return(draw_one())
  lapply(seq_len(nsim), function(s) draw_one())
}

# Pointwise log-likelihood [n_draws x n_sites] over the posterior draws, the
# exact per-site marginal (WAIC / LOO). Reuses .rn_site_loglik.
.tobs_ploglik_royle_nichols <- function(object, n.draws = 1000L, n.threads = 1L) {
  model <- object$model
  draws <- object$draws
  if (!is.null(n.draws) && as.integer(n.draws) < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  p_lam <- model$process_info[[1L]]$p; p_r <- model$process_info[[2L]]$p
  Xl <- model$X_processes[[1L]]; Xr <- model$X_processes[[2L]]
  kk <- model$k_site; nn <- model$n_site; K <- model$K_max
  t(vapply(seq_len(nrow(draws)), function(d) {
    lambda <- exp(as.vector(Xl %*% draws[d, seq_len(p_lam)]))
    r      <- stats::plogis(as.vector(Xr %*% draws[d, p_lam + seq_len(p_r)]))
    r      <- pmin(pmax(r, 1e-8), 1 - 1e-8)
    .rn_site_loglik(lambda, r, kk, nn, K)
  }, numeric(model$n_sites)))
}
