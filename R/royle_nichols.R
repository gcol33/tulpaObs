# royle_nichols.R - Royle-Nichols occupancy (Royle & Nichols 2003; unmarked
# occuRN). Detection heterogeneity is induced by latent abundance:
#
#   N_i           ~ Poisson(lambda_i),   log lambda_i = X_lambda_i . beta_lambda
#   detect at j   ~ Bernoulli(1 - (1 - r_ij)^{N_i}),  logit r_ij = X_r_ij . beta_r
#
# r_ij is the per-individual per-visit detection probability. With r_ij constant
# across a site's visits the visits are exchangeable, so the site sufficient
# statistics are (k_i detections, n_i valid visits) and the latent N marginalises
# in closed form (a Poisson sum to K_max):
#
#   L_i = sum_{N=0}^{K} dpois(N, lambda_i)
#           * [1 - (1 - r_i)^N]^{k_i} * [(1 - r_i)^N]^{n_i - k_i}
#
# When detection varies by visit (a `visits` covariate) the (k_i, n_i) reduction
# no longer holds, so the emission is the full per-visit product inside the same
# Poisson sum:
#
#   L_i = sum_{N=0}^{K} dpois(N, lambda_i)
#           * prod_{j: y_ij=1} [1 - (1 - r_ij)^N]
#           * prod_{j: y_ij=0} [(1 - r_ij)^N]
#
# which reduces to the (k_i, n_i) form when r_ij is constant. The fit maximises
# the exact marginal (optim BFGS) and takes the observed-information vcov from the
# Hessian -- the same closed-form-marginal + observed-Fisher recipe as
# abun()/removal()/distance(), here in R because the marginal is a light Poisson
# sum over the latent count.
#
#   .tobs_build_royle_nichols()   data binder -> model_type = "royle_nichols"
#   .tobs_fit_royle_nichols()     optim over the closed-form marginal
#   .dispatch_royle_nichols()     tobs() entry (bind + fit + assemble)

# ---------------------------------------------------------------------------
# Marginal log-likelihood (per site), the single source of truth reused by the
# fitter, the pointwise log-likelihood, and simulate().
# ---------------------------------------------------------------------------

# Log-sum-exp over the latent-N axis. `term` is [n_sites x (K+1)] (row i, column
# N+1 holds the log joint of N and the site's data); returns length n_sites.
.rn_lse_rows <- function(term) {
  m <- apply(term, 1L, max)
  m + log(rowSums(exp(term - m)))
}

# Site-level marginal: r is one per-individual detection probability per site and
# the data enter through the per-site sufficient statistics (kk detections, nn
# valid visits). Returns length n_sites.
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
  .rn_lse_rows(term)
}

# Visit-level marginal: detection varies by visit. `r_long`, `y_long`, `site_idx`
# are over VALID visits only (one entry per (site, visit) with non-missing y);
# `site_idx` indexes into 1:n_sites. A site with no valid visit contributes 0.
# Returns length n_sites. Reduces to .rn_site_loglik when r_long is constant
# within each site.
.rn_visit_loglik <- function(lambda, r_long, y_long, site_idx, K, n_sites) {
  Ns       <- 0:K
  log1mr   <- log1p(-r_long)                            # log(1 - r_ij), length V
  N_log1mr <- outer(log1mr, Ns, "*")                    # N * log(1 - r_ij), [V x (K+1)]
  logp     <- log1p(-exp(N_log1mr))                     # log(1 - (1-r_ij)^N); -Inf at N=0
  yp       <- y_long * logp                             # detected: log(1 - (1-r)^N)
  yp[is.nan(yp)] <- 0                                   # 0 * -Inf (y=0, N=0) -> 0
  contrib  <- yp + (1 - y_long) * N_log1mr              # per-visit log-emission
  # Sum the visit emissions within each site into a dense [n_sites x (K+1)].
  dc <- matrix(0, n_sites, K + 1L)
  rs <- rowsum(contrib, site_idx)
  dc[as.integer(rownames(rs)), ] <- rs
  logpois <- outer(lambda, Ns, function(l, N) stats::dpois(N, l, log = TRUE))
  .rn_lse_rows(logpois + dc)
}

# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# `y` is an n_sites x max_visits detection matrix (0/1/NA). occ_formula models
# log lambda (abundance); det_formula models logit r at the site level; an
# optional det_visit_formula + det_visit_data (from a `visits` argument) adds
# visit-level detection covariates, unrolled long-form over (site, visit).
.tobs_build_royle_nichols <- function(occ_formula, det_formula, data, y,
                                      K_max = NULL,
                                      det_visit_formula = NULL,
                                      det_visit_data = NULL) {
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

  # Visit-level detection design (site-major grid, intercept dropped, NA -> 0).
  X_r_visit <- .tobs_build_visit_X(det_visit_formula, det_visit_data,
                                   n_sites, max_visits, arm = "detection")
  visit_varying <- !is.null(X_r_visit)

  # Long-form (site, visit) bookkeeping for the visit-varying marginal: the grid
  # is site-major (site 1's visits, then site 2's, ...), matching the flatten in
  # .normalize_visits(); keep only the valid (non-missing y) cells.
  valid_long   <- NULL; vis_site_idx <- NULL; vis_y <- NULL; site_of_cell <- NULL
  if (visit_varying) {
    site_of_cell <- rep(seq_len(n_sites), each = max_visits)
    valid_long   <- as.vector(t(vs))
    vis_site_idx <- site_of_cell[valid_long]
    vis_y        <- as.vector(t(y_int))[valid_long]
  }

  # K_max must cover the Poisson upper tail of the largest plausible lambda; a
  # naive moment guess (below) informs a floor if the caller gives none.
  if (is.null(K_max)) {
    p_hat  <- mean(kk / pmax(nn, 1))
    r0     <- min(max(p_hat, 0.05), 0.7)
    lam0   <- -log(max(1 - p_hat, 1e-3)) / r0
    K_max  <- max(30L, as.integer(stats::qpois(0.9999, max(lam0, 1)) + 10L))
  }

  r_info <- list(name = "r", p = ncol(X_r),
                 coef_names = colnames(X_r), link = "logit")
  if (visit_varying) {
    r_info$visit_coef_names <- colnames(X_r_visit)
    r_info$p_visit          <- ncol(X_r_visit)
  }

  structure(list(
    model_type   = "royle_nichols",
    y            = y_int,
    k_site       = as.integer(kk),
    n_site       = as.integer(nn),
    K_max        = as.integer(K_max),
    X_processes  = list(X_lambda, X_r),
    X_r_visit    = X_r_visit,
    visit_varying = visit_varying,
    site_of_cell = site_of_cell,
    valid_long   = valid_long,
    vis_site_idx = vis_site_idx,
    vis_y        = vis_y,
    formulas     = list(lambda = bind$fe$lambda, r = bind$fe$r,
                        r_visit = det_visit_formula),
    structured_terms = bind$terms,
    data         = data,
    n_sites      = n_sites,
    max_visits   = max_visits,
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda),
           coef_names = colnames(X_lambda), link = "log"),
      r_info
    )
  ), class = "tobs_model")
}

# ---------------------------------------------------------------------------
# Fitter: maximise the closed-form marginal, observed-information vcov
# ---------------------------------------------------------------------------

# Per-site detection linear predictor over the VALID visits, given the coefficient
# blocks. Site-level: r is constant within a site. Visit-varying: the site term is
# broadcast over the site's visits and the visit-level term added.
.rn_r_long <- function(model, beta_r_site, beta_r_visit) {
  X_r      <- model$X_processes[[2L]]
  eta_site <- as.vector(X_r %*% beta_r_site)             # length n_sites
  eta_cell <- eta_site[model$site_of_cell] +
    as.vector(model$X_r_visit %*% beta_r_visit)          # length n_sites * max_visits
  r <- stats::plogis(eta_cell[model$valid_long])
  pmin(pmax(r, 1e-8), 1 - 1e-8)
}

.tobs_fit_royle_nichols <- function(model, verbose = TRUE, ...) {
  X_lambda <- model$X_processes[[1L]]
  X_r      <- model$X_processes[[2L]]
  p_lam    <- ncol(X_lambda); p_rs <- ncol(X_r)
  p_rv     <- if (model$visit_varying) ncol(model$X_r_visit) else 0L
  kk <- model$k_site; nn <- model$n_site; K <- model$K_max

  nll <- function(theta) {
    lambda <- exp(as.vector(X_lambda %*% theta[seq_len(p_lam)]))
    beta_rs <- theta[p_lam + seq_len(p_rs)]
    if (model$visit_varying) {
      beta_rv <- theta[p_lam + p_rs + seq_len(p_rv)]
      r_long  <- .rn_r_long(model, beta_rs, beta_rv)
      ll      <- .rn_visit_loglik(lambda, r_long, model$vis_y,
                                  model$vis_site_idx, K, model$n_sites)
    } else {
      r  <- stats::plogis(as.vector(X_r %*% beta_rs))
      r  <- pmin(pmax(r, 1e-8), 1 - 1e-8)
      ll <- .rn_site_loglik(lambda, r, kk, nn, K)
    }
    val <- -sum(ll)
    if (!is.finite(val)) 1e10 else val
  }

  # Moment initialisation: p_detect = 1 - exp(-lambda * r) at the intercept.
  p_hat <- mean(kk / pmax(nn, 1))
  r0    <- min(max(p_hat, 0.05), 0.7)
  lam0  <- -log(max(1 - p_hat, 1e-3)) / r0
  init  <- c(log(max(lam0, 1e-2)), rep(0, p_lam - 1L),
             stats::qlogis(r0),     rep(0, p_rs - 1L),
             rep(0, p_rv))

  opt <- stats::optim(init, nll, method = "BFGS", hessian = TRUE,
                      control = list(maxit = 500L))
  converged <- opt$convergence == 0L

  r_visit_names <- if (p_rv) paste0("r_", model$process_info[[2L]]$visit_coef_names)
  par_names <- c(paste0("lambda_", model$process_info[[1L]]$coef_names),
                 paste0("r_",      model$process_info[[2L]]$coef_names),
                 r_visit_names)
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
    stop("royle_nichols() requires a `detection` formula (per-individual ",
         "detection r).", call. = FALSE)
  }
  if (is.null(y)) {
    stop("royle_nichols() requires `y` (an N x J 0/1 detection-history ",
         "matrix).", call. = FALSE)
  }
  if (!identical(.map_engine(engine, family = "royle_nichols"), "laplace")) {
    stop("royle_nichols() supports method = \"laplace\" only.", call. = FALSE)
  }
  # Split `detection` into a site-level design and (if `visits` is supplied) a
  # visit-level design, exactly as the occupancy / N-mixture front doors do.
  vd <- .normalize_visits(visits, detection, n_sites = nrow(y),
                          max_visits = ncol(y))
  model <- .tobs_build_royle_nichols(
    occ_formula = formula, det_formula = vd$det_formula, data = data, y = y,
    K_max = family$params$K_max,
    det_visit_formula = vd$det_visit_formula, det_visit_data = vd$visits)
  .tobs_fit_royle_nichols(model, verbose = isTRUE(control$verbose))
}

# ---------------------------------------------------------------------------
# fitted / predict / residuals / pointwise log-likelihood
# ---------------------------------------------------------------------------

# Per-cell (site x visit) detection probability r on the grid, as an
# [n_sites x max_visits] matrix (NA at missing visits). Only meaningful for a
# visit-varying fit.
.rn_r_grid <- function(object) {
  model   <- object$model
  p_lam   <- model$process_info[[1L]]$p
  p_rs    <- model$process_info[[2L]]$p
  beta    <- object$means
  eta_s   <- as.vector(model$X_processes[[2L]] %*% beta[p_lam + seq_len(p_rs)])
  eta_c   <- eta_s[model$site_of_cell] +
    as.vector(model$X_r_visit %*% beta[p_lam + p_rs + seq_len(model$process_info[[2L]]$p_visit)])
  r_mat   <- matrix(stats::plogis(eta_c), model$n_sites, model$max_visits,
                    byrow = TRUE)
  r_mat[is.na(model$y)] <- NA_real_
  r_mat
}

# fitted(): per-site lambda (abundance) and per-individual detection r, plus the
# marginal per-visit detection probability p = 1 - exp(-lambda * r). Site-level
# detection returns r / p as length-n_sites vectors; visit-varying detection
# returns them as [n_sites x max_visits] matrices.
.tobs_fitted_royle_nichols <- function(object) {
  model  <- object$model
  p_lam  <- model$process_info[[1L]]$p
  beta   <- object$means
  lambda <- exp(as.vector(model$X_processes[[1L]] %*% beta[seq_len(p_lam)]))
  if (isTRUE(model$visit_varying)) {
    r <- .rn_r_grid(object)
    list(lambda = lambda, r = r, p = 1 - exp(-lambda * r))
  } else {
    r <- stats::plogis(as.vector(model$X_processes[[2L]] %*%
                                   beta[p_lam + seq_len(model$process_info[[2L]]$p)]))
    list(lambda = lambda, r = r, p = 1 - exp(-lambda * r))
  }
}

# predict(): lambda (abundance, default) or r (detection) at the fitted or new
# design, on the response scale. Visit-varying detection is returned at the fitted
# grid; new visit-level detection covariates are not supported here (use fitted()).
.tobs_predict_royle_nichols <- function(object, newdata = NULL,
                                        type = c("abundance", "detection")) {
  type  <- match.arg(type)
  model <- object$model
  if (identical(type, "detection")) {
    if (isTRUE(model$visit_varying)) {
      if (!is.null(newdata)) {
        stop("royle_nichols() detection prediction at new visit-level ",
             "covariates is not supported; use fitted() for the fitted grid.",
             call. = FALSE)
      }
      return(.rn_r_grid(object))
    }
    off <- model$process_info[[1L]]$p
    X   <- if (is.null(newdata)) model$X_processes[[2L]]
           else stats::model.matrix(model$formulas$r, newdata)
    return(stats::plogis(as.vector(X %*%
             object$means[off + seq_len(model$process_info[[2L]]$p)])))
  }
  X   <- if (is.null(newdata)) model$X_processes[[1L]]
         else stats::model.matrix(model$formulas$lambda, newdata)
  as.vector(exp(X %*% object$means[seq_len(model$process_info[[1L]]$p)]))
}

# residuals(): per-site response / pearson / deviance on the observed per-site
# detection frequency (k / n) against the marginal per-visit detection prob p
# (visit-averaged when detection varies by visit).
.tobs_residuals_royle_nichols <- function(object, type) {
  model <- object$model
  fv    <- .tobs_fitted_royle_nichols(object)
  nn    <- model$n_site; kk <- model$k_site
  obs   <- ifelse(nn > 0L, kk / nn, NA_real_)
  p     <- if (isTRUE(model$visit_varying)) rowMeans(fv$p, na.rm = TRUE) else fv$p
  eps   <- 1e-10
  res <- switch(type,
    response = obs - p,
    pearson  = (obs - p) / sqrt(p * (1 - p) / pmax(nn, 1) + eps),
    deviance = sign(obs - p) * sqrt(2 * pmax(nn, 1) * abs(
      ifelse(obs > 0, obs * log(pmax(obs, eps) / p), 0) +
      ifelse(obs < 1, (1 - obs) * log((1 - obs + eps) / (1 - p + eps)), 0))))
  list(occ = res, det = NULL)
}

# Posterior replicate detection histories: draw a coefficient vector, then per
# site N ~ Poisson(lambda) and y_ij ~ Bernoulli(1 - (1 - r_ij)^N) at the observed
# visit pattern. Returns an n_sites x max_visits matrix (or a list for nsim > 1).
.tobs_simulate_royle_nichols <- function(object, nsim = 1) {
  model <- object$model
  p_lam <- model$process_info[[1L]]$p; p_rs <- model$process_info[[2L]]$p
  Xl <- model$X_processes[[1L]]; Xr <- model$X_processes[[2L]]
  valid <- !is.na(model$y)
  draw_one <- function() {
    idx    <- sample.int(nrow(object$draws), 1L)
    beta   <- object$draws[idx, ]
    lambda <- exp(as.vector(Xl %*% beta[seq_len(p_lam)]))
    Ni     <- stats::rpois(model$n_sites, lambda)
    if (isTRUE(model$visit_varying)) {
      eta_s <- as.vector(Xr %*% beta[p_lam + seq_len(p_rs)])
      eta_c <- eta_s[model$site_of_cell] +
        as.vector(model$X_r_visit %*%
                    beta[p_lam + p_rs + seq_len(model$process_info[[2L]]$p_visit)])
      r_mat <- matrix(stats::plogis(eta_c), model$n_sites, model$max_visits,
                      byrow = TRUE)
      p_mat <- 1 - (1 - r_mat)^Ni
      yy    <- matrix(stats::rbinom(length(p_mat), 1L, as.vector(p_mat)),
                      model$n_sites, model$max_visits)
    } else {
      r   <- stats::plogis(as.vector(Xr %*% beta[p_lam + seq_len(p_rs)]))
      p_i <- 1 - (1 - r)^Ni
      yy  <- matrix(0L, model$n_sites, model$max_visits)
      for (i in seq_len(model$n_sites)) {
        yy[i, ] <- stats::rbinom(model$max_visits, 1L, p_i[i])
      }
    }
    yy[!valid] <- NA_integer_
    yy
  }
  if (nsim == 1L) return(draw_one())
  lapply(seq_len(nsim), function(s) draw_one())
}

# Pointwise log-likelihood [n_draws x n_sites] over the posterior draws, the
# exact per-site marginal (WAIC / LOO). Reuses the site / visit marginal.
.tobs_ploglik_royle_nichols <- function(object, n.draws = 1000L, n.threads = 1L) {
  model <- object$model
  draws <- object$draws
  if (!is.null(n.draws) && as.integer(n.draws) < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  p_lam <- model$process_info[[1L]]$p; p_rs <- model$process_info[[2L]]$p
  Xl <- model$X_processes[[1L]]; Xr <- model$X_processes[[2L]]
  kk <- model$k_site; nn <- model$n_site; K <- model$K_max
  t(vapply(seq_len(nrow(draws)), function(d) {
    lambda <- exp(as.vector(Xl %*% draws[d, seq_len(p_lam)]))
    beta_rs <- draws[d, p_lam + seq_len(p_rs)]
    if (isTRUE(model$visit_varying)) {
      beta_rv <- draws[d, p_lam + p_rs + seq_len(model$process_info[[2L]]$p_visit)]
      r_long  <- .rn_r_long(model, beta_rs, beta_rv)
      .rn_visit_loglik(lambda, r_long, model$vis_y, model$vis_site_idx, K,
                       model$n_sites)
    } else {
      r <- stats::plogis(as.vector(Xr %*% beta_rs))
      r <- pmin(pmax(r, 1e-8), 1 - 1e-8)
      .rn_site_loglik(lambda, r, kk, nn, K)
    }
  }, numeric(model$n_sites)))
}
