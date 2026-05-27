# =============================================================================
# ms_abun.R — community / multispecies N-mixture (spAbundance msNMix)
#
# Per-species N-mixture with Gaussian community hyperpriors on the per-species
# abundance and detection coefficients:
#
#   N_{s,i}        ~ Poisson(lambda_{s,i})            (or NegBin(lambda, r))
#   y_{s,i,j} | N  ~ Binomial(N_{s,i}, p_{s,i,j})
#   log lambda_{s,i} = X_lambda_i . (mu_lambda + b_lambda_s)
#   logit p_{s,i,j}  = X_p_{ij}   . (mu_p      + b_p_s)
#   b_lambda_s ~ N(0, Sigma_lambda),  b_p_s ~ N(0, Sigma_p)     (community RE)
#
# The latent N_{s,i} integrates out per species-site in closed form; the
# per-species coefficient deviations are the random effects, integrated by a
# C++ Laplace-EM in tulpa (`tulpa_nmix_laplace_re`, the fast path for
# gcol33/tulpa#31). This file owns only the family wiring: the data binder, the
# long-form marshalling into tulpa's fitter, and the `tobs_fit` wrapper. tulpa
# owns the marginal math, the per-species mode-finding, the EM covariance
# update, and the marginal observed-information SEs.
#
# Poisson first cut. (A global negative-binomial size is a planned tulpa
# extension; the family already carries the `mixture` flag.)
#
#   .tobs_build_ms_abun()    data binder -> model_type = "ms_nmix"
#   .tobs_ms_nmix_longform() 3D y -> stacked (y, site_idx, species_idx, X_p)
#   .tobs_fit_ms_nmix()      -> tulpa::tulpa_nmix_laplace_re()
#   build_ms_nmix_fit()      wrap into a tobs_fit
#   simulate_ms_abun()       community N-mixture simulator
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind a community N-mixture model. `y` is a 3D array
# [n_sites x max_visits x n_species] or a (named) list of n_sites x max_visits
# count matrices, one per species (mirrors .tobs_build_community). The abundance
# design X_lambda and the site-level detection design X_det_site are shared
# across species (community covariates); the counts differ per species.
.tobs_build_ms_abun <- function(abund_formula, det_formula, data, y, species,
                                det_visit_formula = NULL, det_visit_data = NULL) {
  if (is.list(y) && !is.array(y)) {
    n_species  <- length(y)
    n_sites    <- nrow(y[[1]])
    max_visits <- ncol(y[[1]])
    species_names <- if (is.character(species)) species
                     else if (!is.null(names(y))) names(y)
                     else paste0("sp", seq_len(n_species))
    y_array <- array(NA_integer_, dim = c(n_sites, max_visits, n_species))
    for (s in seq_len(n_species)) y_array[, , s] <- as.integer(round(y[[s]]))
    y <- y_array
  } else {
    if (length(dim(y)) != 3L) {
      stop("y must be a 3D array [n_sites x max_visits x n_species] or a list ",
           "of count matrices.", call. = FALSE)
    }
    species_names <- if (is.character(species)) species
                     else paste0("sp", seq_len(dim(y)[3]))
    storage.mode(y) <- "integer"
  }

  n_sites    <- dim(y)[1]
  max_visits <- dim(y)[2]
  n_species  <- dim(y)[3]
  if (nrow(data) != n_sites) {
    stop(sprintf("y has %d sites but data has %d rows", n_sites, nrow(data)),
         call. = FALSE)
  }
  if (length(species_names) != n_species) {
    stop(sprintf("species has %d names but y has %d species",
                 length(species_names), n_species), call. = FALSE)
  }
  if (any(y < 0L, na.rm = TRUE)) {
    stop("y must contain nonnegative integer counts (or NA).", call. = FALSE)
  }

  bind       <- .tobs_bind_formulas(list(lambda = abund_formula, p = det_formula),
                                    data)
  X_lambda   <- model.matrix(bind$fe$lambda, data)
  X_det_site <- model.matrix(bind$fe$p, data)

  X_det_visit <- NULL
  if (!is.null(det_visit_formula) && !is.null(det_visit_data)) {
    mf <- stats::model.frame(det_visit_formula, det_visit_data,
                             na.action = stats::na.pass)
    X_det_visit <- stats::model.matrix(det_visit_formula, mf)
    X_det_visit[is.na(X_det_visit)] <- 0
    expected_rows <- n_sites * max_visits
    if (nrow(X_det_visit) != expected_rows) {
      stop(sprintf("det_visit_data must have %d rows (n_sites * max_visits), got %d",
                   expected_rows, nrow(X_det_visit)), call. = FALSE)
    }
  }
  det_coef_names <- colnames(X_det_site)
  if (!is.null(X_det_visit)) det_coef_names <- c(det_coef_names, colnames(X_det_visit))

  structure(list(
    model_type   = "ms_nmix",
    y            = y,
    X_processes  = list(X_lambda, X_det_site),
    X_det_visit  = X_det_visit,
    formulas     = list(lambda = bind$fe$lambda, det = bind$fe$p),
    structured_terms = bind$terms,
    data         = data,
    n_sites      = n_sites,
    max_visits   = max_visits,
    n_species    = n_species,
    species_names = species_names,
    process_info = list(
      list(name = "lambda", p = ncol(X_lambda), coef_names = colnames(X_lambda),
           link = "log"),
      list(name = "p",      p = length(det_coef_names), coef_names = det_coef_names,
           link = "logit")
    )
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Long-form marshalling
# ---------------------------------------------------------------------------

# Flatten the 3D count array into the stacked long form tulpa's community
# fitter consumes: one row per observed (species, site, visit). Each species'
# NA visits drop out (site-major order), so a species missing a visit simply
# contributes fewer detection rows. The detection design is the shared site-
# level design replicated over a site's visits, with any visit-level columns
# stacked on.
.tobs_ms_nmix_longform <- function(model) {
  y          <- model$y
  X_det_site <- model$X_processes[[2]]
  X_det_visit <- model$X_det_visit
  n_sites    <- model$n_sites
  max_visits <- model$max_visits
  n_species  <- model$n_species

  site_mat  <- matrix(seq_len(n_sites),  n_sites, max_visits)
  visit_mat <- matrix(seq_len(max_visits), n_sites, max_visits, byrow = TRUE)

  y_long <- integer(0); site_idx <- integer(0); species_idx <- integer(0)
  Xp_list <- vector("list", n_species)
  for (s in seq_len(n_species)) {
    ys      <- y[, , s]
    valid_t <- as.vector(t(!is.na(ys)))                     # site-major mask
    yl      <- as.vector(t(ys))[valid_t]
    si      <- as.vector(t(site_mat))[valid_t]
    vi      <- as.vector(t(visit_mat))[valid_t]
    if (length(yl) == 0L) {
      stop(sprintf("species '%s' has no observed visits (all NA).",
                   model$species_names[s]), call. = FALSE)
    }
    Xp_s <- X_det_site[si, , drop = FALSE]
    if (!is.null(X_det_visit)) {
      Xp_s <- cbind(Xp_s, X_det_visit[(si - 1L) * max_visits + vi, , drop = FALSE])
    }
    y_long      <- c(y_long, as.integer(yl))
    site_idx    <- c(site_idx, as.integer(si))
    species_idx <- c(species_idx, rep.int(s, length(yl)))
    Xp_list[[s]] <- Xp_s
  }
  list(y = y_long, site_idx = site_idx, species_idx = species_idx,
       X_p = do.call(rbind, Xp_list))
}


# ---------------------------------------------------------------------------
# Fitter: delegate to tulpa's C++ community N-mixture Laplace-EM
# ---------------------------------------------------------------------------

.tobs_fit_ms_nmix <- function(model, mixture = "poisson", K_max = NULL,
                              max_iter = 100L, optimizer = "em",
                              n_quad = 1L, lkj_eta = 1, verbose = TRUE) {
  if (!identical(mixture, "poisson")) {
    stop("Community N-mixture currently supports mixture = \"poisson\" only ",
         "(a global negative-binomial size is a planned tulpa extension).",
         call. = FALSE)
  }
  lf  <- .tobs_ms_nmix_longform(model)
  raw <- tulpa::tulpa_nmix_laplace_re(
    y = lf$y, site_idx = lf$site_idx, species_idx = lf$species_idx,
    X_lambda = model$X_processes[[1]], X_p = lf$X_p,
    n_sites = model$n_sites, n_species = model$n_species,
    K_max = K_max, max_iter = as.integer(max_iter),
    optimizer = optimizer, n_quad = as.integer(n_quad), lkj_eta = lkj_eta,
    verbose = isTRUE(verbose))
  build_ms_nmix_fit(raw, model, mixture = mixture)
}


# ---------------------------------------------------------------------------
# Wrap a tulpa community N-mixture fit into a tobs_fit
# ---------------------------------------------------------------------------

# `means` / `vcov` are the community means (mu_lambda, mu_p) and their marginal
# covariance. The community covariances Sigma_lambda / Sigma_p and the per-
# species coefficients (mu + BLUP deviation) are carried as N-mixture community
# structure for ranef() / coef() / simulate().
build_ms_nmix_fit <- function(raw, model, mixture = "poisson") {
  pi_list <- model$process_info
  p_lam   <- pi_list[[1]]$p
  p_p     <- pi_list[[2]]$p
  lam_nm  <- pi_list[[1]]$coef_names
  p_nm    <- pi_list[[2]]$coef_names
  nms     <- c(paste0("lambda_", lam_nm), paste0("p_", p_nm))

  means <- c(as.numeric(raw$mu_lambda), as.numeric(raw$mu_p))
  names(means) <- nms
  vcov  <- as.matrix(raw$vcov)
  rownames(vcov) <- colnames(vcov) <- nms
  sds   <- sqrt(pmax(diag(vcov), 0)); names(sds) <- nms

  Sigma_lambda <- as.matrix(raw$Sigma_lambda)
  Sigma_p      <- as.matrix(raw$Sigma_p)
  dimnames(Sigma_lambda) <- list(lam_nm, lam_nm)
  dimnames(Sigma_p)      <- list(p_nm, p_nm)

  # Per-species coefficients = community mean + BLUP deviation.
  blup_lambda <- as.matrix(raw$b_lambda)
  blup_p      <- as.matrix(raw$b_p)
  coef_lambda <- sweep(blup_lambda, 2, as.numeric(raw$mu_lambda), "+")
  coef_p      <- sweep(blup_p,      2, as.numeric(raw$mu_p), "+")
  rownames(coef_lambda) <- rownames(coef_p) <- model$species_names
  rownames(blup_lambda) <- rownames(blup_p) <- model$species_names
  colnames(coef_lambda) <- colnames(blup_lambda) <- lam_nm
  colnames(coef_p)      <- colnames(blup_p)      <- p_nm

  n_pseudo <- 1000L
  draws <- .rmvn(n_pseudo, means, vcov)
  colnames(draws) <- nms

  structure(list(
    draws = draws, means = means, sds = sds, vcov = vcov,
    n_samples = n_pseudo, n_params = length(means),
    col_names = nms, param_names = nms,
    process_info = pi_list,
    model = model, spatial = NULL,
    method = "laplace",
    mixture = mixture,
    log_lik = raw$log_lik %||% NA_real_,
    ms_community = list(
      Sigma_lambda = Sigma_lambda, Sigma_p = Sigma_p,
      sd_lambda = sqrt(pmax(diag(Sigma_lambda), 0)),
      sd_p      = sqrt(pmax(diag(Sigma_p), 0)),
      coef_lambda = coef_lambda, coef_p = coef_p,
      blup_lambda = blup_lambda, blup_p = blup_p,
      optimizer = raw$optimizer %||% "em",
      n_quad = raw$n_quad %||% 1L, lkj_eta = raw$lkj_eta %||% 1
    ),
    convergence = list(converged = isTRUE(raw$converged),
                       n_iter = raw$n_iter %||% NA_integer_)
  ), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# S3 helpers (routed to from methods.R by model_type == "ms_nmix")
# ---------------------------------------------------------------------------

# ranef(): per-species BLUP deviations (mu + b_s gives the per-species
# coefficients in `coef()`; here we report the deviations b_s). Long form:
# one row per (species, arm, term).
.tobs_ranef_ms_nmix <- function(object) {
  cm <- object$ms_community
  to_long <- function(B, arm) {
    sp <- rownames(B); tm <- colnames(B)
    data.frame(
      species  = rep(sp, times = ncol(B)),
      arm      = arm,
      term     = rep(tm, each = nrow(B)),
      estimate = as.numeric(B),
      stringsAsFactors = FALSE)
  }
  out <- rbind(to_long(cm$blup_lambda, "lambda"), to_long(cm$blup_p, "p"))
  rownames(out) <- NULL
  out
}

# Per-species linear predictors at the posterior-mean (mu + BLUP) coefficients:
# site-level expected abundance lambda_{s,i} and site-level detection p_{s,i},
# each an [n_sites x n_species] matrix.
.tobs_fitted_ms_nmix <- function(object) {
  model <- object$model
  cm    <- object$ms_community
  X_lambda   <- model$X_processes[[1]]
  X_det_site <- model$X_processes[[2]]
  p_p_site   <- ncol(X_det_site)
  lambda <- exp(X_lambda %*% t(cm$coef_lambda))
  # Detection arm: use only the site-level columns of the coefficient vector
  # (visit-level detection covariates, if any, are not site-summarised here).
  p <- plogis(X_det_site %*% t(cm$coef_p[, seq_len(p_p_site), drop = FALSE]))
  dimnames(lambda) <- dimnames(p) <- list(NULL, model$species_names)
  list(lambda = lambda, p = p)
}

# simulate(): draw counts under the fitted per-species coefficients. Each draw
# samples N_{s,i} ~ Poisson(lambda_{s,i}) then y ~ Binomial(N, p) at the
# observed visit pattern, returning a 3D array matching the input y.
.tobs_simulate_ms_nmix <- function(object, nsim = 1) {
  model <- object$model
  fit   <- .tobs_fitted_ms_nmix(object)
  n_sites <- model$n_sites; max_visits <- model$max_visits
  n_species <- model$n_species
  obs_mask <- !is.na(model$y)
  draws <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    y_sim <- array(NA_integer_, dim = c(n_sites, max_visits, n_species),
                   dimnames = list(NULL, NULL, model$species_names))
    for (sp in seq_len(n_species)) {
      N <- stats::rpois(n_sites, fit$lambda[, sp])
      for (i in seq_len(n_sites)) {
        vis <- which(obs_mask[i, , sp])
        if (length(vis))
          y_sim[i, vis, sp] <- stats::rbinom(length(vis), N[i], fit$p[i, sp])
      }
    }
    draws[[s]] <- y_sim
  }
  if (nsim == 1L) draws[[1]] else draws
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate community (multispecies) N-mixture abundance data
#'
#' Per-species Royle (2004) N-mixture with Gaussian community hyperpriors:
#' `beta_lambda_s ~ N(mu_lambda, Sigma_lambda)`,
#' `beta_p_s ~ N(mu_p, Sigma_p)`, then `N_{s,i} ~ Poisson(lambda_{s,i})`
#' (or `NegBin(mu = lambda, size)`) and `y_{s,i,j} ~ Binomial(N_{s,i}, p_{s,i,j})`.
#' The returned `y` is a 3D array `[n_sites x J x n_species]` suitable for
#' [tobs()] with [ms_abun()].
#'
#' @param n_species Number of species (default 12).
#' @param N Number of sites (default 80).
#' @param J Number of replicate visits (default 4).
#' @param n_abund_covs Number of abundance covariates (default 1).
#' @param n_det_covs Number of detection covariates (default 1).
#' @param mu_lambda Community-mean abundance coefficients `c(intercept,
#'   slopes...)` on the log scale. Default `c(log(3), rep(0.4, n_abund_covs))`.
#' @param mu_p Community-mean detection coefficients on the logit scale.
#'   Default `c(0.3, rep(-0.3, n_det_covs))`.
#' @param sd_lambda Per-coefficient community SD for the abundance arm (the
#'   sqrt-diagonal of `Sigma_lambda`). Length 1 (recycled) or `1 + n_abund_covs`.
#'   Default 0.5.
#' @param sd_p Per-coefficient community SD for the detection arm. Default 0.4.
#' @param mixture Abundance mixing distribution: `"poisson"` (default) or
#'   `"negbin"`.
#' @param size Negative-binomial size `r` (ignored under Poisson). Default 3.
#' @param seed Optional random seed.
#' @return A list with `y` (3D count array), `data` (site covariate frame),
#'   `species` (species names), and `truth` (community means / SDs and the
#'   per-species coefficients, lambda, p, latent N).
#' @export
simulate_ms_abun <- function(n_species = 12, N = 80, J = 4,
                             n_abund_covs = 1, n_det_covs = 1,
                             mu_lambda = NULL, mu_p = NULL,
                             sd_lambda = 0.5, sd_p = 0.4,
                             mixture = c("poisson", "negbin"), size = 3,
                             seed = NULL) {
  mixture <- match.arg(mixture)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(mu_lambda)) mu_lambda <- c(log(3), rep(0.4, n_abund_covs))
  if (is.null(mu_p))      mu_p      <- c(0.3, rep(-0.3, n_det_covs))
  p_lam <- length(mu_lambda); p_p <- length(mu_p)
  sd_lambda <- if (length(sd_lambda) == 1L) rep(sd_lambda, p_lam) else sd_lambda
  sd_p      <- if (length(sd_p) == 1L)      rep(sd_p, p_p)        else sd_p

  abund_covs <- data.frame(matrix(stats::rnorm(N * n_abund_covs), N, n_abund_covs))
  names(abund_covs) <- paste0("abund_cov", seq_len(n_abund_covs))
  det_covs <- data.frame(matrix(stats::rnorm(N * n_det_covs), N, n_det_covs))
  names(det_covs) <- paste0("det_cov", seq_len(n_det_covs))
  data <- cbind(abund_covs, det_covs)
  X_lambda <- stats::model.matrix(~ ., abund_covs)
  X_det    <- stats::model.matrix(~ ., det_covs)

  beta_lambda <- matrix(stats::rnorm(n_species * p_lam, 0, rep(sd_lambda, each = n_species)),
                        n_species, p_lam) + matrix(mu_lambda, n_species, p_lam, byrow = TRUE)
  beta_p <- matrix(stats::rnorm(n_species * p_p, 0, rep(sd_p, each = n_species)),
                   n_species, p_p) + matrix(mu_p, n_species, p_p, byrow = TRUE)

  species_names <- paste0("sp", seq_len(n_species))
  y <- array(NA_integer_, dim = c(N, J, n_species),
             dimnames = list(NULL, NULL, species_names))
  lambda <- matrix(NA_real_, N, n_species); p_arr <- matrix(NA_real_, N, n_species)
  Nlat   <- matrix(NA_integer_, N, n_species)
  for (s in seq_len(n_species)) {
    lam <- exp(as.vector(X_lambda %*% beta_lambda[s, ]))
    pp  <- plogis(as.vector(X_det %*% beta_p[s, ]))
    Ns  <- if (identical(mixture, "negbin")) stats::rnbinom(N, size = size, mu = lam)
           else stats::rpois(N, lam)
    for (i in seq_len(N)) y[i, , s] <- stats::rbinom(J, Ns[i], pp[i])
    lambda[, s] <- lam; p_arr[, s] <- pp; Nlat[, s] <- Ns
  }

  list(
    y = y, data = data, species = species_names,
    truth = list(
      mu_lambda = mu_lambda, mu_p = mu_p,
      sd_lambda = sd_lambda, sd_p = sd_p,
      beta_lambda = beta_lambda, beta_p = beta_p,
      lambda = lambda, p = p_arr, N = Nlat,
      mixture = mixture, size = if (identical(mixture, "negbin")) size else NA_real_)
  )
}
