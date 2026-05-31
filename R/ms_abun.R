# =============================================================================
# ms_abun.R — community / multispecies N-mixture (spAbundance msNMix)
#
# Per-species N-mixture with Gaussian community hyperpriors on the per-species
# abundance and detection coefficients:
#
#   N_{s,i}        ~ Poisson(lambda_{s,i})            (or NegBin(lambda, r_s))
#   y_{s,i,j} | N  ~ Binomial(N_{s,i}, p_{s,i,j})
#   log lambda_{s,i} = X_lambda_i . (mu_lambda + b_lambda_s)
#   logit p_{s,i,j}  = X_p_{ij}   . (mu_p      + b_p_s)
#   b_lambda_s ~ N(0, Sigma_lambda),  b_p_s ~ N(0, Sigma_p)     (community RE)
#
# The latent N_{s,i} integrates out per species-site in closed form; the
# per-species coefficient deviations are the random effects, integrated by a
# C++ Laplace-EM in `nmix_laplace_re()`. The community fitter consumes tulpa's
# generic AGHQ RE engine (`tulpa::tulpa_re_aghq()`) through an
# `NMixCommunityOracle` (an `XPtr<tulpa::REGroupOracle>` constructed in
# tulpaObs's src/), so the per-species marginal / score / observed-info
# assembly runs entirely in tulpaObs while the structure-agnostic integration,
# log-Cholesky parametrization and LKJ penalty come from the engine. This file
# owns the family wiring: the data binder, the long-form marshalling into the
# fitter, and the `tobs_fit` wrapper.
#
# Poisson and negative-binomial abundance. Under NB the dispersion is a
# per-species random effect log_r_s ~ N(mu_log_r, sigma_log_r); the oracle widens
# the per-species RE vector with a trailing log_r_s coordinate and the community
# log-dispersion (mu_log_r, sigma_log_r) joins the community hyperparameters.
#
#   .tobs_build_ms_abun()    data binder -> model_type = "ms_nmix"
#   .tobs_ms_nmix_longform() 3D y -> stacked (y, site_idx, species_idx, X_p)
#   .tobs_fit_ms_nmix()      -> nmix_laplace_re()
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
  if (!identical(mixture, "poisson") && !identical(mixture, "negbin")) {
    stop("Community N-mixture supports mixture = \"poisson\" or \"negbin\" ",
         "(got \"", mixture, "\").", call. = FALSE)
  }
  # tulpaObs vocabulary ("poisson" / "negbin") -> nmix_laplace_re's
  # mixing-distribution code ("P" / "NB"). NB makes the dispersion a per-species
  # random effect log_r_s ~ N(mu_log_r, sigma_log_r): mu_log_r joins the community
  # means as a fixed effect and b_logr_s is the trailing per-species RE coordinate,
  # joint-optimised via the analytic-gradient AGHQ path. There is no closed-form EM
  # for NB, so when the user does not pin them, switch the EM defaults to the
  # joint_grad / n_quad = 5 NB defaults (matching nmix_laplace_re()'s own
  # missing()-driven NB defaults).
  mix_code <- if (identical(mixture, "negbin")) "NB" else "P"
  if (mix_code == "NB" && identical(optimizer, "em")) optimizer <- "joint_grad"
  if (mix_code == "NB" && n_quad == 1L)               n_quad    <- 5L
  lf  <- .tobs_ms_nmix_longform(model)
  raw <- nmix_laplace_re(
    y = lf$y, site_idx = lf$site_idx, species_idx = lf$species_idx,
    X_lambda = model$X_processes[[1]], X_p = lf$X_p,
    n_sites = model$n_sites, n_species = model$n_species,
    K_max = K_max, max_iter = as.integer(max_iter),
    mixture = mix_code,
    optimizer = optimizer, n_quad = as.integer(n_quad), lkj_eta = lkj_eta,
    verbose = isTRUE(verbose))
  build_ms_nmix_fit(raw, model, mixture = mixture)
}


# ---------------------------------------------------------------------------
# Spatial fitter: a shared ICAR / BYM2 / proper-CAR field on the abundance arm
# ---------------------------------------------------------------------------

# Fit the spatial community N-mixture (the sfMsNMix analogue, gcol33/tulpaObs#12).
# One spatial unit per site (identity map, as for single-species abun()).
# Dispatches on the abundance-formula spatial term type to the matching nested
# Laplace-EM wrapper in nmix_laplace_re_spatial.R. The NB size r is integrated
# over the outer grid (field-agnostic), so it falls out of the existing r-grid
# plumbing as a global grid hyperparameter (not the per-species log_r RE the
# non-spatial NB path uses).
.tobs_fit_ms_nmix_spatial <- function(model, spatial, mixture = "poisson",
                                      K_max = NULL, max_iter = 100L,
                                      verbose = TRUE) {
  .tobs_reject_weighted_spatial(spatial, "ms_abun() spatial")
  mix_code <- if (identical(mixture, "negbin")) "NB" else "P"
  n_sites <- model$n_sites
  # Continuous Matern field shared across species: thread the SPDE FEM precision
  # Q(range, sigma) and the mesh projection A through the same community
  # Laplace-EM driver, broadcasting the site-level field across species.
  if (identical(spatial$type, "spde")) {
    lf  <- .tobs_ms_nmix_longform(model)
    raw <- nmix_community_laplace_spde(
      lf = lf, X_lambda = model$X_processes[[1]],
      n_sites = n_sites, n_species = model$n_species,
      spatial = spatial, mixture = mix_code, K_max = K_max,
      max_iter = max_iter, verbose = isTRUE(verbose))
    return(build_ms_nmix_fit(raw, model, mixture = mixture, spatial = spatial))
  }
  if (identical(spatial$type, "gp") || identical(spatial$type, "multiscale_gp")) {
    stop(sprintf(
      "ms_abun() spatial supports the continuous field via spde(); the dense ",
      "GP term '%s' is not wired. Use spde() for a mesh-based continuous Matern ",
      "field, or icar() / bym2() / car_proper() for an areal field.",
      spatial$type), call. = FALSE)
  }
  if (!spatial$type %in% c("icar", "bym2", "car_proper")) {
    stop(sprintf(
      "ms_abun() spatial supports the areal terms icar() / bym2() / car_proper() ",
      "and the continuous mesh field spde() under method = \"nested_laplace\"; ",
      "got '%s'. (car() is the improper non-intrinsic CAR; use icar() for the ",
      "intrinsic field.)"), spatial$type, call. = FALSE)
  }
  if ((spatial$n_units %||% n_sites) != n_sites) {
    stop(sprintf("spatial term has %d units but the model has %d sites; one ",
                 "spatial unit per site is required for the community N-mixture.",
                 spatial$n_units, n_sites), call. = FALSE)
  }
  lf  <- .tobs_ms_nmix_longform(model)
  csr <- .nmix_spatial_csr(spatial)
  common <- list(lf = lf, X_lambda = model$X_processes[[1]],
                 n_sites = n_sites, n_species = model$n_species,
                 csr = csr, n_spatial = spatial$n_units %||% n_sites,
                 mixture = mix_code, K_max = K_max,
                 max_iter = as.integer(max_iter), verbose = isTRUE(verbose))
  raw <- switch(
    spatial$type,
    icar = do.call(nmix_community_laplace_icar, common),
    car_proper = do.call(nmix_community_laplace_car_proper,
                         c(common, list(graph = spatial$graph))),
    bym2 = do.call(nmix_community_laplace_bym2,
                   c(common, list(scale_factor = spatial$scale_factor %||%
                                    compute_bym2_scale(spatial$graph))))
  )
  build_ms_nmix_fit(raw, model, mixture = mixture, spatial = spatial)
}


# ---------------------------------------------------------------------------
# Wrap a tulpa community N-mixture fit into a tobs_fit
# ---------------------------------------------------------------------------

# `means` / `vcov` are the community means (mu_lambda, mu_p) and their marginal
# covariance. The community covariances Sigma_lambda / Sigma_p and the per-
# species coefficients (mu + BLUP deviation) are carried as N-mixture community
# structure for ranef() / coef() / simulate().
build_ms_nmix_fit <- function(raw, model, mixture = "poisson", spatial = NULL) {
  pi_list <- model$process_info
  p_lam   <- pi_list[[1]]$p
  p_p     <- pi_list[[2]]$p
  lam_nm  <- pi_list[[1]]$coef_names
  p_nm    <- pi_list[[2]]$coef_names
  nms     <- c(paste0("lambda_", lam_nm), paste0("p_", p_nm))

  # NB makes the dispersion a per-species random effect log_r_s ~ N(mu_log_r,
  # sigma_log_r): the community log-dispersion mu_log_r is the trailing theta
  # coordinate (the joint optimizer estimates it alongside the community means),
  # so append "log_r" to the coefficient surface (means / vcov / sds / draws) and
  # coef() / vcov() / confint() see mu_log_r with its marginal SE. The per-species
  # sizes r_s and the community covariance term sigma_log_r are summarized on
  # `ms_dispersion`.
  # On the SPATIAL path the NB size r is integrated over the outer grid (it is
  # field-agnostic), so it carries no log_r coordinate and is summarized as a
  # grid hyperparameter (raw$dispersion / ms_hyper) rather than the per-species
  # log_r RE the non-spatial NB path appends here.
  is_nb <- is.null(spatial) && identical(mixture, "negbin") &&
           !is.null(raw$mu_log_r) && is.finite(raw$mu_log_r)
  log_r <- if (is_nb) unname(as.numeric(raw$mu_log_r)) else NA_real_

  means <- c(as.numeric(raw$mu_lambda), as.numeric(raw$mu_p))
  if (is_nb) {
    means <- c(means, log_r)
    nms   <- c(nms, "log_r")
  }
  names(means) <- nms
  vcov  <- as.matrix(raw$vcov)
  if (nrow(vcov) != length(nms) || ncol(vcov) != length(nms)) {
    stop(sprintf("build_ms_nmix_fit(): vcov dim %dx%d does not match the %d ",
                 nrow(vcov), ncol(vcov), length(nms)),
         "coefficient names. Expected the joint optimizer to return one ",
         "trailing log_r row/col under NB and a plain (p_lambda + p_p) block ",
         "under Poisson.", call. = FALSE)
  }
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

  ms_dispersion <- if (is_nb) {
    se_logr <- unname(sds["log_r"])
    r_s <- as.numeric(raw$r_s)
    names(r_s) <- model$species_names
    list(r           = exp(log_r),          # community-mean size (LogNormal median)
         log_r       = log_r,               # == mu_log_r
         mu_log_r    = log_r,
         sigma_log_r = unname(as.numeric(raw$sigma_log_r)),
         r_s         = r_s,                 # per-species sizes exp(mu_log_r + b_logr_s)
         r_sd        = exp(log_r) * se_logr)
  } else if (!is.null(spatial)) {
    # Spatial NB: r is grid-integrated, summarized in raw$dispersion.
    raw$dispersion
  } else NULL

  structure(list(
    draws = draws, means = means, sds = sds, vcov = vcov,
    n_samples = n_pseudo, n_params = length(means),
    col_names = nms, param_names = nms,
    n_fixed = length(means), fixed_names = nms,
    process_info = pi_list,
    model = model, spatial = spatial,
    spatial_field = raw$spatial_field,
    ms_hyper = raw$hyper,
    method = if (is.null(spatial)) "laplace" else "nested_laplace",
    mixture = mixture,
    log_lik = raw$log_lik %||% NA_real_,
    ms_community = list(
      Sigma_lambda = Sigma_lambda, Sigma_p = Sigma_p,
      sd_lambda = sqrt(pmax(diag(Sigma_lambda), 0)),
      sd_p      = sqrt(pmax(diag(Sigma_p), 0)),
      coef_lambda = coef_lambda, coef_p = coef_p,
      blup_lambda = blup_lambda, blup_p = blup_p,
      # NB per-species log-dispersion deviation b_logr_s (NULL under Poisson),
      # so ranef() carries the dispersion RE alongside the coefficient REs.
      blup_logr = if (is_nb)
        matrix(as.numeric(raw$b_logr), ncol = 1L,
               dimnames = list(model$species_names, "log_r")) else NULL,
      optimizer = raw$optimizer %||% "em",
      n_quad = raw$n_quad %||% 1L, lkj_eta = raw$lkj_eta %||% 1
    ),
    ms_dispersion = ms_dispersion,
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
  if (!is.null(cm$blup_logr)) out <- rbind(out, to_long(cm$blup_logr, "logr"))
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
  # Abundance log-linear predictor; on the spatial path add the shared field
  # offset f_i (one spatial unit per site) to every species' log lambda.
  eta_lambda <- X_lambda %*% t(cm$coef_lambda)
  if (!is.null(object$spatial_field)) {
    eta_lambda <- sweep(eta_lambda, 1, as.numeric(object$spatial_field), "+")
  }
  lambda <- exp(eta_lambda)
  # Detection arm: use only the site-level columns of the coefficient vector
  # (visit-level detection covariates, if any, are not site-summarised here).
  p <- plogis(X_det_site %*% t(cm$coef_p[, seq_len(p_p_site), drop = FALSE]))
  dimnames(lambda) <- dimnames(p) <- list(NULL, model$species_names)
  list(lambda = lambda, p = p)
}

# simulate(): draw counts under the fitted per-species coefficients. Each draw
# samples N_{s,i} ~ Poisson(lambda_{s,i}) (or NegBin(mu = lambda, size = r_s)
# under mixture = "negbin") then y ~ Binomial(N, p) at the observed visit
# pattern, returning a 3D array matching the input y. Under NB the per-species
# size is the fitted r_s (non-spatial, per-species log_r RE) or the grid-
# integrated community size r (spatial); a non-finite size is the Poisson limit.
.tobs_simulate_ms_nmix <- function(object, nsim = 1) {
  model <- object$model
  fit   <- .tobs_fitted_ms_nmix(object)
  n_sites <- model$n_sites; max_visits <- model$max_visits
  n_species <- model$n_species
  obs_mask <- !is.na(model$y)
  # Per-species NB size: r_s if present (non-spatial NB), else the community r
  # (spatial NB), recycled across species; NA -> Poisson.
  size_s <- rep(NA_real_, n_species)
  if (identical(object$mixture, "negbin") && !is.null(object$ms_dispersion)) {
    rs <- object$ms_dispersion$r_s
    size_s <- if (!is.null(rs) && length(rs) == n_species) as.numeric(rs)
              else rep(as.numeric(object$ms_dispersion$r %||% NA_real_), n_species)
  }
  draws <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    y_sim <- array(NA_integer_, dim = c(n_sites, max_visits, n_species),
                   dimnames = list(NULL, NULL, model$species_names))
    for (sp in seq_len(n_species)) {
      N <- if (is.finite(size_s[sp]))
             stats::rnbinom(n_sites, size = size_s[sp], mu = fit$lambda[, sp])
           else stats::rpois(n_sites, fit$lambda[, sp])
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
#' (or `NegBin(mu = lambda, size = r_s)` with a per-species size
#' `r_s = exp(mu_log_r + b_logr_s)`) and
#' `y_{s,i,j} ~ Binomial(N_{s,i}, p_{s,i,j})`. The returned `y` is a 3D array
#' `[n_sites x J x n_species]` suitable for [tobs()] with [ms_abun()].
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
#' @param size Community-mean negative-binomial size, equal to
#'   \eqn{\exp(\mu_{\log r})} (ignored under Poisson). Default 3. The per-species
#'   sizes are \eqn{r_s = \exp(\mu_{\log r} + b^{\log r}_s)} with
#'   \eqn{b^{\log r}_s \sim N(0, \sigma_{\log r}^2)}.
#' @param sigma_logr Standard deviation of the per-species log-dispersion random
#'   effect `log_r_s` (used only under `mixture = "negbin"`). Default 0.4. Set to
#'   0 for a shared (single-`r`) community.
#' @param graph Optional `N x N` 0/1 adjacency matrix. When supplied, a single
#'   shared spatial field `f` (one value per site) is drawn from a proper GMRF on
#'   the graph, centered and scaled to standard deviation `sigma.field`, and
#'   added to every species' abundance log-linear predictor
#'   (`log lambda_{s,i} = X_lambda_i . beta_lambda_s + f_i`). `N` is taken from
#'   `nrow(graph)`. The truth field is returned in `truth$field`.
#' @param sigma.field Standard deviation of the shared spatial field (used only
#'   when `graph` is supplied). Default 0.6.
#' @param seed Optional random seed.
#' @return A list with `y` (3D count array), `data` (site covariate frame),
#'   `species` (species names), and `truth` (community means / SDs, the
#'   per-species coefficients, lambda, p, latent N, the shared `field` when
#'   `graph` is given, and -- under `negbin` -- `mu_log_r`, `sigma_log_r`, the
#'   per-species `b_logr` / `r_s`).
#' @export
simulate_ms_abun <- function(n_species = 12, N = 80, J = 4,
                             n_abund_covs = 1, n_det_covs = 1,
                             mu_lambda = NULL, mu_p = NULL,
                             sd_lambda = 0.5, sd_p = 0.4,
                             mixture = c("poisson", "negbin"), size = 3,
                             sigma_logr = 0.4,
                             graph = NULL, sigma.field = 0.6,
                             seed = NULL) {
  mixture <- match.arg(mixture)
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(graph)) N <- nrow(graph)
  if (is.null(mu_lambda)) mu_lambda <- c(log(3), rep(0.4, n_abund_covs))
  if (is.null(mu_p))      mu_p      <- c(0.3, rep(-0.3, n_det_covs))
  p_lam <- length(mu_lambda); p_p <- length(mu_p)
  sd_lambda <- if (length(sd_lambda) == 1L) rep(sd_lambda, p_lam) else sd_lambda
  sd_p      <- if (length(sd_p) == 1L)      rep(sd_p, p_p)        else sd_p

  # n_*_covs may be 0 (an intercept-only arm). Build each covariate frame with
  # N rows and the requested number of columns, then stitch a single N-row data
  # frame; a 0-column arm contributes no columns but keeps the row count.
  make_covs <- function(n_covs, prefix) {
    if (n_covs <= 0L) return(data.frame(row.names = seq_len(N)))
    m <- matrix(stats::rnorm(N * n_covs), N, n_covs)
    df <- as.data.frame(m)
    names(df) <- paste0(prefix, seq_len(n_covs))
    df
  }
  abund_covs <- make_covs(n_abund_covs, "abund_cov")
  det_covs   <- make_covs(n_det_covs,   "det_cov")
  data <- data.frame(row.names = seq_len(N))
  if (ncol(abund_covs)) data <- cbind(data, abund_covs)
  if (ncol(det_covs))   data <- cbind(data, det_covs)
  # `~ .` needs at least one column; an intercept-only arm uses `~ 1` with an
  # N-row frame so the design is the N x 1 intercept.
  design_of <- function(df) {
    if (ncol(df)) stats::model.matrix(~ ., df)
    else stats::model.matrix(~ 1, data.frame(row.names = seq_len(N)))
  }
  X_lambda <- design_of(abund_covs)
  X_det    <- design_of(det_covs)

  beta_lambda <- matrix(stats::rnorm(n_species * p_lam, 0, rep(sd_lambda, each = n_species)),
                        n_species, p_lam) + matrix(mu_lambda, n_species, p_lam, byrow = TRUE)
  beta_p <- matrix(stats::rnorm(n_species * p_p, 0, rep(sd_p, each = n_species)),
                   n_species, p_p) + matrix(mu_p, n_species, p_p, byrow = TRUE)

  # Per-species NB size r_s = exp(mu_log_r + b_logr_s), b_logr_s ~ N(0, sigma_logr^2);
  # mu_log_r = log(size). NA-valued under Poisson.
  is_nb    <- identical(mixture, "negbin")
  mu_log_r <- if (is_nb) log(size) else NA_real_
  b_logr   <- if (is_nb) stats::rnorm(n_species, 0, sigma_logr) else rep(0, n_species)
  r_s      <- if (is_nb) exp(mu_log_r + b_logr) else rep(NA_real_, n_species)

  # Shared spatial field f ~ N(0, Q^{-1}) on the graph (Q = D - W + ridge, a
  # proper GMRF), centred and scaled to sigma.field. Zero when no graph.
  field <- rep(0, N)
  if (!is.null(graph)) {
    Q  <- diag(rowSums(graph)) - graph + diag(1e-3, N)
    Lc <- chol(Q)                              # Q = Lc' Lc
    f0 <- backsolve(Lc, stats::rnorm(N))       # ~ N(0, Q^{-1})
    f0 <- f0 - mean(f0)
    field <- as.numeric(sigma.field * f0 / stats::sd(f0))
  }

  species_names <- paste0("sp", seq_len(n_species))
  y <- array(NA_integer_, dim = c(N, J, n_species),
             dimnames = list(NULL, NULL, species_names))
  lambda <- matrix(NA_real_, N, n_species); p_arr <- matrix(NA_real_, N, n_species)
  Nlat   <- matrix(NA_integer_, N, n_species)
  for (s in seq_len(n_species)) {
    lam <- exp(as.vector(X_lambda %*% beta_lambda[s, ]) + field)
    pp  <- plogis(as.vector(X_det %*% beta_p[s, ]))
    Ns  <- if (is_nb) stats::rnbinom(N, size = r_s[s], mu = lam)
           else stats::rpois(N, lam)
    for (i in seq_len(N)) y[i, , s] <- stats::rbinom(J, Ns[i], pp[i])
    lambda[, s] <- lam; p_arr[, s] <- pp; Nlat[, s] <- Ns
  }
  names(r_s) <- names(b_logr) <- species_names

  list(
    y = y, data = data, species = species_names,
    truth = list(
      mu_lambda = mu_lambda, mu_p = mu_p,
      sd_lambda = sd_lambda, sd_p = sd_p,
      beta_lambda = beta_lambda, beta_p = beta_p,
      lambda = lambda, p = p_arr, N = Nlat,
      field = if (!is.null(graph)) field else NULL,
      mixture = mixture,
      size = if (is_nb) size else NA_real_,
      mu_log_r = mu_log_r,
      sigma_log_r = if (is_nb) sigma_logr else NA_real_,
      b_logr = b_logr,
      r_s = r_s)
  )
}
