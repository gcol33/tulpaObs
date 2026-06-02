# =============================================================================
# occu_cover.R - Joint occupancy-detection + cover hurdle family.
#
# A per-cell mixture combining MacKenzie single-season occupancy with a
# positive-cover observation on every detected visit. Plays the role of
# N-mixture, but for vegetation cover instead of counts: cell-level latent
# presence z_i ~ Bernoulli(psi_i), per-visit detection y_ij | z_i = 1
# ~ Bernoulli(p_ij), and per-visit cover y_pos_ij | y_ij = 1 ~ f_pos (beta or
# lognormal) on a third linear predictor. z marginalises out in closed form
# (two states), so the marginal log-likelihood is exact and the fit is a
# direct Laplace on the packed parameter vector
# (beta_occ, beta_p, beta_pos, log_dispersion).
#
# Per-cell likelihood:
#
#   any_det_i : L_i = psi_i * prod_j h_ij
#   no_det_i  : L_i = psi_i * prod_j (1 - p_ij) + (1 - psi_i)
#
#   h_ij = (1 - p_ij) * 1{y_ij = 0}
#        + p_ij       * f_pos(y_pos_ij; eta_pos_ij, dispersion) * 1{y_ij = 1}
#
# Reduces to occu() when f_pos is degenerate and to the plot-level cover
# hurdle when J = 1 and p = 1. v1 covers the non-spatial Laplace path; a
# shared spatial field across the occ and cover arms (the analogue of
# cover()'s nested-Laplace joint engine) is v2.
#
# Files this touches:
#   R/obs_families.R    - occu_cover(positive = ) constructor
#   R/tobs.R            - dispatch switch + .tobs_family_methods entry
#   tests/testthat/     - test-occu-cover.R recovery test
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind the joint occupancy-cover model. The psi predictor is cell-level
# (X_occ, n_sites rows). The detection and positive-cover predictors are
# visit-level (X_p_visit, X_pos_visit, n_sites * max_visits rows in
# site-major order). Visits with NA y are masked out of the likelihood; the
# binder zero-fills NA design rows so matrix algebra stays defined.
.tobs_build_occu_cover <- function(occ_formula, det_formula, pos_formula,
                                   data, y, y_pos, positive,
                                   det_visit_formula = NULL,
                                   det_visit_data   = NULL,
                                   pos_visit_formula = NULL,
                                   pos_visit_data    = NULL) {
  if (!is.matrix(y) || !is.matrix(y_pos)) {
    stop("y and y_pos must be matrices (n_sites x max_visits).", call. = FALSE)
  }
  if (!all(dim(y) == dim(y_pos))) {
    stop("y and y_pos must have identical dimensions.", call. = FALSE)
  }
  if (nrow(y) != nrow(data)) {
    stop(sprintf("y has %d rows but data has %d rows.", nrow(y), nrow(data)),
         call. = FALSE)
  }

  n_sites    <- nrow(y)
  max_visits <- ncol(y)

  # Detection arm carries 0/1 values; NA = visit not observed.
  y_int <- matrix(as.integer(y), n_sites, max_visits)
  valid <- !is.na(y_int)
  if (any(y_int[valid] != 0L & y_int[valid] != 1L)) {
    stop("y must contain only 0, 1, or NA (binary detection per visit).",
         call. = FALSE)
  }
  y_int[!valid] <- 0L

  # Cover arm: meaningful only where y == 1. Zero-fill the rest so matrix
  # algebra stays defined; the likelihood gates evaluation by `y == 1`.
  y_pos_num <- matrix(as.numeric(y_pos), n_sites, max_visits)
  pos_mask  <- valid & (y_int == 1L)
  if (any(!is.finite(y_pos_num[pos_mask]))) {
    stop("y_pos must be finite at every detected visit (y == 1).",
         call. = FALSE)
  }
  if (identical(positive, "beta")) {
    bad <- pos_mask & (y_pos_num <= 0 | y_pos_num >= 1)
    if (any(bad)) {
      stop("Beta positive arm requires 0 < y_pos < 1 at every detected ",
           "visit; clip with pmin(pmax(y_pos, eps), 1 - eps).", call. = FALSE)
    }
  } else {
    bad <- pos_mask & (y_pos_num <= 0)
    if (any(bad)) {
      stop("Lognormal positive arm requires y_pos > 0 at every detected ",
           "visit.", call. = FALSE)
    }
  }
  y_pos_num[!pos_mask] <- 0

  # Reject structured terms in v1 (spatial sharing across arms is v2).
  .occu_cover_reject_structured(occ_formula, "occupancy")
  .occu_cover_reject_structured(det_formula, "detection")
  .occu_cover_reject_structured(pos_formula, "positive cover")

  # Site-level occupancy design.
  X_occ <- stats::model.matrix(occ_formula, data)

  # Site-level detection / positive design (intercept + any site-level covs).
  # The visit-level path adds visit-varying covariates on top.
  X_det_site <- stats::model.matrix(det_formula, data)
  X_pos_site <- stats::model.matrix(pos_formula, data)

  X_det_visit <- .occu_cover_build_visit_X(det_visit_formula, det_visit_data,
                                           n_sites, max_visits, arm = "detection")
  X_pos_visit <- .occu_cover_build_visit_X(pos_visit_formula, pos_visit_data,
                                           n_sites, max_visits, arm = "positive cover")

  det_coef_names <- colnames(X_det_site)
  pos_coef_names <- colnames(X_pos_site)
  if (!is.null(X_det_visit)) det_coef_names <- c(det_coef_names, colnames(X_det_visit))
  if (!is.null(X_pos_visit)) pos_coef_names <- c(pos_coef_names, colnames(X_pos_visit))

  structure(list(
    model_type  = "occu_cover",
    positive    = positive,
    y           = y_int,
    y_pos       = y_pos_num,
    valid       = valid,
    n_sites     = n_sites,
    max_visits  = max_visits,
    X_occ       = X_occ,
    X_det_site  = X_det_site,
    X_pos_site  = X_pos_site,
    X_det_visit = X_det_visit,
    X_pos_visit = X_pos_visit,
    formulas    = list(occ = occ_formula, det = det_formula, pos = pos_formula),
    data        = data,
    process_info = list(
      list(name = "psi", p = ncol(X_occ),
           coef_names = colnames(X_occ), link = "logit"),
      list(name = "p",   p = length(det_coef_names),
           coef_names = det_coef_names, link = "logit"),
      list(name = "pos", p = length(pos_coef_names),
           coef_names = pos_coef_names,
           link = if (positive == "beta") "logit" else "identity")
    )
  ), class = "tobs_model")
}

# Reject formulas containing structured terms (bym2, icar, car, gp, etc.).
# v1 of occu_cover is non-spatial; spatial sharing across the three arms is v2.
.occu_cover_reject_structured <- function(formula, arm) {
  if (is.null(formula)) return(invisible(NULL))
  labs <- attr(stats::terms(formula), "term.labels")
  structured <- c("bym2", "icar", "car", "car_proper", "gp", "spde",
                  "multiscale_gp", "re", "temporal", "svc", "latent", "copy")
  hits <- character(0)
  for (lab in labs) {
    fn <- tryCatch(as.character(as.call(parse(text = lab)[[1]])[[1]]),
                   error = function(e) NA_character_)
    if (!is.na(fn) && fn %in% structured) hits <- c(hits, fn)
  }
  if (length(hits) > 0L) {
    stop(sprintf(paste0(
      "occu_cover() v1 does not support structured terms (%s) on the %s arm. ",
      "Shared spatial / temporal / RE fields across the three arms is v2; ",
      "for now use a plain fixed-effects formula on each."),
      paste(unique(hits), collapse = ", "), arm), call. = FALSE)
  }
  invisible(NULL)
}

# Build the visit-level design matrix from a long-form visit data frame.
# Visit-level data frame is in site-major order (row k = site (k-1) %/% J + 1,
# visit (k-1) %% J + 1). NA-pad cells with fewer than J visits with zero.
.occu_cover_build_visit_X <- function(visit_formula, visit_data,
                                       n_sites, max_visits, arm) {
  if (is.null(visit_formula) || is.null(visit_data)) return(NULL)
  expected_rows <- n_sites * max_visits
  if (nrow(visit_data) != expected_rows) {
    stop(sprintf("%s visit data must have %d rows (n_sites * max_visits), got %d",
                 arm, expected_rows, nrow(visit_data)), call. = FALSE)
  }
  mf <- stats::model.frame(visit_formula, visit_data,
                            na.action = stats::na.pass)
  X <- stats::model.matrix(visit_formula, mf)
  X[is.na(X)] <- 0
  # Drop the intercept column: it would duplicate the site-level intercept.
  int_col <- match("(Intercept)", colnames(X))
  if (!is.na(int_col)) X <- X[, -int_col, drop = FALSE]
  if (ncol(X) == 0L) return(NULL)
  X
}


# ---------------------------------------------------------------------------
# Likelihood
# ---------------------------------------------------------------------------

# Negative log-posterior at packed parameter vector
#   par = c(beta_occ, beta_p, beta_pos, log_dispersion)
# evaluated against the bound model. Gaussian prior with diagonal (mean, prec)
# aligned with par; flat prior when pprec == 0.
.tobs_occu_cover_nlp <- function(par, model, pmean, pprec) {
  pi_list <- model$process_info
  p_occ   <- pi_list[[1L]]$p
  p_p     <- pi_list[[2L]]$p
  p_pos   <- pi_list[[3L]]$p

  bo  <- par[seq_len(p_occ)]
  bp  <- par[p_occ + seq_len(p_p)]
  bpos<- par[p_occ + p_p + seq_len(p_pos)]
  log_disp <- par[length(par)]

  eta <- .occu_cover_eta_from_par(model, bo, bp, bpos)
  ll <- sum(.occu_cover_site_ll(model, eta$psi, eta$p_mat, eta$ep_mat, log_disp))
  penalty <- 0.5 * sum(pprec * (par - pmean)^2)
  -ll + penalty
}

# Build the three arm linear predictors from an arm-split coefficient triple.
# `bo` is the occupancy (psi) coefficient vector; `bp` / `bpos` are the
# detection / cover coefficient vectors, each packed as the site-level block
# followed by the optional visit-level block exactly as the fitter stacks them.
# Returns `psi` (length n_sites), `p_mat` (the per-visit detection probability,
# [n_sites x max_visits]), and `ep_mat` (the per-visit cover linear predictor on
# its link scale, [n_sites x max_visits]). The eta does not depend on the
# response, so the single-species fit and the community per-species marginal
# share one builder (single source of truth for the occu_cover predictors).
.occu_cover_eta_from_par <- function(model, bo, bp, bpos) {
  cl <- function(e) pmin(pmax(e, -30), 30)
  n_sites    <- model$n_sites
  max_visits <- model$max_visits

  psi <- stats::plogis(cl(as.numeric(model$X_occ %*% bo)))

  bp_site <- bp[seq_len(ncol(model$X_det_site))]
  bp_visit <- if (!is.null(model$X_det_visit)) {
    bp[ncol(model$X_det_site) + seq_len(ncol(model$X_det_visit))]
  } else numeric(0)
  bpos_site <- bpos[seq_len(ncol(model$X_pos_site))]
  bpos_visit <- if (!is.null(model$X_pos_visit)) {
    bpos[ncol(model$X_pos_site) + seq_len(ncol(model$X_pos_visit))]
  } else numeric(0)

  # X_*_site is n_sites x ?, X_*_visit is (n_sites * max_visits) x ? in
  # site-major order. Broadcast the site-level eta across visits.
  eta_p_site <- as.numeric(model$X_det_site %*% bp_site)
  p_mat <- matrix(eta_p_site, n_sites, max_visits)
  if (length(bp_visit)) {
    eta_p_visit <- as.numeric(model$X_det_visit %*% bp_visit)
    p_mat <- p_mat + matrix(eta_p_visit, n_sites, max_visits, byrow = TRUE)
  }
  p_mat <- stats::plogis(cl(p_mat))

  eta_pos_site <- as.numeric(model$X_pos_site %*% bpos_site)
  ep_mat <- matrix(eta_pos_site, n_sites, max_visits)
  if (length(bpos_visit)) {
    eta_pos_visit <- as.numeric(model$X_pos_visit %*% bpos_visit)
    ep_mat <- ep_mat + matrix(eta_pos_visit, n_sites, max_visits, byrow = TRUE)
  }

  list(psi = psi, p_mat = p_mat, ep_mat = ep_mat)
}

# Per-site marginal log-likelihood (latent occupancy state z integrated out in
# closed form over its two states), returned as a length-`n_sites` vector. The
# inputs are the per-cell occupancy probability `psi`, the per-visit detection
# probability matrix `p_mat` [n_sites x max_visits], the per-visit cover
# linear-predictor matrix `ep_mat` [n_sites x max_visits], and the scalar
# `log_dispersion`. This is the single source of truth shared by the fit's
# negative-log-posterior and the WAIC / PSIS-LOO pointwise log-likelihood.
.occu_cover_site_ll <- function(model, psi, p_mat, ep_mat, log_disp) {
  cl <- function(e) pmin(pmax(e, -30), 30)
  valid <- model$valid
  y     <- model$y
  y_pos <- model$y_pos
  n_sites    <- model$n_sites
  max_visits <- model$max_visits

  log_p   <- ifelse(valid, log(p_mat),     0)
  log_1mp <- ifelse(valid, log(1 - p_mat), 0)

  # Positive-arm log-density at detected visits.
  pos_mask <- valid & (y == 1L)
  log_f_pos <- matrix(0, n_sites, max_visits)
  if (identical(model$positive, "beta")) {
    phi_d <- exp(log_disp)
    mu_pos <- stats::plogis(cl(ep_mat))
    a <- mu_pos * phi_d
    b <- (1 - mu_pos) * phi_d
    # Beta density (positive support 0 < y < 1).
    dens <- lgamma(phi_d) - lgamma(a) - lgamma(b) +
            (a - 1) * log(y_pos) + (b - 1) * log(1 - y_pos)
    log_f_pos[pos_mask] <- dens[pos_mask]
  } else {
    sigma_pos <- exp(log_disp)
    dens <- -log(y_pos) - log(sigma_pos) - 0.5 * log(2 * pi) -
            0.5 * ((log(y_pos) - ep_mat) / sigma_pos)^2
    log_f_pos[pos_mask] <- dens[pos_mask]
  }

  # Per-visit log-likelihood under z = 1 (occupied).
  log_h <- ifelse(valid,
                  ifelse(y == 1L, log_p + log_f_pos, log_1mp),
                  0)

  any_det <- rowSums(y * valid, na.rm = FALSE) > 0
  log_psi   <- log(pmax(psi, 1e-300))
  log_1mpsi <- log(pmax(1 - psi, 1e-300))

  det_ll <- log_psi + rowSums(log_h)
  # No detection: psi * prod(1-p) + (1-psi). Logsumexp form for stability.
  ln_a <- log_psi   + rowSums(log_1mp)
  ln_b <- log_1mpsi
  m    <- pmax(ln_a, ln_b)
  nodet_ll <- m + log(exp(ln_a - m) + exp(ln_b - m))

  ifelse(any_det, det_ll, nodet_ll)
}


# ---------------------------------------------------------------------------
# Fitter
# ---------------------------------------------------------------------------

.tobs_fit_occu_cover <- function(model,
                                  method   = c("laplace"),
                                  priors   = NULL,
                                  max.iter = 200L,
                                  tol      = 1e-6,
                                  verbose  = TRUE,
                                  sigma.beta = 5,
                                  ...) {
  method <- match.arg(method)

  pi_list <- model$process_info
  p_occ   <- pi_list[[1L]]$p
  p_p     <- pi_list[[2L]]$p
  p_pos   <- pi_list[[3L]]$p
  n_par   <- p_occ + p_p + p_pos + 1L  # +1 for log_dispersion

  par_names <- c(
    paste0("psi_", pi_list[[1L]]$coef_names),
    paste0("p_",   pi_list[[2L]]$coef_names),
    paste0("pos_", pi_list[[3L]]$coef_names),
    if (identical(model$positive, "beta")) "log_phi" else "log_sigma_pos"
  )

  start <- numeric(n_par)
  names(start) <- par_names

  # Warm starts from a separate-fit baseline. Occurrence intercept from
  # empirical detection-any rate; detection intercept at logit(0.5); cover
  # intercept at the marginal positive-mean on its arm's link scale;
  # dispersion at a modestly broad value.
  any_det <- rowSums(model$y * model$valid) > 0
  det_rate <- max(mean(any_det), 1e-3)
  start[1L] <- stats::qlogis(min(max(det_rate, 1e-3), 1 - 1e-3))

  pos_vals <- model$y_pos[model$valid & model$y == 1L]
  if (length(pos_vals) > 0L) {
    if (identical(model$positive, "beta")) {
      start[p_occ + p_p + 1L] <- stats::qlogis(min(max(mean(pos_vals), 1e-3), 1 - 1e-3))
      start[n_par]            <- log(10)   # phi ~ 10 = moderate beta concentration
    } else {
      start[p_occ + p_p + 1L] <- mean(log(pos_vals))
      start[n_par]            <- log(stats::sd(log(pos_vals)) + 0.1)
    }
  } else {
    start[n_par] <- if (identical(model$positive, "beta")) log(10) else log(0.4)
  }

  # Gaussian prior aligned with par. Default sigma.beta on the betas
  # (a weakly-informative N(0, sigma.beta^2)); dispersion stays flat.
  pmean <- numeric(n_par)
  pprec <- numeric(n_par)
  if (isTRUE(is.null(priors)) || !isFALSE(priors)) {
    beta_idx <- seq_len(n_par - 1L)
    pprec[beta_idx] <- 1 / (sigma.beta^2)
  }

  opt <- stats::optim(start, .tobs_occu_cover_nlp,
                       model = model, pmean = pmean, pprec = pprec,
                       method = "BFGS", hessian = TRUE,
                       control = list(maxit = max.iter, reltol = tol,
                                      trace = if (isTRUE(verbose)) 1L else 0L))

  V <- tryCatch(solve(opt$hessian), error = function(e) NULL)
  if (is.null(V)) {
    warning("occu_cover: observed-Fisher Hessian not invertible; SEs unreliable.",
            call. = FALSE)
    V <- matrix(NA_real_, n_par, n_par)
  }
  se <- sqrt(pmax(diag(V), 0))

  means <- opt$par
  names(means) <- par_names
  names(se)    <- par_names
  dimnames(V)  <- list(par_names, par_names)

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = se,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = n_par,
    log_prob     = rep(-opt$value, n_draws),
    log_lik      = -opt$value,
    N            = sum(model$valid)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = NULL,
    method       = "laplace",
    positive     = model$positive,
    convergence  = list(converged = opt$convergence == 0L,
                        n_iter    = opt$counts[1L])
  )), class = c("tobs_fit", "tulpa_fit"))
}

# Draw from MVN via Cholesky; fall back to independent normals if not PD.
.occu_cover_rmvn <- function(n, mu, sigma) {
  p <- length(mu)
  if (any(!is.finite(sigma))) {
    return(matrix(rep(mu, each = n), n, p, byrow = FALSE))
  }
  L <- tryCatch(chol(sigma), error = function(e) NULL)
  z <- matrix(stats::rnorm(n * p), n, p)
  if (is.null(L)) {
    sds <- sqrt(pmax(diag(sigma), 1e-8))
    return(sweep(z * rep(sds, each = n), 2L, mu, "+"))
  }
  sweep(z %*% L, 2L, mu, "+")
}


# ---------------------------------------------------------------------------
# Dispatcher (wired into tobs.R's switch)
# ---------------------------------------------------------------------------

.dispatch_occu_cover <- function(formula, data, family, detection, y, visits,
                                  engine, priors, control,
                                  approx = "gaussian_laplace",
                                  correction = "none", ...) {
  dots <- list(...)

  if (is.null(detection)) {
    stop("occu_cover() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occu_cover() requires `y` (N x J detection-history matrix).",
         call. = FALSE)
  }
  if (is.null(dots$y_pos)) {
    stop("occu_cover() requires `y_pos` (N x J positive-cover matrix; ",
         "values used only where y == 1).", call. = FALSE)
  }

  pos_formula <- dots$positive
  if (is.null(pos_formula)) pos_formula <- detection

  # Detect the coupled spatial field(s) on the psi formula. The spatial path is
  # the joint nested-Laplace engine (shared field(s) across the psi and cover
  # arms); the non-spatial path is plain Laplace on the exact two-state
  # marginal. A weighted areal term adds a second coupled (SVC) field.
  spatial_info <- .occu_cover_spatial_fields(formula, data)
  has_spatial  <- !is.null(spatial_info)

  if (has_spatial && engine == "laplace") {
    stop("occu_cover() found a spatial term (icar/bym2) in the psi formula ",
         "but method = \"laplace\" is non-spatial. Use method = ",
         "\"nested_laplace\" for the spatial v2 path.", call. = FALSE)
  }
  if (!has_spatial && engine == "nested_laplace") {
    stop("occu_cover() with method = \"nested_laplace\" requires a spatial ",
         "term (icar() or bym2()) on the psi formula.", call. = FALSE)
  }

  fe_formula <- if (has_spatial) spatial_info$fe else formula

  # Detection / cover arms never carry a spatial term (the shared field is
  # on the latent state z, not on the observation process). Other structured
  # terms (re, temporal, ...) are not supported on any arm in v1/v2.
  .occu_cover_reject_structured(detection,   "detection")
  .occu_cover_reject_structured(pos_formula, "positive cover")

  vd_det <- .normalize_visits(visits, detection,
                              n_sites = nrow(y), max_visits = ncol(y))
  vd_pos <- .normalize_visits(visits, pos_formula,
                              n_sites = nrow(y), max_visits = ncol(y))

  model <- .tobs_build_occu_cover(
    occ_formula      = fe_formula,
    det_formula      = vd_det$det_formula,
    pos_formula      = vd_pos$det_formula,
    data             = data,
    y                = y,
    y_pos            = dots$y_pos,
    positive         = family$params$positive,
    det_visit_formula = vd_det$det_visit_formula,
    det_visit_data    = vd_det$visits,
    pos_visit_formula = vd_pos$det_visit_formula,
    pos_visit_data    = vd_pos$visits
  )

  if (has_spatial) {
    fields      <- spatial_info$fields
    base_graph  <- fields[[1L]]$graph

    # Resolve the site -> field-node map. With group_var the occupancy units
    # (sites, one per row of `data` / `y`) map onto fewer field nodes (cells),
    # so the same cell field is shared across that cell's sites (e.g. cell-year
    # sites sharing one cell). Without group_var the two coincide 1:1.
    n_cells_field <- nrow(base_graph)
    gv <- spatial_info$group_var
    if (!is.null(gv)) {
      if (!gv %in% names(data)) {
        stop(sprintf("occu_cover() group_var '%s' is not a column of data.", gv),
             call. = FALSE)
      }
      site_cell <- as.integer(data[[gv]])
      if (length(site_cell) != model$n_sites || anyNA(site_cell) ||
          min(site_cell) < 1L || max(site_cell) > n_cells_field) {
        stop(sprintf(paste0(
          "occu_cover() group_var '%s' must be an integer cell index in 1..%d, ",
          "one per site (%d sites)."), gv, n_cells_field, model$n_sites),
          call. = FALSE)
      }
    } else {
      if (model$n_sites != n_cells_field) {
        stop(sprintf(paste0(
          "occu_cover() spatial: %d sites but the graph has %d nodes. Map sites ",
          "to cells with group_var = \"<col>\" on the icar()/bym2() term (e.g. ",
          "site = cell-year), or match the site count to the graph."),
          model$n_sites, n_cells_field), call. = FALSE)
      }
      site_cell <- seq_len(model$n_sites)
    }
    model$site_cell <- site_cell
    model$n_cells   <- n_cells_field

    # joint_coupled (3-arm nested-Laplace via tulpa's cell_coupling spec) is the
    # default: outer-grid integration over (sigma, alpha [, sigma_trend,
    # alpha_trend]) with inner Newton driven by the occu_cover_{lognormal,beta}
    # cell-coupling spec. 150-300x faster than v3 at N=100 and reliably completes
    # at N=200+ where v3 trips on a missing-value compare in its outer BFGS. v3
    # pure-R nested-Laplace and v2's joint Laplace stay reachable via
    # control$engine = "v3_nested" / "v2_joint" as debug escape hatches; both
    # take only the single intercept field.
    engine_pick <- control[["engine"]] %||% "joint_coupled"
    control[["engine"]] <- NULL
    if (engine_pick %in% c("v2_joint", "v3_nested")) {
      if (length(fields) > 1L) {
        stop(sprintf(paste0(
          "occu_cover() engine \"%s\" couples a single shared field; ",
          "weighted SVC field(s) need the default joint_coupled engine."),
          engine_pick), call. = FALSE)
      }
      if (!is.null(gv)) {
        stop(sprintf(paste0(
          "occu_cover() engine \"%s\" binds the field 1:1 to sites and does ",
          "not support group_var; use the default joint_coupled engine."),
          engine_pick), call. = FALSE)
      }
      fit_args <- c(list(model = model, adj = base_graph, priors = priors),
                    control)
      fitter <- if (engine_pick == "v2_joint") .tobs_fit_occu_cover_spatial
                else .tobs_fit_occu_cover_nested
      return(do.call(fitter, fit_args))
    }
    fit_args <- c(list(model = model, fields = fields, priors = priors),
                  control)
    return(do.call(.tobs_fit_occu_cover_joint_coupled, fit_args))
  }

  fit_args <- list(model = model, method = engine, priors = priors)
  fit_args <- c(fit_args, control)
  do.call(.tobs_fit_occu_cover, fit_args)
}


# ---------------------------------------------------------------------------
# Pointwise log-likelihood (WAIC / PSIS-LOO) -- gcol33/tulpaObs#26
# ---------------------------------------------------------------------------

# Pointwise log-likelihood [n_draws x n_sites] for an occu_cover() fit: the
# per-site marginal log-likelihood (latent occupancy state integrated out)
# evaluated at each posterior draw. The spatial nested-Laplace fit samples the
# grid-integrated joint (betas + shared field) via the joint substrate; the
# non-spatial Laplace fit reuses the stored coefficient draws (no field). Both
# feed the same per-draw site-likelihood accumulation.
# Per-draw arm coefficients + dispersion + shared-field contributions for an
# occu_cover() fit. Joint path: sample the grid-integrated posterior and
# accumulate each block's field on the occupancy / cover arms. Non-spatial path:
# read the stored coefficient draws (no field). Returns a list consumed by both
# the pointwise log-likelihood and the posterior-mean plug-in.
.tobs_occu_cover_components <- function(object, n.draws = 1000L) {
  model   <- object$model
  pi_list <- model$process_info
  p_occ   <- pi_list[[1L]]$p
  p_det   <- pi_list[[2L]]$p
  p_pos   <- pi_list[[3L]]$p
  n_sites <- model$n_sites

  if (!is.null(.tobs_joint_fit(object))) {
    bundle <- .tobs_joint_draws(object, n = n.draws)
    S      <- bundle$n
    b_occ  <- bundle$b$occ
    b_det  <- bundle$b$det
    b_pos  <- bundle$b$pos
    disp   <- bundle$disp
    # Each site reads its field node (cell) via site_cell; the field draws
    # blk$z carry one column per node, so blk$z[, units] broadcasts the shared
    # cell field across that cell's sites. Identity when no group_var.
    units  <- model$site_cell %||% seq_len(n_sites)
    wfun   <- function(nm) {
      if (!nm %in% names(model$data)) {
        stop("occu_cover WAIC: trend-field weight column '", nm,
             "' is not in the fitted data.", call. = FALSE)
      }
      as.numeric(model$data[[nm]])
    }
    field_occ <- matrix(0, n_sites, S)
    field_pos <- matrix(0, n_sites, S)
    for (blk in bundle$blocks) {
      z_unit <- blk$z[, units, drop = FALSE]   # [S x n_sites]
      w <- if (is.null(blk$weight)) rep(1, n_sites) else wfun(blk$weight)
      field_occ <- field_occ + t(z_unit * blk$amp_occ) * w
      field_pos <- field_pos + t(z_unit * blk$amp_pos) * w
    }
  } else {
    draws <- object$draws
    if (is.null(draws) || !is.matrix(draws)) {
      stop("occu_cover WAIC: the fit carries no posterior draw matrix.",
           call. = FALSE)
    }
    if (!is.null(n.draws) && n.draws < nrow(draws)) {
      draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
    }
    b_occ <- draws[, seq_len(p_occ), drop = FALSE]
    b_det <- draws[, p_occ + seq_len(p_det), drop = FALSE]
    b_pos <- draws[, p_occ + p_det + seq_len(p_pos), drop = FALSE]
    disp  <- exp(draws[, ncol(draws)])         # trailing log-dispersion column
    S     <- nrow(draws)
    field_occ <- matrix(0, n_sites, S)
    field_pos <- matrix(0, n_sites, S)
  }
  list(b_occ = b_occ, b_det = b_det, b_pos = b_pos, disp = disp,
       field_occ = field_occ, field_pos = field_pos)
}

.tobs_ploglik_occu_cover <- function(object, n.draws = 1000L) {
  c0 <- .tobs_occu_cover_components(object, n.draws)
  .occu_cover_ploglik_core(object$model, c0$b_occ, c0$b_det, c0$b_pos,
                           c0$disp, c0$field_occ, c0$field_pos)
}

# Pointwise log-likelihood at the posterior mean (length n_sites): the per-site
# marginal evaluated at the mean arm coefficients, mean dispersion, and mean
# shared field -- the same core driven by a one-draw mean.
.tobs_occu_cover_loglik_at_mean <- function(object, n.draws = 1000L) {
  c0 <- .tobs_occu_cover_components(object, n.draws)
  as.numeric(.occu_cover_ploglik_core(
    object$model,
    matrix(colMeans(c0$b_occ), nrow = 1L),
    matrix(colMeans(c0$b_det), nrow = 1L),
    matrix(colMeans(c0$b_pos), nrow = 1L),
    mean(c0$disp),
    matrix(rowMeans(c0$field_occ), ncol = 1L),
    matrix(rowMeans(c0$field_pos), ncol = 1L)
  ))
}

# Accumulate the [S x n_sites] pointwise log-likelihood from per-draw arm
# coefficients (`b_occ` / `b_det` / `b_pos`, each [S x p_arm]), per-draw
# dispersion `disp` [S], and per-arm shared-field contributions `field_occ` /
# `field_pos` [n_sites x S] (zero matrices for a non-spatial fit). The detection
# and cover arms split into site-level and visit-level coefficient blocks exactly
# as the fitter packs them.
# Per-draw linear-predictor blocks for an occu_cover() fit: the occupancy
# `eta_psi_all` [n_sites x S] (with the shared field), and the site-level +
# optional visit-level detection / cover blocks the fitter packs. Split out so
# the pointwise log-likelihood, the posterior predictive check, and the PIT all
# build the per-draw psi / p / cover predictors from one source.
.occu_cover_eta_components <- function(model, b_occ, b_det, b_pos,
                                       field_occ, field_pos) {
  p_det_site <- ncol(model$X_det_site)
  p_pos_site <- ncol(model$X_pos_site)
  has_det_visit <- !is.null(model$X_det_visit)
  has_pos_visit <- !is.null(model$X_pos_visit)
  list(
    eta_psi_all = tcrossprod(model$X_occ, b_occ) + field_occ,
    eta_p_site_all = tcrossprod(model$X_det_site,
                                b_det[, seq_len(p_det_site), drop = FALSE]),
    eta_pos_site_all = tcrossprod(model$X_pos_site,
                                  b_pos[, seq_len(p_pos_site), drop = FALSE]) +
                       field_pos,
    eta_p_visit_all = if (has_det_visit) {
      tcrossprod(model$X_det_visit,
                 b_det[, p_det_site + seq_len(ncol(model$X_det_visit)),
                       drop = FALSE])
    } else NULL,
    eta_pos_visit_all = if (has_pos_visit) {
      tcrossprod(model$X_pos_visit,
                 b_pos[, p_pos_site + seq_len(ncol(model$X_pos_visit)),
                       drop = FALSE])
    } else NULL,
    has_det_visit = has_det_visit, has_pos_visit = has_pos_visit
  )
}

# Draw `d`'s occupancy predictor and the [n_sites x max_visits] detection /
# cover linear predictors, folding the visit-level block in site-major order.
.occu_cover_draw_eta <- function(comp, d, n_sites, max_visits) {
  p_eta <- matrix(comp$eta_p_site_all[, d], n_sites, max_visits)
  if (comp$has_det_visit) {
    p_eta <- p_eta + matrix(comp$eta_p_visit_all[, d], n_sites, max_visits,
                            byrow = TRUE)
  }
  ep_mat <- matrix(comp$eta_pos_site_all[, d], n_sites, max_visits)
  if (comp$has_pos_visit) {
    ep_mat <- ep_mat + matrix(comp$eta_pos_visit_all[, d], n_sites, max_visits,
                              byrow = TRUE)
  }
  list(psi_eta = comp$eta_psi_all[, d], p_eta = p_eta, ep_mat = ep_mat)
}

.occu_cover_ploglik_core <- function(model, b_occ, b_det, b_pos, disp,
                                     field_occ, field_pos) {
  n_sites    <- model$n_sites
  max_visits <- model$max_visits
  S <- nrow(b_occ)
  comp <- .occu_cover_eta_components(model, b_occ, b_det, b_pos,
                                     field_occ, field_pos)
  cl <- function(e) pmin(pmax(e, -30), 30)
  ll <- matrix(0, S, n_sites)
  for (d in seq_len(S)) {
    de    <- .occu_cover_draw_eta(comp, d, n_sites, max_visits)
    psi   <- stats::plogis(cl(de$psi_eta))
    p_mat <- stats::plogis(cl(de$p_eta))
    ll[d, ] <- .occu_cover_site_ll(model, psi, p_mat, de$ep_mat, log(disp[d]))
  }
  ll
}


# ---------------------------------------------------------------------------
# Posterior predictive check + PIT (gcol33/tulpaObs#27)
# ---------------------------------------------------------------------------

# Posterior predictive check for an occu_cover() fit. Per draw, the latent
# occupancy z is sampled from its full conditional given the detection history
# (the spOccupancy construction), detection replicates from Bernoulli(z p), and
# cover replicates at the detected visits from the fitted positive part. The
# discrepancy is a Freeman-Tukey (or chi-squared) sum over the detection cells
# plus the positive-part cells, returning a Bayesian p-value.
.tobs_ppc_occu_cover <- function(object,
                                 fit.stat = c("freeman-tukey", "chi-squared"),
                                 n.samples = 500) {
  fit.stat <- match.arg(fit.stat)
  model    <- object$model
  positive <- model$positive %||% "lognormal"
  c0   <- .tobs_occu_cover_components(object, n.samples)
  comp <- .occu_cover_eta_components(model, c0$b_occ, c0$b_det, c0$b_pos,
                                     c0$field_occ, c0$field_pos)
  disp <- c0$disp
  S <- nrow(c0$b_occ)
  n_sites <- model$n_sites; max_visits <- model$max_visits
  y <- model$y; y_pos <- model$y_pos; valid <- model$valid
  cl <- function(e) pmin(pmax(e, -30), 30)
  stat_fn <- if (fit.stat == "freeman-tukey") {
    function(o, e) sum((sqrt(o) - sqrt(e))^2, na.rm = TRUE)
  } else {
    function(o, e) sum((o - e)^2 / (e + 1e-10), na.rm = TRUE)
  }
  any_det <- rowSums(y * valid, na.rm = TRUE) > 0
  n_valid <- rowSums(valid)
  pos_mask <- valid & (y == 1L)

  fit_y <- fit_rep <- numeric(S)
  for (s in seq_len(S)) {
    de    <- .occu_cover_draw_eta(comp, s, n_sites, max_visits)
    psi   <- stats::plogis(cl(de$psi_eta))
    p_mat <- stats::plogis(cl(de$p_eta))
    prod1mp <- exp(rowSums(ifelse(valid, log(1 - p_mat), 0)))
    z_prob <- ifelse(any_det, 1, psi * prod1mp / (psi * prod1mp + (1 - psi)))
    z_prob[n_valid == 0L] <- psi[n_valid == 0L]
    z <- stats::rbinom(n_sites, 1, z_prob)

    exp_det <- z * p_mat
    yrep <- matrix(stats::rbinom(n_sites * max_visits, 1, as.vector(exp_det)),
                   n_sites, max_visits)
    det_obs <- stat_fn(y[valid], exp_det[valid])
    det_rep <- stat_fn(yrep[valid], exp_det[valid])

    if (positive == "lognormal") {
      sg   <- disp[s]
      Epos <- exp(cl(de$ep_mat) + sg^2 / 2)
      ypos_rep <- matrix(exp(stats::rnorm(n_sites * max_visits,
                                          as.vector(de$ep_mat), sg)),
                         n_sites, max_visits)
    } else {
      phi  <- disp[s]
      mu   <- stats::plogis(cl(de$ep_mat))
      Epos <- mu
      ypos_rep <- matrix(stats::rbeta(n_sites * max_visits,
                                      as.vector(mu) * phi,
                                      (1 - as.vector(mu)) * phi),
                         n_sites, max_visits)
    }
    pos_obs <- if (any(pos_mask)) stat_fn(y_pos[pos_mask], Epos[pos_mask]) else 0
    pos_rep <- if (any(pos_mask)) stat_fn(ypos_rep[pos_mask], Epos[pos_mask]) else 0

    fit_y[s]   <- det_obs + pos_obs
    fit_rep[s] <- det_rep + pos_rep
  }
  list(fit.y = fit_y, fit.y.rep = fit_rep,
       bayesian.p = mean(fit_rep > fit_y))
}

# Randomized PIT for an occu_cover() fit, on the per-site detection summary
# (any-detection vs all-zero) marginalized over the latent occupancy state, with
# the shared field projected per site. The detected / non-detected outcome is the
# ordered event; the left and right CDF limits feed the engine's randomized PIT.
.tobs_pit_occu_cover <- function(object, n.samples = 250) {
  model <- object$model
  c0    <- .tobs_occu_cover_components(object, n.samples)
  comp  <- .occu_cover_eta_components(model, c0$b_occ, c0$b_det, c0$b_pos,
                                      c0$field_occ, c0$field_pos)
  S <- nrow(c0$b_occ)
  n_sites <- model$n_sites; max_visits <- model$max_visits
  y <- model$y; valid <- model$valid
  cl <- function(e) pmin(pmax(e, -30), 30)
  any_det <- rowSums(y * valid, na.rm = TRUE) > 0

  Fl <- matrix(0, S, n_sites); Fu <- matrix(0, S, n_sites)
  for (s in seq_len(S)) {
    de    <- .occu_cover_draw_eta(comp, s, n_sites, max_visits)
    psi   <- stats::plogis(cl(de$psi_eta))
    p_mat <- stats::plogis(cl(de$p_eta))
    prod1mp <- exp(rowSums(ifelse(valid, log(1 - p_mat), 0)))
    pdet0 <- psi * prod1mp + (1 - psi)        # P(all-zero detection outcome)
    nd <- !any_det
    Fu[s, nd]      <- pdet0[nd]                # all-zero observed: below the mass
    Fl[s, any_det] <- pdet0[any_det]           # detected observed: above the mass
    Fu[s, any_det] <- 1
  }
  tulpa::tulpa_pit(Fu, cdf_lower = Fl)
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate joint occupancy-detection + cover data
#'
#' Per-cell mixture: latent z_i ~ Bernoulli(psi_i), per-visit detection
#' y_ij | z_i = 1 ~ Bernoulli(p_ij), per-visit cover y_pos_ij | y_ij = 1
#' drawn from `positive` (`"beta"` or `"lognormal"`) on the cover-arm linear
#' predictor. Used by the recovery test and as the generator for the
#' joint occupancy + cover hurdle family (see [occu_cover()]).
#'
#' When `adj` is supplied (a square adjacency matrix), an ICAR field
#' `f[1..N]` is drawn from `MVN(0, Q^-)` (with sum-to-zero constraint),
#' and the linear predictors become
#'
#'     eta_psi_i = X_psi[i, ] %*% beta_occ + sigma * f[i]
#'     eta_pos_ij = X_pos[i, ] %*% beta_pos + alpha * sigma * f[i]
#'
#' matching the v2 nested-Laplace fit's parameterisation.
#'
#' @param N Number of sites (cells).
#' @param J Number of visits per site.
#' @param n_occ_covs,n_det_covs,n_pos_covs Number of covariates on each arm
#'   (drawn IID standard normal).
#' @param beta_occ,beta_p,beta_pos Coefficient vectors c(intercept, slopes).
#'   Defaults pick weakly-informative values: psi intercept at logit(0.4),
#'   p intercept at logit(0.5), cover intercept on the appropriate link.
#' @param positive `"beta"` or `"lognormal"`.
#' @param phi Beta precision when `positive = "beta"` (default 30).
#' @param sigma_pos Lognormal residual SD when `positive = "lognormal"`
#'   (default 0.4).
#' @param adj Optional N x N adjacency matrix. When supplied, generates the
#'   shared ICAR field; when NULL, the simulator is non-spatial (matches v1).
#' @param sigma Spatial field amplitude (used only when `adj` is supplied).
#' @param alpha Cover-arm scaling on the shared field (used only when `adj`
#'   is supplied). 1.0 = arms see the field identically; positive = same sign,
#'   negative = opposite.
#' @param trend Logical; when `TRUE` (and `adj` is supplied) a SECOND shared
#'   ICAR field `f2` (a spatially-varying temporal trend) is generated on the
#'   same graph, weighted by a per-cell covariate `time` drawn IID standard
#'   normal. The trend enters the occupancy and cover predictors as
#'   `sigma_trend * time_i * f2[i]` (occupancy) and
#'   `alpha_trend * sigma_trend * time_i * f2[i]` (cover); the detection
#'   predictor is unaffected. The `time` covariate is per-cell and broadcast
#'   to every visit of that cell.
#' @param sigma_trend Trend-field amplitude (used only when `trend = TRUE`).
#' @param alpha_trend Cover-arm scaling on the trend field (used only when
#'   `trend = TRUE`).
#' @param seed Optional integer seed.
#' @return A list with `y` (N x J detection matrix), `y_pos` (N x J cover
#'   matrix, NA where not detected), `data` (per-site covariate frame, gaining
#'   a `time` column when `trend = TRUE`), `visit_data` (per-visit covariate
#'   frame, N*J rows in site-major order), and `truth` (the coefficients,
#'   dispersion, and field(s) if generated; `f2`, `sigma_trend`, `alpha_trend`,
#'   and `time` when `trend = TRUE`).
#' @export
simulate_occu_cover <- function(N             = 200L,
                                 J             = 4L,
                                 n_occ_covs    = 1L,
                                 n_det_covs    = 1L,
                                 n_pos_covs    = 1L,
                                 beta_occ      = NULL,
                                 beta_p        = NULL,
                                 beta_pos      = NULL,
                                 positive      = c("lognormal", "beta"),
                                 phi           = 30,
                                 sigma_pos     = 0.4,
                                 adj           = NULL,
                                 sigma         = 0.6,
                                 alpha         = 1.0,
                                 trend         = FALSE,
                                 sigma_trend   = 0.6,
                                 alpha_trend   = 1.0,
                                 seed          = NULL) {
  positive <- match.arg(positive)
  if (!is.null(seed)) set.seed(seed)
  N <- as.integer(N); J <- as.integer(J)

  if (is.null(beta_occ)) beta_occ <- c(stats::qlogis(0.4), stats::runif(n_occ_covs, -0.5, 0.5))
  if (is.null(beta_p))   beta_p   <- c(0, stats::runif(n_det_covs, -0.5, 0.5))
  if (is.null(beta_pos)) {
    pos_int <- if (positive == "beta") stats::qlogis(0.3) else log(0.1)
    beta_pos <- c(pos_int, stats::runif(n_pos_covs, -0.5, 0.5))
  }

  # Optional shared ICAR field(s). Draw each f as N(0, Q^-) via the
  # eigendecomposition of Q; the constant (null) component is dropped, giving a
  # zero-mean draw on the constrained space, then divide by sqrt(scale_q) so the
  # field has geo-mean marginal variance 1 (the Sorbye-Rue convention; `sigma *
  # f` then has geo-mean marginal SD sigma, matching INLA's `scale.model = TRUE`
  # and the fitter's parameterisation).
  f  <- numeric(N)
  f2 <- numeric(N)
  time_cov <- numeric(N)
  if (!is.null(adj)) {
    if (!is.matrix(adj) || nrow(adj) != N || ncol(adj) != N) {
      stop("adj must be an N x N adjacency matrix.", call. = FALSE)
    }
    Q       <- .occu_cover_icar_Q(adj)
    scale_q <- .occu_cover_icar_scale(adj)
    eig <- eigen(Q, symmetric = TRUE)
    keep <- eig$values > 1e-8
    draw_field <- function() {
      z_white <- stats::rnorm(sum(keep))
      fk <- as.numeric(eig$vectors[, keep, drop = FALSE] %*%
                         (z_white / sqrt(eig$values[keep])))
      fk <- fk - mean(fk)
      fk / sqrt(scale_q)
    }
    f <- draw_field()
    if (isTRUE(trend)) {
      f2       <- draw_field()
      time_cov <- as.numeric(scale(stats::rnorm(N)))
    }
  }

  # Site-level covariates (psi predictor).
  occ_covs <- data.frame(matrix(stats::rnorm(N * n_occ_covs), N, n_occ_covs))
  names(occ_covs) <- paste0("occ_cov", seq_len(n_occ_covs))
  X_occ <- stats::model.matrix(~ ., occ_covs)
  eta_psi <- as.vector(X_occ %*% beta_occ) + sigma * f
  if (!is.null(adj) && isTRUE(trend)) {
    eta_psi <- eta_psi + sigma_trend * time_cov * f2
  }
  psi <- stats::plogis(eta_psi)
  z_state <- stats::rbinom(N, 1L, psi)

  # Visit-level covariates (p and cover predictors). Same draw used for both
  # arms, mirroring how `tobs_data()`'s `det.covs` matrices feed both formulas.
  det_covs <- data.frame(matrix(stats::rnorm(N * J * n_det_covs), N * J, n_det_covs))
  names(det_covs) <- paste0("det_cov", seq_len(n_det_covs))
  pos_covs <- data.frame(matrix(stats::rnorm(N * J * n_pos_covs), N * J, n_pos_covs))
  names(pos_covs) <- paste0("pos_cov", seq_len(n_pos_covs))
  visit_data <- cbind(det_covs, pos_covs)

  X_p   <- stats::model.matrix(~ ., det_covs)
  X_pos <- stats::model.matrix(~ ., pos_covs)
  eta_p   <- as.vector(X_p   %*% beta_p)
  eta_pos_base <- as.vector(X_pos %*% beta_pos)

  y     <- matrix(0L, N, J)
  y_pos <- matrix(NA_real_, N, J)

  for (i in seq_len(N)) {
    for (j in seq_len(J)) {
      idx <- (i - 1L) * J + j
      if (z_state[i] == 1L) {
        p_ij <- stats::plogis(eta_p[idx])
        d <- stats::rbinom(1L, 1L, p_ij)
        y[i, j] <- d
        if (d == 1L) {
          eta_pos_ij <- eta_pos_base[idx] + alpha * sigma * f[i]
          if (!is.null(adj) && isTRUE(trend)) {
            eta_pos_ij <- eta_pos_ij + alpha_trend * sigma_trend * time_cov[i] * f2[i]
          }
          if (positive == "beta") {
            mu <- stats::plogis(eta_pos_ij)
            y_pos[i, j] <- stats::rbeta(1L, mu * phi, (1 - mu) * phi)
          } else {
            y_pos[i, j] <- exp(stats::rnorm(1L, eta_pos_ij, sigma_pos))
          }
        }
      }
    }
  }

  occ_out <- occ_covs
  has_trend <- !is.null(adj) && isTRUE(trend)
  if (has_trend) occ_out$time <- time_cov

  list(
    y          = y,
    y_pos      = y_pos,
    data       = occ_out,
    visit_data = visit_data,
    adj        = adj,
    truth      = list(
      beta_occ    = beta_occ,
      beta_p      = beta_p,
      beta_pos    = beta_pos,
      psi         = psi,
      z           = z_state,
      positive    = positive,
      phi         = if (positive == "beta")      phi       else NA_real_,
      sigma_pos   = if (positive == "lognormal") sigma_pos else NA_real_,
      f           = f,
      sigma       = if (!is.null(adj)) sigma else NA_real_,
      alpha       = if (!is.null(adj)) alpha else NA_real_,
      f2          = if (has_trend) f2          else NULL,
      time        = if (has_trend) time_cov    else NULL,
      sigma_trend = if (has_trend) sigma_trend else NA_real_,
      alpha_trend = if (has_trend) alpha_trend else NA_real_
    )
  )
}
