# =============================================================================
# occu_multiscale_cover.R - three-level occupancy + cover hurdle
# (gcol33/tulpaObs#29).
#
# Cell-level occupancy (psi) gates plot-level availability (theta), which gates
# per-visit detection (p) and the cover hurdle (pos). The cell is the areal
# field node; plots are the y-rows; visits are the y-columns. Both z (over
# cells) and a (over plots) marginalize in closed form, so the joint marginal
# log-likelihood is exact (the cell-coupling spec in
# src/cell_coupling_occu_multiscale_cover.h). Spatial-only: the fit is the
# joint nested-Laplace cell-coupling path (.tobs_fit_occu_multiscale_cover_*).
#
# Builder, dispatcher, and simulator. The joint fitter lives in
# occu_multiscale_cover_joint_coupled.R.
# =============================================================================


# Build the multiscale model object. `y` / `y_pos` are [n_plots x max_visits];
# `plot_cell` maps each plot (row) to its areal cell (1..n_cells). The psi
# design is per-CELL (cell-level covariates aggregated to the first plot of each
# cell), theta is per-PLOT, p / pos are per-visit (plot-level block + optional
# visit-varying). Detection NAs mask visits out; cover is read only where y == 1.
.tobs_build_occu_multiscale_cover <- function(occ_formula, theta_formula,
                                              det_formula, pos_formula,
                                              data, y, y_pos, plot_cell, n_cells,
                                              positive,
                                              det_visit_formula = NULL,
                                              det_visit_data    = NULL,
                                              pos_visit_formula = NULL,
                                              pos_visit_data    = NULL) {
  if (!is.matrix(y) || !is.matrix(y_pos)) {
    stop("y and y_pos must be matrices (n_plots x max_visits).", call. = FALSE)
  }
  if (!all(dim(y) == dim(y_pos))) {
    stop("y and y_pos must have identical dimensions.", call. = FALSE)
  }
  if (nrow(y) != nrow(data)) {
    stop(sprintf("y has %d rows (plots) but data has %d rows.",
                 nrow(y), nrow(data)), call. = FALSE)
  }
  n_plots    <- nrow(y)
  max_visits <- ncol(y)
  plot_cell  <- as.integer(plot_cell)
  if (length(plot_cell) != n_plots || anyNA(plot_cell) ||
      min(plot_cell) < 1L || max(plot_cell) > n_cells) {
    stop(sprintf(paste0(
      "plot_cell must be an integer cell index in 1..%d, one per plot ",
      "(%d plots)."), n_cells, n_plots), call. = FALSE)
  }

  y_int <- matrix(as.integer(y), n_plots, max_visits)
  valid <- !is.na(y_int)
  if (any(y_int[valid] != 0L & y_int[valid] != 1L)) {
    stop("y must contain only 0, 1, or NA (binary detection per visit).",
         call. = FALSE)
  }
  y_int[!valid] <- 0L

  y_pos_num <- matrix(as.numeric(y_pos), n_plots, max_visits)
  pos_mask  <- valid & (y_int == 1L)
  if (any(!is.finite(y_pos_num[pos_mask]))) {
    stop("y_pos must be finite at every detected visit (y == 1).", call. = FALSE)
  }
  if (identical(positive, "beta")) {
    if (any(pos_mask & (y_pos_num <= 0 | y_pos_num >= 1))) {
      stop("Beta positive arm requires 0 < y_pos < 1 at every detected visit.",
           call. = FALSE)
    }
  } else if (any(pos_mask & (y_pos_num <= 0))) {
    stop("Lognormal positive arm requires y_pos > 0 at every detected visit.",
         call. = FALSE)
  }
  y_pos_num[!pos_mask] <- 0

  # The state-process formula's FE part is cell-level. Build it over plots, then
  # aggregate to one row per cell (the first plot of each cell). Cells with no
  # plot get the intercept + zeroed covariates (prior-mean prediction node).
  X_occ_plot <- stats::model.matrix(occ_formula, data)
  first_plot <- match(seq_len(n_cells), plot_cell)
  X_psi <- matrix(0, n_cells, ncol(X_occ_plot),
                  dimnames = list(NULL, colnames(X_occ_plot)))
  int_col <- match("(Intercept)", colnames(X_occ_plot))
  if (!is.na(int_col)) X_psi[, int_col] <- 1
  present <- !is.na(first_plot)
  X_psi[present, ] <- X_occ_plot[first_plot[present], , drop = FALSE]

  X_theta    <- stats::model.matrix(theta_formula, data)
  X_p_site   <- stats::model.matrix(det_formula, data)
  X_pos_site <- stats::model.matrix(pos_formula, data)
  X_p_visit  <- .occu_cover_build_visit_X(det_visit_formula, det_visit_data,
                                          n_plots, max_visits, arm = "detection")
  X_pos_visit <- .occu_cover_build_visit_X(pos_visit_formula, pos_visit_data,
                                           n_plots, max_visits, arm = "positive cover")

  det_coef_names <- colnames(X_p_site)
  pos_coef_names <- colnames(X_pos_site)
  if (!is.null(X_p_visit))   det_coef_names <- c(det_coef_names, colnames(X_p_visit))
  if (!is.null(X_pos_visit)) pos_coef_names <- c(pos_coef_names, colnames(X_pos_visit))

  structure(list(
    model_type  = "occu_multiscale_cover",
    positive    = positive,
    y           = y_int,
    y_pos       = y_pos_num,
    valid       = valid,
    n_cells     = n_cells,
    n_plots     = n_plots,
    n_sites     = n_plots,   # generic accessors (nobs) read plot-level n
    max_visits  = max_visits,
    plot_cell   = plot_cell,
    X_psi       = X_psi,
    X_theta     = X_theta,
    X_p_site    = X_p_site,
    X_pos_site  = X_pos_site,
    X_p_visit   = X_p_visit,
    X_pos_visit = X_pos_visit,
    formulas    = list(psi = occ_formula, theta = theta_formula,
                       p = det_formula, pos = pos_formula),
    data        = data,
    process_info = list(
      list(name = "psi",   p = ncol(X_psi),
           coef_names = colnames(X_psi),   link = "logit"),
      list(name = "theta", p = ncol(X_theta),
           coef_names = colnames(X_theta), link = "logit"),
      list(name = "p",     p = length(det_coef_names),
           coef_names = det_coef_names,    link = "logit"),
      list(name = "pos",   p = length(pos_coef_names),
           coef_names = pos_coef_names,
           link = if (positive == "beta") "logit" else "identity")
    )
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "occu_multiscale_cover")
# ---------------------------------------------------------------------------

# Posterior-mean arm coefficients, split by process block.
.occu_ms_cover_arm_betas <- function(object) {
  m  <- object$means
  pp <- vapply(object$model$process_info, function(x) as.integer(x$p), integer(1))
  off <- cumsum(c(0L, pp))
  list(psi   = m[off[1] + seq_len(pp[1])],
       theta = m[off[2] + seq_len(pp[2])],
       p     = m[off[3] + seq_len(pp[3])],
       pos   = m[off[4] + seq_len(pp[4])])
}

# fitted(): per-unit posterior-mean predictions for each level of the hierarchy.
# psi is per cell (the areal field is added with coefficient 1); theta / p / cover
# are per plot (site-level designs; the visit-level part of p / cover is constant
# across visits at the posterior mean and so reads off the site design). The cover
# arm carries the shared field with the copy coefficient `alpha`. Lognormal cover
# reports the conditional MEAN E[cover | detected] = exp(eta + sigma_pos^2 / 2);
# beta cover reports the conditional mean plogis(eta). `p_marginal` is the
# per-plot marginal detection probability psi_cell * theta * p (a single visit).
.tobs_fitted_occu_multiscale_cover <- function(object) {
  model <- object$model
  b     <- .occu_ms_cover_arm_betas(object)
  field <- as.numeric(object$spatial_field %||% rep(0, model$n_cells))   # per cell
  alpha <- unname(object$means["alpha"]); if (!is.finite(alpha)) alpha <- 0
  pc    <- model$plot_cell

  eta_psi   <- as.numeric(model$X_psi %*% b$psi) + field                 # per cell
  psi_cell  <- stats::plogis(eta_psi)
  theta     <- stats::plogis(as.numeric(model$X_theta %*% b$theta))      # per plot
  pdet      <- stats::plogis(as.numeric(model$X_p_site %*% b$p))         # per plot
  eta_pos   <- as.numeric(model$X_pos_site %*% b$pos) + alpha * field[pc] # per plot
  if (identical(model$positive, "beta")) {
    cover <- stats::plogis(eta_pos)
  } else {
    sigma_pos <- unname(object$means["phi_pos"]); if (!is.finite(sigma_pos)) sigma_pos <- 0
    cover <- exp(eta_pos + 0.5 * sigma_pos^2)
  }
  list(psi = psi_cell, theta = theta, p = pdet, cover = cover,
       field = field, p_marginal = psi_cell[pc] * theta * pdet)
}

# predict(): in-sample posterior arm predictions. The areal field is tied to the
# cell graph, so prediction at new covariates / cells (X.0 / newdata) is not
# supported (no field at an unseen cell), matching ms_occu_cover_spatial().
.tobs_predict_occu_multiscale_cover <- function(object, type = "state") {
  type <- switch(type,
                 state = "psi", occupancy = "psi", psi = "psi",
                 availability = "theta", theta = "theta",
                 detection = "p", p = "p",
                 cover = "cover", cover_cond = "cover",
                 stop(sprintf(paste0("predict() type '%s' is not supported for ",
                      "occu_multiscale_cover(). Use \"state\"/\"psi\", ",
                      "\"availability\"/\"theta\", \"detection\"/\"p\", or ",
                      "\"cover\"."), type), call. = FALSE))
  fv <- .tobs_fitted_occu_multiscale_cover(object)
  switch(type, psi = fv$psi, theta = fv$theta, p = fv$p, cover = fv$cover)
}


# ---------------------------------------------------------------------------
# Dispatcher (wired into tobs.R's switch). Spatial-only: the psi formula must
# carry an icar()/bym2() term with group_var naming the per-plot cell column.
# ---------------------------------------------------------------------------
.dispatch_occu_multiscale_cover <- function(formula, data, family, detection,
                                            y, visits, engine, priors, control,
                                            approx = "gaussian_laplace",
                                            correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection)) {
    stop("occu_multiscale_cover() requires a `detection` formula (the ",
         "per-visit detection arm p).", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occu_multiscale_cover() requires `y` (an n_plots x max_visits ",
         "detection-history matrix).", call. = FALSE)
  }
  if (is.null(dots$y_pos)) {
    stop("occu_multiscale_cover() requires `y_pos` (n_plots x max_visits; ",
         "values used only where y == 1).", call. = FALSE)
  }

  theta_formula <- dots$availability %||% ~ 1
  pos_formula   <- dots$positive %||% detection

  spatial_info <- .occu_cover_spatial_fields(formula, data)
  if (is.null(spatial_info)) {
    stop("occu_multiscale_cover() is spatial: the state-process formula must ",
         "carry an areal field (icar(graph = adj, group_var = \"<cell>\")) ",
         "naming the per-plot cell column. method must be \"nested_laplace\".",
         call. = FALSE)
  }
  gv <- spatial_info$group_var
  if (is.null(gv)) {
    stop("occu_multiscale_cover() requires group_var on the icar()/bym2() ",
         "term naming the integer cell column (plot -> cell), e.g. ",
         "icar(graph = adj, group_var = \"cell\").", call. = FALSE)
  }
  if (!gv %in% names(data)) {
    stop(sprintf("occu_multiscale_cover() group_var '%s' is not a column of data.",
                 gv), call. = FALSE)
  }

  # Detection / availability / cover arms carry no structured terms.
  .occu_cover_reject_structured(detection,     "detection")
  .occu_cover_reject_structured(theta_formula, "availability")
  .occu_cover_reject_structured(pos_formula,   "positive cover")

  fields  <- spatial_info$fields
  n_cells <- nrow(fields[[1L]]$graph)
  plot_cell <- as.integer(data[[gv]])

  vd_det <- .normalize_visits(visits, detection,
                              n_sites = nrow(y), max_visits = ncol(y))
  vd_pos <- .normalize_visits(visits, pos_formula,
                              n_sites = nrow(y), max_visits = ncol(y))

  model <- .tobs_build_occu_multiscale_cover(
    occ_formula       = spatial_info$fe,
    theta_formula     = theta_formula,
    det_formula       = vd_det$det_formula,
    pos_formula       = vd_pos$det_formula,
    data              = data,
    y                 = y,
    y_pos             = dots$y_pos,
    plot_cell         = plot_cell,
    n_cells           = n_cells,
    positive          = family$params$positive,
    det_visit_formula = vd_det$det_visit_formula,
    det_visit_data    = vd_det$visits,
    pos_visit_formula = vd_pos$det_visit_formula,
    pos_visit_data    = vd_pos$visits
  )

  control[["engine"]] <- NULL
  fit_args <- c(list(model = model, fields = fields, priors = priors), control)
  do.call(.tobs_fit_occu_multiscale_cover_joint_coupled, fit_args)
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate a three-level occupancy + cover hurdle data set
#'
#' Draws cell occupancy `z_c ~ Bernoulli(psi_c)`, plot availability
#' `a_cj | z_c=1 ~ Bernoulli(theta_cj)`, visit detection
#' `y_cjv | a_cj=1 ~ Bernoulli(p_cjv)`, and the cover hurdle
#' `cover_cjv | y_cjv=1 ~ f_pos`, with a shared areal ICAR field on the
#' occupancy (`sigma`) and cover (`alpha * sigma`) arms. The data-generating
#' process for [occu_multiscale_cover()] recovery tests.
#'
#' @param n_cells number of areal cells (graph nodes).
#' @param plots_per_cell plots (the availability units) per cell.
#' @param visits_per_plot detection replicate visits per plot.
#' @param beta_psi,beta_theta,beta_p,beta_pos length-2 `c(intercept, slope)`
#'   coefficients for the four arms (one covariate each).
#' @param positive `"lognormal"` or `"beta"` cover arm.
#' @param phi cover dispersion (lognormal log-scale SD, or beta precision).
#' @param adj `n_cells x n_cells` adjacency; default a 1-D chain.
#' @param sigma areal-field marginal SD on the occupancy arm.
#' @param alpha cover-arm field scaling (`alpha * sigma` on the cover arm).
#' @param seed optional RNG seed.
#' @return A list with `y`, `y_pos` (`[n_plots x visits_per_plot]`), the
#'   plot-level `data` (cell id `cell`, covariates), `adj`, and `truth`.
#' @export
simulate_occu_multiscale_cover <- function(n_cells = 60L,
                                           plots_per_cell = 4L,
                                           visits_per_plot = 2L,
                                           beta_psi   = c(0.4, 0.6),
                                           beta_theta = c(0.2, 0.5),
                                           beta_p     = c(0.0, 0.5),
                                           beta_pos   = c(log(0.10), -0.4),
                                           positive   = c("lognormal", "beta"),
                                           phi        = 0.35,
                                           adj        = NULL,
                                           sigma      = 0.7,
                                           alpha      = 1.0,
                                           seed       = NULL) {
  positive <- match.arg(positive)
  if (!is.null(seed)) set.seed(seed)
  n_cells <- as.integer(n_cells)
  J <- as.integer(visits_per_plot)

  if (is.null(adj)) {
    adj <- matrix(0, n_cells, n_cells)
    for (i in seq_len(n_cells - 1L)) { adj[i, i + 1L] <- 1; adj[i + 1L, i] <- 1 }
  }

  # ICAR field, Sorbye-Rue scaled (geo-mean marginal variance 1).
  Q       <- .occu_cover_icar_Q(adj)
  scale_q <- .occu_cover_icar_scale(adj)
  eig  <- eigen(Q, symmetric = TRUE)
  keep <- eig$values > 1e-8
  z_white <- stats::rnorm(sum(keep))
  f <- as.numeric(eig$vectors[, keep, drop = FALSE] %*% (z_white / sqrt(eig$values[keep])))
  f <- (f - mean(f)) / sqrt(scale_q)

  # Plot -> cell map (balanced).
  plot_cell <- rep(seq_len(n_cells), each = plots_per_cell)
  n_plots   <- length(plot_cell)

  # Covariates: cell-level (psi), plot-level (theta, p, pos).
  x_cell  <- stats::rnorm(n_cells)
  x_plot  <- stats::rnorm(n_plots)
  x_pdet  <- stats::rnorm(n_plots)
  x_cov   <- stats::rnorm(n_plots)

  # Cell occupancy.
  eta_psi <- beta_psi[1L] + beta_psi[2L] * x_cell + sigma * f
  psi     <- stats::plogis(eta_psi)
  z       <- stats::rbinom(n_cells, 1L, psi)

  # Plot availability (only meaningful where the cell is occupied).
  eta_theta <- beta_theta[1L] + beta_theta[2L] * x_plot
  theta     <- stats::plogis(eta_theta)
  a         <- stats::rbinom(n_plots, 1L, theta) * z[plot_cell]

  # Detection + cover per visit.
  y     <- matrix(0L, n_plots, J)
  y_pos <- matrix(NA_real_, n_plots, J)
  eta_p   <- beta_p[1L]   + beta_p[2L]   * x_pdet
  p_plot  <- stats::plogis(eta_p)
  eta_pos <- beta_pos[1L] + beta_pos[2L] * x_cov + alpha * sigma * f[plot_cell]
  for (i in seq_len(n_plots)) {
    if (a[i] == 1L) {
      yi <- stats::rbinom(J, 1L, p_plot[i])
      y[i, ] <- yi
      det <- which(yi == 1L)
      if (length(det)) {
        if (positive == "beta") {
          mu <- stats::plogis(eta_pos[i])
          y_pos[i, det] <- stats::rbeta(length(det), mu * phi, (1 - mu) * phi)
        } else {
          y_pos[i, det] <- stats::rlnorm(length(det), eta_pos[i], phi)
        }
      }
    }
  }

  data <- data.frame(cell = plot_cell, x_cell = x_cell[plot_cell],
                     x_plot = x_plot, x_pdet = x_pdet, x_cov = x_cov)

  list(
    y     = y,
    y_pos = y_pos,
    data  = data,
    adj   = adj,
    truth = list(beta_psi = beta_psi, beta_theta = beta_theta,
                 beta_p = beta_p, beta_pos = beta_pos,
                 positive = positive, phi = phi,
                 sigma = sigma, alpha = alpha,
                 f = f, psi = psi, z = z, theta = theta, a = a,
                 plot_cell = plot_cell)
  )
}
