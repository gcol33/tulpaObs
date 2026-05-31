# predict() for the joint cover-family fits: occu_cover() (3-arm) and the
# cover() hurdle on the nested-Laplace shared-field path (2-arm). One predict
# path serves both, driven by the family-agnostic draw bundle in
# joint_substrate.R (gcol33/tulpaObs#22, #23, #24).
#
# Returns a posterior-draws object with a tidy accessor -- the tulpaObs analogue
# of the hand-rolled INLA pred.inla.joint() collaborators write
# (inla.posterior.sample + control.predictor, then difference the draws). The
# deliverable is "plotting objects, not plots": hand back the per-unit draw
# matrices plus a plot-ready table; the user maps it.
#
# The joint latent (per-arm betas + shared field) is sampled from the grid-
# integrated posterior via tulpa::tulpa_posterior_draws() -- the faithful outer-
# grid mixture sum_k w_k N(m_k, V_k). Every derived quantity (p, mu, expected
# cover, and the change deltas) is computed PER DRAW and only then summarized --
# the "Marginalize Derived Quantities" rule -- so the nonlinear products and
# differences carry the joint posterior's correlation rather than a plug-in of
# posterior means.


# Build the per-arm design matrix at `newdata`, matching the fitted coefficient
# count `p_arm`. occu_cover() arms carry distinct occupancy / cover formulas
# (the cover arm padding any visit-level columns at the reference); the cover()
# hurdle shares one formula across both arms but autoscales each arm separately,
# so the design is rescaled with that arm's stored scaling.
.tobs_joint_arm_design <- function(object, newdata, arm, p_arm) {
  if (inherits(object, "cover_fit")) {
    enc   <- object$encoding
    scale <- if (identical(arm, "occ")) enc$scale_occ else enc$scale_pos
    X <- stats::model.matrix(enc$formula, newdata)
    if (ncol(X) != p_arm) {
      stop(sprintf(paste0(
        "predict(cover): the design from `newdata` has %d columns but the ",
        "fitted model has %d. Check that `newdata` carries every covariate in ",
        "the formula with matching factor levels."), ncol(X), p_arm),
        call. = FALSE)
    }
    return(.apply_scale_to_X(X, scale))
  }
  .occu_cover_arm_design(object$model, newdata, arm, p_arm)
}

# occu_cover() per-arm design at `newdata`. Occupancy ("occ") is cell-level and
# uses its whole design. The cover arm ("pos") is a cell-level (site) block plus,
# when the fit carried visit covariates, trailing visit-level columns held at the
# reference (0) for a per-cell prediction.
.occu_cover_arm_design <- function(model, newdata, arm, p_arm) {
  formula <- if (arm == "occ") model$formulas$occ else model$formulas$pos
  tt <- stats::delete.response(stats::terms(formula))
  X_site <- stats::model.matrix(tt, newdata)
  p_site <- if (arm == "occ") ncol(model$X_occ) else ncol(model$X_pos_site)
  if (ncol(X_site) != p_site) {
    stop(sprintf(paste0(
      "predict(occu_cover): the %s site design from `newdata` has %d ",
      "columns but the fitted model has %d. Check that `newdata` carries ",
      "every cell-level covariate in the %s formula with matching factor ",
      "levels."),
      arm, ncol(X_site), p_site, arm), call. = FALSE)
  }
  if (p_site == p_arm) return(X_site)
  cbind(X_site, matrix(0.0, nrow(X_site), p_arm - p_site))
}

# Posterior-summary table for a single quantity: per-unit mean + level CI.
.occu_cover_summ <- function(mat, cell, level) {
  a <- (1 - level) / 2
  qs <- t(apply(mat, 1L, stats::quantile, probs = c(a, 1 - a), names = FALSE))
  data.frame(cell = cell,
             mean = rowMeans(mat),
             sd   = apply(mat, 1L, stats::sd),
             lwr  = qs[, 1L],
             upr  = qs[, 2L])
}

# Core predict handler for the joint cover-family fits. `object` is an
# occu_cover() fit (3-arm) or a cover() hurdle fit on the nested-Laplace path
# (2-arm); both expose a joint nested-Laplace object via `.tobs_joint_fit()`.
.tobs_predict_joint <- function(object, newdata = NULL,
                                type = "occurrence", times = NULL,
                                level = 0.95, nsim = 1000L,
                                draws = TRUE, time_col = NULL) {
  type <- match.arg(type, c("occurrence", "cover_cond", "cover_exp", "change"))
  if (is.null(.tobs_joint_fit(object))) {
    stop("predict() requires a joint nested-Laplace fit (method = ",
         "\"nested_laplace\"); this fit carries no joint object.",
         call. = FALSE)
  }

  bundle  <- .tobs_joint_draws(object, n = nsim)
  n_cells <- bundle$n_cells
  n_field <- length(bundle$blocks)

  # Trend (time-varying field) fits weight blocks 2.. by a per-cell covariate.
  # A single `time_col` drives every trend block (and is set as the time column
  # in the change map); require it when the fit has any trend field.
  if (n_field > 1L) {
    if (is.null(time_col)) time_col <- object$trend_weight
    if (is.null(time_col)) {
      stop("predict(): this fit has ", n_field - 1L,
           " time-varying (trend) field(s); pass `time_col = \"<column>\"`, ",
           "the per-cell covariate that weights the trend field(s) (the same ",
           "column used at fit time via control$trend).", call. = FALSE)
    }
    for (b in 2:n_field) bundle$blocks[[b]]$weight <- time_col
  }

  if (is.null(newdata)) {
    newdata <- object$model$data
    if (is.null(newdata)) {
      stop("predict(): `newdata` is required (one row per spatial unit, or a ",
           "`cell` column indexing the field cells).", call. = FALSE)
    }
  }

  # cell map: explicit `cell` column, else row i -> field cell i.
  if (!is.null(newdata$cell)) {
    cell <- as.integer(newdata$cell)
  } else {
    cell <- seq_len(nrow(newdata))
  }
  if (anyNA(cell) || any(cell < 1L) || any(cell > n_cells)) {
    stop("predict(): `cell` must index field cells 1..", n_cells,
         " (add a `cell` column to `newdata`, or pass one row per field cell ",
         "in cell order).", call. = FALSE)
  }

  state <- function(nd) {
    X_occ <- .tobs_joint_arm_design(object, nd, "occ", ncol(bundle$b$occ))
    X_pos <- .tobs_joint_arm_design(object, nd, "pos", ncol(bundle$b$pos))
    wf <- function(nm) {
      if (!nm %in% names(nd)) {
        stop("predict(): trend-field weight column '", nm,
             "' is not in `newdata`.", call. = FALSE)
      }
      as.numeric(nd[[nm]])
    }
    eta_occ <- .tobs_joint_arm_eta(bundle, X_occ, "occ", cell, wf)
    eta_pos <- .tobs_joint_arm_eta(bundle, X_pos, "pos", cell, wf)
    p <- stats::plogis(eta_occ)
    mu <- if (identical(bundle$positive, "beta")) {
      stats::plogis(eta_pos)
    } else {
      exp(sweep(eta_pos, 2L, bundle$disp^2 / 2, "+"))
    }
    list(p = p, mu = mu, E = p * mu)
  }

  # --- single-time quantities ---------------------------------------------
  if (type != "change") {
    st <- state(newdata)
    mat <- switch(type,
                  occurrence = st$p,
                  cover_cond = st$mu,
                  cover_exp  = st$E)
    tbl <- .occu_cover_summ(mat, cell, level)
    attr(tbl, "quantity") <- type
    if (isTRUE(draws)) attr(tbl, "draws") <- stats::setNames(list(mat), type)
    class(tbl) <- c("tobs_prediction", "data.frame")
    return(tbl)
  }

  # --- change between two values of the time covariate ---------------------
  if (is.null(times) || length(times) != 2L) {
    stop("predict(type = \"change\") needs `times = c(t1, t2)`.", call. = FALSE)
  }
  if (is.null(time_col)) time_col <- object$trend_weight
  if (is.null(time_col)) {
    stop("predict(type = \"change\") needs the name of the time covariate. ",
         "Pass `time_col = \"<column>\"` (the covariate whose change between ",
         "`times[1]` and `times[2]` drives the prediction).", call. = FALSE)
  }
  if (!time_col %in% names(newdata)) {
    stop("predict(): `time_col = \"", time_col, "\"` is not a column of ",
         "`newdata`.", call. = FALSE)
  }

  nd1 <- newdata; nd1[[time_col]] <- times[1L]
  nd2 <- newdata; nd2[[time_col]] <- times[2L]
  s1 <- state(nd1)
  s2 <- state(nd2)

  # Per-draw deltas + the exact additive decomposition
  #   delta_exp = p2 mu2 - p1 mu1 = (p2 - p1) mu1 + p2 (mu2 - mu1).
  d_p    <- s2$p  - s1$p
  d_cond <- s2$mu - s1$mu
  d_exp  <- s2$E  - s1$E
  d_occ  <- (s2$p - s1$p) * s1$mu
  d_ab   <- s2$p * (s2$mu - s1$mu)

  a <- (1 - level) / 2
  qlo <- function(m) apply(m, 1L, stats::quantile, probs = a,     names = FALSE)
  qhi <- function(m) apply(m, 1L, stats::quantile, probs = 1 - a, names = FALSE)

  tbl <- data.frame(
    cell             = cell,
    p_T1             = rowMeans(s1$p),
    p_T2             = rowMeans(s2$p),
    delta_p          = rowMeans(d_p),
    cover_cond_T1    = rowMeans(s1$mu),
    cover_cond_T2    = rowMeans(s2$mu),
    delta_cover_cond = rowMeans(d_cond),
    cover_exp_T1     = rowMeans(s1$E),
    cover_exp_T2     = rowMeans(s2$E),
    delta_cover_exp  = rowMeans(d_exp),
    delta_cover_from_occ = rowMeans(d_occ),
    delta_cover_from_ab  = rowMeans(d_ab)
  )
  # .lwr / .upr for each delta at `level`.
  tbl$delta_p.lwr <- qlo(d_p);   tbl$delta_p.upr <- qhi(d_p)
  tbl$delta_cover_cond.lwr <- qlo(d_cond); tbl$delta_cover_cond.upr <- qhi(d_cond)
  tbl$delta_cover_exp.lwr  <- qlo(d_exp);  tbl$delta_cover_exp.upr  <- qhi(d_exp)
  tbl$delta_cover_from_occ.lwr <- qlo(d_occ); tbl$delta_cover_from_occ.upr <- qhi(d_occ)
  tbl$delta_cover_from_ab.lwr  <- qlo(d_ab);  tbl$delta_cover_from_ab.upr  <- qhi(d_ab)

  attr(tbl, "quantity") <- "change"
  attr(tbl, "times")    <- times
  if (isTRUE(draws)) {
    attr(tbl, "draws") <- list(
      p_T1 = s1$p, p_T2 = s2$p, delta_p = d_p,
      cover_cond_T1 = s1$mu, cover_cond_T2 = s2$mu, delta_cover_cond = d_cond,
      cover_exp_T1 = s1$E, cover_exp_T2 = s2$E, delta_cover_exp = d_exp,
      delta_cover_from_occ = d_occ, delta_cover_from_ab = d_ab
    )
  }
  class(tbl) <- c("tobs_prediction", "data.frame")
  tbl
}
