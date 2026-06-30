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
  formula <- switch(arm, occ = model$formulas$occ,
                         det = model$formulas$det,
                         pos = model$formulas$pos)
  tt <- stats::delete.response(stats::terms(formula))
  X_site <- stats::model.matrix(tt, newdata)
  p_site <- switch(arm, occ = ncol(model$X_occ),
                        det = ncol(model$X_det_site),
                        pos = ncol(model$X_pos_site))
  if (ncol(X_site) != p_site) {
    stop(sprintf(paste0(
      "predict(occu_cover): the %s site design from `newdata` has %d ",
      "columns but the fitted model has %d. Check that `newdata` carries ",
      "every cell-level covariate in the %s formula with matching factor ",
      "levels."),
      arm, ncol(X_site), p_site, arm), call. = FALSE)
  }
  if (p_site == p_arm) return(X_site)
  # The remaining coefficients are visit-level detection / positive covariates
  # (e.g. the time axis). Rebuild their design from `newdata` -- one prediction
  # row per cell -- with the same builder, reference (k - 1) coding, and column
  # order as the fit, so a covariate supplied in `newdata` enters the linear
  # predictor (the change map varies the time covariate this way). A covariate
  # ABSENT from `newdata` is held at its reference (numeric 0 / factor base level
  # via the zero columns), the long-standing per-cell default.
  vf <- switch(arm, det = model$formulas$det_visit,
                    pos = model$formulas$pos_visit, NULL)
  if (!is.null(vf)) {
    nd <- newdata
    for (v in setdiff(all.vars(vf), names(nd))) nd[[v]] <- 0
    X_visit <- tryCatch(
      .tobs_build_visit_X(vf, nd, nrow(nd), 1L, arm = arm),
      error = function(e) NULL)
    if (!is.null(X_visit) && ncol(X_site) + ncol(X_visit) == p_arm)
      return(cbind(X_site, X_visit))
  }
  cbind(X_site, matrix(0.0, nrow(X_site), p_arm - p_site))
}

# Per-arm RE BLUP offset at `newdata` (gcol33/tulpaObs#102, #103): an
# [n_rows x n_draws] matrix added to the arm's linear predictor, SUMMED over
# every RE term on that arm (crossed / nested groupings each contribute). For
# each term, a row's grouping level (the term's grouping variable, e.g.
# `habitat`) is matched against the fitted levels and the matching group's latent
# draws are added; an UNSEEN level (no match) or a `newdata` without the grouping
# column shrinks that term to the population mean (offset 0). The draws come from
# the grid-integrated posterior so the offset is marginalised over the joint, not
# a plug-in of the BLUP mean.
.occu_cover_re_offset <- function(bundle, arm, nd, n_rows, n_draws) {
  off   <- matrix(0, n_rows, n_draws)
  terms <- Filter(function(r) identical(r$arm, arm), bundle$re %||% list())
  for (re in terms) {
    if (is.null(re$var) || !(re$var %in% names(nd))) next
    codes <- match(as.character(nd[[re$var]]), re$levels)
    seen  <- which(!is.na(codes))
    if (!length(seen)) next
    nc <- re$n_coefs %||% 1L
    if (nc == 1L) {
      # Random intercept: add the group's per-draw offset (draws is [n x n_groups]).
      off[seen, ] <- off[seen, ] + t(re$draws[, codes[seen], drop = FALSE])
    } else {
      # Random slope: draws is coefficient-major [n x (n_coefs * n_groups)]; the
      # row offset is sum_c w_c[row] * b_draw[group, c], the intercept weight 1 and
      # each slope weight the covariate value from `newdata` (0 / absent column ->
      # that coefficient drops out, the same shrink as an unseen level).
      ng <- re$n_groups
      for (cc in seq_len(nc)) {
        nm <- re$coef_names[cc]
        w  <- if (identical(nm, "(Intercept)")) rep(1, length(seen))
              else if (nm %in% names(nd)) as.numeric(nd[[nm]][seen])
              else next
        cols <- (cc - 1L) * ng + codes[seen]
        off[seen, ] <- off[seen, ] + w * t(re$draws[, cols, drop = FALSE])
      }
    }
  }
  off
}

# Cover mean on the natural scale from the cover linear predictor `eta_pos`
# ([nrow x n]). `bundle$disp` is the per-draw positive-arm dispersion read off
# the `phi_pos` grid axis: the residual SD (lognormal) / precision (beta) for the
# non-latent paths, the integrated cover-latent SD sigma_u for the latent path.
#
# Non-latent: the per-visit cover mean. Lognormal (and ordinal interval-censored
# Gaussian, its log-cover sibling) log-normal mean exp(eta + sigma^2/2); beta
# mean plogis(eta) (the beta mean given the linear predictor).
#
# Latent (cover_aggregate == "latent"): the cover mean marginalized over the
# per-unit cover latent u ~ N(0, sigma_u^2), with the within-unit dispersion
# (sigma_eps lognormal / precision beta) held fixed at model$cover_latent_disp2.
# Lognormal: total log-scale variance sigma_eps^2 + sigma_u^2, so the marginal
# mean is exp(eta + (sigma_eps^2 + sigma_u^2)/2). Beta: the marginal mean is
# E_u[plogis(eta + u)] (no closed form), integrated by Gauss-Hermite. sigma_u
# rides the grid per draw, so the marginal cover is integrated over the joint
# outer grid -- the marginalize-derived-quantities rule, not a plug-in of the
# posterior-mean sigma_u.
.tobs_cover_mu <- function(eta_pos, bundle, object) {
  latent <- identical(object$model$cover_aggregate, "latent")
  if (identical(bundle$positive, "beta_oi")) {
    # One-inflated Beta (gcol33/tulpaObs#108): the conditional cover mixes the
    # constant ceiling mass (cover = 1) with the interior Beta mean,
    # E[cover | y > 0] = pi + (1 - pi) * plogis(eta_pos). cover()-only, so no
    # latent cover-aggregate variant.
    pi1 <- object$pi_one %||% 0
    return(pi1 + (1 - pi1) * stats::plogis(eta_pos))
  }
  if (identical(bundle$positive, "beta")) {
    if (!latent) return(stats::plogis(eta_pos))
    sigma_u <- bundle$disp
    gh <- .gauss_hermite(15L)
    c0 <- 1 / sqrt(pi)
    mu <- matrix(0, nrow(eta_pos), ncol(eta_pos))
    for (g in seq_along(gh$x)) {
      shift <- sqrt(2) * sigma_u * gh$x[g]                # per-draw u node
      mu <- mu + (gh$w[g] * c0) *
        stats::plogis(sweep(eta_pos, 2L, shift, "+"))
    }
    return(mu)
  }
  if (identical(bundle$positive, "lognormal_trunc")) {
    # Truncated-lognormal conditional cover mean: E[exp(t) | t <= u] with
    # t ~ N(eta, sg^2) and u = log(1) = 0 (cover <= 1) =
    # exp(eta + sg^2/2) * Phi((u - eta - sg^2)/sg) / Phi((u - eta)/sg). Per-draw sg
    # sweeps over columns. lognormal_trunc has no `latent` cover-aggregate variant.
    sg <- bundle$disp
    za <- sweep(0 - eta_pos, 2L, sg, "/")           # (u - eta)/sg, u = 0
    za_m <- sweep(za, 2L, sg, "-")                  # za - sg
    mean_log <- exp(sweep(eta_pos, 2L, sg^2 / 2, "+"))
    return(mean_log * stats::pnorm(za_m) / stats::pnorm(za))
  }
  if (!latent) return(exp(sweep(eta_pos, 2L, bundle$disp^2 / 2, "+")))
  sigma_eps <- as.numeric(object$model$cover_latent_disp2)
  log_var   <- sigma_eps^2 + bundle$disp^2               # per-draw total variance
  exp(sweep(eta_pos, 2L, log_var / 2, "+"))
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
  type <- match.arg(type, c("occurrence", "detection", "cover_cond",
                            "cover_exp", "change"))
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
  # in the change map); require it when the fit has any trend field. Arm-specific
  # fits (gcol33/tulpaObs#65) instead carry their own per-block weight column name
  # and may stack several intercept-only fields, so the positional 2.. convention
  # does not apply -- their weights resolve directly from newdata via `wf`, and a
  # purely intercept arm-specific fit needs no time_col (gcol33/tulpaObs#95).
  if (n_field > 1L && !isTRUE(object$armspecific)) {
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
    # Per-arm RE BLUP offset (gcol33/tulpaObs#102): added when the fit carries an
    # RE on that arm AND `newdata` supplies the grouping column; otherwise the
    # prediction is at the population mean (offset 0), the field-only behaviour.
    eta_occ <- .tobs_joint_arm_eta(bundle, X_occ, "occ", cell, wf) +
               .occu_cover_re_offset(bundle, "psi", nd, nrow(X_occ), bundle$n)
    eta_pos <- .tobs_joint_arm_eta(bundle, X_pos, "pos", cell, wf) +
               .occu_cover_re_offset(bundle, "pos", nd, nrow(X_pos), bundle$n)
    p <- stats::plogis(eta_occ)
    mu <- .tobs_cover_mu(eta_pos, bundle, object)
    out <- list(p = p, mu = mu, E = p * mu)
    if (!is.null(bundle$b$det)) {
      X_det   <- .tobs_joint_arm_design(object, nd, "det", ncol(bundle$b$det))
      eta_det <- .tobs_joint_arm_eta(bundle, X_det, "det", cell) +
                 .occu_cover_re_offset(bundle, "p", nd, nrow(X_det), bundle$n)
      out$p_det <- stats::plogis(eta_det)
    }
    out
  }

  # --- single-time quantities ---------------------------------------------
  if (type != "change") {
    st <- state(newdata)
    if (identical(type, "detection") && is.null(st$p_det)) {
      stop("predict(type = \"detection\") needs a detection arm; this is a ",
           "cover() hurdle fit (occupancy + cover only). Use occu_cover() for a ",
           "detection prediction.", call. = FALSE)
    }
    mat <- switch(type,
                  occurrence = st$p,
                  detection  = st$p_det,
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

  # Start / end occupancy uncertainty (sd + CI at `level`) and the directional
  # posterior probability P(delta > 0) per cell -- the certainty that the quantity
  # increased. Both are taken over draws, so they carry the joint posterior rather
  # than a plug-in of the means (the marginalize-derived-quantities rule). These
  # are the per-cell change-certainty quantities a spatially-varying-trend
  # occupancy fit reports; they are pure additions to the table.
  qsd <- function(m) apply(m, 1L, stats::sd)
  tbl$p_T1.sd <- qsd(s1$p); tbl$p_T1.lwr <- qlo(s1$p); tbl$p_T1.upr <- qhi(s1$p)
  tbl$p_T2.sd <- qsd(s2$p); tbl$p_T2.lwr <- qlo(s2$p); tbl$p_T2.upr <- qhi(s2$p)
  tbl$delta_p.prob_pos          <- rowMeans(d_p    > 0)
  tbl$delta_cover_cond.prob_pos <- rowMeans(d_cond > 0)
  tbl$delta_cover_exp.prob_pos  <- rowMeans(d_exp  > 0)

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


# Occupancy ("occ") arm design at `newdata` for the rerouted standalone occu()
# SVC joint fit. The single-season occu() model carries the occupancy formula at
# `model$formulas$occ` (the fixed-effects psi formula the spatial term was
# stripped from) and the fitted design at `model$X_processes[[1]]`; the field is
# added separately by `.tobs_joint_arm_eta`, not by this design. Detection
# ("det") uses the detection formula / `model$X_processes[[2]]`; visit-level
# columns (if any) are held at the reference (0) for a per-cell prediction.
.occu_joint_arm_design <- function(model, newdata, arm) {
  if (identical(arm, "occ")) {
    formula <- model$formulas$occ
    p_arm   <- ncol(model$X_processes[[1L]])
  } else {
    formula <- model$formulas$det
    p_arm   <- ncol(model$X_processes[[2L]]) +
               (if (is.null(model$X_det_visit)) 0L else ncol(model$X_det_visit))
  }
  tt <- stats::delete.response(stats::terms(formula))
  X_site <- stats::model.matrix(tt, newdata)
  p_site <- if (identical(arm, "occ")) ncol(model$X_processes[[1L]])
            else ncol(model$X_processes[[2L]])
  if (ncol(X_site) != p_site) {
    stop(sprintf(paste0(
      "predict(occu): the %s design from `newdata` has %d columns but the ",
      "fitted model has %d. Check that `newdata` carries every covariate in ",
      "the %s formula with matching factor levels."),
      arm, ncol(X_site), p_site, arm), call. = FALSE)
  }
  if (p_site == p_arm) return(X_site)
  cbind(X_site, matrix(0.0, nrow(X_site), p_arm - p_site))
}

# Core predict handler for the rerouted standalone occu() SVC joint fit
# (gcol33/tulpaObs#81): the occupancy-only twin of `.tobs_predict_joint`. The
# occupancy psi and detection p are computed PER DRAW from the grid-integrated
# joint posterior (betas + shared field) and only then summarized -- the
# marginalize-derived-quantities rule, so the per-cell psi carries the joint
# posterior's correlation, including the spatial field at each cell. Supports
# `type = "occupancy" | "detection" | "both" | "change"`; "change" reads the
# occupancy difference between two values of the trend-field weight covariate.
.tobs_predict_occu_joint <- function(object, newdata = NULL,
                                     type = "occupancy", times = NULL,
                                     level = 0.95, nsim = 1000L,
                                     draws = TRUE, time_col = NULL) {
  type <- match.arg(type, c("occupancy", "detection", "both", "change"))
  if (is.null(.tobs_joint_fit(object))) {
    stop("predict() requires a joint nested-Laplace fit; this occu() fit ",
         "carries no joint object.", call. = FALSE)
  }

  bundle  <- .tobs_joint_draws(object, n = nsim)
  n_cells <- bundle$n_cells
  n_field <- length(bundle$blocks)

  if (n_field > 1L) {
    if (is.null(time_col)) time_col <- object$trend_weight
    if (is.null(time_col)) {
      stop("predict(): this fit has ", n_field - 1L,
           " time-varying (trend) field(s); pass `time_col = \"<column>\"`, ",
           "the per-cell covariate that weights the trend field(s).",
           call. = FALSE)
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

  occ_state <- function(nd) {
    X_occ <- .occu_joint_arm_design(object$model, nd, "occ")
    wf <- function(nm) {
      if (!nm %in% names(nd)) {
        stop("predict(): trend-field weight column '", nm,
             "' is not in `newdata`.", call. = FALSE)
      }
      as.numeric(nd[[nm]])
    }
    eta_occ <- .tobs_joint_arm_eta(bundle, X_occ, "occ", cell, wf)
    stats::plogis(eta_occ)
  }
  det_state <- function(nd) {
    X_det <- .occu_joint_arm_design(object$model, nd, "det")
    eta_det <- .tobs_joint_arm_eta(bundle, X_det, "det", cell)
    stats::plogis(eta_det)
  }

  if (type %in% c("occupancy", "detection", "both")) {
    out <- list()
    if (type %in% c("occupancy", "both")) {
      psi <- occ_state(newdata)
      tbl <- .occu_cover_summ(psi, cell, level)
      attr(tbl, "quantity") <- "occupancy"
      if (isTRUE(draws)) attr(tbl, "draws") <- list(occupancy = psi)
      class(tbl) <- c("tobs_prediction", "data.frame")
      if (identical(type, "occupancy")) return(tbl)
      out$occupancy <- tbl
    }
    if (type %in% c("detection", "both")) {
      pp  <- det_state(newdata)
      tbl <- .occu_cover_summ(pp, cell, level)
      attr(tbl, "quantity") <- "detection"
      if (isTRUE(draws)) attr(tbl, "draws") <- list(detection = pp)
      class(tbl) <- c("tobs_prediction", "data.frame")
      if (identical(type, "detection")) return(tbl)
      out$detection <- tbl
    }
    return(out)
  }

  # type == "change": occupancy difference between two trend-weight values.
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
  p1 <- occ_state(nd1)
  p2 <- occ_state(nd2)
  d_p <- p2 - p1
  a <- (1 - level) / 2
  qlo <- function(m) apply(m, 1L, stats::quantile, probs = a,     names = FALSE)
  qhi <- function(m) apply(m, 1L, stats::quantile, probs = 1 - a, names = FALSE)
  qsd <- function(m) apply(m, 1L, stats::sd)
  # Start / end occupancy uncertainty and the directional posterior probability
  # P(delta > 0), both over draws (marginalize-derived-quantities); pure additions
  # alongside delta_psi.
  tbl <- data.frame(
    cell      = cell,
    psi_T1    = rowMeans(p1),
    psi_T1.sd = qsd(p1), psi_T1.lwr = qlo(p1), psi_T1.upr = qhi(p1),
    psi_T2    = rowMeans(p2),
    psi_T2.sd = qsd(p2), psi_T2.lwr = qlo(p2), psi_T2.upr = qhi(p2),
    delta_psi = rowMeans(d_p),
    delta_psi.lwr = qlo(d_p),
    delta_psi.upr = qhi(d_p),
    delta_psi.prob_pos = rowMeans(d_p > 0))
  attr(tbl, "quantity") <- "change"
  attr(tbl, "times")    <- times
  if (isTRUE(draws)) {
    attr(tbl, "draws") <- list(psi_T1 = p1, psi_T2 = p2, delta_psi = d_p)
  }
  class(tbl) <- c("tobs_prediction", "data.frame")
  tbl
}
