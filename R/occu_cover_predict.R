# predict() for the joint cover-family fits: occu_cover() (3-arm) and the
# cover() hurdle on the nested-Laplace shared-field path (2-arm). One predict
# path serves both, driven by the family-agnostic draw bundle in
# joint_substrate.R.
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


# Where one arm's fitted design lives on the fit: the arm's formula, the
# cell-level (site) coefficient count the fit carried, the visit-level formula
# whose columns trail the site block (NULL when the arm carries none), and the
# front-door name the width-mismatch message reports. occu_cover() stores a
# per-arm site design plus its visit formulas; the standalone occu() SVC fit
# stores the occupancy / detection designs in `X_processes` and holds any
# visit-level detection columns at their reference.
.tobs_joint_arm_slots <- function(object, arm) {
  model <- object$model
  if (isTRUE(object$occu_only_joint)) {
    occ <- identical(arm, "occ")
    return(list(
      formula = if (occ) model$formulas$occ else model$formulas$det,
      p_site  = ncol(model$X_processes[[if (occ) 1L else 2L]]),
      visit   = NULL,
      label   = "occu"))
  }
  list(
    formula = switch(arm, occ = model$formulas$occ,
                          det = model$formulas$det,
                          pos = model$formulas$pos),
    p_site  = switch(arm, occ = ncol(model$X_occ),
                          det = ncol(model$X_det_site),
                          pos = ncol(model$X_pos_site)),
    visit   = switch(arm, det = model$formulas$det_visit,
                          pos = model$formulas$pos_visit, NULL),
    label   = "occu_cover")
}

# Build the per-arm design matrix at `newdata`, matching the fitted coefficient
# count `p_arm`. The cover() hurdle shares one formula across both arms but
# autoscales each arm separately, so its design is rescaled with that arm's
# stored scaling. Every other arm takes its own formula: the cell-level block is
# checked against the width the fit carried, and any remaining coefficients are
# visit-level covariates (e.g. the time axis). When the arm carries a visit
# formula those are rebuilt from `newdata` -- one prediction row per cell -- with
# the same builder, reference (k - 1) coding, and column order as the fit, so a
# covariate supplied in `newdata` enters the linear predictor (the change map
# varies the time covariate this way). A covariate ABSENT from `newdata`, and
# every trailing column on an arm carrying no visit formula, is held at its
# reference (numeric 0 / factor base level via the zero columns).
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
  slots  <- .tobs_joint_arm_slots(object, arm)
  X_site <- .tobs_arm_site_design(slots$formula, newdata, slots$p_site, arm,
                                  slots$label)
  if (slots$p_site == p_arm) return(X_site)
  if (!is.null(slots$visit)) {
    nd <- newdata
    for (v in setdiff(all.vars(slots$visit), names(nd))) nd[[v]] <- 0
    X_visit <- tryCatch(
      .tobs_build_visit_X(slots$visit, nd, nrow(nd), 1L, arm = arm),
      error = function(e) NULL)
    if (!is.null(X_visit) && ncol(X_site) + ncol(X_visit) == p_arm)
      return(cbind(X_site, X_visit))
  }
  cbind(X_site, matrix(0.0, nrow(X_site), p_arm - slots$p_site))
}

# Per-arm RE BLUP offset at `newdata`: an [n_rows x n_draws] matrix added to the
# arm's linear predictor, SUMMED over every RE term on that arm (crossed / nested
# groupings each contribute). For each term, a row's grouping level (the term's
# grouping variable, e.g. `habitat`) is matched against the fitted levels and the
# matching group's latent draws are added; an UNSEEN level (no match) or a
# `newdata` without the grouping column shrinks that term to the population mean
# (offset 0). The draws come from the grid-integrated posterior so the offset is
# marginalised over the joint, not a plug-in of the BLUP mean.
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
    # One-inflated Beta: the conditional cover mixes the constant ceiling mass
    # (cover = 1) with the interior Beta mean, E[cover | y > 0] = pi + (1 -
    # pi) * plogis(eta_pos). cover()-only, so no latent cover-aggregate
    # variant.
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
  if (identical(bundle$positive, "gaussian")) {
    # Identity-Gaussian arm: the response IS the mean, so the conditional cover
    # mean is the linear predictor itself, mu = eta. No latent cover-aggregate
    # variant (the latent path is lognormal / beta only).
    return(eta_pos)
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

# RE arm key per bundle arm slot: a fit records its random effects under the
# predictor name ("psi" occupancy, "p" detection, "pos" cover) while the draw
# bundle keys the coefficient blocks by arm slot.
.TOBS_JOINT_RE_ARM <- c(occ = "psi", det = "p", pos = "pos")

# Point the fit's shared trend (time-varying) field blocks at the prediction's
# time column and return the bundle. A shared trend block is weighted by a
# per-cell covariate and a single `time_col` drives every one of them (it is also
# the change map's time column), so it is required as soon as the fit carries
# one. `field_specs` labels each block's arm and weight; a fit that carries none
# (the standalone occu() route, older cover fits) follows the ordering the field
# assembly guarantees, block 1 the unweighted intercept field and blocks 2.. the
# trend fields. cover() arm-specific fits and occu_cover() arm-specific cover
# fields carry their own per-block weight column and keep it -- an intercept-only
# arm-specific field needs no time_col; a pos-arm trend field resolves its weight
# from newdata.
.tobs_joint_trend_blocks <- function(object, bundle, time_col) {
  n_field <- length(bundle$blocks)
  if (n_field <= 1L || isTRUE(object$armspecific)) return(bundle)
  specs <- object$field_specs
  shared_trend <- if (!is.null(specs) && length(specs) >= n_field) {
    vapply(seq_len(n_field), function(b)
      identical(specs[[b]]$arm, "shared") && !is.null(specs[[b]]$weight),
      logical(1))
  } else {
    c(FALSE, rep(TRUE, n_field - 1L))
  }
  if (!any(shared_trend)) return(bundle)
  if (is.null(time_col)) time_col <- object$trend_weight
  if (is.null(time_col)) {
    stop("predict(): this fit has ", sum(shared_trend),
         " shared time-varying (trend) field(s); pass `time_col = ",
         "\"<column>\"`, the per-cell covariate that weights the trend ",
         "field(s) (the same column used at fit time).", call. = FALSE)
  }
  for (b in which(shared_trend)) bundle$blocks[[b]]$weight <- time_col
  bundle
}

# Prediction frame and its cell map. `newdata` defaults to the fit's own data
# (one row per spatial unit); an explicit `cell` column indexes the field cells,
# otherwise row i is field cell i.
.tobs_joint_predict_cells <- function(object, newdata, n_cells) {
  if (is.null(newdata)) {
    newdata <- object$model$data
    if (is.null(newdata)) {
      stop("predict(): `newdata` is required (one row per spatial unit, or a ",
           "`cell` column indexing the field cells).", call. = FALSE)
    }
  }
  cell <- if (!is.null(newdata$cell)) as.integer(newdata$cell)
          else seq_len(nrow(newdata))
  if (anyNA(cell) || any(cell < 1L) || any(cell > n_cells)) {
    stop("predict(): `cell` must index field cells 1..", n_cells,
         " (add a `cell` column to `newdata`, or pass one row per field cell ",
         "in cell order).", call. = FALSE)
  }
  list(newdata = newdata, cell = cell)
}

# Per-draw quantities at `nd` for the arms the fit carries. `bundle$b` is the
# roster: the occupancy arm ("occ") is always present, the detection arm ("det")
# on occu() and occu_cover() fits, the cover arm ("pos") on occu_cover() and
# cover() fits. `arms` narrows the evaluation to the arms the caller reads, so a
# quantity is never charged the design of an arm it does not use. Returns `p`
# (occupancy), `p_det` (detection), `mu` (conditional cover) and `E` (expected
# cover) for whichever arms were evaluated.
#
# Each arm's linear predictor is the fitted betas, plus the shared-field
# contribution, plus that arm's RE BLUP offset: the group's latent draws are
# added when the fit carries an RE on the arm AND `newdata` supplies the grouping
# column; otherwise the term shrinks to the population mean (offset 0), the
# field-only behaviour.
.tobs_joint_arm_states <- function(object, bundle, nd, cell, arms = NULL) {
  present <- names(bundle$b)[!vapply(bundle$b, is.null, logical(1))]
  arms    <- if (is.null(arms)) present else intersect(arms, present)
  wf <- function(nm) {
    if (!nm %in% names(nd)) {
      stop("predict(): trend-field weight column '", nm,
           "' is not in `newdata`.", call. = FALSE)
    }
    as.numeric(nd[[nm]])
  }
  eta <- function(arm) {
    X <- .tobs_joint_arm_design(object, nd, arm, ncol(bundle$b[[arm]]))
    .tobs_joint_arm_eta(bundle, X, arm, cell, wf) +
      .occu_cover_re_offset(bundle, .TOBS_JOINT_RE_ARM[[arm]], nd, nrow(X),
                            bundle$n)
  }
  out <- list()
  if ("occ" %in% arms) out$p     <- stats::plogis(eta("occ"))
  if ("det" %in% arms) out$p_det <- stats::plogis(eta("det"))
  if ("pos" %in% arms) {
    out$mu <- .tobs_cover_mu(eta("pos"), bundle, object)
    if (!is.null(out$p)) out$E <- out$p * out$mu
  }
  out
}

# The prediction frames a change map differences: `newdata` with the time
# covariate held at each of `times`. Two times give one difference; K times give
# a trajectory, every step differenced against `times[1]`.
.tobs_joint_change_frames <- function(object, newdata, times, time_col) {
  if (is.null(times) || length(times) < 2L) {
    stop("predict(type = \"change\") needs `times = c(t1, t2)` (one change) or ",
         "`times = c(t1, ..., tK)` (a trajectory, each step against t1).",
         call. = FALSE)
  }
  if (!is.numeric(times) || anyNA(times) || !all(is.finite(times))) {
    stop("predict(type = \"change\"): `times` must be finite numeric values.",
         call. = FALSE)
  }
  if (is.null(time_col)) time_col <- object$trend_weight
  if (is.null(time_col)) {
    stop("predict(type = \"change\") needs the name of the time covariate. ",
         "Pass `time_col = \"<column>\"` (the covariate the `times` are held ",
         "at, whose movement drives the prediction).", call. = FALSE)
  }
  if (!time_col %in% names(newdata)) {
    stop("predict(): `time_col = \"", time_col, "\"` is not a column of ",
         "`newdata`.", call. = FALSE)
  }
  lapply(times, function(tk) { nd <- newdata; nd[[time_col]] <- tk; nd })
}

# Per-row summaries over a draw matrix: the `level` central-interval endpoints
# and the posterior SD.
.tobs_draw_lwr <- function(m, level)
  apply(m, 1L, stats::quantile, probs = (1 - level) / 2, names = FALSE)
.tobs_draw_upr <- function(m, level)
  apply(m, 1L, stats::quantile, probs = 1 - (1 - level) / 2, names = FALSE)
.tobs_draw_sd <- function(m) apply(m, 1L, stats::sd)

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

  bundle <- .tobs_joint_trend_blocks(
    object, .tobs_joint_draws(object, n = nsim), time_col)
  nc      <- .tobs_joint_predict_cells(object, newdata, bundle$n_cells)
  newdata <- nc$newdata
  cell    <- nc$cell

  state <- function(nd) .tobs_joint_arm_states(object, bundle, nd, cell)

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

  # --- the time covariate held at each of `times` --------------------------
  # Two times give one difference; K times give a trajectory, every step
  # differenced against the first. One draw set serves the whole table, so the
  # steps share a posterior and their deltas are jointly valid.
  nds <- .tobs_joint_change_frames(object, newdata, times, time_col)
  st  <- lapply(nds, state)
  K   <- length(st)
  s1  <- st[[1L]]

  # Per-draw deltas against the baseline step + the exact additive decomposition
  #   delta_exp = pk muk - p1 mu1 = (pk - p1) mu1 + pk (muk - mu1).
  d <- lapply(st[-1L], function(sk) list(
    p    = sk$p  - s1$p,
    cond = sk$mu - s1$mu,
    exp  = sk$E  - s1$E,
    occ  = (sk$p - s1$p) * s1$mu,
    ab   = sk$p * (sk$mu - s1$mu)))

  qlo <- function(m) .tobs_draw_lwr(m, level)
  qhi <- function(m) .tobs_draw_upr(m, level)
  qsd <- .tobs_draw_sd

  # A delta belongs to the step it was taken at, and with two times there is
  # only one step -- so the delta columns carry a `_T<k>` suffix only when the
  # table holds a trajectory. Two times therefore keep the column names (and
  # their order) a change map has always had.
  sfx <- if (K == 2L) rep("", K) else paste0("_T", seq_len(K))
  lvl <- function(nm, k) sprintf("%s_T%d", nm, k)
  dl  <- function(nm, k) paste0(nm, sfx[k])

  tbl <- data.frame(cell = cell)
  add <- function(tbl, nm, v) { tbl[[nm]] <- v; tbl }

  # Quantity-major: every step's level, then that quantity's deltas.
  for (k in seq_len(K)) tbl <- add(tbl, lvl("p", k), rowMeans(st[[k]]$p))
  for (k in 2:K)        tbl <- add(tbl, dl("delta_p", k), rowMeans(d[[k - 1L]]$p))
  for (k in seq_len(K)) tbl <- add(tbl, lvl("cover_cond", k), rowMeans(st[[k]]$mu))
  for (k in 2:K)        tbl <- add(tbl, dl("delta_cover_cond", k),
                                   rowMeans(d[[k - 1L]]$cond))
  for (k in seq_len(K)) tbl <- add(tbl, lvl("cover_exp", k), rowMeans(st[[k]]$E))
  for (k in 2:K)        tbl <- add(tbl, dl("delta_cover_exp", k),
                                   rowMeans(d[[k - 1L]]$exp))
  for (k in 2:K)        tbl <- add(tbl, dl("delta_cover_from_occ", k),
                                   rowMeans(d[[k - 1L]]$occ))
  for (k in 2:K)        tbl <- add(tbl, dl("delta_cover_from_ab", k),
                                   rowMeans(d[[k - 1L]]$ab))

  # .lwr / .upr for each delta at `level`.
  for (q in c("p", "cover_cond", "cover_exp", "cover_from_occ", "cover_from_ab")) {
    slot <- switch(q, p = "p", cover_cond = "cond", cover_exp = "exp",
                   cover_from_occ = "occ", cover_from_ab = "ab")
    for (k in 2:K) {
      nm <- dl(paste0("delta_", q), k)
      m  <- d[[k - 1L]][[slot]]
      tbl <- add(tbl, paste0(nm, ".lwr"), qlo(m))
      tbl <- add(tbl, paste0(nm, ".upr"), qhi(m))
    }
  }

  # Per-step occupancy uncertainty (sd + CI at `level`) and the directional
  # posterior probability P(delta > 0) per cell -- the certainty that the quantity
  # increased. Both are taken over draws, so they carry the joint posterior rather
  # than a plug-in of the means (the marginalize-derived-quantities rule). These
  # are the per-cell change-certainty quantities a spatially-varying-trend
  # occupancy fit reports; they are pure additions to the table.
  for (k in seq_len(K)) {
    nm <- lvl("p", k)
    tbl <- add(tbl, paste0(nm, ".sd"),  qsd(st[[k]]$p))
    tbl <- add(tbl, paste0(nm, ".lwr"), qlo(st[[k]]$p))
    tbl <- add(tbl, paste0(nm, ".upr"), qhi(st[[k]]$p))
  }
  for (q in c("p", "cover_cond", "cover_exp")) {
    slot <- switch(q, p = "p", cover_cond = "cond", cover_exp = "exp")
    for (k in 2:K) {
      tbl <- add(tbl, paste0(dl(paste0("delta_", q), k), ".prob_pos"),
                 rowMeans(d[[k - 1L]][[slot]] > 0))
    }
  }

  attr(tbl, "quantity") <- "change"
  attr(tbl, "times")    <- times
  if (isTRUE(draws)) {
    dr <- list()
    for (k in seq_len(K)) dr[[lvl("p", k)]] <- st[[k]]$p
    for (k in 2:K)        dr[[dl("delta_p", k)]] <- d[[k - 1L]]$p
    for (k in seq_len(K)) dr[[lvl("cover_cond", k)]] <- st[[k]]$mu
    for (k in 2:K)        dr[[dl("delta_cover_cond", k)]] <- d[[k - 1L]]$cond
    for (k in seq_len(K)) dr[[lvl("cover_exp", k)]] <- st[[k]]$E
    for (k in 2:K)        dr[[dl("delta_cover_exp", k)]] <- d[[k - 1L]]$exp
    for (k in 2:K)        dr[[dl("delta_cover_from_occ", k)]] <- d[[k - 1L]]$occ
    for (k in 2:K)        dr[[dl("delta_cover_from_ab", k)]] <- d[[k - 1L]]$ab
    attr(tbl, "draws") <- dr
  }
  class(tbl) <- c("tobs_prediction", "data.frame")
  tbl
}


# Site-level design at `newdata` for one arm, checked against the width the fit
# carried. A mismatch means `newdata` is missing a covariate, or carries a factor
# with different levels, which would otherwise shift every coefficient onto the
# wrong column silently -- so it errors naming the arm and both widths. `label`
# names the front door in the message.
.tobs_arm_site_design <- function(formula, newdata, p_site, arm, label) {
  tt     <- stats::delete.response(stats::terms(formula))
  X_site <- stats::model.matrix(tt, newdata)
  if (ncol(X_site) != p_site) {
    stop(sprintf(paste0(
      "predict(%s): the %s site design from `newdata` has %d columns but the ",
      "fitted model has %d. Check that `newdata` carries every cell-level ",
      "covariate in the %s formula with matching factor levels."),
      label, arm, ncol(X_site), p_site, arm), call. = FALSE)
  }
  X_site
}

# Core predict handler for the rerouted standalone occu() SVC joint fit. The
# occupancy psi and detection p are computed PER DRAW from the grid-integrated
# joint posterior (betas + shared field + the arm's RE BLUP offset) and only
# then summarized -- the marginalize-derived-quantities rule, so the per-cell
# psi carries the joint posterior's correlation, including the spatial field at
# each cell. Supports `type = "occupancy" | "detection" | "both" | "change"`;
# "change" reads the occupancy difference between two values of the trend-field
# weight covariate.
.tobs_predict_occu_joint <- function(object, newdata = NULL,
                                     type = "occupancy", times = NULL,
                                     level = 0.95, nsim = 1000L,
                                     draws = TRUE, time_col = NULL) {
  type <- match.arg(type, c("occupancy", "detection", "both", "change"))
  if (is.null(.tobs_joint_fit(object))) {
    stop("predict() requires a joint nested-Laplace fit; this occu() fit ",
         "carries no joint object.", call. = FALSE)
  }

  bundle <- .tobs_joint_trend_blocks(
    object, .tobs_joint_draws(object, n = nsim), time_col)
  nc      <- .tobs_joint_predict_cells(object, newdata, bundle$n_cells)
  newdata <- nc$newdata
  cell    <- nc$cell

  occ_state <- function(nd)
    .tobs_joint_arm_states(object, bundle, nd, cell, "occ")$p
  det_state <- function(nd)
    .tobs_joint_arm_states(object, bundle, nd, cell, "det")$p_det

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
  nds <- .tobs_joint_change_frames(object, newdata, times, time_col)
  p1  <- occ_state(nds[[1L]])
  p2  <- occ_state(nds[[2L]])
  d_p <- p2 - p1
  qlo <- function(m) .tobs_draw_lwr(m, level)
  qhi <- function(m) .tobs_draw_upr(m, level)
  qsd <- .tobs_draw_sd
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
