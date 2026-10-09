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
         " shared time-varying (trend) field(s); pass `time.col = ",
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

# Row pooling for predict(weights = ). `newdata` may carry several rows per
# prediction cell, one per level of a covariate the cell mixes (a habitat
# composition, a set of plot sizes); `weights` gives each row's share of its cell.
# Resolves the weights (a numeric vector or the name of a `newdata` column),
# normalises them within each cell and returns the pooling map: `w` the
# normalised per-row weight, `g` the row's output-row index, `cell` the output
# cells in order of first appearance. NULL weights -> NULL (no pooling).
.tobs_predict_pool <- function(weights, newdata, cell) {
  if (is.null(weights)) return(NULL)
  if (is.character(weights)) {
    if (length(weights) != 1L || !weights %in% names(newdata)) {
      stop("predict(): `weights = \"", paste(weights, collapse = "\", \""),
           "\"` must name one column of `newdata`.", call. = FALSE)
    }
    weights <- newdata[[weights]]
  }
  w <- as.numeric(weights)
  if (length(w) != nrow(newdata)) {
    stop("predict(): `weights` has length ", length(w), " but `newdata` has ",
         nrow(newdata), " rows (one weight per row).", call. = FALSE)
  }
  if (anyNA(w) || !all(is.finite(w)) || any(w < 0)) {
    stop("predict(): `weights` must be finite and non-negative.", call. = FALSE)
  }
  key <- unique(cell)
  g   <- match(cell, key)
  tot <- as.vector(rowsum(w, g, reorder = FALSE))
  if (any(tot <= 0)) {
    stop("predict(): the `weights` of cell(s) ",
         paste(utils::head(key[tot <= 0], 5L), collapse = ", "),
         " sum to zero; every cell needs a positive total weight.", call. = FALSE)
  }
  list(w = w / tot[g], g = g, cell = key)
}

# Pool per-draw arm states over each cell's rows (`pool` from
# .tobs_predict_pool()). Occupancy, detection and expected cover are the
# weighted means of their rows; conditional cover is the pooled expected cover
# over the pooled occupancy, i.e. E[cover | present] under the cell's mixture,
# which reduces to the weighted mean of the rows' conditional cover when
# occupancy does not vary across them. Pooling happens per draw, before any
# summary or difference, so the pooled table and its change share one draw set.
.tobs_pool_states <- function(st, pool) {
  if (is.null(pool)) return(st)
  agg <- function(m) {
    out <- rowsum(pool$w * m, pool$g, reorder = FALSE)
    dimnames(out) <- NULL
    out
  }
  out <- list()
  if (!is.null(st$p))     out$p     <- agg(st$p)
  if (!is.null(st$p_det)) out$p_det <- agg(st$p_det)
  if (!is.null(st$mu)) {
    if (!is.null(st$E) && !is.null(out$p)) {
      out$E  <- agg(st$E)
      mu_w   <- agg(st$mu)
      pos    <- out$p > 0
      out$mu <- mu_w
      out$mu[pos] <- out$E[pos] / out$p[pos]
    } else {
      out$mu <- agg(st$mu)
    }
  }
  out
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

# The prediction frames of a change map or trajectory: `newdata` with the time
# covariate held at each of `times`. For a change map, two times give one
# difference and K times a sequence of steps, each differenced against
# `times[1]`; a trajectory takes one or more times.
.tobs_joint_change_frames <- function(object, newdata, times, time_col,
                                      type = "change") {
  if (identical(type, "trajectory")) {
    if (is.null(times) || !length(times)) {
      stop("predict(type = \"trajectory\") needs `times`: the values of the ",
           "time covariate to predict at.", call. = FALSE)
    }
  } else if (is.null(times) || length(times) < 2L) {
    stop("predict(type = \"change\") needs `times = c(t1, t2)` (one change) or ",
         "`times = c(t1, ..., tK)` (a trajectory, each step against t1).",
         call. = FALSE)
  }
  if (!is.numeric(times) || anyNA(times) || !all(is.finite(times))) {
    stop("predict(type = \"", type, "\"): `times` must be finite numeric values.",
         call. = FALSE)
  }
  if (is.null(time_col)) time_col <- object$trend_weight
  if (is.null(time_col)) {
    stop("predict(type = \"", type, "\") needs the name of the time covariate. ",
         "Pass `time.col = \"<column>\"` (the covariate the `times` are held ",
         "at, whose movement drives the prediction).", call. = FALSE)
  }
  if (!time_col %in% names(newdata)) {
    stop("predict(): `time.col = \"", time_col, "\"` is not a column of ",
         "`newdata`.", call. = FALSE)
  }
  lapply(times, function(tk) { nd <- newdata; nd[[time_col]] <- tk; nd })
}

# Monte Carlo precision of draw-based quantiles. A sample quantile at
# probability p from n draws has standard error sqrt(p (1 - p) / n) s(p), with
# s(p) = dQ/dp the sparsity (inverse density at the quantile). s(p) is the
# difference quotient of the sample quantiles at p - h and p + h, h the
# Hall-Sheather (1988) bandwidth, the standard choice for quantile standard
# errors (quantreg's `hs = TRUE`); a bracket clipped at 0 or 1 divides by the
# span it kept. Dividing by the row's posterior SD states the error in units of
# posterior spread, comparable across quantities on different scales. The SD is
# floored at `floor` on the quantity's own scale: where a probability saturates
# near 0 or 1 its SD goes to ~0 and an error relative to it is unstable. `m` is
# [rows x draws]; returns the [length(probs) x rows] errors. A row with no
# spread and no floor has no MC error (0).
.tobs_hall_sheather <- function(n, p, alpha = 0.05) {
  zp <- stats::qnorm(p)
  n^(-1 / 3) * stats::qnorm(1 - alpha / 2)^(2 / 3) *
    (1.5 * stats::dnorm(zp)^2 / (2 * zp^2 + 1))^(1 / 3)
}
.tobs_mc_quantile_error <- function(m, probs, floor = 0) {
  n  <- ncol(m)
  h  <- .tobs_hall_sheather(n, probs)
  lo <- pmax(probs - h, 0)
  hi <- pmin(probs + h, 1)
  q  <- matrix(apply(m, 1L, stats::quantile, probs = c(lo, hi), names = FALSE),
               nrow = 2L * length(probs))
  np <- length(probs)
  se <- (q[np + seq_len(np), , drop = FALSE] - q[seq_len(np), , drop = FALSE]) *
    (sqrt(probs * (1 - probs) / n) / (hi - lo))
  mu <- rowMeans(m)
  scale <- pmax(sqrt(rowSums((m - mu)^2) / (n - 1L)), floor)
  err <- sweep(se, 2L, scale, "/")
  err[, !(scale > 0)] <- 0
  err
}

# Long trajectory table over the K prediction frames `nds` (one per time). For
# every quantity `state_at(nd)` returns -- a named list of [n x nsim] draw
# matrices -- it reports `<q>_mean`, `<q>_median` and the `level` interval
# `<q>_lwr` / `<q>_upr` at each time; one row per unit x time, units fastest.
# Every time is evaluated on the same draw set, so the rows are jointly
# consistent across units, quantities and times. `aggregate = TRUE` first
# averages each quantity over the units per draw, giving one row per time. A
# time's draws are summarised and released before the next is evaluated, so
# memory stays flat in K unless `draws = TRUE` keeps them. `attr(, "mc_se_max")`
# is the largest MC standard error of any reported median or bound, over every
# row, quantity and time, in units of that row's posterior SD (floored at
# `mc.floor`); `attr(, ".mc_se")` carries every one of those errors for
# .tobs_mc_trajectory() to pool, bound fastest, then row, quantity and time, and
# `attr(, ".mc_layout")` the quantity names and rows per time that index them.
.tobs_trajectory_table <- function(state_at, nds, times, cell, level,
                                   aggregate = FALSE, draws = FALSE,
                                   mc.floor = 0) {
  a <- (1 - level) / 2
  probs <- c(0.5, a, 1 - a)
  mc_se <- list()
  summ <- function(m, q) {
    qs <- matrix(apply(m, 1L, stats::quantile, probs = probs, names = FALSE),
                 nrow = 3L)
    mc_se[[length(mc_se) + 1L]] <<-
      as.vector(.tobs_mc_quantile_error(m, probs, mc.floor))
    stats::setNames(list(rowMeans(m), qs[1L, ], qs[2L, ], qs[3L, ]),
                    paste0(q, c("_mean", "_median", "_lwr", "_upr")))
  }
  rows <- vector("list", length(nds))
  dr   <- list()
  n_draws <- NA_integer_
  for (k in seq_along(nds)) {
    st <- state_at(nds[[k]])
    n_draws <- ncol(st[[1L]])
    if (aggregate) st <- lapply(st, function(m) matrix(colMeans(m), nrow = 1L))
    if (isTRUE(draws))
      for (q in names(st)) dr[[sprintf("%s_T%d", q, k)]] <- st[[q]]
    key <- if (aggregate) list(time = times[k])
           else list(cell = cell, time = rep(times[k], nrow(st[[1L]])))
    rows[[k]] <- as.data.frame(
      c(key, do.call(c, lapply(names(st), function(q) summ(st[[q]], q)))),
      check.names = FALSE)
  }
  tbl <- do.call(rbind, rows)
  rownames(tbl) <- NULL
  attr(tbl, "quantity")  <- "trajectory"
  attr(tbl, "times")     <- times
  attr(tbl, "aggregate") <- aggregate
  mc_se <- unlist(mc_se, use.names = FALSE)
  attr(tbl, "mc_se_max") <- max(mc_se)
  attr(tbl, ".mc_se")    <- mc_se
  attr(tbl, ".mc_layout") <- list(quantities = names(st), n_rows = nrow(st[[1L]]),
                                  bounds = c("median", "lwr", "upr"))
  attr(tbl, "nsim_used") <- n_draws
  if (isTRUE(draws)) attr(tbl, "draws") <- dr
  class(tbl) <- c("tobs_prediction", "data.frame")
  tbl
}

# `nsim`: a positive whole number of draws, or "auto" (trajectory only) to draw
# until every reported quantile meets `mc.tol`. Returns the integer count or
# "auto".
.tobs_check_nsim <- function(nsim, type, mc.tol, nsim.max, mc.floor = 0) {
  if (!is.numeric(mc.floor) || length(mc.floor) != 1L || !is.finite(mc.floor) ||
      mc.floor < 0) {
    stop("predict(): `mc.floor` must be one non-negative number (the smallest ",
         "posterior SD the Monte Carlo error is measured against, on the ",
         "quantity's own scale).", call. = FALSE)
  }
  if (identical(nsim, "auto")) {
    if (!identical(type, "trajectory")) {
      stop("predict(nsim = \"auto\") applies to type = \"trajectory\" only; ",
           "pass a whole number of draws.", call. = FALSE)
    }
    if (!is.numeric(mc.tol) || length(mc.tol) != 1L || !is.finite(mc.tol) ||
        mc.tol <= 0) {
      stop("predict(): `mc.tol` must be one positive number (the Monte Carlo ",
           "error allowed on each reported bound, as a fraction of its ",
           "posterior SD).", call. = FALSE)
    }
    if (!is.numeric(nsim.max) || length(nsim.max) != 1L ||
        !is.finite(nsim.max) || nsim.max < 2 || nsim.max != round(nsim.max)) {
      stop("predict(): `nsim.max` must be a whole number of at least 2.",
           call. = FALSE)
    }
    return("auto")
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || !is.finite(nsim) ||
      nsim < 2 || nsim != round(nsim)) {
    stop("predict(): `nsim` must be a whole number of draws (at least 2) or ",
         "\"auto\".", call. = FALSE)
  }
  as.integer(nsim)
}

# Draw count of the first "auto" round, the step a round grows by at least, and
# the factor it grows by at most.
.TOBS_MC_BATCH  <- 250L
.TOBS_MC_GROWTH <- 4L

# Trajectory table with the draw count chosen by Monte Carlo precision.
# `state_of(bundle)` returns the per-frame state function over that draw
# bundle; `redraw(n)` returns a fresh bundle of `n` joint posterior draws.
# With an integer `nsim` the table is built once on `bundle`. With "auto", the
# table is rebuilt on fresh draws until the worst MC error meets `mc.tol`.
# A quantile's MC error times sqrt(n) is a property of the posterior alone, so
# every round is an independent estimate of the same per-bound constant c; the
# rounds' estimates are pooled with inverse-variance weights n^(2/3) (the
# Hall-Sheather estimate's CV falls as n^(-1/3)). Pooling keeps the noise in
# one round's largest estimate from setting the draw count. While
# max(c) / sqrt(n) > mc.tol the next round draws the count that projects to the
# target, at least one batch more, at most `.TOBS_MC_GROWTH` times the current
# count (a small round's largest estimate is noisy upward, and projecting from
# it alone overshoots) and at most `nsim.max`. Every round evaluates
# one draw set end to end, so peak memory is that of a single run at the final
# count and all quantities keep sharing one draw set. `mc_se_max` reports the
# pooled worst error at the final count and `mc_binding` the cell, time,
# quantity and bound it belongs to, so a run the cap stopped shows where.
.tobs_mc_trajectory <- function(state_of, bundle, redraw, nds, times, cell,
                                level, aggregate, draws, nsim, mc.tol,
                                nsim.max, mc.floor = 0) {
  auto <- identical(nsim, "auto")
  c_sum <- 0
  w_sum <- 0
  repeat {
    tbl <- .tobs_trajectory_table(state_of(bundle), nds, times, cell, level,
                                  aggregate, draws, mc.floor)
    n     <- attr(tbl, "nsim_used")
    w     <- n^(2 / 3)
    c_sum <- c_sum + w * attr(tbl, ".mc_se") * sqrt(n)
    w_sum <- w_sum + w
    c_bar <- c_sum / w_sum
    c_max <- max(c_bar)
    err   <- c_max / sqrt(n)
    if (!auto || err <= mc.tol || n >= nsim.max) break
    proj <- ceiling(1.1 * (c_max / mc.tol)^2 / .TOBS_MC_BATCH) * .TOBS_MC_BATCH
    bundle <- redraw(as.integer(min(nsim.max, .TOBS_MC_GROWTH * n,
                                    max(n + .TOBS_MC_BATCH, proj))))
  }
  attr(tbl, "mc_binding")  <- .tobs_mc_binding(which.max(c_bar), err,
                                               attr(tbl, ".mc_layout"), times,
                                               cell, aggregate)
  attr(tbl, ".mc_se")      <- NULL
  attr(tbl, ".mc_layout")  <- NULL
  attr(tbl, "mc_se_max")   <- err
  attr(tbl, "mc_tol")      <- if (auto) mc.tol else NA_real_
  attr(tbl, "nsim_capped") <- auto && err > mc.tol
  tbl
}


# The cell, time, quantity and bound of error element `i` of a trajectory
# table's `.mc_se` (bound fastest, then row, quantity and time), as a one-row
# data frame carrying its error `mc_se`. `cell` is NA on an aggregated table.
.tobs_mc_binding <- function(i, mc_se, layout, times, cell, aggregate) {
  nb <- length(layout$bounds)
  nr <- layout$n_rows
  nq <- length(layout$quantities)
  i0 <- i - 1L
  data.frame(
    cell     = if (aggregate) NA_integer_ else cell[(i0 %/% nb) %% nr + 1L],
    time     = times[i0 %/% (nb * nr * nq) + 1L],
    quantity = layout$quantities[(i0 %/% (nb * nr)) %% nq + 1L],
    bound    = layout$bounds[i0 %% nb + 1L],
    mc_se    = mc_se,
    stringsAsFactors = FALSE)
}
# `aggregate` is a single TRUE / FALSE and applies to a trajectory only.
.tobs_check_aggregate <- function(aggregate, type) {
  if (!is.logical(aggregate) || length(aggregate) != 1L || is.na(aggregate)) {
    stop("predict(): `aggregate` must be TRUE or FALSE.", call. = FALSE)
  }
  if (aggregate && !identical(type, "trajectory")) {
    stop("predict(aggregate = TRUE) applies to type = \"trajectory\" only.",
         call. = FALSE)
  }
  aggregate
}

# Per-row summaries over a draw matrix: the `level` central-interval endpoints
# and the posterior SD.
.tobs_draw_lwr <- function(m, level)
  apply(m, 1L, stats::quantile, probs = (1 - level) / 2, names = FALSE)
.tobs_draw_upr <- function(m, level)
  apply(m, 1L, stats::quantile, probs = 1 - (1 - level) / 2, names = FALSE)
.tobs_draw_sd <- function(m) apply(m, 1L, stats::sd)

# Coefficient-level predict() for a non-joint occu_cover() fit (laplace / nuts,
# `.tobs_joint_fit(object)` is NULL). Reuses `.tobs_occu_cover_components()`,
# the per-draw arm coefficients / field / random-effect offsets fitted() and
# WAIC / LOO already share (`R/occu_cover_diag.R`), and builds the new design
# at `newdata` with `.tobs_joint_arm_design()` -- which reads the fit's own
# formulas rather than the joint substrate, so it works unchanged here. A
# shared spatial field or a detection/cover-arm random effect is tied to the
# fitted cell graph / grouping, so newdata prediction (kriging to an unseen
# cell or an unseen level) is refused for those; the in-sample surface (no
# `newdata`) is unaffected -- predict.tobs_fit() routes that to fitted()
# instead of here, which already folds in the field and every random effect
# (#351).
.tobs_predict_occu_cover_coef <- function(object, newdata, type, level, nsim,
                                          weights = NULL) {
  type <- match.arg(type, c("occurrence", "detection", "cover_cond",
                            "cover_exp", "change", "trajectory"))
  if (identical(nsim, "auto")) {
    stop("predict(nsim = \"auto\") applies to type = \"trajectory\" on a joint ",
         "nested-Laplace fit (method = \"nested_laplace\"); pass a whole ",
         "number of draws.", call. = FALSE)
  }
  if (type %in% c("change", "trajectory")) {
    stop("predict(type = \"", type, "\") needs a joint nested-Laplace fit ",
         "(method = \"nested_laplace\"); this fit carries no joint object.",
         call. = FALSE)
  }
  cmp <- .tobs_occu_cover_components(object, n.draws = nsim)
  if (any(cmp$field_occ != 0) || any(cmp$field_pos != 0)) {
    stop("predict(newdata = ) is not supported for this occu_cover() fit: it ",
         "carries a shared spatial field, which is tied to the fitted cell ",
         "graph and cannot be evaluated at a new location. Call predict() ",
         "with no `newdata` for the in-sample fit, or refit with method = ",
         "\"nested_laplace\" for a field-aware newdata / change-map ",
         "predictor.", call. = FALSE)
  }
  if ((!is.null(cmp$off_det) && ncol(cmp$off_det) > 0L) ||
      (!is.null(cmp$off_pos) && ncol(cmp$off_pos) > 0L)) {
    stop("predict(newdata = ) is not supported for this occu_cover() fit: it ",
         "carries a detection / cover-arm random effect, which is tied to the ",
         "fitted grouping. Call predict() with no `newdata` for the in-sample ",
         "fit.", call. = FALSE)
  }
  p_occ <- ncol(cmp$b_occ); p_det <- ncol(cmp$b_det); p_pos <- ncol(cmp$b_pos)
  X_occ   <- .tobs_joint_arm_design(object, newdata, "occ", p_occ)
  eta_occ <- tcrossprod(X_occ, cmp$b_occ)
  st <- list(p = stats::plogis(eta_occ))
  if (identical(type, "detection")) {
    X_det    <- .tobs_joint_arm_design(object, newdata, "det", p_det)
    st$p_det <- stats::plogis(tcrossprod(X_det, cmp$b_det))
  }
  if (type %in% c("cover_cond", "cover_exp")) {
    X_pos   <- .tobs_joint_arm_design(object, newdata, "pos", p_pos)
    eta_pos <- tcrossprod(X_pos, cmp$b_pos)
    st$mu   <- .occu_cover_mu_from_eta(eta_pos, cmp$disp, object$model$positive)
    st$E    <- st$p * st$mu
  }
  unit <- if (!is.null(newdata$cell)) as.integer(newdata$cell)
          else seq_len(nrow(newdata))
  pool <- .tobs_predict_pool(weights, newdata, unit)
  st   <- .tobs_pool_states(st, pool)

  mat <- switch(type,
                occurrence = st$p,
                detection  = st$p_det,
                cover_cond = st$mu,
                cover_exp  = st$E)

  tbl <- .occu_cover_summ(mat, if (is.null(pool)) seq_len(nrow(newdata))
                               else pool$cell, level)
  attr(tbl, "quantity") <- type
  attr(tbl, "draws") <- stats::setNames(list(mat), type)
  class(tbl) <- c("tobs_prediction", "data.frame")
  tbl
}

# Core predict handler for the joint cover-family fits. `object` is an
# occu_cover() fit (3-arm) or a cover() hurdle fit on the nested-Laplace path
# (2-arm); both expose a joint nested-Laplace object via `.tobs_joint_fit()`.
.tobs_predict_joint <- function(object, newdata = NULL,
                                type = "occurrence", times = NULL,
                                level = 0.95, nsim = 1000L,
                                draws = TRUE, time_col = NULL,
                                weights = NULL, aggregate = FALSE,
                                mc.tol = 0.05, nsim.max = 10000L,
                                mc.floor = 0.001) {
  type <- match.arg(type, c("occurrence", "detection", "cover_cond",
                            "cover_exp", "change", "trajectory"))
  aggregate <- .tobs_check_aggregate(aggregate, type)
  nsim <- .tobs_check_nsim(nsim, type, mc.tol, nsim.max, mc.floor)
  if (is.null(.tobs_joint_fit(object))) {
    stop("predict() requires a joint nested-Laplace fit (method = ",
         "\"nested_laplace\"); this fit carries no joint object.",
         call. = FALSE)
  }

  redraw <- function(n) .tobs_joint_trend_blocks(
    object, .tobs_joint_draws(object, n = n), time_col)
  bundle <- redraw(if (identical(nsim, "auto"))
                     min(.TOBS_MC_BATCH, as.integer(nsim.max)) else nsim)
  nc      <- .tobs_joint_predict_cells(object, newdata, bundle$n_cells)
  newdata <- nc$newdata
  cell    <- nc$cell
  pool    <- .tobs_predict_pool(weights, newdata, cell)
  out_cell <- if (is.null(pool)) cell else pool$cell

  state <- function(nd, b = bundle)
    .tobs_pool_states(.tobs_joint_arm_states(object, b, nd, cell), pool)

  # --- trajectory: every quantity at each of `times` -----------------------
  if (identical(type, "trajectory")) {
    nds <- .tobs_joint_change_frames(object, newdata, times, time_col, type)
    return(.tobs_mc_trajectory(function(b) function(nd) {
      s <- state(nd, b)
      list(psi = s$p, cover_cond = s$mu, cover_exp = s$E)
    }, bundle, redraw, nds, times, out_cell, level, aggregate, draws, nsim,
    mc.tol, nsim.max, mc.floor))
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
    tbl <- .occu_cover_summ(mat, out_cell, level)
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
  s1  <- st[[1L]]
  steps <- lapply(st, function(sk)
    list(psi = sk$p, cover_cond = sk$mu, cover_exp = sk$E))
  # Exact additive split of the expected-cover change against the first step:
  #   delta_exp = pk muk - p1 mu1 = (pk - p1) mu1 + pk (muk - mu1).
  parts <- lapply(st[-1L], function(sk) list(
    cover_from_occ = (sk$p - s1$p) * s1$mu,
    cover_from_ab  = sk$p * (sk$mu - s1$mu)))
  .tobs_change_table(steps, out_cell, times, level, draws, parts)
}

# Per-cell change table over K time steps, shared by the joint predict handlers.
# `steps` is a list over the K times, each a named list of [n x nsim] draw
# matrices carrying the same quantities at every step. Every level `<q>_T<k>`
# reports its posterior mean with `.sd` / `.lwr` / `.upr`; every delta against
# the first step `delta_<q>` reports its mean with `.lwr` / `.upr` and the
# posterior probability `.prob_pos` = P(delta > 0). All summaries are taken over
# one draw set, so levels and deltas are jointly consistent. `parts` optionally
# adds further per-step delta draw matrices (a list over steps 2..K). With two
# times the delta columns carry no step suffix; with K > 2 they carry `_T<k>`.
.tobs_change_table <- function(steps, cell, times, level, draws, parts = NULL) {
  K   <- length(steps)
  sfx <- if (K == 2L) rep("", K) else paste0("_T", seq_len(K))
  summ_level <- function(nm, m) stats::setNames(
    list(rowMeans(m), .tobs_draw_sd(m), .tobs_draw_lwr(m, level),
         .tobs_draw_upr(m, level)),
    paste0(nm, c("", ".sd", ".lwr", ".upr")))
  summ_delta <- function(nm, m) stats::setNames(
    list(rowMeans(m), .tobs_draw_lwr(m, level), .tobs_draw_upr(m, level),
         rowMeans(m > 0)),
    paste0(nm, c("", ".lwr", ".upr", ".prob_pos")))

  cols <- list(cell = cell)
  dr   <- list()
  for (q in names(steps[[1L]])) {
    for (k in seq_len(K)) {
      nm <- sprintf("%s_T%d", q, k)
      m  <- steps[[k]][[q]]
      cols <- c(cols, summ_level(nm, m)); dr[[nm]] <- m
    }
    for (k in seq_len(K)[-1L]) {
      nm <- paste0("delta_", q, sfx[k])
      m  <- steps[[k]][[q]] - steps[[1L]][[q]]
      cols <- c(cols, summ_delta(nm, m)); dr[[nm]] <- m
    }
  }
  for (q in names(parts[[1L]])) {
    for (k in seq_len(K)[-1L]) {
      nm <- paste0("delta_", q, sfx[k])
      m  <- parts[[k - 1L]][[q]]
      cols <- c(cols, summ_delta(nm, m)); dr[[nm]] <- m
    }
  }

  tbl <- as.data.frame(cols, check.names = FALSE)
  attr(tbl, "quantity") <- "change"
  attr(tbl, "times")    <- times
  if (isTRUE(draws)) attr(tbl, "draws") <- dr
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
                                     draws = TRUE, time_col = NULL,
                                     weights = NULL, aggregate = FALSE,
                                     mc.tol = 0.05, nsim.max = 10000L,
                                     mc.floor = 0.001) {
  type <- match.arg(type, c("occupancy", "detection", "both", "change",
                            "trajectory"))
  aggregate <- .tobs_check_aggregate(aggregate, type)
  nsim <- .tobs_check_nsim(nsim, type, mc.tol, nsim.max, mc.floor)
  if (is.null(.tobs_joint_fit(object))) {
    stop("predict() requires a joint nested-Laplace fit; this occu() fit ",
         "carries no joint object.", call. = FALSE)
  }

  redraw <- function(n) .tobs_joint_trend_blocks(
    object, .tobs_joint_draws(object, n = n), time_col)
  bundle <- redraw(if (identical(nsim, "auto"))
                     min(.TOBS_MC_BATCH, as.integer(nsim.max)) else nsim)
  nc      <- .tobs_joint_predict_cells(object, newdata, bundle$n_cells)
  newdata <- nc$newdata
  cell    <- nc$cell
  pool    <- .tobs_predict_pool(weights, newdata, cell)
  if (!is.null(pool)) cell <- pool$cell

  occ_state <- function(nd, b = bundle) .tobs_pool_states(
    .tobs_joint_arm_states(object, b, nd, nc$cell, "occ"), pool)$p
  det_state <- function(nd) .tobs_pool_states(
    .tobs_joint_arm_states(object, bundle, nd, nc$cell, "det"), pool)$p_det

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

  if (identical(type, "trajectory")) {
    nds <- .tobs_joint_change_frames(object, newdata, times, time_col, type)
    return(.tobs_mc_trajectory(
      function(b) function(nd) list(psi = occ_state(nd, b)), bundle, redraw,
      nds, times, cell, level, aggregate, draws, nsim, mc.tol, nsim.max,
      mc.floor))
  }

  # type == "change": occupancy at each of `times`, differenced against the first.
  nds <- .tobs_joint_change_frames(object, newdata, times, time_col)
  steps <- lapply(nds, function(nd) list(psi = occ_state(nd)))
  .tobs_change_table(steps, cell, times, level, draws)
}
