# predict() for the joint occu_cover fit (gcol33/tulpaObs#22).
#
# Returns a posterior-draws object with a tidy accessor -- the tulpaObs analogue
# of the hand-rolled INLA pred.inla.joint() collaborators write
# (inla.posterior.sample + control.predictor, then difference the draws). The
# deliverable is "plotting objects, not plots": hand back the per-unit draw
# matrices plus a plot-ready table; the user maps it.
#
# Draws source: the joint latent (beta_psi, beta_pos, field) is sampled from the
# grid-integrated posterior via tulpa::tulpa_posterior_draws() -- the faithful
# outer-grid mixture sum_k w_k N(m_k, V_k), with the field amplitude restored per
# draw from the grid cell the draw came from (attr "cells"). Every derived
# quantity (p, mu, expected cover, and the change deltas) is computed PER DRAW
# and only then summarized -- the "Marginalize Derived Quantities" rule -- so the
# nonlinear products and differences carry the joint posterior's correlation
# rather than a plug-in of posterior means.

# Reconstruct the per-arm linear predictors at `newdata` for every posterior
# draw, returning per-unit x nsim matrices of p (occupancy), mu (conditional
# cover) and E (expected cover). Shared across the single-time and the change
# (two-time) paths so the same posterior draws back both states.
#
# `blocks` is a list of shared-field blocks, each:
#   $z       : nsim x n_cells latent draws (unit-scale ICAR field)
#   $sigma   : length-nsim per-draw field amplitude (the donor-arm scaling)
#   $alpha   : length-nsim per-draw copy coefficient (the cover-arm scaling)
#   $weight  : NULL for the intercept field (weight 1), else the name of the
#              per-cell covariate weighting this (trend / SVC) field. The field's
#              contribution to a unit is weight_value * amplitude * z[cell].
# This loop is the single source of truth for both the single-shared-field fit
# (one block, no weight) and the trend fit (intercept block + N weighted trend
# blocks); nothing is special-cased on block count.
.occu_cover_state_draws <- function(object, newdata, b_psi, b_pos, blocks,
                                    cell) {
    model    <- object$model
    positive <- object$positive %||% "lognormal"

    # Occupancy is purely cell-level: the occ design is the full psi design.
    X_occ0 <- .occu_cover_arm_design(model, newdata, arm = "occ",
                                     p_arm = ncol(b_psi))
    # The cover arm carries cell-level (site) coefficients plus, when the fit
    # used visit covariates, trailing visit-level ones. A per-cell prediction
    # has no visit, so the visit columns are held at the reference (0) -- the
    # prediction is at the average visit condition.
    X_pos0 <- .occu_cover_arm_design(model, newdata, arm = "pos",
                                     p_arm = ncol(b_pos))

    # Fixed-effect part. tcrossprod(X, B) is [n_unit x nsim].
    eta_occ <- tcrossprod(X_occ0, b_psi)
    eta_pos <- tcrossprod(X_pos0, b_pos)

    # Shared-field part. For each block, the per-unit contribution is
    #   w[unit] * amplitude[draw] * z[draw, cell[unit]].
    # `z[, cell] * sigma` scales each draw-row by its amplitude (sigma has
    # length nsim == nrow, so the column-wise recycling hits row d with
    # sigma[d]); the transpose to [n_unit x nsim] is then row-scaled by the
    # per-unit weight w (length n_unit == nrow after transpose).
    for (blk in blocks) {
        z_cell <- blk$z[, cell, drop = FALSE]            # nsim x n_unit
        w <- if (is.null(blk$weight)) rep(1, length(cell)) else {
            if (!blk$weight %in% names(newdata)) {
                stop("predict(occu_cover): trend-field weight column '",
                     blk$weight, "' is not in `newdata`.", call. = FALSE)
            }
            as.numeric(newdata[[blk$weight]])
        }
        eta_occ <- eta_occ + t(z_cell * blk$sigma) * w
        eta_pos <- eta_pos + t(z_cell * (blk$alpha * blk$sigma)) * w
    }

    p <- stats::plogis(eta_occ)
    if (identical(positive, "beta")) {
        mu <- stats::plogis(eta_pos)
    } else {
        sigma_pos <- object$sigma_pos %||% 1
        mu <- exp(eta_pos + sigma_pos^2 / 2)
    }
    list(p = p, mu = mu, E = p * mu)
}

# Build the per-cell design matrix for one arm at `newdata`, matching the fitted
# coefficient count `p_arm`. Occupancy ("occ") is cell-level and uses its whole
# design. The cover arm ("pos") splits into a cell-level (site) block built from
# the stored site formula and, when the fit carried visit covariates, trailing
# visit-level columns; per-cell prediction holds those at the reference (0).
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
    # Pad visit-level columns (held at the reference) so the design matches the
    # fitted coefficient vector length.
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

# Core predict handler for occu_cover joint fits.
.tobs_predict_occu_cover <- function(object, newdata = NULL,
                                     type = "occurrence", times = NULL,
                                     level = 0.95, nsim = 1000L,
                                     draws = TRUE, time_col = NULL) {
    type <- match.arg(type, c("occurrence", "cover_cond", "cover_exp", "change"))
    jf <- object$joint_fit
    if (is.null(jf)) {
        stop("predict(occu_cover) requires a joint spatial fit (method = ",
             "\"joint_coupled\"); this fit carries no `joint_fit`.",
             call. = FALSE)
    }
    model   <- object$model
    layout  <- jf$arm_layout
    n_cells <- model$n_sites

    # Shared field blocks: block 1 is the time-invariant intercept field;
    # blocks 2.. are time-varying (trend / SVC) fields. `field_starts` (multi-
    # block fit) or `phi_start` (single-block fit) carries the 0-based offset
    # of each field; every field is n_cells long.
    starts  <- layout$field_starts %||% layout$phi_start
    n_field <- length(starts)

    # For a trend fit the change map (and any single-time prediction) needs the
    # per-cell time covariate that weights each trend field. Require its column
    # name so newdata can be evaluated at the right time(s).
    if (n_field > 1L) {
        if (is.null(time_col)) time_col <- object$trend_weight
        if (is.null(time_col)) {
            stop("predict(occu_cover): this fit has ", n_field - 1L,
                 " time-varying (trend) field(s); pass `time_col = \"<column>\"`, ",
                 "the per-cell covariate that weights the trend field(s) (the ",
                 "same column used at fit time via control$trend).",
                 call. = FALSE)
        }
    }

    if (is.null(newdata)) newdata <- model$data

    # cell map: explicit `cell` column, else row i -> field cell i.
    if (!is.null(newdata$cell)) {
        cell <- as.integer(newdata$cell)
    } else {
        cell <- seq_len(nrow(newdata))
    }
    if (anyNA(cell) || any(cell < 1L) || any(cell > n_cells)) {
        stop("predict(occu_cover): `cell` must index field cells 1..", n_cells,
             " (add a `cell` column to `newdata`, or pass one row per field cell ",
             "in cell order).", call. = FALSE)
    }

    # --- joint latent draws: (beta_psi, beta_pos, every field block) ---------
    p_psi <- layout$p[1L]
    p_pos <- layout$p[3L]
    idx_psi <- layout$beta_start[1L] + seq_len(p_psi)
    idx_pos <- layout$beta_start[3L] + seq_len(p_pos)
    field_idx_list <- lapply(starts, function(s0) s0 + seq_len(n_cells))
    idx <- c(idx_psi, idx_pos, unlist(field_idx_list))

    D <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = nsim)
    cells_drawn <- attr(D, "cells")
    b_psi <- D[, seq_len(p_psi), drop = FALSE]
    b_pos <- D[, p_psi + seq_len(p_pos), drop = FALSE]

    # Per-block latent draws + per-draw (sigma, alpha) amplitudes. Amplitude
    # axes are bare ("sigma"/"alpha") on a single-block fit and block-prefixed
    # ("b<b>.sigma"/"b<b>.alpha") on a multi-block fit; resolve by trying the
    # prefixed name first, then the bare name.
    tg <- jf$theta_grid
    cn <- colnames(tg)
    amp_axis <- function(b, name) {
        j <- match(paste0("b", b, ".", name), cn)
        if (is.na(j)) j <- match(name, cn)
        if (is.na(j)) return(rep(if (name == "alpha") 1 else 1, nsim))
        as.numeric(tg[cells_drawn, j])
    }
    base <- p_psi + p_pos
    blocks <- lapply(seq_len(n_field), function(b) {
        zb <- D[, base + (b - 1L) * n_cells + seq_len(n_cells), drop = FALSE]
        list(z      = zb,
             sigma  = amp_axis(b, "sigma"),
             alpha  = amp_axis(b, "alpha"),
             weight = if (b == 1L) NULL else time_col)
    })

    state <- function(nd) {
        .occu_cover_state_draws(object, nd, b_psi, b_pos, blocks, cell)
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
        stop("predict(occu_cover, type = \"change\") needs `times = c(t1, t2)`.",
             call. = FALSE)
    }
    if (is.null(time_col)) time_col <- object$trend_weight
    if (is.null(time_col)) {
        stop("predict(occu_cover, type = \"change\") needs the name of the time ",
             "covariate. Pass `time_col = \"<column>\"` (the covariate whose ",
             "change between `times[1]` and `times[2]` drives the prediction).",
             call. = FALSE)
    }
    if (!time_col %in% names(newdata)) {
        stop("predict(occu_cover): `time_col = \"", time_col, "\"` is not a ",
             "column of `newdata`.", call. = FALSE)
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
