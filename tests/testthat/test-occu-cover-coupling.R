# FD-checks for every closed-form derivative of the OccuCoverCoupling
# CellCouplingSpec (tulpaObs consumer of gcol33/tulpa#32 Layer B.2), across its
# three positive-arm policies -- lognormal, beta and identity-Gaussian -- and
# both cover granularities (per-visit and cell-aggregated).
#
# The spec lives in src/cell_coupling_occu_cover.h; the direct evaluators
# `cpp_eval_occu_cover_{lognormal,beta,gaussian}[_agg]_cell` return the cell
# log-density at given etas plus every nonzero derivative buffer the kernel
# reads (gradients, diagonal neg-Hess, cross-Hess (psi, p_v), (p_v, p_w)). Each
# is compared against a central finite difference of cell_ll under an explicit
# tolerance, so a failure reports the size of the disagreement.
#
# Curvature modes. The det branch is complete-data (z is known once a visit
# detects), so Observed and Expected coincide there and both have to clear the
# FD checks on every arm, cover arm included. The nodet branch's Expected
# curvature is the complete-data Fisher information, which is not the second
# derivative of the mixture, so it is checked against its closed form instead.

fd1 <- function(f, h = 1e-5) (f(+h) - f(-h)) / (2 * h)
fd2 <- function(f, h = 1e-4) (f(+h) - 2 * f(0) + f(-h)) / (h * h)
fd_cross <- function(g, h = 1e-3) {
    (g(+h, +h) - g(+h, -h) - g(-h, +h) + g(-h, -h)) / (4 * h * h)
}

bump <- function(vec, i, h) { vec[i] <- vec[i] + h; vec }

# Per-family FD settings. `tol1` bounds the gradient checks and `tol2` the
# curvature checks: the beta arm carries a digamma / trigamma evaluation and
# holds a decade looser than the two arms whose curvature is the closed-form
# 1 / sigma^2. `seed` roots that family's synthetic cells.
oc_fam <- list(
    lognormal = list(tol1 = 1e-6, tol2 = 1e-4, seed = 101L),
    beta      = list(tol1 = 1e-5, tol2 = 1e-3, seed = 707L),
    gaussian  = list(tol1 = 1e-6, tol2 = 1e-4, seed = 313L)
)
oc_fams <- names(oc_fam)

# Family-agnostic evaluator. `disp` is the pos arm's dispersion: the SD on the
# log scale for lognormal, the SD on the response scale for gaussian, the
# precision for beta. `agg = TRUE` selects the cell-aggregated spec, whose
# eta_pos / y_pos are length 1. The dispersion formal is named `phi_pos` on the
# beta evaluators and `sigma_pos` on the other four, so arguments go
# positionally and this stays one call site.
eval_oc <- function(fam, eta_psi, eta_p, eta_pos, y_det, y_pos, disp,
                    curvature = "observed", agg = FALSE) {
    f <- switch(paste0(fam, if (agg) "_agg" else ""),
        lognormal     = cpp_eval_occu_cover_lognormal_cell,
        beta          = cpp_eval_occu_cover_beta_cell,
        gaussian      = cpp_eval_occu_cover_gaussian_cell,
        lognormal_agg = cpp_eval_occu_cover_lognormal_agg_cell,
        beta_agg      = cpp_eval_occu_cover_beta_agg_cell,
        gaussian_agg  = cpp_eval_occu_cover_gaussian_agg_cell)
    f(eta_psi, eta_p, eta_pos, y_det, y_pos, disp, curvature)
}

cell_ll_oc <- function(fam, d, eta_psi = d$eta_psi, eta_p = d$eta_p,
                       eta_pos = d$eta_pos, agg = FALSE) {
    eval_oc(fam, eta_psi, eta_p, eta_pos, d$y_det, d$y_pos, d$disp,
            agg = agg)$cell_ll
}

# R twin of PosPolicy::log_density, for the value checks.
oc_logdens <- function(fam, y, eta, disp) {
    switch(fam,
        beta     = dbeta(y, plogis(eta) * disp, (1 - plogis(eta)) * disp,
                         log = TRUE),
        gaussian = dnorm(y, eta, disp, log = TRUE),
        dlnorm(y, eta, disp, log = TRUE))
}

# Response-scale centre of the positive arm at eta: lognormal median, beta
# mean, gaussian mean.
oc_centre <- function(fam, eta) {
    switch(fam, beta = plogis(eta), gaussian = eta, exp(eta))
}

# Build a synthetic cell. `any_det = FALSE` forces the all-undetected branch.
# y_pos is drawn on the family's response scale at detected visits and left at
# 0 elsewhere, where the spec never reads it.
setup_cell <- function(seed, Jc, any_det, fam) {
    set.seed(seed)
    eta_psi <- rnorm(1, 0, 0.5)
    eta_p   <- rnorm(Jc, -0.3, 0.6)
    eta_pos <- if (fam == "beta") rnorm(Jc, 0.4, 0.4) else rnorm(Jc, 1.2, 0.4)
    if (any_det) {
        y_det <- rbinom(Jc, 1, 0.6)
        if (all(y_det == 0)) y_det[1] <- 1L
    } else {
        y_det <- rep(0L, Jc)
    }
    disp  <- if (fam == "beta") 12 else 0.35
    y_pos <- switch(fam,
        beta = vapply(seq_len(Jc), function(v) {
            if (y_det[v] == 1L) {
                mu <- plogis(eta_pos[v])
                rbeta(1, mu * disp, (1 - mu) * disp)
            } else 0
        }, numeric(1)),
        gaussian = rnorm(Jc, eta_pos, 0.3) * ifelse(y_det == 1, 1, 0),
        rlnorm(Jc, eta_pos, 0.3) * ifelse(y_det == 1, 1, 0))
    list(eta_psi = eta_psi, eta_p = eta_p, eta_pos = eta_pos,
         y_det = as.integer(y_det), y_pos = y_pos, disp = disp)
}

# Score and diagonal negative Hessian of every arm against central differences
# of the same cell density. Drives both granularities: the pos loop runs over
# the pos arm's own row count, which is Jc per-visit and 1 aggregated.
check_grad_diag_oc <- function(fam, d, label, curvature = "observed",
                               agg = FALSE) {
    tol1 <- oc_fam[[fam]]$tol1
    tol2 <- oc_fam[[fam]]$tol2
    res  <- eval_oc(fam, d$eta_psi, d$eta_p, d$eta_pos, d$y_det, d$y_pos,
                    d$disp, curvature = curvature, agg = agg)

    expect_equal(res$grad_psi,
                 fd1(function(h) cell_ll_oc(fam, d, eta_psi = d$eta_psi + h,
                                            agg = agg)),
                 tolerance = tol1, info = paste0(label, ": grad_psi"))
    expect_equal(res$neg_hess_psi,
                 -fd2(function(h) cell_ll_oc(fam, d, eta_psi = d$eta_psi + h,
                                             agg = agg)),
                 tolerance = tol2, info = paste0(label, ": neg_hess_psi"))

    for (v in seq_along(d$eta_p)) {
        expect_equal(res$grad_p[v],
                     fd1(function(h) cell_ll_oc(fam, d,
                                                eta_p = bump(d$eta_p, v, h),
                                                agg = agg)),
                     tolerance = tol1, info = paste0(label, ": grad_p[", v, "]"))
        expect_equal(res$neg_hess_p[v],
                     -fd2(function(h) cell_ll_oc(fam, d,
                                                 eta_p = bump(d$eta_p, v, h),
                                                 agg = agg)),
                     tolerance = tol2,
                     info = paste0(label, ": neg_hess_p[", v, "]"))
    }

    for (v in seq_along(d$eta_pos)) {
        expect_equal(res$grad_pos[v],
                     fd1(function(h) cell_ll_oc(fam, d,
                                                eta_pos = bump(d$eta_pos, v, h),
                                                agg = agg)),
                     tolerance = tol1,
                     info = paste0(label, ": grad_pos[", v, "]"))
        expect_equal(res$neg_hess_pos[v],
                     -fd2(function(h) cell_ll_oc(fam, d,
                                                 eta_pos = bump(d$eta_pos, v, h),
                                                 agg = agg)),
                     tolerance = tol2,
                     info = paste0(label, ": neg_hess_pos[", v, "]"))
    }

    res
}

# Cross-Hessian checks (-d2 cell_ll / d eta_a d eta_b). Nonzero only in the
# nodet branch, where the occupancy mixture couples psi to every p_v and the
# p_v to one another.
check_cross_oc <- function(fam, d, res, label, tol = 1e-3) {
    Jc <- length(d$eta_p)
    for (v in seq_len(Jc)) {
        fdv <- -fd_cross(function(a, b) cell_ll_oc(
            fam, d, eta_psi = d$eta_psi + a, eta_p = bump(d$eta_p, v, b)))
        expect_equal(res$cross_psi_p[v], fdv, tolerance = tol,
                     info = paste0(label, ": cross_psi_p[", v, "]"))
    }
    for (v in seq_len(Jc)) for (w in seq_len(Jc)) if (v != w) {
        fdv <- -fd_cross(function(a, b) cell_ll_oc(
            fam, d, eta_p = bump(bump(d$eta_p, v, a), w, b)))
        expect_equal(res$cross_p_p[v, w], fdv, tolerance = tol,
                     info = paste0(label, ": cross_p_p[", v, ",", w, "]"))
    }
}


# ---------------------------------------------------------------------------
# Per-visit cover.
# ---------------------------------------------------------------------------

for (fam in oc_fams) {

    seed0 <- oc_fam[[fam]]$seed

    test_that(paste0(fam, " det case: value, gradients + diagonal neg-Hess match FD"), {
        d   <- setup_cell(seed0 + 0L, Jc = 4L, any_det = TRUE, fam = fam)
        res <- check_grad_diag_oc(fam, d, label = paste0(fam, "-det"))

        # cell_ll against the closed-form det density: log psi + the detection
        # terms + one log f_pos per detected visit.
        det <- d$y_det == 1L
        p_v <- plogis(d$eta_p)
        expect_equal(res$cell_ll,
                     log(plogis(d$eta_psi)) +
                         sum(log(ifelse(det, p_v, 1 - p_v))) +
                         sum(oc_logdens(fam, d$y_pos[det], d$eta_pos[det], d$disp)),
                     tolerance = 1e-10, info = paste0(fam, "-det: cell_ll"))

        # The det branch factorises: every cross-Hessian buffer stays zero.
        expect_true(all(res$cross_psi_p == 0))
        expect_true(all(res$cross_p_p == 0))

        # z is known here, so the Fisher curvature is the observed Hessian and
        # has to clear the same FD checks on every arm, cover arm included.
        res_e <- check_grad_diag_oc(fam, d, label = paste0(fam, "-det-fisher"),
                                    curvature = "expected")
        expect_equal(res_e$neg_hess_psi, res$neg_hess_psi, tolerance = 1e-12)
        expect_equal(res_e$neg_hess_p,   res$neg_hess_p,   tolerance = 1e-12)
        expect_equal(res_e$neg_hess_pos, res$neg_hess_pos, tolerance = 1e-12)
        expect_true(all(res_e$cross_psi_p == 0))
        expect_true(all(res_e$cross_p_p == 0))
    })

    test_that(paste0(fam, " nodet case: value, gradients, diagonal + cross neg-Hess match FD"), {
        d   <- setup_cell(seed0 + 1L, Jc = 4L, any_det = FALSE, fam = fam)
        res <- check_grad_diag_oc(fam, d, label = paste0(fam, "-nodet"))
        check_cross_oc(fam, d, res, label = paste0(fam, "-nodet"))

        psi <- plogis(d$eta_psi)
        P0  <- prod(1 - plogis(d$eta_p))
        expect_equal(res$cell_ll, log(psi * P0 + (1 - psi)), tolerance = 1e-10,
                     info = paste0(fam, "-nodet: cell_ll"))
    })

    test_that(paste0(fam, " nodet Fisher curvature: block-diagonal PSD, shares the gradient"), {
        d     <- setup_cell(seed0 + 1L, Jc = 4L, any_det = FALSE, fam = fam)
        res_o <- eval_oc(fam, d$eta_psi, d$eta_p, d$eta_pos, d$y_det, d$y_pos,
                         d$disp, curvature = "observed")
        res_e <- eval_oc(fam, d$eta_psi, d$eta_p, d$eta_pos, d$y_det, d$y_pos,
                         d$disp, curvature = "expected")

        # Fisher scoring changes the curvature, not the score.
        expect_equal(res_e$grad_psi, res_o$grad_psi, tolerance = 1e-12)
        expect_equal(res_e$grad_p,   res_o$grad_p,   tolerance = 1e-12)

        # Complete-data Fisher diagonals: psi(1-psi) and gamma_c p_v(1-p_v),
        # gamma_c = P(z = 1 | all undetected) = psi P0 / L.
        psi     <- plogis(d$eta_psi)
        p       <- plogis(d$eta_p)
        P0      <- prod(1 - p)
        L       <- psi * P0 + (1 - psi)
        gamma_c <- psi * P0 / L
        expect_equal(res_e$neg_hess_psi, psi * (1 - psi), tolerance = 1e-10)
        for (v in seq_along(p)) {
            expect_equal(res_e$neg_hess_p[v], gamma_c * p[v] * (1 - p[v]),
                         tolerance = 1e-10,
                         info = paste0(fam, "-fisher: neg_hess_p[", v, "]"))
        }

        # Block-diagonal (no cross-Hessian) and PSD (non-negative diagonals).
        expect_true(all(res_e$cross_psi_p == 0))
        expect_true(all(res_e$cross_p_p == 0))
        expect_true(res_e$neg_hess_psi >= 0)
        expect_true(all(res_e$neg_hess_p >= 0))
    })

    test_that(paste0(fam, " nodet case with a single visit: derivs match FD"), {
        d   <- setup_cell(seed0 + 2L, Jc = 1L, any_det = FALSE, fam = fam)
        res <- check_grad_diag_oc(fam, d, label = paste0(fam, "-nodet-J1"))
        check_cross_oc(fam, d, res, label = paste0(fam, "-nodet-J1"))
    })

    test_that(paste0(fam, " det case with one detection out of five visits: derivs match FD"), {
        d <- setup_cell(seed0 + 3L, Jc = 5L, any_det = TRUE, fam = fam)
        d$y_det <- c(0L, 0L, 1L, 0L, 0L)
        # Cover observed exactly at the arm's response-scale centre.
        d$y_pos <- oc_centre(fam, d$eta_pos) * d$y_det
        check_grad_diag_oc(fam, d, label = paste0(fam, "-det-1of5"))
    })
}


# ---------------------------------------------------------------------------
# Cell-aggregated cover (tulpaObs#33). The pos arm carries ONE row per cell (the
# mean / median cover over the cell's detected visits), so the det branch adds a
# single log f_pos(ybar; eta_pos_cell). eta_pos / y_pos are length 1; eta_p and
# y_det stay length J.
# ---------------------------------------------------------------------------

setup_cell_agg <- function(seed, Jc, y_det, fam) {
    set.seed(seed)
    eta_psi <- rnorm(1, 0, 0.5)
    eta_p   <- rnorm(Jc, -0.2, 0.6)
    eta_pos <- if (fam == "beta") rnorm(1, 0.3, 0.4) else rnorm(1, 1.0, 0.4)
    disp    <- if (fam == "beta") 14 else 0.35
    y_pos   <- switch(fam,
        beta = {
            mu <- plogis(eta_pos)
            rbeta(1, mu * disp, (1 - mu) * disp)
        },
        gaussian = eta_pos + rnorm(1, 0, 0.2),
        exp(eta_pos + rnorm(1, 0, 0.2)))
    list(eta_psi = eta_psi, eta_p = eta_p, eta_pos = eta_pos,
         y_det = as.integer(y_det), y_pos = y_pos, disp = disp)
}

for (fam in oc_fams) {

    seed0 <- oc_fam[[fam]]$seed

    test_that(paste0("aggregated ", fam, " det case: pos derivs match FD, psi/p match the per-visit spec"), {
        d   <- setup_cell_agg(seed0 + 4L, Jc = 5L,
                              y_det = c(1L, 0L, 1L, 1L, 0L), fam = fam)
        res <- check_grad_diag_oc(fam, d, label = paste0(fam, "-agg-det"),
                                  agg = TRUE)
        expect_length(res$grad_pos, 1L)
        expect_length(res$neg_hess_pos, 1L)

        # Complete data again, so Fisher has to clear the same FD checks.
        res_e <- check_grad_diag_oc(fam, d, agg = TRUE, curvature = "expected",
                                    label = paste0(fam, "-agg-det-fisher"))
        expect_equal(res_e$neg_hess_pos, res$neg_hess_pos, tolerance = 1e-12)

        # psi / p block is independent of the cover granularity: equal to the
        # per-visit spec fed the same psi / p etas and the aggregate at each
        # detected visit. The densities then differ by the (n_det - 1) repeated
        # cover factors the per-visit spec adds.
        Jc <- length(d$eta_p)
        pv <- eval_oc(fam, d$eta_psi, d$eta_p, rep(d$eta_pos, Jc), d$y_det,
                      ifelse(d$y_det == 1L, d$y_pos, 0), d$disp)
        expect_equal(res$grad_psi,     pv$grad_psi,     tolerance = 1e-12)
        expect_equal(res$grad_p,       pv$grad_p,       tolerance = 1e-12)
        expect_equal(res$neg_hess_psi, pv$neg_hess_psi, tolerance = 1e-12)
        expect_equal(res$neg_hess_p,   pv$neg_hess_p,   tolerance = 1e-12)

        n_det <- sum(d$y_det == 1L)
        lf    <- oc_logdens(fam, d$y_pos, d$eta_pos, d$disp)
        expect_equal(res$cell_ll, pv$cell_ll - (n_det - 1) * lf,
                     tolerance = 1e-9, info = paste0(fam, "-agg: cell_ll"))
    })

    test_that(paste0("aggregated ", fam, " nodet case: pos arm contributes nothing"), {
        d   <- setup_cell_agg(seed0 + 5L, Jc = 4L, y_det = rep(0L, 4L), fam = fam)
        res <- eval_oc(fam, d$eta_psi, d$eta_p, d$eta_pos, d$y_det, d$y_pos,
                       d$disp, agg = TRUE)
        expect_equal(res$grad_pos[1], 0, tolerance = 1e-12)
        expect_equal(res$neg_hess_pos[1], 0, tolerance = 1e-12)

        # cell_ll equals the per-visit spec's nodet density (pos ignored either way).
        Jc <- length(d$eta_p)
        pv <- eval_oc(fam, d$eta_psi, d$eta_p, rep(d$eta_pos, Jc), d$y_det,
                      rep(0, Jc), d$disp)
        expect_equal(res$cell_ll, pv$cell_ll, tolerance = 1e-12)
    })
}


test_that("nodet derivatives are family-independent (pos arm idle)", {
    d   <- setup_cell(oc_fam$lognormal$seed + 1L, Jc = 4L, any_det = FALSE,
                      fam = "lognormal")
    ref <- eval_oc("lognormal", d$eta_psi, d$eta_p, d$eta_pos, d$y_det,
                   d$y_pos, d$disp)
    for (fam in setdiff(oc_fams, "lognormal")) {
        alt <- eval_oc(fam, d$eta_psi, d$eta_p, d$eta_pos, d$y_det, d$y_pos,
                       if (fam == "beta") 12 else d$disp)
        for (nm in c("cell_ll", "grad_psi", "grad_p", "neg_hess_psi",
                     "neg_hess_p", "cross_psi_p", "cross_p_p")) {
            expect_equal(alt[[nm]], ref[[nm]], tolerance = 1e-12,
                         info = paste0(fam, " vs lognormal: ", nm))
        }
    }
})

test_that("every occu_cover coupling spec is registered at package load", {
    for (nm in c("occu_cover_lognormal", "occu_cover_beta", "occu_cover_gaussian",
                 "occu_cover_lognormal_agg", "occu_cover_beta_agg",
                 "occu_cover_gaussian_agg")) {
        expect_true(tulpa:::cpp_cell_coupling_registry_has(nm), info = nm)
    }
})
