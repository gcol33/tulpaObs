# FD-checks for every closed-form derivative of the
# OccuMultiscaleCover{Lognormal,Beta}Coupling CellCouplingSpec
# (three-level occupancy + cover hurdle, gcol33/tulpaObs#29).
#
# The spec lives in src/cell_coupling_occu_multiscale_cover.h; the direct
# evaluators `cpp_eval_occu_multiscale_cover_{lognormal,beta}_cell` return the
# cell log-density at given etas + every nonzero derivative buffer the kernel
# reads (gradients; diagonal neg-Hess; cross-Hess (psi,theta), (psi,p),
# (theta,theta), (theta,p), (p,p)). We FD-check each against numerical 1st / 2nd
# derivatives of cell_ll, in three cell regimes:
#   * branch A, mixed: some plots detect, some do not (cell occupied)
#   * branch A, all plots detect (fully factorized -> zero cross)
#   * branch B: no detection anywhere (nested z-over-cell, a-over-plot mixture)

fd1 <- function(f, h = 1e-5) (f(+h) - f(-h)) / (2 * h)
fd2 <- function(f, h = 1e-4) (f(+h) - 2 * f(0) + f(-h)) / (h * h)
fd_cross <- function(g, h = 1e-3) {
    (g(+h, +h) - g(+h, -h) - g(-h, +h) + g(-h, -h)) / (4 * h * h)
}

# Family-agnostic evaluator wrappers. `fam` is "lognormal" or "beta".
eval_ms <- function(fam, eta_psi, eta_theta, eta_p, eta_pos,
                    y_det, y_pos, plot_sizes, disp, curvature = "observed") {
    if (fam == "beta") {
        cpp_eval_occu_multiscale_cover_beta_cell(
            eta_psi, eta_theta, eta_p, eta_pos, y_det, y_pos,
            plot_sizes, phi_pos = disp, curvature = curvature)
    } else {
        cpp_eval_occu_multiscale_cover_lognormal_cell(
            eta_psi, eta_theta, eta_p, eta_pos, y_det, y_pos,
            plot_sizes, sigma_pos = disp, curvature = curvature)
    }
}
cell_ll_ms <- function(fam, d, eta_psi = d$eta_psi, eta_theta = d$eta_theta,
                       eta_p = d$eta_p, eta_pos = d$eta_pos) {
    eval_ms(fam, eta_psi, eta_theta, eta_p, eta_pos,
            d$y_det, d$y_pos, d$plot_sizes, d$disp)$cell_ll
}

# Build a synthetic cell. `det_pattern` is a list (one vector per plot) of 0/1
# detection flags per visit; plot_sizes and Jc derive from it. y_pos drawn
# positive at detected visits (0 elsewhere) under the requested family.
setup_cell_ms <- function(seed, det_pattern, fam = "lognormal") {
    set.seed(seed)
    M          <- length(det_pattern)
    plot_sizes <- vapply(det_pattern, length, integer(1))
    Jc         <- sum(plot_sizes)
    y_det      <- as.integer(unlist(det_pattern))
    eta_psi    <- rnorm(1, 0.2, 0.5)
    eta_theta  <- rnorm(M, 0.1, 0.6)
    eta_p      <- rnorm(Jc, -0.2, 0.6)
    if (fam == "beta") {
        eta_pos <- rnorm(Jc, 0.3, 0.4)
        disp    <- 12
        y_pos   <- vapply(seq_len(Jc), function(v) {
            if (y_det[v] == 1L) {
                mu <- plogis(eta_pos[v]); rbeta(1, mu * disp, (1 - mu) * disp)
            } else 0
        }, numeric(1))
    } else {
        eta_pos <- rnorm(Jc, 1.1, 0.4)
        disp    <- 0.35
        y_pos   <- rlnorm(Jc, eta_pos, 0.3) * ifelse(y_det == 1, 1, 0)
    }
    list(eta_psi = eta_psi, eta_theta = eta_theta, eta_p = eta_p,
         eta_pos = eta_pos, y_det = y_det, y_pos = y_pos,
         plot_sizes = as.integer(plot_sizes), disp = disp, M = M, Jc = Jc)
}

# Perturb-one-eta helpers returning modified eta vectors.
bump <- function(vec, i, h) { vec[i] <- vec[i] + h; vec }

check_grad_diag_ms <- function(fam, d, label, tol1 = 1e-5, tol2 = 1e-3) {
    res <- eval_ms(fam, d$eta_psi, d$eta_theta, d$eta_p, d$eta_pos,
                   d$y_det, d$y_pos, d$plot_sizes, d$disp)

    # psi
    expect_equal(res$grad_psi,
                 fd1(function(h) cell_ll_ms(fam, d, eta_psi = d$eta_psi + h)),
                 tolerance = tol1, info = paste0(label, ": grad_psi"))
    expect_equal(res$neg_hess_psi,
                 -fd2(function(h) cell_ll_ms(fam, d, eta_psi = d$eta_psi + h)),
                 tolerance = tol2, info = paste0(label, ": neg_hess_psi"))

    # theta
    for (j in seq_len(d$M)) {
        expect_equal(res$grad_theta[j],
                     fd1(function(h) cell_ll_ms(fam, d, eta_theta = bump(d$eta_theta, j, h))),
                     tolerance = tol1, info = paste0(label, ": grad_theta[", j, "]"))
        expect_equal(res$neg_hess_theta[j],
                     -fd2(function(h) cell_ll_ms(fam, d, eta_theta = bump(d$eta_theta, j, h))),
                     tolerance = tol2, info = paste0(label, ": neg_hess_theta[", j, "]"))
    }

    # p and pos
    for (v in seq_len(d$Jc)) {
        expect_equal(res$grad_p[v],
                     fd1(function(h) cell_ll_ms(fam, d, eta_p = bump(d$eta_p, v, h))),
                     tolerance = tol1, info = paste0(label, ": grad_p[", v, "]"))
        expect_equal(res$neg_hess_p[v],
                     -fd2(function(h) cell_ll_ms(fam, d, eta_p = bump(d$eta_p, v, h))),
                     tolerance = tol2, info = paste0(label, ": neg_hess_p[", v, "]"))
        expect_equal(res$grad_pos[v],
                     fd1(function(h) cell_ll_ms(fam, d, eta_pos = bump(d$eta_pos, v, h))),
                     tolerance = tol1, info = paste0(label, ": grad_pos[", v, "]"))
        expect_equal(res$neg_hess_pos[v],
                     -fd2(function(h) cell_ll_ms(fam, d, eta_pos = bump(d$eta_pos, v, h))),
                     tolerance = tol2, info = paste0(label, ": neg_hess_pos[", v, "]"))
    }
    res
}

# Cross-Hessian FD checks (-d2 cell_ll / d eta_a d eta_b).
check_cross_ms <- function(fam, d, res, label, tol = 2e-3) {
    # (psi, theta_j)
    for (j in seq_len(d$M)) {
        fdv <- -fd_cross(function(a, b) cell_ll_ms(
            fam, d, eta_psi = d$eta_psi + a, eta_theta = bump(d$eta_theta, j, b)))
        expect_equal(res$cross_psi_theta[j], fdv, tolerance = tol,
                     info = paste0(label, ": cross_psi_theta[", j, "]"))
    }
    # (psi, p_v)
    for (v in seq_len(d$Jc)) {
        fdv <- -fd_cross(function(a, b) cell_ll_ms(
            fam, d, eta_psi = d$eta_psi + a, eta_p = bump(d$eta_p, v, b)))
        expect_equal(res$cross_psi_p[v], fdv, tolerance = tol,
                     info = paste0(label, ": cross_psi_p[", v, "]"))
    }
    # (theta_j, theta_k), j != k
    for (j in seq_len(d$M)) for (k in seq_len(d$M)) if (j != k) {
        fdv <- -fd_cross(function(a, b) cell_ll_ms(
            fam, d, eta_theta = bump(bump(d$eta_theta, j, a), k, b)))
        expect_equal(res$cross_theta_theta[j, k], fdv, tolerance = tol,
                     info = paste0(label, ": cross_theta_theta[", j, ",", k, "]"))
    }
    # (theta_j, p_v)
    for (j in seq_len(d$M)) for (v in seq_len(d$Jc)) {
        fdv <- -fd_cross(function(a, b) cell_ll_ms(
            fam, d, eta_theta = bump(d$eta_theta, j, a), eta_p = bump(d$eta_p, v, b)))
        expect_equal(res$cross_theta_p[j, v], fdv, tolerance = tol,
                     info = paste0(label, ": cross_theta_p[", j, ",", v, "]"))
    }
    # (p_v, p_w), v != w
    for (v in seq_len(d$Jc)) for (w in seq_len(d$Jc)) if (v != w) {
        fdv <- -fd_cross(function(a, b) cell_ll_ms(
            fam, d, eta_p = bump(bump(d$eta_p, v, a), w, b)))
        expect_equal(res$cross_p_p[v, w], fdv, tolerance = tol,
                     info = paste0(label, ": cross_p_p[", v, ",", w, "]"))
    }
}

for (fam in c("lognormal", "beta")) {

    test_that(paste0(fam, " branch B (no detection): grad/diag/cross match FD"), {
        # 3 plots, 2/1/2 visits, none detect -> nested z/a mixture, all cross
        # blocks active.
        d <- setup_cell_ms(seed = 11L,
                           det_pattern = list(c(0, 0), c(0), c(0, 0)), fam = fam)
        res <- check_grad_diag_ms(fam, d, label = paste0(fam, "-B"))
        check_cross_ms(fam, d, res, label = paste0(fam, "-B"))
    })

    test_that(paste0(fam, " branch A mixed (some plots detect): grad/diag/cross match FD"), {
        # plot 1 detects, plot 2 does not, plot 3 detects: log psi separates,
        # the nodet plot contributes a within-plot (theta, p) mixture.
        d <- setup_cell_ms(seed = 22L,
                           det_pattern = list(c(1, 0), c(0, 0), c(0, 1)), fam = fam)
        res <- check_grad_diag_ms(fam, d, label = paste0(fam, "-Amix"))
        check_cross_ms(fam, d, res, label = paste0(fam, "-Amix"))
    })

    test_that(paste0(fam, " branch A all-detect: factorized, every cross is zero"), {
        d <- setup_cell_ms(seed = 33L,
                           det_pattern = list(c(1, 1), c(1), c(1, 0)), fam = fam)
        # plot 3 visit 2 undetected but plot 3 still detects (visit 1) -> det plot
        res <- check_grad_diag_ms(fam, d, label = paste0(fam, "-Aall"))
        expect_true(all(res$cross_psi_theta == 0))
        expect_true(all(res$cross_psi_p == 0))
        expect_true(all(res$cross_theta_theta == 0))
        expect_true(all(res$cross_theta_p == 0))
        expect_true(all(res$cross_p_p == 0))
    })

    test_that(paste0(fam, " single plot single visit: derivs match FD"), {
        d <- setup_cell_ms(seed = 44L, det_pattern = list(c(0)), fam = fam)
        res <- check_grad_diag_ms(fam, d, label = paste0(fam, "-1x1"))
        check_cross_ms(fam, d, res, label = paste0(fam, "-1x1"))
    })
}

test_that("Expected (Fisher) curvature: block-diagonal, PSD, shares the gradient (branch B)", {
    d <- setup_cell_ms(seed = 55L,
                       det_pattern = list(c(0, 0), c(0), c(0, 0)), fam = "beta")
    res_o <- eval_ms("beta", d$eta_psi, d$eta_theta, d$eta_p, d$eta_pos,
                     d$y_det, d$y_pos, d$plot_sizes, d$disp, curvature = "observed")
    res_e <- eval_ms("beta", d$eta_psi, d$eta_theta, d$eta_p, d$eta_pos,
                     d$y_det, d$y_pos, d$plot_sizes, d$disp, curvature = "expected")

    # Scores identical (Fisher changes curvature, not the gradient).
    expect_equal(res_e$grad_psi,   res_o$grad_psi,   tolerance = 1e-12)
    expect_equal(res_e$grad_theta, res_o$grad_theta, tolerance = 1e-12)
    expect_equal(res_e$grad_p,     res_o$grad_p,     tolerance = 1e-12)

    # Complete-data Fisher diagonals: psi(1-psi); gamma_c theta_j(1-theta_j);
    # gamma_c gamma_j p(1-p), gamma_c = psi M / L, gamma_j = theta_j P0_j / m_j.
    psi   <- plogis(d$eta_psi)
    theta <- plogis(d$eta_theta)
    p     <- plogis(d$eta_p)
    off   <- c(0, cumsum(d$plot_sizes))
    P0    <- vapply(seq_len(d$M), function(j)
        prod(1 - p[(off[j] + 1):off[j + 1]]), numeric(1))
    m     <- theta * P0 + (1 - theta)
    M     <- prod(m)
    L     <- psi * M + (1 - psi)
    gamma_c <- psi * M / L
    expect_equal(res_e$neg_hess_psi, psi * (1 - psi), tolerance = 1e-9)
    for (j in seq_len(d$M)) {
        expect_equal(res_e$neg_hess_theta[j], gamma_c * theta[j] * (1 - theta[j]),
                     tolerance = 1e-9, info = paste0("fisher theta[", j, "]"))
        gamma_j <- theta[j] * P0[j] / m[j]
        for (v in (off[j] + 1):off[j + 1]) {
            expect_equal(res_e$neg_hess_p[v],
                         gamma_c * gamma_j * p[v] * (1 - p[v]),
                         tolerance = 1e-9, info = paste0("fisher p[", v, "]"))
        }
    }
    # Block-diagonal (no cross) and PSD (non-negative diagonals).
    expect_true(all(res_e$cross_psi_theta == 0))
    expect_true(all(res_e$cross_psi_p == 0))
    expect_true(all(res_e$cross_theta_theta == 0))
    expect_true(all(res_e$cross_theta_p == 0))
    expect_true(all(res_e$cross_p_p == 0))
    expect_true(res_e$neg_hess_psi >= 0)
    expect_true(all(res_e$neg_hess_theta >= 0))
    expect_true(all(res_e$neg_hess_p >= 0))
})

test_that("branch B nodet derivatives are family-independent (pos arm idle)", {
    d_ln <- setup_cell_ms(seed = 66L,
                          det_pattern = list(c(0, 0), c(0)), fam = "lognormal")
    res_ln <- eval_ms("lognormal", d_ln$eta_psi, d_ln$eta_theta, d_ln$eta_p,
                      d_ln$eta_pos, d_ln$y_det, d_ln$y_pos, d_ln$plot_sizes, d_ln$disp)
    res_bt <- eval_ms("beta", d_ln$eta_psi, d_ln$eta_theta, d_ln$eta_p,
                      d_ln$eta_pos, d_ln$y_det, d_ln$y_pos, d_ln$plot_sizes, 12)
    expect_equal(res_ln$cell_ll,       res_bt$cell_ll,       tolerance = 1e-12)
    expect_equal(res_ln$grad_psi,      res_bt$grad_psi,      tolerance = 1e-12)
    expect_equal(res_ln$grad_theta,    res_bt$grad_theta,    tolerance = 1e-12)
    expect_equal(res_ln$grad_p,        res_bt$grad_p,        tolerance = 1e-12)
    expect_equal(res_ln$neg_hess_psi,  res_bt$neg_hess_psi,  tolerance = 1e-12)
    expect_equal(res_ln$cross_p_p,     res_bt$cross_p_p,     tolerance = 1e-12)
})

test_that("multiscale specs register and reduce to 2-level when theta -> 1, one plot", {
    # One plot holding all visits, theta -> +Inf (availability = 1): the plot
    # layer collapses and the cell density matches the 2-level occu_cover.
    set.seed(77L)
    Jc <- 4L
    eta_psi <- 0.3; eta_p <- rnorm(Jc, -0.2, 0.5); eta_pos <- rnorm(Jc, 1.0, 0.3)
    y_det <- c(1L, 0L, 1L, 0L)
    y_pos <- rlnorm(Jc, eta_pos, 0.3) * ifelse(y_det == 1, 1, 0)
    ms <- cpp_eval_occu_multiscale_cover_lognormal_cell(
        eta_psi, eta_theta = 30, eta_p = eta_p, eta_pos = eta_pos,
        y_det = y_det, y_pos = y_pos, plot_sizes = Jc, sigma_pos = 0.35)
    oc <- cpp_eval_occu_cover_lognormal_cell(
        eta_psi, eta_p = eta_p, eta_pos = eta_pos,
        y_det = y_det, y_pos = y_pos, sigma_pos = 0.35)
    expect_equal(ms$cell_ll, oc$cell_ll, tolerance = 1e-8)
    expect_equal(ms$grad_p,  oc$grad_p,  tolerance = 1e-6)
    expect_equal(ms$grad_pos, oc$grad_pos, tolerance = 1e-6)
})
