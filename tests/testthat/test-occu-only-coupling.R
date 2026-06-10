# FD-checks for every closed-form derivative of the OccuOnlyCoupling
# CellCouplingSpec (gcol33/tulpaObs#81): the single-arm occupancy mixture (psi/p
# arms, no cover arm) the standalone occu() SVC bar fits through the joint
# direct-grid engine.
#
# The spec lives in src/cell_coupling_occu_only.h; the direct evaluator
# `cpp_eval_occu_only_cell` returns the cell log-density at given etas + every
# nonzero derivative buffer the kernel reads (grad_psi, grad_p, diagonal
# neg-Hess, and the nodet (psi, p_v) / (p_v, p_w) cross-Hessians). We FD-check
# each against numerical 1st / 2nd derivatives of cell_ll. The occupancy /
# detection derivatives are the SAME shared blocks the occu_cover spec uses
# (occu_det_psi_p_block / occu_nodet_block); this test confirms the cover arm
# was dropped without altering them.

occu_cell_ll_at <- function(eta_psi, eta_p, y_det) {
    cpp_eval_occu_only_cell(eta_psi = eta_psi, eta_p = eta_p, y_det = y_det)$cell_ll
}

ofd1 <- function(f, h = 1e-5) (f(+h) - f(-h)) / (2 * h)
ofd2 <- function(f, h = 1e-4) (f(+h) - 2 * f(0) + f(-h)) / (h * h)
ofd_cross <- function(g, h = 1e-3) {
    (g(+h, +h) - g(+h, -h) - g(-h, +h) + g(-h, -h)) / (4 * h * h)
}

occu_setup_cell <- function(seed, Jc, any_det) {
    set.seed(seed)
    eta_psi <- rnorm(1, 0, 0.5)
    eta_p   <- rnorm(Jc, -0.3, 0.6)
    if (any_det) {
        y_det <- rbinom(Jc, 1, 0.6)
        if (all(y_det == 0)) y_det[1] <- 1L
    } else {
        y_det <- rep(0L, Jc)
    }
    list(eta_psi = eta_psi, eta_p = eta_p, y_det = as.integer(y_det))
}

occu_check_grad_diag <- function(d, label, tol1 = 1e-6, tol2 = 1e-4) {
    res <- cpp_eval_occu_only_cell(d$eta_psi, d$eta_p, d$y_det)

    grad_psi_fd <- ofd1(function(h)
        occu_cell_ll_at(d$eta_psi + h, d$eta_p, d$y_det))
    expect_equal(res$grad_psi, grad_psi_fd, tolerance = tol1,
                 info = paste0(label, ": grad_psi"))

    nh_psi_fd <- -ofd2(function(h)
        occu_cell_ll_at(d$eta_psi + h, d$eta_p, d$y_det))
    expect_equal(res$neg_hess_psi, nh_psi_fd, tolerance = tol2,
                 info = paste0(label, ": neg_hess_psi"))

    Jc <- length(d$eta_p)
    for (v in seq_len(Jc)) {
        grad_p_fd <- ofd1(function(h) {
            ep <- d$eta_p; ep[v] <- ep[v] + h
            occu_cell_ll_at(d$eta_psi, ep, d$y_det)
        })
        expect_equal(res$grad_p[v], grad_p_fd, tolerance = tol1,
                     info = paste0(label, ": grad_p[", v, "]"))

        nh_p_fd <- -ofd2(function(h) {
            ep <- d$eta_p; ep[v] <- ep[v] + h
            occu_cell_ll_at(d$eta_psi, ep, d$y_det)
        })
        expect_equal(res$neg_hess_p[v], nh_p_fd, tolerance = tol2,
                     info = paste0(label, ": neg_hess_p[", v, "]"))
    }
    res
}

test_that("det case: gradients + diagonal neg-Hess match FD; cross-Hess zero", {
    d   <- occu_setup_cell(seed = 101L, Jc = 4L, any_det = TRUE)
    res <- occu_check_grad_diag(d, label = "det")
    expect_true(all(res$cross_psi_p == 0))
    expect_true(all(res$cross_p_p == 0))
})

test_that("nodet case: gradients + diagonal neg-Hess + cross-Hess match FD", {
    d   <- occu_setup_cell(seed = 202L, Jc = 4L, any_det = FALSE)
    res <- occu_check_grad_diag(d, label = "nodet")
    Jc  <- length(d$eta_p)

    for (v in seq_len(Jc)) {
        cross_fd <- -ofd_cross(function(a, b) {
            ep <- d$eta_p; ep[v] <- ep[v] + b
            occu_cell_ll_at(d$eta_psi + a, ep, d$y_det)
        })
        expect_equal(res$cross_psi_p[v], cross_fd, tolerance = 1e-3,
                     info = paste0("nodet: cross_psi_p[", v, "]"))
    }
    for (v in seq_len(Jc)) {
        for (w in seq_len(Jc)) {
            if (v == w) next
            cross_pp_fd <- -ofd_cross(function(a, b) {
                ep <- d$eta_p; ep[v] <- ep[v] + a; ep[w] <- ep[w] + b
                occu_cell_ll_at(d$eta_psi, ep, d$y_det)
            })
            expect_equal(res$cross_p_p[v, w], cross_pp_fd, tolerance = 1e-3,
                         info = paste0("nodet: cross_p_p[", v, ",", w, "]"))
        }
    }
})

test_that("nodet Fisher curvature is block-diagonal PSD and shares the gradient", {
    d <- occu_setup_cell(seed = 202L, Jc = 4L, any_det = FALSE)
    res_o <- cpp_eval_occu_only_cell(d$eta_psi, d$eta_p, d$y_det,
                                     curvature = "observed")
    res_e <- cpp_eval_occu_only_cell(d$eta_psi, d$eta_p, d$y_det,
                                     curvature = "expected")
    expect_equal(res_e$grad_psi, res_o$grad_psi, tolerance = 1e-12)
    expect_equal(res_e$grad_p,   res_o$grad_p,   tolerance = 1e-12)

    psi     <- plogis(d$eta_psi)
    p       <- plogis(d$eta_p)
    P0      <- prod(1 - p)
    L       <- psi * P0 + (1 - psi)
    gamma_c <- psi * P0 / L
    expect_equal(res_e$neg_hess_psi, psi * (1 - psi), tolerance = 1e-10)
    for (v in seq_along(p)) {
        expect_equal(res_e$neg_hess_p[v], gamma_c * p[v] * (1 - p[v]),
                     tolerance = 1e-10, info = paste0("fisher: neg_hess_p[", v, "]"))
    }
    expect_true(all(res_e$cross_psi_p == 0))
    expect_true(all(res_e$cross_p_p == 0))
    expect_true(res_e$neg_hess_psi >= 0)
    expect_true(all(res_e$neg_hess_p >= 0))
})

test_that("det Fisher equals observed (complete-data, z known); single-visit nodet", {
    d <- occu_setup_cell(seed = 101L, Jc = 4L, any_det = TRUE)
    res_o <- cpp_eval_occu_only_cell(d$eta_psi, d$eta_p, d$y_det,
                                     curvature = "observed")
    res_e <- cpp_eval_occu_only_cell(d$eta_psi, d$eta_p, d$y_det,
                                     curvature = "expected")
    expect_equal(res_e$neg_hess_psi, res_o$neg_hess_psi, tolerance = 1e-12)
    expect_equal(res_e$neg_hess_p,   res_o$neg_hess_p,   tolerance = 1e-12)

    d1 <- occu_setup_cell(seed = 303L, Jc = 1L, any_det = FALSE)
    occu_check_grad_diag(d1, label = "nodet-J1")
})

test_that("occu_only cell density matches the occu_cover spec with the cover arm muted", {
    # With NO detected visit (any_det = FALSE) the occu_cover cell density has no
    # cover contribution, so its psi/p log-density must equal the occu_only one
    # at the same occupancy / detection etas (the shared nodet block).
    d <- occu_setup_cell(seed = 404L, Jc = 5L, any_det = FALSE)
    ll_occu  <- cpp_eval_occu_only_cell(d$eta_psi, d$eta_p, d$y_det)$cell_ll
    ll_cover <- cpp_eval_occu_cover_lognormal_cell(
        eta_psi = d$eta_psi, eta_p = d$eta_p,
        eta_pos = rep(0, length(d$eta_p)),
        y_det = d$y_det, y_pos = rep(0, length(d$eta_p)), sigma_pos = 0.4)$cell_ll
    expect_equal(ll_occu, ll_cover, tolerance = 1e-12)
})
