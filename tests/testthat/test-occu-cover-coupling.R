# FD-checks for every closed-form derivative of the
# OccuCoverLognormalCoupling CellCouplingSpec
# (tulpaObs consumer of gcol33/tulpa#32 Layer B.2).
#
# The spec lives in src/cell_coupling_occu_cover.h; the direct evaluator
# `cpp_eval_occu_cover_lognormal_cell` returns the cell log-density at
# given etas + every nonzero derivative buffer the kernel would read
# (gradients, diagonal neg-Hess, cross-Hess (psi, p_v), (p_v, p_w)). We
# FD-check each against numerical 1st / 2nd derivatives of cell_ll.

cell_ll_at <- function(eta_psi, eta_p, eta_pos,
                       y_det, y_pos, sigma_pos) {
    cpp_eval_occu_cover_lognormal_cell(
        eta_psi  = eta_psi,
        eta_p    = eta_p,
        eta_pos  = eta_pos,
        y_det    = y_det,
        y_pos    = y_pos,
        sigma_pos = sigma_pos
    )$cell_ll
}

fd1 <- function(f, h = 1e-5) {
    (f(+h) - f(-h)) / (2 * h)
}

fd2 <- function(f, h = 1e-4) {
    (f(+h) - 2 * f(0) + f(-h)) / (h * h)
}

fd_cross <- function(g, h = 1e-3) {
    # g(a, b) returns cell_ll at offset (a, b)
    (g(+h, +h) - g(+h, -h) - g(-h, +h) + g(-h, -h)) / (4 * h * h)
}

setup_cell <- function(seed, Jc, any_det) {
    set.seed(seed)
    eta_psi <- rnorm(1, 0, 0.5)
    eta_p   <- rnorm(Jc, -0.3, 0.6)
    eta_pos <- rnorm(Jc, 1.2, 0.4)
    if (any_det) {
        y_det <- rbinom(Jc, 1, 0.6)
        if (all(y_det == 0)) y_det[1] <- 1L
    } else {
        y_det <- rep(0L, Jc)
    }
    y_pos <- rlnorm(Jc, eta_pos, 0.3) * ifelse(y_det == 1, 1, 0)
    sigma_pos <- 0.35
    list(eta_psi = eta_psi,
         eta_p   = eta_p,
         eta_pos = eta_pos,
         y_det   = as.integer(y_det),
         y_pos   = y_pos,
         sigma_pos = sigma_pos)
}

check_grad_diag <- function(d, label, tol1 = 1e-6, tol2 = 1e-4) {
    res <- cpp_eval_occu_cover_lognormal_cell(
        eta_psi   = d$eta_psi,
        eta_p     = d$eta_p,
        eta_pos   = d$eta_pos,
        y_det     = d$y_det,
        y_pos     = d$y_pos,
        sigma_pos = d$sigma_pos
    )

    grad_psi_fd <- fd1(function(h) cell_ll_at(
        d$eta_psi + h, d$eta_p, d$eta_pos,
        d$y_det, d$y_pos, d$sigma_pos
    ))
    expect_equal(res$grad_psi, grad_psi_fd,
                 tolerance = tol1, info = paste0(label, ": grad_psi"))

    nh_psi_fd <- -fd2(function(h) cell_ll_at(
        d$eta_psi + h, d$eta_p, d$eta_pos,
        d$y_det, d$y_pos, d$sigma_pos
    ))
    expect_equal(res$neg_hess_psi, nh_psi_fd,
                 tolerance = tol2, info = paste0(label, ": neg_hess_psi"))

    Jc <- length(d$eta_p)
    for (v in seq_len(Jc)) {
        grad_p_fd <- fd1(function(h) {
            eta_p_perturb <- d$eta_p
            eta_p_perturb[v] <- eta_p_perturb[v] + h
            cell_ll_at(d$eta_psi, eta_p_perturb, d$eta_pos,
                       d$y_det, d$y_pos, d$sigma_pos)
        })
        expect_equal(res$grad_p[v], grad_p_fd,
                     tolerance = tol1,
                     info = paste0(label, ": grad_p[", v, "]"))

        nh_p_fd <- -fd2(function(h) {
            eta_p_perturb <- d$eta_p
            eta_p_perturb[v] <- eta_p_perturb[v] + h
            cell_ll_at(d$eta_psi, eta_p_perturb, d$eta_pos,
                       d$y_det, d$y_pos, d$sigma_pos)
        })
        expect_equal(res$neg_hess_p[v], nh_p_fd,
                     tolerance = tol2,
                     info = paste0(label, ": neg_hess_p[", v, "]"))

        grad_pos_fd <- fd1(function(h) {
            eta_pos_perturb <- d$eta_pos
            eta_pos_perturb[v] <- eta_pos_perturb[v] + h
            cell_ll_at(d$eta_psi, d$eta_p, eta_pos_perturb,
                       d$y_det, d$y_pos, d$sigma_pos)
        })
        expect_equal(res$grad_pos[v], grad_pos_fd,
                     tolerance = tol1,
                     info = paste0(label, ": grad_pos[", v, "]"))

        nh_pos_fd <- -fd2(function(h) {
            eta_pos_perturb <- d$eta_pos
            eta_pos_perturb[v] <- eta_pos_perturb[v] + h
            cell_ll_at(d$eta_psi, d$eta_p, eta_pos_perturb,
                       d$y_det, d$y_pos, d$sigma_pos)
        })
        expect_equal(res$neg_hess_pos[v], nh_pos_fd,
                     tolerance = tol2,
                     info = paste0(label, ": neg_hess_pos[", v, "]"))
    }

    res
}

test_that("det case: gradients + diagonal neg-Hess match FD", {
    d <- setup_cell(seed = 101L, Jc = 4L, any_det = TRUE)
    res <- check_grad_diag(d, label = "det")
    # det case: every cross-Hessian buffer is zero
    expect_true(all(res$cross_psi_p == 0))
    expect_true(all(res$cross_p_p == 0))
})

test_that("nodet case: gradients + diagonal neg-Hess match FD", {
    d <- setup_cell(seed = 202L, Jc = 4L, any_det = FALSE)
    res <- check_grad_diag(d, label = "nodet")

    # Cross-Hess (psi, p_v) FD vs closed form
    Jc <- length(d$eta_p)
    for (v in seq_len(Jc)) {
        nh_cross_psi_pv_fd <- -fd_cross(function(a, b) {
            eta_p_perturb <- d$eta_p
            eta_p_perturb[v] <- eta_p_perturb[v] + b
            cell_ll_at(d$eta_psi + a, eta_p_perturb, d$eta_pos,
                       d$y_det, d$y_pos, d$sigma_pos)
        })
        expect_equal(res$cross_psi_p[v], nh_cross_psi_pv_fd,
                     tolerance = 1e-3,
                     info = paste0("nodet: cross_psi_p[", v, "]"))
    }

    # Cross-Hess (p_v, p_w) FD vs closed form
    for (v in seq_len(Jc)) {
        for (w in seq_len(Jc)) {
            if (v == w) {
                next  # diagonal is in neg_hess_p, checked above
            }
            nh_cross_pp_fd <- -fd_cross(function(a, b) {
                eta_p_perturb <- d$eta_p
                eta_p_perturb[v] <- eta_p_perturb[v] + a
                eta_p_perturb[w] <- eta_p_perturb[w] + b
                cell_ll_at(d$eta_psi, eta_p_perturb, d$eta_pos,
                           d$y_det, d$y_pos, d$sigma_pos)
            })
            expect_equal(res$cross_p_p[v, w], nh_cross_pp_fd,
                         tolerance = 1e-3,
                         info = paste0("nodet: cross_p_p[", v, ",", w, "]"))
        }
    }
})

test_that("nodet case with single visit: derivs match FD", {
    d <- setup_cell(seed = 303L, Jc = 1L, any_det = FALSE)
    check_grad_diag(d, label = "nodet-J1")
})

test_that("det case with single detection out of many visits: derivs match FD", {
    d <- setup_cell(seed = 404L, Jc = 5L, any_det = TRUE)
    d$y_det <- c(0L, 0L, 1L, 0L, 0L)
    d$y_pos <- ifelse(d$y_det == 1, exp(d$eta_pos), 0)
    check_grad_diag(d, label = "det-1of5")
})

test_that("spec is registered under occu_cover_lognormal at package load", {
    expect_true(tulpa:::cpp_cell_coupling_registry_has("occu_cover_lognormal"))
})
