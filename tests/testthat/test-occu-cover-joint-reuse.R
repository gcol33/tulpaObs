# Factor reuse + grad-only + parallel outer grid through the real occu_cover
# beta cell-coupling spec (gcol33/tulpa#46). Exercises the sparse joint path
# (force.sparse = TRUE) where the inner_refresh / grad-only / parallel levers
# live: the beta arm honours CellDerivs::grad_only (skips its trigamma on a
# reuse step), so the converged joint fit must be invariant to inner.refresh and
# to n.threads.outer.

.fit_joint_beta <- function(seed, inner.refresh = 1L, n.threads.outer = 1L) {
    N <- 30L; J <- 4L
    adj <- matrix(0L, N, N)
    for (s in seq_len(N)) {
        if (s > 1L) adj[s, s - 1L] <- 1L
        if (s < N)  adj[s, s + 1L] <- 1L
    }
    sim <- simulate_occu_cover(N = N, J = J, positive = "beta", phi = 25,
                               adj = adj, sigma = 0.8, alpha = 1.0, seed = seed)
    long <- data.frame(
        site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
        y = as.vector(t(sim$y)),
        det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1)
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                    det.covs = c("det_cov1", "pos_cov1"))
    cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
    y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
    suppressWarnings(tobs(
        formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
        family = occu_cover("beta"),
        detection = ~ det_cov1, positive = ~ pos_cov1,
        y = od$y, y_pos = y_pos, visits = od$det.covs,
        method = "nested_laplace",
        control = list(verbose = FALSE, max.iter = 500L, engine = "joint_coupled",
                       force.sparse = TRUE, adaptive.grid = FALSE,
                       inner.refresh = inner.refresh,
                       n.threads.outer = n.threads.outer)))
}

.expect_joint_equiv <- function(a, b, lm_tol = 1e-5, mode_tol = 1e-4) {
    ja <- a$joint_fit; jb <- b$joint_fit
    expect_equal(as.numeric(jb$log_marginal), as.numeric(ja$log_marginal),
                 tolerance = lm_tol)
    expect_equal(as.numeric(jb$modes), as.numeric(ja$modes), tolerance = mode_tol)
}

test_that("beta cover joint fit is invariant to inner.refresh (grad-only reuse)", {
    skip_on_cran()
    skip_if_fast()
    fit1 <- .fit_joint_beta(6789L, inner.refresh = 1L)
    fit3 <- .fit_joint_beta(6789L, inner.refresh = 3L)
    .expect_joint_equiv(fit1, fit3)
})

test_that("beta cover joint fit is invariant to n.threads.outer (parallel sparse)", {
    skip_on_cran()
    skip_if_fast()
    skip_if_not(parallel::detectCores() >= 2L, "needs multi-core")
    n_outer <- max(2L, min(4L, parallel::detectCores() - 1L))
    fit_s <- .fit_joint_beta(6789L, n.threads.outer = 1L)
    fit_p <- .fit_joint_beta(6789L, n.threads.outer = n_outer)
    .expect_joint_equiv(fit_s, fit_p)
})
