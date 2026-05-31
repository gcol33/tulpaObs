# predict() for the joint occu_cover fit (gcol33/tulpaObs#22).
#
# The joint_coupled engine builds an ICAR shared field internally, couples it
# onto the cover arm with the alpha axis, and integrates (sigma, alpha) on the
# outer grid. predict() samples the joint latent via tulpa::tulpa_posterior_draws
# and marginalizes every derived quantity per draw.
#
# Fixture mirrors test-occu-cover-joint-coupled.R: occupancy is cell-level
# (formula on cell_dat), detection + cover are visit-level (od$det.covs).

.ocp_chain_adj <- function(N) {
    adj <- matrix(0L, N, N)
    for (s in seq_len(N)) {
        if (s > 1L) adj[s, s - 1L] <- 1L
        if (s < N)  adj[s, s + 1L] <- 1L
    }
    adj
}

# Build + fit a spatial occu_cover through the joint_coupled engine. A cell-level
# `year` covariate is added to the occupancy arm so type = "change" between two
# years moves occupancy (and, through the shared field, expected cover).
.ocp_build_fit <- function(N = 40L, J = 5L, seed = 202L) {
    adj <- .ocp_chain_adj(N)
    sim <- simulate_occu_cover(
        N = N, J = J, positive = "lognormal",
        adj = adj, sigma = 0.8, alpha = 0.6, seed = seed
    )
    long <- data.frame(
        site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
        y = as.vector(t(sim$y)),
        det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
    )
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                    det.covs = c("det_cov1", "pos_cov1"))
    # cell-level covariate frame: simulator's occ_cov1 plus a `year` covariate.
    year <- as.numeric(scale(seq_len(N)))
    cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data, year = year)
    y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

    fit <- suppressWarnings(tobs(
        formula = ~ occ_cov1 + year + bym2(graph = adj), data = cell_dat,
        family = occu_cover("lognormal"),
        detection = ~ det_cov1, positive = ~ pos_cov1,
        y = od$y, y_pos = y_pos, visits = od$det.covs,
        method = "nested_laplace",
        control = list(verbose = FALSE, max.iter = 300L,
                       engine = "joint_coupled",
                       sigma.grid = c(0.5, 1.0),
                       alpha.grid = c(0, 0.5, 1.0))
    ))
    list(fit = fit, cell_dat = cell_dat, adj = adj, N = N, f = sim$truth$f)
}

test_that("the joint_coupled occu_cover fit carries a sampleable joint_fit", {
    skip_on_cran()
    f <- .ocp_build_fit()
    expect_s3_class(f$fit$joint_fit, "tulpa_nested_laplace_joint")
    expect_false(is.null(f$fit$joint_fit$Q_csc_p_per_grid))   # store_Q forced on
    jf <- f$fit$joint_fit
    layout <- jf$arm_layout
    idx <- layout$beta_start[1L] + seq_len(layout$p[1L])
    D <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = 50L)
    expect_equal(dim(D), c(50L, layout$p[1L]))
})

test_that("single-quantity predictions: cell + mean + CI + draw matrices", {
    skip_on_cran()
    f <- .ocp_build_fit()
    for (ty in c("occurrence", "cover_cond", "cover_exp")) {
        pr <- predict(f$fit, newdata = f$cell_dat, type = ty, nsim = 400L)
        expect_s3_class(pr, "tobs_prediction")
        expect_true(all(c("cell", "mean", "sd", "lwr", "upr") %in% names(pr)))
        expect_equal(nrow(pr), f$N)
        expect_true(all(pr$lwr <= pr$upr))
        dm <- attr(pr, "draws")[[ty]]
        expect_equal(dim(dm), c(f$N, 400L))
        expect_true(all(pr$mean >= apply(dm, 1L, min) - 1e-9 &
                        pr$mean <= apply(dm, 1L, max) + 1e-9))
    }
    pr <- predict(f$fit, newdata = f$cell_dat, type = "occurrence", nsim = 400L)
    expect_true(all(pr$mean >= 0 & pr$mean <= 1))           # p is a probability
    pr <- predict(f$fit, newdata = f$cell_dat, type = "cover_exp", nsim = 400L)
    expect_true(all(pr$mean >= 0))                          # expected cover >= 0
})

test_that("change decomposition identity holds per draw and in the summary", {
    skip_on_cran()
    f <- .ocp_build_fit()
    pr <- predict(f$fit, newdata = f$cell_dat, type = "change",
                  times = c(0, 1), time_col = "year", nsim = 600L)
    dr <- attr(pr, "draws")
    expect_equal(dr$delta_cover_from_occ + dr$delta_cover_from_ab,
                 dr$delta_cover_exp, tolerance = 1e-6)
    expect_equal(dr$delta_cover_exp, dr$cover_exp_T2 - dr$cover_exp_T1,
                 tolerance = 1e-6)
    expect_equal(pr$delta_cover_from_occ + pr$delta_cover_from_ab,
                 pr$delta_cover_exp, tolerance = 1e-6)
})

test_that("type = change emits the exact column contract keyed by cell", {
    skip_on_cran()
    f <- .ocp_build_fit()
    pr <- predict(f$fit, newdata = f$cell_dat, type = "change",
                  times = c(0, 1), time_col = "year", nsim = 200L)
    point_cols <- c("p_T1", "p_T2", "delta_p",
                    "cover_cond_T1", "cover_cond_T2", "delta_cover_cond",
                    "cover_exp_T1", "cover_exp_T2", "delta_cover_exp",
                    "delta_cover_from_occ", "delta_cover_from_ab")
    delta_ci <- as.vector(outer(
        c("delta_p", "delta_cover_cond", "delta_cover_exp",
          "delta_cover_from_occ", "delta_cover_from_ab"),
        c(".lwr", ".upr"), paste0))
    expect_true(all(c("cell", point_cols, delta_ci) %in% names(pr)))
    expect_equal(nrow(pr), f$N)
    expect_s3_class(pr, "tobs_prediction")
    # year drives occupancy -> delta_p is non-trivial across cells.
    expect_gt(stats::sd(pr$delta_p), 0)
})

test_that("in-sample occurrence tracks the plug-in predictor", {
    skip_on_cran()
    f <- .ocp_build_fit()
    fit <- f$fit
    pr <- predict(fit, type = "occurrence", nsim = 2000L)   # defaults to training
    X <- model.matrix(~ occ_cov1 + year, f$cell_dat)
    b <- fit$means[c("psi_(Intercept)", "psi_occ_cov1", "psi_year")]
    eta <- as.numeric(X %*% b) + fit$spatial$sigma_mean * fit$spatial_field
    p_plug <- plogis(eta)
    # high correlation pins the field amplitude + betas (a wrong sigma factor
    # would collapse it); the level gap is Jensen on plogis.
    expect_gt(stats::cor(pr$mean, p_plug), 0.95)
    expect_lt(max(abs(pr$mean - p_plug)), 0.2)
})

# Trend (time-varying field) fit: the change map moves cover over time THROUGH
# the spatially-varying trend field, the collaborator's headline use case.
.ocp_build_trend_fit <- function(N = 36L, J = 5L, seed = 303L) {
    adj <- .ocp_chain_adj(N)
    sim <- simulate_occu_cover(
        N = N, J = J, positive = "lognormal", adj = adj,
        sigma = 0.8, alpha = 0.6, trend = TRUE,
        sigma_trend = 0.7, alpha_trend = 0.5, seed = seed
    )
    long <- data.frame(
        site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
        y = as.vector(t(sim$y)),
        det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
    )
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                    det.covs = c("det_cov1", "pos_cov1"))
    cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)   # carries `time`
    y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
    fit <- suppressWarnings(tobs(
        formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
        family = occu_cover("lognormal"),
        detection = ~ det_cov1, positive = ~ pos_cov1,
        y = od$y, y_pos = y_pos, visits = od$det.covs,
        method = "nested_laplace",
        control = list(verbose = FALSE, max.iter = 250L,
                       engine = "joint_coupled",
                       trend = list(weight = "time"),
                       sigma.grid = c(0.5, 1.0),
                       alpha.grid = c(0, 0.5))
    ))
    list(fit = fit, cell_dat = cell_dat, N = N)
}

test_that("trend fit: change map runs through the time-varying field", {
    skip_on_cran()
    skip_if_fast()
    f <- .ocp_build_trend_fit()
    # this fit has >1 coupled field block
    expect_gt(length(f$fit$joint_fit$arm_layout$field_starts), 1L)
    # auto-resolve time_col from the fit's stored trend weight
    pr <- predict(f$fit, newdata = f$cell_dat, type = "change",
                  times = c(-1, 1), nsim = 400L)
    expect_s3_class(pr, "tobs_prediction")
    expect_equal(nrow(pr), f$N)
    # decomposition identity still exact per draw
    dr <- attr(pr, "draws")
    expect_equal(dr$delta_cover_from_occ + dr$delta_cover_from_ab,
                 dr$delta_cover_exp, tolerance = 1e-6)
    # the trend field genuinely moves the change across cells (not all equal)
    expect_gt(stats::sd(pr$delta_cover_exp), 0)
    # single-time prediction also works on the trend fit
    p1 <- predict(f$fit, newdata = f$cell_dat, type = "cover_exp",
                  time_col = "time", nsim = 300L)
    expect_equal(nrow(p1), f$N)
    expect_true(all(p1$mean >= 0))
})

test_that("trend fit errors clearly when time_col is unavailable", {
    skip_on_cran()
    skip_if_fast()
    f <- .ocp_build_trend_fit()
    # drop the stored trend weight and don't pass time_col -> clear error
    fit2 <- f$fit
    fit2$trend_weight <- NULL
    expect_error(
        predict(fit2, newdata = f$cell_dat, type = "occurrence"),
        "time_col"
    )
})

# Calibration: the per-cell change CI must cover the KNOWN delta in occupancy
# induced by moving the time covariate from t1 to t2 through the trend field.
# With the field held at its truth, the true per-cell change is
#   delta_p[i] = plogis(eta_base[i] + sigma_trend t2 f2[i])
#              - plogis(eta_base[i] + sigma_trend t1 f2[i]),
# eta_base[i] = X_occ[i] beta_occ + sigma f[i]. A broken field amplitude
# (sigma / alpha) factor or a plug-in of posterior means into the nonlinear
# plogis would collapse this coverage. Returns truth alongside the fit so the
# delta is computed from the generative parameters, not re-estimated.
.ocp_build_trend_fit_truth <- function(N = 40L, J = 6L, seed = 202L) {
    adj <- .ocp_chain_adj(N)
    sim <- simulate_occu_cover(
        N = N, J = J, positive = "lognormal", adj = adj,
        sigma = 0.8, alpha = 0.6, trend = TRUE,
        sigma_trend = 0.7, alpha_trend = 0.5, seed = seed
    )
    long <- data.frame(
        site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
        y = as.vector(t(sim$y)),
        det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
    )
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                    det.covs = c("det_cov1", "pos_cov1"))
    cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
    y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
    fit <- suppressWarnings(tobs(
        formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
        family = occu_cover("lognormal"),
        detection = ~ det_cov1, positive = ~ pos_cov1,
        y = od$y, y_pos = y_pos, visits = od$det.covs,
        method = "nested_laplace",
        control = list(verbose = FALSE, max.iter = 250L, engine = "joint_coupled",
                       trend = list(weight = "time"),
                       sigma.grid = c(0.5, 1.0), alpha.grid = c(0, 0.5))
    ))
    list(fit = fit, cell_dat = cell_dat, truth = sim$truth, N = N)
}

test_that("change CI covers the known per-cell occupancy change", {
    skip_on_cran()
    skip_if_fast()
    f  <- .ocp_build_trend_fit_truth()
    tr <- f$truth
    t1 <- -1; t2 <- 1

    X_occ    <- stats::model.matrix(~ occ_cov1, f$cell_dat)
    eta_base <- as.vector(X_occ %*% tr$beta_occ) + tr$sigma * tr$f
    p_at     <- function(t) stats::plogis(eta_base + tr$sigma_trend * t * tr$f2)
    delta_true <- p_at(t2) - p_at(t1)

    pr <- predict(f$fit, newdata = f$cell_dat, type = "change",
                  times = c(t1, t2), time_col = "time", nsim = 1500L)

    covered <- pr$delta_p.lwr <= delta_true & delta_true <= pr$delta_p.upr
    # 95% CIs over the N cells of one realization; conservative floor (observed
    # coverage runs ~0.93-1.0, a broken amplitude factor would crater it).
    expect_gt(mean(covered), 0.8)
    # the point estimate tracks the true change in sign and magnitude.
    expect_gt(stats::cor(pr$delta_p, delta_true), 0.5)
})
