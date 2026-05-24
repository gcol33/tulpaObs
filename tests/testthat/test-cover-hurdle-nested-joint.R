# Joint nested-Laplace cover_hurdle smoke test (Phase 1c first cut).

simulate_joint_lognormal_cover <- function(N = 250, n_s = 30,
                                           sigma = 0.6, rho = 0.7,
                                           alpha = 1.0,
                                           beta_occ = c(-0.3, 0.7),
                                           beta_pos = c(0.4, -0.5),
                                           sd_pos = 0.4, seed = 11) {
    set.seed(seed)
    spatial_idx <- sample.int(n_s, N, replace = TRUE)
    phi   <- rnorm(n_s, 0, 1)
    theta <- rnorm(n_s, 0, 1)
    w_s   <- sigma * (sqrt(rho) * phi + sqrt(1 - rho) * theta)

    x <- rnorm(N)
    eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
    occur   <- rbinom(N, 1, plogis(eta_occ))

    is_pos  <- occur == 1L
    eta_pos <- beta_pos[1] + beta_pos[2] * x + alpha * w_s[spatial_idx]
    log_y   <- rnorm(N, eta_pos, sd_pos)
    y       <- ifelse(is_pos, exp(log_y), 0)
    y       <- pmin(y, 1 - 1e-6)

    list(
        data = data.frame(x = x, region = factor(spatial_idx)),
        y    = y,
        truth = list(beta_occ = beta_occ, beta_pos = beta_pos,
                     sigma = sigma, rho = rho, alpha = alpha,
                     sd_pos = sd_pos)
    )
}

test_that("cover() with engine='nested_laplace' returns a cover_fit shape", {
    sim <- simulate_joint_lognormal_cover(N = 200, n_s = 25, seed = 13)

    n_s <- nlevels(sim$data$region)
    nbr <- lapply(seq_len(n_s),
                  function(s) setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L)))
    adj <- matrix(0L, n_s, n_s)
    for (s in seq_len(n_s)) for (j in nbr[[s]]) adj[s, j] <- 1L

    fit <- tobs(
        formula  = ~ x + bym2(graph = adj, group_var = "region"),
        data     = sim$data,
        family   = cover("lognormal"),
        y        = sim$y,
        method   = "nested_laplace",
        control  = list(
            sigma_grid     = c(0.4, 0.8),
            rho_grid       = c(0.5, 0.9),
            sigma_pos_grid = c(0.0, 0.6)
        )
    )

    expect_s3_class(fit, "cover_fit")
    expect_equal(fit$positive, "lognormal")
    expect_true(fit$converged)
    expect_true(all(is.finite(fit$beta_occ)))
    expect_true(all(is.finite(fit$beta_pos)))
    expect_named(fit$hyperpar,
                 c("spatial", "engine", "sigma_pos", "sigma_pos_sd"))
    expect_identical(fit$hyperpar$engine, "nested_laplace")

    # Sanity: posterior-weighted slope estimates land on the right side of zero.
    expect_gt(fit$beta_occ[2], 0)   # true 0.7 > 0
    expect_lt(fit$beta_pos[2], 0)   # true -0.5 < 0
})

# ---- ICAR + CAR_proper backends through cover() -------------------------- #

test_that("cover(engine='nested_laplace') accepts ICAR spatial spec", {
    sim <- simulate_joint_lognormal_cover(N = 200, n_s = 25, seed = 17)

    n_s <- nlevels(sim$data$region)
    nbr <- lapply(seq_len(n_s),
                  function(s) setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L)))
    adj <- matrix(0L, n_s, n_s)
    for (s in seq_len(n_s)) for (j in nbr[[s]]) adj[s, j] <- 1L

    fit <- tobs(
        formula  = ~ x + car(graph = adj, group_var = "region"),
        data     = sim$data,
        family   = cover("lognormal"),
        y        = sim$y,
        method   = "nested_laplace",
        control  = list(
            tau_grid       = c(1.0, 4.0),
            sigma_pos_grid = c(0.0, 0.75)
        )
    )

    expect_s3_class(fit, "cover_fit")
    expect_true(fit$converged)
    expect_true(all(is.finite(fit$beta_occ)))
    expect_true(all(is.finite(fit$beta_pos)))
})

# ---- Beta-positive variant through cover(engine='nested_laplace') ------ #

simulate_joint_beta_cover <- function(N = 250, n_s = 30,
                                      sigma = 0.5, rho = 0.7,
                                      alpha = 1.0, phi = 30,
                                      beta_occ = c(-0.3, 0.7),
                                      beta_pos = c(0.4, -0.5),
                                      seed = 23) {
    set.seed(seed)
    spatial_idx <- sample.int(n_s, N, replace = TRUE)
    phi_field   <- rnorm(n_s, 0, 1)
    theta_field <- rnorm(n_s, 0, 1)
    w_s         <- sigma * (sqrt(rho) * phi_field + sqrt(1 - rho) * theta_field)

    x <- rnorm(N)
    eta_occ <- beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]
    occur   <- rbinom(N, 1, plogis(eta_occ))

    eta_pos <- beta_pos[1] + beta_pos[2] * x + alpha * w_s[spatial_idx]
    mu_pos  <- plogis(eta_pos)
    y       <- numeric(N)
    is_pos  <- occur == 1L
    y[is_pos]  <- rbeta(sum(is_pos),
                       mu_pos[is_pos] * phi,
                       (1 - mu_pos[is_pos]) * phi)
    y[!is_pos] <- 0
    y <- pmin(pmax(y, 0), 1 - 1e-6)

    list(
        data = data.frame(x = x, region = factor(spatial_idx)),
        y    = y,
        truth = list(beta_occ = beta_occ, beta_pos = beta_pos,
                     sigma = sigma, rho = rho, alpha = alpha, phi = phi)
    )
}

test_that("cover('beta', engine='nested_laplace') BYM2 returns cover_fit", {
    sim <- simulate_joint_beta_cover(N = 220, n_s = 25, seed = 27)

    n_s <- nlevels(sim$data$region)
    nbr <- lapply(seq_len(n_s),
                  function(s) setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L)))
    adj <- matrix(0L, n_s, n_s)
    for (s in seq_len(n_s)) for (j in nbr[[s]]) adj[s, j] <- 1L

    fit <- tobs(
        formula  = ~ x + bym2(graph = adj, group_var = "region"),
        data     = sim$data,
        family   = cover("beta"),
        y        = sim$y,
        method   = "nested_laplace",
        control  = list(
            sigma_grid     = c(0.4, 0.8),
            rho_grid       = c(0.5, 0.9),
            sigma_pos_grid = c(0.0, 0.6)
        )
    )

    expect_s3_class(fit, "cover_fit")
    expect_equal(fit$positive, "beta")
    expect_true(fit$converged)
    expect_true(all(is.finite(fit$beta_occ)))
    expect_true(all(is.finite(fit$beta_pos)))
    expect_true(is.finite(fit$phi_pos) && fit$phi_pos > 0)
    expect_named(fit$hyperpar,
                 c("spatial", "engine", "phi_pos", "phi_pos_sd"))

    # Posterior-weighted slope estimates land on the right side of zero.
    expect_gt(fit$beta_occ[2], 0)
    expect_lt(fit$beta_pos[2], 0)
})

test_that("cover('beta', engine='nested_laplace') accepts ICAR spatial spec", {
    sim <- simulate_joint_beta_cover(N = 220, n_s = 25, seed = 31)

    n_s <- nlevels(sim$data$region)
    nbr <- lapply(seq_len(n_s),
                  function(s) setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L)))
    adj <- matrix(0L, n_s, n_s)
    for (s in seq_len(n_s)) for (j in nbr[[s]]) adj[s, j] <- 1L

    fit <- tobs(
        formula  = ~ x + car(graph = adj, group_var = "region"),
        data     = sim$data,
        family   = cover("beta"),
        y        = sim$y,
        method   = "nested_laplace",
        control  = list(
            tau_grid       = c(1.0, 4.0),
            sigma_pos_grid = c(0.0, 0.75)
        )
    )

    expect_s3_class(fit, "cover_fit")
    expect_equal(fit$positive, "beta")
    expect_true(fit$converged)
    expect_true(is.finite(fit$phi_pos) && fit$phi_pos > 0)
    expect_true(all(is.finite(fit$beta_pos)))
})

test_that("cover(engine='nested_laplace') accepts CAR_proper spatial spec", {
    sim <- simulate_joint_lognormal_cover(N = 200, n_s = 25, seed = 19)

    n_s <- nlevels(sim$data$region)
    nbr <- lapply(seq_len(n_s),
                  function(s) setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L)))
    adj <- matrix(0L, n_s, n_s)
    for (s in seq_len(n_s)) for (j in nbr[[s]]) adj[s, j] <- 1L

    fit <- tobs(
        formula  = ~ x + car_proper(graph = adj, group_var = "region"),
        data     = sim$data,
        family   = cover("lognormal"),
        y        = sim$y,
        method   = "nested_laplace",
        control  = list(
            tau_grid       = c(1.0, 4.0),
            rho_car_grid   = c(0.7, 0.95),
            sigma_pos_grid = c(0.0, 0.75)
        )
    )

    expect_s3_class(fit, "cover_fit")
    expect_true(fit$converged)
    expect_true(all(is.finite(fit$beta_occ)))
    expect_true(all(is.finite(fit$beta_pos)))
})
