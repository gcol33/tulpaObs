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

    spatial <- tulpa::spatial_bym2(adj, level = "group", group_var = "region")

    fit <- tobs(
        formula  = ~ x,
        data     = sim$data,
        family   = cover("lognormal"),
        y        = sim$y,
        spatial  = spatial,
        engine   = "nested_laplace",
        control  = list(
            sigma_grid = c(0.4, 0.8),
            rho_grid   = c(0.5, 0.9),
            alpha_grid = c(0.0, 1.0)
        )
    )

    expect_s3_class(fit, "cover_fit")
    expect_equal(fit$positive, "lognormal")
    expect_true(fit$converged)
    expect_true(all(is.finite(fit$beta_occ)))
    expect_true(all(is.finite(fit$beta_pos)))
    expect_named(fit$hyperpar, c("spatial", "sigma_pos", "engine"))
    expect_identical(fit$hyperpar$engine, "nested_laplace")

    # Sanity: posterior-weighted slope estimates land on the right side of zero.
    expect_gt(fit$beta_occ[2], 0)   # true 0.7 > 0
    expect_lt(fit$beta_pos[2], 0)   # true -0.5 < 0
})
