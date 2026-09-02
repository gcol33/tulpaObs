# cover() spatial: control$prior.phi re-weights the cover-arm dispersion grid.
#
# The beta-precision phi on the positive arm rides its own outer-grid axis
# (control$phi.grid). control$prior.phi forwards a regularizing hyperprior to
# tulpa's prior_phi ( + tulpaObs side), so the phi axis is weighted by the
# chosen density instead of an implicit flat prior. A sharp half-normal on the
# precision pulls the posterior phi_pos below the flat fit.

simulate_beta_cover_pp <- function(N = 240, n_s = 25, sigma = 0.5, rho = 0.7,
                                   phi = 30, beta_occ = c(-0.3, 0.7),
                                   beta_pos = c(0.4, -0.5), seed = 23) {
    set.seed(seed)
    spatial_idx <- sample.int(n_s, N, replace = TRUE)
    w_s <- sigma * (sqrt(rho) * rnorm(n_s) + sqrt(1 - rho) * rnorm(n_s))
    x <- rnorm(N)
    occur <- rbinom(N, 1, plogis(beta_occ[1] + beta_occ[2] * x + w_s[spatial_idx]))
    mu_pos <- plogis(beta_pos[1] + beta_pos[2] * x + w_s[spatial_idx])
    y <- numeric(N); is_pos <- occur == 1L
    y[is_pos] <- rbeta(sum(is_pos), mu_pos[is_pos] * phi, (1 - mu_pos[is_pos]) * phi)
    y <- pmin(pmax(y, 0), 1 - 1e-6)
    list(data = data.frame(x = x, region = factor(spatial_idx)), y = y)
}

test_that("control$prior.phi shrinks the cover-arm precision toward zero", {
  skip_if_fast()
    sim <- simulate_beta_cover_pp(N = 240, n_s = 25, phi = 30, seed = 27)
    adj <- chain_adj(nlevels(sim$data$region))
    phi_grid <- c(5, 15, 40, 100, 250)

    # The simulator puts w_s in the occurrence logit AND in mu_pos, so the arms
    # share the field and the formula has to say so. Uncoupled, w_s stays in the
    # cover arm's residual and phi_pos pins against the bottom of the axis
    # (8.9 against a truth of 30) -- and a prior that shrinks TOWARD zero then
    # has nowhere left to pull, so the comparison below turns on numerical noise
    # and can land either way. Coupled, phi_pos recovers to 23.1 and the
    # half-normal moves it to 22.8, which is the effect this block is about.
    fit_one <- function(prior_phi) tobs(
        formula = ~ x + bym2(graph = adj, group_var = "region") +
            share(spatial(), alpha = grid(c(0.5, 1.0, 1.5))),
        data = sim$data, family = cover("beta"), y = sim$y,
        method = "nested_laplace",
        control = list(sigma.grid = c(0.4, 0.8), rho.grid = c(0.5, 0.9),
                       phi.grid = phi_grid,
                       prior.phi = prior_phi))

    flat <- fit_one(NULL)
    # A half-normal on the precision down-weights the large-phi cells.
    hn   <- fit_one(list("half_normal", 20))

    expect_s3_class(flat, "cover_fit")
    expect_true(flat$converged && hn$converged)
    expect_true(is.finite(flat$phi_pos) && is.finite(hn$phi_pos))
    expect_lt(hn$phi_pos, flat$phi_pos)
})
