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

.chain_adj_cover <- function(n_s) {
    nbr <- lapply(seq_len(n_s),
                  function(s) setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L)))
    adj <- matrix(0L, n_s, n_s)
    for (s in seq_len(n_s)) for (j in nbr[[s]]) adj[s, j] <- 1L
    adj
}

test_that("control$prior.phi shrinks the cover-arm precision toward zero", {
  skip_if_fast()
    sim <- simulate_beta_cover_pp(N = 240, n_s = 25, phi = 30, seed = 27)
    adj <- .chain_adj_cover(nlevels(sim$data$region))
    phi_grid <- c(5, 15, 40, 100, 250)

    fit_one <- function(prior_phi) tobs(
        formula = ~ x + bym2(graph = adj, group_var = "region"),
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
