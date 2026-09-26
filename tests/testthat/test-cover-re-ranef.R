# cover() + a grouped random intercept on the nested_laplace joint path reports
# its per-group BLUPs under the factor labels and its SD by name (#373).

cover_re_fixture <- function() {
  chain_adj <- function(n) {
    a <- matrix(0L, n, n)
    for (i in seq_len(n - 1L)) { a[i, i + 1L] <- 1L; a[i + 1L, i] <- 1L }
    a
  }
  set.seed(12); n_s <- 12; ng <- 6; reps <- 6; adj <- chain_adj(n_s)
  combos <- expand.grid(cell = seq_len(n_s), obs = seq_len(ng))
  combos <- combos[rep(seq_len(nrow(combos)), each = reps), ]
  cell <- combos$cell; N <- length(cell); x <- rnorm(n_s)[cell]
  lab <- sprintf("obs%02d", sample.int(ng)); btrue <- rnorm(ng, sd = 0.9)
  b <- btrue[combos$obs]
  occur <- rbinom(N, 1, plogis(-0.2 + 0.6 * x + b))
  mu <- plogis(0.2 - 0.4 * x + b)
  y <- numeric(N); pos <- occur == 1
  y[pos] <- pmin(pmax(rbeta(sum(pos), mu[pos] * 18, (1 - mu[pos]) * 18), 1e-6),
                 1 - 1e-6)
  d <- data.frame(x = x, region = factor(cell), obs = factor(lab[combos$obs]))
  list(d = d, y = y, adj = adj, truth = stats::setNames(btrue, lab))
}

test_that("cover() (1|g) exposes labelled BLUPs and a named sigma_re", {
  skip_if_fast()
  skip_on_cran()
  fx <- cover_re_fixture()
  adj <- fx$adj
  ctl <- list(verbose = FALSE, progress = FALSE, aggregate.occ = FALSE,
              sigma.grid = c(0.5, 1.0), rho.grid = 0.5,
              sigma.re.grid = c(0.5, 1.0), phi.grid = c(8, 18, 40),
              adaptive.grid = FALSE, max.iter = 300L)
  fit <- tobs(~ x + bym2(graph = adj, group_var = "region") + (1 | obs),
              data = fx$d, family = cover("beta"), y = fx$y,
              method = "nested_laplace", control = ctl)
  rf <- ranef(fit)
  expect_setequal(rf$level, levels(fx$d$obs))
  expect_identical(nrow(rf), 6L)
  expect_true(all(is.finite(rf$estimate)) && all(rf$std.error > 0))
  # The BLUPs follow the simulated group effects.
  est <- stats::setNames(rf$estimate, rf$level)[names(fx$truth)]
  expect_gt(stats::cor(est, fx$truth), 0.5)
  expect_named(fit$sigma_re, "sigma_re")
  expect_true(fit$sigma_re >= 0.5 && fit$sigma_re <= 1.0)
})
