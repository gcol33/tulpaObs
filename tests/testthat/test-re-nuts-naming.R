# A random-effect term reports the same parameters on either engine (#381):
# sigma_* / cor_* on the natural scale and re_*[k] group effects, with the
# group effects kept out of summary() / coef().

re_slope_fixture <- function() {
  set.seed(5); G <- 12; per <- 20; N <- G * per; J <- 4
  g <- rep(seq_len(G), each = per); w <- rnorm(N); x <- rnorm(N)
  S <- matrix(c(0.64, 0.4, 0.4, 1), 2)
  B <- matrix(rnorm(G * 2), G) %*% chol(S)
  z <- rbinom(N, 1, plogis(0.2 + 0.5 * x + B[g, 1] + B[g, 2] * w))
  y <- matrix(rbinom(N * J, 1, z * 0.6), N, J)
  list(d = data.frame(x = x, w = w, g = factor(g)), y = y)
}

test_that("the NUTS sampler coordinates are reported as Laplace names them", {
  design <- list(list(n_coefs = 2L, n_groups = 3L, group_label = "g1",
                      coef_names = c("(Intercept)", "w"), correlated = TRUE,
                      levels = c("a", "b", "c")))
  nc <- 2L
  n_raw <- length(.tobs_re_nuts_param_names(design))
  draws <- matrix(rnorm(5 * (2 + n_raw), sd = 0.3), 5)
  nat <- .tobs_re_nuts_natural(draws, design, 2L)
  expect_identical(ncol(nat), n_raw)
  expect_equal(unname(nat[, 1:nc]), exp(draws[, 2 + 1:nc]))
  cr <- nat[, grep("^cor_", colnames(nat))]
  expect_true(all(abs(cr) <= 1))
  expect_false(any(grepl("^(log_sigma|chol|z)_", colnames(nat))))
})

test_that("occu() (w | g) on NUTS reports sigma / cor like Laplace", {
  skip_if_fast()
  skip_on_cran()
  fx <- re_slope_fixture()
  ctl <- list(progress = FALSE, verbose = FALSE)
  fl <- suppressWarnings(tobs(~ x + (w | g), data = fx$d, y = fx$y,
                              detection = ~ 1, family = occu(),
                              method = "laplace", control = ctl))
  fn <- suppressWarnings(tobs(~ x + (w | g), data = fx$d, y = fx$y,
                              detection = ~ 1, family = occu(), method = "nuts",
                              control = c(ctl, list(n.iter = 150L,
                                                    n.warmup = 150L))))
  expect_identical(names(fn$means), names(fl$means))
  expect_identical(rownames(summary(fn)), rownames(summary(fl)))
  expect_identical(names(coef(fn)), names(coef(fl)))
  expect_gt(fn$means[["sigma_g1_w"]], 0)
  expect_lte(abs(fn$means[["cor_g1_(Intercept)_w"]]), 1)
  expect_identical(nrow(ranef(fn)), nrow(ranef(fl)))
})
