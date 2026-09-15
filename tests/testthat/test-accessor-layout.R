# One layout for the stats accessors on every tobs fit (#332): coef() is a flat
# vector named `<arm>_<term>` as vcov() / confint() / summary() name it, tidy()
# splits that name into `arm` and `term`, glance() has one column set, and
# simulate() returns a list of `nsim` responses carrying the "seed" attribute.

test_that("the term splitter reads an arm only where the arm declares the term", {
  arms <- list(f_sp1 = c("(Intercept)", "x"), f_sp1_sp2 = c("(Intercept)", "x"),
               mu = c("(Intercept)", "x"))
  sp <- tulpaObs:::.tobs_split_terms(
    c("f_sp1_x", "f_sp1_sp2_(Intercept)", "mu_x", "mu_log_r", "log_r"), arms)
  expect_identical(sp$arm, c("f_sp1", "f_sp1_sp2", "mu", NA, NA))
  expect_identical(sp$term, c("x", "(Intercept)", "x", "mu_log_r", "log_r"))

  # An arm whose process records no coefficient names matches on its prefix.
  sp0 <- tulpaObs:::.tobs_split_terms(c("psi_a", "p_b"), list(psi = NULL, p = NULL))
  expect_identical(sp0$arm, c("psi", "p"))
})

.al_occu <- function() {
  sim <- simulate_occu(N = 40L, J = 3L, seed = 21L)
  suppressWarnings(tobs(~ occ_cov1, data = sim$data, family = occu(),
                        detection = ~ det_cov1, y = sim$y,
                        control = list(verbose = FALSE, progress = FALSE)))
}

test_that("coef() is flat and coef(arm =) reads one arm", {
  fit <- .al_occu()
  cf <- coef(fit)
  expect_type(cf, "double")
  expect_identical(names(cf), rownames(vcov(fit)))
  expect_identical(names(cf), rownames(confint(fit)))
  expect_identical(names(cf), rownames(summary(fit)))
  expect_identical(coef(fit, arm = "psi"),
                   stats::setNames(unname(cf[c("psi_(Intercept)", "psi_occ_cov1")]),
                                   c("(Intercept)", "occ_cov1")))
  expect_error(coef(fit, arm = "lambda"), "must be one of \"psi\", \"p\"")

  td <- tidy(fit)
  expect_identical(names(td)[1:6], c("arm", "term", "estimate", "std.error",
                                     "conf.low", "conf.high"))
  expect_identical(paste(td$arm, td$term, sep = "_"), names(cf))
})

test_that("simulate() follows the stats::simulate() contract", {
  fit <- .al_occu()
  s <- simulate(fit, nsim = 3L, seed = 9L)
  expect_type(s, "list")
  expect_length(s, 3L)
  expect_named(s, c("sim_1", "sim_2", "sim_3"))
  expect_equal(dim(s[[1L]]), dim(fit$model$y))
  expect_identical(as.integer(attr(s, "seed")), 9L)
  expect_identical(attr(attr(s, "seed"), "kind"), as.list(RNGkind()))

  # One replicate is still a list, and a supplied seed reproduces it.
  s1 <- simulate(fit, nsim = 1L, seed = 9L)
  expect_length(s1, 1L)
  expect_identical(s1[[1L]], simulate(fit, nsim = 1L, seed = 9L)[[1L]])

  # A supplied seed leaves the caller's stream where it was.
  set.seed(4L); before <- stats::runif(1L)
  set.seed(4L); simulate(fit, seed = 123L); after <- stats::runif(1L)
  expect_identical(after, before)

  # Without a seed the attribute is the stream the call started from.
  set.seed(5L); state <- .Random.seed
  expect_identical(attr(simulate(fit), "seed"), state)

  expect_error(simulate(fit, nsim = 0), "single positive integer")
  expect_error(simulate(fit, nsim = 1, bogus = TRUE), "no further arguments")
})

test_that("every family reports the one layout", {
  skip_on_cran()
  skip_if_fast()
  glance_cols <- c("nobs", "df", "logLik", "n_fixed", "n_samples",
                   "n_divergent", "mean_accept", "converged")
  for (nm in names(.SBC_REG_FIXTURES)) {
    fit <- suppressWarnings(suppressMessages(.SBC_REG_FIXTURES[[nm]]()))
    cf <- coef(fit)
    expect_type(cf, "double")
    expect_identical(names(cf), rownames(vcov(fit)), info = nm)
    expect_identical(names(cf), rownames(summary(fit)), info = nm)
    expect_identical(names(summary(fit))[1:4],
                     c("estimate", "std.error", colnames(confint(fit))), info = nm)

    td <- tidy(fit)
    expect_identical(names(td)[1:6], c("arm", "term", "estimate", "std.error",
                                       "conf.low", "conf.high"), info = nm)
    on_arm <- !is.na(td$arm)
    expect_identical(paste(td$arm, td$term, sep = "_")[on_arm], names(cf)[on_arm],
                     info = nm)
    expect_identical(td$term[!on_arm], names(cf)[!on_arm], info = nm)
    for (a in unique(td$arm[on_arm])) {
      expect_equal(unname(coef(fit, arm = a)), td$estimate[on_arm & td$arm == a],
                   info = paste(nm, a))
    }

    expect_identical(names(glance(fit))[seq_along(glance_cols)], glance_cols,
                     info = nm)

    s <- simulate(fit, nsim = 2L, seed = 3L)
    expect_type(s, "list")
    expect_named(s, c("sim_1", "sim_2"), info = nm)
    expect_false(is.null(attr(s, "seed")), info = nm)
    expect_identical(class(s[[1L]]), class(s[[2L]]), info = nm)
  }
})

test_that("simulate() on the families that refused before draws their response", {
  skip_on_cran()
  skip_if_fast()
  # cover("lognormal"): absent plots are exactly zero, present plots carry a
  # positive cover whose log scatters like the fitted arm.
  cv <- suppressWarnings(.SBC_REG_FIXTURES$cover())
  y <- simulate(cv, nsim = 1L, seed = 2L)[[1L]]
  expect_length(y, length(cv$encoding$y))
  expect_true(all(y >= 0))
  expect_lt(abs(stats::sd(log(y[y > 0])) - stats::sd(log(cv$encoding$y[cv$encoding$y > 0]))),
            0.25)
  expect_gt(mean(y > 0), 0.1)
  expect_lt(abs(mean(y > 0) - mean(cv$encoding$y > 0)), 0.2)

  # dyn_occu(): a detection array in y's layout, unsurveyed cells kept.
  dy <- suppressWarnings(.SBC_REG_FIXTURES$dyn_occu())
  yd <- simulate(dy, nsim = 1L, seed = 2L)[[1L]]
  expect_identical(dim(yd), dim(dy$model$y))
  expect_identical(is.na(yd), is.na(dy$model$y))
  expect_true(all(yd[!is.na(yd)] %in% c(0L, 1L)))

  # int_occu(): one detection matrix per source, named as `y` was.
  io <- suppressWarnings(.SBC_REG_FIXTURES$int_occu())
  yi <- simulate(io, nsim = 1L, seed = 2L)[[1L]]
  expect_identical(names(yi), names(io$model$y_sources))
  expect_identical(lapply(yi, dim), lapply(io$model$y_sources, dim))

  # occu_categorical(): codes 0..K, as the response was given.
  oc <- suppressWarnings(.SBC_REG_FIXTURES$occu_categorical())
  yc <- simulate(oc, nsim = 1L, seed = 2L)[[1L]]
  expect_true(all(yc %in% 0:oc$K))

  # occu_multiscale_cover() and t_occu(): the response grid, cover only where
  # a visit detected.
  om <- suppressWarnings(.SBC_REG_FIXTURES$occu_multiscale_cover())
  ym <- simulate(om, nsim = 1L, seed = 2L)[[1L]]
  expect_identical(dim(ym$y), dim(om$model$y))
  expect_true(all(is.na(ym$y_pos[!is.na(ym$y) & ym$y == 0L])))
  to <- suppressWarnings(.SBC_REG_FIXTURES$t_occu())
  expect_identical(dim(simulate(to, nsim = 1L, seed = 2L)[[1L]]), dim(to$model$y))
})
