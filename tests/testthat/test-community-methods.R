# test-community-methods.R - predict() / residuals() for the community-occupancy
# families and fitted() / predict() / residuals() for jsdm(), which previously
# errored (a whole family cluster with no post-fit surface).

test_that("predict() / residuals() work on ms_occu", {
  skip_on_cran()
  sim <- simulate_ms_occu(N = 50, J = 3, n_species = 5, seed = 1)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:5), method = "laplace",
              control = list(verbose = FALSE))
  psi <- predict(fit)
  expect_equal(dim(psi), c(50L, 5L))
  expect_true(all(psi >= 0 & psi <= 1))
  expect_equal(colnames(psi), fit$model$species_names)
  # newdata recomputes from the per-species coefficients.
  pn <- predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
  expect_equal(dim(pn), c(3L, 5L))
  # detection prediction.
  pd <- predict(fit, type = "detection")
  expect_equal(dim(pd), c(50L, 5L))
  rr <- residuals(fit)
  expect_equal(dim(rr$occ), c(50L, 5L))
  expect_true(all(is.finite(rr$occ)))
})

test_that("predict() / residuals() work on ms_dyn_occu and ms_int_occu", {
  skip_on_cran()
  skip_if_fast()
  sd_ <- simulate_ms_dyn_occu(N = 40, J = 3, n_species = 5, n_seasons = 3,
                              gamma = 0.2, epsilon = 0.1, seed = 2)
  fd <- tobs(~ 1, data = sd_$data, family = ms_dyn_occu(), detection = ~ 1,
             y = sd_$y, species = paste0("sp", 1:5), method = "laplace",
             control = list(verbose = FALSE))
  expect_equal(dim(predict(fd)), c(40L, 5L))
  expect_true(all(is.finite(residuals(fd)$occ)))

  si <- simulate_ms_int_occu(N = 60, J = c(3, 4), n_species = 5, seed = 3)
  fi <- tobs(~ 1, data = si$data, family = ms_int_occu(), detection = ~ 1,
             y = si$y, species = paste0("sp", 1:5), method = "laplace",
             control = list(verbose = FALSE))
  expect_equal(dim(predict(fi)), c(60L, 5L))
  expect_true(all(is.finite(residuals(fi)$occ)))
})

test_that("fitted() / predict() / residuals() work on jsdm", {
  skip_on_cran()
  # Since jsdm() IS the community GLMM with a logit link, so it shares the
  # ms_count() post-fit surface: fitted()/predict() return the per-(site,
  # species) mean on the response scale (a probability here) under `$mu`, and
  # residuals() returns the unit-level series in `$occ`.
  sim <- simulate_ms_occu(N = 40, J = 1, n_species = 5, seed = 4)
  yj  <- apply(sim$y, c(1, 3),
               function(v) as.integer(any(v[!is.na(v)] == 1)))
  fit <- tobs(~ x, data = sim$data, family = jsdm(), y = yj,
              species = paste0("sp", 1:5), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  f <- fitted(fit)
  expect_equal(dim(f$mu), c(40L, 5L))
  expect_true(all(f$mu >= 0 & f$mu <= 1))
  expect_equal(dim(predict(fit)), c(40L, 5L))
  expect_equal(dim(predict(fit, newdata = data.frame(x = c(-1, 0, 1)))),
               c(3L, 5L))
  rr <- residuals(fit)
  expect_equal(dim(rr$occ), c(40L, 5L))
  expect_true(all(is.finite(rr$occ)))
})


# ---------------------------------------------------------------------------
# The community predict() arm registry (gcol33/tulpaObs#256)
# ---------------------------------------------------------------------------

test_that("every registered community arm is wired to the shared predictor", {
  arms <- tulpaObs:::.TOBS_MS_PREDICT_ARMS
  expect_true(length(arms) >= 6L)

  for (mt in names(arms)) {
    # the family reaches the shared handler, and reports its first arm by
    # default -- both derived from this table, not restated per family
    expect_identical(unname(tulpaObs:::.TOBS_MS_PREDICT_ALIAS[[mt]]),
                     "ms_community")
    expect_identical(unname(tulpaObs:::.TOBS_PREDICT_NEWDATA_TYPE[[mt]]),
                     names(arms[[mt]])[[1L]])
    # and `terms` is refused for it rather than answered from the wrong design
    expect_true(mt %in% tulpaObs:::.TOBS_PREDICT_NO_TERMS)

    for (nm in names(arms[[mt]])) {
      a <- arms[[mt]][[nm]]
      if (!is.null(a$product)) {
        expect_true(all(a$product %in% names(arms[[mt]])))
        next
      }
      expect_true(a$inv %in% c("logit", "log", "cover"))
      expect_true(is.character(a$fitted) && nzchar(a$fitted))
      expect_true(is.character(a$coef) && nzchar(a$coef))
      # an arm with no formula slot has to say why it cannot take newdata
      if (is.null(a$formula)) expect_true(nzchar(a$newdata_reason))
    }
  }
  # the abundance / scale arms are log-linked, which is the whole of #256
  expect_identical(arms$ms_nmix$abundance$inv, "log")
  expect_identical(arms$ms_distance$lambda$inv, "log")
  expect_identical(arms$ms_distance$sigma$inv, "log")
})


test_that("predict(newdata = ) on ms_abun uses the log link, not plogis", {
  skip_on_cran()
  sim <- simulate_ms_abun(n_species = 4, N = 40, J = 3, seed = 1)
  fit <- suppressWarnings(
    tobs(~ abund_cov1, data = sim$data, family = ms_abun(mixture = "poisson"),
         detection = ~ 1, y = sim$y, species = paste0("sp", 1:4),
         method = "laplace", control = list(verbose = FALSE, max.iter = 15L)))

  nd <- data.frame(abund_cov1 = c(-1, 0, 1))
  X  <- stats::model.matrix(~ abund_cov1, nd)
  b  <- fit$ms_community$coef_lambda

  pr <- predict(fit, newdata = nd)
  expect_equal(dim(pr), c(3L, 4L))
  expect_equal(colnames(pr), fit$model$species_names)
  expect_equal(unname(pr), unname(exp(X %*% t(b))))
  # the occupancy fallback this replaces reported plogis() of the same eta
  expect_false(isTRUE(all.equal(unname(pr),
                                unname(stats::plogis(X %*% t(b))))))
  # `newdata` is consulted: rows differ, and the in-sample call is fitted()
  expect_false(isTRUE(all.equal(pr[1L, ], pr[3L, ])))
  expect_identical(predict(fit), fitted(fit)$lambda)
  # the detection arm keeps its own link
  expect_true(all(predict(fit, newdata = nd, type = "detection") > 0 &
                    predict(fit, newdata = nd, type = "detection") < 1))
  # a type this family does not carry is named, not silently substituted
  expect_error(predict(fit, newdata = nd, type = "occupancy"),
               "not a response")
})


test_that("predict(newdata = ) on ms_distance reports lambda and sigma", {
  skip_on_cran()
  cutp <- c(0, 25, 50, 75, 100)
  sim <- simulate_ms_distance(n_species = 3, N = 40, cutpoints = cutp,
                              transect = "line", key = "halfnorm", seed = 0L)
  fit <- suppressWarnings(
    tobs(~ abund_cov1, data = sim$data,
         family = ms_distance(key = "halfnorm", transect = "line",
                              cutpoints = cutp),
         detection = ~ 1, y = sim$y, species = paste0("sp", 1:3),
         method = "laplace", control = list(verbose = FALSE, max.iter = 15L)))

  nd <- data.frame(abund_cov1 = c(-0.5, 0.5))
  X  <- stats::model.matrix(~ abund_cov1, nd)
  expect_equal(unname(predict(fit, newdata = nd)),
               unname(exp(X %*% t(fit$ms_community$coef_lambda))))
  sg <- predict(fit, newdata = nd, type = "sigma")
  expect_true(all(sg > 0))
  expect_equal(dim(sg), c(2L, 3L))
  expect_identical(predict(fit), fitted(fit)$lambda)
  # "abundance" resolves onto this family's own name for that arm
  expect_equal(predict(fit, newdata = nd, type = "abundance"),
               predict(fit, newdata = nd, type = "lambda"))
})


test_that("predict(newdata = ) on ms_occu_cover reports all three arms", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_occu_cover(n_species = 3, N = 30, J = 3,
                                positive = "lognormal", seed = 1)
  fit <- suppressWarnings(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
         species = sim$species, method = "laplace",
         control = list(verbose = FALSE, max.iter = 15L)))

  fv <- fitted(fit)
  expect_identical(predict(fit), fv$psi)
  expect_identical(predict(fit, type = "cover_cond"), fv$cover)
  expect_identical(predict(fit, type = "cover"), fv$cover)
  # the unconditional cover is the product the spatial twin reports under the
  # same name
  expect_equal(predict(fit, type = "cover_exp"), fv$psi * fv$cover)

  # at the training rows the newdata predictor reproduces fitted()
  nd <- data.frame(occ_cov1 = sim$data$occ_cov1[1:5])
  expect_equal(unname(predict(fit, newdata = nd)), unname(fv$psi[1:5, ]))
  ce <- predict(fit, newdata = nd, type = "cover_exp")
  expect_equal(dim(ce), c(5L, 3L))
  expect_true(all(ce > 0))
  expect_equal(unname(ce),
               unname(predict(fit, newdata = nd) *
                        predict(fit, newdata = nd, type = "cover_cond")))
})


test_that("a community fit carrying a latent field refuses newdata", {
  skip_on_cran()
  sim <- simulate_ms_occu(N = 50, J = 3, n_species = 5, seed = 1)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:5), method = "laplace",
              control = list(verbose = FALSE))
  expect_null(tulpaObs:::.tobs_ms_field_bound(fit))

  # a fit whose sites carry an in-sample factor contribution has no value for
  # it at an unseen row, so the predictor declines rather than dropping it
  fit$model$occu_factor_offset <- matrix(0, 50L, 5L)
  expect_identical(tulpaObs:::.tobs_ms_field_bound(fit), "occu_factor_offset")
  expect_error(predict(fit, newdata = data.frame(x = 0)),
               "tied to the in-sample sites")
})
