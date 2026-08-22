.fit_simple <- function(formula = ~ elev, det = ~ 1, n = 50, seed = 42,
                       method = "laplace") {
  set.seed(seed)
  d <- data.frame(elev = rnorm(n))
  psi <- plogis(0.5 + 0.5 * d$elev)
  z <- rbinom(n, 1, psi)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  fit <- tobs(formula = formula, data = d, family = occu(),
              detection = det, y = y, method = method,
              control = list(verbose = FALSE))
  list(fit = fit, y = y, d = d, n = n)
}

test_that("S3 methods work on single-season fit", {
  res <- .fit_simple()
  fit <- res$fit; y <- res$y; n <- res$n

  cf <- coef(fit)
  expect_type(cf, "list")
  expect_length(cf$psi, 2)
  expect_length(cf$p, 1)

  ci <- confint(fit)
  expect_true(nrow(ci) >= 3)

  V <- vcov(fit)
  expect_true(all(diag(V) > 0))

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")

  expect_equal(nobs(fit), sum(y >= 0, na.rm = TRUE))

  fv <- fitted(fit)
  expect_named(fv, c("psi", "p", "z"))
  expect_length(fv$psi, n)
  expect_true(all(fv$psi >= 0 & fv$psi <= 1))

  r <- residuals(fit)
  expect_named(r, c("occ", "det"))
  expect_length(r$occ, n)
  expect_equal(dim(r$det), c(n, 3))

  y_sim <- simulate(fit, nsim = 1, seed = 1)
  expect_equal(dim(y_sim), dim(y))

  pred <- predict(fit)
  expect_named(pred, c("psi", "p", "z"))

  X0 <- model.matrix(~ elev, data.frame(elev = c(-1, 0, 1)))
  pred_dm <- predict(fit, X.0 = X0)
  expect_equal(nrow(pred_dm), 3)
  expect_true(all(pred_dm$mean >= 0 & pred_dm$mean <= 1))

  td <- tulpa::tidy(fit)
  expect_s3_class(td, "data.frame")

  gl <- tulpa::glance(fit)
  expect_s3_class(gl, "data.frame")

  re <- tulpa::ranef(fit)
  expect_true(is.data.frame(re) || is.list(re))
})

test_that("convergence()/converged() read one record across families (tulpaObs#88)", {
  # occu(): the verdict lives at fit$convergence (the package-wide convention).
  occ <- .fit_simple()$fit
  rec <- convergence(occ)
  expect_type(rec, "list")
  expect_true(is.logical(rec$converged))
  expect_true("n_iter" %in% names(rec))
  expect_identical(converged(occ), isTRUE(occ$convergence$converged))

  # cover(): historically only fit$converged was set, so a consumer reading
  # fit$convergence$converged got NULL/NA (the bug). It must now carry the same
  # unified record, and the accessor must agree with the top-level flag.
  sim <- simulate_cover(N = 200, seed = 7)
  cov_fit <- tobs(formula = ~ x, data = sim$data, family = cover("beta"),
                  y = sim$y)
  expect_s3_class(cov_fit, "cover_fit")
  expect_false(is.null(cov_fit$convergence))
  expect_true(is.logical(cov_fit$convergence$converged))
  expect_identical(cov_fit$convergence$converged, cov_fit$converged)
  expect_true("sla_status" %in% names(cov_fit$convergence))

  # One accessor, same shape and meaning for both families.
  expect_identical(converged(cov_fit), isTRUE(cov_fit$converged))
  expect_named(convergence(cov_fit)[c("converged", "n_iter")],
               c("converged", "n_iter"))
})

test_that("non-NUTS fits report NA sampler diagnostics, NUTS reports numeric", {
  skip_if_fast()
  # tulpaObs#17: a Laplace / nested-Laplace fit ran no HMC trajectory, so the
  # NUTS-only sampler-health fields (acceptance, divergence, tree depth, step
  # size) must be NA rather than the constants 1 / 0 / 0 / 0 -- otherwise a user
  # checking sampler health reads "no sampler ran" as "sampler ran cleanly".
  fit_lap <- .fit_simple(method = "laplace")$fit
  expect_identical(fit_lap$method, "laplace")
  expect_true(all(is.na(fit_lap$accept_prob)))
  expect_true(all(is.na(fit_lap$divergent)))
  expect_true(all(is.na(fit_lap$treedepth)))
  expect_true(is.na(fit_lap$epsilon))

  # Areal nested-Laplace N-mixture build path (abun.R) -- the same rule.
  set.seed(11)
  side <- 5L; ng <- side * side
  co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i != j && abs(co$x[i]-co$x[j]) + abs(co$y[i]-co$y[j]) == 1L) adj[i, j] <- 1L
  phi <- as.numeric(scale(rnorm(ng))) * 0.5; phi <- phi - mean(phi)
  x_ab <- rnorm(ng)
  N <- rpois(ng, exp(log(5) + 0.5 * x_ab + phi))
  yA <- matrix(NA_integer_, ng, 5L)
  for (i in seq_len(ng)) yA[i, ] <- rbinom(5L, N[i], plogis(0.4))
  fit_nl <- tobs(~ abund_cov1 + icar(graph = adj),
                 data = data.frame(abund_cov1 = x_ab), family = abun(),
                 detection = ~ 1, y = yA, method = "nested_laplace",
                 control = list(verbose = FALSE))
  expect_identical(fit_nl$method, "nested_laplace")
  expect_true(all(is.na(fit_nl$accept_prob)))
  expect_true(all(is.na(fit_nl$divergent)))
  expect_true(all(is.na(fit_nl$treedepth)))
  expect_true(is.na(fit_nl$epsilon))

  # A NUTS fit, by contrast, carries real numeric diagnostics.
  fit_nuts <- .fit_simple(method = "nuts")$fit
  expect_identical(fit_nuts$method, "nuts")
  expect_true(any(is.finite(fit_nuts$accept_prob)))
  expect_false(all(is.na(fit_nuts$divergent)))
})

test_that("WAIC works on single-season fit", {
  res <- .fit_simple(formula = ~ elev, n = 30, seed = 42)
  w <- waic(res$fit)
  expect_true(is.finite(w$waic))
  expect_true(is.finite(w$elpd))
  expect_true(w$p_waic >= 0)
})

test_that("PPC works on single-season fit", {
  res <- .fit_simple(formula = ~ 1, n = 30, seed = 42)
  ppc <- ppc(res$fit, n.samples = 50)
  expect_length(ppc$fit.y, 50)
  expect_length(ppc$fit.y.rep, 50)
  expect_true(ppc$bayesian.p >= 0 && ppc$bayesian.p <= 1)
})

test_that("compare_models works", {
  skip_if_fast()
  set.seed(42)
  n <- 30
  d <- data.frame(x = rnorm(n))
  z <- rbinom(n, 1, 0.5)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)

  fit1 <- tobs(~ 1, d, family = occu(), detection = ~ 1, y = y,
               method = "nuts",
               control = list(n.iter = 200, n.warmup = 100, seed = 42, verbose = FALSE))
  fit2 <- tobs(~ x, d, family = occu(), detection = ~ 1, y = y,
               method = "nuts",
               control = list(n.iter = 200, n.warmup = 100, seed = 42, verbose = FALSE))

  comp <- tulpa::compare_models(null = fit1, elev = fit2)
  expect_s3_class(comp, "data.frame")
  expect_equal(nrow(comp), 2)
})

test_that("simulation functions work", {
  sim <- simulate_occu(N = 20, J = 3, seed = 42)
  expect_equal(dim(sim$y), c(20, 3))
  expect_equal(nrow(sim$data), 20)

  sim_ms <- simulate_ms_occu(N = 10, J = 3, n_species = 3, seed = 42)
  expect_equal(dim(sim_ms$y), c(10, 3, 3))

  sim_t <- simulate_dyn_occu(N = 10, J = 3, n_seasons = 4, seed = 42)
  expect_equal(dim(sim_t$y), c(10, 3, 4))
})

test_that("predict(terms=) varies one term and rejects a longer vector", {
  set.seed(4)
  n <- 60
  d <- data.frame(elev = stats::rnorm(n), slope = stats::rnorm(n))
  psi <- plogis(0.4 + 0.6 * d$elev - 0.5 * d$slope)
  z <- rbinom(n, 1, psi)
  y <- matrix(rbinom(n * 3, 1, z * 0.5), n, 3)
  fit <- tobs(~ elev + slope, data = d, family = occu(), detection = ~ 1,
              y = y, control = list(verbose = FALSE))

  # The documented mode: a grid over `elev`, every other design column at its
  # mean, on the occupancy probability scale.
  pr <- predict(fit, terms = "elev", n_points = 12L)
  expect_s3_class(pr, "tobs_prediction")
  expect_equal(nrow(pr), 12L)
  expect_identical(attr(pr, "term"), "elev")
  expect_true(all(pr$estimate >= 0 & pr$estimate <= 1))
  expect_true(all(pr$lower <= pr$estimate & pr$estimate <= pr$upper))

  # A second term is not a grouping variable here: it would be held at its
  # mean, answering a different question under a plausible-looking result. It
  # is an error, not element one of the vector.
  expect_error(predict(fit, terms = c("elev", "slope")), "one term")
  expect_error(predict(fit, terms = character(0)), "one term")
  expect_error(predict(fit, terms = 1), "one term")
  expect_error(predict(fit, terms = NA_character_), "one term")

  # An unknown single term still reports the coefficient names it looked in.
  expect_error(predict(fit, terms = "nope"), "not found")

  # The detection arm takes the same route and the same contract.
  expect_error(predict(fit, terms = c("(Intercept)", "elev"),
                       type = "detection"), "one term")

  # tobs_marginal_effect() calls the same fitter, so it inherits the contract.
  me <- tobs_marginal_effect(fit, "elev", n_points = 8L)
  expect_equal(nrow(me), 8L)
  expect_error(tobs_marginal_effect(fit, c("elev", "slope")), "one term")

  # plot() reads its input shape from this table alone (x / estimate / lower /
  # upper plus the term and process attributes) and returns it invisibly. The
  # middle level travels beside them under its own name.
  expect_named(as.data.frame(pr),
               c("x", "estimate", "lower", "upper", "q50"))
  expect_identical(attr(pr, "process"), "psi")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(plot(pr), "tobs_prediction")
})

test_that("predict(quantiles=) drives the levels AND the column names", {
  res <- .fit_simple()
  fit <- res$fit
  X0 <- model.matrix(~ elev, data.frame(elev = c(-1, 0, 1)))

  # Default layout, with the levels recorded on the table so a downstream
  # consumer can read what the interval columns mean.
  pd <- predict(fit, X.0 = X0)
  expect_named(pd, c("mean", "sd", "q2.5", "q50", "q97.5"))
  expect_equal(attr(pd, "quantiles"), c(0.025, 0.5, 0.975))

  # A non-default request names its own levels and fills them with those
  # levels' values, checked against quantile() on the same draws.
  pr <- predict(fit, X.0 = X0, quantiles = c(0.1, 0.5, 0.9))
  expect_named(pr, c("mean", "sd", "q10", "q50", "q90"))
  expect_equal(attr(pr, "quantiles"), c(0.1, 0.5, 0.9))
  psi <- tulpaObs:::.tobs_psi_draws(fit$draws, X0,
                                    fit$model$process_info[[1]]$p)
  expect_equal(pr$q10, unname(apply(psi, 2, stats::quantile, 0.1)))
  expect_equal(pr$q50, unname(apply(psi, 2, stats::quantile, 0.5)))
  expect_equal(pr$q90, unname(apply(psi, 2, stats::quantile, 0.9)))

  # Narrowing the interval moves the endpoints, it does not only relabel them
  # (the failure was a table whose columns stated the default level whatever
  # was asked for, gcol33/tulpaObs#242).
  expect_true(all(pr$q10 > pd$q2.5))
  expect_true(all(pr$q90 < pd$q97.5))

  # A length-2 vector used to put the 10% point in q2.5 and leave q97.5 all NA.
  expect_error(predict(fit, X.0 = X0, quantiles = c(0.1, 0.9)), "length 3")
  expect_error(predict(fit, X.0 = X0, quantiles = c(0.9, 0.5, 0.1)),
               "strictly increasing")
  expect_error(predict(fit, X.0 = X0, quantiles = c(0, 0.5, 1)),
               "strictly inside")
  expect_error(predict(fit, X.0 = X0, quantiles = c(0.025, NA, 0.975)),
               "no NA")

  # The terms route reports the same levels: lower / upper at the requested
  # ends, and the middle level as its own named column rather than dropped.
  tm <- predict(fit, terms = "elev", n_points = 8L,
                quantiles = c(0.1, 0.5, 0.9))
  expect_named(as.data.frame(tm),
               c("x", "estimate", "lower", "upper", "q50"))
  expect_equal(attr(tm, "quantiles"), c(0.1, 0.5, 0.9))
  expect_true(all(tm$lower <= tm$q50 & tm$q50 <= tm$upper))
  wide <- predict(fit, terms = "elev", n_points = 8L)
  expect_true(all(tm$lower > wide$lower))
  expect_true(all(tm$upper < wide$upper))
  expect_error(predict(fit, terms = "elev", quantiles = c(0.1, 0.9)),
               "length 3")
})


test_that("tobs_predict_spatial builds newocc.covs from the fitted formula", {
  # Minimal object carrying what the predictor reads plus the fitted state
  # formula and its data. Identity link and an all-zero field, so the
  # prediction is the linear predictor exactly.
  fit_data <- data.frame(x = c(0, 1, 0, 1), z = c(0, 0, 1, 1))
  d <- cbind(rep(0.5, 3), rep(2, 3), rep(-2, 3))
  colnames(d) <- c("psi_(Intercept)", "psi_x", "psi_z")
  obj <- structure(list(
    draws = d, spatial = list(type = "icar"), spatial_field = rep(0, 3),
    model = list(
      process_info = list(list(name = "psi", p = 3L, link = "identity",
                               coef_names = c("(Intercept)", "x", "z"))),
      formulas = list(occ = ~ x + z), data = fit_data)),
    class = c("tobs_fit", "tulpa_fit"))
  nodes <- cbind(0:2, 0)
  nc <- cbind(c(0, 1), 0)

  right <- data.frame(x = c(1, 0), z = c(0, 1))
  wrong <- right[, c("z", "x")]           # same data, the user's column order

  pr <- tobs_predict_spatial(obj, nc, newocc.covs = right, node.coords = nodes)
  expect_equal(pr$mean, c(0.5 + 2, 0.5 - 2))

  # `newocc.covs` is a data.frame, so its column order is the user's. The
  # covariate has to meet the coefficient of the same NAME (tulpaObs#243).
  expect_equal(tobs_predict_spatial(obj, nc, newocc.covs = wrong,
                                    node.coords = nodes)$mean, pr$mean)

  # The quantile columns are named from the levels, the same spelling the
  # design-matrix predictor uses (tulpaObs#242).
  expect_named(pr, c("mean", "sd", "q2.5", "q50", "q97.5"))
  expect_named(tobs_predict_spatial(obj, nc, newocc.covs = right,
                                    node.coords = nodes,
                                    quantiles = c(0.1, 0.9)),
               c("mean", "sd", "q10", "q90"))

  # Fewer covariate rows than coordinates used to recycle silently.
  expect_error(
    tobs_predict_spatial(obj, cbind(c(0, 0.3, 0.6, 0.9), 0),
                         newocc.covs = right, node.coords = nodes),
    "one covariate row per prediction location")

  # A frame missing a fitted term fails rather than shifting the pairing.
  expect_error(
    tobs_predict_spatial(obj, nc, newocc.covs = data.frame(x = c(1, 0)),
                         node.coords = nodes),
    "not found")

  # newocc.covs = NULL predicts at the fitted covariate means, which the
  # intercept-plus-zeros design expresses only for a centred fit.
  expect_error(tobs_predict_spatial(obj, nc, node.coords = nodes),
               "not centred")
  cobj <- obj
  cobj$model$data <- data.frame(x = c(-1, 1, -1, 1), z = c(-1, -1, 1, 1))
  expect_equal(tobs_predict_spatial(cobj, nc, node.coords = nodes)$mean,
               rep(0.5, 2))

  # A factor reaches model.matrix() on the fit's own levels, so a newdata frame
  # holding a subset of them still fills the right columns. The raw columns as
  # a design matrix could not represent this at all.
  fdat <- data.frame(g = factor(c("a", "b", "c", "a")))
  df <- cbind(rep(0.1, 3), rep(1, 3), rep(-1, 3))
  colnames(df) <- c("psi_(Intercept)", "psi_gb", "psi_gc")
  fobj <- structure(list(
    draws = df, spatial = list(type = "icar"), spatial_field = rep(0, 3),
    model = list(
      process_info = list(list(name = "psi", p = 3L, link = "identity",
                               coef_names = c("(Intercept)", "gb", "gc"))),
      formulas = list(occ = ~ g), data = fdat)),
    class = c("tobs_fit", "tulpa_fit"))
  expect_equal(
    tobs_predict_spatial(fobj, nc, node.coords = nodes,
                         newocc.covs = data.frame(g = factor(c("c", "a"))))$mean,
    c(0.1 - 1, 0.1))

  # A poly() term re-evaluates on the basis the fit used, so predicting at
  # three of the fitted rows reproduces those rows of the fitted design.
  # Rebuilding the basis from those rows alone is a different answer.
  pdat <- data.frame(x = c(0, 1, 2, 3))
  X_fit <- model.matrix(~ poly(x, 2), pdat)
  beta <- c(0.1, 1, -1)
  dp <- matrix(rep(beta, each = 3L), nrow = 3L)
  colnames(dp) <- paste0("psi_", colnames(X_fit))
  pobj <- structure(list(
    draws = dp, spatial = list(type = "icar"), spatial_field = rep(0, 3),
    model = list(
      process_info = list(list(name = "psi", p = 3L, link = "identity",
                               coef_names = colnames(X_fit))),
      formulas = list(occ = ~ poly(x, 2)), data = pdat)),
    class = c("tobs_fit", "tulpa_fit"))
  rows <- c(1L, 2L, 4L)
  pp <- tobs_predict_spatial(pobj, cbind(c(0, 1, 2), 0), node.coords = nodes,
                             newocc.covs = pdat[rows, , drop = FALSE])
  expect_equal(pp$mean, as.numeric(X_fit[rows, ] %*% beta))
  naive <- model.matrix(~ poly(x, 2), pdat[rows, , drop = FALSE])
  expect_false(isTRUE(all.equal(pp$mean, as.numeric(naive %*% beta))))
})


test_that("nobs() resolves a per-family handler and refuses an unknown type", {
  bare <- function(mt) structure(list(model = list(model_type = mt)),
                                 class = c("tobs_fit", "tulpa_fit"))

  # An unregistered model type used to return NA, which AIC() / BIC() then
  # carried with no indication why (gcol33/tulpaObs#245).
  expect_error(nobs(bare("not_a_family")), "no observation count registered")
  expect_error(nobs(bare("not_a_family")), "not_a_family")

  # Every alias points at a model type that does carry a handler, so the table
  # cannot drift into naming a target that was renamed or removed.
  for (target in unique(unname(tulpaObs:::.TOBS_NOBS_ALIAS))) {
    expect_type(tulpaObs:::.tobs_s3_handler("nobs", target), "closure")
  }

  # Every model type the package builds resolves one, directly or by alias.
  built <- c("single", "dynamic", "integrated", "nmix", "removal", "fp_occu",
             "count", "distance", "dyn_abun", "ms_nmix", "ms_distance",
             "ms_occu_cover", "ms_occu_cover_spatial",
             "occu_multiscale_cover", "ms_occu", "ms_dyn_occu", "ms_count",
             "ms_int_occu", "dyn_int_occu", "occu_cover", "occu_multi",
             "occu_ttd", "royle_nichols", "double_observer", "gdistremoval",
             "distsamp_open", "t_occu")
  for (mt in built) {
    expect_type(
      tulpaObs:::.tobs_s3_handler("nobs", mt, tulpaObs:::.TOBS_NOBS_ALIAS),
      "closure")
  }
})


test_that("predict(terms=) on the N-mixture route holds the same contract", {
  sim <- simulate_abun(N = 40, J = 3, n_abund_covs = 2, n_det_covs = 1,
                       seed = 3)
  fit <- tobs(~ abund_cov1 + abund_cov2, data = sim$data, family = abun(),
              detection = ~ 1, y = sim$y, control = list(verbose = FALSE))

  pr <- predict(fit, terms = "abund_cov1", n_points = 10L)
  expect_s3_class(pr, "tobs_prediction")
  expect_equal(nrow(pr), 10L)
  expect_identical(attr(pr, "term"), "abund_cov1")
  expect_true(all(pr$estimate > 0))   # abundance intensity, not a probability

  expect_error(predict(fit, terms = c("abund_cov1", "abund_cov2")), "one term")
})

test_that("families whose predictor has no terms= argument say so", {
  # `terms` is read by the classic occupancy route and the N-mixture route
  # above. Every other family returns from predict() before reaching it, so a
  # supplied `terms` used to come back as an in-sample fitted() vector that had
  # quietly ignored it -- the same half-answer as reading terms[1] of a vector.
  # Each is told which argument its own predictor varies covariates through.
  sim_rn <- simulate_royle_nichols(N = 60L, J = 4L, beta_lambda = c(0.3, 0.4),
                                   beta_r = -0.8, seed = 11L)
  fit_rn <- tobs(~ x, detection = ~ 1, data = sim_rn$data, y = sim_rn$y,
                 family = royle_nichols(),
                 control = list(verbose = FALSE, progress = FALSE))
  expect_error(predict(fit_rn, terms = "x"), "newdata")

  sim_fp <- simulate_fp_occu(N = 80L, J = 5L, beta_psi = c(0.2, 0.5),
                             p11 = 0.6, p10 = 0.05, b = 0.5, seed = 12L)
  fit_fp <- tobs(~ occ_cov1, detection = ~ 1, data = sim_fp$data, y = sim_fp$y,
                 family = fp_occu(),
                 control = list(verbose = FALSE, progress = FALSE))
  expect_error(predict(fit_fp, terms = "occ_cov1"), "X\\.0")

  # The guard is keyed on the model type, so it must not disturb a predict()
  # that passed no `terms` at all -- both still return their in-sample fit.
  expect_length(predict(fit_rn), nrow(sim_rn$y))
  expect_length(predict(fit_fp), nrow(sim_fp$y))

  # Every model type named in the no-terms set is a real dispatch target, so
  # the set cannot drift into naming a family that does honour `terms`.
  expect_true(all(.TOBS_PREDICT_NO_TERMS %in% names(.tobs_family_methods)))
  expect_false(any(c("nmix", "removal") %in% .TOBS_PREDICT_NO_TERMS))
})

test_that("tobs_data long format conversion works", {
  df <- expand.grid(site = 1:5, visit = 1:3)
  df$detected <- rbinom(15, 1, 0.3)
  df$effort <- rnorm(15)
  df$habitat <- rep(c("forest", "grass", "forest", "grass", "forest"), each = 3)

  od <- tobs_data(df, y = "detected", site = "site", visit = "visit",
                  occ.covs = "habitat", det.covs = "effort")
  expect_s3_class(od, "tobs_data")
  expect_equal(nrow(od$y), 5)
  expect_equal(ncol(od$y), 3)
  expect_equal(nrow(od$occ.covs), 5)
})
