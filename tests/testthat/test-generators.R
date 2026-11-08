# One generator per family: simulate_<family>(), simulate(fit) and sbc() all
# draw through `.tobs_gen_<key>(model, theta)`. Each block below takes the
# layout and theta the exported simulator hands the generator, feeds the same
# theta to the fit route (a fit object whose draw is that theta) and to the
# SBC route (the registry's `simulate` slot at that theta), and asserts that
# each route hands the generator that theta and returns exactly what the
# generator draws from the same stream position.

# Every call a route makes to family `key`'s generator inside `expr`: the
# layout, the theta and the RNG state at the call. The generator still runs.
gen_calls <- function(key, expr) {
  name <- paste0(".tobs_gen_", key)
  real <- utils::getFromNamespace(name, "tulpaObs")
  calls <- list()
  fake <- function(model, theta) {
    calls[[length(calls) + 1L]] <<- list(model = model, theta = theta,
                                         rng = get(".Random.seed", globalenv()))
    real(model, theta)
  }
  args <- stats::setNames(list(fake), name)
  do.call(testthat::local_mocked_bindings,
          c(args, list(.package = "tulpaObs", .env = environment())))
  value <- force(expr)
  list(value = value, calls = calls)
}

# The generator's draw for a recorded call, replayed from the stream position
# the route was at.
gen_replay <- function(key, call) {
  assign(".Random.seed", call$rng, envir = globalenv())
  tulpaObs:::.tobs_gen_strip(
    tulpaObs:::.tobs_generator(key)(call$model, call$theta))
}

# A route made exactly one generator call, handed it theta `theta`, and
# returned what the generator draws from that stream position.
expect_route <- function(key, run, theta, response = identity) {
  expect_length(run$calls, 1L)
  call <- run$calls[[1L]]
  expect_equal(call$theta, theta, ignore_attr = TRUE, tolerance = 1e-12)
  expect_identical(response(run$value), gen_replay(key, call))
}

set_rng <- function(seed) {
  set.seed(seed)
  get(".Random.seed", globalenv())
}


test_that("gdistremoval: the three routes draw through one generator", {
  a <- gen_calls("gdistremoval", simulate_gdistremoval(N = 20, seed = 1))
  ca <- a$calls[[1L]]
  rep_a <- gen_replay("gdistremoval", ca)
  expect_identical(unname(a$value$y), unname(rep_a$yDist))
  expect_identical(unname(a$value$y_rem), unname(rep_a$yRem))

  fit <- list(model = ca$model, draws = tulpaObs:::.tobs_theta_row(ca$theta))
  set.seed(5)
  b <- gen_calls("gdistremoval",
                 tulpaObs:::.tobs_simulate_gdistremoval(fit, 1L)[[1L]])
  expect_route("gdistremoval", b, ca$theta)

  spec <- list(model = c(ca$model, list(data = a$value$data)))
  c_ <- gen_calls("gdistremoval",
                  tulpaObs:::.TOBS_SBC_REGISTRY$gdistremoval$simulate(
                    spec, ca$theta, 9L))
  expect_route("gdistremoval", c_, ca$theta, function(d) d$y)
})


test_that("occu_multi: the three routes draw through one generator", {
  a <- gen_calls("occu_multi", simulate_occu_multi(S = 2, N = 15, J = 3, seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("occu_multi", ca))

  fit <- list(model = ca$model, draws = tulpaObs:::.tobs_theta_row(ca$theta))
  set.seed(5)
  b <- gen_calls("occu_multi", tulpaObs:::.tobs_simulate_occu_multi(fit, 1L)[[1L]])
  expect_route("occu_multi", b, ca$theta)

  spec <- list(model = ca$model,
               data_obs = list(cells = a$value$data, site = seq_len(15)))
  c_ <- gen_calls("occu_multi",
                  tulpaObs:::.TOBS_SBC_REGISTRY$occu_multi$simulate(
                    spec, ca$theta, 9L))
  expect_route("occu_multi", c_, ca$theta, function(d) d$y)
})


test_that("occu_categorical: the three routes draw through one generator", {
  a <- gen_calls("occu_categorical", simulate_occu_categorical(N = 40, seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("occu_categorical", ca))

  # A non-finite covariance makes the coefficient draw the mean itself, so the
  # draw that reaches the generator is theta.
  Beta <- ca$theta$Beta
  p <- nrow(Beta); K <- ncol(Beta) + 1L
  fit <- list(encoding = list(X_occ = ca$model$X, K = K),
              beta_occ = ca$theta$beta_occ,
              vcov_occ = matrix(NA_real_, p, p),
              beta_class = Beta,
              vcov_class = matrix(NA_real_, p * (K - 1L), p * (K - 1L)))
  set.seed(5)
  b <- gen_calls("occu_categorical",
                 tulpaObs:::.tobs_simulate_occu_categorical(fit, 1L)[[1L]])
  expect_route("occu_categorical", b, ca$theta)

  dimnames(Beta) <- list(c("(Intercept)", "x"), paste0("k", seq_len(K - 1L)))
  obs_fit <- list(K = K, beta_occ = c(`(Intercept)` = 0, x = 0),
                  beta_class = Beta)
  th <- c(stats::setNames(ca$theta$beta_occ, c("occ_(Intercept)", "occ_x")),
          stats::setNames(as.vector(Beta),
                          as.vector(outer(rownames(Beta), colnames(Beta),
                                          function(r, k) paste0("class_", k, "_", r)))))
  spec <- list(args = list(formula = ~ x), fit_obs = obs_fit,
               data_obs = list(cells = a$value$data, site = seq_len(40)))
  c_ <- gen_calls("occu_categorical",
                  tulpaObs:::.TOBS_SBC_REGISTRY$occu_categorical$simulate(
                    spec, th, 9L))
  expect_route("occu_categorical", c_, ca$theta, function(d) d$y)
})


test_that("t_occu: the three routes draw through one generator", {
  a <- gen_calls("t_occu", simulate_t_occu(N = 20, T.seasons = 3, J = 2,
                                           beta.occ = c(0.2, 0.5), seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("t_occu", ca))

  model <- c(ca$model, list(process_info = list(
    list(p = 2L, coef_names = c("(Intercept)", "x")),
    list(p = 1L, coef_names = "(Intercept)"))))
  fit <- list(model = model, temporal_field = ca$theta$eta,
              draws = tulpaObs:::.tobs_theta_row(
                c(ca$theta$beta_psi, ca$theta$beta_p)))
  set.seed(5)
  b <- gen_calls("t_occu", tulpaObs:::.tobs_simulate_t_occu(fit, 1L)[[1L]])
  expect_route("t_occu", b, ca$theta)

  # The SBC theta carries the AR1 hyperparameters; the year effect it hands
  # the generator is the AR1 draw at the head of the replicate's stream.
  th <- c(`psi_(Intercept)` = 0.2, psi_x = 0.5,
          `p_(Intercept)` = ca$theta$beta_p, log_sigma_ar1 = log(0.7),
          rho_ar1 = 0.6)
  spec <- list(model = model,
               data_obs = list(cells = a$value$data, site = seq_len(20)))
  c_ <- gen_calls("t_occu",
                  tulpaObs:::.TOBS_SBC_REGISTRY$t_occu$simulate(spec, th, 9L))
  set.seed(9L)
  eta <- tulpaObs:::.t_occu_draw_ar1(3L, 0.7, 0.6)
  expect_route("t_occu", c_,
               list(beta_psi = c(0.2, 0.5), beta_p = ca$theta$beta_p, eta = eta),
               function(d) d$y)
})


test_that("dyn_int_occu: the three routes draw through one generator", {
  a <- gen_calls("dyn_int_occu",
                 simulate_dyn_int_occu(N = 20, T.seasons = 3, seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("dyn_int_occu", ca))

  fit <- list(model = ca$model,
              draws = tulpaObs:::.tobs_theta_row(ca$theta$beta))
  set.seed(5)
  b <- gen_calls("dyn_int_occu",
                 tulpaObs:::.tobs_simulate_dyn_int_occu(fit, 1L)[[1L]])
  th <- list(beta = ca$theta$beta)
  expect_route("dyn_int_occu", b, th)

  spec <- list(model = ca$model,
               data_obs = list(cells = a$value$data, site = seq_len(20)))
  c_ <- gen_calls("dyn_int_occu",
                  tulpaObs:::.TOBS_SBC_REGISTRY$dyn_int_occu$simulate(
                    spec, ca$theta$beta, 9L))
  expect_route("dyn_int_occu", c_, th, function(d) d$y)
})


test_that("cover: the three routes draw through one generator", {
  a <- gen_calls("cover", simulate_cover(N = 40, seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("cover", ca))

  # The fit route reads one posterior draw of both arms' predictors; here that
  # draw is the predictor theta implies.
  X <- ca$model$X_occ
  eta_occ <- as.vector(X %*% ca$theta$b_occ) + ca$theta$off_occ
  eta_pos <- as.vector(X %*% ca$theta$b_pos) + ca$theta$off_pos
  testthat::local_mocked_bindings(
    .tobs_cover_eta_draws = function(object, n.draws, rows) {
      list(eta_occ = rbind(eta_occ), eta_pos = rbind(eta_pos),
           disp = ca$theta$disp)
    }, .package = "tulpaObs")
  fit <- list(encoding = list(), positive = "lognormal")
  set.seed(5)
  b <- gen_calls("cover", tulpaObs:::.tobs_simulate_cover(fit, 1L)[[1L]])
  expect_length(b$calls, 1L)
  expect_equal(b$calls[[1L]]$theta,
               list(off_occ = eta_occ, off_pos = eta_pos, disp = ca$theta$disp),
               ignore_attr = TRUE)
  expect_identical(b$value, gen_replay("cover", b$calls[[1L]]))
  set_rng(5)
  expect_identical(b$value, tulpaObs:::.tobs_gen_strip(
    tulpaObs:::.tobs_gen_cover(ca$model, ca$theta)))

  th <- c(`occ_(Intercept)` = -0.5, occ_x = 0.8,
          `pos_(Intercept)` = -1.0, pos_x = 0.3)
  spec <- list(args = list(presence = ~ x, positive = ~ x), sigma_pos = 0.4,
               data_obs = list(cells = a$value$data, site = seq_len(40)))
  c_ <- gen_calls("cover",
                  tulpaObs:::.TOBS_SBC_REGISTRY$cover$simulate(spec, th, 9L))
  expect_route("cover", c_,
               list(b_occ = c(-0.5, 0.8), b_pos = c(-1.0, 0.3), disp = 0.4),
               function(d) d$y)
  set_rng(9L)
  expect_identical(c_$value$y, tulpaObs:::.tobs_gen_strip(
    tulpaObs:::.tobs_gen_cover(ca$model, ca$theta)))
})


test_that("occu_multiscale_cover: the three routes draw through one generator", {
  a <- gen_calls("occu_multiscale_cover",
                 simulate_occu_multiscale_cover(n.cells = 6, plots.per.cell = 2,
                                                visits.per.plot = 2, alpha = 1,
                                                seed = 1))
  ca <- a$calls[[1L]]
  rep_a <- gen_replay("occu_multiscale_cover", ca)
  expect_identical(a$value$y, rep_a$y)
  expect_identical(a$value$y_pos, rep_a$y_pos)

  pinfo <- lapply(1:4, function(k) list(p = 2L))
  model <- c(ca$model, list(process_info = pinfo))
  th <- ca$theta
  par <- c(th$b_psi, th$b_theta, th$b_p_site, th$b_pos_site, th$disp)
  names(par) <- c("psi_a", "psi_b", "theta_a", "theta_b", "p_a", "p_b",
                  "pos_a", "pos_b", "phi_pos")
  fit <- list(model = model, draws = tulpaObs:::.tobs_theta_row(par),
              spatial_field = th$field_psi, means = c(alpha = 1))
  set.seed(5)
  b <- gen_calls("occu_multiscale_cover",
                 tulpaObs:::.tobs_simulate_occu_multiscale_cover(fit, 1L)[[1L]])
  expect_route("occu_multiscale_cover", b, th)

  # The SBC registration is the non-spatial engine's: the field is off.
  names(par)[9L] <- "log_sigma_pos"
  par[9L] <- log(0.35)
  spec <- list(model = model,
               data_obs = list(cells = a$value$data, site = model$plot_cell,
                               n_cells = 6L))
  c_ <- gen_calls("occu_multiscale_cover",
                  tulpaObs:::.TOBS_SBC_REGISTRY$occu_multiscale_cover$simulate(
                    spec, par, 9L))
  th0 <- th
  th0$disp <- exp(log(0.35))
  th0$field_psi <- rep(0, 6L); th0$field_pos <- rep(0, 12L)
  th0$b_p_visit <- par[integer(0)]; th0$b_pos_visit <- par[integer(0)]
  expect_route("occu_multiscale_cover", c_, th0,
               function(d) list(y = d$y, y_pos = d$y_pos))
})


test_that("occu_cover: the three routes draw through one generator", {
  a <- gen_calls("occu_cover", simulate_occu_cover(N = 12, J = 3, seed = 1))
  ca <- a$calls[[1L]]
  grid_a <- tulpaObs:::.occu_cover_gen_shape(ca$model,
                                             gen_replay("occu_cover", ca))
  expect_identical(a$value$y, grid_a$y)
  expect_identical(a$value$y_pos, grid_a$y_pos)

  lt <- ca$model; th <- ca$theta
  model <- list(n_sites = 12L, max_visits = 3L, valid = lt$valid,
                y = a$value$y, y_pos = a$value$y_pos, X_occ = lt$X_occ,
                X_det_site = lt$X_det_site, X_pos_site = lt$X_pos_site,
                X_det_visit = lt$X_det_visit, X_pos_visit = lt$X_pos_visit,
                positive = "lognormal", cover_aggregate = "none")
  testthat::local_mocked_bindings(
    .tobs_occu_cover_components = function(object, n.draws) {
      list(b_occ = rbind(th$b_occ), b_det = rbind(th$b_det),
           b_pos = rbind(th$b_pos), disp = th$disp,
           field_occ = matrix(th$field_occ, ncol = 1L),
           field_pos = matrix(th$field_pos, ncol = 1L),
           off_det = matrix(th$off_det, ncol = 1L),
           off_pos = matrix(th$off_pos, ncol = 1L))
    }, .package = "tulpaObs")
  set.seed(5)
  b <- gen_calls("occu_cover",
                 tulpaObs:::.tobs_simulate_occu_cover(list(model = model), 1L)[[1L]])
  expect_length(b$calls, 1L)
  expect_equal(b$calls[[1L]]$theta, th, ignore_attr = TRUE)
  expect_identical(b$value, tulpaObs:::.occu_cover_gen_shape(
    b$calls[[1L]]$model, gen_replay("occu_cover", b$calls[[1L]])))
  set_rng(5)
  expect_identical(b$value, tulpaObs:::.occu_cover_gen_shape(
    lt, tulpaObs:::.tobs_gen_occu_cover(lt, th)))

  # Field off, so the SBC theta map draws nothing ahead of the generator.
  par <- c(th$b_occ, th$b_det, th$b_pos, disp = th$disp)
  attr(par, "blocks") <- list(occ = 1:2, det = 3:4, pos = 5:6)
  spec <- list(model = model, positive = "lognormal", has_field = FALSE,
               site_cell = seq_len(12L),
               data_obs = list(cells = a$value$data, site = seq_len(12L),
                               visits = NULL))
  c_ <- gen_calls("occu_cover",
                  tulpaObs:::.TOBS_SBC_REGISTRY$occu_cover$simulate(spec, par, 9L))
  expect_length(c_$calls, 1L)
  expect_identical(list(y = c_$value$y, y_pos = c_$value$y_pos),
                   tulpaObs:::.occu_cover_gen_shape(
                     lt, gen_replay("occu_cover", c_$calls[[1L]])))
  set_rng(9L)
  expect_identical(c_$value$y, tulpaObs:::.occu_cover_gen_shape(
    lt, tulpaObs:::.tobs_gen_occu_cover(lt, th))$y)
})


test_that("dyn_occu: the three routes draw through one generator", {
  a <- gen_calls("dynamic", simulate_dyn_occu(N = 15, J = 2, n.seasons = 3,
                                              seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("dynamic", ca))

  X1 <- matrix(1, 15, 1L, dimnames = list(NULL, "(Intercept)"))
  model <- list(y = ca$model$y, n_sites = 15L, n_seasons = 3L,
                X_processes = list(X1, X1, X1, X1),
                process_info = lapply(1:4, function(k) list(p = 1L)),
                data = a$value$data)
  par <- stats::qlogis(c(ca$theta$psi1[1L], ca$theta$p[1L], ca$theta$gam[1L],
                         ca$theta$eps[1L]))
  fit <- list(model = model, draws = tulpaObs:::.tobs_theta_row(par))
  set.seed(5)
  b <- gen_calls("dynamic", tulpaObs:::.tobs_simulate_dynamic(fit, 1L)[[1L]])
  expect_route("dynamic", b, ca$theta)

  c_ <- gen_calls("dynamic",
                  tulpaObs:::.TOBS_SBC_REGISTRY$dyn_occu$simulate(
                    list(model = model), par, 9L))
  expect_route("dynamic", c_, ca$theta, function(d) d$y)
})


test_that("int_occu: the three routes draw through one generator", {
  a <- gen_calls("integrated", simulate_int_occu(N.total = 30, n.shared = 6,
                                                 seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("integrated", ca))

  X_occ <- stats::model.matrix(~ x, a$value$data)
  X1 <- matrix(1, 30, 1L, dimnames = list(NULL, "(Intercept)"))
  model <- c(ca$model, list(
    X_processes = list(X_occ, X1, X1), data = a$value$data,
    process_info = list(list(p = 2L), list(p = 1L, name = "s1"),
                        list(p = 1L, name = "s2"))))
  par <- c(0.5, 0.3, 0.2, -0.1)
  fit <- list(model = model, draws = tulpaObs:::.tobs_theta_row(par))
  set.seed(5)
  b <- gen_calls("integrated", tulpaObs:::.tobs_simulate_integrated(fit, 1L)[[1L]])
  expect_route("integrated", b, ca$theta)

  c_ <- gen_calls("integrated",
                  tulpaObs:::.TOBS_SBC_REGISTRY$int_occu$simulate(
                    list(model = model), par, 9L))
  expect_route("integrated", c_, ca$theta,
               function(d) stats::setNames(d$y, NULL))
})


# A fit object for the families SBC reaches through `simulate(fit)`: the
# replicate at theta is that fit's own simulate() handler on a one-row draw
# matrix holding theta.
as_tobs_fit <- function(model, type, ...) {
  structure(c(list(model = c(model, list(model_type = type))), list(...)),
            class = "tobs_fit")
}

# The shared SBC replicate stores a count response as double; the counts
# themselves are what the generator drew.
sbc_counts <- function(d) {
  y <- d$y
  storage.mode(y) <- "integer"
  y
}

# The SBC route of a family registered on `.tobs_sbc_simple_entry()` or the
# shared default replicate.
sbc_simple_route <- function(key, family, fit, theta, seed = 9L) {
  spec <- list(fit_obs = fit,
               data_obs = list(cells = data.frame(row.names = seq_len(
                 fit$model$n_sites)), site = seq_len(fit$model$n_sites)))
  gen_calls(key, tulpaObs:::.TOBS_SBC_REGISTRY[[family]]$simulate(
    spec, theta, seed))
}


test_that("occu_ttd: the three routes draw through one generator", {
  a <- gen_calls("occu_ttd", simulate_occu_ttd(N = 15, J = 3, seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("occu_ttd", ca))

  fit <- as_tobs_fit(ca$model, "occu_ttd",
                     draws = tulpaObs:::.tobs_theta_row(ca$theta))
  set.seed(5)
  b <- gen_calls("occu_ttd", tulpaObs:::.tobs_simulate_occu_ttd(fit, 1L)[[1L]])
  expect_route("occu_ttd", b, ca$theta)

  c_ <- sbc_simple_route("occu_ttd", "occu_ttd", fit, ca$theta)
  expect_route("occu_ttd", c_, ca$theta, function(d) d$y)
})


test_that("royle_nichols: the three routes draw through one generator", {
  for (rv in list(NULL, 0.6)) {
    a <- gen_calls("royle_nichols",
                   simulate_royle_nichols(N = 20, J = 4, beta.r.visit = rv,
                                          seed = 1))
    ca <- a$calls[[1L]]
    expect_identical(a$value$y, gen_replay("royle_nichols", ca))

    fit <- as_tobs_fit(ca$model, "royle_nichols",
                       draws = tulpaObs:::.tobs_theta_row(ca$theta))
    set.seed(5)
    b <- gen_calls("royle_nichols",
                   tulpaObs:::.tobs_simulate_royle_nichols(fit, 1L)[[1L]])
    expect_route("royle_nichols", b, ca$theta)
    if (is.null(rv)) {
      c_ <- sbc_simple_route("royle_nichols", "royle_nichols", fit, ca$theta)
      expect_route("royle_nichols", c_, ca$theta, sbc_counts)
    }
  }
})


test_that("double_observer: the three routes draw through one generator", {
  for (type in c("independent", "dependent")) {
    a <- gen_calls("double_observer",
                   simulate_double_observer(N = 20, type = type, seed = 1))
    ca <- a$calls[[1L]]
    expect_identical(unname(a$value$y), gen_replay("double_observer", ca))

    fit <- as_tobs_fit(ca$model, "double_observer",
                       draws = tulpaObs:::.tobs_theta_row(ca$theta))
    set.seed(5)
    b <- gen_calls("double_observer",
                   tulpaObs:::.tobs_simulate_double_observer(fit, 1L)[[1L]])
    expect_route("double_observer", b, ca$theta)

    c_ <- sbc_simple_route("double_observer", "double_observer", fit, ca$theta)
    expect_route("double_observer", c_, ca$theta, sbc_counts)
  }
})


test_that("distsamp_open: the three routes draw through one generator", {
  a <- gen_calls("distsamp_open",
                 simulate_distsamp_open(N = 15, n.seasons = 3, seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(unname(a$value$y), gen_replay("distsamp_open", ca))

  X_l <- stats::model.matrix(~ abund_cov1, a$value$data)
  X_s <- stats::model.matrix(~ det_cov1, a$value$data)
  testthat::local_mocked_bindings(
    .dso_unpack = function(theta, model) {
      list(lambda = exp(as.vector(X_l %*% theta[1:2])),
           sigma = exp(as.vector(X_s %*% theta[3:4])),
           omega = ca$theta$omega, gamma = ca$theta$gamma)
    }, .package = "tulpaObs")
  par <- c(a$value$truth$beta_lambda, a$value$truth$beta_sigma)
  fit <- as_tobs_fit(c(ca$model, list(n_sites = 15L)), "distsamp_open",
                     draws = tulpaObs:::.tobs_theta_row(par))
  set.seed(5)
  b <- gen_calls("distsamp_open",
                 tulpaObs:::.tobs_simulate_distsamp_open(fit, 1L)[[1L]])
  expect_route("distsamp_open", b, ca$theta)
})


test_that("count: the three routes draw through one generator", {
  for (resp in c("poisson", "negbin", "gaussian", "binomial")) {
    a <- gen_calls("count", simulate_count(N = 30, response = resp, seed = 1))
    ca <- a$calls[[1L]]
    expect_identical(a$value$y, gen_replay("count", ca))

    phi <- if (identical(resp, "gaussian")) ca$theta$disp^2 else ca$theta$disp
    fit <- as_tobs_fit(c(ca$model, list(process_info = list(list(p = 2L)))),
                       "count", draws = tulpaObs:::.tobs_theta_row(ca$theta$beta),
                       count_dispersion = list(phi = phi))
    set.seed(5)
    b <- gen_calls("count", tulpaObs:::.tobs_simulate_count(fit, 1L)[[1L]])
    expect_route("count", b, ca$theta)
  }
})


test_that("ms_count / jsdm: the routes draw through one generator", {
  for (resp in c("poisson", "negbin", "gaussian")) {
    a <- gen_calls("ms_count", simulate_ms_count(N = 25, n.species = 4,
                                                 response = resp, seed = 1))
    ca <- a$calls[[1L]]
    expect_identical(a$value$y, gen_replay("ms_count", ca))

    disp <- switch(resp, negbin = list(r_s = ca$theta$disp),
                   gaussian = list(variance = ca$theta$disp^2), list())
    fit <- list(model = ca$model,
                ms_community = list(coef_mu = ca$theta$B), ms_dispersion = disp)
    set.seed(5)
    b <- gen_calls("ms_count", tulpaObs:::.tobs_simulate_ms_count(fit, 1L)[[1L]])
    expect_route("ms_count", b, ca$theta)
  }
  a <- gen_calls("ms_count", simulate_jsdm(N = 25, n.species = 4, seed = 1))
  ca <- a$calls[[1L]]
  rep_a <- gen_replay("ms_count", ca)
  storage.mode(rep_a) <- "integer"
  expect_identical(a$value$y, rep_a)
})


# The count / multistate families: the generator is a C++ kernel taking theta;
# simulate(fit) picks the posterior draw on the R side and hands it over.

test_that("abun: the three routes draw through one generator", {
  for (mix in c("poisson", "zip")) {
    a <- gen_calls("nmix", simulate_abun(N = 15, J = 3, mixture = mix, seed = 1))
    ca <- a$calls[[1L]]
    expect_identical(a$value$y, gen_replay("nmix", ca))

    zi <- identical(mix, "zip")
    fit <- as_tobs_fit(ca$model, "nmix",
                       draws = tulpaObs:::.tobs_theta_row(ca$theta$beta),
                       nmix_dispersion = list(r = NULL), zero_inflated = zi,
                       zi_omega = if (zi) 0.3 else NA_real_)
    set.seed(5)
    b <- gen_calls("nmix", tulpaObs:::.tobs_simulate_nmix(fit, 1L)[[1L]])
    expect_route("nmix", b, ca$theta)

    c_ <- sbc_simple_route("nmix", "abun", fit, ca$theta$beta)
    expect_route("nmix", c_, ca$theta, sbc_counts)
  }
})


test_that("removal: the three routes draw through one generator", {
  a <- gen_calls("removal", simulate_removal(N = 15, K = 3, seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("removal", ca))

  fit <- as_tobs_fit(ca$model, "removal",
                     draws = tulpaObs:::.tobs_theta_row(ca$theta$beta),
                     nmix_dispersion = list(r = NULL))
  set.seed(5)
  b <- gen_calls("removal", tulpaObs:::.tobs_simulate_removal(fit, 1L)[[1L]])
  expect_route("removal", b, ca$theta)

  c_ <- sbc_simple_route("removal", "removal", fit, ca$theta$beta)
  expect_route("removal", c_, ca$theta, sbc_counts)
})


test_that("distance: the three routes draw through one generator", {
  for (key in c("halfnorm", "hazard")) {
    a <- gen_calls("distance", simulate_distance(N = 15, key = key, seed = 1))
    ca <- a$calls[[1L]]
    expect_identical(a$value$y, gen_replay("distance", ca))

    fit <- as_tobs_fit(ca$model, "distance",
                       draws = tulpaObs:::.tobs_theta_row(ca$theta$beta),
                       nmix_dispersion = list(r = NULL),
                       distance_shape = list(shape = ca$theta$shape))
    set.seed(5)
    b <- gen_calls("distance", tulpaObs:::.tobs_simulate_distance(fit, 1L)[[1L]])
    expect_route("distance", b, ca$theta)

    c_ <- sbc_simple_route("distance", "distance", fit, ca$theta$beta)
    expect_route("distance", c_, ca$theta, sbc_counts)
  }
})


test_that("fp_occu: the three routes draw through one generator", {
  a <- gen_calls("fp_occu", simulate_fp_occu(N = 20, J = 3, seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("fp_occu", ca))

  fit <- as_tobs_fit(ca$model, "fp_occu",
                     draws = tulpaObs:::.tobs_theta_row(ca$theta$beta))
  set.seed(5)
  b <- gen_calls("fp_occu", tulpaObs:::.tobs_simulate_fp_occu(fit, 1L)[[1L]])
  expect_route("fp_occu", b, ca$theta)

  c_ <- sbc_simple_route("fp_occu", "fp_occu", fit, ca$theta$beta)
  expect_route("fp_occu", c_, ca$theta, sbc_counts)
})


test_that("dyn_abun: the three routes draw through one generator", {
  a <- gen_calls("dyn_abun", simulate_dyn_abun(N = 12, n.seasons = 3, J = 2,
                                               seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("dyn_abun", ca))

  fit <- as_tobs_fit(ca$model, "dyn_abun",
                     draws = tulpaObs:::.tobs_theta_row(ca$theta$beta),
                     mixture = "poisson")
  set.seed(5)
  b <- gen_calls("dyn_abun", tulpaObs:::.tobs_simulate_dyn_abun(fit, 1L)[[1L]])
  expect_route("dyn_abun", b, ca$theta)

  c_ <- sbc_simple_route("dyn_abun", "dyn_abun", fit, ca$theta$beta)
  expect_route("dyn_abun", c_, ca$theta, sbc_counts)

  # A zero-inflated draw reaches the generator as its structural-zero share.
  zfit <- fit
  zfit$draws <- tulpaObs:::.tobs_theta_row(c(ca$theta$beta, zi_logit = 0.4))
  zfit$mixture <- "zip"
  set.seed(5)
  z <- gen_calls("dyn_abun", tulpaObs:::.tobs_simulate_dyn_abun(zfit, 1L)[[1L]])
  expect_equal(z$calls[[1L]]$theta$zi, stats::plogis(0.4))
})


test_that("occu: the three routes draw through one generator", {
  a <- gen_calls("single", simulate_occu(N = 20, J = 3, seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("single", ca))

  fit <- as_tobs_fit(c(ca$model, list(n_sites = 20L)), "single",
                     draws = tulpaObs:::.tobs_theta_row(ca$theta$beta))
  set.seed(5)
  b <- gen_calls("single", tulpaObs:::.tobs_simulate_single(fit, 1L)[[1L]])
  expect_route("single", b, ca$theta)

  c_ <- sbc_simple_route("single", "occu", fit, ca$theta$beta)
  expect_route("single", c_, ca$theta, sbc_counts)

  # The per-visit detection group effect is the generator's detection offset.
  g <- gen_calls("single", simulate_occu(N = 20, J = 3, n.visit.groups = 4,
                                         seed = 1))
  cg <- g$calls[[1L]]
  expect_identical(g$value$y, gen_replay("single", cg))
  grp <- as.integer(g$value$visits$visit_group)
  expect_equal(cg$theta$off_det,
               matrix(g$value$truth$b_visit[grp], 20L, 3L, byrow = TRUE))
})


# The community families: the generator is a C++ kernel taking per-species,
# per-site parameters; simulate(fit) hands it the fitted community values.

test_that("ms_occu: the three routes draw through one generator", {
  a <- gen_calls("ms_occu", simulate_ms_occu(N = 20, J = 3, n.species = 3,
                                             seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("ms_occu", ca))

  tr <- a$value$truth
  X1 <- matrix(1, 20, 1L, dimnames = list(NULL, "(Intercept)"))
  sp <- paste0("sp", 1:3)
  model <- list(n_sites = 20L, max_visits = 3L, n_species = 3L,
                species_names = sp, valid = ca$model$valid,
                X_occ = stats::model.matrix(~ x, a$value$data), X_det = X1,
                process_info = list(list(p = 2L, coef_names = c("(Intercept)", "x")),
                                    list(p = 1L, coef_names = "(Intercept)")))
  fit <- as_tobs_fit(model, "ms_occu",
                     ms_community = list(coef_psi = tr$beta_species,
                                         coef_p = tr$alpha_species))
  set.seed(5)
  b <- gen_calls("ms_occu", tulpaObs:::.tobs_simulate_ms_occu(fit, 1L)[[1L]])
  expect_route("ms_occu", b, ca$theta)

  th <- stats::setNames(
    as.vector(rbind(t(tr$beta_species), t(tr$alpha_species))),
    as.vector(outer(c("psi_(Intercept)", "psi_x", "p_(Intercept)"), sp,
                    function(p, s) paste0(s, "_", p))))
  spec <- list(model = fit$model, fit_obs = fit,
               data_obs = list(cells = a$value$data, site = seq_len(20)))
  c_ <- gen_calls("ms_occu",
                  tulpaObs:::.TOBS_SBC_REGISTRY$ms_occu$simulate(spec, th, 9L))
  expect_route("ms_occu", c_, ca$theta, function(d) d$y)
})


test_that("ms_abun: the exported and fit routes draw through one generator", {
  for (mix in c("negbin", "zip")) {
    a <- gen_calls("ms_nmix", simulate_ms_abun(n.species = 3, N = 15, J = 3,
                                               mixture = mix, seed = 1))
    ca <- a$calls[[1L]]
    expect_identical(a$value$y, gen_replay("ms_nmix", ca))

    tr <- a$value$truth
    model <- list(n_sites = 15L, max_visits = 3L, n_species = 3L,
                  species_names = a$value$species, y = ca$model$y,
                  X_processes = list(
                    stats::model.matrix(~ abund_cov1, a$value$data),
                    stats::model.matrix(~ det_cov1, a$value$data)))
    fit <- as_tobs_fit(model, "ms_nmix", mixture = mix,
                       ms_community = list(coef_lambda = tr$beta_lambda,
                                           coef_p = tr$beta_p),
                       ms_dispersion = if (mix == "negbin") list(r_s = tr$r_s),
                       ms_zi = if (mix == "zip") list(omega_s = tr$omega_s))
    set.seed(5)
    b <- gen_calls("ms_nmix", tulpaObs:::.tobs_simulate_ms_nmix(fit, 1L)[[1L]])
    expect_route("ms_nmix", b, ca$theta)
  }
})


test_that("ms_dyn_occu: the exported and fit routes draw through one generator", {
  a <- gen_calls("ms_dyn_occu", simulate_ms_dyn_occu(N = 12, J = 2, n.species = 3,
                                                     n.seasons = 3, seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("ms_dyn_occu", ca))

  tr <- a$value$truth
  X1 <- matrix(1, 12, 1L, dimnames = list(NULL, "(Intercept)"))
  model <- list(n_sites = 12L, max_visits = 2L, n_seasons = 3L, n_species = 3L,
                species_names = NULL, valid = ca$model$valid,
                X_psi1 = X1, X_p = X1, X_gamma = X1, X_eps = X1,
                process_info = list(NULL, NULL,
                                    list(coef_names = "(Intercept)"),
                                    list(coef_names = "(Intercept)")))
  fit <- as_tobs_fit(model, "ms_dyn_occu",
                     means = c(`gamma_(Intercept)` = stats::qlogis(tr$gamma),
                               `eps_(Intercept)` = stats::qlogis(tr$epsilon)),
                     ms_community = list(
                       coef_psi1 = matrix(stats::qlogis(tr$psi1_species), ncol = 1L),
                       coef_p = matrix(stats::qlogis(tr$p_species), ncol = 1L)))
  set.seed(5)
  b <- gen_calls("ms_dyn_occu",
                 tulpaObs:::.tobs_simulate_ms_dyn_occu(fit, 1L)[[1L]])
  expect_route("ms_dyn_occu", b, ca$theta)
})


test_that("ms_int_occu: the exported and fit routes draw through one generator", {
  a <- gen_calls("ms_int_occu", simulate_ms_int_occu(N = 12, n.species = 3,
                                                     seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("ms_int_occu", ca))

  tr <- a$value$truth
  X1 <- matrix(1, 12, 1L, dimnames = list(NULL, "(Intercept)"))
  model <- c(ca$model[c("n_sources", "n_visits", "n_sites", "n_species", "valid")],
             list(process_names = c("a", "b"), species_names = NULL,
                  X_psi = X1, X_p = list(X1, X1)))
  fit <- as_tobs_fit(model, "ms_int_occu", ms_community = list(
    coef_psi = matrix(stats::qlogis(tr$psi_species), ncol = 1L),
    coef_a = matrix(stats::qlogis(tr$p_det[[1L]]), ncol = 1L),
    coef_b = matrix(stats::qlogis(tr$p_det[[2L]]), ncol = 1L)))
  set.seed(5)
  b <- gen_calls("ms_int_occu",
                 tulpaObs:::.tobs_simulate_ms_int_occu(fit, 1L)[[1L]])
  expect_route("ms_int_occu", b, ca$theta)
})


test_that("ms_occu_cover: the exported and fit routes draw through one generator", {
  a <- gen_calls("ms_occu_cover", simulate_ms_occu_cover(n.species = 3, N = 12,
                                                         J = 3, seed = 1))
  ca <- a$calls[[1L]]
  rep_a <- gen_replay("ms_occu_cover", ca)
  expect_identical(a$value$y, rep_a$y)
  expect_identical(a$value$y_pos, rep_a$y_pos)

  tr <- a$value$truth; vd <- a$value$visit_data
  model <- list(n_sites = 12L, max_visits = 3L, n_species = 3L,
                species_names = a$value$species, positive = "lognormal",
                valid = array(TRUE, c(12L, 3L, 3L)),
                X_occ = stats::model.matrix(~ occ_cov1, a$value$data),
                X_det_site = matrix(0, 12L, 0L), X_pos_site = matrix(0, 12L, 0L),
                X_det_visit = stats::model.matrix(~ det_cov1, vd),
                X_pos_visit = stats::model.matrix(~ pos_cov1, vd))
  fit <- as_tobs_fit(model, "ms_occu_cover", means = c(log_sigma_pos = log(0.4)),
                     ms_community = list(coef_occ = tr$beta_occ,
                                         coef_p = tr$beta_p,
                                         coef_pos = tr$beta_pos))
  set.seed(5)
  b <- gen_calls("ms_occu_cover",
                 tulpaObs:::.tobs_simulate_ms_occu_cover(fit, 1L)[[1L]])
  expect_length(b$calls, 1L)
  expect_equal(b$calls[[1L]]$theta, ca$theta, ignore_attr = TRUE,
               tolerance = 1e-12)
  expect_identical(b$value, gen_replay("ms_occu_cover", b$calls[[1L]]))
})


test_that("ms_occu_cover_spatial: the exported route draws through the generator", {
  adj <- (abs(outer(1:10, 1:10, "-")) == 1) * 1
  a <- gen_calls("ms_occu_cover",
                 simulate_ms_occu_cover_spatial(adj, n.species = 3, seed = 1))
  ca <- a$calls[[1L]]
  rep_a <- gen_replay("ms_occu_cover", ca)
  expect_identical(a$value$y, rep_a$y)
  expect_identical(a$value$y_pos, rep_a$y_pos)
})


test_that("ms_distance: the exported and fit routes draw through one generator", {
  a <- gen_calls("ms_distance", simulate_ms_distance(n.species = 3, N = 12,
                                                     seed = 1))
  ca <- a$calls[[1L]]
  expect_identical(a$value$y, gen_replay("ms_distance", ca))

  tr <- a$value$truth
  fit <- as_tobs_fit(ca$model, "ms_distance", ms_community = list(
    coef_lambda = tr$beta_lambda, coef_sigma = tr$beta_sigma))
  set.seed(5)
  b <- gen_calls("ms_distance",
                 tulpaObs:::.tobs_simulate_ms_distance(fit, 1L)[[1L]])
  expect_route("ms_distance", b, ca$theta)
})
