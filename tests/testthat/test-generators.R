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
