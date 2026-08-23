# The fitted latent field on the five count / observation families (abun,
# removal, distance, fp_occu, dyn_abun) is part of the arm's linear predictor,
# so every post-fit door that rebuilds eta from `model$X_processes` has to add it
# back: fitted() / residuals() / predict() and the pointwise log-likelihood
# behind WAIC / LOO / DIC / CPO (#254). The offset is recorded per process on
# `model$field_eta_offset` (R/field_offset.R).
#
# A field written in the `detection =` formula is arm-tagged, and a fitter that
# loads its field on the state arm regardless would silently fit a
# spatially-varying abundance where a spatially-varying detection was asked for
# (#255); every such fitter rejects it through one shared check.

.cffo_adj <- function(side) {
  ng <- side * side; co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i != j && abs(co$x[i]-co$x[j]) + abs(co$y[i]-co$y[j]) == 1L) adj[i, j] <- 1L
  adj
}
.cffo_field <- function(adj, sd_phi = 0.7, seed = 1L) {
  set.seed(seed); ng <- nrow(adj); phi <- as.numeric(scale(stats::rnorm(ng)))
  for (r in 1:3) { pn <- phi
    for (i in seq_len(ng)) { nb <- which(adj[i, ] == 1L); pn[i] <- 0.5*phi[i] + 0.5*mean(phi[nb]) }
    phi <- pn }
  phi <- sd_phi * as.numeric(scale(phi)); phi - mean(phi)
}
.cffo_nmix_data <- function(adj, phi, J = 5L, seed = 3L) {
  ng <- nrow(adj)
  set.seed(seed)
  x_ab <- stats::rnorm(ng); x_det <- stats::rnorm(ng)
  lambda <- exp(log(6) + 0.5 * x_ab + phi)
  p <- stats::plogis(0.4 + 0.4 * x_det)
  N <- stats::rpois(ng, lambda)
  y <- matrix(NA_integer_, ng, J)
  for (i in seq_len(ng)) y[i, ] <- stats::rbinom(J, N[i], p[i])
  list(y = y, data = data.frame(abund_cov1 = x_ab, det_cov1 = x_det),
       x_ab = x_ab, x_det = x_det)
}


test_that("a detection-formula areal field is rejected, not fitted on the state arm (#255)", {
  adj <- .cffo_adj(4L)
  d   <- .cffo_nmix_data(adj, .cffo_field(adj, seed = 11L), J = 4L)

  # abun() routes no arm on EITHER method, so both reject rather than loading a
  # spatially-varying abundance where a spatially-varying detection was asked for.
  for (m in c("nested_laplace", "nuts")) {
    expect_error(
      tobs(~ abund_cov1, data = d$data, family = abun(),
           detection = ~ det_cov1 + icar(graph = adj), y = d$y, method = m,
           control = list(verbose = FALSE, progress = FALSE)),
      "detection-arm field")
  }

  # The check itself, on the shared helper, for the arm tagging the parser sets.
  spec <- tulpaObs:::.tobs_term_icar(graph = adj)
  spec$shared <- c(FALSE, TRUE)
  expect_true(tulpaObs:::.tobs_det_arm_spatial(spec))
  expect_error(
    tulpaObs:::.tobs_reject_det_arm_spatial(spec, "x()", "abundance", "a logit",
                                            "is not wired."),
    "detection-arm field")
  spec$shared <- c(TRUE, FALSE)
  expect_false(tulpaObs:::.tobs_det_arm_spatial(spec))
  expect_silent(
    tulpaObs:::.tobs_reject_det_arm_spatial(spec, "x()", "abundance", "a logit",
                                            "is not wired."))
})


test_that("abun() areal field enters fitted() and the pointwise log-likelihood (#254)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .cffo_adj(5L)
  phi <- .cffo_field(adj, seed = 11L)
  d   <- .cffo_nmix_data(adj, phi)

  fit <- tobs(~ abund_cov1 + icar(graph = adj), data = d$data, family = abun(),
              detection = ~ det_cov1, y = d$y, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))

  off <- fit$model$field_eta_offset
  expect_length(off, 2L)
  expect_equal(off[[1L]], as.numeric(fit$spatial_field))
  expect_null(off[[2L]])                              # no detection-arm field

  # fitted() is the fit's own predictor, field included.
  fv <- fitted(fit)
  expect_equal(log(fv$lambda),
               as.vector(cbind(1, d$x_ab) %*% fit$means[1:2]) + fit$spatial_field)
  # The detection arm carries no field and is untouched.
  expect_equal(as.numeric(fv$p),
               stats::plogis(as.vector(fit$model$X_processes[[2]] %*% fit$means[3:4])))
  expect_equal(predict(fit), fv)

  # WAIC / LOO score the same predictor. Scoring the coefficients alone (the
  # pre-#254 behaviour, reproduced by dropping the slot) is markedly worse, so
  # the assertion cannot pass on a fit that ignores the field.
  # Pareto-k warns on a fixture this small; the comparison is the point here.
  elpd <- function(f) suppressWarnings(loo::loo(f))$estimates["elpd_loo", "Estimate"]
  blind <- fit; blind$model$field_eta_offset <- NULL
  expect_gt(elpd(fit), elpd(blind))

  # simulate() reaches the kernel through the same offset, as a design column
  # whose coefficient is pinned at 1.
  ab <- tulpaObs:::.tobs_sim_arm_block(fit$model, fit$draws, 2L)
  expect_identical(as.integer(ab$p[1L]),
                   as.integer(fit$model$process_info[[1L]]$p + 1L))
  expect_equal(as.numeric(ab$X[[1L]][, ncol(ab$X[[1L]])]), off[[1L]])
  expect_true(all(ab$draws[, ab$p[1L]] == 1))
  expect_identical(as.integer(ab$p[2L]),
                   as.integer(fit$model$process_info[[2L]]$p))
})


test_that("a detection-arm field lands on the detection arm's own row layout (#254)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .cffo_adj(6L); ng <- nrow(adj)
  phi <- .cffo_field(adj, sd_phi = 0.7, seed = 41L)
  set.seed(41L); x <- as.numeric(scale(stats::rnorm(ng)))
  lambda <- exp(log(9) + 0.5 * x)                    # abundance: no field
  p <- stats::plogis(0.2 + phi)                       # detection: spatial field
  K <- 4L; Nn <- stats::rpois(ng, lambda); y <- matrix(0L, ng, K); rem <- Nn
  for (k in 1:K) { y[, k] <- stats::rbinom(ng, rem, p); rem <- rem - y[, k] }

  # Control: the same term in the STATE formula labels and loads the state arm.
  fit <- tobs(~ x + icar(graph = adj), data = data.frame(x = x), family = removal(),
              detection = ~ 1, y = y, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$spatial_field_arm, "abundance")
  expect_null(fit$model$field_eta_offset[[2L]])

  fit <- tobs(~ x, data = data.frame(x = x), family = removal(),
              detection = ~ icar(graph = adj), y = y, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$spatial_field_arm, "detection")

  off <- fit$model$field_eta_offset
  expect_null(off[[1L]])
  # removal's detection design is one row per PASS, so the per-site field is
  # expanded through the model's site index -- the same expansion the fitter
  # applies inside the likelihood.
  expect_length(off[[2L]], length(fit$model$y_long))
  expect_equal(off[[2L]], as.numeric(fit$spatial_field)[fit$model$site_idx])
  expect_equal(stats::qlogis(as.numeric(fitted(fit)$p)),
               as.vector(fit$model$X_processes[[2]] %*%
                           fit$means[grep("^p_", names(fit$means))]) + off[[2L]])
  expect_true(is.finite(
    suppressWarnings(loo::loo(fit))$estimates["elpd_loo", "Estimate"]))
})


test_that("a sampled (NUTS) field enters fitted() and the criteria too (#254)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .cffo_adj(5L)
  phi <- .cffo_field(adj, seed = 11L)
  d   <- .cffo_nmix_data(adj, phi)

  fit <- tobs(~ abund_cov1 + car_proper(graph = adj), data = d$data,
              family = abun(), detection = ~ det_cov1, y = d$y, method = "nuts",
              control = list(verbose = FALSE, progress = FALSE, n.chains = 1L,
                             n.iter = 150L, n.warmup = 150L))
  expect_identical(fit$method, "nuts")
  off <- fit$model$field_eta_offset
  expect_equal(off[[1L]], as.numeric(fit$spatial_field))
  expect_equal(log(fitted(fit)$lambda),
               as.vector(cbind(1, d$x_ab) %*% fit$means[1:2]) + fit$spatial_field)
  expect_true(is.finite(
    suppressWarnings(loo::loo(fit))$estimates["elpd_loo", "Estimate"]))

  # logLik() reads the data log-likelihood at the posterior mean. The sampled
  # field is not one of the coefficient draw columns, so the family's own
  # evaluation runs it at offset 0; it is re-evaluated through the kernel the
  # criteria use, or logLik() / AIC() / BIC() would describe a model WAIC does
  # not score.
  par   <- stats::setNames(colMeans(fit$draws), colnames(fit$draws))
  blind <- fit$model; blind$field_eta_offset <- NULL
  expect_equal(as.numeric(stats::logLik(fit)),
               tulpaObs:::.tobs_laplace_marginal_loglik(fit$model, par)$loglik)
  expect_gt(as.numeric(stats::logLik(fit)),
            tulpaObs:::.tobs_laplace_marginal_loglik(blind, par)$loglik)
})


test_that("predict() adds the field in sample and not at a new design (#254)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .cffo_adj(6L); ng <- nrow(adj)
  phi <- .cffo_field(adj, sd_phi = 0.7, seed = 80L)
  set.seed(80L); x <- as.numeric(scale(stats::rnorm(ng)))
  psi <- stats::plogis(0.4 + 0.6 * x)                 # occupancy: field
  p11 <- stats::plogis(0.5)
  J <- 10L; z <- stats::rbinom(ng, 1, stats::plogis(0.4 + 0.6 * x + phi))
  y <- matrix(0L, ng, J)
  for (i in seq_len(ng)) for (j in seq_len(J)) {
    y[i, j] <- if (z[i] == 1L)
      sample(0:2, 1, prob = c(1 - p11, 0.2 * p11, 0.8 * p11))
    else sample(0:1, 1, prob = c(0.92, 0.08))
  }
  fit <- tobs(~ x + icar(graph = adj), data = data.frame(x = x), family = fp_occu(),
              detection = ~ 1, y = y, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  off <- fit$model$field_eta_offset[[1L]]
  expect_equal(off, as.numeric(fit$spatial_field))

  # `.tobs_predict_fp_occu()` rebuilds eta rather than routing to fitted(), so
  # the in-sample call has to carry the field itself.
  expect_equal(predict(fit, type = "psi"), fitted(fit)$psi)
  expect_equal(stats::qlogis(predict(fit, type = "psi")),
               as.vector(cbind(1, x) %*% fit$means[1:2]) + off)
  # A caller-supplied design is new rows with no field node, so it does not.
  X0 <- cbind(1, x)
  expect_equal(stats::qlogis(predict(fit, X.0 = X0, type = "psi")),
               as.vector(X0 %*% fit$means[1:2]))
})


test_that("a temporal-only field rides the same offset, mapped by period (#254)", {
  skip_on_cran()
  skip_if_fast()
  Tt <- 6L; per_t <- 12L; N <- Tt * per_t
  set.seed(101L)
  period <- rep(seq_len(Tt), each = per_t)
  rho <- 0.7; sig <- 0.5; u <- numeric(Tt)
  u[1] <- stats::rnorm(1, 0, sig / sqrt(1 - rho^2))
  for (t in 2:Tt) u[t] <- rho * u[t - 1] + stats::rnorm(1, 0, sig)
  u <- u - mean(u)
  x <- stats::rnorm(N)
  lambda <- exp(log(8) + 0.5 * x + u[period])
  K <- 4L; Nn <- stats::rpois(N, lambda); y <- matrix(0L, N, K); rem <- Nn
  for (k in 1:K) { y[, k] <- stats::rbinom(N, rem, 0.45); rem <- rem - y[, k] }

  fit <- tobs(~ x + temporal(period, type = "ar1"),
              data = data.frame(x = x, period = period), family = removal(),
              detection = ~ 1, y = y, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  off <- fit$model$field_eta_offset
  expect_null(fit$spatial_field)
  expect_length(off[[1L]], N)
  expect_equal(off[[1L]], as.numeric(fit$temporal_field)[period])
  expect_equal(log(fitted(fit)$lambda),
               as.vector(cbind(1, x) %*% fit$means[1:2]) + off[[1L]])
})


test_that("a fit with no field carries no offset and is unchanged (#254)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .cffo_adj(5L)
  d   <- .cffo_nmix_data(adj, .cffo_field(adj, seed = 11L))
  fit <- tobs(~ abund_cov1, data = d$data, family = abun(),
              detection = ~ det_cov1, y = d$y, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_null(fit$model$field_eta_offset)
  expect_null(tulpaObs:::.tobs_eta_offset(fit$model, 1L))
  expect_identical(log(fitted(fit)$lambda),
                   as.vector(cbind(1, d$x_ab) %*% fit$means[1:2]))

  # The simulation kernels take a field as a design COLUMN with its coefficient
  # pinned at 1, so with no field the arm block is exactly what they were handed
  # before and simulate() runs the same stream.
  pk <- vapply(fit$model$process_info, function(pp) pp$p, integer(1))
  ab <- tulpaObs:::.tobs_sim_arm_block(fit$model, fit$draws, 2L)
  expect_identical(ab$X[[1L]], fit$model$X_processes[[1L]])
  expect_identical(ab$X[[2L]], fit$model$X_processes[[2L]])
  expect_identical(ab$draws, fit$draws[, seq_len(sum(pk)), drop = FALSE])
  expect_identical(as.integer(ab$p), as.integer(pk))
})
