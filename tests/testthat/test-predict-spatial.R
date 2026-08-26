# =============================================================================
# test-predict-spatial.R - tobs_predict_spatial(), the documented IDW-on-the-
# field spatial predictor.
#
# WHAT IS ASSERTED, AND AGAINST WHAT
#
# The interpolation is inverse-distance weighting over the k = 5 nearest fitted
# field nodes with weights 1 / (d + 1e-10), so every quantity below has a closed
# form and is asserted against it rather than against a stored number:
#
#   * at a node, the weight on that node is 1e10 against O(1) elsewhere, so the
#     prediction is that node's own value;
#   * with two nodes, w1 / (w1 + w2) = 1 - x along the connecting segment, so
#     p = 1 IDW is exactly linear interpolation there -- the midpoint is the
#     mean and the quarter point is the 3:1 weighted mean;
#   * with more than five nodes the sixth-nearest carries exactly zero weight,
#     and the fifth-nearest carries (1/d5) / sum_{j<=5} (1/dj).
#
# The state link is applied last, so the same eta is checked through logit, log
# and identity.
#
# SCOPE. Most blocks drive the interpolation through a minimal object carrying
# only what the function reads (`draws` with its field columns, `spatial$coords`,
# `model$process_info[[1]]`). That is a unit test of the interpolation kernel: it
# proves the weighting, the truncation and the link, and it would catch a change
# to any of them. It does NOT prove that a fitted model reaches that kernel --
# the end-to-end blocks at the bottom carry that, and one of them is currently
# skipped (see its message).
# =============================================================================

# Minimal object exposing exactly the fields tobs_predict_spatial() reads. The
# field columns are named `gp_w[i]`, one of the prefixes the function matches.
.tps_obj <- function(node_coords, node_values, link = "identity", beta = 0,
                     n_draws = 4L) {
  n_f <- nrow(node_coords)
  d <- matrix(0, n_draws, 1L + n_f)
  d[, 1L] <- beta
  for (k in seq_len(n_f)) d[, 1L + k] <- node_values[[k]]
  colnames(d) <- c("psi_(Intercept)", paste0("gp_w[", seq_len(n_f), "]"))
  pi1 <- list(name = "psi", p = 1L)
  if (!is.null(link)) pi1$link <- link
  structure(list(draws = d,
                 spatial = list(type = "gp", coords = node_coords),
                 model = list(process_info = list(pi1))),
            class = c("tobs_fit", "tulpa_fit"))
}

# IDW over the k nearest nodes, written independently of the implementation.
.tps_idw <- function(node_coords, node_values, pt, k = 5L) {
  d <- sqrt((node_coords[, 1] - pt[1])^2 + (node_coords[, 2] - pt[2])^2)
  k <- min(k, length(d))
  nn <- order(d)[seq_len(k)]
  w <- 1 / (d[nn] + 1e-10)
  sum(w / sum(w) * node_values[nn])
}


test_that("a prediction at a fitted node returns that node's field value", {
  nodes <- rbind(c(0, 0), c(1, 0), c(0, 1))
  vals  <- c(1, 3, -2)
  pr <- tobs_predict_spatial(.tps_obj(nodes, vals), nodes)

  # The 1e-10 distance offset leaves ~1e-10 of weight on the other nodes, so
  # this is exact to well within 1e-8, not merely close.
  expect_equal(pr$mean, vals, tolerance = 1e-8)
})


test_that("IDW between two nodes is their distance-weighted mean", {
  nodes <- rbind(c(0, 0), c(1, 0))
  vals  <- c(1, 3)
  obj <- .tps_obj(nodes, vals)

  # Equidistant -> the plain mean.
  expect_equal(tobs_predict_spatial(obj, matrix(c(0.5, 0), 1, 2))$mean, 2)

  # Quarter point: weights 1/0.25 and 1/0.75 normalise to 0.75 / 0.25.
  expect_equal(tobs_predict_spatial(obj, matrix(c(0.25, 0), 1, 2))$mean,
               0.75 * 1 + 0.25 * 3)

  # With two nodes, p = 1 IDW along the segment reduces to linear interpolation.
  xs <- seq(0, 1, by = 0.1)
  expect_equal(tobs_predict_spatial(obj, cbind(xs, 0))$mean, 1 + 2 * xs,
               tolerance = 1e-8)
})


test_that("IDW is monotone in distance to the higher-valued node", {
  nodes <- rbind(c(0, 0), c(1, 0))
  obj <- .tps_obj(nodes, c(1, 3))
  m <- tobs_predict_spatial(obj, cbind(seq(0, 1, by = 0.05), 0))$mean

  expect_true(all(diff(m) > 0))
  expect_true(all(m >= 1 - 1e-8 & m <= 3 + 1e-8))

  # Moving off the segment toward neither node leaves the value between them.
  off <- tobs_predict_spatial(obj, cbind(0.5, seq(0, 5, by = 0.5)))$mean
  expect_true(all(off >= 1 & off <= 3))
})


test_that("only the k = 5 nearest nodes contribute", {
  # Nodes on a line at distances 1..6 from the prediction point. The value 100
  # sits on the sixth-nearest, which the k = 5 truncation must exclude outright.
  nodes <- cbind(seq_len(6), 0)
  origin <- matrix(c(0, 0), 1, 2)

  far <- .tps_obj(nodes, c(0, 0, 0, 0, 0, 100))
  expect_equal(tobs_predict_spatial(far, origin)$mean, 0)

  # Moved onto the fifth-nearest it contributes exactly its normalised weight,
  # which pins the weight formula and not just the truncation rule.
  vals <- c(0, 0, 0, 0, 100, 0)
  w <- 1 / (seq_len(5) + 1e-10)
  inside <- .tps_obj(nodes, vals)
  expect_equal(tobs_predict_spatial(inside, origin)$mean,
               100 * w[5] / sum(w), tolerance = 1e-8)

  # Fewer nodes than k: every node contributes.
  two <- .tps_obj(nodes[1:2, , drop = FALSE], c(0, 100))
  expect_equal(tobs_predict_spatial(two, origin)$mean,
               .tps_idw(nodes[1:2, , drop = FALSE], c(0, 100), c(0, 0)),
               tolerance = 1e-8)
})


test_that("the state link is applied to the interpolated linear predictor", {
  nodes <- rbind(c(0, 0), c(1, 0))
  mid <- matrix(c(0.5, 0), 1, 2)
  eta <- 2  # the IDW value at the midpoint of nodes valued 1 and 3

  expect_equal(tobs_predict_spatial(.tps_obj(nodes, c(1, 3), "identity"), mid)$mean,
               eta)
  expect_equal(tobs_predict_spatial(.tps_obj(nodes, c(1, 3), "logit"), mid)$mean,
               stats::plogis(eta))
  expect_equal(tobs_predict_spatial(.tps_obj(nodes, c(1, 3), "log"), mid)$mean,
               exp(eta))

  # A count family returning logit-of-log-lambda is the failure this guards.
  expect_false(isTRUE(all.equal(
    tobs_predict_spatial(.tps_obj(nodes, c(1, 3), "log"), mid)$mean,
    tobs_predict_spatial(.tps_obj(nodes, c(1, 3), "logit"), mid)$mean)))

  # A process carrying no link falls back to logit (the occupancy default).
  expect_equal(tobs_predict_spatial(.tps_obj(nodes, c(1, 3), NULL), mid)$mean,
               stats::plogis(eta))

  # An unrecognised link is refused rather than silently treated as identity.
  expect_error(tobs_predict_spatial(.tps_obj(nodes, c(1, 3), "probit"), mid),
               "unsupported state-process link")
})


test_that("the summary columns are computed over the draws", {
  # Field node 2 varies across the five draws, node 1 is fixed at 1. At the
  # midpoint the prediction is the per-draw mean of the two, so mean / sd /
  # quantiles have exact closed forms.
  nodes <- rbind(c(0, 0), c(1, 0))
  node2 <- c(1, 2, 3, 4, 5)
  d <- cbind(0, 1, node2)
  colnames(d) <- c("psi_(Intercept)", "gp_w[1]", "gp_w[2]")
  obj <- structure(list(
    draws = d, spatial = list(type = "gp", coords = nodes),
    model = list(process_info = list(list(name = "psi", p = 1L,
                                          link = "identity")))),
    class = c("tobs_fit", "tulpa_fit"))

  truth <- (1 + node2) / 2
  pr <- tobs_predict_spatial(obj, matrix(c(0.5, 0), 1, 2))
  expect_equal(pr$mean, mean(truth))
  expect_equal(pr$sd, stats::sd(truth))
  expect_equal(pr$q50, unname(stats::quantile(truth, 0.5)))

  # The quantile columns are named from the requested probabilities as a
  # percentage, on stats::quantile()'s own convention: 0.025 -> q2.5, 0.5 ->
  # q50. predict() and predict_terms() name theirs the same way, so one level
  # reads as one column name across the three predictors (#242).
  expect_named(pr, c("mean", "sd", "q2.5", "q50", "q97.5"))
  expect_equal(attr(pr, "quantiles"), c(0.025, 0.5, 0.975))

  # A non-default `quantiles` names its columns for the levels it was given,
  # and those columns hold the values of those levels.
  pr2 <- tobs_predict_spatial(obj, matrix(c(0.5, 0), 1, 2),
                              quantiles = c(0.1, 0.9))
  expect_named(pr2, c("mean", "sd", "q10", "q90"))
  expect_equal(pr2$q10, unname(stats::quantile(truth, 0.1)))
  expect_equal(pr2$q90, unname(stats::quantile(truth, 0.9)))
  expect_equal(attr(pr2, "quantiles"), c(0.1, 0.9))

  # A malformed `quantiles` is named as the bad argument rather than reaching
  # quantile() as an NA.
  bad <- function(q) tobs_predict_spatial(obj, matrix(c(0.5, 0), 1, 2),
                                          quantiles = q)
  expect_error(bad(c(0.9, 0.1)), "strictly increasing")
  expect_error(bad(c(0, 0.5)), "strictly inside")
  expect_error(bad(c(0.1, NA)), "no NA")
  expect_error(bad("0.5"), "numeric vector")

  # One row per requested location.
  expect_s3_class(pr, "data.frame")
  expect_equal(nrow(tobs_predict_spatial(obj, cbind(c(0, 0.5, 1), 0))), 3L)
})


# A spatial fit whose fixed effects are held at `betas` across draws, with an
# identically-zero field: whatever the interpolation does it adds 0, so the
# fixed-effect half is what is being read. Carries the state formula and the
# fitted data, which is where the design is expanded from (#243).
.tps_fe_obj <- function(formula, data, coef_names, betas,
                        n_nodes = 3L, n_draws = 3L) {
  d <- matrix(rep(betas, each = n_draws), n_draws, length(betas))
  colnames(d) <- paste0("psi_", coef_names)
  structure(list(
    draws = d, spatial = list(type = "icar"),
    spatial_field = rep(0, n_nodes),
    model = list(
      formulas = list(occ = formula, det = ~ 1),
      data = data,
      process_info = list(list(name = "psi", p = length(betas),
                               link = "identity",
                               coef_names = coef_names)))),
    class = c("tobs_fit", "tulpa_fit"))
}

test_that("fixed effects enter through newocc.covs", {
  x <- c(0, 1, -0.5)
  obj <- .tps_fe_obj(~ x, data.frame(x = x - mean(x)),
                     c("(Intercept)", "x"), c(0.5, 2))
  nodes <- cbind(0:2, 0)

  pr <- tobs_predict_spatial(obj, cbind(x, 0), newocc.covs = data.frame(x = x),
                             node.coords = nodes)
  expect_equal(pr$mean, 0.5 + 2 * x)

  # Without covariates the design is intercept-only, so every location shares
  # the fixed-effect value (the field, when it enters, is what separates them).
  flat <- tobs_predict_spatial(obj, cbind(x, 0), node.coords = nodes)
  expect_equal(flat$mean, rep(0.5, 3))

  # A fit that declares a spatial component but carries no field at all is a
  # broken object, not an intercept-only prediction.
  no_field <- obj
  no_field$spatial_field <- NULL
  expect_error(tobs_predict_spatial(no_field, cbind(x, 0), node.coords = nodes),
               "no fitted field")
})


test_that("newocc.covs is paired with the coefficients by name, not position", {
  # beta_x = 2 and beta_z = -2, so swapping the two columns is visible in the
  # prediction if they are read positionally (#243).
  fit_d <- data.frame(x = c(-1, 0, 1), z = c(1, 0, -1))
  obj <- .tps_fe_obj(~ x + z, fit_d, c("(Intercept)", "x", "z"),
                     c(0.5, 2, -2))
  nodes <- cbind(0:2, 0)

  right <- data.frame(x = c(1, 0, 0.5), z = c(0, 1, -0.5))
  wrong <- right[, c("z", "x")]          # same data, user's column order

  pr_right <- tobs_predict_spatial(obj, cbind(0:2, 0), newocc.covs = right,
                                   node.coords = nodes)
  pr_wrong <- tobs_predict_spatial(obj, cbind(0:2, 0), newocc.covs = wrong,
                                   node.coords = nodes)
  expect_equal(pr_right$mean, 0.5 + 2 * right$x - 2 * right$z)
  expect_equal(pr_wrong$mean, pr_right$mean)
})


test_that("factor and poly terms expand on the basis the fit used", {
  nodes <- cbind(0:2, 0)

  # A factor expands through the fitted levels, so a character column and the
  # dummy the fit coded agree.
  fit_g <- data.frame(g = factor(c("a", "b", "a", "b")))
  obj_g <- .tps_fe_obj(~ g, fit_g, c("(Intercept)", "gb"), c(0.5, 1.5))
  g_new <- c("a", "b", "a")
  pr_g <- tobs_predict_spatial(obj_g, cbind(0:2, 0),
                               newocc.covs = data.frame(g = g_new),
                               node.coords = nodes)
  expect_equal(pr_g$mean, 0.5 + 1.5 * (g_new == "b"))

  # An unseen factor level cannot be placed on the fitted contrast.
  expect_error(tobs_predict_spatial(obj_g, cbind(0:2, 0),
                                    newocc.covs = data.frame(g = c("a", "c", "b")),
                                    node.coords = nodes),
               "new levels")

  # poly() is orthogonal on the FIT data, so its basis must come from there. If
  # it were rebuilt from whatever frame is passed, predicting at a subset would
  # not agree with the same rows of a prediction over the whole frame.
  fit_x <- data.frame(x = c(-2, -1, 0, 1, 2))
  obj_p <- .tps_fe_obj(~ poly(x, 2), fit_x,
                       c("(Intercept)", "poly(x, 2)1", "poly(x, 2)2"),
                       c(0.5, 2, -1), n_nodes = 5L)
  n5 <- cbind(0:4, 0)
  full <- tobs_predict_spatial(obj_p, cbind(0:4, 0), newocc.covs = fit_x,
                               node.coords = n5)
  sub  <- tobs_predict_spatial(obj_p, cbind(0:2, 0),
                               newocc.covs = fit_x[1:3, , drop = FALSE],
                               node.coords = n5)
  expect_equal(sub$mean, full$mean[1:3])
})


test_that("newocc.covs = NULL predicts at the fitted means, or says it cannot", {
  nodes <- cbind(0:2, 0)

  # Centred covariate: the intercept-only design IS the covariate mean.
  centred <- .tps_fe_obj(~ x, data.frame(x = c(-1, 0, 1)),
                         c("(Intercept)", "x"), c(0.5, 2))
  expect_equal(
    tobs_predict_spatial(centred, cbind(0:2, 0), node.coords = nodes)$mean,
    rep(0.5, 3))

  # Uncentred: zeros are not the mean, so predicting there would report at a
  # covariate value the fit never saw. That is refused rather than assumed.
  uncentred <- .tps_fe_obj(~ x, data.frame(x = c(1, 2, 3)),
                           c("(Intercept)", "x"), c(0.5, 2))
  expect_error(
    tobs_predict_spatial(uncentred, cbind(0:2, 0), node.coords = nodes),
    "not centred")
})


test_that("a fit without a spatial component is refused", {
  obj <- structure(list(
    draws = matrix(0, 2, 1, dimnames = list(NULL, "psi_(Intercept)")),
    spatial = NULL,
    model = list(process_info = list(list(name = "psi", p = 1L)))),
    class = c("tobs_fit", "tulpa_fit"))
  expect_error(tobs_predict_spatial(obj, matrix(c(0, 0), 1, 2)),
               "requires a model fitted with a spatial component")
})


# --- end to end, on a fitted model -------------------------------------------

# 20 chain-linked cells x 3 sites, J = 3, with a smooth field on the cells.
.tps_fixture <- function(seed = 11L, n_cells = 20L, reps = 3L, J = 3L) {
  set.seed(seed)
  adj <- chain_adj(n_cells)
  f0 <- sin(2 * pi * seq_len(n_cells) / n_cells)
  f0 <- f0 - mean(f0)
  n_sites <- n_cells * reps
  cell <- rep(seq_len(n_cells), each = reps)
  x <- as.numeric(scale(stats::rnorm(n_sites)))
  z <- stats::rbinom(n_sites, 1L,
                     stats::plogis(stats::qlogis(0.4) + 0.5 * x + f0[cell]))
  y <- matrix(0L, n_sites, J)
  for (i in seq_len(n_sites)) {
    y[i, ] <- stats::rbinom(J, 1L, z[i] * stats::plogis(0.5))
  }
  list(adj = adj, y = y, f0 = f0, data = data.frame(x = x, cell = cell))
}

.tps_fit_icar <- function(s) {
  tobs(~ x + icar(graph = s$adj, group_var = "cell"), detection = ~ 1,
       data = s$data, family = occu(), y = s$y, method = "nested_laplace",
       control = list(verbose = FALSE, progress = FALSE))
}


test_that("the fixed-effect half matches a recomputation from the draws", {
  skip_on_cran()
  skip_if_fast()

  s <- .tps_fixture()
  fit <- .tps_fit_icar(s)

  newx  <- c(-1, 0, 1)
  cells <- c(2, 8, 15)
  nodes <- cbind(seq_along(fit$spatial_field), 0)
  # Predicting AT node coordinates makes the interpolated field that node's own
  # value (IDW gives the coincident node all the weight), so the field half is
  # known in closed form and the fixed-effect half can be checked against it.
  pr <- tobs_predict_spatial(fit, cbind(cells, 0),
                             newocc.covs = data.frame(x = newx),
                             node.coords = nodes)

  # Independent recomputation of the occupancy-arm contribution from the same
  # draws. This pins WHICH draw columns the predictor reads: the psi block is
  # the leading process_info[[1]]$p columns, and reading one column further
  # would pull in the detection intercept.
  p_occ <- fit$model$process_info[[1]]$p
  X0 <- cbind(1, newx)
  eta <- fit$draws[, seq_len(p_occ), drop = FALSE] %*% t(X0)
  # This backend reports a posterior-mean field, so it enters every draw as the
  # same offset. Tolerance covers the 1e-10 IDW distance epsilon, which leaves
  # the coincident node a weight of 1 only to ~1e-10.
  eta <- sweep(eta, 2L, as.numeric(fit$spatial_field[cells]), "+")
  expect_equal(pr$mean, colMeans(stats::plogis(eta)), tolerance = 1e-6)
  expect_equal(pr$sd, apply(stats::plogis(eta), 2, stats::sd),
               tolerance = 1e-6)

  # Occupancy probabilities, and intervals in the right order.
  expect_true(all(pr$mean > 0 & pr$mean < 1))
  expect_true(all(pr$q25 <= pr$q500 & pr$q500 <= pr$q975))
})


test_that("predictions vary across locations once a field is fitted", {
  skip_on_cran()
  skip_if_fast()

  s <- .tps_fixture()
  fit <- .tps_fit_icar(s)
  # This fixture's graph is a chain, so node i is placed at (i, 0). An areal
  # field's nodes are graph vertices and carry no geometry of their own, so
  # this placement is the caller's to make.
  nodes <- cbind(seq_along(fit$spatial_field), 0)

  # Without it the call must refuse: an intercept-only prediction is
  # indistinguishable from a fit whose field happens to be flat.
  expect_error(tobs_predict_spatial(fit, cbind(c(2, 8, 15), 0)),
               "node\\.coords")

  # The truth on this fixture is a full sine wave across the cell chain, so
  # cells 2, 8 and 15 sit at genuinely different field values.
  pr <- tobs_predict_spatial(fit, cbind(c(2, 8, 15), 0), node.coords = nodes)
  expect_gt(diff(range(pr$mean)), 0.05)

  # At a cell centre the interpolated field must reproduce that cell's own
  # fitted value, which is what makes this the field predictor rather than an
  # intercept-only one.
  at_nodes <- tobs_predict_spatial(fit, nodes, node.coords = nodes)
  b0 <- fit$means[["psi_(Intercept)"]]
  expect_gt(stats::cor(stats::qlogis(at_nodes$mean) - b0,
                       as.numeric(fit$spatial_field)), 0.95)
})
