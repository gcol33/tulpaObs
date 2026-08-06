# =============================================================================
# test-predict-spatial.R - tobs_predict_spatial(), the documented IDW-on-the-
# field spatial predictor (gcol33/tulpaObs#179).
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
  expect_equal(pr$q500, unname(stats::quantile(truth, 0.5)))

  # The quantile columns are named from the requested probabilities: the
  # percentage with its decimal point deleted, so 0.025 -> q25, 0.5 -> q500.
  expect_named(pr, c("mean", "sd", "q25", "q500", "q975"))
  pr2 <- tobs_predict_spatial(obj, matrix(c(0.5, 0), 1, 2),
                              quantiles = c(0.1, 0.9))
  expect_named(pr2, c("mean", "sd", "q100", "q900"))
  expect_equal(pr2$q100, unname(stats::quantile(truth, 0.1)))

  # One row per requested location.
  expect_s3_class(pr, "data.frame")
  expect_equal(nrow(tobs_predict_spatial(obj, cbind(c(0, 0.5, 1), 0))), 3L)
})


test_that("fixed effects enter through newocc.covs", {
  # Intercept 0.5 and slope 2 held fixed across draws, no field columns, so the
  # prediction is beta0 + beta1 * x exactly.
  d <- cbind(rep(0.5, 3), rep(2, 3))
  colnames(d) <- c("psi_(Intercept)", "psi_x")
  obj <- structure(list(
    draws = d, spatial = list(type = "icar"),
    model = list(process_info = list(list(name = "psi", p = 2L,
                                          link = "identity")))),
    class = c("tobs_fit", "tulpa_fit"))

  x <- c(0, 1, -0.5)
  pr <- tobs_predict_spatial(obj, cbind(x, 0), newocc.covs = data.frame(x = x))
  expect_equal(pr$mean, 0.5 + 2 * x)

  # Without covariates the design is intercept-only, so every location shares
  # the fixed-effect value (the field, when it enters, is what separates them).
  flat <- tobs_predict_spatial(obj, cbind(x, 0))
  expect_equal(flat$mean, rep(0.5, 3))
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

.tps_chain_adj <- function(N) {
  a <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) a[s, s - 1L] <- 1L
    if (s < N)  a[s, s + 1L] <- 1L
  }
  a
}

# 20 chain-linked cells x 3 sites, J = 3, with a smooth field on the cells.
.tps_fixture <- function(seed = 11L, n_cells = 20L, reps = 3L, J = 3L) {
  set.seed(seed)
  adj <- .tps_chain_adj(n_cells)
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

  newx <- c(-1, 0, 1)
  pr <- tobs_predict_spatial(fit, cbind(c(2, 8, 15), 0),
                             newocc.covs = data.frame(x = newx))

  # Independent recomputation of the occupancy-arm contribution from the same
  # draws. This pins WHICH draw columns the predictor reads: the psi block is
  # the leading process_info[[1]]$p columns, and reading one column further
  # would pull in the detection intercept.
  p_occ <- fit$model$process_info[[1]]$p
  X0 <- cbind(1, newx)
  eta <- fit$draws[, seq_len(p_occ), drop = FALSE] %*% t(X0)
  expect_equal(pr$mean, colMeans(stats::plogis(eta)), tolerance = 1e-8)
  expect_equal(pr$sd, apply(stats::plogis(eta), 2, stats::sd),
               tolerance = 1e-8)

  # Occupancy probabilities, and intervals in the right order.
  expect_true(all(pr$mean > 0 & pr$mean < 1))
  expect_true(all(pr$q25 <= pr$q500 & pr$q500 <= pr$q975))
})


test_that("predictions vary across locations once a field is fitted", {
  skip_on_cran()
  skip_if_fast()
  skip(paste(
    "tobs_predict_spatial() does not reach the fitted field on any route that",
    "produces one. Under laplace / nested_laplace `fit$draws` carries only the",
    "fixed-effect columns and an areal term's `fit$spatial$coords` is NULL, so",
    "the interpolation is skipped and every location gets the same prediction",
    "(measured on this fixture: the spread over three cells is exactly 0",
    "against a fitted field of sd 1.44). Under NUTS an areal field's columns",
    "are named spatial_field[i], which the branch's phi_spatial / w_gp",
    "patterns do not match, so it is skipped there too; the one pattern that",
    "does match a real column, gp_w, belongs to a gp() term whose coords are",
    "stored flattened, so that route errors on fit_coords[, 1] instead.",
    "Reported from gcol33/tulpaObs#179; unskip with the fix."))

  s <- .tps_fixture()
  fit <- .tps_fit_icar(s)

  # The truth on this fixture is a full sine wave across the cell chain, so
  # cells 2, 8 and 15 sit at genuinely different field values.
  pr <- tobs_predict_spatial(fit, cbind(c(2, 8, 15), 0))
  expect_gt(diff(range(pr$mean)), 0.05)

  # At a cell centre the interpolated field must reproduce that cell's own
  # fitted value, which is what makes this the field predictor rather than an
  # intercept-only one.
  at_nodes <- tobs_predict_spatial(fit, cbind(seq_along(fit$spatial_field), 0))
  b0 <- fit$means[["psi_(Intercept)"]]
  expect_gt(stats::cor(stats::qlogis(at_nodes$mean) - b0,
                       as.numeric(fit$spatial_field)), 0.95)
})
