# =============================================================================
# test-occu-cover-batch.R
# - batched multi-response occu_cover.
#
# The batched-independent path fits B species in one tobs() call, each with the
# per-species model. The defining gate: a B-species batch must be PER-SPECIES
# bit-identical to B independent single-species fits at the same grid -- the
# batch only reorganises the work, it does not change the statistics.
#
# This test is written against the public API, so it gates BOTH the Stage-1
# looped backend (bit-identical by construction) and the future fused
# block-diagonal C++ backend (which must match to floating tolerance).
# =============================================================================


# ---- fast-tier unit tests on the multi-response plumbing (no model fits) ----

test_that("multi-response detection distinguishes single vs batched y", {
  # Single response: a plain matrix or a data.frame -> NULL (single-species).
  expect_null(tulpaObs:::.tobs_multiresponse_n(matrix(0, 4, 3)))
  expect_null(tulpaObs:::.tobs_multiresponse_n(NULL))
  expect_null(tulpaObs:::.tobs_multiresponse_n(
    as.data.frame(matrix(0, 4, 3))))

  # A list of >= 2 matrices, or a 3D array, is multi-response.
  expect_identical(tulpaObs:::.tobs_multiresponse_n(
    list(matrix(0, 4, 3), matrix(0, 4, 3))), 2L)
  expect_identical(tulpaObs:::.tobs_multiresponse_n(
    array(0, dim = c(4, 3, 5))), 5L)

  # A list with one matrix is degenerate (treated as a 1-batch list).
  expect_identical(tulpaObs:::.tobs_multiresponse_n(list(matrix(0, 4, 3))), 1L)
})


test_that("response slicing pulls one species out of list / 3D array", {
  m1 <- matrix(1, 4, 3); m2 <- matrix(2, 4, 3)
  expect_identical(tulpaObs:::.tobs_response_slice(list(m1, m2), 1L), m1)
  expect_identical(tulpaObs:::.tobs_response_slice(list(m1, m2), 2L), m2)

  arr <- array(0, dim = c(4, 3, 2))
  arr[, , 1L] <- 1; arr[, , 2L] <- 2
  expect_equal(tulpaObs:::.tobs_response_slice(arr, 1L), matrix(1, 4, 3))
  expect_equal(tulpaObs:::.tobs_response_slice(arr, 2L), matrix(2, 4, 3))
})


test_that("species labels resolve from species / names / default; validated", {
  y <- list(a = matrix(0, 2, 2), b = matrix(0, 2, 2))
  expect_identical(tulpaObs:::.tobs_batch_species_labels(NULL, y, 2L),
                   c("a", "b"))
  expect_identical(
    tulpaObs:::.tobs_batch_species_labels(c("x", "y"), y, 2L), c("x", "y"))
  expect_identical(
    tulpaObs:::.tobs_batch_species_labels(NULL,
      list(matrix(0, 2, 2), matrix(0, 2, 2)), 2L), c("sp1", "sp2"))
  expect_error(
    tulpaObs:::.tobs_batch_species_labels(c("a"), y, 2L), "length")
  expect_error(
    tulpaObs:::.tobs_batch_species_labels(c("a", "a"), y, 2L), "unique")
})


# ---- the bit-identity gate (heavy: real nested-Laplace joint fits) ----

test_that("2-species batch is per-species bit-identical to 2 independent fits", {
  skip_on_cran()
  skip_if_fast()

  N <- 24L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }

  # One shared design (covariates + graph), two response realisations. The
  # batch shares data/visits/formula across species; only y / y_pos differ.
  sim1 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 0.8, alpha = 1.0, seed = 101L)
  sim2 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 0.8, alpha = 1.0, seed = 202L)

  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim1$y)),
    det_cov1 = sim1$visit_data$det_cov1, pos_cov1 = sim1$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim1$data)

  y1 <- od$y;             yp1 <- sim1$y_pos; yp1[is.na(yp1)] <- 0
  y2 <- sim2$y;           yp2 <- sim2$y_pos; yp2[is.na(yp2)] <- 0

  # The fused batch integrates ONE shared fixed outer grid across species
  # (per-species refinement is not shareable). For the bit-identity comparison
  # both sides must integrate the same fixed tensor, which needs BOTH refinement
  # passes off: adaptive.grid = FALSE turns off edge-driven refinement, and
  # var.of.means.consistency = FALSE turns off the post-integration consistency
  # pass (joint-path default ON, independent of adaptive.grid -- it otherwise
  # refines a peaked axis on the single fits but not the fused batch; see). The
  # default batch backend is "looped" (correct + fastest); this gate opts into
  # the FUSED backend (control$batch.backend = "fused") to validate it against
  # independent fits. The default fit integrates the cover-arm dispersion on an
  # outer axis about each species' own pre-fit, so the two species reach the
  # fused driver with different dispersion nodes over one shared cell layout.
  ctrl <- list(verbose = FALSE, max.iter = 200L, engine = "joint",
               adaptive.grid = FALSE, var.of.means.consistency = FALSE,
               diagnose.k = FALSE)
  ctrl_fused <- c(ctrl, list(batch.backend = "fused"))

  fit_one <- function(yy, ypp) {
    suppressWarnings(tobs(
      formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
      family = occu_cover("lognormal"),
      detection = ~ det_cov1, positive = ~ pos_cov1,
      y = yy, y_pos = ypp, visits = od$det.covs,
      method = "nested_laplace", control = ctrl
    ))
  }

  batch <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = list(a = y1, b = y2), y_pos = list(yp1, yp2),
    visits = od$det.covs,
    method = "nested_laplace", control = ctrl_fused
  ))

  expect_s3_class(batch, "tobs_batch")
  expect_identical(batch$n_species, 2L)
  expect_identical(batch$species, c("a", "b"))
  expect_identical(batch$backend, "fused")
  phi_nodes <- lapply(batch$fits, function(f)
    sort(unique(f$joint_fit$theta_grid[, "phi_pos"])))
  expect_true(all(lengths(phi_nodes) > 1L))
  expect_false(isTRUE(all.equal(phi_nodes[["a"]], phi_nodes[["b"]])))

  # The DEFAULT backend (no batch.backend) is the looped path.
  batch_default <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = list(a = y1, b = y2), y_pos = list(yp1, yp2),
    visits = od$det.covs, method = "nested_laplace", control = ctrl
  ))
  expect_identical(batch_default$backend, "looped")
  expect_equal(batch_default$fits[["a"]]$means, batch$fits[["a"]]$means,
               tolerance = 1e-7)

  ind1 <- fit_one(y1, yp1)
  ind2 <- fit_one(y2, yp2)

  # Per-species equivalence to an independent fit. The fused backend only
  # reorganises the work, so the posterior summaries match. means / sds /
  # spatial_field are weighted sums over the outer grid (order-invariant); the
  # raw per-cell log_marginal vector is compared as a set (the fused cpp_grid
  # orders cells differently from the single-block grid, but it is the same set
  # of (sigma, alpha) cells with identical values).
  for (pair in list(list(batch$fits[["a"]], ind1),
                    list(batch$fits[["b"]], ind2))) {
    fb <- pair[[1L]]; fi <- pair[[2L]]
    expect_equal(fb$means, fi$means, tolerance = 1e-7)
    expect_equal(fb$sds,   fi$sds,   tolerance = 1e-7)
    expect_equal(fb$spatial_field, fi$spatial_field, tolerance = 1e-7)
    expect_equal(sort(fb$joint_fit$log_marginal),
                 sort(fi$joint_fit$log_marginal), tolerance = 1e-7)
    expect_equal(sort(unique(fb$joint_fit$theta_grid[, "phi_pos"])),
                 sort(unique(fi$joint_fit$theta_grid[, "phi_pos"])))
    expect_equal(sort(fb$joint_fit$weights), sort(fi$joint_fit$weights),
                 tolerance = 1e-7)
  }

  # The two species are genuinely different fits (not an accidental alias).
  expect_false(isTRUE(all.equal(batch$fits[["a"]]$means,
                                batch$fits[["b"]]$means)))

  # S3 surface.
  expect_type(coef(batch), "list")
  expect_named(coef(batch), c("a", "b"))
  expect_identical(tobs_get(batch, "b"), batch$fits[["b"]])
})


# The DEFAULT (looped) backend's own correctness, independent of the fused
# driver gated above -- this is the path every user who has not opted into
# control$batch.backend = "fused" actually runs, so it stays asserted directly
# against independent fits rather than transitively through a comparison to
# the fused output.
test_that("the default looped batch backend is per-species bit-identical to independent fits", {
  skip_on_cran()
  skip_if_fast()

  N <- 24L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim1 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 0.8, alpha = 1.0, seed = 101L)
  sim2 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 0.8, alpha = 1.0, seed = 202L)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim1$y)),
    det_cov1 = sim1$visit_data$det_cov1, pos_cov1 = sim1$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim1$data)
  y1 <- od$y;   yp1 <- sim1$y_pos; yp1[is.na(yp1)] <- 0
  y2 <- sim2$y; yp2 <- sim2$y_pos; yp2[is.na(yp2)] <- 0

  ctrl <- list(verbose = FALSE, max.iter = 200L, engine = "joint",
               adaptive.grid = FALSE, var.of.means.consistency = FALSE,
               diagnose.k = FALSE)
  fit_one <- function(yy, ypp) {
    suppressWarnings(tobs(
      formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
      family = occu_cover("lognormal"),
      detection = ~ det_cov1, positive = ~ pos_cov1,
      y = yy, y_pos = ypp, visits = od$det.covs,
      method = "nested_laplace", control = ctrl
    ))
  }

  batch_default <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = list(a = y1, b = y2), y_pos = list(yp1, yp2),
    visits = od$det.covs, method = "nested_laplace", control = ctrl
  ))
  expect_identical(batch_default$backend, "looped")

  ind1 <- fit_one(y1, yp1)
  ind2 <- fit_one(y2, yp2)

  for (pair in list(list(batch_default$fits[["a"]], ind1),
                    list(batch_default$fits[["b"]], ind2))) {
    fb <- pair[[1L]]; fi <- pair[[2L]]
    expect_equal(fb$means, fi$means, tolerance = 1e-7)
    expect_equal(fb$sds,   fi$sds,   tolerance = 1e-7)
    expect_equal(fb$spatial_field, fi$spatial_field, tolerance = 1e-7)
  }
})


# A one-node `phi.grid.pos` pins the cover dispersion. The fused driver takes
# each species' held dispersion from its own prepared arm, so a pinned batch
# holds the stated node for every species rather than each species' pre-fit, and
# carries no dispersion axis.
test_that("a fused batch with a one-node phi.grid.pos holds the node for every species", {
  skip_on_cran()
  skip_if_fast()

  N <- 24L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim1 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 0.8, alpha = 1.0, seed = 101L)
  sim2 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 0.8, alpha = 1.0, seed = 202L)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim1$y)),
    det_cov1 = sim1$visit_data$det_cov1, pos_cov1 = sim1$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim1$data)
  y1 <- od$y;   yp1 <- sim1$y_pos; yp1[is.na(yp1)] <- 0
  y2 <- sim2$y; yp2 <- sim2$y_pos; yp2[is.na(yp2)] <- 0

  ctrl_at <- function(node) list(
    verbose = FALSE, max.iter = 200L, engine = "joint", adaptive.grid = FALSE,
    var.of.means.consistency = FALSE, diagnose.k = FALSE, phi.grid.pos = node)
  fit <- function(yy, ypp, ctrl) suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = yy, y_pos = ypp, visits = od$det.covs,
    method = "nested_laplace", control = ctrl))

  batch <- fit(list(a = y1, b = y2), list(yp1, yp2),
               c(ctrl_at(0.4), list(batch.backend = "fused")))
  expect_identical(batch$backend, "fused")
  ind <- list(a = fit(y1, yp1, ctrl_at(0.4)), b = fit(y2, yp2, ctrl_at(0.4)))

  # The two species' pre-fits differ, so a batch that fell back to them could
  # not hold one value for both.
  pre <- vapply(list(y1 = list(y1, yp1), y2 = list(y2, yp2)), function(r)
    tulpaObs:::.tobs_joint_fit(fit(r[[1L]], r[[2L]], ctrl_at(NULL)))$responses$pos$phi,
    numeric(1))
  expect_false(isTRUE(all.equal(pre[[1L]], pre[[2L]])))

  for (sp in c("a", "b")) {
    fb <- batch$fits[[sp]]; fi <- ind[[sp]]
    expect_false("phi_pos" %in% colnames(fb$joint_fit$theta_grid))
    # Posterior draws read the fused fit like any other, and hold the node.
    d_b <- tulpaObs:::.tobs_joint_draws(fb, n = 100L)
    d_i <- tulpaObs:::.tobs_joint_draws(fi, n = 100L)
    expect_equal(d_b$disp, rep(0.4, 100L), tolerance = 1e-12)
    expect_equal(d_i$disp, rep(0.4, 100L), tolerance = 1e-12)
    expect_equal(fb$joint_fit$responses$pos$phi, 0.4^2, tolerance = 1e-12)
    expect_equal(fi$joint_fit$responses$pos$phi, 0.4^2, tolerance = 1e-12)
    expect_equal(fb$model$cover_pos_disp, 0.4, tolerance = 1e-12)
    expect_equal(fb$model$cover_pos_disp, fi$model$cover_pos_disp)
    expect_equal(fb$means, fi$means, tolerance = 1e-7)
    expect_equal(fb$sds,   fi$sds,   tolerance = 1e-7)
    expect_equal(fb$spatial_field, fi$spatial_field, tolerance = 1e-7)
    expect_equal(sort(fb$joint_fit$log_marginal),
                 sort(fi$joint_fit$log_marginal), tolerance = 1e-7)
  }

  # The node reaches the fused kernel: another node is another posterior.
  other <- fit(list(a = y1, b = y2), list(yp1, yp2),
               c(ctrl_at(0.7), list(batch.backend = "fused")))
  expect_false(isTRUE(all.equal(other$fits[["a"]]$means, batch$fits[["a"]]$means)))
})


# ---- a fused species fit is the object its own tobs() call returns ----

# The fit fields that record wall-clock time, which no two runs share, at any
# depth of a fit.
.batch_timing_fields <- c("timing", "diagnose_cost_ratio", "elapsed",
                          "run_time", "runtime")

.batch_drop_timing <- function(x) {
  if (!is.list(x)) return(x)
  keep <- !(names(x) %in% .batch_timing_fields)
  if (is.null(names(x))) keep <- rep(TRUE, length(x))
  at <- attributes(x)
  out <- lapply(x[keep], .batch_drop_timing)
  at$names <- names(x)[keep]
  attributes(out) <- at
  out
}

# The inner-layer Pareto k-hats a joint fit reports, and the fit without them.
.batch_inner_k <- function(fit) {
  jf <- fit$joint_fit
  list(inner = jf$inner_pareto_k, correction = jf$skew_correction$pareto_k)
}
.batch_drop_inner_k <- function(fit) {
  fit$joint_fit$inner_pareto_k <- NULL
  fit$joint_fit$skew_correction$pareto_k <- NULL
  fit
}

# The class and field names of `x`, recursively, as one flat record.
.batch_structure <- function(x, path = "fit") {
  here <- paste(path, paste(class(x), collapse = "/"))
  if (!is.list(x)) return(here)
  nm <- names(x) %||% rep("", length(x))
  c(here, unlist(lapply(seq_along(x), function(i)
    .batch_structure(x[[i]], paste0(path, "$", nm[[i]], "[", i, "]")))))
}

test_that("a fused species fit is the fit, draws, prediction and summary of its own tobs() call", {
  skip_on_cran()
  skip_if_fast()

  N <- 24L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim1 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 0.8, alpha = 1.0, seed = 101L)
  sim2 <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                              adj = adj, sigma = 0.8, alpha = 1.0, seed = 202L)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim1$y)),
    det_cov1 = sim1$visit_data$det_cov1, pos_cov1 = sim1$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim1$data)
  y1 <- od$y;   yp1 <- sim1$y_pos; yp1[is.na(yp1)] <- 0
  y2 <- sim2$y; yp2 <- sim2$y_pos; yp2[is.na(yp2)] <- 0

  # The engine defaults, refinement and placement included: every step of a
  # fused species' fit other than its main grid solve runs as it does alone.
  ctrl <- list(verbose = FALSE, max.iter = 200L, engine = "joint",
               diagnose.k = FALSE, progress = FALSE)
  fit <- function(yy, ypp, control) suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = yy, y_pos = ypp, visits = od$det.covs,
    method = "nested_laplace", control = control))

  # Each fit draws posterior samples, so the fused species and the independent
  # fits are run from one seed, in species order.
  set.seed(31L)
  batch <- fit(list(a = y1, b = y2), list(yp1, yp2),
               c(ctrl, list(batch.backend = "fused")))
  expect_identical(batch$backend, "fused")

  # The fused species share one design, which the per-species no-detection
  # compression would break, so the independent fits held to them here are
  # built uncompressed too.
  op <- options(tulpaObs.compress_nodet = FALSE)
  on.exit(options(op), add = TRUE)
  set.seed(31L)
  ind <- list(a = fit(y1, yp1, ctrl), b = fit(y2, yp2, ctrl))

  for (sp in c("a", "b")) {
    fb <- batch$fits[[sp]]; fi <- ind[[sp]]
    expect_identical(class(fb), class(fi), info = sp)
    expect_identical(class(fb$joint_fit), class(fi$joint_fit), info = sp)
    expect_identical(.batch_structure(.batch_drop_timing(fb)),
                     .batch_structure(.batch_drop_timing(fi)), info = sp)
    # A lone occu_cover cell hands the kernel its all-undetected (p, p)
    # curvature as a rank-1 term, a batched cell as the dense block, so the
    # fused kernel agrees with the lone one to the last few bits rather than bit
    # for bit (test-occu-cover-batch-fused.R); every number is held to that.
    # The inner Pareto k-hat is a tail-shape fit to a handful of points on the
    # probed curve and moves by thousandths under that perturbation, so it is
    # held at its own resolution. A formula's environment is the frame it was
    # written in, which for the batch is the batch call's.
    k_b <- .batch_inner_k(fb); k_i <- .batch_inner_k(fi)
    expect_equal(k_b, k_i, tolerance = 0.05, info = paste(sp, "inner k-hat"))
    expect_equal(.batch_drop_timing(.batch_drop_inner_k(fb)),
                 .batch_drop_timing(.batch_drop_inner_k(fi)),
                 tolerance = 1e-9, ignore_formula_env = TRUE, info = sp)

    set.seed(11L); d_b <- tulpaObs:::.tobs_joint_draws(fb, n = 200L)
    set.seed(11L); d_i <- tulpaObs:::.tobs_joint_draws(fi, n = 200L)
    expect_identical(.batch_structure(d_b), .batch_structure(d_i))
    expect_equal(d_b, d_i, tolerance = 1e-8, info = paste(sp, "draws"))

    set.seed(12L); p_b <- suppressWarnings(predict(fb))
    set.seed(12L); p_i <- suppressWarnings(predict(fi))
    expect_identical(.batch_structure(p_b), .batch_structure(p_i))
    expect_equal(p_b, p_i, tolerance = 1e-8, info = paste(sp, "predict"))

    set.seed(14L); f_b <- fitted(fb)
    set.seed(14L); f_i <- fitted(fi)
    expect_equal(f_b, f_i, tolerance = 1e-9, info = paste(sp, "fitted"))
    set.seed(15L); s_b <- .batch_drop_timing(summary(fb))
    set.seed(15L); s_i <- .batch_drop_timing(summary(fi))
    expect_identical(.batch_structure(s_b), .batch_structure(s_i))
    expect_equal(s_b, s_i, tolerance = 1e-9, ignore_formula_env = TRUE,
                 info = paste(sp, "summary"))
  }
  expect_false(isTRUE(all.equal(batch$fits[["a"]]$means,
                                batch$fits[["b"]]$means)))

  # Against the default (compressed) independent fits the posterior agrees to
  # the solver tolerance.
  options(op)
  set.seed(31L)
  ind_c <- list(a = fit(y1, yp1, ctrl), b = fit(y2, yp2, ctrl))
  for (sp in c("a", "b")) {
    fb <- batch$fits[[sp]]; fi <- ind_c[[sp]]
    expect_equal(fb$means, fi$means, tolerance = 1e-7)
    expect_equal(fb$sds,   fi$sds,   tolerance = 1e-7)
    set.seed(16L); f_b <- fitted(fb)
    set.seed(16L); f_i <- fitted(fi)
    expect_equal(f_b, f_i, tolerance = 1e-7)
    set.seed(13L); d_b <- tulpaObs:::.tobs_joint_draws(fb, n = 200L)
    set.seed(13L); d_i <- tulpaObs:::.tobs_joint_draws(fi, n = 200L)
    expect_equal(d_b$b, d_i$b, tolerance = 1e-6)
    expect_equal(d_b$disp, d_i$disp, tolerance = 1e-6)
  }
})
