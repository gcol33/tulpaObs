# =============================================================================
# test-occu-cover-batch.R - batched multi-response occu_cover (gcol33/tulpa#66).
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
  # refines a peaked axis on the single fits but not the fused batch; see
  # gcol33/tulpaObs#58, gcol33/tulpa#69). The default batch backend is "looped"
  # (correct + fastest); this gate opts into the FUSED backend
  # (control$batch.backend = "fused") to validate it against independent fits.
  ctrl <- list(verbose = FALSE, max.iter = 200L, engine = "joint_coupled",
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
  }

  # The two species are genuinely different fits (not an accidental alias).
  expect_false(isTRUE(all.equal(batch$fits[["a"]]$means,
                                batch$fits[["b"]]$means)))

  # S3 surface.
  expect_type(coef(batch), "list")
  expect_named(coef(batch), c("a", "b"))
  expect_identical(tobs_get(batch, "b"), batch$fits[["b"]])
})
