# =============================================================================
# test-occu-cover-checkpoint.R - grid-cell checkpoint/resume forwarded from
# occu_cover() (and the cover hurdle) into tulpa_nested_laplace_joint()
# (gcol33/tulpa#50 consumer wiring).
#
# control$checkpoint = list(path =, resume =) is passed verbatim to the joint
# engine, whose outer grid appends each completed cell to `path`. A checkpointed
# fit must equal an un-checkpointed one, a resume must load the finished cells
# (re-appending nothing) and reproduce the fit, and a torn tail from a killed
# run must be discarded and re-solved to the same answer.
# =============================================================================

skip_on_cran()
skip_if_fast()

.cc_build <- function(N = 30L, J = 4L, positive = "lognormal", seed = 12345L) {
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(
    N = N, J = J, positive = positive,
    adj = adj, sigma = 0.8, alpha = 1.0, seed = seed
  )
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  list(adj = adj, cell_dat = cell_dat, od = od, y_pos = y_pos, positive = positive)
}

.cc_fit <- function(d, checkpoint = NULL) {
  ctrl <- list(verbose = FALSE, max.iter = 500L, engine = "joint")
  if (!is.null(checkpoint)) ctrl$checkpoint <- checkpoint
  suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = d$adj), data = d$cell_dat,
    family = occu_cover(d$positive),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
    method = "nested_laplace", control = ctrl
  ))
}

.cc_expect_equiv <- function(fit_a, fit_b, tol = 1e-8) {
  expect_equal(as.numeric(fit_b$joint_fit$log_marginal),
               as.numeric(fit_a$joint_fit$log_marginal),
               tolerance = tol, info = "log_marginal differs")
  expect_equal(as.numeric(fit_b$means), as.numeric(fit_a$means),
               tolerance = tol, info = "posterior means differ")
  expect_equal(as.numeric(fit_b$spatial_field),
               as.numeric(fit_a$spatial_field),
               tolerance = tol, info = "spatial field differs")
}

test_that("occu_cover checkpoint fit equals an un-checkpointed fit", {
  d <- .cc_build(seed = 21L)
  path <- tempfile(fileext = ".ckpt")
  on.exit(unlink(path), add = TRUE)

  fit_plain <- .cc_fit(d)
  fit_ckpt  <- .cc_fit(d, checkpoint = list(path = path, resume = FALSE))
  .cc_expect_equiv(fit_plain, fit_ckpt)
  expect_true(file.exists(path))
  expect_gt(file.size(path), 16)
})

test_that("occu_cover resume loads completed cells and re-appends nothing", {
  d <- .cc_build(seed = 22L)
  path <- tempfile(fileext = ".ckpt")
  on.exit(unlink(path), add = TRUE)

  fit1 <- .cc_fit(d, checkpoint = list(path = path, resume = FALSE))
  size_after_1 <- file.size(path)
  fit2 <- .cc_fit(d, checkpoint = list(path = path, resume = TRUE))
  .cc_expect_equiv(fit1, fit2)
  expect_equal(file.size(path), size_after_1,
               info = "resume re-appended cells it should have loaded")
})

test_that("occu_cover resume after a torn tail re-solves to the same fit", {
  d <- .cc_build(seed = 23L)
  path <- tempfile(fileext = ".ckpt")
  on.exit(unlink(path), add = TRUE)

  fit_full <- .cc_fit(d, checkpoint = list(path = path, resume = FALSE))
  full_bytes <- readBin(path, "raw", n = file.size(path))
  keep <- as.integer(length(full_bytes) * 0.6)
  writeBin(full_bytes[seq_len(keep)], path)

  fit_resumed <- .cc_fit(d, checkpoint = list(path = path, resume = TRUE))
  .cc_expect_equiv(fit_full, fit_resumed)
})

test_that("cover hurdle forwards checkpoint into the joint engine", {
  # Single-y beta cover hurdle on a chain BYM2 field; the checkpoint flows
  # through the .dispatch_cover() -> fit_cover_hurdle_joint_nested() path.
  set.seed(24L)
  n_s <- 12L; N <- 120L
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) {
    for (j in setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L))) adj[s, j] <- 1L
  }
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  w_s   <- 0.6 * rnorm(n_s)
  x     <- rnorm(N)
  occur <- rbinom(N, 1, plogis(-0.3 + 0.7 * x + w_s[spatial_idx]))
  mu    <- plogis(0.4 - 0.5 * x + 1.0 * w_s[spatial_idx])
  y     <- numeric(N)
  ip    <- occur == 1L
  y[ip] <- rbeta(sum(ip), mu[ip] * 30, (1 - mu[ip]) * 30)
  y     <- pmin(pmax(y, 0), 1 - 1e-6)
  dat   <- data.frame(x = x, region = factor(spatial_idx))

  path <- tempfile(fileext = ".ckpt")
  on.exit(unlink(path), add = TRUE)

  run <- function(checkpoint = NULL) {
    ctrl <- list(sigma.grid = c(0.4, 0.7), rho.grid = c(0.5, 0.8),
                 adaptive.grid = FALSE)
    if (!is.null(checkpoint)) ctrl$checkpoint <- checkpoint
    suppressWarnings(tobs(
      formula = ~ x + bym2(graph = adj, group_var = "region"),
      data = dat, family = cover("beta"), y = y,
      method = "nested_laplace", control = ctrl
    ))
  }

  fit_plain <- run()
  fit_ckpt  <- run(checkpoint = list(path = path, resume = FALSE))
  expect_s3_class(fit_ckpt, "cover_fit")
  expect_true(file.exists(path))
  expect_gt(file.size(path), 16)
  expect_equal(as.numeric(fit_ckpt$joint$log_marginal),
               as.numeric(fit_plain$joint$log_marginal), tolerance = 1e-8,
               info = "cover-hurdle checkpoint changed the fit")

  size1 <- file.size(path)
  fit_resume <- run(checkpoint = list(path = path, resume = TRUE))
  expect_equal(file.size(path), size1,
               info = "cover-hurdle resume re-appended loaded cells")
  expect_equal(as.numeric(fit_resume$joint$log_marginal),
               as.numeric(fit_plain$joint$log_marginal), tolerance = 1e-8)
})

test_that("occu_cover checkpoint validates its argument via the engine", {
  d <- .cc_build(seed = 25L)
  expect_error(
    .cc_fit(d, checkpoint = list(resume = TRUE)),
    "path"
  )
  expect_error(
    .cc_fit(d, checkpoint = list(path = tempfile(), resume = "yes")),
    "resume"
  )
})
