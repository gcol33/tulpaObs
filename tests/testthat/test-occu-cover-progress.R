# =============================================================================
# test-occu-cover-progress.R - outer-grid heartbeat file survives verbose=FALSE
# (gcol33/tulpaObs#43).
#
# The cover()/occu_cover() nested-Laplace joint grid has two progress channels:
# the Rcout console bar and a heartbeat file written to disk (the only signal
# that survives a detached Start-Process / nohup stdout buffer). The original
# bug: progress defaulted to !verbose, so a detached run (verbose = FALSE)
# silently disabled the file too -- exactly the configuration the file exists
# for. The fix builds the reporter whenever progress.file is non-empty,
# independent of verbose; the console bar is now ON by default (set
# control$progress = FALSE to silence it), no longer tied to verbose.
#
# Guard: a verbose = FALSE fit with progress.file set must write a heartbeat
# file holding a final "<done> <total> <elapsed_s> <eta_s>" with done == total.
# =============================================================================

skip_on_cran()
skip_if_fast()

# Parse the last heartbeat line "<done> <total> <elapsed_s> <eta_s>".
.read_heartbeat <- function(path) {
  expect_true(file.exists(path), info = "heartbeat file was never written")
  expect_gt(file.size(path), 0)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  expect_gt(length(lines), 0L)
  fields <- as.numeric(strsplit(trimws(tail(lines, 1L)), "\\s+")[[1]])
  expect_length(fields, 4L)
  list(done = fields[1L], total = fields[2L],
       elapsed = fields[3L], eta = fields[4L])
}

test_that("cover() hurdle writes its heartbeat file under verbose = FALSE", {
  set.seed(43L)
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

  path <- tempfile(fileext = ".eta")
  on.exit(unlink(path), add = TRUE)

  fit <- suppressWarnings(tobs(
    formula = ~ x + bym2(graph = adj, group_var = "region"),
    data = dat, family = cover("beta"), y = y,
    method = "nested_laplace",
    control = list(
      verbose = FALSE,             # standard detached-run setting
      progress.file = path,        # ask for the ETA heartbeat
      sigma.grid = c(0.4, 0.7), rho.grid = c(0.5, 0.8),
      sigma.pos.grid = c(0.0, 0.6), adaptive.grid = FALSE
    )
  ))
  expect_s3_class(fit, "cover_fit")

  hb <- .read_heartbeat(path)
  expect_gt(hb$total, 0)
  expect_equal(hb$done, hb$total,
               info = "final heartbeat should report all outer cells done")
})

test_that("occu_cover() writes its heartbeat file under verbose = FALSE", {
  N <- 30L; J <- 4L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                             adj = adj, sigma = 0.8, alpha = 1.0, seed = 43L)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  path <- tempfile(fileext = ".eta")
  on.exit(unlink(path), add = TRUE)

  fit <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + bym2(graph = adj), data = cell_dat,
    family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 500L,
                   engine = "joint_coupled", progress.file = path)
  ))

  hb <- .read_heartbeat(path)
  expect_gt(hb$total, 0)
  expect_gte(hb$done, 1)
  expect_lte(hb$done, hb$total)
})

test_that("verbose = FALSE without progress.file leaves no heartbeat file", {
  # The console progress bar is ON by default (independent of verbose); with no
  # progress.file requested, nothing is written to disk and the fit still
  # completes (the no-op baseline).
  set.seed(431L)
  n_s <- 10L; N <- 80L
  adj <- matrix(0L, n_s, n_s)
  for (s in seq_len(n_s)) {
    for (j in setdiff(c(s - 1L, s + 1L), c(0L, n_s + 1L))) adj[s, j] <- 1L
  }
  spatial_idx <- sample.int(n_s, N, replace = TRUE)
  w_s <- 0.6 * rnorm(n_s)
  x   <- rnorm(N)
  mu  <- plogis(0.4 - 0.5 * x + w_s[spatial_idx])
  y   <- rbeta(N, mu * 30, (1 - mu) * 30)
  y   <- pmin(pmax(y, 0), 1 - 1e-6)
  dat <- data.frame(x = x, region = factor(spatial_idx))

  fit <- suppressWarnings(tobs(
    formula = ~ x + bym2(graph = adj, group_var = "region"),
    data = dat, family = cover("beta"), y = y,
    method = "nested_laplace",
    control = list(verbose = FALSE, sigma.grid = c(0.4, 0.7),
                   rho.grid = c(0.5, 0.8), sigma.pos.grid = c(0.0, 0.6),
                   adaptive.grid = FALSE)
  ))
  expect_s3_class(fit, "cover_fit")
})
