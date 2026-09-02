# =============================================================================
# test-progress-all-variants.R - the outer-grid heartbeat file is emitted under
# verbose = FALSE for EVERY fitting backend, not just the nested-Laplace cover
# grid. Each backend reads the same scoped `tulpa.nl_progress` option tobs()
# sets from control$progress[.file]:
#   * EM-Laplace occupancy            -> tulpa's em_laplace R loop
#   * community EM (ms_occu)          -> community_em.R R loop
#   * NUTS (any family)               -> tulpa's HMC sampler (C++ GridProgress)
#   * count-marginal Laplace (abun)   -> in-tree Newton (C++ GridProgress)
#   * community N-mixture EM (ms_abun)-> in-tree EM (C++ GridProgress)
# The console bar is ON by default and independent of verbose; the heartbeat
# file is the detached-run channel and must appear whenever progress.file is set.
# =============================================================================

# Final heartbeat line "<done> <total> <elapsed_s> <eta_s>" parses to 4 numbers
# with done >= 1 and total > 0.
.expect_heartbeat <- function(path, label) {
  expect_true(file.exists(path), info = paste(label, "- no heartbeat file"))
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  expect_gt(length(ln), 0L)
  f <- as.numeric(strsplit(trimws(tail(ln, 1L)), "[[:space:]]+")[[1]])
  expect_length(f, 4L)
  expect_gte(f[1], 1)          # done
  expect_gt(f[2], 0)           # total
  expect_true(all(is.finite(f)))
}

test_that("occu laplace (EM) writes a heartbeat under verbose = FALSE", {
  skip_on_cran()
  skip_if_fast()
  set.seed(1)
  sim <- simulate_occu(N = 80, J = 3, seed = 1)
  path <- tempfile(fileext = ".eta"); on.exit(unlink(path), add = TRUE)
  tobs(~ occ_cov1, data = sim$data, family = occu(), detection = ~ det_cov1,
       y = sim$y, method = "laplace",
       control = list(verbose = FALSE, progress.file = path))
  .expect_heartbeat(path, "occu-laplace")
})

test_that("ms_occu (community EM) writes a heartbeat under verbose = FALSE", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_occu(N = 90, J = 3, n_species = 8, seed = 2)
  path <- tempfile(fileext = ".eta"); on.exit(unlink(path), add = TRUE)
  tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1, y = sim$y,
       species = paste0("sp", seq_len(8)), method = "laplace",
       control = list(verbose = FALSE, progress.file = path))
  .expect_heartbeat(path, "ms_occu-cem")
})

test_that("abun NUTS writes a heartbeat under verbose = FALSE", {
  skip_on_cran()
  skip_if_fast()
  set.seed(3)
  sim <- simulate_abun(N = 60, J = 4, n_abund_covs = 1, n_det_covs = 1, seed = 3)
  path <- tempfile(fileext = ".eta"); on.exit(unlink(path), add = TRUE)
  tobs(~ abund_cov1, data = sim$data, y = sim$y, family = abun(),
       detection = ~ det_cov1, method = "nuts",
       control = list(n.iter = 300L, n.warmup = 200L, seed = 7L,
                      verbose = FALSE, progress.file = path))
  .expect_heartbeat(path, "abun-nuts")
})

test_that("abun laplace (count-marginal Newton) writes a heartbeat", {
  skip_on_cran()
  skip_if_fast()
  set.seed(4)
  sim <- simulate_abun(N = 200, J = 4, n_abund_covs = 2, n_det_covs = 1, seed = 4)
  path <- tempfile(fileext = ".eta"); on.exit(unlink(path), add = TRUE)
  tobs(~ abund_cov1 + abund_cov2, data = sim$data, y = sim$y, family = abun(),
       detection = ~ det_cov1, method = "laplace",
       control = list(verbose = FALSE, progress.file = path))
  .expect_heartbeat(path, "abun-laplace")
})

test_that("ms_abun (community N-mixture EM) writes a heartbeat", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_abun(n_species = 10, N = 70, J = 4,
                          n_abund_covs = 1, n_det_covs = 1, seed = 5)
  path <- tempfile(fileext = ".eta"); on.exit(unlink(path), add = TRUE)
  tobs(~ abund_cov1, data = sim$data, y = sim$y, family = ms_abun(),
       detection = ~ det_cov1, species = sim$species, method = "laplace",
       control = list(verbose = FALSE, progress.file = path))
  .expect_heartbeat(path, "ms_abun-cem")
})

test_that("NUTS progress does not perturb the sampler (byte-exact draws)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(6)
  sim <- simulate_abun(N = 50, J = 4, n_abund_covs = 1, n_det_covs = 1, seed = 6)
  ctl <- list(n.iter = 300L, n.warmup = 200L, n.chains = 1L, seed = 9L,
              verbose = FALSE)
  off <- tobs(~ abund_cov1, data = sim$data, y = sim$y, family = abun(),
              detection = ~ det_cov1, method = "nuts", control = ctl)
  on  <- tobs(~ abund_cov1, data = sim$data, y = sim$y, family = abun(),
              detection = ~ det_cov1, method = "nuts",
              control = c(ctl, list(progress.file = tempfile(fileext = ".eta"))))
  expect_equal(off$draws, on$draws, tolerance = 0)
})
