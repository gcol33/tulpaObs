# =============================================================================
# test-cover-perarm.R - cover() per-arm formulas and share() placement.
#
# `cover(presence = ~ ..., positive = ~ ...)` gives each hurdle arm its own fixed
# effects (two independent designs); the single shared `formula` stays the
# back-compat spelling. A bare spatial() in one per-arm formula is that arm's own
# field. share(spatial()) in the positive formula reuses the presence field: with
# the default amplitude it reproduces the shared (both-arm) field byte-for-byte,
# while share(spatial(), alpha = ) sets one amplitude grid for the whole field and
# share(spatial(), terms = list(intercept = , trend = )) sets per-component grids
# so the intercept and trend blocks can decouple.
# =============================================================================

.cp_sim <- function(seed = 1L, N = 3000L) {
  set.seed(seed)
  x1 <- stats::rnorm(N); x2 <- stats::rnorm(N)
  pres <- stats::rbinom(N, 1L, stats::plogis(-0.2 + 1.0 * x1))
  mu   <- stats::plogis(-0.5 + 0.8 * x2)
  y <- numeric(N); pos <- pres == 1L
  y[pos] <- stats::rbeta(sum(pos), mu[pos] * 20, (1 - mu[pos]) * 20)
  y[y >= 1] <- 1 - 1e-6
  list(data = data.frame(x1 = x1, x2 = x2), y = y)
}

test_that("per-arm formulas give each arm its own fixed effects and recover", {
  skip_on_cran()
  s <- .cp_sim()
  fit <- tobs(presence = ~ x1, positive = ~ x2,
              family = cover(response = "beta"), data = s$data, y = s$y,
              method = "laplace")

  # Each arm carries ONLY its own covariate.
  expect_true("x1" %in% names(fit$beta_occ) && !("x2" %in% names(fit$beta_occ)))
  expect_true("x2" %in% names(fit$beta_pos) && !("x1" %in% names(fit$beta_pos)))
  # ... and recovers its truth.
  expect_lt(abs(fit$beta_occ[["x1"]] - 1.0), 0.15)
  expect_lt(abs(fit$beta_pos[["x2"]] - 0.8), 0.15)
})

test_that("the shared single formula stays back-compat (both arms share the FE)", {
  skip_on_cran()
  s <- .cp_sim()
  fit <- tobs(~ x1 + x2, family = cover(response = "beta"),
              data = s$data, y = s$y, method = "laplace")
  expect_true(all(c("x1", "x2") %in% names(fit$beta_occ)))
  expect_true(all(c("x1", "x2") %in% names(fit$beta_pos)))
})

# ---------------------------------------------------------------------------
# share() in the SINGLE shared formula (#298)
#
# The shared spelling has no positive formula to place a share() in, so a
# share() written in `formula` is the only way that door can state a coupling.
# It used to die inside the strip with `get1index`, naming no term the user
# wrote.
# ---------------------------------------------------------------------------

.cp_shared_sim <- function(seed = 7L, g = 4L, N = 400L) {
  set.seed(seed)
  n_cells <- g * g
  co  <- expand.grid(r = seq_len(g), c = seq_len(g))
  adj <- matrix(0L, n_cells, n_cells)
  for (i in seq_len(n_cells)) for (j in seq_len(n_cells)) {
    if (i < j && abs(co$r[i] - co$r[j]) + abs(co$c[i] - co$c[j]) == 1L) {
      adj[i, j] <- 1L; adj[j, i] <- 1L
    }
  }
  cell <- sample.int(n_cells, N, replace = TRUE)
  time <- as.numeric(scale(stats::rnorm(N)))
  z    <- as.numeric(scale(stats::rnorm(n_cells)))
  occ  <- stats::rbinom(N, 1L, stats::plogis(0.2 + 0.3 * time + 0.8 * z[cell]))
  cover <- ifelse(occ == 1L,
                  pmin(exp(-1.6 + 0.2 * time + 0.8 * z[cell] +
                             stats::rnorm(N, 0, 0.4)), 1 - 1e-6), 0)
  list(adj = adj, data = data.frame(cell = cell, time = time, cover = cover))
}

.cp_shared_ctrl <- list(verbose = FALSE, n.threads = 1L,
                        sigma.grid = exp(seq(log(0.3), log(1.5),
                                             length.out = 3)))

test_that("share() in the single shared cover() formula is accepted (#298)", {
  skip_on_cran()
  skip_if_fast()
  s   <- .cp_shared_sim()
  adj <- s$adj
  fit_at <- function(fml, ctrl = .cp_shared_ctrl) {
    suppressWarnings(suppressMessages(tobs(
      formula = fml, data = s$data, family = cover("lognormal"),
      y = s$data$cover, method = "nested_laplace", control = ctrl)))
  }
  alpha_of <- function(f) {
    h <- f$hyperpar$spatial
    unname(unlist(h[grep("alpha", names(h))]))[1L]
  }

  base_f <- ~ time + icar(graph = adj, group_var = "cell")

  # It fits at all -- this is what used to be `get1index`.
  f_bare <- fit_at(update(base_f, ~ . + share(spatial())))
  expect_s3_class(f_bare, "tobs_fit")

  # A bare share() asks for the engine's own axis, so it reproduces the fit
  # that writes no coupling at all. (When cover()'s no-coupling default flips
  # to pinned-at-zero per #297, THIS is the pair that has to be re-read: the
  # bare share() keeps the estimated axis and the bare formula loses it.)
  expect_equal(alpha_of(f_bare), alpha_of(fit_at(base_f)), tolerance = 1e-8)

  # A stated amplitude reaches the fit.
  expect_equal(alpha_of(fit_at(update(base_f, ~ . + share(spatial(),
                                                          alpha = 0.5)))),
               0.5, tolerance = 1e-8)

  # A stated grid is integrated rather than pinned: the reported amplitude is
  # inside the stated nodes.
  a_grid <- alpha_of(fit_at(update(base_f, ~ . + share(spatial(),
                                                       alpha = grid(c(0, 1))))))
  expect_gte(a_grid, 0); expect_lte(a_grid, 1)
})

test_that("shared-formula share() refuses the same conflicts the per-arm one does", {
  skip_on_cran()
  skip_if_fast()
  s   <- .cp_shared_sim()
  adj <- s$adj
  fit_at <- function(fml, ctrl = .cp_shared_ctrl) {
    suppressWarnings(suppressMessages(tobs(
      formula = fml, data = s$data, family = cover("lognormal"),
      y = s$data$cover, method = "nested_laplace", control = ctrl)))
  }
  # Both spellings of the amplitude at once.
  expect_error(
    fit_at(~ time + icar(graph = adj, group_var = "cell") + share(spatial()),
           ctrl = c(.cp_shared_ctrl, list(alpha.grid = c(0, 1)))),
    "not both")
  # A coupling with no field to couple names the field, not an index error.
  expect_error(fit_at(~ time + share(spatial())),
               "needs a spatial field")
})

test_that("a spatial field placed in the positive formula is an arm-specific field", {
  skip_on_cran()
  side <- 6L; nc <- side * side
  adj <- matrix(0L, nc, nc)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    if (r > 1L)   adj[idx(r, c), idx(r - 1L, c)] <- 1L
    if (r < side) adj[idx(r, c), idx(r + 1L, c)] <- 1L
    if (c > 1L)   adj[idx(r, c), idx(r, c - 1L)] <- 1L
    if (c < side) adj[idx(r, c), idx(r, c + 1L)] <- 1L
  }
  set.seed(3); reps <- 10L; cell <- rep(seq_len(nc), each = reps); N <- length(cell)
  x <- stats::rnorm(N); tt <- stats::rnorm(N)
  pres <- stats::rbinom(N, 1L, stats::plogis(0.2 + 0.5 * x))
  mu <- stats::plogis(0.3 * tt)
  y <- numeric(N); pos <- pres == 1L
  y[pos] <- stats::rbeta(sum(pos), mu[pos] * 15, (1 - mu[pos]) * 15); y[y >= 1] <- 1 - 1e-6
  dat <- data.frame(x = x, t = tt, cell = cell)
  ctrl <- list(progress = FALSE, integration = "ccd")

  fit_pa <- suppressWarnings(tobs(
    presence = ~ x + t, positive = ~ x + t + spatial(~ 1 || cell, graph = adj),
    family = cover(response = "beta"), data = dat, y = y,
    method = "nested_laplace", control = ctrl))

  expect_s3_class(fit_pa, "cover_fit")
  expect_true(isTRUE(fit_pa$armspecific))
  expect_identical(fit_pa$armspec_blocks[[1L]]$arm, "positive")
})

test_that("share(spatial()) in the positive formula == the shared-formula field (both arms)", {
  skip_on_cran()
  side <- 6L; nc <- side * side
  adj <- matrix(0L, nc, nc)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    if (r > 1L)   adj[idx(r, c), idx(r - 1L, c)] <- 1L
    if (r < side) adj[idx(r, c), idx(r + 1L, c)] <- 1L
    if (c > 1L)   adj[idx(r, c), idx(r, c - 1L)] <- 1L
    if (c < side) adj[idx(r, c), idx(r, c + 1L)] <- 1L
  }
  set.seed(7); reps <- 10L; cell <- rep(seq_len(nc), each = reps); N <- length(cell)
  x <- stats::rnorm(N); tt <- stats::rnorm(N)
  u <- stats::rnorm(nc, sd = 0.6)[cell]          # a genuine shared cell field
  pres <- stats::rbinom(N, 1L, stats::plogis(0.2 + 0.5 * x + u))
  mu <- stats::plogis(0.3 * tt + 0.7 * u)        # same field drives the cover arm
  y <- numeric(N); pos <- pres == 1L
  y[pos] <- stats::rbeta(sum(pos), mu[pos] * 15, (1 - mu[pos]) * 15); y[y >= 1] <- 1 - 1e-6
  dat <- data.frame(x = x, t = tt, cell = cell)
  ctrl <- list(progress = FALSE, integration = "ccd")

  # Canonical shared-field spelling: field on the presence arm, share() couples it.
  fit_copy <- suppressWarnings(tobs(
    presence = ~ x + t + spatial(~ 1 || cell, graph = adj),
    positive = ~ x + t + share(spatial()),
    family = cover(response = "beta"), data = dat, y = y,
    method = "nested_laplace", control = ctrl))
  # Shared-formula spelling: one bar in the shared formula reaches both arms.
  fit_to <- suppressWarnings(tobs(
    ~ x + t + spatial(~ 1 || cell, graph = adj),
    family = cover(response = "beta"), data = dat, y = y,
    method = "nested_laplace", control = ctrl))

  # Byte-identical across every numeric fit field (timings excluded).
  common <- intersect(names(fit_copy), names(fit_to))
  varying <- grepl("time|elapsed|second|runtime", common, ignore.case = TRUE)
  checked <- 0L
  for (nm in common[!varying]) {
    a <- fit_copy[[nm]]; b <- fit_to[[nm]]
    if (is.numeric(a) && is.numeric(b) && length(a) == length(b) && length(a) > 0L) {
      expect_equal(unname(a), unname(b), tolerance = 1e-8, info = nm)
      checked <- checked + 1L
    }
  }
  expect_gt(checked, 3L)   # guard: the loop actually compared coefficients + field
})

test_that("share() on the presence formula, or without a presence field, errors", {
  s <- .cp_sim(N = 200L)
  # share() on presence (wrong arm)
  expect_error(
    tobs(presence = ~ x1 + share(spatial()), positive = ~ x2,
         family = cover(response = "beta"), data = s$data, y = s$y,
         method = "nested_laplace"),
    "positive")
  # share() with no presence field to copy
  expect_error(
    tobs(presence = ~ x1, positive = ~ x2 + share(spatial()),
         family = cover(response = "beta"), data = s$data, y = s$y,
         method = "nested_laplace"),
    "needs a spatial field")
})

test_that("a temporal() / re() term in a per-arm formula is rejected with a pointer", {
  s <- .cp_sim(N = 200L)
  expect_error(
    tobs(presence = ~ x1 + (1 | x2), positive = ~ x2,
         family = cover(response = "beta"), data = s$data, y = s$y,
         method = "laplace"),
    "not routed per-arm")
})

test_that("only one per-arm formula errors (need both, or the shared one)", {
  s <- .cp_sim(N = 200L)
  expect_error(
    tobs(presence = ~ x1, family = cover(response = "beta"),
         data = s$data, y = s$y, method = "laplace"),
    "BOTH")
})


# ---------------------------------------------------------------------------
# share() coupling-amplitude handoff (whole-field + per-component decoupling).
# Guards the alpha -> control$alpha.grid[.trend] handoff the roadmap flagged as
# untested, and the per-component grammar that lets cover()'s intercept and
# trend blocks decouple (parity with occu_cover()).
# ---------------------------------------------------------------------------

# 6x6 grid with a shared intercept field u0 and a shared trend field u1
# (weighted by t), both driving occupancy and cover. The presence formula's
# `spatial(~ 1 + t || cell)` bar carries both blocks; share() transfers them.
.cov_trend_sim <- function(seed = 11L, side = 6L, reps = 12L) {
  nc <- side * side
  adj <- matrix(0L, nc, nc)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    if (r > 1L)   adj[idx(r, c), idx(r - 1L, c)] <- 1L
    if (r < side) adj[idx(r, c), idx(r + 1L, c)] <- 1L
    if (c > 1L)   adj[idx(r, c), idx(r, c - 1L)] <- 1L
    if (c < side) adj[idx(r, c), idx(r, c + 1L)] <- 1L
  }
  set.seed(seed); cell <- rep(seq_len(nc), each = reps); N <- length(cell)
  x <- stats::rnorm(N); tt <- stats::rnorm(N)
  u0 <- stats::rnorm(nc, sd = 0.6)[cell]; u1 <- stats::rnorm(nc, sd = 0.4)[cell]
  pres <- stats::rbinom(N, 1L, stats::plogis(0.2 + 0.5 * x + u0 + tt * u1))
  mu   <- stats::plogis(0.3 * x + 0.7 * u0 + 0.5 * tt * u1)
  y <- numeric(N); pos <- pres == 1L
  y[pos] <- stats::rbeta(sum(pos), mu[pos] * 15, (1 - mu[pos]) * 15); y[y >= 1] <- 1 - 1e-6
  list(data = data.frame(x = x, t = tt, cell = cell), y = y, adj = adj)
}

# Byte-equality across every numeric fit field (NA-in-both and timings excluded).
.cover_fit_equal <- function(a, b, tol = 1e-8) {
  common  <- intersect(names(a), names(b))
  varying <- grepl("time|elapsed|second|runtime", common, ignore.case = TRUE)
  checked <- 0L
  for (nm in common[!varying]) {
    av <- a[[nm]]; bv <- b[[nm]]
    if (is.numeric(av) && is.numeric(bv) && length(av) == length(bv) && length(av) > 0L) {
      keep <- !(is.na(av) & is.na(bv))
      if (any(keep)) {
        expect_equal(unname(av[keep]), unname(bv[keep]), tolerance = tol, info = nm)
        checked <- checked + 1L
      }
    }
  }
  checked
}

test_that("cover() whole-field share(alpha = grid()) == the shared-formula + control handoff", {
  skip_on_cran()
  s <- .cov_trend_sim(); g <- c(0, 0.5, 1)
  ctrl <- list(progress = FALSE, integration = "ccd")

  fit_copy <- suppressWarnings(tobs(
    presence = ~ x + spatial(~ 1 + t || cell, graph = s$adj),
    positive = ~ x + share(spatial(), alpha = grid(g)),
    family = cover(response = "beta"), data = s$data, y = s$y,
    method = "nested_laplace", control = ctrl))
  fit_ctrl <- suppressWarnings(tobs(
    ~ x + spatial(~ 1 + t || cell, graph = s$adj),
    family = cover(response = "beta"), data = s$data, y = s$y,
    method = "nested_laplace",
    control = c(ctrl, list(alpha.grid = g, alpha.grid.trend = g))))

  expect_s3_class(fit_copy, "cover_fit")
  expect_gt(.cover_fit_equal(fit_copy, fit_ctrl), 3L)
})

test_that("cover() per-component share(terms=) with equal grids == the whole-field share(alpha=)", {
  skip_on_cran()
  s <- .cov_trend_sim(); g <- c(0, 0.5, 1)
  ctrl <- list(progress = FALSE, integration = "ccd")

  fit_terms <- suppressWarnings(tobs(
    presence = ~ x + spatial(~ 1 + t || cell, graph = s$adj),
    positive = ~ x + share(spatial(), terms = list(intercept = grid(g), trend = grid(g))),
    family = cover(response = "beta"), data = s$data, y = s$y,
    method = "nested_laplace", control = ctrl))
  fit_whole <- suppressWarnings(tobs(
    presence = ~ x + spatial(~ 1 + t || cell, graph = s$adj),
    positive = ~ x + share(spatial(), alpha = grid(g)),
    family = cover(response = "beta"), data = s$data, y = s$y,
    method = "nested_laplace", control = ctrl))

  expect_gt(.cover_fit_equal(fit_terms, fit_whole), 3L)
})

test_that("cover() per-component share(terms=) decouples the intercept and trend blocks", {
  skip_on_cran()
  s <- .cov_trend_sim(); g <- c(0, 0.5, 1)
  ctrl <- list(progress = FALSE, integration = "ccd")

  # intercept coupled on grid g, trend pinned at alpha = 0 (decoupled).
  fit_dec <- suppressWarnings(tobs(
    presence = ~ x + spatial(~ 1 + t || cell, graph = s$adj),
    positive = ~ x + share(spatial(), terms = list(intercept = grid(g), trend = 0)),
    family = cover(response = "beta"), data = s$data, y = s$y,
    method = "nested_laplace", control = ctrl))
  # The low-level control spelling of the same coupling.
  fit_ctrl <- suppressWarnings(tobs(
    ~ x + spatial(~ 1 + t || cell, graph = s$adj),
    family = cover(response = "beta"), data = s$data, y = s$y,
    method = "nested_laplace",
    control = c(ctrl, list(alpha.grid = g, alpha.grid.trend = 0))))
  # Whole-field coupling of BOTH blocks -- must differ from the decoupled fit.
  fit_whole <- suppressWarnings(tobs(
    presence = ~ x + spatial(~ 1 + t || cell, graph = s$adj),
    positive = ~ x + share(spatial(), alpha = grid(g)),
    family = cover(response = "beta"), data = s$data, y = s$y,
    method = "nested_laplace", control = ctrl))

  # Per-component grids map exactly onto control$alpha.grid[.trend].
  expect_gt(.cover_fit_equal(fit_dec, fit_ctrl), 3L)
  # Decoupling the trend genuinely changes the fit.
  lm_dec   <- fit_dec$log_marginal   %||% fit_dec$logLik
  lm_whole <- fit_whole$log_marginal %||% fit_whole$logLik
  if (!is.null(lm_dec) && !is.null(lm_whole)) expect_false(isTRUE(all.equal(lm_dec, lm_whole)))
})

test_that("cover() share(terms=) rejects unknown, incomplete, and no-trend component specs", {
  s <- .cov_trend_sim(reps = 4L)
  ctrl <- list(progress = FALSE, integration = "ccd")

  # Unknown component name.
  expect_error(suppressWarnings(tobs(
    presence = ~ x + spatial(~ 1 + t || cell, graph = s$adj),
    positive = ~ x + share(spatial(), terms = list(intercept = grid(c(0, 1)), slope = 0)),
    family = cover(response = "beta"), data = s$data, y = s$y,
    method = "nested_laplace", control = ctrl)),
    "unknown field component")

  # Incomplete: trend block left unaddressed.
  expect_error(suppressWarnings(tobs(
    presence = ~ x + spatial(~ 1 + t || cell, graph = s$adj),
    positive = ~ x + share(spatial(), terms = list(intercept = grid(c(0, 1)))),
    family = cover(response = "beta"), data = s$data, y = s$y,
    method = "nested_laplace", control = ctrl)),
    "every field block")

  # A "trend" component named when the presence field is intercept-only.
  expect_error(suppressWarnings(tobs(
    presence = ~ x + spatial(~ 1 || cell, graph = s$adj),
    positive = ~ x + share(spatial(), terms = list(intercept = grid(c(0, 1)), trend = 0)),
    family = cover(response = "beta"), data = s$data, y = s$y,
    method = "nested_laplace", control = ctrl)),
    "no weighted trend block")
})
