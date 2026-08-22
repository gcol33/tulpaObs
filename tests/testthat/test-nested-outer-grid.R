# =============================================================================
# test-nested-outer-grid.R -- the outer-grid controls reach the multi-block
# nested-Laplace route.
#
# `.tobs_block_from_spatial()` / `_temporal()` / `_re()` read their grid
# overrides off the TERM object (`spatial$sigma_grid`, `spatial$tau_grid`, ...),
# and no constructor or binder ever wrote those fields: `.tobs_term()` adds
# `id` and `label` and nothing else. So twelve guards were dead and this route
# had no outer-grid override at all, while `cover()` and `occu_cover()` drove
# the same engine and exposed `control$sigma.grid` / `rho.grid` / `tau.grid` --
# and the control validator ADMITTED those names here, so a user setting one
# got it silently dropped.
#
# The grids now arrive as an argument threaded from `control`.
# =============================================================================

.nog_fixture <- function(n = 40L, seed = 4L) {
  sim <- simulate_occu(N = n, J = 3, seed = seed)
  adj <- matrix(0L, n, n)
  adj[cbind(seq_len(n - 1L), 2:n)] <- 1L
  adj[cbind(2:n, seq_len(n - 1L))] <- 1L
  list(sim = sim, adj = adj,
       data = transform(sim$data, cell = seq_len(nrow(sim$data))))
}

test_that("a supplied grid is the one the latent block carries", {
  fx <- .nog_fixture()
  m <- tobs(~ occ_cov1 + icar(graph = fx$adj, group_var = "cell"),
            data = fx$data, family = occu(), detection = ~ 1, y = fx$sim$y,
            method = "nested_laplace",
            control = list(verbose = FALSE, progress = FALSE))
  sp <- .tobs_structures_from_model(m$model)$spatial
  sf <- .tobs_resolve_occu_spatial_fields(sp, m$model)[[1L]]
  base <- .tobs_block_from_spatial(sf, 40L, seq_len(40L), m$model)
  set  <- .tobs_block_from_spatial(sf, 40L, seq_len(40L), m$model,
                                   grids = list(tau = c(0.5, 1, 2)))
  expect_equal(as.numeric(set$tau_grid), c(0.5, 1, 2))
  expect_false(isTRUE(all.equal(as.numeric(set$tau_grid),
                                as.numeric(base$tau_grid))))
  # No override leaves the block's own default untouched.
  expect_identical(.tobs_block_from_spatial(sf, 40L, seq_len(40L), m$model,
                                            grids = list())$tau_grid,
                   base$tau_grid)
})

test_that("control$tau.grid changes the fitted field", {
  skip_if_fast()
  skip_on_cran()
  fx <- .nog_fixture()
  f <- function(...) suppressWarnings(tobs(
    ~ occ_cov1 + icar(graph = fx$adj, group_var = "cell"),
    data = fx$data, family = occu(), detection = ~ 1, y = fx$sim$y,
    method = "nested_laplace",
    control = c(list(verbose = FALSE, progress = FALSE), list(...))))
  a <- f()
  b <- f(tau.grid = c(0.5, 1, 2))
  # A grid pinned well below the default's precision range leaves a wider
  # field; before the fix the two fits were bit-identical because the knob
  # never reached the block.
  expect_false(identical(a$spatial_field, b$spatial_field))
  expect_gt(stats::sd(b$spatial_field), stats::sd(a$spatial_field))
  # The same call with no override reproduces the default fit exactly.
  expect_identical(f()$spatial_field, a$spatial_field)
})

test_that("the outer-grid control names are admitted on this route", {
  for (nm in c("sigma.grid", "rho.grid", "tau.grid", "range.grid"))
    expect_true(nm %in% .tobs_control_groups$nested_laplace_joint, info = nm)
  expect_true(all(c("laplace_em", "nested_laplace_joint") %in%
                  .tobs_control_allow("nested_laplace", "none")))
  # The block builders take the grids as an argument; nothing reads them off
  # the term object, which is where they used to be looked for and never set.
  for (fn in c(".tobs_block_from_spatial", ".tobs_block_from_temporal",
               ".tobs_block_from_re")) {
    f <- get(fn, envir = asNamespace("tulpaObs"))
    expect_true("grids" %in% names(formals(f)), info = fn)
    src <- paste(deparse(body(f)), collapse = " ")
    for (term in c("spatial", "temporal", "re"))
      for (ax in c("sigma", "rho", "tau", "range"))
        expect_false(grepl(paste0(term, "$", ax, "_grid"), src, fixed = TRUE),
                     info = paste(fn, term, ax))
  }
})
