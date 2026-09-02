# The multiscale cover model is built around a copy coefficient: `alpha` is in
# `means`, in `draws`, in .occu_mscale_cover_surface_at(), and the simulator
# draws the cover arm as `... + alpha * sigma * f[plot_cell]`. Its cover formula
# could not name it until gcol33/tulpaObs#294 -- the raw formula reached
# .occu_cover_reject_structured(), whose list ends in "copy", because this door
# skipped the strip occu_cover() does. These tests pin the route open.

# What reaches the joint fitter, without paying for a fit.
mscopy_control <- function(pos, method = "nested_laplace",
                           sim = mscopy_sim(), control = list()) {
  ns   <- asNamespace("tulpaObs")
  nm   <- ".tobs_fit_occu_multiscale_cover_joint"
  orig <- get(nm, envir = ns)
  cap  <- new.env(parent = emptyenv())
  unlockBinding(nm, ns)
  assign(nm, function(...) { cap$args <- list(...); structure(list(), class = "stub") },
         envir = ns)
  on.exit({ assign(nm, orig, envir = ns); lockBinding(nm, ns) }, add = TRUE)
  tobs(formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
       data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
       detection = ~ x_pdet, availability = ~ x_plot, positive = pos,
       y = sim$y, y_pos = sim$y_pos, method = method, control = control)
  cap$args
}

mscopy_sim <- function() {
  simulate_occu_multiscale_cover(n_cells = 12L, plots_per_cell = 3L,
                                 visits_per_plot = 2L, seed = 1L)
}


test_that("occu_multiscale_cover() accepts copy() on the cover formula", {
  # A bare copy() asks for the engine's own axis, so it states neither key.
  a <- mscopy_control(~ x_cov + copy(spatial()))
  expect_null(a$alpha.grid)
  expect_null(a$alpha.n)
})

test_that("occu_multiscale_cover(): copy() states nodes, a resolution, or a pin", {
  a <- mscopy_control(~ x_cov + copy(spatial(), alpha = grid(c(0.25, 0.5, 1))))
  expect_equal(as.numeric(a$alpha.grid), c(0.25, 0.5, 1))
  expect_null(a$alpha.n)

  a <- mscopy_control(~ x_cov + copy(spatial(), alpha = grid(n = 9)))
  expect_identical(a$alpha.n, 9L)
  expect_null(a$alpha.grid)

  # A scalar pins the amplitude: a one-node axis.
  a <- mscopy_control(~ x_cov + copy(spatial(), alpha = 0.5))
  expect_equal(as.numeric(a$alpha.grid), 0.5)
})

test_that("occu_multiscale_cover(): a fit naming no copy is left alone", {
  # This door's no-copy meaning is the engine's DEFAULT amplitude axis, not
  # occu_cover()'s pinned alpha = 0, so adding the copy() route must not touch
  # a fit that writes none. Measured against HEAD before the route landed: the
  # translation run unconditionally pins alpha.grid = 0 and silently decouples
  # every such fit (gcol33/tulpaObs#297 tracks the inconsistency itself).
  a <- mscopy_control(~ x_cov)
  expect_null(a$alpha.grid)
  expect_null(a$alpha.n)
})

test_that("occu_multiscale_cover(): copy() is refused where no field enters", {
  # method = "laplace" / "nuts" are the non-spatial paths (iid cells, field
  # fixed at 0), so a copy amplitude has nothing to scale.
  for (m in c("laplace", "nuts")) {
    expect_error(mscopy_control(~ x_cov + copy(spatial()), method = m),
                 "non-spatial path")
  }
})

test_that("occu_multiscale_cover(): copy() and the control spelling are exclusive", {
  expect_error(
    mscopy_control(~ x_cov + copy(spatial(), alpha = grid(c(0.5, 1))),
                   control = list(alpha.grid = c(0.5, 1))),
    "not both")
})

test_that("occu_multiscale_cover(): a copy() with no field to copy is refused", {
  sim <- mscopy_sim()
  expect_error(
    tobs(formula = ~ x_cell, data = sim$data,
         family = occu_multiscale_cover(response = "lognormal"),
         detection = ~ x_pdet, availability = ~ x_plot,
         positive = ~ x_cov + copy(spatial()),
         y = sim$y, y_pos = sim$y_pos, method = "nested_laplace"),
    "areal")
})
