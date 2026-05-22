# End-to-end reproducer for tulpa#25: spatial_car(level='group') with a
# data set covering only a subset of adjacency cells. Pre-fix this errored
# at validate_spatial(); post-fix it must run to convergence and report a
# spatial field of length n_spatial_units (= 9), not n_obs (= 5).

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(tulpa)
})

set.seed(123)

# 3x3 raster: 9 cells, rook-or-diagonal adjacency
adj <- matrix(0L, 9, 9)
for (i in 1:3) for (j in 1:3) {
  c1 <- (i - 1) * 3 + j
  for (di in -1:1) for (dj in -1:1) {
    if (di == 0 && dj == 0) next
    ii <- i + di; jj <- j + dj
    if (ii >= 1 && ii <= 3 && jj >= 1 && jj <= 3) {
      c2 <- (ii - 1) * 3 + jj
      adj[c1, c2] <- 1L
    }
  }
}

# Observations in only 5 of the 9 cells
dat <- data.frame(
  cell_idx = c(1L, 2L, 4L, 5L, 8L),
  x        = rnorm(5),
  cover    = runif(5, 0.1, 0.9)
)

spec <- tulpa::spatial_car(adj, level = "group", group_var = "cell_idx")

cat("Spec built. n_spatial =", spec$n_spatial, "\n")
cat("Calling validate_spatial directly... ")
tulpa:::validate_spatial(spec, dat)
cat("OK\n")

cat("Calling prior_from_spec to inspect spatial_idx... ")
prior <- tulpa::prior_from_spec(spec, dat)
cat("OK\n")
cat("  spatial_idx =", prior$spatial_idx, "\n")
cat("  n_spatial_units =", prior$n_spatial_units, "\n")

stopifnot(identical(prior$spatial_idx, c(1L, 2L, 4L, 5L, 8L)))
stopifnot(prior$n_spatial_units == 9L)

# Now try the cover-hurdle nested_laplace fit as in the issue body.
fit <- tryCatch(
  tulpaObs::tobs(
    formula = ~ x,
    data    = dat,
    family  = tulpaObs::cover(positive = "beta"),
    y       = dat$cover,
    spatial = spec,
    engine  = "nested_laplace",
    control = list(verbose = FALSE)
  ),
  error = function(e) e
)

if (inherits(fit, "error")) {
  cat("tobs() errored:\n  ", conditionMessage(fit), "\n", sep = "")
} else {
  cat("tobs() returned a fit (class:", paste(class(fit), collapse = "/"), ")\n")
}
