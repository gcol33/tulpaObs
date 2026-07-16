# Dumps the outer hyperparameter grid for the single-season areal path so the
# recovery test's comments quote the real axis instead of assuming its
# direction. The note recorded tau in [0.3, 30] over 9 cells; this reads the
# axis, its ordering, and where the weight lands for a null vs an interior field.
devtools::load_all(".")  # run from the repo root

chain_adj <- function(N) {
  a <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) a[s, s - 1L] <- 1L
    if (s < N)  a[s, s + 1L] <- 1L
  }
  a
}
smooth_field <- function(N, sd_target, phase) {
  f <- sin(2 * pi * (seq_len(N) / N) + phase)
  f <- f - mean(f)
  f * (sd_target / stats::sd(f))
}
sim <- function(seed, sigma_f, n_cells = 40L, reps = 6L, J = 4L, b_x = 0.5) {
  set.seed(seed)
  adj <- chain_adj(n_cells)
  f0 <- if (sigma_f > 0) smooth_field(n_cells, sigma_f, phase = 0.7) else rep(0, n_cells)
  n_sites <- n_cells * reps
  cell <- rep(seq_len(n_cells), each = reps)
  x <- as.numeric(scale(stats::rnorm(n_sites)))
  z <- stats::rbinom(n_sites, 1L, plogis(stats::qlogis(0.35) + b_x * x + f0[cell]))
  y <- matrix(0L, n_sites, J)
  for (i in seq_len(n_sites)) y[i, ] <- stats::rbinom(J, 1L, z[i] * plogis(0.4))
  list(adj = adj, y = y, f0 = f0, data = data.frame(x = x, cell = cell))
}

for (sigma_f in c(0, 1.0)) {
  s <- sim(1L, sigma_f)
  f <- tobs(~ x + icar(graph = s$adj, group_var = "cell"),
            detection = ~ 1, data = s$data, family = occu(), y = s$y,
            method = "nested_laplace",
            control = list(verbose = FALSE, progress = FALSE))
  of <- f$nested_laplace$occ_fit
  tg <- of$theta_grid
  w  <- of$weights; w <- w / sum(w)
  cat("\n=== sigma_f =", sigma_f, "===\n")
  cat("theta_names:", paste(of$theta_names, collapse = ","), "\n")
  cat("grid class:", class(tg), " dim:", paste(dim(tg), collapse = "x"), "\n")
  tgv <- if (is.matrix(tg)) tg[, 1] else as.numeric(tg)
  cat("axis values:", paste(sprintf("%.4f", tgv), collapse = ", "), "\n")
  cat("implied sigma=1/sqrt(tau):", paste(sprintf("%.4f", 1 / sqrt(tgv)), collapse = ", "), "\n")
  cat("weights    :", paste(sprintf("%.4f", w), collapse = ", "), "\n")
  cat("peak cell  :", which.max(w), " tau =", sprintf("%.4f", tgv[which.max(w)]), "\n")
  cat("reported field sd:", sprintf("%.4f", stats::sd(as.numeric(f$spatial_field))), "\n")
  nm <- grep("sigma", names(f$means), value = TRUE)
  if (length(nm)) for (n in nm) cat("means[", n, "] =", sprintf("%.4f", f$means[[n]]), "\n")
}
