# Probe: int_occu() + icar() under nested_laplace.
#
# The finding note recorded int_occu as "inherits the M fix, untested". It did
# not inherit it: build_integrated_callbacks' m_step_encode had no latent_prior
# branch (M = 1000 whenever spatial_occ was NULL, i.e. exactly the nested case),
# .tobs_laplace_nested never passed latent_prior to it (so the field-aware
# E-step was dead code), and use_louis excluded "integrated" (so state SEs were
# NA). Measures all three post-fix.
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

sim <- function(seed, field, n_cells = 40L, reps = 6L, J = c(4L, 3L),
                sigma_f = 1.0, b_x = 0.5) {
  set.seed(seed)
  adj <- chain_adj(n_cells)
  f0 <- if (field) smooth_field(n_cells, sigma_f, phase = 0.7) else rep(0, n_cells)
  n_sites <- n_cells * reps
  cell <- rep(seq_len(n_cells), each = reps)
  x <- as.numeric(scale(stats::rnorm(n_sites)))
  z <- stats::rbinom(n_sites, 1L,
                     plogis(stats::qlogis(0.35) + b_x * x + f0[cell]))
  y <- lapply(J, function(j) {
    m <- matrix(0L, n_sites, j)
    for (i in seq_len(n_sites)) m[i, ] <- stats::rbinom(j, 1L, z[i] * plogis(0.4))
    m
  })
  names(y) <- c("s1", "s2")
  list(adj = adj, y = y, f0 = f0, data = data.frame(x = x, cell = cell))
}

fit_one <- function(s) {
  tobs(~ x + icar(graph = s$adj, group_var = "cell"),
       detection = ~ 1, data = s$data, family = int_occu(), y = s$y,
       method = "nested_laplace",
       control = list(verbose = FALSE, progress = FALSE))
}

n_seed <- 8L
for (field in c(FALSE, TRUE)) {
  est <- se <- fsd <- peak <- ncell <- rep(NA_real_, n_seed)
  for (i in seq_len(n_seed)) {
    s <- sim(i, field)
    f <- try(fit_one(s), silent = TRUE)
    if (inherits(f, "try-error")) { cat("seed", i, "ERROR:", conditionMessage(attr(f, "condition")), "\n"); next }
    est[i] <- f$means[["psi_x"]]
    se[i]  <- if ("psi_x" %in% names(f$sds)) f$sds[["psi_x"]] else NA_real_
    fsd[i] <- stats::sd(as.numeric(f$spatial_field))
    w <- f$nested_laplace$occ_fit$weights
    if (!is.null(w)) {
      peak[i]  <- which.max(w)
      ncell[i] <- length(w)
    }
  }
  cat("\n=== field =", field, "(truth fsd =", if (field) 1.0 else 0, ") ===\n")
  cat("slope   median", sprintf("%.4f", stats::median(est, na.rm = TRUE)),
      " truth 0.5\n")
  cat("se      finite ", sum(is.finite(se)), "/", n_seed, "\n")
  cat("field sd median", sprintf("%.4f", stats::median(fsd, na.rm = TRUE)), "\n")
  cat("peak grid cell ", paste(peak, collapse = ","), " of ",
      paste(unique(ncell[!is.na(ncell)]), collapse = ","), "\n")
  if (all(is.finite(se))) {
    hit <- abs(est - 0.5) <= 1.96 * se
    cat("coverage      ", sprintf("%.2f", mean(hit, na.rm = TRUE)), "\n")
  }
}
