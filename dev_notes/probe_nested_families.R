# Probe: does method = "nested_laplace" now run for integrated / community /
# dynamic occupancy with a spatial (bym2) latent block? Checks the pipeline
# returns a sensibly-shaped tobs_fit with the multi-block prior attached and a
# recovered latent field of the right length (one value per site).
suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE))

chain_adj <- function(n) {
  adj <- matrix(0, n, n)
  for (i in seq_len(n - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
  adj
}

ok <- function(label, expr) {
  res <- tryCatch(expr, error = function(e) e)
  if (inherits(res, "error")) {
    cat(sprintf("[FAIL] %s: %s\n", label, conditionMessage(res)))
    return(invisible(NULL))
  }
  cat(sprintf("[ OK ] %s\n", label))
  invisible(res)
}

# ---- integrated ----------------------------------------------------------
set.seed(1)
n <- 40; adj <- chain_adj(n)
elev <- rnorm(n)
d_int <- data.frame(elev = elev)
z <- rbinom(n, 1, plogis(0.2 + 0.5 * elev))
J <- c(3, 4)
y_int <- lapply(J, function(j) {
  y <- matrix(0L, n, j)
  for (i in seq_len(n)) if (z[i]) y[i, ] <- rbinom(j, 1, 0.4)
  y
})
names(y_int) <- c("s1", "s2")
fit_int <- ok("int_occu + bym2", tobs(
  ~ elev + bym2(graph = adj), data = d_int, family = int_occu(),
  detection = ~ 1, y = y_int, method = "nested_laplace",
  control = list(max.iter = 8L, verbose = FALSE)))
if (!is.null(fit_int)) {
  cat("   field length:", length(fit_int$spatial_field),
      " (expect", n, ")\n")
  cat("   prior type:", fit_int$nested_laplace$multi_prior$type, "\n")
}

# ---- community -----------------------------------------------------------
set.seed(2)
n_sites <- 30; n_species <- 3; adj_c <- chain_adj(n_sites)
d_ms <- data.frame(x = rnorm(n_sites))
y_ms <- list()
for (s in seq_len(n_species)) {
  zz <- rbinom(n_sites, 1, plogis(rnorm(1, 0, 0.3)))
  ys <- matrix(0L, n_sites, 3)
  for (i in seq_len(n_sites)) if (zz[i]) ys[i, ] <- rbinom(3, 1, 0.4)
  y_ms[[paste0("sp", s)]] <- ys
}
fit_ms <- ok("ms_occu + bym2", tobs(
  ~ 1 + bym2(graph = adj_c), data = d_ms, family = ms_occu(),
  detection = ~ 1, y = y_ms, species = TRUE, method = "nested_laplace",
  control = list(max.iter = 8L, verbose = FALSE)))
if (!is.null(fit_ms)) {
  cat("   field length:", length(fit_ms$spatial_field),
      " (expect", n_sites, ")\n")
}

# ---- dynamic -------------------------------------------------------------
set.seed(3)
n_sites <- 40; n_seasons <- 3; J <- 3; adj_d <- chain_adj(n_sites)
d_dyn <- data.frame(elev = rnorm(n_sites))
psi1 <- plogis(0.2 + 0.4 * d_dyn$elev)
zmat <- matrix(0L, n_sites, n_seasons)
zmat[, 1] <- rbinom(n_sites, 1, psi1)
for (t in 2:n_seasons) {
  surv <- zmat[, t - 1] * (1 - rbinom(n_sites, 1, plogis(-1)))
  col  <- (1 - zmat[, t - 1]) * rbinom(n_sites, 1, plogis(-1.5))
  zmat[, t] <- surv + col
}
y_dyn <- array(NA_integer_, dim = c(n_sites, J, n_seasons))
for (i in seq_len(n_sites)) for (t in seq_len(n_seasons))
  y_dyn[i, , t] <- if (zmat[i, t]) rbinom(J, 1, 0.4) else 0L
fit_dyn <- ok("dyn_occu + bym2", tobs(
  ~ elev + bym2(graph = adj_d), data = d_dyn, family = dyn_occu(),
  detection = ~ 1, y = y_dyn, col_formula = ~ 1, ext_formula = ~ 1,
  method = "nested_laplace", control = list(max.iter = 8L, verbose = FALSE)))
if (!is.null(fit_dyn)) {
  cat("   field length:", length(fit_dyn$spatial_field),
      " (expect", n_sites, ")\n")
}

cat("\nDONE\n")
