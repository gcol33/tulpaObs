# Repro: areal-spatial N-mixture path fails under devtools::load_all() because
# R/nmix_laplace_spatial.R calls the tulpa-internal `.nl_normalise_weights_safe`
# unqualified. It is defined in tulpa's namespace, not tulpaObs's, and is not
# imported, so it does not resolve from inside the tulpaObs namespace.
#
# Affected: tobs(~ x + icar(graph=adj) | bym2(...) | car_proper(...),
#                family = abun(), method = "nested_laplace")
# Call sites: R/nmix_laplace_spatial.R lines 167, 360, 595.
# Fix: qualify as tulpa:::.nl_normalise_weights_safe(...) at all three sites, or
#      add it to the tulpaObs imports (it is a tulpa-internal, so the qualified
#      call or a tiny re-export in tulpa is the clean route).
#
#   "/c/Program Files/R/R-4.6.0/bin/Rscript.exe" dev_notes/repro_nmix_spatial_unqualified_helper.R

suppressMessages(devtools::load_all("."))

cat("tulpa  has .nl_normalise_weights_safe:",
    exists(".nl_normalise_weights_safe", where = asNamespace("tulpa")), "\n")
cat("tulpaObs has .nl_normalise_weights_safe:",
    exists(".nl_normalise_weights_safe", where = asNamespace("tulpaObs")), "\n")

grid_adj <- function(side) {
  ng <- side * side
  co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i != j && abs(co$x[i] - co$x[j]) + abs(co$y[i] - co$y[j]) == 1L)
      adj[i, j] <- 1L
  adj
}
adj <- grid_adj(5)
set.seed(1)
ng     <- nrow(adj)
x_ab   <- rnorm(ng)
x_det  <- rnorm(ng)
lambda <- exp(log(5) + 0.5 * x_ab)
p      <- plogis(0.3 + 0.4 * x_det)
N      <- rpois(ng, lambda)
y      <- matrix(NA_integer_, ng, 5)
for (i in seq_len(ng)) y[i, ] <- rbinom(5, N[i], p[i])
d <- data.frame(abund_cov1 = x_ab, det_cov1 = x_det)

res <- tryCatch(
  tobs(~ abund_cov1 + icar(graph = adj), data = d, family = abun(),
       detection = ~ det_cov1, y = y, method = "nested_laplace",
       control = list(verbose = FALSE)),
  error = function(e) conditionMessage(e)
)
cat("result:", if (is.character(res)) paste("ERROR:", res) else "fit OK", "\n")
