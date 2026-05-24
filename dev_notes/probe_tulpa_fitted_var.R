# Does the rebuilt tulpa engine return fitted_eta_var directly?
#   "/c/Program Files/R/R-4.6.0/bin/Rscript.exe" dev_notes/probe_tulpa_fitted_var.R

library(tulpa)

set.seed(1)
n <- 25                     # 5x5 grid
gx <- gy <- 5
adj <- matrix(0, n, n)
idx_of <- function(i, j) (j - 1) * gx + i
for (i in seq_len(gx)) for (j in seq_len(gy)) {
  a <- idx_of(i, j)
  if (i < gx) { b <- idx_of(i + 1, j); adj[a, b] <- 1; adj[b, a] <- 1 }
  if (j < gy) { b <- idx_of(i, j + 1); adj[a, b] <- 1; adj[b, a] <- 1 }
}
deg <- rowSums(adj)
rp  <- as.integer(c(0, cumsum(deg)))
ci  <- as.integer(unlist(lapply(seq_len(n), function(i) which(adj[i, ] == 1)) ) - 1)

y <- rbinom(n, 1, 0.5)
prior <- list(type = "icar", spatial_idx = seq_len(n), n_spatial_units = n,
              adj_row_ptr = rp, adj_col_idx = ci, n_neighbors = as.integer(deg))

res <- tulpa::tulpa_nested_laplace(
  y = y, n_trials = rep(1L, n), X = matrix(1, n, 1),
  prior = prior)

cat("names(res):\n"); print(names(res))
cat("\nhas fitted_eta:     ", !is.null(res$fitted_eta), "\n")
cat("has fitted_eta_var: ", !is.null(res$fitted_eta_var), "\n")
if (!is.null(res$fitted_eta_var)) {
  cat("dim fitted_eta_var:", paste(dim(res$fitted_eta_var), collapse = " x "), "\n")
  cat("range:", paste(round(range(res$fitted_eta_var), 4), collapse = " .. "), "\n")
  cat("col1 (site1 across grid):", paste(round(res$fitted_eta_var[, 1], 4), collapse = " "), "\n")
}
