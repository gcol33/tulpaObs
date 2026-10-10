# gcol33/tulpaObs#414: temporal() in the occu_cover() occurrence formula on the
# joint nested-Laplace route. 10 cells on a chain graph, 5 intervals, 6 plots per
# site; runs in seconds. Before the fix the rw1 fit carried no time structure on
# psi (hyper axes `sigma phi_pos`, no psi time coefficient at all).
suppressPackageStartupMessages({library(tulpa); library(tulpaObs)})
cat("tulpa", as.character(packageVersion("tulpa")), "\n")
cat("tulpaObs", as.character(packageVersion("tulpaObs")),
    packageDescription("tulpaObs")$RemoteSha %||% "", "\n")
set.seed(2)
nc <- 10L; n_int <- 5L
adj <- matrix(0L, nc, nc)
for (i in 1:(nc - 1)) { adj[i, i + 1] <- 1L; adj[i + 1, i] <- 1L }
rows <- list()
for (cell in seq_len(nc)) for (ti in seq_len(n_int)) {
  J <- 6L; tsc <- ti - 3
  psi <- plogis(-0.3 + 0.6 * sin(ti) + 0.15 * (cell - nc / 2))
  z <- rbinom(1, 1, psi)
  x1 <- rnorm(J); p <- plogis(0.2 + 0.4 * x1)
  occur <- if (z) rbinom(J, 1, p) else rep(0L, J)
  cover <- ifelse(occur == 1, pmin(pmax(rbeta(J, 2, 4), 1e-4), 1 - 1e-4), 0)
  rows[[length(rows) + 1]] <- data.frame(
    cell_idx = cell, yr.interval = ti, site_key = paste(cell, ti, sep = "_"),
    visit = seq_len(J), time.sc = tsc, x1 = x1, occur = as.integer(occur),
    cover.flat = cover)
}
dd <- do.call(rbind, rows)
ctrl <- list(engine = "joint", n.threads = 1L, n.threads.outer = 1L,
             max.iter = 200L, sigma.grid = c(0.5, 1.5), phi.grid.pos = c(2, 10),
             integration = "grid", verbose = FALSE, progress = FALSE)
go <- function(f) tryCatch({
  fit <- tobs(occurrence = f, data = dd,
              family = occu_cover(response = "beta", cover.aggregate = "none"),
              detection = ~ x1, positive = ~ time.sc, response = "occur",
              y.pos = "cover.flat", site = "site_key", visit = "visit",
              det.covs = c("x1", "time.sc"),
              occ.covs = c("cell_idx", "time.sc", "yr.interval"),
              method = "nested_laplace", control = ctrl)
  cat("OK; coef names:", paste(names(fit$means), collapse = " "), "\n")
  cat("   hyper axes:", paste(colnames(fit$joint_fit$theta_grid), collapse = " "),
      "\n")
  fit
}, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
cat("--- rw1 on yr.interval\n")
f1 <- go(~ temporal(yr.interval, type = "rw1") +
           spatial(~ 1 || cell_idx, graph = adj))
if (!is.null(f1$re$temporal)) {
  cat("   per-interval effect:", round(f1$re$temporal$blup, 3), "\n")
}
cat("--- linear\n")
f0 <- go(~ time.sc + spatial(~ 1 || cell_idx, graph = adj))
