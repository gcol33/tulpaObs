# gcol33/tulpaObs#415: occu_cover joint route, detection random intercept
# (1 | habitat). Measures ranef() std.error against the posterior SD of the same
# latent columns drawn from the fit's own grid mixture, and p_waic against the
# fixed-habitat fit, at a railed RE sigma grid and at one spanning the truth.
suppressPackageStartupMessages({library(tulpa); library(tulpaObs)})
cat("tulpa", as.character(packageVersion("tulpa")), "| tulpaObs",
    as.character(packageVersion("tulpaObs")),
    packageDescription("tulpaObs")$RemoteSha %||% "", "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a

set.seed(415)
nc <- 40L; J <- 40L
adj <- matrix(0L, nc, nc); for (i in 1:(nc - 1)) { adj[i, i + 1] <- 1L; adj[i + 1, i] <- 1L }
hab_eff <- c(-3.9, -2.0, -1.0, 0, 0.5, 1.0, 1.5, 2.1)   # truth, logit detection
hab_lev <- paste0("H", 1:8)
rows <- list()
for (cell in seq_len(nc)) {
  psi <- plogis(0.3 + 0.1 * (cell - nc / 2)); z <- rbinom(1, 1, psi)
  x1 <- rnorm(J); h <- sample.int(8L, J, replace = TRUE)
  p <- plogis(0.2 + 0.4 * x1 + hab_eff[h])
  occur <- if (z) rbinom(J, 1, p) else rep(0L, J)
  cover <- ifelse(occur == 1, pmin(pmax(rbeta(J, 2, 4), 1e-4), 1 - 1e-4), 0)
  rows[[cell]] <- data.frame(cell_idx = cell, site_key = as.character(cell),
    visit = seq_len(J), x1 = x1, habitat = factor(hab_lev[h], levels = hab_lev),
    occur = as.integer(occur), cover.flat = cover)
}
dd <- do.call(rbind, rows)
cat("sites", nc, "visits", nrow(dd), "detections", sum(dd$occur), "\n")

base_ctrl <- list(engine = "joint", n.threads = 1L, n.threads.outer = 1L,
                  max.iter = 200L, integration = "grid", verbose = FALSE,
                  progress = FALSE)
fit_one <- function(det, ctrl_extra = list()) {
  tobs(occurrence = ~ 1 + spatial(~ 1 || cell_idx, graph = adj), data = dd,
       family = occu_cover(response = "beta", cover.aggregate = "none"),
       detection = det, positive = ~ 1, response = "occur", y.pos = "cover.flat",
       site = "site_key", visit = "visit", det.covs = c("x1", "habitat"),
       occ.covs = c("cell_idx"), method = "nested_laplace",
       control = c(base_ctrl, ctrl_extra))
}

report <- function(lbl, fit) {
  cat("\n==== ", lbl, "\n")
  w <- tryCatch(waic(fit), error = function(e) { cat("waic error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(w)) {
    est <- if (!is.null(w$estimates)) w$estimates else NULL
    print(est)
  }
  jf <- fit$joint_fit
  tg <- jf$theta_grid
  for (a in colnames(tg)) {
    agg <- tapply(jf$weights, tg[, a], sum)
    cat(sprintf("  axis %-10s ", a), paste(sprintf("%.3g:%.3f", as.numeric(names(agg)), agg), collapse = " "), "\n")
  }
  if (!is.null(fit$re) && length(fit$re)) {
    re <- Filter(function(t) is.list(t) && !is.null(t$latent_idx), fit$re)
    for (t in re) {
      D <- tulpa::tulpa_posterior_draws(jf, idx = t$latent_idx, n = 4000)
      draw_sd <- apply(D, 2, sd)
      # within-cell (conditional) SD alone, from the stored precision at the
      # heaviest cell: diag of the inverse restricted to these columns
      cat("  ranef blup:      ", sprintf("%7.3f", t$blup), "\n")
      cat("  ranef blup_sd:   ", sprintf("%7.1e", t$blup_sd), "\n")
      cat("  draw SD (mixture):", sprintf("%7.3f", draw_sd), "\n")
      cat("  draw mean - mean :", sprintf("%7.3f", colMeans(D) - mean(colMeans(D))), "\n")
    }
    print(ranef(fit))
  }
  cf <- fit$means; sdv <- fit$sds
  sel <- grepl("^p_|sigma_re", names(cf))
  print(data.frame(est = round(cf[sel], 4), sd = round(sdv[sel], 4)))
}

# FITS selects which fits run (fix, rail, span); COMPRESS = 0 turns the
# detection-pattern compression off (option tulpaObs.compress_nodet).
options(tulpaObs.compress_nodet = Sys.getenv("COMPRESS", "1") != "0")
cat("compress_nodet:", getOption("tulpaObs.compress_nodet"), "
")
which_fits <- strsplit(Sys.getenv("FITS", "fix,rail,span"), ",")[[1]]
specs <- list(
  fix  = list(lbl = "fixed habitat", det = ~ x1 + habitat, ctrl = list()),
  rail = list(lbl = "RE, sigma grid c(0.2, 0.6, 1.5)", det = ~ x1 + (1 | habitat),
              ctrl = list(re.sigma.grid.p = c(0.2, 0.6, 1.5))),
  span = list(lbl = "RE, sigma grid c(0.5 .. 4)", det = ~ x1 + (1 | habitat),
              ctrl = list(re.sigma.grid.p = c(0.5, 1, 1.5, 2, 3, 4))),
  wide = list(lbl = "RE, sigma pinned at 20 (near-flat prior)", det = ~ x1 + (1 | habitat),
              ctrl = list(re.sigma.grid.p = 20)))
fits <- list()
for (k in which_fits) {
  fits[[k]] <- fit_one(specs[[k]]$det, specs[[k]]$ctrl)
  report(specs[[k]]$lbl, fits[[k]])
}
saveRDS(fits, file.path(Sys.getenv("W415", "."),
                        sprintf("fits415_c%s.rds", Sys.getenv("COMPRESS", "1"))))
