`%||%` <- function(a, b) if (is.null(a)) b else a
# gcol33/tulpaObs#415: the joint engine's `iid` block is unit-precision with the
# SD on d_fac (tulpa src/nested_laplace_joint_multi.cpp, `type == "iid"`:
# eta += sigma * u, u ~ N(0, 1)), so the latent columns tulpaObs reads as BLUPs
# are u, not the effect. This rescales the RE draws by each draw's cell sigma
# and re-scores WAIC, and computes the RE posterior SD of the centred effect
# from the same draws. Reads fits saved by fixture415.R (W415 dir).
suppressPackageStartupMessages({library(tulpa); library(tulpaObs)})
w <- Sys.getenv("W415", ".")
fits <- do.call(c, lapply(list.files(w, "^fits415.*[.]rds$", full.names = TRUE), readRDS))
fits <- fits[!duplicated(names(fits))]
orig <- tulpaObs:::.tobs_joint_draws
scaled <- function(object, n, ...) {
  b <- orig(object, n = n, ...)
  if (length(b$re)) {
    tg <- object$joint_fit$theta_grid
    re_axes <- grep("^b[0-9]+\\.sigma$", colnames(tg), value = TRUE)
    # the RE blocks trail the field blocks; one sigma axis per RE block, in order
    re_axes <- tail(re_axes, length(b$re))
    for (i in seq_along(b$re)) {
      sig <- tg[b$cells, re_axes[i]]
      b$re[[i]]$draws <- b$re[[i]]$draws * sig
    }
  }
  b
}
for (nm in c("fix", "rail", "span", "wide")) {
  f <- fits[[nm]]; if (is.null(f)) next
  set.seed(7)
  w0 <- suppressWarnings(waic(f))$estimates
  assignInNamespace(".tobs_joint_draws", scaled, "tulpaObs")
  set.seed(7)
  w1 <- suppressWarnings(waic(f))$estimates
  b <- tulpaObs:::.tobs_joint_draws(f, n = 4000)
  assignInNamespace(".tobs_joint_draws", orig, "tulpaObs")
  cat(sprintf("%-5s p_waic as shipped %8.2f -> sigma-scaled %7.2f | elpd %9.2f -> %9.2f\n",
              nm, w0["p_waic", 1], w1["p_waic", 1], w0["elpd_waic", 1], w1["elpd_waic", 1]))
  if (length(b$re)) {
    d <- b$re[[1]]$draws
    dc <- d - rowMeans(d)
    cat("   centred effect mean:", sprintf("%6.2f", colMeans(dc)), "\n")
    cat("   centred effect SD:  ", sprintf("%6.3f", apply(dc, 2, sd)), "\n")
    cat("   ranef() std.error:  ", sprintf("%6.0e", ranef(f)$std.error), "\n")
    lev <- b$b$det[, 1] + d
    cat("   intercept + b_g:    ", sprintf("%6.2f", colMeans(lev)), "\n")
  } else {
    cf <- f$means
    lev <- cf["p_(Intercept)"] + c(0, cf[grep("^p_habitat", names(cf))])
    cat("   fixed level means:  ", sprintf("%6.2f", lev), "\n")
  }
}
cat("truth levels (0.2 + hab):", sprintf("%6.2f", 0.2 + c(-3.9, -2.0, -1.0, 0, 0.5, 1.0, 1.5, 2.1)), "\n")
