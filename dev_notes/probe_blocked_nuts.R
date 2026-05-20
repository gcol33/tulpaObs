# probe_blocked_nuts.R
# Smoke-test the four components flagged in CLAUDE.md as "wired but blocked by
# upstream tulpa NUTS bugs":
#   - tobs_temporal
#   - tobs_re (multi-term)
#   - tobs_svc
#   - tobs_latent
#
# Goal: find out whether each still crashes, errors cleanly, or works.
# Tiny N / few warmup so we fail fast.

options(error = NULL)
devtools::load_all(quiet = TRUE)

set.seed(1)
N <- 40; J <- 3
sim <- simulate_occu(N = N, J = J, seed = 1)
coords <- matrix(runif(N * 2), N, 2)
sim$data$grp  <- sample.int(5, N, replace = TRUE)
sim$data$time <- sample.int(4, N, replace = TRUE)

ctl_base <- list(iter = 50, warmup = 25, seed = 1, verbose = FALSE)

run <- function(label, expr) {
  cat(sprintf("\n=== %s ===\n", label))
  t0 <- Sys.time()
  res <- tryCatch(eval(expr),
                  error    = function(e) list(status = "error",   msg = conditionMessage(e)),
                  warning  = function(w) list(status = "warning", msg = conditionMessage(w)))
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(res, "tobs_fit")) {
    cat(sprintf("OK (%.2fs) — class: %s\n", dt, paste(class(res), collapse = "/")))
  } else if (is.list(res) && !is.null(res$status)) {
    cat(sprintf("%s (%.2fs): %s\n", toupper(res$status), dt, res$msg))
  } else {
    cat(sprintf("Returned non-fit (%.2fs): class %s\n", dt, paste(class(res), collapse = "/")))
  }
  invisible(res)
}

# (1) temporal AR1
r1 <- run("temporal AR1 / NUTS", quote(
  tobs(~ occ_cov1, data = sim$data, family = occu(),
       detection = ~ det_cov1, y = sim$y,
       temporal = tobs_temporal(type = "ar1", time = "time"),
       engine = "nuts", control = ctl_base)
))

# (2) multi-term RE (two REs in a list)
r2 <- run("multi-term RE / NUTS", quote(
  tobs(~ occ_cov1, data = sim$data, family = occu(),
       detection = ~ det_cov1, y = sim$y,
       re = list(
         tobs_re(group = "grp",  type = "intercept"),
         tobs_re(group = "time", type = "intercept")
       ),
       engine = "nuts", control = ctl_base)
))

# (3) SVC — passed via control$svc (.tobs_fit_model arg, not user-facing in tobs())
r3 <- run("SVC / NUTS", quote(
  tobs(~ occ_cov1, data = sim$data, family = occu(),
       detection = ~ det_cov1, y = sim$y,
       engine = "nuts",
       control = c(ctl_base, list(svc = tobs_svc(indices = 1L, coords = coords, nn = 8))))
))

# (4) latent factors via control$latent on community model
n_sp <- 4
ms <- simulate_ms_occu(N = N, J = J, n_species = n_sp, seed = 1)
sp_names <- paste0("sp", seq_len(n_sp))
r4 <- run("latent factors / NUTS (ms_occu)", quote(
  tobs(~ x, data = ms$data, family = ms_occu(),
       detection = ~ 1, y = ms$y, species = sp_names,
       engine = "nuts",
       control = c(ctl_base, list(latent = tobs_latent(n_factors = 2))))
))

cat("\n--- summary ---\n")
for (nm in c("r1", "r2", "r3", "r4")) {
  v <- get(nm)
  if (inherits(v, "tobs_fit")) cat(sprintf("  %s: OK\n", nm))
  else if (is.list(v) && !is.null(v$status))
    cat(sprintf("  %s: %s — %s\n", nm, v$status, substr(v$msg, 1, 120)))
  else cat(sprintf("  %s: unexpected (%s)\n", nm, paste(class(v), collapse = "/")))
}
