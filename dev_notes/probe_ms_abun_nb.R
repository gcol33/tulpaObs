# probe_ms_abun_nb.R - verify NB community wiring after lifting the
# .tobs_fit_ms_nmix Poisson-only gate.

suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})

set.seed(31L)
sim <- simulate_ms_abun(
  n_species   = 8L,
  N           = 60L,
  J           = 4L,
  n_abund_covs = 1L,
  n_det_covs  = 1L,
  mu_lambda   = c(log(3), 0.4),
  mu_p        = c(0.3, -0.3),
  sd_lambda   = 0.5,
  sd_p        = 0.4,
  mixture     = "negbin",
  size        = 4,
  seed        = 31L
)

cat("== ms_abun(mixture = \"negbin\") smoke ==\n")
cat(sprintf("  truth: mu_lambda = (%.2f, %.2f), mu_p = (%.2f, %.2f), r = %g\n",
            sim$truth$mu_lambda[1L], sim$truth$mu_lambda[2L],
            sim$truth$mu_p[1L],      sim$truth$mu_p[2L],
            sim$truth$size))

t0 <- proc.time()[["elapsed"]]
fit <- tryCatch(
  tobs(
    formula   = ~ abund_cov1,
    detection = ~ det_cov1,
    data      = sim$data,
    family    = ms_abun(mixture = "negbin"),
    y         = sim$y,
    species   = sim$species,
    control   = list(verbose = FALSE)
  ),
  error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
)
dt <- proc.time()[["elapsed"]] - t0
if (is.null(fit)) quit(save = "no", status = 1L)
cat(sprintf("  fit OK in %.2fs\n", dt))
cat(sprintf("  class:    %s\n", paste(class(fit), collapse = "/")))
cat(sprintf("  mixture:  %s\n", fit$mixture))
cat(sprintf("  estimated r: %s\n",
            if (!is.null(fit$ms_dispersion)) signif(fit$ms_dispersion$r, 4)
            else "(missing)"))
cat("\n  community means (estimate vs truth):\n")
m <- fit$means; s <- fit$sds
# Drop log_r from the comparison frame; it's NOT a community mean.
keep <- setdiff(names(m), "log_r")
m_k  <- m[keep]; s_k <- s[keep]
truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
names(truth) <- keep
print(round(rbind(estimate = m_k, truth = truth, sd = s_k,
                  z = (m_k - truth) / s_k), 3))
cat("\n  community SDs:\n")
cm <- fit$ms_community
print(round(rbind(estimate = c(cm$sd_lambda, cm$sd_p),
                  truth   = c(sim$truth$sd_lambda, sim$truth$sd_p)), 3))
