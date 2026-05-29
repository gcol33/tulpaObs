# probe_ms_abun_profile.R - profile a single ms_abun() fit and break down
# where the time goes (warm-start vs EM mode-find vs covariance M-step).

suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})

# Two fixtures: tiny and the test's coverage-block size, both Poisson so the
# fast EM path is exercised.
fixtures <- list(
  tiny     = list(n_species = 5L,  N = 20L, J = 3L, seed = 41L),
  coverage = list(n_species = 12L, N = 60L, J = 4L, seed = 42L)
)

for (nm in names(fixtures)) {
  cfg <- fixtures[[nm]]
  cat(sprintf("\n=== %s: S=%d, N=%d, J=%d ===\n",
              nm, cfg$n_species, cfg$N, cfg$J))
  sim <- simulate_ms_abun(
    n_species   = cfg$n_species, N = cfg$N, J = cfg$J,
    n_abund_covs = 1L, n_det_covs = 1L,
    mu_lambda   = c(log(3), 0.4),
    mu_p        = c(0.3, -0.3),
    sd_lambda   = 0.5, sd_p = 0.4,
    mixture     = "poisson",
    seed        = cfg$seed
  )

  pf <- tempfile(fileext = ".out")
  Rprof(pf, interval = 0.02, line.profiling = FALSE)
  t0 <- proc.time()[["elapsed"]]
  fit <- tobs(
    formula   = ~ abund_cov1,
    detection = ~ det_cov1,
    data      = sim$data,
    family    = ms_abun(mixture = "poisson"),
    y         = sim$y,
    species   = sim$species,
    control   = list(verbose = FALSE)
  )
  dt <- proc.time()[["elapsed"]] - t0
  Rprof(NULL)
  cat(sprintf("  total fit: %.2fs\n", dt))
  prof <- summaryRprof(pf)
  cat("\n  by.self top 12 (self.pct >= 1):\n")
  bs <- prof$by.self
  bs <- bs[bs$self.pct >= 1, ]
  print(utils::head(bs, 12L))
  cat("\n  by.total top 10 (total.pct >= 5):\n")
  bt <- prof$by.total
  bt <- bt[bt$total.pct >= 5, ]
  print(utils::head(bt, 10L))
}
