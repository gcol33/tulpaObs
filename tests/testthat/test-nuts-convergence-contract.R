# =============================================================================
# test-nuts-convergence-contract.R -- every family that advertises "nuts" must
# return a per-parameter convergence record.
#
# `n.chains` is a documented control on every NUTS path, so every NUTS path owes
# the diagnostic that makes a chain count worth setting. The contract is the one
# `summary.tobs_fit` and `print.tobs_fit` read:
#
#   fit$convergence$parameter / $rhat / $ess_bulk / $ess_tail
#
# and the names in `$parameter` have to be the names `summary()` puts on its
# rows, or the record is present but unreadable. Both halves are asserted here:
# the record exists and is finite, AND summary() resolves an Rhat for every row
# it prints. Sampler settings are the smallest that still mix -- this file tests
# plumbing, not recovery, so a wide interval is fine.
# =============================================================================

.nconv_ctl <- function(...) {
  utils::modifyList(
    list(n.iter = 150L, n.warmup = 150L, n.chains = 2L, seed = 1L,
         verbose = FALSE, progress = FALSE),
    list(...))
}

# One tiny fit per family that lists "nuts" in .tobs_family_methods.
.nconv_cases <- list(
  occu = function() {
    set.seed(3); N <- 60L; J <- 4L
    z <- stats::rbinom(N, 1, 0.6)
    y <- matrix(stats::rbinom(N * J, 1, rep(z, times = J) * 0.45), N, J)
    tobs(~ 1, data = data.frame(x = stats::rnorm(N)), y = y, detection = ~ 1,
         family = occu(), method = "nuts", control = .nconv_ctl())
  },
  dyn_occu = function() {
    set.seed(4); n <- 60L; Tn <- 3L; J <- 3L
    z <- matrix(NA_integer_, n, Tn); z[, 1] <- stats::rbinom(n, 1, 0.5)
    for (t in 2:Tn) z[, t] <- ifelse(z[, t - 1] == 1,
                                     stats::rbinom(n, 1, 0.8),
                                     stats::rbinom(n, 1, 0.3))
    y <- array(0L, c(n, J, Tn))
    for (i in seq_len(n)) for (t in seq_len(Tn))
      y[i, , t] <- if (z[i, t]) stats::rbinom(J, 1, 0.5) else 0L
    tobs(~ 1, data = data.frame(idx = seq_len(n)), family = dyn_occu(), y = y,
         detection = ~ 1, colonization = ~ 1, extinction = ~ 1,
         method = "nuts", control = .nconv_ctl())
  },
  int_occu = function() {
    sim <- simulate_int_occu(N_total = 150, n_data = 1L, J = 4L,
                             beta_occ = c(0, 0.4), beta_det = list(c(0, -0.3)),
                             seed = 6)
    tobs(~ x, data = sim$data, family = int_occu(), detection = ~ 1, y = sim$y,
         method = "nuts", control = .nconv_ctl())
  },
  ms_occu = function() {
    sim <- simulate_ms_occu(N = 40, J = 3, n_species = 4,
                            beta_comm_mean = c(0, 0.6),
                            beta_comm_sd = c(0.6, 0.3),
                            alpha_comm_mean = c(0.2), alpha_comm_sd = c(0.5),
                            seed = 7)
    tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1, y = sim$y,
         species = paste0("sp", seq_len(4)), method = "nuts",
         control = .nconv_ctl(n.iter = 80L, n.warmup = 80L))
  },
  ms_dyn_occu = function() {
    sim <- simulate_ms_dyn_occu(N = 45, J = 3, n_species = 5, n_seasons = 3,
                                beta_comm_mean = c(0.3), beta_comm_sd = c(0.7),
                                gamma = 0.2, epsilon = 0.1, seed = 12)
    tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
         y = sim$y, species = paste0("sp", seq_len(5)), method = "nuts",
         control = .nconv_ctl())
  },
  ms_int_occu = function() {
    sim <- simulate_ms_int_occu(N = 60, J = c(3, 3), n_species = 5, n_data = 2,
                                seed = 23)
    tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
         y = sim$y, species = paste0("sp", seq_len(5)), method = "nuts",
         control = .nconv_ctl())
  },
  jsdm = function() {
    set.seed(9); N <- 80L; S <- 5L
    d <- data.frame(x = stats::rnorm(N))
    X <- stats::model.matrix(~ x, d)
    bs <- vapply(1:2, function(j)
      stats::rnorm(S, c(0.2, 0.8)[j], c(0.4, 0.3)[j]), numeric(S))
    y <- matrix(stats::rbinom(N * S, 1, stats::plogis(X %*% t(bs))), N, S,
                dimnames = list(NULL, paste0("sp", seq_len(S))))
    tobs(~ x, data = d, family = jsdm(), y = y, species = colnames(y),
         method = "nuts", control = .nconv_ctl())
  },
  ms_count = function() {
    sim <- simulate_ms_count(N = 40, n_species = 4, beta_comm_mean = c(1, 0.5),
                             response = "poisson", seed = 5)
    tobs(~ x, data = sim$data, family = ms_count(), y = sim$y,
         species = colnames(sim$y), method = "nuts",
         control = .nconv_ctl(n.iter = 80L, n.warmup = 80L))
  },
  abun = function() {
    sim <- simulate_abun(N = 50, J = 3, n_abund_covs = 1, n_det_covs = 1,
                         seed = 11)
    tobs(~ abund_cov1, data = sim$data, family = abun(), detection = ~ det_cov1,
         y = sim$y, method = "nuts", control = .nconv_ctl())
  },
  ms_abun = function() {
    sim <- simulate_ms_abun(n_species = 3, N = 25, J = 3, n_abund_covs = 1,
                            n_det_covs = 1, mu_lambda = c(log(3), 0.4),
                            mu_p = c(0.3, -0.3), sd_lambda = 0.5, sd_p = 0.4,
                            seed = 7)
    tobs(~ abund_cov1, data = sim$data, family = ms_abun(K_max = 60L),
         detection = ~ det_cov1, y = sim$y, species = sim$species,
         method = "nuts", control = .nconv_ctl(n.iter = 80L, n.warmup = 80L))
  },
  removal = function() {
    sim <- simulate_removal(N = 50, K = 4, n_abund_covs = 1, n_det_covs = 1,
                            beta_lambda = c(log(6), 0.4), beta_p = c(0.3, -0.3),
                            seed = 12)
    tobs(~ abund_cov1, data = sim$data, family = removal(K_max = 40L),
         detection = ~ det_cov1, y = sim$y, method = "nuts",
         control = .nconv_ctl())
  },
  distance = function() {
    cuts <- seq(0, 1, length.out = 5)
    sim <- simulate_distance(N = 60, cutpoints = cuts, key = "halfnorm",
                             transect = "line", n_abund_covs = 1,
                             n_sigma_covs = 1, beta_lambda = c(log(20), 0.3),
                             beta_sigma = c(log(0.45), 0.2), seed = 13)
    tobs(~ abund_cov1, data = sim$data,
         family = distance(key = "halfnorm", transect = "line",
                           cutpoints = sim$cutpoints),
         detection = ~ sigma_cov1, y = sim$y, method = "nuts",
         control = .nconv_ctl())
  },
  fp_occu = function() {
    sim <- simulate_fp_occu(N = 120, J = 5, n_occ_covs = 1,
                            beta_psi = c(stats::qlogis(0.5), 0.6), p11 = 0.6,
                            p10 = 0.05, b = 0.5, seed = 15)
    tobs(~ occ_cov1, data = sim$data, family = fp_occu(), detection = ~ 1,
         y = sim$y, method = "nuts", control = .nconv_ctl())
  },
  dyn_abun = function() {
    # The forward recursion is cubic in K_max, so the truncation (not the site
    # count) sets this fixture's cost.
    sim <- simulate_dyn_abun(N = 25, T = 2, J = 3, n_abund_covs = 1,
                             beta_lambda = c(log(3), 0.4), p = 0.5, omega = 0.6,
                             gamma = 1.2, seed = 14)
    tobs(~ abund_cov1, data = sim$data, family = dyn_abun(K_max = 12),
         detection = ~ 1, y = sim$y, method = "nuts",
         control = .nconv_ctl(n.iter = 80L, n.warmup = 80L))
  },
  cover = function() {
    sim <- simulate_cover(N = 120, beta_occ = c(-0.4, 0.8),
                          beta_pos = c(-1, 0.3), sigma_pos = 0.4, seed = 1)
    tobs(~ x, data = sim$data, family = cover("lognormal"), y = sim$y,
         method = "nuts", control = .nconv_ctl())
  },
  occu_cover = function() {
    N <- 60L; J <- 3L
    sim <- simulate_occu_cover(N = N, J = J, n_occ_covs = 1L, n_det_covs = 1L,
                               n_pos_covs = 1L, sigma_pos = 0.4,
                               positive = "lognormal", seed = 2)
    long <- data.frame(site_id = rep(seq_len(N), each = J),
                       visit = rep(seq_len(J), times = N),
                       y = as.vector(t(sim$y)),
                       det_cov1 = sim$visit_data$det_cov1,
                       pos_cov1 = sim$visit_data$pos_cov1)
    od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                    det.covs = c("det_cov1", "pos_cov1"))
    y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
    tobs(~ occ_cov1, data = cbind(data.frame(site_id = seq_len(N)), sim$data),
         family = occu_cover("lognormal"), detection = ~ det_cov1,
         positive = ~ pos_cov1, y = od$y, y_pos = y_pos,
         visits = od$det.covs, method = "nuts", control = .nconv_ctl())
  },
  occu_multiscale_cover = function() {
    sim <- simulate_occu_multiscale_cover(
      n_cells = 30L, plots_per_cell = 3L, visits_per_plot = 3L,
      positive = "lognormal", phi = 0.35, sigma = 0, alpha = 0, seed = 909L)
    suppressWarnings(tobs(
      formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
      data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
      detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
      y = sim$y, y_pos = sim$y_pos, method = "nuts", control = .nconv_ctl()))
  },
  ms_occu_cover = function() {
    N <- 60L; J <- 3L; S <- 5L
    sim <- simulate_ms_occu_cover(
      n_species = S, N = N, J = J,
      mu_occ = c(stats::qlogis(0.45), 0.7), mu_p = c(0.2, -0.4),
      mu_pos = c(log(0.12), 0.5), sd_occ = 0.5, sd_p = 0.4, sd_pos = 0.4,
      positive = "lognormal", sigma_pos = 0.4, seed = 21)
    long <- data.frame(site_id = rep(seq_len(N), each = J),
                       visit = rep(seq_len(J), times = N), yy = 0L,
                       det_cov1 = sim$visit_data$det_cov1,
                       pos_cov1 = sim$visit_data$pos_cov1)
    vis <- tobs_data(long, y = "yy", site = "site_id", visit = "visit",
                     det.covs = c("det_cov1", "pos_cov1"))$det.covs
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1, y = sim$y,
         y_pos = sim$y_pos, visits = vis, species = sim$species,
         method = "nuts", control = .nconv_ctl())
  })


test_that("every family advertising nuts has a convergence-contract case", {
  advertised <- names(Filter(function(m) "nuts" %in% m,
                             tulpaObs:::.tobs_family_methods))
  expect_setequal(advertised, names(.nconv_cases))
})


test_that("every NUTS family reports per-parameter Rhat / ESS", {
  skip_on_cran()
  skip_if_fast()
  for (nm in names(.nconv_cases)) {
    fit <- .nconv_cases[[nm]]()
    cv  <- fit$convergence

    expect_false(is.null(cv$rhat), info = nm)
    expect_false(is.null(cv$ess_bulk), info = nm)
    expect_false(is.null(cv$ess_tail), info = nm)
    expect_true(all(is.finite(cv$rhat)), info = nm)
    expect_true(all(cv$ess_bulk > 0), info = nm)
    expect_true(all(cv$ess_tail > 0), info = nm)
    # Tail ESS is the 5% / 95% indicator ESS, a different quantity from the bulk
    # ESS -- not the bulk value re-reported under the tail name.
    expect_false(identical(unname(cv$ess_bulk), unname(cv$ess_tail)), info = nm)

    # The record has to be readable: summary() matches its rows against
    # cv$parameter, so every printed row must resolve an Rhat.
    s <- summary(fit)
    expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(s)), info = nm)
    expect_true(all(is.finite(s$rhat)), info = nm)

    # `converged` on a sampled fit reports chain mixing, not the warm-start
    # optimiser: TRUE only when every reported split-Rhat is below 1.01.
    expect_identical(cv$converged, all(cv$rhat < 1.01), info = nm)
    expect_equal(fit$max_rhat, max(cv$rhat), info = nm)
  }
})
