# =============================================================================
# helper-sbc-registry.R
# -- fixtures and the means accessor the SBC registry tests share
#
# Read by test-sbc-registry.R, which covers the registry contract, and by each
# test-sbc-acceptance-<family>.R, which carries one family's calibration
# measurement.
# =============================================================================

.sbc_reg_ctl <- list(verbose = FALSE, progress = FALSE)

# One fixture per registered non-coupled family, small enough that a refit costs
# a fraction of a second. Each returns the fitted model SBC then rechecks.
.SBC_REG_FIXTURES <- list(
  occu = function(N = 60L) {
    sim <- simulate_occu(N = N, J = 4L, seed = 11L)
    suppressWarnings(tobs(~ occ_cov1, data = sim$data, family = occu(),
                          detection = ~ det_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  count = function(N = 120L) {
    sim <- simulate_count(N = N, beta = c(1, 0.5), seed = 12L)
    suppressWarnings(tobs(~ x, data = sim$data, family = count("poisson"),
                          y = sim$y, method = "laplace",
                          control = .sbc_reg_ctl))
  },
  abun = function(N = 60L) {
    sim <- simulate_abun(N = N, J = 4L, seed = 13L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data, family = abun(),
                          detection = ~ det_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  royle_nichols = function(N = 100L) {
    sim <- simulate_royle_nichols(N = N, J = 5L, seed = 14L)
    suppressWarnings(tobs(~ x, data = sim$data, family = royle_nichols(),
                          detection = ~ x, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  occu_ttd = function(N = 100L) {
    sim <- simulate_occu_ttd(N = N, J = 4L, seed = 15L)
    suppressWarnings(tobs(~ psi_cov1, data = sim$data,
                          family = occu_ttd(surveyLength = sim$Tmax),
                          detection = ~ rate_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  fp_occu = function(N = 120L) {
    sim <- simulate_fp_occu(N = N, J = 5L, seed = 16L)
    suppressWarnings(tobs(~ occ_cov1, data = sim$data, family = fp_occu(),
                          detection = ~ occ_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  removal = function(N = 80L) {
    sim <- simulate_removal(N = N, K = 4L, seed = 17L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data, family = removal(),
                          detection = ~ det_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  distance = function(N = 80L) {
    sim <- simulate_distance(N = N, seed = 18L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data,
                          family = distance(cutpoints = sim$cutpoints,
                                            key = "halfnorm",
                                            transect = "line"),
                          detection = ~ sigma_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  double_observer = function(N = 150L) {
    sim <- simulate_double_observer(N = N, beta_lambda = c(log(8), 0.4),
                                    beta_p1 = c(stats::qlogis(0.5), 0.2),
                                    beta_p2 = c(stats::qlogis(0.45), -0.1),
                                    seed = 19L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data,
                          family = double_observer(),
                          detection = ~ det_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  dyn_occu = function(N = 80L) {
    sim <- simulate_dyn_occu(N = N, J = 4L, n_seasons = 5L,
                             beta_occ = c(0.2, 0.6), beta_det = c(0.4),
                             gamma = 0.25, epsilon = 0.15, seed = 21L)
    suppressWarnings(tobs(~ x, data = sim$data, family = dyn_occu(),
                          detection = ~ 1, colonization = ~ 1,
                          extinction = ~ 1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  int_occu = function(N = 60L) {
    set.seed(22L)
    x_cov <- rnorm(N); det_cov <- rnorm(N)
    z <- rbinom(N, 1, plogis(0.2 + 0.7 * x_cov))
    mk <- function(J, p0) {
      p <- plogis(p0 + 0.4 * det_cov)
      y <- matrix(0L, N, J)
      for (i in seq_len(N)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p[i])
      y
    }
    dat <- data.frame(occ_cov = x_cov, det_cov = det_cov)
    yy <- list(src1 = mk(4L, -0.2), src2 = mk(3L, -0.5))
    suppressWarnings(tobs(~ occ_cov, data = dat, family = int_occu(),
                          detection = ~ det_cov, y = yy,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  dyn_abun = function(N = 80L) {
    sim <- simulate_dyn_abun(N = N, T = 4L, J = 3L,
                             beta_lambda = c(log(5), 0.3), p = 0.5,
                             omega = 0.6, gamma = 1.0, seed = 23L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data, family = dyn_abun(),
                          detection = ~ 1, omega = ~ 1, gamma = ~ 1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  gdistremoval = function(N = 150L) {
    cutp <- c(0, 25, 50, 75, 100)
    sim <- simulate_gdistremoval(N = N, cutpoints = cutp, n_periods = 4L,
                                 beta_lambda = c(log(30), 0.3),
                                 beta_sigma = c(log(18), 0.1),
                                 beta_r = c(stats::qlogis(0.4), -0.2), seed = 11L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data, y = sim$y,
                          y_rem = sim$y_rem,
                          family = gdistremoval(cutpoints = cutp),
                          detection = ~ det_cov1, removal = ~ rem_cov1,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  occu_categorical = function(N = 300L) {
    sim <- simulate_occu_categorical(N = N, seed = 11L)
    suppressWarnings(tobs(~ x, data = sim$data, family = occu_categorical(),
                          y = sim$y, method = "laplace", control = .sbc_reg_ctl))
  },
  distsamp_open = function(N = 100L) {
    cutp <- c(0, 10, 20, 30, 40)
    sim <- simulate_distsamp_open(N = N, cutpoints = cutp, n_seasons = 4L,
                                  beta_lambda = c(log(15), 0.3),
                                  beta_sigma = c(log(15), 0.1),
                                  omega = 0.7, gamma = 2.5, seed = 7L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data,
                          family = distsamp_open(cutpoints = cutp),
                          detection = ~ det_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  occu_multi = function(N = 150L) {
    sim <- simulate_occu_multi(S = 2L, N = N, J = 4L, seed = 5L)
    suppressWarnings(tobs(~ scov1, data = sim$data, family = occu_multi(),
                          detection = ~ 1, y = sim$y, species = sim$species,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  dyn_int_occu = function(N = 60L) {
    sim <- simulate_dyn_int_occu(N = N, T_seasons = 4L, S = 2L, J = 3L,
                                 psi1 = 0.5, gamma = 0.3, eps = 0.2,
                                 p = c(0.4, 0.6), seed = 25L)
    suppressWarnings(tobs(~ 1, data = sim$data, family = dyn_int_occu(),
                          detection = ~ 1, colonization = ~ 1,
                          extinction = ~ 1, y = sim$y, sources = sim$sources,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  t_occu = function(N = 100L) {
    sim <- simulate_t_occu(N = N, T_seasons = 6L, J = 3L, beta_occ = c(0.2, 0.6),
                           p = 0.4, rho = 0.6, sigma = 0.7, seed = 31L)
    suppressWarnings(tobs(~ x, data = sim$data, family = t_occu(),
                          detection = ~ 1, y = sim$y,
                          method = "pg_gibbs", control = .sbc_reg_ctl))
  },
  # ms_occu: S=20, not a smaller/faster count (R/sbc.R section 6j).
  # At S=5 this family's Laplace-EM posterior is measurably
  # non-Gaussian (validated against method="nuts", the exact reference:
  # Vf/Cinv 3-15x too narrow) and posterior SBC fails hard; at S=20 the
  # same plain Laplace-EM calibrates cleanly. Shrinking N/S here for
  # fixture speed would silently resurrect that failure.
  ms_occu = function(N = 80L, n_species = 20L) {
    sim <- simulate_ms_occu(N = N, J = 4L, n_species = n_species,
                            beta_comm_mean = c(0.2, 0.5), beta_comm_sd = c(0.6, 0.3),
                            alpha_comm_mean = c(0, 0.3), alpha_comm_sd = c(0.4, 0.2),
                            seed = 0L)
    suppressWarnings(tobs(~ x, data = sim$data, family = ms_occu(),
                          detection = ~ x, y = sim$y,
                          species = paste0("sp", seq_len(n_species)),
                          method = "laplace", control = .sbc_reg_ctl))
  },
  # ms_int_occu: S=14 (matching this family's own recovery-test fixture),
  # not a smaller/faster count (R/sbc.R section 6l). Shares ms_occu's
  # species-count-dependent non-Gaussianity; a smaller fixture silently
  # resurrects the original S~3 failure.
  ms_int_occu = function(N = 140L, n_species = 14L) {
    sim <- simulate_ms_int_occu(N = N, J = c(3, 4), n_species = n_species,
                                n_data = 2, seed = 0L)
    suppressWarnings(tobs(~ 1, data = sim$data, family = ms_int_occu(),
                          detection = ~ 1, y = sim$y,
                          species = paste0("sp", seq_len(n_species)),
                          method = "laplace", control = .sbc_reg_ctl))
  },
  ms_occu_cover = function(N = 60L) {
    sim <- simulate_ms_occu_cover(n_species = 4L, N = N, J = 5L,
                                  positive = "lognormal",
                                  sd_occ = 0.6, sd_p = 0.4, sd_pos = 0.3,
                                  seed = 3L)
    suppressWarnings(tobs(~ occ_cov1, data = sim$data,
                          family = ms_occu_cover("lognormal"),
                          detection = ~ 1, positive = ~ 1,
                          y = sim$y, y_pos = sim$y_pos,
                          species = paste0("sp", seq_len(4L)),
                          method = "laplace", control = .sbc_reg_ctl))
  },
  cover = function(N = 200L) {
    sim <- simulate_cover(N = N, beta_occ = c(-0.5, 0.8), beta_pos = c(-1.0, 0.3),
                          sigma_pos = 0.4, response = "lognormal", seed = 51L)
    suppressWarnings(tobs(~ x, data = sim$data, family = cover("lognormal"),
                          y = sim$y, method = "laplace", control = .sbc_reg_ctl))
  },
  # ms_int_occu: not registered (R/sbc.R section 6l). Confirmed (not just
  # suspected) to share ms_occu's failure mode.
  occu_multiscale_cover = function(n_cells = 40L) {
    sim <- simulate_occu_multiscale_cover(
      n_cells = n_cells, plots_per_cell = 4L, visits_per_plot = 2L,
      beta_psi = c(0.4, 0.6), beta_theta = c(0.2, 0.5),
      beta_p = c(0.0, 0.5), beta_pos = c(log(0.10), -0.4),
      positive = "lognormal", phi = 0.35, seed = 61L)
    suppressWarnings(tobs(~ x_cell + icar(graph = sim$adj, group_var = "cell"),
                          data = sim$data,
                          family = occu_multiscale_cover(response = "lognormal"),
                          detection = ~ x_pdet, availability = ~ x_plot,
                          positive = ~ x_cov, y = sim$y, y_pos = sim$y_pos,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  # ms_count: S=20 (matching ms_occu's resolved scale), not a smaller/faster
  # count (R/sbc.R section 6n). Multi-seed (0, 1, 2) originally confirmed the same Cov(mu, b_s)-adjacent bias, worse than
  # ms_occu/ms_int_occu on a small fixture; resolved the same species-count
  # way. A smaller fixture silently resurrects that failure.
  ms_count = function(N = 150L, n_species = 20L) {
    sim <- simulate_ms_count(N = N, n_species = n_species, response = "poisson",
                             seed = 0L)
    suppressWarnings(tobs(~ x, data = sim$data, family = ms_count("poisson"),
                          y = sim$y, species = paste0("sp", seq_len(n_species)),
                          method = "laplace", control = .sbc_reg_ctl))
  },
  # jsdm() shares ms_count()'s exact community Laplace-EM -- same S=20 scale
  # for the same reason.
  jsdm = function(N = 150L, n_species = 20L) {
    sim <- simulate_jsdm(N = N, n_species = n_species, seed = 0L)
    suppressWarnings(tobs(~ x, data = sim$data, family = jsdm(),
                          y = sim$y, species = paste0("sp", seq_len(n_species)),
                          method = "laplace", control = .sbc_reg_ctl))
  },
  ms_distance = function(N = 150L, n_species = 20L) {
    cutp <- c(0, 25, 50, 75, 100)
    sim <- simulate_ms_distance(n_species = n_species, N = N, cutpoints = cutp,
                                transect = "line", key = "halfnorm", seed = 0L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data,
                          family = ms_distance(key = "halfnorm", transect = "line",
                                               cutpoints = cutp),
                          detection = ~ 1, y = sim$y,
                          species = paste0("sp", seq_len(n_species)),
                          method = "laplace", control = .sbc_reg_ctl))
  },
  ms_dyn_occu = function(N = 150L, n_species = 20L, n_seasons = 4L, J = 3L) {
    sim <- simulate_ms_dyn_occu(N = N, J = J, n_species = n_species,
                                n_seasons = n_seasons, seed = 0L)
    suppressWarnings(tobs(~ 1, data = sim$data, family = ms_dyn_occu(),
                          detection = ~ 1, y = sim$y,
                          species = paste0("sp", seq_len(n_species)),
                          method = "laplace", control = .sbc_reg_ctl))
  },
  # ms_abun() needs the AGHQ/joint_fd engine (n.quad > 1) for Cinv/Bf to be
  # available at all -- the default optimizer = "em" (n.quad = 1) uses a
  # different engine that does not expose them.
  ms_abun = function(N = 150L, n_species = 20L, J = 3L) {
    sim <- simulate_ms_abun(n_species = n_species, N = N, J = J, seed = 0L)
    suppressWarnings(tobs(~ 1, data = sim$data,
                          family = ms_abun(mixture = "poisson"),
                          detection = ~ 1, y = sim$y,
                          species = paste0("sp", seq_len(n_species)),
                          method = "laplace",
                          control = utils::modifyList(
                            .sbc_reg_ctl, list(optimizer = "joint_fd", n.quad = 3L))))
  }
)


# occu_categorical carries no fit$means -- two independent Laplace-Gaussian
# blocks (presence, class), not tulpa's usual joint MVN summary -- so this
# flattens it with the SAME arm-prefix scheme .tobs_sbc_draws_occu_categorical()
# uses, letting the generic checks below read one name per reported
# coefficient regardless of family. A COMMUNITY family's fit$means is the
# community MEAN alone (P-length), never what its own SBC design ranks --
# a community family's theta is the per-species REALIZED coefficient (mu + b_s,
# S x P, species-major) plus any shared blocks, so it must be checked BEFORE
# the fit$means short-circuit below or it silently returns the wrong
# (right-sized-for-the-wrong-quantity) vector. The registry entry exposes the
# layout, so this reads it rather than restating each family's arm list.
.sbc_reg_means <- function(fit) {
  fam   <- attr(fit, "tobs_family")$name
  entry <- tulpaObs:::.TOBS_SBC_REGISTRY[[fam]]
  if (!is.null(entry$names)) {
    nm <- entry$names(fit$model)
    cm <- fit$ms_community
    theta <- as.vector(t(do.call(
      cbind, lapply(nm$arms, function(a) cm[[a$coef]]))))
    vals <- c(theta, fit$means[nm$global_cols])
    names(vals) <- nm$cols
    return(vals)
  }
  if (!is.null(fit$means)) return(fit$means)
  if (inherits(fit, "occu_categorical_fit")) {
    occ <- fit$beta_occ
    names(occ) <- paste0("occ_", names(occ))
    cls <- as.numeric(fit$beta_class)
    names(cls) <- as.vector(outer(rownames(fit$beta_class),
                                  colnames(fit$beta_class),
                                  function(r, cl) paste0("class_", cl, "_", r)))
    return(c(occ, cls))
  }
  if (inherits(fit, "cover_fit")) {
    occ <- fit$beta_occ; names(occ) <- paste0("occ_", names(occ))
    pos <- fit$beta_pos; names(pos) <- paste0("pos_", names(pos))
    return(c(occ, pos))
  }
  stop("no means accessor for this fit class (", paste(class(fit), collapse = "/"),
       ")", call. = FALSE)
}
