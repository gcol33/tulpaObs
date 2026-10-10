# ranef() returns one layout on every fit: `arm, group, level, term, estimate,
# std.error`, whichever branch of ranef.tobs_fit() a fit reaches. Each branch is
# exercised on a hand-built fit object, so the contract holds without fitting.

.ranef_cols <- c("arm", "group", "level", "term", "estimate", "std.error")

.ranef_fake_fit <- function(model_type, ...) {
  structure(c(list(model = list(
    model_type = model_type,
    process_info = list(list(name = "psi", coef_names = c("(Intercept)", "x")),
                        list(name = "p",   coef_names = "(Intercept)")))),
    list(...)), class = c("tobs_fit", "tulpa_fit"))
}

test_that("ranef(): the re_effects branch (single-species Laplace / NUTS)", {
  d <- list(group_label = "g1", n_groups = 2L, levels = c("a", "b"),
            coef_names = c("(Intercept)", "x"), group_expr = quote(g),
            process = 1L)
  tab <- .tobs_re_effects_table(d, estimate = 1:4 / 10, std.error = rep(0.2, 4))
  dv <- list(group_label = "p1", n_groups = 2L, levels = c("o1", "o2"),
             coef_names = "(Intercept)", group_expr = quote(obs),
             process = integer(0))
  tabv <- .tobs_re_effects_table(dv, estimate = c(-1, 1), std.error = c(0.3, 0.4))
  fit <- .ranef_fake_fit("occu", re_effects = list(g1 = tab, p1 = tabv))
  re <- ranef(fit)
  expect_identical(names(re), .ranef_cols)
  expect_identical(re$arm, c(rep("psi", 4L), "p", "p"))
  expect_identical(re$group, c(rep("g", 4L), "obs", "obs"))
  expect_identical(re$level, c("a", "b", "a", "b", "o1", "o2"))
  expect_identical(re$term, c("(Intercept)", "(Intercept)", "x", "x",
                              "(Intercept)", "(Intercept)"))
  expect_equal(re$estimate, c(1:4 / 10, -1, 1))
  expect_equal(re$std.error, c(rep(0.2, 4), 0.3, 0.4))
})

test_that("ranef(): an effect shared across arms names both", {
  d <- list(group_label = "g1", n_groups = 1L, levels = "a",
            coef_names = "(Intercept)", group_expr = NULL, process = 1:2)
  tab <- .tobs_re_effects_table(d, estimate = 0.5, std.error = 0.1)
  re <- ranef(.ranef_fake_fit("occu", re_effects = list(g1 = tab)))
  expect_identical(re$arm, "psi+p")
  expect_identical(re$group, "g1")
})

test_that("ranef(): the per-arm BLUP blocks of a joint occu_cover fit", {
  fit <- .ranef_fake_fit(
    "occu_cover",
    re = list(
      "p:habitat" = list(arm = "p", var = "habitat", levels = c("h1", "h2"),
                         blup = c(0.1, -0.1), blup_sd = c(0.5, 0.6)),
      pos = list(arm = "pos", var = "site", levels = c("s1", "s2", "s3"),
                 coef_names = c("(Intercept)", "z"),
                 blup = matrix(1:6, 3L, 2L, dimnames = list(NULL, c("(Intercept)", "z"))),
                 blup_sd = matrix(0.1, 3L, 2L))))
  re <- ranef(fit)
  expect_identical(names(re), .ranef_cols)
  expect_identical(re$arm, c("p", "p", rep("pos", 6L)))
  expect_identical(re$group, c("habitat", "habitat", rep("site", 6L)))
  expect_identical(re$level, c("h1", "h2", rep(c("s1", "s2", "s3"), 2L)))
  expect_identical(re$term, c("(Intercept)", "(Intercept)",
                              rep("(Intercept)", 3L), rep("z", 3L)))
  expect_equal(re$estimate, c(0.1, -0.1, 1:6))
  expect_equal(re$std.error, c(0.5, 0.6, rep(0.1, 6L)))
})

test_that("ranef(): the single BLUP block of a sampled count fit", {
  fit <- .ranef_fake_fit(
    "nmix",
    re = list(arm = "lambda", group_label = "g1", var = "g", n_groups = 3L,
              levels = c("a", "b", "c"), sigma = 0.5,
              blup = c(a = 0.2, b = 0, c = -0.2), blup_sd = c(0.3, 0.3, 0.3)))
  re <- ranef(fit)
  expect_identical(names(re), .ranef_cols)
  expect_identical(re$arm, rep("lambda", 3L))
  expect_identical(re$group, rep("g", 3L))
  expect_identical(re$level, c("a", "b", "c"))
  expect_identical(re$term, rep("(Intercept)", 3L))
  expect_equal(unname(re$estimate), c(0.2, 0, -0.2))
  # A block that records only the positional label reports that label.
  fit$re$var <- NULL
  expect_identical(ranef(fit)$group, rep("g1", 3L))
})

test_that("ranef(): a community fit's per-species deviations", {
  sp <- paste0("sp", 1:3)
  B_psi <- matrix(seq_len(6) / 10, 3L, 2L,
                  dimnames = list(sp, c("(Intercept)", "x")))
  B_p   <- matrix(c(-1, 0, 1), 3L, 1L, dimnames = list(sp, "(Intercept)"))
  fit <- .ranef_fake_fit("ms_occu", ms_community = list(blup_psi = B_psi,
                                                        blup_p = B_p))
  re <- ranef(fit)
  expect_identical(names(re), .ranef_cols)
  expect_identical(re$arm, c(rep("psi", 6L), rep("p", 3L)))
  expect_true(all(re$group == "species"))
  expect_identical(re$level, c(rep(sp, 2L), sp))
  expect_identical(re$term, c(rep("(Intercept)", 3L), rep("x", 3L),
                              rep("(Intercept)", 3L)))
  expect_equal(re$estimate, c(seq_len(6) / 10, -1, 0, 1))
  expect_true(all(is.na(re$std.error)))
  # A block the family does not carry (the NB log_r deviation under Poisson)
  # is skipped, not a row of NAs.
  fit2 <- .ranef_fake_fit("ms_nmix", ms_community = list(blup_lambda = B_psi,
                                                         blup_p = B_p))
  expect_identical(nrow(ranef(fit2)), 9L)
  expect_setequal(unique(ranef(fit2)$arm), c("lambda", "p"))
})

.ranef_spd <- function(d, seed) {
  set.seed(seed)
  A <- matrix(stats::rnorm(d * d), d, d)
  crossprod(A) + diag(0.1, d)
}

test_that("ranef(): community std.error is sqrt(diag(Cinv)) through the index map", {
  sp <- paste0("sp", 1:3)
  B_psi <- matrix(seq_len(6) / 10, 3L, 2L,
                  dimnames = list(sp, c("(Intercept)", "x")))
  B_p   <- matrix(c(-1, 0, 1), 3L, 1L, dimnames = list(sp, "(Intercept)"))
  Cinv  <- lapply(1:3, function(s) .ranef_spd(3L, s))
  # A non-contiguous, out-of-order map: psi sits at positions 3 and 1 of the
  # stacked vector, p at position 2.
  idx <- list(blup_psi = c(3L, 1L), blup_p = 2L)
  cm <- list(blup_psi = B_psi, blup_p = B_p, Cinv = Cinv, blup_idx = idx)
  re <- ranef(.ranef_fake_fit("ms_occu", ms_community = cm))
  sd_of <- function(j) vapply(Cinv, function(C) sqrt(C[j, j]), numeric(1))
  expect_equal(re$std.error, c(sd_of(3L), sd_of(1L), sd_of(2L)))
  expect_false(anyNA(re$std.error))

  # A species without a covariance (a failed per-species solve) reads NA on
  # its own rows only.
  cm$Cinv[2L] <- list(NULL)
  re2 <- ranef(.ranef_fake_fit("ms_occu", ms_community = cm))
  expect_identical(is.na(re2$std.error), rep(c(FALSE, TRUE, FALSE), 3L))

  # A map that addresses past the covariance, or a covariance with no map for
  # a block the fit reports, is an error rather than a column of NAs.
  bad <- cm; bad$Cinv <- Cinv; bad$blup_idx$blup_p <- 4L
  expect_error(ranef(.ranef_fake_fit("ms_occu", ms_community = bad)),
               "addresses position 4")
  bad <- cm; bad$Cinv <- Cinv; bad$blup_idx$blup_p <- NULL
  expect_error(ranef(.ranef_fake_fit("ms_occu", ms_community = bad)),
               "no index map")
  bad <- cm; bad$Cinv <- Cinv; bad$blup_idx$blup_psi <- 1L
  expect_error(ranef(.ranef_fake_fit("ms_occu", ms_community = bad)),
               "does not match the BLUP layout")
})

test_that("ranef(): every community handler carries the covariance through", {
  sp <- paste0("sp", 1:2)
  blk <- function(k) matrix(seq_len(2 * k) / 10, 2L, k,
                            dimnames = list(sp, paste0("c", seq_len(k))))
  # model_type -> the blup fields its ranef() handler reports.
  fields <- list(
    ms_occu       = c("blup_psi", "blup_p"),
    ms_dyn_occu   = c("blup_psi1", "blup_p"),
    ms_int_occu   = c("blup_psi", "blup_p1", "blup_p2"),
    ms_count      = "blup_mu",
    ms_occu_cover = c("blup_occ", "blup_p", "blup_pos"),
    ms_distance   = c("blup_lambda", "blup_sigma"),
    ms_nmix       = c("blup_lambda", "blup_p", "blup_logr", "blup_omega"))
  for (mt in names(fields)) {
    f <- fields[[mt]]
    width <- seq_along(f)
    off <- cumsum(c(0L, width))
    cm <- stats::setNames(lapply(width, blk), f)
    cm$blup_idx <- stats::setNames(lapply(seq_along(f), function(k)
      off[k] + seq_len(width[k])), f)
    d <- sum(width)
    cm$Cinv <- lapply(1:2, function(s) .ranef_spd(d, 10L * s))
    fit <- .ranef_fake_fit(mt, ms_community = cm)
    fit$model$process_names <- c("p1", "p2")
    re <- ranef(fit)
    expect_identical(nrow(re), 2L * d, info = mt)
    expect_false(anyNA(re$std.error), info = mt)
    ref <- unlist(lapply(f, function(fld) vapply(
      cm$blup_idx[[fld]], function(j)
        vapply(cm$Cinv, function(C) sqrt(C[j, j]), numeric(1)),
      numeric(2))), use.names = FALSE)
    expect_equal(re$std.error, ref, info = mt)
  }
})

# Fitted routes: each builder records the index map it slices b_s with, so the
# reported SD is checked against the fit's own covariance and, on a sampler,
# against the draws themselves.
.ranef_sd_from_cinv <- function(fit, field) {
  cm <- fit$ms_community
  do.call(rbind, lapply(cm$Cinv, function(C)
    sqrt(diag(as.matrix(C)))[cm$blup_idx[[field]]]))
}

.ranef_expect_cinv_sd <- function(fit, arms) {
  re <- ranef(fit)
  expect_false(anyNA(re$std.error))
  expect_length(fit$ms_community$Cinv, nrow(fit$ms_community[[arms[[1L]]]]))
  for (a in names(arms)) {
    expect_equal(re$std.error[re$arm == a],
                 as.numeric(.ranef_sd_from_cinv(fit, arms[[a]])), info = a)
  }
}

test_that("ranef(): fitted community Laplace routes report sqrt(diag(Cinv))", {
  skip_on_cran()
  sim <- simulate_ms_occu(N = 60, J = 4, n.species = 6, seed = 7)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:6), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  .ranef_expect_cinv_sd(fit, c(psi = "blup_psi", p = "blup_p"))

  s4 <- simulate_ms_count(N = 60, n.species = 5, response = "poisson", seed = 2)
  fit <- tobs(~ x, data = s4$data, family = ms_count(), y = s4$y,
              species = colnames(s4$y), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  .ranef_expect_cinv_sd(fit, c(mu = "blup_mu"))
})

test_that("ranef(): ms_abun EM and AGHQ routes carry the per-species covariance", {
  skip_on_cran()
  skip_if_fast()
  # Default n.quad = 1 EM: the covariance comes from cpp_nmix_community_em().
  s6 <- simulate_ms_abun(n.species = 5, N = 30, J = 3, n.abund.covs = 1,
                         n.det.covs = 1, seed = 3)
  fit <- tobs(~ abund_cov1, data = s6$data, y = s6$y, family = ms_abun(),
              detection = ~ det_cov1, species = s6$species, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$ms_community$optimizer, "em")
  expect_equal(dim(fit$ms_community$Cinv[[1L]]), c(4L, 4L))
  .ranef_expect_cinv_sd(fit, c(lambda = "blup_lambda", p = "blup_p"))

  # Negative binomial (joint AGHQ): the log_r deviation is the trailing
  # coordinate of the same covariance.
  s7 <- simulate_ms_abun(n.species = 4, N = 20, J = 3, n.abund.covs = 1,
                         n.det.covs = 1, mixture = "negbin", size = 5,
                         sigma.logr = 0.4, seed = 1)
  fit <- suppressWarnings(tobs(
    ~ abund_cov1, data = s7$data, y = s7$y,
    family = ms_abun(mixture = "negbin"), detection = ~ det_cov1,
    species = s7$species, method = "laplace",
    control = list(verbose = FALSE, progress = FALSE, n.quad = 2L,
                   max.iter = 20L)))
  expect_equal(dim(fit$ms_community$Cinv[[1L]]), c(5L, 5L))
  .ranef_expect_cinv_sd(fit, c(lambda = "blup_lambda", p = "blup_p",
                               logr = "blup_logr"))
})

test_that("ranef(): sampler routes report the SD of the sampled deviations", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_occu(N = 60, J = 4, n.species = 6, seed = 7)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:6), method = "nuts",
              control = list(verbose = FALSE, progress = FALSE, n.iter = 100L,
                             n.warmup = 100L, n.chains = 1L, seed = 1))
  .ranef_expect_cinv_sd(fit, c(psi = "blup_psi", p = "blup_p"))
  lay <- fit$nuts$layout; dr <- fit$nuts$draws
  Bd <- vapply(seq_len(nrow(dr)), function(i)
    tulpaObs:::.ms_ocs_b_from_z(dr[i, ], lay), matrix(0, 6L, lay$P))
  sd_draws <- apply(Bd, c(1L, 2L), stats::sd)
  idx <- fit$ms_community$blup_idx
  expect_equal(ranef(fit)$std.error,
               c(as.numeric(sd_draws[, idx$blup_psi]),
                 as.numeric(sd_draws[, idx$blup_p])))

  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:6), method = "pg_gibbs",
              control = list(verbose = FALSE, n.iter = 300L, n.warmup = 100L,
                             n.chains = 2L, seed = 1))
  .ranef_expect_cinv_sd(fit, c(psi = "blup_psi", p = "blup_p"))
})

test_that("ranef(): grid BLUP covariance is the law of total covariance", {
  # Two species, d = 2, three outer-grid nodes (one carrying zero weight).
  w <- c(0.5, 0.5, 0)
  blup_nodes <- list(matrix(c(1, 2, 3, 4), 2L), matrix(c(2, 0, 1, 4), 2L),
                     matrix(99, 2L, 2L))
  cov_nodes <- lapply(1:3, function(k)
    lapply(1:2, function(s) .ranef_spd(2L, 100L * k + s)))
  center <- 0.5 * blup_nodes[[1L]] + 0.5 * blup_nodes[[2L]]
  got <- tulpaObs:::.tobs_grid_blup_cov(w, blup_nodes, cov_nodes, center)
  for (s in 1:2) {
    ref <- matrix(0, 2L, 2L)
    for (k in 1:2) {
      dk <- blup_nodes[[k]][s, ] - center[s, ]
      ref <- ref + w[k] * (cov_nodes[[k]][[s]] + tcrossprod(dk))
    }
    expect_equal(got[[s]], (ref + t(ref)) / 2)
  }
})

test_that("ranef(): spatial community routes carry the per-species covariance", {
  skip_on_cran()
  skip_if_fast()
  adj <- rook_adj(4L)
  so <- simulate_ms_occu(N = 16, J = 4, n.species = 5, seed = 3)
  fit <- tobs(~ x + icar(graph = adj), data = so$data, family = ms_occu(),
              detection = ~ 1, y = so$y, species = paste0("sp", 1:5),
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  .ranef_expect_cinv_sd(fit, c(psi = "blup_psi", p = "blup_p"))

  ss <- simulate_ms_occu_cover_spatial(adj, n.species = 4L, J = 3L, K = 1L,
                                       sd.occ = 0.5, sd.load = 1.1,
                                       sigma.pos = 0.4, seed = 3L)
  fit <- tobs(~ occ_cov1 + icar(graph = adj), data = ss$data,
              family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
              positive = ~ pos_cov1, y = ss$y, y.pos = ss$y_pos,
              species = ss$species, method = "laplace",
              control = list(verbose = FALSE, n.factors = 1L, sd.load = 1.1))
  expect_identical(nrow(ranef(fit)), 4L * 6L)
  .ranef_expect_cinv_sd(fit, c(psi = "blup_occ", p = "blup_p",
                               pos = "blup_pos"))
})

test_that("ranef(): a fit without random effects is the zero-row table", {
  fit <- structure(list(model = list()), class = c("tobs_fit", "tulpa_fit"))
  re <- ranef(fit)
  expect_identical(names(re), .ranef_cols)
  expect_identical(nrow(re), 0L)
})
