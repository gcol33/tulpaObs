# =============================================================================
# test-sbc-registry.R
# -- the SBC family registry
#
# `test-sbc.R` measures the coupled `occu_cover` route. This file covers the
# OTHER half of the registry: the families whose site marginals multiply, and
# the registry contract itself.
#
# Two tiers:
#
#   STRUCTURE (always on, no fitting). Every entry supplies the callbacks the
#   driver reads, and every callback slot holds a function. A family cannot be
#   added half-registered and discovered at the first run.
#
#   CONTRACT (skip_on_cran). Per family: the callbacks compose, the pooled data
#   set carries both blocks under disjoint site labels, and refitting the
#   observed data ALONE reproduces the observed fit -- which is what says the
#   rebuilt formulas, family and method are the ones it came from. A fit
#   carrying a structured term is refused rather than scored on its
#   coefficients alone.
#
# The ACCEPTANCE tier -- the calibration measurement, the reported posterior
# uniform and a deliberately mis-scaled control not -- is one file per family,
# test-sbc-acceptance-<family>.R. Fixtures and the means accessor are in
# helper-sbc-registry.R.
# =============================================================================


# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------

test_that("every registry entry supplies the callbacks the driver reads", {
  reg <- tulpaObs:::.TOBS_SBC_REGISTRY
  req <- tulpaObs:::.TOBS_SBC_REQUIRED
  expect_gt(length(reg), 0L)

  for (fam in names(reg)) {
    e <- reg[[fam]]
    info <- paste("family", fam)
    expect_true(all(req %in% names(e)), info = info)
    for (slot in req) expect_true(is.function(e[[slot]]), info = info)
    # Optional slots, each a function when present. A rank arm needs one of the
    # two log-likelihood readers; the rest fall back to the shared driver's.
    for (slot in c("spec", "data", "pool", "loglik", "loglik_many")) {
      if (!is.null(e[[slot]])) expect_true(is.function(e[[slot]]), info = info)
    }
    expect_true(!is.null(e$loglik) || !is.null(e$loglik_many), info = info)
    # The key IS the family name `tobs()` dispatches on, so an entry filed under
    # a name no fit reports would be unreachable.
    expect_true(exists(fam, envir = asNamespace("tulpaObs"),
                       mode = "function"), info = info)
  }

  expect_identical(tulpaObs:::.tobs_sbc_registered(), sort(names(reg)))
  # The fixtures below cover every registered family bar the coupled one, which
  # test-sbc.R owns. A new entry lands here without a fixture and this fails.
  expect_setequal(names(reg),
                  c("occu_cover", names(.SBC_REG_FIXTURES)))
})


test_that("the roster in the error message names what is registered", {
  expect_error(sbc(list()), "No sbc method is registered")
  fake <- structure(list(model = list()), class = "tobs_fit",
                    tobs_family = list(name = "not_a_family"))
  expect_error(sbc(fake), "not registered for family")
  expect_error(sbc(fake), "occu_cover")
})


# ---------------------------------------------------------------------------
# Contract
# ---------------------------------------------------------------------------

test_that("each registered family composes its callbacks end to end", {
  skip_on_cran()

  for (fam in names(.SBC_REG_FIXTURES)) {
    fit <- .SBC_REG_FIXTURES[[fam]]()
    expect_identical(attr(fit, "tobs_family")$name, fam)

    m <- sbc(fit, model.only = TRUE)
    expect_true(all(c("data_obs", "fit", "draw_theta", "simulate", "pool",
                      "arms", "group_ids") %in% names(m)), info = fam)
    # Every fixed effect the fit reports is scored; nothing is silently fixed.
    expect_setequal(attr(m, "quantities"), names(.sbc_reg_means(fit)))
    expect_length(attr(m, "fixed"), 0L)

    # Refitting the observed data ALONE reproduces the observed fit: the
    # rebuilt formula, family, method and control are the ones it came from.
    again <- m$fit(m$data_obs)
    expect_equal(.sbc_reg_means(again), .sbc_reg_means(fit), tolerance = 1e-8,
                info = fam)

    # `draw_theta` is reproducible from its seed alone and moves with it, which
    # is what the driver's RNG split relies on.
    t1 <- m$draw_theta(fit, 11L)
    expect_identical(t1, m$draw_theta(fit, 11L), info = fam)
    expect_false(isTRUE(all.equal(unname(t1),
                                  unname(m$draw_theta(fit, 12L)))), info = fam)

    # Pooling keeps BOTH data sets under disjoint site labels: fitting the
    # replicate alone would be ordinary SBC under a hand-made prior, not the
    # posterior experiment, and shared labels would deny the premise
    # `tulpa::sbc()` verifies.
    rep <- m$simulate(t1, 4400L)
    pl  <- m$pool(m$data_obs, rep)
    # Row count (one row per unit in `cells`/`y`'s first axis) and GROUP count
    # (unique exchangeable-unit labels) coincide for every family except
    # occu_multiscale_cover, where several plot-rows share one cell group id
    # -- the premise check is about groups, the slicing below is about rows,
    # and conflating them under one n_o/n_r silently passed for every family
    # that happens to have one row per group.
    n_row_o <- nrow(m$data_obs$cells); n_row_r <- nrow(rep$cells)
    n_grp_o <- length(unique(m$group_ids(m$data_obs)))
    n_grp_r <- length(unique(m$group_ids(rep)))
    expect_length(unique(m$group_ids(pl)), n_grp_o + n_grp_r)
    expect_identical(nrow(pl$cells), n_row_o + n_row_r)
    # `y` is a plain 2D matrix for most families, a 3D [site x visit x season]
    # array for the multi-season group (pooled on the site axis alone), a 4D
    # [site x visit x season x species] array for the community dynamic group,
    # a list of one matrix (or, for the product shape, one 3D array) per source
    # for the multi-source group (each pooled on the site axis independently),
    # or a plain length-N vector for a family with one observation per unit and
    # no visit/season axis (occu_categorical). Every one of those pools on the
    # first axis and leaves the rest alone, so the slice indexes axis 1 at
    # whatever rank the family brings rather than enumerating the ranks.
    slice_site <- function(z) {
      d <- dim(z)
      if (is.null(d)) return(z[seq_len(n_row_o)])
      do.call(`[`, c(list(z), list(seq_len(n_row_o)),
                     rep(list(bquote()), length(d) - 1L), list(drop = FALSE)))
    }
    if (is.list(pl$y) && !is.data.frame(pl$y)) {
      expect_identical(vapply(pl$y, nrow, integer(1)),
                       stats::setNames(rep(n_row_o + n_row_r, length(pl$y)), names(pl$y)),
                       info = fam)
      obs_slice <- lapply(pl$y, slice_site)
    } else if (is.null(dim(pl$y))) {
      expect_identical(length(pl$y), n_row_o + n_row_r)
      obs_slice <- pl$y[seq_len(n_row_o)]
    } else {
      expect_identical(nrow(pl$y), n_row_o + n_row_r)
      obs_slice <- slice_site(pl$y)
    }
    # Values only -- a family's own y may carry dimnames (site/bin/season
    # labels) the freshly-built pooled array never does; that carries no
    # information this check is about.
    strip_dn <- function(z) if (is.list(z)) lapply(z, unname) else unname(z)
    expect_equal(strip_dn(obs_slice), strip_dn(m$data_obs$y), info = fam)

    pooled <- m$fit(pl)
    expect_s3_class(pooled, "tobs_fit")
    expect_true(all(is.finite(.sbc_reg_means(pooled))), info = fam)

    # The rank arm reads the family's own marginal, so a reference set and the
    # truth score through one kernel and come back finite.
    a <- m$arms(pooled, list(theta = t1))
    expect_true("log_lik" %in% names(a$posterior), info = fam)
    expect_setequal(setdiff(names(a$posterior), "log_lik"),
                    attr(m, "quantities"))
  }
})


test_that("a structured term is refused, not scored on the coefficients", {
  skip_on_cran()
  fit <- .SBC_REG_FIXTURES$occu(N = 40L)

  # The field is a latent quantity shared across sites that theta does not
  # hold; scoring the coefficients alone would report a calibration the
  # experiment did not measure.
  spatial <- fit
  spatial$spatial <- list(type = "icar")
  expect_error(sbc(spatial, model.only = TRUE), "structured term")
  expect_error(sbc(spatial, model.only = TRUE), "fresh cells")

  re <- fit
  re$re_effects <- list(g = 1)
  expect_error(sbc(re, model.only = TRUE), "structured term")
})


test_that("a visit-level observation design is refused, not rebuilt", {
  skip_on_cran()

  # The replicate comes from the family's own simulate() kernel, which builds
  # detection from the SITE-level design; a visit-level column would be scored
  # by the refit and missing from the data it sees, and the ranks would not say
  # so. Refused at build time instead.
  N <- 40L; J <- 3L
  sim <- simulate_occu(N = N, J = J, seed = 31L)
  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)),
                     det_v = stats::rnorm(N * J))
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = "det_v")
  fit <- suppressWarnings(tobs(~ occ_cov1, data = sim$data, family = occu(),
                               detection = ~ det_v, y = od$y,
                               visits = od$det.covs, method = "laplace",
                               control = .sbc_reg_ctl))
  expect_true(ncol(fit$model$X_det_visit) > 0L)
  expect_error(sbc(fit, model.only = TRUE), "visit-level")
})


test_that("the replicate generator draws at the theta it is handed", {
  skip_on_cran()

  # count(): the mean response has a closed form in theta alone, so the
  # generator can be checked against it rather than against itself. The band is
  # absolute and sized on the estimator's own Monte-Carlo error: at 200
  # replicates of 120 Poisson sites the standard error of the pooled mean is
  # under 0.02, so 0.08 is four of them.
  fit <- .SBC_REG_FIXTURES$count()
  m <- sbc(fit, model.only = TRUE)
  th <- m$draw_theta(fit, 3L)
  mu <- exp(as.vector(fit$model$X_occ %*% th[seq_len(ncol(fit$model$X_occ))]))

  reps <- vapply(seq_len(200L),
                 function(i) mean(m$simulate(th, 77000L + i)$y), numeric(1))
  expect_lt(abs(mean(reps) - mean(mu)), 0.08)
})

