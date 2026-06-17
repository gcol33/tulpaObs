# occu_cover() random intercept on the observation arms (detection / positive
# cover), composed with the shared spatial field on the nested-Laplace joint
# engine (gcol33/tulpaObs#102). A per-visit categorical grouping (e.g. an EUNIS
# habitat class) enters the detection (or cover) linear predictor as an iid RE
# block whose variance integrates on the outer grid and whose per-group BLUPs are
# partially pooled. The block rides ONE arm: the detection arm carries field_coef
# = 1 (so the iid block scatters) while the shared field is still skipped on
# detection by its 0-node sentinel.

.ocor_grid_adj <- function(side) {
  ng <- side * side
  co <- expand.grid(x = seq_len(side), y = seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i != j && abs(co$x[i] - co$x[j]) + abs(co$y[i] - co$y[j]) == 1L) adj[i, j] <- 1L
  adj
}

.ocor_sim <- function(seed, side = 7L, J = 6L, n_g = 6L, sigma_re_p = 0.8) {
  adj <- .ocor_grid_adj(side)
  simulate_occu_cover(
    N = nrow(adj), J = J, n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    positive = "lognormal", adj = adj, sigma = 0.6, alpha = 0.6,
    re_det_groups = n_g, sigma_re_p = sigma_re_p, seed = seed)
}

test_that("occu_cover() detection RE: fit runs and reports the RE on the p arm", {
  skip_on_cran()
  sim <- .ocor_sim(1L, side = 6L, n_g = 6L)
  adj <- .ocor_grid_adj(6L)
  fit <- tobs(occurrence = ~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = occu_cover("lognormal"),
              detection = ~ det_cov1 + (1 | habitat), positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nested_laplace")
  # Detection-arm RE hyperparameter + per-group BLUPs, keyed by arm "p".
  expect_true("sigma_re_p" %in% names(fit$means))
  expect_gt(fit$means[["sigma_re_p"]], 0)
  expect_false(is.null(fit$re$p))
  expect_identical(fit$re$p$arm, "p")
  expect_identical(fit$re$p$var, "habitat")
  expect_length(fit$re$p$blup, 6L)
  expect_identical(fit$re$p$levels, paste0("hab", 1:6))
  # ranef() surfaces the BLUP table with an `arm` column and the level labels.
  rf <- ranef(fit)
  expect_s3_class(rf, "data.frame")
  expect_true(all(c("arm", "group", "blup", "blup_sd") %in% names(rf)))
  expect_true(all(rf$arm == "p"))
  expect_equal(nrow(rf), 6L)
})

test_that("occu_cover() detection RE: a no-RE fit is unchanged (no fit$re)", {
  skip_on_cran()
  sim <- .ocor_sim(2L, side = 6L, n_g = 6L)
  adj <- .ocor_grid_adj(6L)
  fit0 <- tobs(occurrence = ~ occ_cov1 + icar(graph = adj), data = sim$data,
               family = occu_cover("lognormal"),
               detection = ~ det_cov1, positive = ~ pos_cov1,
               y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
               method = "nested_laplace",
               control = list(verbose = FALSE, progress = FALSE))
  expect_null(fit0$re)
  expect_false(any(grepl("sigma_re", names(fit0$means))))
})

test_that("occu_cover() detection RE: predict() handles seen and unseen levels", {
  skip_on_cran()
  sim <- .ocor_sim(3L, side = 6L, n_g = 6L)
  adj <- .ocor_grid_adj(6L)
  N   <- nrow(adj)
  fit <- tobs(occurrence = ~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = occu_cover("lognormal"),
              detection = ~ det_cov1 + (1 | habitat), positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  nd <- data.frame(occ_cov1 = sim$data$occ_cov1)
  nd$habitat <- factor(c(rep(paste0("hab", 1:6), length.out = N - 2L),
                         "newA", "newB"))
  pd <- predict(fit, newdata = nd, type = "detection", draws = FALSE)
  expect_s3_class(pd, "tobs_prediction")
  expect_true(all(c("cell", "mean", "lwr", "upr") %in% names(pd)))
  expect_true(all(pd$mean > 0 & pd$mean < 1))
  # An unseen habitat shrinks to the population-mean detection (offset 0).
  p_pop  <- stats::plogis(fit$means[["p_(Intercept)"]])
  unseen <- which(nd$habitat %in% c("newA", "newB"))
  expect_equal(pd$mean[unseen], rep(p_pop, length(unseen)), tolerance = 0.02)
  # A seen-level habitat differs from the population mean for at least one group.
  expect_gt(max(abs(pd$mean[-unseen] - p_pop)), 0.02)
})

test_that("occu_cover() positive-cover RE: fit runs and reports the RE on the pos arm", {
  skip_on_cran()
  sim <- .ocor_sim(4L, side = 6L, n_g = 6L)
  adj <- .ocor_grid_adj(6L)
  fit <- tobs(occurrence = ~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = occu_cover("lognormal", cover_aggregate = "none"),
              detection = ~ det_cov1, positive = ~ pos_cov1 + (1 | habitat),
              y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_false(is.null(fit$re$pos))
  expect_identical(fit$re$pos$arm, "pos")
  expect_true("sigma_re_pos" %in% names(fit$means))
})

test_that("occu_cover() detection RE recovers BLUPs and detects the variance", {
  skip_on_cran()
  skip_if_fast()
  seeds <- 1:6
  res <- t(vapply(seeds, function(s) {
    sim <- .ocor_sim(s, side = 9L, J = 8L, n_g = 6L, sigma_re_p = 0.8)
    adj <- .ocor_grid_adj(9L)
    fit <- tobs(occurrence = ~ occ_cov1 + icar(graph = adj), data = sim$data,
                family = occu_cover("lognormal"),
                detection = ~ det_cov1 + (1 | habitat),
                positive  = ~ pos_cov1 + copy(spatial()),
                y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
                method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE,
                               integration = "ccd"))
    c(blup_cor = stats::cor(fit$re$p$blup, sim$truth$b_p_re),
      sd_re    = unname(fit$means[["sigma_re_p"]]))
  }, numeric(2)))
  mn <- colMeans(res)
  # BLUPs recover the per-group structure (strong rank/scale agreement).
  expect_gt(mn[["blup_cor"]], 0.7)
  # The variance is detected (positive). Binary-RE inner-Laplace attenuation makes
  # the grid-integrated sigma a lower bound on the small-cluster truth, so it is
  # checked as detected rather than exact (matches the occupancy-arm RE path).
  expect_gt(mn[["sd_re"]], 0.3)
})

test_that("occu_cover() observation RE gates the unsupported configurations", {
  skip_on_cran()
  sim <- .ocor_sim(5L, side = 5L, n_g = 5L)
  adj <- .ocor_grid_adj(5L)
  base <- list(data = sim$data, family = occu_cover("lognormal"),
               y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
               control = list(progress = FALSE))
  # A random slope on the detection arm needs the joint-engine weighted-iid /
  # multivariate free-Sigma blocks tracked in gcol33/tulpa#114; gated until then.
  expect_error(
    do.call(tobs, c(list(occurrence = ~ occ_cov1 + icar(graph = adj),
                         detection = ~ (det_cov1 | habitat), positive = ~ pos_cov1,
                         method = "nested_laplace"), base)),
    "tulpa#114")
  # An observation-arm RE needs the joint nested-Laplace engine.
  expect_error(
    do.call(tobs, c(list(occurrence = ~ occ_cov1,
                         detection = ~ det_cov1 + (1 | habitat), positive = ~ pos_cov1,
                         method = "laplace"), base)),
    "nested_laplace")
  # A positive-cover RE needs per-visit cover; an explicit cell-aggregated cover
  # arm (cover_aggregate = "mean") has one row per unit and cannot carry it.
  expect_error(
    tobs(occurrence = ~ occ_cov1 + icar(graph = adj),
         detection = ~ det_cov1, positive = ~ occ_cov1 + (1 | habitat),
         data = sim$data,
         family = occu_cover("lognormal", cover_aggregate = "mean"),
         y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
         method = "nested_laplace", control = list(progress = FALSE)),
    "per-visit cover")
})

# --- crossed / nested random intercepts on the detection arm (#103) -----------

test_that("occu_cover() crossed detection RE: two groupings fit and report", {
  skip_on_cran()
  adj <- .ocor_grid_adj(6L)
  sim <- simulate_occu_cover(
    N = nrow(adj), J = 6L, n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    positive = "lognormal", adj = adj, sigma = 0.6, alpha = 0.6,
    re_det_groups = 6L, sigma_re_p = 0.8,
    re_det = list(observer = list(K = 4L, sigma = 0.6, prefix = "obs")),
    seed = 11L)
  fit <- tobs(occurrence = ~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = occu_cover("lognormal"),
              detection = ~ det_cov1 + (1 | habitat) + (1 | observer),
              positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE,
                             integration = "ccd",
                             re.sigma.grid.p = c(0.4, 1.0)))
  expect_s3_class(fit, "tobs_fit")
  # Two p-arm RE blocks, keyed "p:<var>", each with its own variance component.
  expect_false(is.null(fit$re[["p:habitat"]]))
  expect_false(is.null(fit$re[["p:observer"]]))
  expect_identical(fit$re[["p:habitat"]]$arm, "p")
  expect_identical(fit$re[["p:observer"]]$var, "observer")
  expect_length(fit$re[["p:habitat"]]$blup, 6L)
  expect_length(fit$re[["p:observer"]]$blup, 4L)
  expect_true(all(c("sigma_re_p_habitat", "sigma_re_p_observer") %in%
                    names(fit$means)))
  # ranef() stacks both groupings with an `arm` + `var` column.
  rf <- ranef(fit)
  expect_true(all(c("arm", "var", "group", "blup", "blup_sd") %in% names(rf)))
  expect_setequal(unique(rf$var), c("habitat", "observer"))
  expect_equal(nrow(rf), 10L)
})

test_that("occu_cover() nested detection RE: (1 | region/site) fits two terms", {
  skip_on_cran()
  adj <- .ocor_grid_adj(6L)
  sim <- simulate_occu_cover(
    N = nrow(adj), J = 6L, n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    positive = "lognormal", adj = adj, sigma = 0.6, alpha = 0.6,
    re_det = list(region = list(K = 3L, sigma = 0.7, prefix = "region"),
                  site   = list(K = 2L, sigma = 0.5, prefix = "s",
                                nested_in = "region")),
    seed = 12L)
  fit <- tobs(occurrence = ~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = occu_cover("lognormal"),
              detection = ~ det_cov1 + (1 | region/site), positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE,
                             integration = "ccd",
                             re.sigma.grid.p = c(0.4, 1.0)))
  expect_s3_class(fit, "tobs_fit")
  # Desugars to re(region) + re(region:site): two p-arm blocks (3 + 6 groups).
  expect_false(is.null(fit$re[["p:region"]]))
  expect_false(is.null(fit$re[["p:region:site"]]))
  expect_equal(fit$re[["p:region"]]$n_groups, 3L)
  expect_equal(fit$re[["p:region:site"]]$n_groups, 6L)
})

test_that("occu_cover() crossed RE: predict sums both groupings, unseen shrinks", {
  skip_on_cran()
  adj <- .ocor_grid_adj(6L)
  N   <- nrow(adj)
  sim <- simulate_occu_cover(
    N = N, J = 6L, n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    positive = "lognormal", adj = adj, sigma = 0.6, alpha = 0.6,
    re_det_groups = 6L, sigma_re_p = 0.8,
    re_det = list(observer = list(K = 4L, sigma = 0.6, prefix = "obs")),
    seed = 13L)
  fit <- tobs(occurrence = ~ occ_cov1 + icar(graph = adj), data = sim$data,
              family = occu_cover("lognormal"),
              detection = ~ det_cov1 + (1 | habitat) + (1 | observer),
              positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE,
                             integration = "ccd",
                             re.sigma.grid.p = c(0.4, 1.0)))
  nd <- data.frame(occ_cov1 = sim$data$occ_cov1)
  nd$habitat  <- factor(rep(paste0("hab", 1:6), length.out = N))
  nd$observer <- factor(c(rep(paste0("obs", 1:4), length.out = N - 1L), "newObs"))
  pd <- predict(fit, newdata = nd, type = "detection", draws = FALSE)
  expect_true(all(pd$mean > 0 & pd$mean < 1))
  # An unseen observer level shrinks that term to 0, so the row's detection is
  # driven only by the (seen) habitat term + fixed effect -- finite, in (0, 1).
  unseen <- which(nd$observer == "newObs")
  expect_true(all(is.finite(pd$mean[unseen])))
})

test_that(".occu_cover_obs_re_design builds slope Z (scaffolding for tulpa#114)", {
  # Engine-independent: a random-slope bar parses + designs into a 2-coef
  # correlated block with intercept + slope columns, even though the joint fit is
  # gated on gcol33/tulpa#114. Drives the design builder directly.
  pp <- tulpaObs:::.occu_cover_obs_re_parse(~ x + (area | habitat), "detection")
  expect_true(pp$has_slope)
  expect_length(pp$terms, 1L)
  n_sites <- 4L; max_visits <- 2L
  vdf <- data.frame(
    habitat = factor(rep(c("a", "b"), length.out = n_sites * max_visits)),
    area    = as.numeric(seq_len(n_sites * max_visits)))
  valid <- matrix(TRUE, n_sites, max_visits)
  des <- tulpaObs:::.occu_cover_obs_re_design(
    pp, data = data.frame(row = seq_len(n_sites)), visit_df = vdf,
    valid = valid, n_sites = n_sites, max_visits = max_visits, arm = "detection")
  d <- des[[1L]]
  expect_identical(d$type, "slope")
  expect_true(d$correlated)
  expect_equal(d$n_coefs, 2L)
  expect_identical(d$coef_names, c("(Intercept)", "area"))
  expect_equal(ncol(d$Z), 2L)
  expect_equal(nrow(d$Z), n_sites * max_visits)
  expect_true(all(d$Z[, 1L] == 1))            # intercept column
})

test_that("occu_cover() crossed detection RE recovers both variances + BLUPs", {
  skip_on_cran()
  skip_if_fast()
  seeds <- 1:5
  res <- t(vapply(seeds, function(s) {
    adj <- .ocor_grid_adj(9L)
    sim <- simulate_occu_cover(
      N = nrow(adj), J = 8L, n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
      positive = "lognormal", adj = adj, sigma = 0.6, alpha = 0.6,
      re_det_groups = 6L, sigma_re_p = 0.8,
      re_det = list(observer = list(K = 5L, sigma = 0.6, prefix = "obs")),
      seed = s)
    fit <- tobs(occurrence = ~ occ_cov1 + icar(graph = adj), data = sim$data,
                family = occu_cover("lognormal"),
                detection = ~ det_cov1 + (1 | habitat) + (1 | observer),
                positive  = ~ pos_cov1 + copy(spatial()),
                y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
                method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE,
                               integration = "ccd"))
    re_h <- fit$re[["p:habitat"]]; re_o <- fit$re[["p:observer"]]
    th_h <- sim$truth$re_det$habitat$b[re_h$levels]   # align by level label
    th_o <- sim$truth$re_det$observer$b[re_o$levels]
    c(cor_h = stats::cor(re_h$blup, th_h),
      cor_o = stats::cor(re_o$blup, th_o),
      sd_h  = unname(fit$means[["sigma_re_p_habitat"]]),
      sd_o  = unname(fit$means[["sigma_re_p_observer"]]))
  }, numeric(4)))
  mn <- colMeans(res)
  # Both crossed groupings recover their per-group structure and detect variance
  # (binary-RE inner-Laplace attenuation makes the integrated sigma a lower bound,
  # so it is checked as detected, matching the single-grouping path).
  expect_gt(mn[["cor_h"]], 0.6)
  expect_gt(mn[["cor_o"]], 0.6)
  expect_gt(mn[["sd_h"]], 0.25)
  expect_gt(mn[["sd_o"]], 0.2)
})
