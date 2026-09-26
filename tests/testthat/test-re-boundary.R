# Grouped random-effect components tested against their boundary on the
# single-species AGHQ fits (#375 SD at zero, #382 correlation at +-1): one
# record per log-scale SD coordinate on `convergence(fit)$re_boundary`, one
# warning per fit.

# A 2 x 2 full block in tulpa's log-Cholesky packing (log L11, L21, log L22),
# with the SE of each coordinate supplied.
full_ref <- function(L11, L21, L22, se = c(0.2, 0.3, 0.2)) {
  L <- matrix(c(L11, L21, 0, L22), 2L, 2L)
  list(re_par = c(log(L11), L21, log(L22)),
       re_par_se = stats::setNames(se, c("log_L11", "L21", "log_L22")),
       re_par_layout = list(list(label = "re1", nc = 2L, full = TRUE,
                                 index = 1:3,
                                 coord = c("log_L11", "L21", "log_L22"))),
       Sigma_list = list(L %*% t(L)))
}
slope_design <- list(list(group_label = "g1", coef_names = c("(Intercept)", "w")))

test_that("a full block yields an SD record and a correlation record", {
  recs <- .tobs_re_boundary(full_ref(0.8, 0.6, 0.5), slope_design)
  expect_named(recs, c("sigma_g1_(Intercept)", "cor_g1_(Intercept)_w"))
  expect_identical(recs[[1]]$kind, "sd")
  expect_identical(recs[[2]]$kind, "correlation")
  expect_equal(recs[[1]]$estimate, 0.8)
  # cor = L21 / sqrt(L21^2 + L22^2)
  expect_equal(recs[[2]]$estimate, 0.6 / sqrt(0.6^2 + 0.5^2))
  expect_true(all(vapply(recs, `[[`, logical(1), "distinguishable")))
})

test_that("a correlation at +-1 is the conditional SD at zero", {
  # L22 -> 0 with a flat log L22: the curvature says nothing separates the
  # block from a singular one.
  recs <- .tobs_re_boundary(full_ref(0.8, 0.9, 1e-3, se = c(0.2, 0.3, 25)),
                            slope_design)
  cr <- recs[["cor_g1_(Intercept)_w"]]
  expect_gt(abs(cr$estimate), 0.999)
  expect_false(cr$distinguishable)
  expect_warning(.tobs_warn_re_boundary(recs),
                 "cor_g1_\\(Intercept\\)_w.*re.lkj")
})

test_that("a diagonal block tests each log-SD", {
  ref <- list(re_par = log(c(0.5, 0.01)),
              re_par_se = c(log_sd_1 = 0.2, log_sd_2 = 30),
              re_par_layout = list(list(label = "re1", nc = 2L, full = FALSE,
                                        index = 1:2,
                                        coord = c("log_sd_1", "log_sd_2"))),
              Sigma_list = list(diag(c(0.25, 1e-4))))
  recs <- .tobs_re_boundary(ref, slope_design)
  expect_named(recs, c("sigma_g1_(Intercept)", "sigma_g1_w"))
  expect_true(recs[[1]]$distinguishable)
  expect_false(recs[[2]]$distinguishable)
})

test_that("the fit tail warns on a boundary record and not on a clear one", {
  bad <- .tobs_re_boundary(full_ref(0.8, 0.9, 1e-3, se = c(0.2, 0.3, 25)),
                           slope_design)
  fit <- list(method = "laplace",
              convergence = list(converged = TRUE, n_iter = 5L,
                                 re_boundary = bad))
  expect_warning(.tobs_finalize_convergence(fit, list(name = "occu")),
                 "at the boundary")
  fit$convergence$re_boundary <- .tobs_re_boundary(full_ref(0.8, 0.6, 0.5),
                                                   slope_design)
  expect_no_warning(.tobs_finalize_convergence(fit, list(name = "occu")))
})

sim_fp_psi_re <- function(n_groups, per_group, J, seed, beta0 = qlogis(0.45),
                          sigma_re = 0.7, p11 = 0.65, p10 = 0.05, b = 0.5) {
  set.seed(seed)
  N <- n_groups * per_group; g <- rep(seq_len(n_groups), each = per_group)
  bg <- rnorm(n_groups, 0, sigma_re)
  z <- rbinom(N, 1L, plogis(beta0 + bg[g])); y <- matrix(0L, N, J)
  for (i in seq_len(N)) if (z[i] == 1L) {
    det <- rbinom(J, 1L, p11); cert <- rbinom(J, 1L, b)
    y[i, ] <- ifelse(det == 1L, ifelse(cert == 1L, 2L, 1L), 0L)
  } else y[i, ] <- rbinom(J, 1L, p10)
  list(y = y, data = data.frame(g = factor(g)))
}

fit_fp_re <- function(seed) {
  s <- sim_fp_psi_re(16L, 10L, 5L, seed)
  w <- character()
  fit <- withCallingHandlers(
    tobs(~ 1 + (1 | g), data = s$data, family = fp_occu(), detection = ~ 1,
         y = s$y, method = "laplace",
         control = list(n.quad = 5L, progress = FALSE, verbose = FALSE)),
    warning = function(x) { w <<- c(w, conditionMessage(x))
                            invokeRestart("muffleWarning") })
  list(fit = fit, w = w)
}

test_that("fp_occu() psi (1|g) collapsing to zero is recorded and warned (#375)", {
  skip_on_cran()
  r <- fit_fp_re(21L)
  rec <- convergence(r$fit)$re_boundary[["sigma_g1_(Intercept)"]]
  expect_lt(rec$estimate, 0.05)
  expect_false(rec$distinguishable)
  expect_length(grep("at the boundary", r$w), 1L)
})

test_that("fp_occu() psi (1|g) away from zero stays silent", {
  skip_on_cran()
  r <- fit_fp_re(22L)
  rec <- convergence(r$fit)$re_boundary[["sigma_g1_(Intercept)"]]
  expect_true(rec$distinguishable)
  expect_length(grep("at the boundary", r$w), 0L)
})
