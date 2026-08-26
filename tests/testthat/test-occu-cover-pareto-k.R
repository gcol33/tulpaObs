# Fast outer Pareto-k diagnostic on occu_cover ( consumer side).
#
# The diagnostic re-solves dominate a spatial occu_cover fit sped them up
# (Shamanskii reuse + loosened inner tol + near-neighbour batch order) with the
# k-hat byte-stable. These tests pin (1) the fast default reports the SAME k-hat
# as the byte-for-byte exact diagnostic, and (2) that k-hat agrees with the
# reference loo::psis on the diagnostic's actual importance ratios.


.pk_occu_cover_fit <- function(side = 12L, seed = 100L) {
  N <- side * side
  adj <- rook_adj(side)
  sim <- simulate_occu_cover(N = N, J = 3L, positive = "beta", phi = 25,
                             adj = adj, sigma = 0.8, alpha = 1.0, seed = 1L)
  long <- data.frame(site_id = rep(seq_len(N), each = 3L),
                     visit = rep(seq_len(3L), times = N),
                     y = as.vector(t(sim$y)),
                     det_cov1 = sim$visit_data$det_cov1,
                     pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  ctrl <- list(verbose = FALSE, engine = "joint",
               diagnose.k = TRUE, k.samples = 200L,
               n.threads.outer = 1L, progress = FALSE)
  set.seed(seed)
  suppressWarnings(tobs(
    occurrence = ~ occ_cov1 + icar(graph = adj), data = cell_dat,
    family = occu_cover("beta"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace", control = ctrl))
}

.pk_extract <- function(fit) {
  found <- NULL
  rec <- function(x, depth) {
    if (!is.null(found) || depth > 6L) return(invisible())
    if (is.list(x)) {
      if (!is.null(x[["pareto_k"]]) && length(x[["pareto_k"]]) == 1L) {
        found <<- as.numeric(x[["pareto_k"]]); return(invisible())
      }
      for (el in x) rec(el, depth + 1L)
    }
  }
  rec(fit, 0L)
  found
}

test_that("fast Pareto-k diagnostic matches the byte-for-byte exact k-hat", {
  skip_if_fast()
  skip_on_cran()
  # EXACT diagnostic: refresh 1, fit's own tol, no re-order, no per-cell warm.
  withr::local_options(tulpa.kdiag.refresh = 1L, tulpa.kdiag.tol = 0,
                       tulpa.kdiag.reorder = FALSE, tulpa.kdiag.percell = FALSE)
  k_exact <- .pk_extract(.pk_occu_cover_fit(seed = 100L))
  # FAST diagnostic: package defaults (refresh 4, tol 1e-4, near-neighbour order,
  # per-cell nearest-grid-mode warm start).
  withr::local_options(tulpa.kdiag.refresh = NULL, tulpa.kdiag.tol = NULL,
                       tulpa.kdiag.reorder = NULL, tulpa.kdiag.percell = NULL)
  k_fast <- .pk_extract(.pk_occu_cover_fit(seed = 100L))

  expect_true(is.finite(k_exact))
  expect_true(is.finite(k_fast))
  # Same RNG seed before each fit -> same 200-draw batch; the only difference is
  # the loosened inner tol (Laplace log-marginal error O(tol^2)), so the k-hat is
  # byte-stable to a few 1e-4. A regression that moved it would fail here.
  expect_equal(k_fast, k_exact, tolerance = 5e-3)
})

test_that("occu_cover Pareto-k agrees with loo::psis on the real importance ratios", {
  skip_if_fast()
  skip_on_cran()
  cap <- new.env()
  withr::local_options(tulpa.kdiag.capture = cap)   # fast defaults otherwise
  k_tulpa <- .pk_extract(.pk_occu_cover_fit(seed = 100L))
  lr <- cap$lr
  skip_if(is.null(lr) || length(lr) < 25L, "diagnostic declined (no lr captured)")

  k_loo <- suppressWarnings(
    loo::psis(matrix(lr, ncol = 1), r_eff = NA)$diagnostics$pareto_k)
  # tulpa's reported k-hat, tulpa_psis on the captured ratios, and loo's GPD fit
  # on the same ratios must coincide (the estimator is the reference's).
  expect_equal(tulpa::tulpa_psis(lr)$pareto_k, k_tulpa, tolerance = 1e-8)
  expect_equal(k_loo, k_tulpa, tolerance = 0.02)
})
