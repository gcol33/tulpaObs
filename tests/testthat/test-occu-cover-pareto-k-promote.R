# Promote the joint nested-Laplace outer Pareto-k diagnostic to the tobs_fit top
# level + glance(), including its mode-Hessian-vs-grid-moment proposal source. The
# raw joint engine attaches pareto_k / pareto_k_is_ess / pareto_k_scope /
# pareto_k_proposal_source to the object the postprocess wrappers nest at
# $joint_fit; a user reading fit$pareto_k directly should not have to reach into
# $joint_fit.

# --------------------------------------------------------------------------- #
# Extractor: surface a ran diagnostic, stay inert when diagnose.k was off       #
# --------------------------------------------------------------------------- #

# pareto_k_is_ess is the importance-sampling ESS on the PSIS-smoothed weights
# (numeric; the engine prints it as "IS-ESS = ..."), NOT a boolean flag.
.pkp_jf <- function(k = 0.42, is_ess = 180, src = "mode_hessian") {
  list(pareto_k = k, pareto_k_is_ess = is_ess,
       pareto_k_scope = "outer (hyperparameter) Gaussian proposal",
       pareto_k_proposal_source = src)
}

test_that(".tobs_promote_pareto_k surfaces a diagnostic that ran", {
  pk <- tulpaObs:::.tobs_promote_pareto_k(.pkp_jf(0.51, 182, "mode_hessian"))
  expect_equal(pk$pareto_k, 0.51)
  expect_equal(pk$pareto_k_is_ess, 182)          # numeric IS-ESS, carried as-is
  expect_identical(pk$pareto_k_proposal_source, "mode_hessian")
  expect_identical(pk$pareto_k_scope, "outer (hyperparameter) Gaussian proposal")
})

test_that(".tobs_promote_pareto_k surfaces on a finite IS-ESS alone (gate OR branch)", {
  # The `ran` gate accepts a usable number from either field: a finite IS-ESS
  # surfaces the diagnostic even if the k-hat itself came back NA.
  pk <- tulpaObs:::.tobs_promote_pareto_k(.pkp_jf(NA_real_, 41, NA_character_))
  expect_true(is.na(pk$pareto_k))
  expect_equal(pk$pareto_k_is_ess, 41)
})

test_that(".tobs_promote_pareto_k is inert when diagnose.k was off (all NA)", {
  # diagnose.k OFF: the engine leaves every field NA -> nothing to promote.
  off <- list(pareto_k = NA_real_, pareto_k_is_ess = NA_real_,
              pareto_k_scope = "outer (hyperparameter) Gaussian proposal",
              pareto_k_proposal_source = NA_character_)
  expect_null(tulpaObs:::.tobs_promote_pareto_k(off))
  # A non-joint fit carries none of the fields.
  expect_null(tulpaObs:::.tobs_promote_pareto_k(list(foo = 1)))
  expect_null(tulpaObs:::.tobs_promote_pareto_k(NULL))
})

# --------------------------------------------------------------------------- #
# glance.tobs_fit: extra columns only when the diagnostic ran                   #
# --------------------------------------------------------------------------- #

.pkp_fit <- function(jf = NULL) {
  structure(list(
    n_fixed = 4L, n_samples = 1000L, log_prob = -120.5, converged = TRUE,
    joint_fit = jf
  ), class = c("tobs_fit", "tulpa_fit"))
}

test_that("glance.tobs_fit adds pareto-k columns for a joint-coupled fit", {
  g <- glance(.pkp_fit(.pkp_jf(0.33, 176, "mode_hessian")))
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 1L)
  expect_true(all(c("pareto_k", "pareto_k_is_ess",
                    "pareto_k_proposal_source") %in% names(g)))
  expect_equal(g$pareto_k, 0.33)
  expect_equal(g$pareto_k_is_ess, 176)           # carried numeric, not coerced to logical
  expect_identical(g$pareto_k_proposal_source, "mode_hessian")
})

test_that("glance.tobs_fit reads the promoted top-level fields too", {
  # A fit whose postprocess already promoted the fields to the top level (the
  # production path) glances identically to one that only carries $joint_fit.
  fit <- .pkp_fit(NULL)
  fit$pareto_k <- 0.61; fit$pareto_k_is_ess <- 158
  fit$pareto_k_proposal_source <- "grid_moment"
  g <- glance(fit)
  expect_equal(g$pareto_k, 0.61)
  expect_equal(g$pareto_k_is_ess, 158)
  expect_identical(g$pareto_k_proposal_source, "grid_moment")
})

test_that("glance.tobs_fit adds no pareto-k columns when the diagnostic was off", {
  g_off  <- glance(.pkp_fit(NULL))                                   # no joint fit
  g_na   <- glance(.pkp_fit(list(pareto_k = NA_real_,
                                 pareto_k_is_ess = NA_real_,
                                 pareto_k_proposal_source = NA_character_)))
  for (g in list(g_off, g_na)) {
    expect_s3_class(g, "data.frame")
    expect_false(any(grepl("^pareto_k", names(g))))
  }
})

# --------------------------------------------------------------------------- #
# End-to-end: a real spatial occu_cover fit promotes the diagnostic            #
# --------------------------------------------------------------------------- #

test_that("occu_cover() spatial fit surfaces pareto_k at the top level + glance", {
  skip_if_fast()
  skip_on_cran()

  N <- 40L; J <- 5L
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  sim <- simulate_occu_cover(N = N, J = J, positive = "lognormal",
                             adj = adj, sigma = 1, alpha = 1, seed = 104L)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

  fit_args <- function(diag_k) {
    suppressWarnings(tobs(
      occurrence = ~ occ_cov1 + icar(graph = adj), data = cell_dat,
      family = occu_cover("lognormal"),
      detection = ~ det_cov1, positive = ~ pos_cov1,
      y = od$y, y_pos = y_pos, visits = od$det.covs,
      method = "nested_laplace",
      control = list(verbose = FALSE, diagnose.k = diag_k, k.samples = 200L)
    ))
  }

  fit <- fit_args(TRUE)
  # Promoted to the top level: a user reading fit$pareto_k directly works.
  expect_true("pareto_k" %in% names(fit))
  expect_true(is.numeric(fit$pareto_k) && is.finite(fit$pareto_k))
  # The IS-ESS surfaces as a finite number, not a boolean flag.
  expect_true(is.numeric(fit$pareto_k_is_ess) && is.finite(fit$pareto_k_is_ess))
  expect_true("pareto_k_proposal_source" %in% names(fit))
  # The full set of outer-proposal sources the joint engine reports: the
  # single-Gaussian grid-moment proposal, its moment-matching refinement, the
  # grid-mixture proposal, the FD mode-Hessian delta-collapse fallback, or NA
  # when diagnose.k stayed off.
  expect_true(fit$pareto_k_proposal_source %in%
                c("mode_hessian", "grid_moment", "moment_matched",
                  "grid_mixture", NA_character_))
  expect_identical(fit$pareto_k_scope, "outer (hyperparameter) Gaussian proposal")
  # Same values as the nested raw object it was promoted from.
  expect_equal(fit$pareto_k, fit$joint_fit$pareto_k)
  expect_identical(fit$pareto_k_proposal_source,
                   fit$joint_fit$pareto_k_proposal_source)

  # glance() surfaces them, carrying the IS-ESS as the same finite number as the
  # top-level field (not coerced to logical).
  g <- glance(fit)
  expect_true(all(c("pareto_k", "pareto_k_proposal_source") %in% names(g)))
  expect_equal(g$pareto_k, fit$pareto_k)
  expect_equal(g$pareto_k_is_ess, fit$pareto_k_is_ess)

  # diagnose.k OFF (the default): inert -- no top-level field, no glance column.
  fit0 <- fit_args(FALSE)
  expect_false("pareto_k" %in% names(fit0))
  expect_false(any(grepl("^pareto_k", names(glance(fit0)))))
})
