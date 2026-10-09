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

test_that("ranef(): a fit without random effects is the zero-row table", {
  fit <- structure(list(model = list()), class = c("tobs_fit", "tulpa_fit"))
  re <- ranef(fit)
  expect_identical(names(re), .ranef_cols)
  expect_identical(nrow(re), 0L)
})
