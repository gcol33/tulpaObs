# =============================================================================
# test-occu-cover-residual.R -- share(spatial(), residual = ) on occu_cover().
#
# The cover arm carries the occurrence field as `alpha * w` plus its own
# DEVIATION: `residual = "full"` gives that deviation one latent per node (the
# arm-specific block, in different coordinates), `residual = r` gives it r basis
# functions built orthogonal to the shared field.
#
# What is asserted here: the syntax resolves, the gates fire, the basis is
# orthogonal and correctly capped, `residual = "full"` is the SAME fit as the
# placed-bar spelling, and the rank sweep behaves as measured
# (NOTES_measurements.md). What is NOT asserted is the plan's original
# expectation that the sweep converges on a BETTER fit as r grows -- it converges
# on the full-rank one, and in a confounded fixture that one is the degenerate
# end, which is why the rank floor is a guard here rather than a target.
# =============================================================================

.res_fit <- function(sim, adj, positive, ...) {
  suppressWarnings(tobs(
    occurrence = ~ occ_cov1 + icar(graph = adj, group_var = "cell"),
    detection  = ~ 1,
    positive   = positive,
    family     = occu_cover(response = "lognormal"),
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    method = "nested_laplace",
    control = list(progress = FALSE, integration = "ccd"), ...))
}

# The fixture the rank sweep is measured on: a shared field the cover arm copies
# (alpha = 1) AND an independent cover-arm field, which is the configuration the
# amplitude and the deviation are confounded in. With alpha = 0 there is nothing
# to be confounded with and the question does not arise.
.res_sim <- function(adj, seed) {
  simulate_occu_cover(
    N = nrow(adj), J = 6L, positive = "lognormal",
    beta_occ = c(qlogis(0.7), 0.3), beta_p = c(qlogis(0.65), 0.1),
    beta_pos = c(log(0.25), 0.0), sigma_pos = 0.3, adj = adj,
    sigma = 0.5, alpha = 1.0,
    pos_field = TRUE, sigma_pos_int = 0.6, sigma_pos_trend = 0.0, seed = seed)
}

# The truth's OWN orthogonal decomposition, per seed. The simulator draws the
# shared field and the cover-arm field independently, but a finite draw of two
# independent fields is not orthogonal, so the amplitude the identified
# parameterization targets is the projection of the seed's realized cover
# surface onto its realized occurrence surface -- not the nominal `alpha = 1`.
# Scoring against the constant would spend most of the budget on draw noise
# (the estimand rule in CLAUDE.md, #155).
.res_truth <- function(sim) {
  tr    <- sim$truth
  s_occ <- tr$sigma * tr$f
  s_cov <- tr$alpha * s_occ + tr$sigma_pos_int * tr$g0
  a     <- sum(s_cov * s_occ) / sum(s_occ * s_occ)
  list(alpha = a, delta = s_cov - a * s_occ, s_occ = s_occ)
}


# --- syntax ------------------------------------------------------------------

test_that("residual = takes a whole number or \"full\", and nothing else", {
  cp <- tulpaObs:::.tobs_term_copy(spatial(), residual = "full")
  expect_identical(cp$residual$rank, "full")
  expect_identical(tulpaObs:::.tobs_term_copy(spatial(), residual = 50)$residual$rank,
                   50L)
  expect_null(tulpaObs:::.tobs_term_copy(spatial())$residual)

  expect_error(tulpaObs:::.tobs_term_copy(spatial(), residual = "half"),
               "the only word is \"full\"", fixed = TRUE)
  expect_error(tulpaObs:::.tobs_term_copy(spatial(), residual = 2.5),
               "whole number of basis functions", fixed = TRUE)
  expect_error(tulpaObs:::.tobs_term_copy(spatial(), residual = 0),
               "whole number of basis functions", fixed = TRUE)
})

test_that("print() reports the deviation the share() asked for", {
  expect_output(print(tulpaObs:::.tobs_term_copy(spatial(), residual = "full")),
                "full (one per node)", fixed = TRUE)
  expect_output(print(tulpaObs:::.tobs_term_copy(spatial(), residual = 7)),
                "rank 7", fixed = TRUE)
})


# --- the basis ---------------------------------------------------------------

test_that("the residual basis is orthogonal to the reference field it is built off", {
  adj <- rook_adj(6L)
  set.seed(1)
  w <- as.numeric(scale(stats::rnorm(nrow(adj))))
  for (r in c(1L, 2L, 5L, 15L)) {
    B <- tulpaObs:::.tobs_residual_basis(adj, w, r)
    expect_identical(dim(B), c(nrow(adj), r))
    # Exact, and by construction rather than on average: the projection is
    # applied to the basis, so EVERY coefficient vector gives a deviation
    # orthogonal to the shared field.
    expect_lt(max(abs(crossprod(B, w))), 1e-10)
    # And to the constant, so the deviation carries no level (the intercept
    # owns that).
    expect_lt(max(abs(colSums(B))), 1e-10)
    # Sorbye-Rue: unit geometric-mean marginal variance, so the coefficient SD
    # means the same thing here as a field SD does on an areal block.
    expect_equal(exp(mean(log(rowSums(B * B)))), 1, tolerance = 1e-8)
  }
})

test_that("the basis rank is capped by the graph, counting components", {
  adj <- rook_adj(4L)                       # 16 nodes, connected -> 15
  set.seed(2); w <- stats::rnorm(nrow(adj))
  expect_identical(ncol(tulpaObs:::.tobs_residual_basis(adj, w, 15L)), 15L)
  expect_error(tulpaObs:::.tobs_residual_basis(adj, w, 16L),
               "at most 15 basis functions")

  # Two disjoint components -> two null directions, so the cap drops by two.
  b <- rook_adj(3L); k <- nrow(b)
  adj2 <- matrix(0L, 2L * k, 2L * k)
  adj2[seq_len(k), seq_len(k)] <- b
  adj2[k + seq_len(k), k + seq_len(k)] <- b
  set.seed(3); w2 <- stats::rnorm(2L * k)
  expect_identical(ncol(tulpaObs:::.tobs_residual_basis(adj2, w2, 16L)), 16L)
  expect_error(tulpaObs:::.tobs_residual_basis(adj2, w2, 17L),
               "2 connected component")
})


# --- gates -------------------------------------------------------------------

test_that("a residual is refused where its block cannot go", {
  adj <- rook_adj(4L)
  sim <- .res_sim(adj, 1L)

  # Not a grid-integrated engine.
  expect_error(
    suppressWarnings(tobs(
      occurrence = ~ occ_cov1, detection = ~ 1,
      positive = ~ 1 + share(spatial(), residual = "full"),
      family = occu_cover(response = "lognormal"),
      data = sim$data, y = sim$y, y_pos = sim$y_pos, method = "laplace",
      control = list(progress = FALSE))),
    "rides the joint nested-Laplace engine", fixed = TRUE)

  # Two spellings of one thing: a placed cover-arm field AND a residual.
  expect_error(
    .res_fit(sim, adj, ~ 1 + spatial(~ 1 || cell, graph = adj) +
                        share(spatial(), residual = "full")),
    "already carries a field of its own", fixed = TRUE)

  # Two share()s, two answers.
  expect_error(
    .res_fit(sim, adj, ~ 1 + share(spatial(), residual = "full") +
                        share(spatial(), residual = 4)),
    "one share() names it", fixed = TRUE)

  # Nothing to deviate from.
  expect_error(
    .res_fit(sim, adj, ~ 1 + share(spatial(), alpha = 0, residual = "full")),
    "an amplitude pinned at 0 shares nothing", fixed = TRUE)

  # No shared field at all. The share() itself has nothing to select, so the
  # existing selector error is what a reader gets -- the residual gate does not
  # need to restate it, and would only bury the more specific message.
  expect_error(
    suppressWarnings(tobs(
      occurrence = ~ occ_cov1, detection = ~ 1,
      positive = ~ 1 + share(spatial(), residual = 4),
      family = occu_cover(response = "lognormal"),
      data = sim$data, y = sim$y, y_pos = sim$y_pos,
      method = "nested_laplace", control = list(progress = FALSE))),
    "needs a spatial field on the occurrence formula", fixed = TRUE)
})


# --- residual = "full" is the placed-bar fit, re-reported --------------------

test_that("residual = \"full\" fits the same model as the placed cover-arm bar", {
  skip_if_fast()
  skip_on_cran()

  adj <- rook_adj(5L)
  sim <- .res_sim(adj, 101L)
  a <- .res_fit(sim, adj, ~ 1 + share(spatial(), residual = "full"))
  b <- .res_fit(sim, adj, ~ 1 + spatial(~ 1 || cell, graph = adj) +
                            share(spatial()))

  # The same blocks, so the same fit -- to the bit, not merely close. The two
  # spellings compile to one prior list; only the reporting differs.
  expect_identical(unname(a$means), unname(b$means))
  expect_identical(a$spatial_field, b$spatial_field)
  expect_identical(a$pos_field,     b$pos_field)
  expect_equal(a$log_lik, b$log_lik)

  # What `residual =` adds is the identified reading, and only it carries one.
  expect_null(b$residual)
  expect_identical(a$residual$rank, "full")
  expect_identical(a$residual$n_basis, length(a$spatial_field))
  # The raw amplitude and the identified one are different numbers: the raw one
  # is what the outer grid integrated with the deviation free to absorb the
  # shared direction, the identified one is read off the surface they sum to.
  expect_false(isTRUE(all.equal(a$residual$alpha, a$residual$alpha_raw)))
})


# --- the deviation is orthogonal to the shared component, at every rank ------

test_that("the reported deviation is orthogonal to the occurrence surface", {
  skip_if_fast()
  skip_on_cran()

  adj <- rook_adj(5L)
  sim <- .res_sim(adj, 101L)
  for (spec in list(quote(share(spatial(), residual = "full")),
                    quote(share(spatial(), residual = 3)),
                    quote(share(spatial(), residual = 12)))) {
    f <- .res_fit(sim, adj, stats::reformulate(c("1", deparse(spec))))
    d <- f$residual
    expect_length(d$field, nrow(adj))
    expect_true(all(is.finite(d$field)))
    # Exact by construction: alpha is DEFINED as the projection, so what is left
    # has no component along it.
    expect_lt(abs(sum(d$field * d$occ_field)) /
                max(sum(d$occ_field^2), 1e-12), 1e-8)
    # The surface the two halves sum to is the one the fit actually carries.
    expect_equal(d$cover_field, d$alpha * d$occ_field + d$field,
                 tolerance = 1e-10)
  }
})

test_that("a rank above the graph's capacity is refused at the fit", {
  skip_if_fast()
  skip_on_cran()
  adj <- rook_adj(4L)
  sim <- .res_sim(adj, 5L)
  expect_error(.res_fit(sim, adj, ~ 1 + share(spatial(), residual = 40)),
               "at most 15 basis functions")
})


# --- composition: the deviation's blocks trail everything else ---------------
#
# The basis-coefficient blocks are appended LAST, after the field blocks and
# after the RE blocks, and both of those are read by position: the field copy
# indices name field blocks, and each RE descriptor is a position relative to
# the field count. Appending last is what keeps both readings valid, so a fit
# carrying an RE or a second (trend) field alongside a deviation is the case
# that would catch it if that ever stopped being true.

test_that("a deviation composes with a psi RE and with a trend field", {
  skip_if_fast()
  skip_on_cran()

  adj <- rook_adj(5L)
  sim <- .res_sim(adj, 3L)
  sim$data$grp <- factor(rep(1:5, length.out = nrow(adj)))

  ortho_of <- function(f) abs(sum(f$residual$field * f$residual$occ_field)) /
                          sum(f$residual$occ_field^2)

  for (rr in list("\"full\"", "4")) {
    f <- suppressWarnings(tobs(
      occurrence = ~ occ_cov1 + icar(graph = adj, group_var = "cell") + re(grp),
      detection  = ~ 1,
      positive   = stats::reformulate(
        c("1", sprintf("share(spatial(), residual = %s)", rr))),
      family = occu_cover(response = "lognormal"),
      data = sim$data, y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
      control = list(progress = FALSE, integration = "ccd")))
    # The RE is still read off its own block, at its own position.
    expect_true("sigma_re" %in% names(f$means))
    expect_true(is.finite(f$means[["sigma_re"]]) && f$means[["sigma_re"]] > 0)
    expect_length(f$residual$field, nrow(adj))
    expect_lt(ortho_of(f), 1e-8)
  }

  # And with a coupled trend field, where the deviation's blocks trail two
  # field blocks rather than one.
  f <- suppressWarnings(tobs(
    occurrence = ~ occ_cov1 + spatial(~ 1 + time || cell, graph = adj),
    detection  = ~ 1,
    positive   = ~ 1 + share(spatial(), residual = 4),
    family = occu_cover(response = "lognormal"),
    data = sim$data, y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
    control = list(progress = FALSE, integration = "ccd")))
  expect_false(is.null(f$trend_field))
  expect_length(f$residual$field, nrow(adj))
  expect_lt(ortho_of(f), 1e-8)
})


# --- step 6: the rank sweep, against the full-rank formulation ---------------
#
# Measured over the FULL seed set, never a subset: a three-seed read of the
# neighbouring #299 fixture inverted its own verdict, and this fixture family is
# exposed to exactly that.
#
# The measurement (12 seeds, NOTES_measurements.md) separates two effects that
# the plan had bundled together:
#
#   * The AMPLITUDE's identification comes from the deviation's SD being pinned,
#     not from the rank. Bias against the seed's own orthogonal amplitude is
#     +1.04 for the full-rank block with a free SD, +0.36 for the same block
#     held at the pinned SD, and +0.28 / +0.42 / +0.49 at ranks 2 / 8 / 32.
#   * The DEVIATION's recovery comes from the rank. Even pinned, the full-rank
#     block returns a deviation of SD 0.028 against a truth of 0.62 -- it puts
#     essentially all of its mass along the shared direction, which is the #110
#     confounding -- while ranks 2 / 8 / 32 return 0.37 / 0.29 / 0.25.
#
# So the sweep DOES converge on the full-rank fit as r grows, and converging on
# it is a loss, not a gain. The assertions below are gross-regression guards on
# the measured behaviour, not reproductions of it.

test_that("the rank-r deviation recovers structure the full-rank block does not", {
  skip_if_fast()
  skip_on_cran()

  adj   <- rook_adj(8L)
  seeds <- 1:12
  ranks <- c(2L, 8L, 32L)

  rec <- lapply(seeds, function(s) {
    sim <- .res_sim(adj, s)
    tt  <- .res_truth(sim)
    one <- function(f) {
      d <- f$residual
      c(alpha_err = d$alpha - tt$alpha,
        sd_dev    = stats::sd(d$field),
        cor_dev   = stats::cor(d$field, tt$delta),
        ortho     = abs(sum(d$field * d$occ_field)) / sum(d$occ_field^2))
    }
    out <- list(full = one(.res_fit(sim, adj,
                                    ~ 1 + share(spatial(), residual = "full"))))
    for (r in ranks) {
      out[[paste0("r", r)]] <- one(.res_fit(
        sim, adj, stats::reformulate(
          c("1", sprintf("share(spatial(), residual = %d)", r)))))
    }
    out
  })
  grab <- function(tag, what) vapply(rec, function(x) unname(x[[tag]][what]),
                                     numeric(1))

  # Orthogonality is exact at every rank and every seed -- the one thing here
  # that is a property of the construction rather than of the data.
  for (tag in c("full", paste0("r", ranks))) {
    expect_lt(max(grab(tag, "ortho")), 1e-8)
  }

  sd_full <- stats::median(grab("full", "sd_dev"))
  sd_r8   <- stats::median(grab("r8",   "sd_dev"))
  cor_r8  <- stats::median(grab("r8",   "cor_dev"))

  # The full-rank block collapses the orthogonal deviation (measured 0.075
  # against a truth of 0.62); a rank-8 deviation keeps a real one (0.287). The
  # ratio is ~3.8x, so a 2x floor is a regression guard, not the measurement.
  expect_gt(sd_r8, 2 * sd_full)
  # ... and it is structure, not noise: median correlation with the seed's own
  # orthogonal truth is 0.44 (worst seed 0.17).
  expect_gt(cor_r8, 0.20)

  # The identified amplitude is recovered to within a fraction of its own scale.
  # Median absolute error is 0.42 at rank 8; the guard is wide enough that only
  # a collapse or a runaway reaches it.
  expect_lt(stats::median(abs(grab("r8", "alpha_err"))), 0.9)
  # Per-seed gross-regression guard, budgeted as one.
  expect_lt(max(abs(grab("r8", "alpha_err"))), 3)

  # Rank is a REGULARIZER here, not only a cost knob: the deviation shrinks
  # toward the full-rank block's collapsed one as r grows. Asserted as an
  # ordering of the extremes over the seed set, which is what the 12-seed read
  # supports (0.367 -> 0.251); the intermediate ordering is reported, not
  # asserted.
  expect_gt(stats::median(grab("r2", "sd_dev")),
            stats::median(grab("r32", "sd_dev")))
})
