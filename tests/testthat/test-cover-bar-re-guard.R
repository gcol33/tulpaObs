# Soft guard for the cover()/occu_cover() bar-vs-spatial-field papercut. A bare `|
# / ||` lme4 bar in a cover formula is a legitimate random effect, NOT a spatial
# field. RE bars must keep fitting as REs (behavior unchanged). The guard only
# INFORMS (message, not warning/error) in the strong-signal confusion case: the
# bar's grouping factor is also an areal term's graph-node group_var. These run at
# parse level (.encode_cover_terms / .occu_cover_spatial_fields) -- no model fit,
# so they run always.

# Small chain adjacency on `n` nodes (self-contained).
.guard_chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (i in seq_len(n - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  adj
}

.guard_data <- function(n_cells = 16L, N = 64L, seed = 1) {
  set.seed(seed)
  data.frame(
    cell = rep(seq_len(n_cells), length.out = N),
    site = rep(seq_len(n_cells %/% 2L), length.out = N),
    time = rnorm(N)
  )
}

# ---- 1. same factor as graph node -> message -------------------------------

test_that("a bar sharing the areal group_var emits the guidance message", {
  adj <- .guard_chain_adj(16L)
  dat <- .guard_data()
  f <- ~ time + (1 + time | cell) + icar(graph = adj, group_var = "cell")

  expect_message(
    tulpaObs:::.encode_cover_terms(f, dat),
    "random effect|spatial"
  )
})

test_that("the same-factor message also fires through occu_cover's parser", {
  adj <- .guard_chain_adj(16L)
  dat <- .guard_data()
  f <- ~ time + (1 | cell) + icar(graph = adj, group_var = "cell")

  expect_message(
    tulpaObs:::.occu_cover_spatial_fields(f, dat),
    "random effect|spatial"
  )
})

# ---- 2. unrelated factor -> silent -----------------------------------------

test_that("a bar on an unrelated factor does not message", {
  adj <- .guard_chain_adj(16L)
  dat <- .guard_data()
  # bar groups by `site`; the field's graph-node group_var is `cell`.
  f <- ~ time + (1 + time | site) + icar(graph = adj, group_var = "cell")

  expect_no_message(tulpaObs:::.encode_cover_terms(f, dat))
})

test_that("a bar with no spatial term in the formula does not message", {
  dat <- .guard_data()
  f <- ~ time + (1 + time | cell)

  expect_no_message(tulpaObs:::.encode_cover_terms(f, dat))
})

# ---- 3. behavior unchanged: the bar still parses as a tobs_re ---------------

test_that("the same-factor bar still fits as a random effect (re populated)", {
  adj <- .guard_chain_adj(16L)
  dat <- .guard_data()
  f <- ~ time + (1 + time | cell) + icar(graph = adj, group_var = "cell")

  enc <- suppressMessages(tulpaObs:::.encode_cover_terms(f, dat))
  expect_false(is.null(enc$re))
  expect_length(enc$re, 1L)
  expect_s3_class(enc$re[[1L]], "tobs_re")
  # The field is still resolved from the areal term (the bar did not become one).
  expect_false(is.null(enc$spatial))
})

test_that("the no-spatial bar parses as a tobs_re with spatial NULL", {
  dat <- .guard_data()
  f <- ~ time + (1 + time | cell)

  enc <- tulpaObs:::.encode_cover_terms(f, dat)
  expect_false(is.null(enc$re))
  expect_length(enc$re, 1L)
  expect_s3_class(enc$re[[1L]], "tobs_re")
  expect_null(enc$spatial)
})
