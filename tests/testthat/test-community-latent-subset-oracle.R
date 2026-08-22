# The site-subset contract every community-latent oracle's `ll_cell(eta, idx)`
# implements: with `idx` given, the result must equal `ll_cell(eta)[idx, ,
# drop = FALSE]` -- evaluating only the requested sites is an optimisation,
# never a different answer. These are direct, ungated unit tests of that
# contract (not a full fit), so they stay fast and run on every push.
#
# The case worth a dedicated regression test is `length(idx) == 1`: `vapply`
# only returns a matrix when its `FUN.VALUE` has length > 1, so a single
# pending site silently collapsed the result to a plain vector and broke
# `rowSums()` downstream in `hval()` (R/community_latent.R) -- caught by the
# ms_occu / ms_distance bit-identical fit checks in dev_notes, not by a smoke
# test, which is why it is pinned here at the oracle level.

.subset_contract <- function(oracle, eta, Ns) {
  full <- oracle$ll_cell(eta)
  expect_equal(dim(full), c(Ns, oracle$n_species))
  for (idx in list(1L, Ns, c(1L, Ns), seq_len(Ns), sample(Ns, min(3L, Ns)))) {
    got <- oracle$ll_cell(eta, idx = idx)
    expect_equal(dim(got), c(length(idx), oracle$n_species))
    expect_equal(got, full[idx, , drop = FALSE])
  }
}

test_that(".tobs_ms_count_oracle() ll_cell(idx=) matches ll_cell()[idx, ]", {
  set.seed(1)
  Ns <- 7L; S <- 4L
  eta <- matrix(rnorm(Ns * S), Ns, S)
  y_pois <- matrix(rpois(Ns * S, 3), Ns, S)
  .subset_contract(tulpaObs:::.tobs_ms_count_oracle(y_pois, link = "log"), eta, Ns)
  y_bern <- matrix(rbinom(Ns * S, 1, 0.5), Ns, S)
  .subset_contract(tulpaObs:::.tobs_ms_count_oracle(y_bern, link = "logit"), eta, Ns)
})

test_that(".tobs_ms_distance_oracle() ll_cell(idx=) matches ll_cell()[idx, ]", {
  set.seed(2)
  Ns <- 7L; S <- 3L
  cutpoints <- c(0, 25, 50, 75, 100)
  y <- array(rpois(Ns * 4L * S, 2), c(Ns, 4L, S))
  model <- list(n_species = S, n_sites = Ns, n_bins = 4L, y = y,
               cutpoints = cutpoints, transect = "line", key = "halfnorm",
               quad_order = 32L)
  eng <- tulpaObs:::.tobs_ms_distance_engine(model)
  eta_sig_list <- lapply(seq_len(S), function(s) rep(log(20), Ns))
  or <- tulpaObs:::.tobs_ms_distance_oracle(eng, eta_sig_list, 0, Ns, S)
  eta <- matrix(rnorm(Ns * S, log(5), 0.3), Ns, S)
  .subset_contract(or, eta, Ns)
})

test_that(".tobs_ms_occu_oracle() ll_cell(idx=) matches ll_cell()[idx, ]", {
  set.seed(3)
  Ns <- 7L; S <- 3L
  eta_p <- matrix(rnorm(Ns * S, 0, 0.5), Ns, S)
  su <- lapply(seq_len(S), function(s) {
    y <- matrix(rbinom(Ns * 4L, 1, 0.3), Ns, 4L)
    valid <- matrix(TRUE, Ns, 4L)
    tulpaObs:::.ms_int_occu_sp_summary(list(y), list(valid))
  })
  or <- tulpaObs:::.tobs_ms_occu_oracle(su, eta_p)
  eta <- matrix(rnorm(Ns * S, 0, 0.5), Ns, S)
  .subset_contract(or, eta, Ns)
})

# ms_abun's oracle defers the actual site-subset optimisation (deferred:
# visit-indexed detection design), but must still honour the contract by
# slicing its full result -- the same test applies with no speedup expected.
test_that(".tobs_ms_abun_oracle() ll_cell(idx=) matches ll_cell()[idx, ]", {
  set.seed(4)
  Ns <- 7L; S <- 3L
  d <- simulate_ms_abun(n_species = S, N = Ns, J = 3L, seed = 4L)
  model <- .tobs_build_ms_abun(~ 1, ~ 1, d$data, d$y, colnames(d$y))
  ms <- tulpaObs:::.tobs_ms_abun_marginals(model)
  eta_p_list <- lapply(seq_len(S), function(s)
    as.numeric(ms$X_p[[s]] %*% rep(0, ncol(ms$X_p[[s]]))))
  or <- tulpaObs:::.tobs_ms_abun_oracle(ms$marg, eta_p_list, Ns, S)
  eta <- matrix(rnorm(Ns * S, log(3), 0.3), Ns, S)
  .subset_contract(or, eta, Ns)
})
