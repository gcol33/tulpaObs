# test-ms-ocs-theta-length.R - the ms_ocs C++ entry points hand theta's buffer to
# kernels that index it by a layout derived from the spec, so a vector of any
# other length is read past its end (gcol33/tulpaObs#236). Each entry point is
# pinned against the R mirror of that layout (.ms_ocs_npar_inner for the packed
# inner latent, .ms_ocs_nuts_layout for the full NUTS coordinate vector).

.ocs_len_grid_adj <- function(nr, nc) {
  N <- nr * nc
  adj <- matrix(0L, N, N)
  idx <- function(r, c) (c - 1L) * nr + r
  for (r in seq_len(nr)) for (c in seq_len(nc)) {
    k <- idx(r, c)
    if (r > 1L)  adj[k, idx(r - 1L, c)] <- 1L
    if (r < nr)  adj[k, idx(r + 1L, c)] <- 1L
    if (c > 1L)  adj[k, idx(r, c - 1L)] <- 1L
    if (c < nc)  adj[k, idx(r, c + 1L)] <- 1L
  }
  adj
}

# A small icar fixture plus the two vectors the entry points take, packed from
# truth: `inner` is the packed inner latent, `theta` the full NUTS coordinate
# vector at the requested loading parameterisation.
.ocs_len_fixture <- function(K = 1L, constrain = FALSE, seed = 321L) {
  adj <- .ocs_len_grid_adj(4L, 4L)              # N = 16 cells
  S <- 4L
  sim <- simulate_ms_occu_cover_spatial(adj, n_species = S, K = K, J = 3L,
                                        seed = seed)
  model <- tulpaObs:::.tobs_build_ms_occu_cover_spatial(
    occ_formula = ~ occ_cov1, det_formula = ~ det_cov1, pos_formula = ~ pos_cov1,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    positive = "lognormal", species = sim$species, adj = adj, K = K)
  d  <- tulpaObs:::.ms_ocs_dims(model)
  tr <- sim$truth
  Lmat <- if (K == 1L) matrix(tr$L, S, 1L) else tr$L
  Lblock <- if (constrain) tulpaObs:::.ms_ocs_L_to_lfree(Lmat, S, K)
            else as.numeric(Lmat)
  mu <- c(tr$mu_occ, tr$mu_p, tr$mu_pos)
  b  <- as.numeric(t(cbind(tr$b_occ, tr$b_p, tr$b_pos)))
  inner <- c(mu, b, Lblock, as.numeric(tr$w), log(tr$sigma_pos))
  chol_v <- function(Sig) tulpaObs:::.ms_ocs_chol_pack(t(chol(Sig)))
  theta <- c(inner, chol_v(diag(0.4^2, d$P_occ)),
             chol_v(diag(0.35^2, d$P_p)), chol_v(diag(0.3^2, d$P_pos)),
             log(rep(1.3, K)))
  list(d = d, spec = tulpaObs:::.ms_ocs_nuts_spec(model),
       pri = tulpaObs:::.ms_ocs_nuts_priors(), inner = inner, theta = theta)
}


test_that("the marginal entry points reject a theta_inner off the packed layout (#236)", {
  fx <- .ocs_len_fixture(K = 1L)
  n_inner <- as.integer(tulpaObs:::.ms_ocs_npar_inner(fx$d, FALSE))
  expect_identical(length(fx$inner), n_inner)

  # The accepted length is exactly the R mirror's, and the gradient comes back
  # over the same layout.
  expect_true(is.finite(tulpaObs:::cpp_ms_ocs_marginal_ll(fx$spec, fx$inner)))
  expect_length(tulpaObs:::cpp_ms_ocs_marginal_grad(fx$spec, fx$inner), n_inner)

  for (bad in list(fx$inner[-n_inner], c(fx$inner, 0))) {
    expect_error(tulpaObs:::cpp_ms_ocs_marginal_ll(fx$spec, bad),
                 "theta_inner length")
    expect_error(tulpaObs:::cpp_ms_ocs_marginal_grad(fx$spec, bad),
                 "theta_inner length")
  }
})

test_that("cpp_ms_ocs_joint_logpost rejects a theta off the NUTS layout (#236)", {
  fx  <- .ocs_len_fixture(K = 1L)
  lay <- tulpaObs:::.ms_ocs_nuts_layout(fx$d, FALSE)
  expect_identical(length(fx$theta), as.integer(lay$total))

  ok <- tulpaObs:::cpp_ms_ocs_joint_logpost(fx$spec, fx$theta, fx$pri, 5, 1.0,
                                            FALSE)
  expect_true(is.finite(ok$lp))
  expect_length(ok$grad, as.integer(lay$total))

  for (bad in list(fx$theta[-lay$total], c(fx$theta, 0))) {
    expect_error(
      tulpaObs:::cpp_ms_ocs_joint_logpost(fx$spec, bad, fx$pri, 5, 1.0, FALSE),
      "theta length")
  }
})

test_that("cpp_ms_ocs_joint_logpost rejects a theta packed for the other `constrain` (#236)", {
  # The loading block is S*K unconstrained and K*S - K(K-1)/2 constrained, so
  # from K = 2 up the two packings differ in length. Passing a constrained theta
  # with constrain = FALSE is the realistic short read, not a typo.
  free <- .ocs_len_fixture(K = 2L, constrain = FALSE)
  con  <- .ocs_len_fixture(K = 2L, constrain = TRUE)
  lay_free <- tulpaObs:::.ms_ocs_nuts_layout(free$d, FALSE)
  lay_con  <- tulpaObs:::.ms_ocs_nuts_layout(con$d,  TRUE)
  expect_gt(lay_free$total, lay_con$total)

  expect_error(
    tulpaObs:::cpp_ms_ocs_joint_logpost(free$spec, con$theta, free$pri, 5, 1.0,
                                        FALSE),
    "theta length")
  expect_error(
    tulpaObs:::cpp_ms_ocs_joint_logpost(free$spec, free$theta, free$pri, 5, 1.0,
                                        TRUE),
    "theta length")

  # Each packing is accepted under its own flag.
  expect_true(is.finite(tulpaObs:::cpp_ms_ocs_joint_logpost(
    free$spec, free$theta, free$pri, 5, 1.0, FALSE)$lp))
  expect_true(is.finite(tulpaObs:::cpp_ms_ocs_joint_logpost(
    con$spec, con$theta, con$pri, 5, 1.0, TRUE)$lp))
})

test_that("cpp_ms_ocs_nuts rejects a theta0 off the NUTS layout (#236)", {
  # The sampler's guard runs before any draw, so this returns without sampling.
  fx  <- .ocs_len_fixture(K = 1L)
  lay <- tulpaObs:::.ms_ocs_nuts_layout(fx$d, FALSE)
  expect_error(
    tulpaObs:::cpp_ms_ocs_nuts(fx$spec, fx$theta[-lay$total], fx$pri, 5, 1.0,
                               NULL, n_iter = 2L, n_warmup = 1L,
                               max_treedepth = 3L, adapt_delta = 0.8,
                               seed = 1L, verbose = FALSE, constrain = FALSE),
    "theta0 length")
})
