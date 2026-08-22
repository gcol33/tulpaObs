# The block-offset walk shared by the community NUTS layouts. All six
# pack to one recipe --
#
#   mu [P] | {z_s} species-major [S*P] | per-arm chol blocks | trailing blocks
#
# -- and each used to write that arithmetic out itself, including four private
# spellings of the triangular count .ms_ocs_chol_dim() already had. These pin
# the shared walk and the family aliases layered on it; each family's own
# byte-exactness oracle test is what pins its layout against the C++ twin.

test_that("the walk lays mu, the species blocks, the chols and the trailing", {
  lay <- .ms_ocs_layout(list(list(name = "a", width = 2L),
                             list(name = "b", width = 3L)),
                        n_species = 4L,
                        trailing = list(list(name = "g", size = 2L)))
  expect_identical(lay$P, 5L)
  expect_identical(lay$mu, 1:5)
  expect_identical(lay$b_off, 5L)
  expect_identical(lay$idx$a, 1:2)
  expect_identical(lay$idx$b, 3:5)
  # chol section starts after mu and the four species blocks of width 5.
  expect_identical(lay$chol$a, 5L + 20L + seq_len(3L))
  expect_identical(lay$chol$b, 5L + 20L + 3L + seq_len(6L))
  expect_identical(lay$trailing$g, 5L + 20L + 9L + seq_len(2L))
  expect_identical(lay$total, 36L)
  # Nothing overlaps and nothing is skipped.
  all_idx <- c(lay$mu,
               unlist(lapply(seq_len(4L), function(s) .ms_ocs_b_idx(lay, s))),
               lay$chol$a, lay$chol$b, lay$trailing$g)
  expect_identical(sort(all_idx), seq_len(lay$total))
})

test_that("the arms declaration is what .ms_ocs_b_from_z reads", {
  lay <- .ms_ocs_layout(list(list(name = "a", width = 2L),
                             list(name = "b", width = 1L)), n_species = 3L)
  expect_length(lay$arms, 2L)
  expect_identical(lay$arms[[1L]], .ms_ocs_arm(lay$idx$a, lay$chol$a, 2L))
  expect_identical(lay$arms[[2L]], .ms_ocs_arm(lay$idx$b, lay$chol$b, 1L))
  # A one-dimensional arm needs no branch: its 1 x 1 log-Cholesky is exp().
  theta <- numeric(lay$total)
  theta[lay$chol$b] <- log(2)
  theta[.ms_ocs_b_idx(lay, 2L)[3L]] <- 1.5
  B <- .ms_ocs_b_from_z(theta, lay)
  expect_identical(dim(B), c(3L, 3L))
  expect_equal(B[2L, 3L], 2 * 1.5)
})

test_that("the triangular count comes from .ms_ocs_chol_dim alone", {
  for (p in 0:6) expect_identical(.ms_ocs_chol_dim(p), as.integer(p * (p + 1) / 2))
  lay <- .ms_ocs_layout(list(list(name = "a", width = 4L)), n_species = 2L)
  expect_identical(lay$q[["a"]], .ms_ocs_chol_dim(4L))
  # No community NUTS layout writes the count itself any more.
  layouts <- list(.tobs_ms_occu_nuts_layout, .tobs_ms_int_occu_nuts_layout,
                  .tobs_ms_occu_cover_nuts_layout, .tobs_ms_abun_nuts_layout,
                  .tobs_ms_dyn_occu_nuts_layout, .tobs_ms_count_nuts_layout)
  for (f in layouts) {
    src <- paste(deparse(body(f)), collapse = " ")
    expect_false(grepl("+ 1L)/2L", gsub(" ", "", src), fixed = TRUE))
    expect_true(grepl(".ms_ocs_layout(", src, fixed = TRUE))
  }
})

test_that("a zero-size trailing block claims no coordinate", {
  # ms_count's gaussian log_phi block is absent under poisson / negbin.
  lay <- .ms_ocs_layout(list(list(name = "a", width = 2L)), n_species = 3L,
                        trailing = list(list(name = "phi", size = 0L)))
  expect_identical(lay$trailing$phi, integer(0))
  expect_identical(lay$total,
                   .ms_ocs_layout(list(list(name = "a", width = 2L)),
                                  n_species = 3L)$total)
})

test_that("each family exposes the names its own logpost reads", {
  S <- 3L
  lo <- .tobs_ms_occu_nuts_layout(2L, 1L, S)
  expect_identical(lo$psi, 1:2); expect_identical(lo$p, 3L)
  expect_identical(lo$q_psi, 3L); expect_identical(lo$q_p, 1L)

  lc <- .tobs_ms_occu_cover_nuts_layout(2L, 1L, 2L, S)
  expect_identical(lc$occ, 1:2); expect_identical(lc$pos, 4:5)
  expect_identical(lc$log_disp, lc$total)

  ld <- .tobs_ms_dyn_occu_nuts_layout(2L, 1L, 2L, 3L, S)
  expect_identical(ld$G, 5L)
  expect_identical(ld$global, ld$global_off + 1:5)
  expect_identical(ld$total, ld$global_off + ld$G)
  expect_identical(ld$gam, 1:2); expect_identical(ld$eps, 3:5)

  li <- .tobs_ms_int_occu_nuts_layout(2L, c(1L, 3L), S)
  expect_identical(li$D, 2L)
  expect_identical(li$P, 6L)
  expect_length(li$p, 2L); expect_length(li$chol_p, 2L)
  expect_identical(li$p[[1L]], 3L); expect_identical(li$p[[2L]], 4:6)
  expect_identical(li$q_p, c(1L, 6L))

  la <- .tobs_ms_abun_nuts_layout(2L, 1L, S, TRUE)
  expect_identical(la$P, 4L)
  expect_identical(la$logr, 4L)
  expect_identical(la$q_logr, 1L)
  la0 <- .tobs_ms_abun_nuts_layout(2L, 1L, S, FALSE)
  expect_identical(la0$logr, integer(0))
  expect_identical(la0$chol_logr, integer(0))
  expect_identical(la0$q_logr, 0L)

  lg <- .tobs_ms_count_nuts_layout(2L, S, "gaussian")
  expect_identical(lg$chol, lg$chol_beta)   # the Poisson cross-check reads this
  expect_length(lg$logphi, S)
  expect_identical(lg$total, max(lg$logphi))
  expect_identical(.tobs_ms_count_nuts_layout(2L, S, "poisson")$logphi, integer(0))
})
