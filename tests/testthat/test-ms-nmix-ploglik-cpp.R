# Batched pointwise log-likelihood for the community N-mixture (ms_abun) family.
# The former R loop reconstructed each species' deviation b = C z from the
# non-centered NUTS draw and called the per-species Royle marginal; the C++
# kernel does the log-Cholesky reconstruction and the per-(species, site)
# marginal (compute_nmix_site) internally, parallel over draws. Because both call
# the identical per-site kernel, the batched result is BYTE-IDENTICAL to the R
# loop (the oracle here, via .tobs_ms_abun_nuts_b_from_z) and thread invariant.

.mk_ms_nmix <- function(is_nb, n_species = 4L, n_sites = 12L, p_lam = 2L,
                        p_p = 2L, M = 20L, seed = 101) {
  set.seed(seed)
  visits <- 2L; rows <- list()
  for (s in seq_len(n_species)) for (si in seq_len(n_sites)) for (v in seq_len(visits))
    rows[[length(rows) + 1L]] <- c(s, si)
  mat <- do.call(rbind, rows); n_obs <- nrow(mat)
  lf <- list(species_idx = mat[, 1L], site_idx = mat[, 2L],
             y = rpois(n_obs, 1.5), X_p = cbind(1, rnorm(n_obs)))
  lay <- .tobs_ms_abun_nuts_layout(p_lam, p_p, n_species, is_nb)
  list(lf = lf, X_lambda = cbind(1, rnorm(n_sites)), lay = lay,
       draws = matrix(rnorm(M * lay$total, 0, 0.4), M, lay$total),
       n_sites = n_sites, is_nb = is_nb, K_max = as.integer(max(lf$y) + 100L))
}

.ms_nmix_oracle <- function(d) {
  lay <- d$lay; S <- lay$n_species; M <- nrow(d$draws)
  margs <- .tobs_ms_abun_nuts_marginals(d$lf, d$X_lambda, d$n_sites,
                                        if (d$is_nb) "NB" else "P", d$K_max)
  out <- matrix(0, M, S * d$n_sites)
  for (m in seq_len(M)) {
    mu <- d$draws[m, lay$mu]; B <- .tobs_ms_abun_nuts_b_from_z(d$draws[m, ], lay)
    for (s in seq_len(S)) {
      b_s <- B[s, ]; r <- if (d$is_nb) exp(mu[lay$logr] + b_s[lay$logr]) else Inf
      out[m, (s - 1L) * d$n_sites + seq_len(d$n_sites)] <-
        margs[[s]]$eval_beta(mu[lay$lambda] + b_s[lay$lambda],
                             mu[lay$p] + b_s[lay$p], r = r)$log_lik_site
    }
  }
  out
}

.ms_nmix_new <- function(d, nt) {
  lay <- d$lay; clogr <- if (d$is_nb) as.integer(lay$chol_logr[1L]) - 1L else 0L
  cpp_ms_nmix_ploglik_batch(as.integer(d$lf$y), as.integer(d$lf$species_idx),
    as.integer(d$lf$site_idx), d$lf$X_p, d$X_lambda, d$draws,
    as.integer(lay$mu[1L]) - 1L, as.integer(lay$b_off),
    as.integer(lay$chol_lam[1L]) - 1L, as.integer(lay$chol_p[1L]) - 1L, clogr,
    as.integer(lay$p_lam), as.integer(lay$p_p), as.integer(lay$n_species),
    as.integer(d$n_sites), d$is_nb, as.integer(d$K_max), nt)
}

test_that("community N-mixture batched ploglik == R reconstruction loop", {
  for (nb in c(FALSE, TRUE)) {
    d <- .mk_ms_nmix(nb)
    R <- .ms_nmix_oracle(d); C1 <- .ms_nmix_new(d, 1L)
    expect_equal(C1, R, tolerance = 1e-10, info = if (nb) "NB" else "P")
    expect_identical(.ms_nmix_new(d, 4L), C1)
  }
})
