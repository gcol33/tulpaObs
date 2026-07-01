# Batched per-(cell, species) pointwise log-likelihood for the spatial-factor
# community occupancy + cover family (ms_occu_cover_spatial). The former R loop
# unpacked the NUTS draw (mu, per-species b, shared fields W, occupancy loadings
# L, optional cover loadings Lpos, log-dispersion), assembled each species'
# predictors with the shared-factor offset W L[s,] (and W Lpos[s,] on cover), and
# evaluated the dense occu_cover per-cell marginal. The C++ kernel does all of
# that, parallel over draws, mirroring the R oracle (~1e-15) and thread invariant;
# both cover families and both cover-factor settings.

.mk_ms_ocs <- function(positive, cf, N = 10L, J = 2L, S = 3L, K = 2L, seed = 111) {
  set.seed(seed)
  X_det_visit <- cbind(rnorm(N * J))
  y <- array(0L, c(N, J, S)); ypos <- array(0, c(N, J, S)); valid <- array(FALSE, c(N, J, S))
  for (s in seq_len(S)) for (c in seq_len(N)) for (j in seq_len(J)) {
    valid[c, j, s] <- runif(1) < 0.85
    if (valid[c, j, s]) { y[c, j, s] <- rbinom(1L, 1L, 0.3)
      if (y[c, j, s] == 1L)
        ypos[c, j, s] <- if (positive == "beta") runif(1, 1e-3, 1 - 1e-3) else rlnorm(1) }
  }
  model <- list(model_type = "ms_occu_cover_spatial", positive = positive,
                cover_aggregate = "none", n_sites = N, n_species = S, max_visits = J,
                K = K, cover_factor = cf, X_occ = cbind(1, rnorm(N)),
                X_det_site = cbind(1, rnorm(N)), X_det_visit = X_det_visit,
                X_pos_site = cbind(1, rnorm(N)), X_pos_visit = NULL,
                y = y, y_pos = ypos, valid = valid,
                process_info = list(list(p = 2L), list(p = 3L), list(p = 2L)))
  d <- .ms_ocs_dims(model)
  list(model = model, d = d,
       total = d$P + d$S * d$P + d$S * d$K +
               (if (cf) d$S * d$K else 0L) + d$N * d$K + 1L)
}

.ms_ocs_oracle <- function(model, d, draws) {
  cl <- .tobs_clamp_eta; M <- nrow(draws); out <- matrix(0, M, d$N * d$S)
  for (i in seq_len(M)) {
    up <- .ms_ocs_unpack(draws[i, ], d); LL <- matrix(0, d$N, d$S)
    for (s in seq_len(d$S)) {
      v <- .ms_occu_cover_species_view(model, s); th <- up$mu + up$b[[s]]
      eta <- .occu_cover_eta_from_par(v, th[d$occ_idx], th[d$p_idx], th[d$pos_idx])
      eta$psi <- stats::plogis(cl(as.numeric(v$X_occ %*% th[d$occ_idx]) +
                                    as.numeric(up$W %*% up$L[s, ])))
      if (d$cover_factor) eta$ep_mat <- eta$ep_mat + as.numeric(up$W %*% up$Lpos[s, ])
      LL[, s] <- .occu_cover_site_ll(v, eta$psi, eta$p_mat, eta$ep_mat, up$ld)
    }
    out[i, ] <- as.numeric(LL)
  }
  out
}

.ms_ocs_new <- function(model, d, draws, nt) {
  cpp_ms_ocs_ploglik(draws, model$X_occ, model$X_det_site,
    if (is.null(model$X_det_visit)) matrix(0, d$N * model$max_visits, 0L) else model$X_det_visit,
    model$X_pos_site, matrix(0, d$N * model$max_visits, 0L),
    as.integer(model$y), as.numeric(model$y_pos), as.integer(model$valid),
    d$N, model$max_visits, d$S, d$K, d$P_occ, d$P_p, d$P_pos,
    ncol(model$X_det_site), ncol(model$X_pos_site), d$cover_factor,
    identical(model$positive, "beta"), nt)
}

test_that("spatial-factor community occu+cover ploglik: C++ == R oracle", {
  for (positive in c("beta", "lognormal")) for (cf in c(FALSE, TRUE)) {
    x <- .mk_ms_ocs(positive, cf); M <- 20L
    draws <- matrix(rnorm(M * x$total, 0, 0.4), M, x$total)
    R <- .ms_ocs_oracle(x$model, x$d, draws)
    C1 <- .ms_ocs_new(x$model, x$d, draws, 1L)
    expect_equal(C1, R, tolerance = 1e-8, info = sprintf("%s cf=%s", positive, cf))
    expect_identical(.ms_ocs_new(x$model, x$d, draws, 4L), C1)
  }
})
