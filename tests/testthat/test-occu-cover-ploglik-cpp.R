# Parallel C++ pointwise log-likelihood for the compact (ragged) occu_cover()
# path (cpp_occu_cover_ploglik_ragged, the WAIC / PSIS-LOO input). The R oracle
# is the per-draw ragged marginal (.occu_cover_site_ll_ragged folded through
# .occu_cover_eta_components / .occu_cover_draw_eta_ragged) -- exactly what
# .occu_cover_ploglik_core ran before the port. The C++ kernel accumulates each
# site's visit sums in visit (index) order, matching R's rowsum(), so the two
# agree to libm rounding (~1e-13); and because the draw loop has no shared
# writes, the result is thread-count invariant (a divergence there is a bug, not
# a rounding artefact). Both cover families and both visit-block layouts.

# A synthetic ragged occu_cover model: site + optional visit designs, a random
# detection/cover history, everything .occu_cover_ploglik_core reads directly.
# Bypasses fitting -- the kernel is a pure function of (model, draws, field).
.mk_ragged_model <- function(positive, has_dv, has_pv, n_sites = 40L,
                             mean_visits = 6, seed = 1) {
  set.seed(seed)
  nv   <- pmax(1L, rpois(n_sites, mean_visits))
  site <- rep(seq_len(n_sites), nv); V <- length(site)
  ydet <- rbinom(V, 1L, 0.35)
  ypos <- numeric(V); det <- ydet == 1L
  ypos[det] <- if (identical(positive, "beta")) runif(sum(det), 1e-3, 1 - 1e-3)
               else rlnorm(sum(det))
  structure(list(
    model_type = "occu_cover", positive = positive, ragged = TRUE,
    cover_aggregate = "none", site_of_visit = site, y_det_visit = ydet,
    y_pos_visit = ypos, n_visits_valid = V, n_sites = n_sites,
    max_visits = max(nv),
    X_occ = cbind(1, rnorm(n_sites), rnorm(n_sites)),
    X_det_site = cbind(1, rnorm(n_sites)),
    X_pos_site = cbind(1, rnorm(n_sites), rnorm(n_sites)),
    X_det_visit = if (has_dv) cbind(rnorm(V), rnorm(V)) else NULL,
    X_pos_visit = if (has_pv) cbind(rnorm(V)) else NULL
  ), class = "tobs_model")
}

.ploglik_oracle <- function(model, b_occ, b_det, b_pos, disp,
                            field_occ, field_pos) {
  S  <- nrow(b_occ)
  cl <- .tobs_clamp_eta
  comp <- .occu_cover_eta_components(model, b_occ, b_det, b_pos,
                                     field_occ, field_pos)
  ll <- matrix(0, S, model$n_sites)
  for (j in seq_len(S)) {
    de    <- .occu_cover_draw_eta_ragged(comp, j, model$site_of_visit)
    psi   <- stats::plogis(cl(de$psi_eta))
    p_vec <- stats::plogis(cl(de$p_eta))
    ll[j, ] <- .occu_cover_site_ll_ragged(model, psi, p_vec, de$ep, log(disp[j]))
  }
  ll
}

.mk_draws <- function(model, S = 40L, seed = 7) {
  set.seed(seed)
  p_occ <- ncol(model$X_occ)
  p_det <- ncol(model$X_det_site) +
    (if (!is.null(model$X_det_visit)) ncol(model$X_det_visit) else 0L)
  p_pos <- ncol(model$X_pos_site) +
    (if (!is.null(model$X_pos_visit)) ncol(model$X_pos_visit) else 0L)
  list(
    b_occ = matrix(rnorm(S * p_occ, 0, 0.7), S, p_occ),
    b_det = matrix(rnorm(S * p_det, 0, 0.7), S, p_det),
    b_pos = matrix(rnorm(S * p_pos, 0, 0.7), S, p_pos),
    disp  = if (identical(model$positive, "beta")) runif(S, 2, 20)
            else runif(S, 0.3, 1.5),
    field_occ = matrix(rnorm(model$n_sites * S, 0, 0.5), model$n_sites, S),
    field_pos = matrix(rnorm(model$n_sites * S, 0, 0.5), model$n_sites, S)
  )
}

test_that("C++ ragged pointwise loglik matches the R oracle (both families, visit blocks)", {
  grid <- expand.grid(positive = c("beta", "lognormal"),
                      has_dv = c(FALSE, TRUE), has_pv = c(FALSE, TRUE),
                      stringsAsFactors = FALSE)
  for (r in seq_len(nrow(grid))) {
    m  <- .mk_ragged_model(grid$positive[r], grid$has_dv[r], grid$has_pv[r])
    dr <- .mk_draws(m)
    R  <- .ploglik_oracle(m, dr$b_occ, dr$b_det, dr$b_pos, dr$disp,
                          dr$field_occ, dr$field_pos)
    C  <- .occu_cover_ploglik_core(m, dr$b_occ, dr$b_det, dr$b_pos, dr$disp,
                                   dr$field_occ, dr$field_pos, n_threads = 1L)
    expect_equal(C, R, tolerance = 1e-8,
                 info = sprintf("%s dv=%s pv=%s", grid$positive[r],
                                grid$has_dv[r], grid$has_pv[r]))
  }
})

# A dense (padded [n_sites x max_visits]) no-aggregation occu_cover model, flat
# enough to drive .occu_cover_ploglik_core through the dense->ragged adaptor.
.mk_dense_model <- function(positive, has_dv, has_pv, n_sites = 30L, mv = 8L,
                            seed = 21) {
  set.seed(seed)
  nv <- sample(2:mv, n_sites, replace = TRUE)
  valid <- matrix(FALSE, n_sites, mv)
  for (i in seq_len(n_sites)) valid[i, seq_len(nv[i])] <- TRUE
  y <- matrix(-1L, n_sites, mv); ypos <- matrix(0, n_sites, mv)
  for (i in seq_len(n_sites)) for (v in seq_len(nv[i])) {
    y[i, v] <- rbinom(1L, 1L, 0.35)
    if (y[i, v] == 1L)
      ypos[i, v] <- if (positive == "beta") runif(1, 1e-3, 1 - 1e-3) else rlnorm(1)
  }
  R <- n_sites * mv
  structure(list(
    model_type = "occu_cover", positive = positive, cover_aggregate = "none",
    n_sites = n_sites, max_visits = mv, valid = valid, y = y, y_pos = ypos,
    X_occ = cbind(1, rnorm(n_sites), rnorm(n_sites)),
    X_det_site = cbind(1, rnorm(n_sites)),
    X_pos_site = cbind(1, rnorm(n_sites), rnorm(n_sites)),
    X_det_visit = if (has_dv) cbind(rnorm(R), rnorm(R)) else NULL,
    X_pos_visit = if (has_pv) cbind(rnorm(R)) else NULL
  ), class = "tobs_model")
}

# Dense R oracle: loop draws through .occu_cover_site_ll (the padded-grid path).
.ploglik_oracle_dense <- function(model, b_occ, b_det, b_pos, disp,
                                  field_occ, field_pos) {
  S  <- nrow(b_occ); cl <- .tobs_clamp_eta
  comp <- .occu_cover_eta_components(model, b_occ, b_det, b_pos,
                                     field_occ, field_pos)
  ll <- matrix(0, S, model$n_sites)
  for (j in seq_len(S)) {
    de    <- .occu_cover_draw_eta(comp, j, model$n_sites, model$max_visits)
    psi   <- stats::plogis(cl(de$psi_eta))
    p_mat <- stats::plogis(cl(de$p_eta))
    ll[j, ] <- .occu_cover_site_ll(model, psi, p_mat, de$ep_mat, log(disp[j]),
                                   units = NULL)
  }
  ll
}

test_that("C++ dense (padded) occu_cover pointwise loglik matches the R oracle", {
  grid <- expand.grid(positive = c("beta", "lognormal"),
                      has_dv = c(FALSE, TRUE), has_pv = c(FALSE, TRUE),
                      stringsAsFactors = FALSE)
  for (r in seq_len(nrow(grid))) {
    m  <- .mk_dense_model(grid$positive[r], grid$has_dv[r], grid$has_pv[r])
    dr <- .mk_draws(m)
    R  <- .ploglik_oracle_dense(m, dr$b_occ, dr$b_det, dr$b_pos, dr$disp,
                                dr$field_occ, dr$field_pos)
    C  <- .occu_cover_ploglik_core(m, dr$b_occ, dr$b_det, dr$b_pos, dr$disp,
                                   dr$field_occ, dr$field_pos, n_threads = 1L)
    expect_equal(C, R, tolerance = 1e-8,
                 info = sprintf("%s dv=%s pv=%s", grid$positive[r],
                                grid$has_dv[r], grid$has_pv[r]))
  }
})

test_that("C++ ragged pointwise loglik is thread-count invariant", {
  m  <- .mk_ragged_model("beta", TRUE, TRUE, n_sites = 60L, mean_visits = 10)
  dr <- .mk_draws(m, S = 64L)
  C1 <- .occu_cover_ploglik_core(m, dr$b_occ, dr$b_det, dr$b_pos, dr$disp,
                                 dr$field_occ, dr$field_pos, n_threads = 1L)
  C8 <- .occu_cover_ploglik_core(m, dr$b_occ, dr$b_det, dr$b_pos, dr$disp,
                                 dr$field_occ, dr$field_pos, n_threads = 8L)
  # No shared writes across draws -> identical to the bit, not merely close.
  expect_identical(C8, C1)
})
