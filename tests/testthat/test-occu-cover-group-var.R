# =============================================================================
# test-occu-cover-group-var.R - occu_cover() with group_var: occupancy units
# (sites) decoupled from field nodes (cells). Many sites share one cell field
# node (site = cell x period), so a per-site trend weight gives an occupancy
# time trend on a shared areal field. The icar()/bym2() term carries
# group_var = "<col>" mapping each site row to a field node.
#
# Covers: the fit runs with n_sites > n_cells, the fields stay length n_cells
# (not n_sites), both coupled fields are exposed, and (recovery) the per-cell
# field shapes + arm coefficients recover the generative truth.
# =============================================================================

.gv_chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (s in seq_len(n)) { if (s > 1L) adj[s, s-1L] <- 1L; if (s < n) adj[s, s+1L] <- 1L }
  adj
}

# Simulate sites = cells x periods sharing a per-cell ICAR field, with a per-site
# time weight driving the trend on psi (and, via alpha_trend, on cover).
.gv_sim <- function(n_cells, n_per, J, adj, seed,
                    b_psi = c(0, 0.6), b_p = c(stats::qlogis(0.55), 0.5),
                    b_pos = c(stats::qlogis(0.3), 0.4),
                    sigma = 0.9, alpha = 1.2, sigma_trend = 0.7,
                    alpha_trend = 0.8, phi = 25) {
  set.seed(seed)
  n_sites <- n_cells * n_per
  Q  <- tulpaObs:::.occu_cover_icar_Q(adj)
  sq <- tulpaObs:::.occu_cover_icar_scale(adj)
  eig <- eigen(Q, symmetric = TRUE); keep <- eig$values > 1e-8
  draw_f <- function() {
    zk <- stats::rnorm(sum(keep))
    fk <- as.numeric(eig$vectors[, keep, drop = FALSE] %*% (zk / sqrt(eig$values[keep])))
    (fk - mean(fk)) / sqrt(sq)
  }
  f <- draw_f(); f2 <- draw_f()
  site <- data.frame(cell_idx = rep(seq_len(n_cells), times = n_per),
                     period   = rep(seq_len(n_per), each = n_cells))
  site$time.sc <- as.numeric(scale(site$period))
  eta_psi <- b_psi[1] + b_psi[2]*site$time.sc +
             sigma*f[site$cell_idx] + sigma_trend*site$time.sc*f2[site$cell_idx]
  z <- stats::rbinom(n_sites, 1L, stats::plogis(eta_psi))
  det_cov <- matrix(stats::rnorm(n_sites*J), n_sites, J)
  pos_cov <- matrix(stats::rnorm(n_sites*J), n_sites, J)
  Y <- matrix(0L, n_sites, J); Ypos <- matrix(0, n_sites, J)
  for (i in seq_len(n_sites)) for (j in seq_len(J)) {
    if (z[i] == 1L) {
      d <- stats::rbinom(1L, 1L, stats::plogis(b_p[1] + b_p[2]*det_cov[i,j]))
      Y[i,j] <- d
      if (d == 1L) {
        mu <- stats::plogis(b_pos[1] + b_pos[2]*pos_cov[i,j] +
                            alpha*sigma*f[site$cell_idx[i]] +
                            alpha_trend*sigma_trend*site$time.sc[i]*f2[site$cell_idx[i]])
        Ypos[i,j] <- stats::rbeta(1L, mu*phi, (1-mu)*phi)
      }
    }
  }
  Ypos[Ypos <= 0] <- 0; Ypos[Ypos >= 1] <- 1 - 1e-6
  list(site = site, Y = Y, Ypos = Ypos,
       vd = data.frame(det_cov = as.vector(t(det_cov)),
                       pos_cov = as.vector(t(pos_cov))),
       f = f - mean(f), f2 = f2 - mean(f2),
       truth = c(b_psi, b_p, b_pos),
       n_sites = n_sites, n_cells = n_cells)
}

.gv_fit <- function(sim, adj, max.iter = 300L) {
  suppressWarnings(tobs(
    formula = ~ time.sc + icar(graph = adj, group_var = "cell_idx") +
                icar(graph = adj, weight = time.sc, group_var = "cell_idx"),
    data = sim$site, family = occu_cover("beta"),
    detection = ~ det_cov, positive = ~ pos_cov,
    y = sim$Y, y_pos = sim$Ypos, visits = sim$vd,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = max.iter, engine = "joint_coupled",
                   sigma.grid = exp(seq(log(0.4), log(2.0), length.out = 4)),
                   alpha.grid = c(0, 0.6, 1.5),
                   alpha.grid.trend = c(0, 0.6, 1.5),
                   adaptive.grid = TRUE)
  ))
}


test_that("occu_cover group_var runs with more sites than field nodes", {
  n_cells <- 16L; n_per <- 4L; J <- 12L
  adj <- .gv_chain_adj(n_cells)
  sim <- .gv_sim(n_cells, n_per, J, adj, seed = 4242L)
  expect_gt(sim$n_sites, sim$n_cells)            # 64 sites, 16 cells

  fit <- .gv_fit(sim, adj)
  expect_s3_class(fit, "tobs_fit")
  expect_identical(attr(fit, "tobs_family")$name, "occu_cover")

  # Fields are sized by the GRAPH (cells), not the occupancy units (sites).
  expect_length(fit$spatial_field, n_cells)
  expect_length(fit$trend_field,   n_cells)
  expect_lt(abs(mean(fit$spatial_field)), 1e-6)
  expect_lt(abs(mean(fit$trend_field)),   1e-6)

  # The site -> cell map is carried on the model and the psi design is per site.
  expect_length(fit$model$site_cell, sim$n_sites)
  expect_equal(nrow(fit$model$X_occ), sim$n_sites)
  expect_identical(fit$model$n_cells, n_cells)

  # All arm intercepts + both coupling scales surface and are finite.
  expect_true(all(c("sigma","alpha","sigma_trend","alpha_trend") %in% names(fit$means)))
  expect_true(all(is.finite(fit$means[c("psi_(Intercept)","p_(Intercept)",
                                         "pos_(Intercept)","psi_time.sc")])))

  # Issue-16 joint vcov spans betas + BOTH per-cell fields (2 * n_cells field
  # columns), one sum-to-zero row per ICAR block.
  p_beta <- sum(vapply(fit$process_info, function(p) length(p$coef_names), 1L))
  expect_equal(dim(fit$joint_vcov), c(p_beta + 2L*n_cells, p_beta + 2L*n_cells))
})


test_that("occu_cover group_var recovers fields and slopes (multi-seed)", {
  skip_on_cran()
  n_cells <- 20L; n_per <- 6L; J <- 15L
  adj <- .gv_chain_adj(n_cells)
  seeds <- c(101L, 202L, 303L, 404L, 505L)

  fcor_i <- fcor_t <- numeric(0)
  psi_slope <- p_slope <- numeric(0)
  for (sd in seeds) {
    sim <- .gv_sim(n_cells, n_per, J, adj, seed = sd)
    fit <- .gv_fit(sim, adj)
    fcor_i <- c(fcor_i, cor(fit$spatial_field, sim$f))
    fcor_t <- c(fcor_t, cor(fit$trend_field,   sim$f2))
    psi_slope <- c(psi_slope, unname(fit$means[["psi_time.sc"]]))
    p_slope   <- c(p_slope,   unname(fit$means[["p_det_cov"]]))
  }

  # The shared per-cell fields reconstruct despite many sites per node.
  expect_gt(median(fcor_i), 0.9)
  expect_gt(median(fcor_t), 0.85)
  # The per-site occupancy time slope and detection slope recover their sign and
  # rough magnitude (truth 0.6 and 0.5).
  expect_gt(median(psi_slope), 0.2)
  expect_gt(median(p_slope),   0.3)
})
