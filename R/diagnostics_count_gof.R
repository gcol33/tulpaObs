# diagnostics_count_gof.R - dispersion / zero-inflation / outlier goodness-of-fit
# for the count families (nmix / removal / distance / dyn_abun). The single-
# season versions use detection-history semantics (0/1 per visit); these use the
# per-site TOTAL count, the natural overdispersion / excess-zero unit for an
# N-mixture-type model. Each compares the observed statistic to its posterior
# predictive distribution (simulate() replicates) and returns a tail p-value,
# mirroring the single-season tobs_test_* return shape.

.tobs_count_gof_families <- c("nmix", "removal", "distance", "dyn_abun")

# Per-site total counts: observed vector [n_sites] and simulated [n_sites x nsim].
# nmix / removal store long-form counts (y_long, site_idx); distance stores an
# [n_sites x n_bins] matrix; dyn_abun an [n_sites x visits x seasons] array.
# simulate() returns replicates in the family's native shape, reduced to a
# per-site total the same way.
.tobs_count_gof_totals <- function(object, n.samples) {
  model <- object$model; mt <- model$model_type; n_sites <- model$n_sites
  obs <- switch(mt,
    nmix    = ,
    removal = as.numeric(tapply(model$y_long,
                factor(model$site_idx, levels = seq_len(n_sites)), sum)),
    distance = rowSums(model$y, na.rm = TRUE),
    dyn_abun = apply(model$y, 1L, function(a) sum(a, na.rm = TRUE)),
    stop("Count goodness-of-fit is not defined for model_type = '", mt, "'.",
         call. = FALSE))
  obs[is.na(obs)] <- 0

  sims <- simulate(object, nsim = n.samples)
  if (n.samples == 1L) sims <- list(sims)
  site_total <- function(a) {
    if (length(dim(a)) <= 2L) rowSums(a, na.rm = TRUE)
    else apply(a, 1L, function(x) sum(x, na.rm = TRUE))
  }
  list(obs = as.numeric(obs),
       sim = vapply(sims, site_total, numeric(n_sites)))
}

.tobs_test_dispersion_count <- function(object, n.samples) {
  tt <- .tobs_count_gof_totals(object, n.samples)
  obs_var <- stats::var(tt$obs)
  sim_var <- apply(tt$sim, 2L, stats::var)
  list(observed = obs_var, expected = mean(sim_var),
       ratio = obs_var / mean(sim_var), p.value = mean(sim_var >= obs_var))
}

.tobs_test_zero_inflation_count <- function(object, n.samples) {
  tt <- .tobs_count_gof_totals(object, n.samples)
  obs <- sum(tt$obs == 0); sim <- colSums(tt$sim == 0)
  list(observed = obs, expected = mean(sim),
       ratio = obs / max(mean(sim), 1), p.value = mean(sim >= obs))
}

.tobs_test_outliers_count <- function(object, n.samples) {
  tt <- .tobs_count_gof_totals(object, n.samples)
  lower <- apply(tt$sim, 1L, stats::quantile, 0.025)
  upper <- apply(tt$sim, 1L, stats::quantile, 0.975)
  n_out   <- sum(tt$obs < lower | tt$obs > upper)
  sim_out <- colSums(tt$sim < lower | tt$sim > upper)
  list(observed = n_out, expected = mean(sim_out),
       ratio = n_out / max(mean(sim_out), 1), p.value = mean(sim_out >= n_out))
}
