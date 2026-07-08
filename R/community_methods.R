# community_methods.R - predict() / residuals() for the community-occupancy
# families (ms_occu / ms_dyn_occu / ms_int_occu). These share the per-species
# coefficient structure (community mean + BLUP), exposed by fitted(); predict
# reads the per-species coefficients, residuals compare the fitted occupancy
# probability to the per-species ever-detected indicator.

# Per-(site, species) ever-detected indicator (1 if the species was detected at
# least once at the site). For the dynamic family this is the season-1 indicator,
# to pair with the season-1 occupancy psi1 that fitted() exposes.
.tobs_community_ever_detected <- function(model) {
  mt <- model$model_type; ns <- model$n_sites; nsp <- model$n_species
  z <- matrix(0L, ns, nsp)
  if (identical(mt, "ms_occu")) {
    for (s in seq_len(nsp))
      z[, s] <- as.integer(rowSums(model$valid[, , s] & model$y[, , s] == 1L) > 0L)
  } else if (identical(mt, "ms_int_occu")) {
    for (s in seq_len(nsp)) {
      det <- Reduce(`+`, lapply(seq_len(model$n_sources), function(d)
        rowSums(model$valid[[d]][, , s] & model$y[[d]][, , s] == 1L)))
      z[, s] <- as.integer(det > 0L)
    }
  } else {                                    # ms_dyn_occu: season 1
    for (s in seq_len(nsp))
      z[, s] <- as.integer(rowSums(model$valid[, , 1L, s] &
                                     model$y[, , 1L, s] == 1L) > 0L)
  }
  colnames(z) <- model$species_names
  z
}

# Binary residual (deviance / pearson / response) of an observed 0/1 matrix
# against a fitted-probability matrix of the same shape.
.tobs_resid_binary <- function(obs, phat, type) {
  eps <- 1e-10
  switch(type,
    response = obs - phat,
    pearson  = (obs - phat) / sqrt(phat * (1 - phat) + eps),
    deviance = sign(obs - phat) * sqrt(2 * abs(
      ifelse(obs == 1, -log(phat + eps), -log(1 - phat + eps)))))
}

# predict() for the community-occupancy families. In-sample (newdata = NULL)
# returns the per-species fitted probability matrix [n_sites x n_species];
# `newdata` recomputes it from the per-species coefficients (non-spatial only).
.tobs_predict_ms_community <- function(object, newdata = NULL,
                                       type = c("occupancy", "detection")) {
  type  <- match.arg(type)
  model <- object$model; mt <- model$model_type

  if (is.null(newdata)) {
    fv <- fitted(object)
    if (identical(type, "detection")) return(fv$p)
    return(if (identical(mt, "ms_dyn_occu")) fv$psi1 else fv$psi)
  }
  if (!is.null(object$spatial_field)) {
    stop("predict(newdata = ) is not supported for a spatial community fit: the ",
         "shared areal field is tied to the in-sample cells. Call predict() ",
         "without `newdata` for the in-sample per-species posterior.",
         call. = FALSE)
  }
  cm <- object$ms_community
  if (identical(type, "occupancy")) {
    f    <- model$formulas$occ
    coef <- if (identical(mt, "ms_dyn_occu")) cm$coef_psi1 else cm$coef_psi
  } else {
    if (identical(mt, "ms_int_occu")) {
      stop("predict(type = \"detection\", newdata = ) is not defined for the ",
           "multi-source ms_int_occu() (each source carries its own detection ",
           "design). Call predict() without `newdata`, or fitted().",
           call. = FALSE)
    }
    f    <- model$formulas$det
    coef <- cm$coef_p
  }
  X   <- stats::model.matrix(f, newdata)
  out <- stats::plogis(X %*% t(coef))
  colnames(out) <- model$species_names
  out
}

# residuals() for the community-occupancy families: per-species occupancy
# residual of the fitted occupancy probability (smoothed z where available, else
# the season-1 psi1) against the ever-detected indicator. Detection-level
# residuals are per-species and left NULL (the smoothed per-visit state is not
# stored), matching residuals.tobs_fit()'s det = NULL for non-single fits.
.tobs_residuals_ms_community <- function(object, type) {
  model <- object$model
  fv    <- fitted(object)
  occ_p <- if (identical(model$model_type, "ms_dyn_occu")) fv$psi1 else fv$z
  z_obs <- .tobs_community_ever_detected(model)
  list(occ = .tobs_resid_binary(z_obs, occ_p, type), det = NULL)
}
