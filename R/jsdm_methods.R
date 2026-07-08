# jsdm_methods.R - fitted() / predict() / residuals() for the joint species
# distribution model. jsdm observes presence / absence directly (no detection
# process): the fitted quantity is the per-species occupancy probability psi,
# and the residual is the binary presence residual against it.
#
# eta_{i,s} = X_i beta + field_i + b_s. The shared fixed effects always enter; a
# shared areal field (spatial jsdm) and per-species intercept BLUPs enter when
# the fit surfaces them (object$spatial_field / object$jsdm_re$blup). The
# non-spatial Laplace jsdm integrates the per-species intercept out without
# storing BLUPs, so there psi is the community-level occupancy (identical across
# species except through covariates), consistent with the WAIC scored on the
# fixed-effect predictor (.tobs_ploglik_jsdm).

# Per-species occupancy eta at a site-level design X [n x p]: returns [n x S].
.tobs_jsdm_eta <- function(object, X) {
  model <- object$model
  beta  <- object$means[seq_len(model$process_info[[1L]]$p)]
  eta   <- matrix(as.numeric(X %*% beta), nrow(X), model$n_species)
  blup  <- object$jsdm_re$blup
  if (!is.null(blup)) eta <- sweep(eta, 2L, as.numeric(blup), "+")
  eta
}

.tobs_fitted_jsdm <- function(object) {
  model <- object$model
  eta   <- .tobs_jsdm_eta(object, model$X_occ)
  if (!is.null(object$spatial_field))
    eta <- sweep(eta, 1L, as.numeric(object$spatial_field), "+")
  psi <- stats::plogis(eta)
  colnames(psi) <- model$species_names
  list(psi = psi)
}

# predict(): in-sample returns fitted()$psi; newdata recomputes psi from the
# occupancy design (non-spatial - a shared field cannot be evaluated off-grid).
.tobs_predict_jsdm <- function(object, newdata = NULL) {
  if (is.null(newdata)) return(fitted(object)$psi)
  if (!is.null(object$spatial_field)) {
    stop("predict(newdata = ) is not supported for a spatial jsdm() fit: the ",
         "shared areal field is tied to the in-sample cells. Call predict() ",
         "without `newdata`.", call. = FALSE)
  }
  X   <- stats::model.matrix(object$model$formulas$occ, newdata)
  psi <- stats::plogis(.tobs_jsdm_eta(object, X))
  colnames(psi) <- object$model$species_names
  psi
}

# residuals(): binary presence residual y - psi, per (site, species).
.tobs_residuals_jsdm <- function(object, type) {
  psi <- fitted(object)$psi
  list(occ = .tobs_resid_binary(object$model$y_mat, psi, type), det = NULL)
}
