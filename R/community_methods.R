# community_methods.R - predict() / residuals() for the community families.
# They share the per-species coefficient structure (community mean + BLUP),
# exposed by fitted(); predict reads the per-species coefficients, residuals
# compare the fitted occupancy probability to the per-species ever-detected
# indicator.

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

# ---------------------------------------------------------------------------
# predict() for the community families
# ---------------------------------------------------------------------------

# The response arms of every community family: for each arm the `model$formulas`
# slot its design is built from, the `ms_community` coefficient block, the
# `fitted()` slot it is reported in, and the inverse link taking eta to the
# response scale. A family gains a real `newdata` predictor -- and with it the
# default response type `.TOBS_PREDICT_NEWDATA_TYPE` reports and the `terms`
# refusal `.TOBS_PREDICT_NO_TERMS` carries -- by gaining a row here. The FIRST
# arm of a family is the type it reports when the caller names none.
#
# A `product` arm is the elementwise product of two other arms of the same
# family rather than a design of its own; a `newdata_reason` arm is one whose
# design cannot be rebuilt from a data frame, and says why.
.TOBS_MS_PREDICT_ARMS <- list(
  ms_occu = list(
    occupancy = list(formula = "occ", coef = "coef_psi",  fitted = "psi",
                     inv = "logit"),
    detection = list(formula = "det", coef = "coef_p",    fitted = "p",
                     inv = "logit")),
  ms_dyn_occu = list(
    occupancy = list(formula = "occ", coef = "coef_psi1", fitted = "psi1",
                     inv = "logit"),
    detection = list(formula = "det", coef = "coef_p",    fitted = "p",
                     inv = "logit")),
  ms_int_occu = list(
    occupancy = list(formula = "occ", coef = "coef_psi",  fitted = "psi",
                     inv = "logit"),
    detection = list(coef = "coef_p", fitted = "p", inv = "logit",
                     newdata_reason = paste0(
                       "the multi-source ms_int_occu() carries one detection ",
                       "design per source"))),
  # The abundance arm is log-linked: reaching the occupancy fallback below
  # reported plogis() of a log lambda as an occupancy probability
  # (gcol33/tulpaObs#256).
  ms_nmix = list(
    abundance = list(formula = "lambda", coef = "coef_lambda",
                     fitted = "lambda", inv = "log"),
    detection = list(formula = "det",    coef = "coef_p",
                     fitted = "p",      inv = "logit")),
  ms_distance = list(
    lambda = list(formula = "lambda", coef = "coef_lambda", fitted = "lambda",
                  inv = "log"),
    sigma  = list(formula = "sigma",  coef = "coef_sigma",  fitted = "sigma",
                  inv = "log")),
  # Same response vocabulary as the spatial-factor twin
  # (.tobs_ms_ocs_predict_state), so the two halves of one family answer to the
  # same `type`.
  ms_occu_cover = list(
    occupancy  = list(formula = "occ", coef = "coef_occ", fitted = "psi",
                      inv = "logit"),
    detection  = list(formula = "det", coef = "coef_p",   fitted = "p",
                      inv = "logit"),
    cover_cond = list(formula = "pos", coef = "coef_pos", fitted = "cover",
                      inv = "cover"),
    cover_exp  = list(product = c("occupancy", "cover_cond"))))

# One handler serves every family in the table above.
.TOBS_MS_PREDICT_ALIAS <- stats::setNames(
  rep("ms_community", length(.TOBS_MS_PREDICT_ARMS)),
  names(.TOBS_MS_PREDICT_ARMS))

# `type` spellings accepted beside a family's own arm names, applied only when
# the family carries the arm they resolve to.
.TOBS_MS_PREDICT_TYPE_ALIAS <- c(state = "occupancy", occurrence = "occupancy",
                                 cover = "cover_cond",
                                 abundance = "lambda", lambda = "abundance")

# Model slots holding an in-sample latent contribution to a community linear
# predictor (an areal field or a latent-factor loading, both indexed by the
# fitted sites). A fit carrying one cannot be predicted at new rows: there is no
# field at an unseen site.
.TOBS_MS_FIELD_SLOTS <- c("occu_field_offset", "occu_factor_offset",
                          "nmix_field_offset", "nmix_factor_offset",
                          "distance_field_offset", "distance_factor_offset")

.tobs_ms_field_bound <- function(object) {
  if (!is.null(object$spatial_field)) return("shared areal field")
  hit <- .TOBS_MS_FIELD_SLOTS[vapply(.TOBS_MS_FIELD_SLOTS, function(k)
    !is.null(object$model[[k]]), logical(1))]
  if (length(hit)) hit[[1L]] else NULL
}

# Resolve a requested response type against one family's arm table.
.tobs_ms_predict_type <- function(type, arms, mt) {
  if (is.null(type) || length(type) != 1L || is.na(type)) return(names(arms)[[1L]])
  if (type %in% names(arms)) return(type)
  alias <- unname(.TOBS_MS_PREDICT_TYPE_ALIAS[type])
  if (!is.na(alias) && alias %in% names(arms)) return(alias)
  stop(sprintf("predict(type = \"%s\") is not a response of a %s() fit; use %s.",
               type, mt, paste(sprintf('"%s"', names(arms)), collapse = ", ")),
       call. = FALSE)
}

# The arm's design at `newdata` and the coefficient columns it multiplies,
# matched by column name where the block carries names and by leading width
# otherwise.
.tobs_ms_arm_design <- function(f, newdata, coef, arm, mt) {
  X  <- stats::model.matrix(f, newdata)
  cn <- colnames(coef)
  # A registered arm formula is the arm's SITE-level design. A detection or
  # cover block may carry visit-level columns beyond it, which are not
  # site-summarised; the arm reports the site-level predictor at them, the same
  # columns fitted() reads for the same rows.
  if (!is.null(cn) && all(colnames(X) %in% cn))
    return(list(X = X, coef = coef[, colnames(X), drop = FALSE]))
  if (ncol(coef) >= ncol(X))
    return(list(X = X, coef = coef[, seq_len(ncol(X)), drop = FALSE]))
  stop(sprintf(paste0("predict(newdata = , type = \"%s\") on a %s() fit built a ",
                      "%d-column design for a %d-column coefficient block. ",
                      "`newdata` needs the columns the arm was fit on: %s."),
               arm, mt, ncol(X), ncol(coef),
               paste(cn %||% "(unnamed)", collapse = ", ")), call. = FALSE)
}

# eta -> response. The cover arm's link is the fit's own positive family, read
# through the same helper fitted() reads it through.
.tobs_ms_inv_link <- function(eta, inv, object) {
  switch(inv,
         logit = stats::plogis(eta),
         log   = exp(eta),
         cover = .tobs_ms_cover_response(eta, object),
         stop(sprintf("unregistered inverse link '%s'", inv), call. = FALSE))
}

# predict() for the community families. In-sample (newdata = NULL) returns the
# per-species fitted matrix [n_sites x n_species] of the requested arm;
# `newdata` recomputes it from the per-species coefficients (non-spatial only).
.tobs_predict_ms_community <- function(object, newdata = NULL, type = NULL) {
  model <- object$model
  mt    <- model$model_type %||% ""
  arms  <- .TOBS_MS_PREDICT_ARMS[[mt]]
  if (is.null(arms)) {
    stop(sprintf(paste0("predict() has no registered response arms for model ",
                        "type '%s'. Register it in .TOBS_MS_PREDICT_ARMS."), mt),
         call. = FALSE)
  }
  type <- .tobs_ms_predict_type(type, arms, mt)

  if (is.null(newdata)) {
    fv <- fitted(object)
    insample <- function(k) {
      a <- arms[[k]]
      if (!is.null(a$product)) return(Reduce(`*`, lapply(a$product, insample)))
      fv[[a$fitted]]
    }
    return(insample(type))
  }

  bound <- .tobs_ms_field_bound(object)
  if (!is.null(bound)) {
    stop("predict(newdata = ) is not supported for this community fit: its ",
         bound, " is tied to the in-sample sites, so there is no value for it ",
         "at a new row. Call predict() without `newdata` for the in-sample ",
         "per-species posterior.", call. = FALSE)
  }

  cm <- object$ms_community
  arm_at <- function(k) {
    a <- arms[[k]]
    if (!is.null(a$product)) return(Reduce(`*`, lapply(a$product, arm_at)))
    if (!is.null(a$newdata_reason)) {
      stop(sprintf(paste0("predict(type = \"%s\", newdata = ) is not defined ",
                          "for a %s() fit: %s. Call predict() without ",
                          "`newdata`, or fitted()."),
                   k, mt, a$newdata_reason), call. = FALSE)
    }
    d <- .tobs_ms_arm_design(model$formulas[[a$formula]], newdata,
                             cm[[a$coef]], k, mt)
    .tobs_ms_inv_link(d$X %*% t(d$coef), a$inv, object)
  }
  out <- arm_at(type)
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
