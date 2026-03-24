#' Define a single-season occupancy model
#'
#' @param occ_formula Formula for occupancy probability (e.g., `~ elevation + forest`)
#' @param det_formula Formula for detection probability (e.g., `~ effort`)
#' @param data A data frame with site-level covariates
#' @param y Detection history matrix (n_sites x max_visits). Use NA for missing visits.
#' @param det_visit_formula Optional formula for visit-level detection covariates
#' @param det_visit_data Optional data frame with visit-level covariates
#'   (n_sites * max_visits rows, one per site-visit combination)
#'
#' @return A `tulpaOcc_model` object
#' @export
occ <- function(occ_formula, det_formula, data, y,
                det_visit_formula = NULL, det_visit_data = NULL) {

  if (!is.matrix(y)) {
    stop("y must be a matrix (n_sites x max_visits)")
  }
  if (nrow(y) != nrow(data)) {
    stop(sprintf("y has %d rows but data has %d rows", nrow(y), nrow(data)))
  }

  # Build design matrices
  X_occ <- model.matrix(occ_formula, data)
  X_det <- model.matrix(det_formula, data)

  # Handle visit-level detection covariates
  X_det_visit <- NULL
  if (!is.null(det_visit_formula) && !is.null(det_visit_data)) {
    X_det_visit <- model.matrix(det_visit_formula, det_visit_data)
    expected_rows <- nrow(y) * ncol(y)
    if (nrow(X_det_visit) != expected_rows) {
      stop(sprintf("det_visit_data must have %d rows (n_sites * max_visits), got %d",
                   expected_rows, nrow(X_det_visit)))
    }
  }

  # Replace NA with -1 in detection history (C++ convention)
  y_int <- matrix(as.integer(y), nrow = nrow(y), ncol = ncol(y))
  y_int[is.na(y_int)] <- -1L

  structure(
    list(
      y = y_int,
      X_occ = X_occ,
      X_det = X_det,
      X_det_visit = X_det_visit,
      occ_formula = occ_formula,
      det_formula = det_formula,
      n_sites = nrow(y),
      max_visits = ncol(y),
      p_occ = ncol(X_occ),
      p_det = ncol(X_det),
      p_det_visit = if (!is.null(X_det_visit)) ncol(X_det_visit) else 0L,
      occ_names = colnames(X_occ),
      det_names = colnames(X_det),
      det_visit_names = if (!is.null(X_det_visit)) colnames(X_det_visit) else character(0)
    ),
    class = "tulpaOcc_model"
  )
}

#' @export
print.tulpaOcc_model <- function(x, ...) {
  cat("Single-season occupancy model\n")
  cat(sprintf("  Sites: %d, Max visits: %d\n", x$n_sites, x$max_visits))
  cat(sprintf("  Occupancy covariates (%d): %s\n",
              x$p_occ, paste(x$occ_names, collapse = ", ")))
  cat(sprintf("  Detection covariates (%d): %s\n",
              x$p_det, paste(x$det_names, collapse = ", ")))
  if (x$p_det_visit > 0) {
    cat(sprintf("  Visit-level detection (%d): %s\n",
                x$p_det_visit, paste(x$det_visit_names, collapse = ", ")))
  }
  n_detected <- sum(apply(x$y, 1, function(row) any(row[row >= 0] == 1)))
  cat(sprintf("  Naive occupancy: %.1f%%\n", 100 * n_detected / x$n_sites))
  invisible(x)
}
