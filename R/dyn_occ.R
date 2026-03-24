#' Define a multi-season dynamic occupancy model
#'
#' @param occ_formula Formula for initial occupancy (season 1)
#' @param det_formula Formula for detection probability
#' @param col_formula Formula for colonization probability (default `~ 1`)
#' @param ext_formula Formula for extinction probability (default `~ 1`)
#' @param data A data frame with site-level covariates
#' @param y 3D array of detection histories `[n_sites x max_visits x n_seasons]`,
#'   or a list of matrices (one per season). Use NA for missing visits.
#'
#' @return A `tulpaOcc_dynmodel` object
#' @export
dynOcc <- function(occ_formula, det_formula, data, y,
                   col_formula = ~ 1, ext_formula = ~ 1) {

  # Handle y as 3D array or list of matrices
  if (is.list(y) && !is.array(y)) {
    n_seasons <- length(y)
    n_sites <- nrow(y[[1]])
    max_visits <- ncol(y[[1]])
    y_array <- array(NA_integer_, dim = c(n_sites, max_visits, n_seasons))
    for (t in seq_len(n_seasons)) {
      y_array[, , t] <- as.integer(y[[t]])
    }
    y <- y_array
  }

  if (length(dim(y)) != 3) {
    stop("y must be a 3D array [n_sites x max_visits x n_seasons] or a list of matrices")
  }

  n_sites <- dim(y)[1]
  max_visits <- dim(y)[2]
  n_seasons <- dim(y)[3]

  if (nrow(data) != n_sites) {
    stop(sprintf("y has %d sites but data has %d rows", n_sites, nrow(data)))
  }

  # Build design matrices
  X_occ <- model.matrix(occ_formula, data)
  X_det <- model.matrix(det_formula, data)
  X_col <- model.matrix(col_formula, data)
  X_ext <- model.matrix(ext_formula, data)

  # Flatten detection history: [site * n_seasons * max_visits]
  y_int <- as.integer(y)
  y_int[is.na(y_int)] <- -1L

  # Compute n_visits and any_detected per site-season
  n_visits <- integer(n_sites * n_seasons)
  any_detected <- logical(n_sites * n_seasons)

  for (i in seq_len(n_sites)) {
    for (t in seq_len(n_seasons)) {
      idx <- (i - 1) * n_seasons + (t - 1)
      season_data <- y[i, , t]
      valid <- !is.na(season_data) & season_data >= 0
      # Use original y (before NA replacement) for valid check
      raw <- y[i, , t]
      raw[is.na(raw)] <- -1L
      valid_visits <- raw >= 0
      n_visits[idx + 1] <- sum(valid_visits)
      any_detected[idx + 1] <- any(raw[valid_visits] == 1)
    }
  }

  structure(
    list(
      y_flat = y_int,
      n_visits = n_visits,
      any_detected = any_detected,
      X_occ = X_occ,
      X_det = X_det,
      X_col = X_col,
      X_ext = X_ext,
      occ_formula = occ_formula,
      det_formula = det_formula,
      col_formula = col_formula,
      ext_formula = ext_formula,
      n_sites = n_sites,
      n_seasons = n_seasons,
      max_visits = max_visits,
      p_occ = ncol(X_occ),
      p_det = ncol(X_det),
      p_col = ncol(X_col),
      p_ext = ncol(X_ext),
      occ_names = colnames(X_occ),
      det_names = colnames(X_det),
      col_names = colnames(X_col),
      ext_names = colnames(X_ext)
    ),
    class = "tulpaOcc_dynmodel"
  )
}

#' @export
print.tulpaOcc_dynmodel <- function(x, ...) {
  cat("Multi-season dynamic occupancy model\n")
  cat(sprintf("  Sites: %d, Seasons: %d, Max visits: %d\n",
              x$n_sites, x$n_seasons, x$max_visits))
  cat(sprintf("  Initial occupancy (%d): %s\n",
              x$p_occ, paste(x$occ_names, collapse = ", ")))
  cat(sprintf("  Detection (%d): %s\n",
              x$p_det, paste(x$det_names, collapse = ", ")))
  cat(sprintf("  Colonization (%d): %s\n",
              x$p_col, paste(x$col_names, collapse = ", ")))
  cat(sprintf("  Extinction (%d): %s\n",
              x$p_ext, paste(x$ext_names, collapse = ", ")))
  invisible(x)
}
