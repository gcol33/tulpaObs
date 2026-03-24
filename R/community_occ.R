#' Define a community (multi-species) occupancy model
#'
#' @param occ_formula Formula for occupancy probability
#' @param det_formula Formula for detection probability
#' @param data A data frame with site-level covariates
#' @param y 3D array of detection histories `[n_sites x max_visits x n_species]`,
#'   or a named list of matrices (one per species)
#' @param species_names Optional character vector of species names
#'
#' @return A `tulpaOcc_community` object
#' @export
communityOcc <- function(occ_formula, det_formula, data, y,
                         species_names = NULL) {

  # Handle y as 3D array or list of matrices
  if (is.list(y) && !is.array(y)) {
    n_species <- length(y)
    n_sites <- nrow(y[[1]])
    max_visits <- ncol(y[[1]])
    if (is.null(species_names) && !is.null(names(y))) {
      species_names <- names(y)
    }
    y_array <- array(NA_integer_, dim = c(n_sites, max_visits, n_species))
    for (s in seq_len(n_species)) {
      y_array[, , s] <- as.integer(y[[s]])
    }
    y <- y_array
  }

  if (length(dim(y)) != 3) {
    stop("y must be a 3D array [n_sites x max_visits x n_species] or a list of matrices")
  }

  n_sites <- dim(y)[1]
  max_visits <- dim(y)[2]
  n_species <- dim(y)[3]

  if (nrow(data) != n_sites) {
    stop(sprintf("y has %d sites but data has %d rows", n_sites, nrow(data)))
  }

  if (is.null(species_names)) {
    species_names <- paste0("sp", seq_len(n_species))
  }

  # Build design matrices (site-level, replicated per species)
  X_occ <- model.matrix(occ_formula, data)
  X_det <- model.matrix(det_formula, data)

  # Expand to site-species observations: N = n_sites * n_species
  # Order: all species for site 1, then all species for site 2, ...
  N <- n_sites * n_species
  X_occ_expanded <- X_occ[rep(seq_len(n_sites), each = n_species), , drop = FALSE]
  X_det_expanded <- X_det[rep(seq_len(n_sites), each = n_species), , drop = FALSE]

  # Flatten detection history to [N x max_visits]
  # y_expanded[obs, visit] where obs = (site-1)*n_species + (species-1) + 1
  y_expanded <- matrix(NA_integer_, nrow = N, ncol = max_visits)
  for (i in seq_len(n_sites)) {
    for (s in seq_len(n_species)) {
      obs <- (i - 1) * n_species + s
      y_expanded[obs, ] <- as.integer(y[i, , s])
    }
  }
  y_expanded[is.na(y_expanded)] <- -1L

  # Species grouping (1-based, for RE)
  species_group <- rep(seq_len(n_species), times = n_sites)

  structure(
    list(
      y = y_expanded,
      X_occ = X_occ_expanded,
      X_det = X_det_expanded,
      species_group = as.integer(species_group),
      occ_formula = occ_formula,
      det_formula = det_formula,
      n_sites = n_sites,
      n_species = n_species,
      max_visits = max_visits,
      N = N,
      p_occ = ncol(X_occ),
      p_det = ncol(X_det),
      occ_names = colnames(X_occ),
      det_names = colnames(X_det),
      species_names = species_names
    ),
    class = "tulpaOcc_community"
  )
}

#' @export
print.tulpaOcc_community <- function(x, ...) {
  cat("Community occupancy model\n")
  cat(sprintf("  Sites: %d, Species: %d, Max visits: %d\n",
              x$n_sites, x$n_species, x$max_visits))
  cat(sprintf("  Observations (site x species): %d\n", x$N))
  cat(sprintf("  Occupancy covariates (%d): %s\n",
              x$p_occ, paste(x$occ_names, collapse = ", ")))
  cat(sprintf("  Detection covariates (%d): %s\n",
              x$p_det, paste(x$det_names, collapse = ", ")))
  cat(sprintf("  Species RE: intercept on psi and p\n"))
  invisible(x)
}
