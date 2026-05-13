#' Define an occupancy model
#'
#' Unified entry point for single-season, dynamic, and community occupancy models.
#' The model type is inferred from arguments:
#' - **Single-season**: no `col_formula`, no `species`
#' - **Dynamic**: `col_formula` and/or `ext_formula` provided
#' - **Community**: `species` provided (3D array or named list for y)
#'
#' @param occ_formula Formula for occupancy probability (e.g., `~ elevation + forest`).
#'   For dynamic models this is the initial occupancy formula.
#' @param det_formula Formula for detection probability (e.g., `~ effort`).
#' @param data A data frame with site-level covariates.
#' @param y Detection history. One of:
#'   - **Matrix** (n_sites x max_visits): single-season
#'   - **3D array** or list of matrices: dynamic (n_sites x max_visits x n_seasons)
#'     or community (n_sites x max_visits x n_species)
#' @param col_formula Formula for colonization (dynamic models). Default NULL.
#' @param ext_formula Formula for extinction (dynamic models). Default NULL.
#' @param species Character vector of species names, or TRUE to auto-name.
#'   Signals a community model. Default NULL.
#' @param det_visit_formula Optional formula for visit-level detection covariates.
#' @param det_visit_data Optional data frame with visit-level covariates.
#'
#' @return A `TulpaObs` object
#' @export
occu <- function(occ_formula, det_formula = NULL, data, y,
                 col_formula = NULL, ext_formula = NULL,
                 species = NULL, integrated = FALSE, jsdm = FALSE,
                 det_visit_formula = NULL, det_visit_data = NULL) {

  # ---- Determine model type ----
  is_dynamic <- !is.null(col_formula) || !is.null(ext_formula)
  is_community <- !is.null(species) && !isTRUE(jsdm)
  is_integrated <- isTRUE(integrated)
  is_jsdm <- isTRUE(jsdm)

  if (is_dynamic && is_community) {
    stop("Dynamic community models are not yet supported. ",
         "Use col_formula/ext_formula OR species, not both.")
  }

  if (is_jsdm) {
    return(build_jsdm_model(occ_formula, data, y, species))
  }

  if (is_integrated) {
    return(build_integrated_model(occ_formula, det_formula, data, y))
  }

  if (is_dynamic) {
    return(build_dynamic_model(occ_formula, det_formula, data, y,
                                col_formula, ext_formula))
  }

  if (is_community) {
    return(build_community_model(occ_formula, det_formula, data, y, species))
  }

  # Single-season
  if (is.null(det_formula)) stop("det_formula required for non-JSDM models")
  build_single_model(occ_formula, det_formula, data, y,
                     det_visit_formula, det_visit_data)
}

# ============================================================================
# Internal builders
# ============================================================================

build_single_model <- function(occ_formula, det_formula, data, y,
                                det_visit_formula = NULL, det_visit_data = NULL) {
  if (!is.matrix(y)) {
    stop("y must be a matrix (n_sites x max_visits)")
  }
  if (nrow(y) != nrow(data)) {
    stop(sprintf("y has %d rows but data has %d rows", nrow(y), nrow(data)))
  }

  X_occ <- model.matrix(occ_formula, data)
  X_det <- model.matrix(det_formula, data)

  X_det_visit <- NULL
  if (!is.null(det_visit_formula) && !is.null(det_visit_data)) {
    X_det_visit <- model.matrix(det_visit_formula, det_visit_data)
    expected_rows <- nrow(y) * ncol(y)
    if (nrow(X_det_visit) != expected_rows) {
      stop(sprintf("det_visit_data must have %d rows (n_sites * max_visits), got %d",
                   expected_rows, nrow(X_det_visit)))
    }
  }

  y_int <- matrix(as.integer(y), nrow = nrow(y), ncol = ncol(y))
  y_int[is.na(y_int)] <- -1L

  n_detected <- sum(apply(y_int, 1, function(row) any(row[row >= 0] == 1)))

  structure(list(
    model_type = "single",
    y = y_int,
    X_processes = list(X_occ, X_det),
    X_det_visit = X_det_visit,
    formulas = list(occ = occ_formula, det = det_formula),
    n_sites = nrow(y),
    max_visits = ncol(y),
    process_info = list(
      list(name = "psi", p = ncol(X_occ), coef_names = colnames(X_occ)),
      list(name = "p",   p = ncol(X_det), coef_names = colnames(X_det))
    ),
    det_visit_names = if (!is.null(X_det_visit)) colnames(X_det_visit) else character(0),
    naive_occ = n_detected / nrow(y)
  ), class = "TulpaObs")
}

build_dynamic_model <- function(occ_formula, det_formula, data, y,
                                 col_formula, ext_formula) {
  if (is.null(col_formula)) col_formula <- ~ 1
  if (is.null(ext_formula)) ext_formula <- ~ 1

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

  X_occ <- model.matrix(occ_formula, data)
  X_det <- model.matrix(det_formula, data)
  X_col <- model.matrix(col_formula, data)
  X_ext <- model.matrix(ext_formula, data)

  # Flatten detection history and compute metadata
  y_int <- as.integer(y)
  y_int[is.na(y_int)] <- -1L

  n_visits <- integer(n_sites * n_seasons)
  any_detected <- logical(n_sites * n_seasons)

  for (i in seq_len(n_sites)) {
    for (t in seq_len(n_seasons)) {
      idx <- (i - 1) * n_seasons + (t - 1)
      raw <- y[i, , t]
      raw[is.na(raw)] <- -1L
      valid <- raw >= 0
      n_visits[idx + 1] <- sum(valid)
      any_detected[idx + 1] <- any(raw[valid] == 1)
    }
  }

  structure(list(
    model_type = "dynamic",
    y_flat = y_int,
    n_visits = n_visits,
    any_detected = any_detected,
    X_processes = list(X_occ, X_det, X_col, X_ext),
    formulas = list(occ = occ_formula, det = det_formula,
                    col = col_formula, ext = ext_formula),
    n_sites = n_sites,
    n_seasons = n_seasons,
    max_visits = max_visits,
    process_info = list(
      list(name = "psi1",    p = ncol(X_occ), coef_names = colnames(X_occ)),
      list(name = "p",       p = ncol(X_det), coef_names = colnames(X_det)),
      list(name = "gamma",   p = ncol(X_col), coef_names = colnames(X_col)),
      list(name = "epsilon", p = ncol(X_ext), coef_names = colnames(X_ext))
    )
  ), class = "TulpaObs")
}

build_community_model <- function(occ_formula, det_formula, data, y, species) {
  # Handle y as 3D array or list of matrices
  if (is.list(y) && !is.array(y)) {
    n_species <- length(y)
    n_sites <- nrow(y[[1]])
    max_visits <- ncol(y[[1]])
    species_names <- if (is.character(species)) species
                     else if (!is.null(names(y))) names(y)
                     else paste0("sp", seq_len(n_species))
    y_array <- array(NA_integer_, dim = c(n_sites, max_visits, n_species))
    for (s in seq_len(n_species)) {
      y_array[, , s] <- as.integer(y[[s]])
    }
    y <- y_array
  } else {
    species_names <- if (is.character(species)) species
                     else paste0("sp", seq_len(dim(y)[3]))
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

  if (length(species_names) != n_species) {
    stop(sprintf("species has %d names but y has %d species",
                 length(species_names), n_species))
  }

  X_occ <- model.matrix(occ_formula, data)
  X_det <- model.matrix(det_formula, data)

  # Expand to site-species: N = n_sites * n_species
  N <- n_sites * n_species
  X_occ_expanded <- X_occ[rep(seq_len(n_sites), each = n_species), , drop = FALSE]
  X_det_expanded <- X_det[rep(seq_len(n_sites), each = n_species), , drop = FALSE]

  # Flatten detection history to [N x max_visits]
  y_expanded <- matrix(NA_integer_, nrow = N, ncol = max_visits)
  for (i in seq_len(n_sites)) {
    for (s in seq_len(n_species)) {
      obs <- (i - 1) * n_species + s
      y_expanded[obs, ] <- as.integer(y[i, , s])
    }
  }
  y_expanded[is.na(y_expanded)] <- -1L

  # Species grouping (1-based)
  species_group <- rep(seq_len(n_species), times = n_sites)

  structure(list(
    model_type = "community",
    y = y_expanded,
    X_processes = list(X_occ_expanded, X_det_expanded),
    formulas = list(occ = occ_formula, det = det_formula),
    n_sites = n_sites,
    n_species = n_species,
    max_visits = max_visits,
    N = N,
    species_group = as.integer(species_group),
    species_names = species_names,
    process_info = list(
      list(name = "psi", p = ncol(X_occ), coef_names = colnames(X_occ)),
      list(name = "p",   p = ncol(X_det), coef_names = colnames(X_det))
    )
  ), class = "TulpaObs")
}

build_integrated_model <- function(occ_formula, det_formula, data, y) {
  if (!is.list(y) || is.array(y)) {
    stop("For integrated models, y must be a list of detection matrices (one per source)")
  }

  n_sources <- length(y)
  if (inherits(det_formula, "formula")) {
    det_formulas <- rep(list(det_formula), n_sources)
  } else if (is.list(det_formula)) {
    if (length(det_formula) != n_sources) {
      stop(sprintf("det_formula list has %d elements but y has %d sources",
                   length(det_formula), n_sources))
    }
    det_formulas <- det_formula
  } else {
    stop("det_formula must be a formula or list of formulas for integrated models")
  }

  n_sites <- nrow(data)
  X_occ <- model.matrix(occ_formula, data)

  y_sources <- vector("list", n_sources)
  X_det_list <- vector("list", n_sources)
  site_maps <- vector("list", n_sources)

  for (s in seq_len(n_sources)) {
    ys <- y[[s]]
    if (!is.matrix(ys)) stop(sprintf("y[[%d]] must be a matrix", s))

    if (nrow(ys) == n_sites) {
      site_maps[[s]] <- as.integer(seq_len(n_sites) - 1L)
    } else if (!is.null(rownames(ys))) {
      site_maps[[s]] <- as.integer(match(rownames(ys), rownames(data)) - 1L)
    } else {
      site_maps[[s]] <- as.integer(seq_len(nrow(ys)) - 1L)
    }

    y_int <- matrix(as.integer(ys), nrow = nrow(ys), ncol = ncol(ys))
    y_int[is.na(y_int)] <- -1L
    y_sources[[s]] <- y_int

    src_rows <- site_maps[[s]] + 1L
    X_det_list[[s]] <- model.matrix(det_formulas[[s]], data[src_rows, , drop = FALSE])
  }

  X_processes <- vector("list", 1 + n_sources)
  X_processes[[1]] <- X_occ

  process_info <- list(
    list(name = "psi", p = ncol(X_occ), coef_names = colnames(X_occ))
  )

  for (s in seq_len(n_sources)) {
    X_det_s <- matrix(0, n_sites, ncol(X_det_list[[s]]))
    src_rows <- site_maps[[s]] + 1L
    X_det_s[src_rows, ] <- X_det_list[[s]]
    X_processes[[1 + s]] <- X_det_s

    src_name <- if (!is.null(names(y))) names(y)[s] else paste0("p", s)
    process_info[[1 + s]] <- list(
      name = src_name, p = ncol(X_det_list[[s]]),
      coef_names = colnames(X_det_list[[s]])
    )
  }

  structure(list(
    model_type = "integrated",
    y_sources = y_sources,
    site_maps = site_maps,
    X_processes = X_processes,
    formulas = c(list(occ = occ_formula), setNames(det_formulas, names(y))),
    n_sites = n_sites,
    n_sources = n_sources,
    process_info = process_info
  ), class = "TulpaObs")
}

build_jsdm_model <- function(occ_formula, data, y, species) {
  # y: 3D array [n_sites x 1 x n_species] or named list of vectors/single-col matrices
  if (is.list(y) && !is.array(y)) {
    n_species <- length(y)
    n_sites <- length(y[[1]])
    species_names <- if (is.character(species)) species
                     else if (!is.null(names(y))) names(y)
                     else paste0("sp", seq_len(n_species))
    y_mat <- matrix(NA_integer_, n_sites, n_species)
    for (s in seq_len(n_species)) {
      y_mat[, s] <- as.integer(if (is.matrix(y[[s]])) y[[s]][, 1] else y[[s]])
    }
  } else if (is.matrix(y)) {
    # y: n_sites x n_species matrix of 0/1
    n_sites <- nrow(y)
    n_species <- ncol(y)
    species_names <- if (is.character(species)) species
                     else if (!is.null(colnames(y))) colnames(y)
                     else paste0("sp", seq_len(n_species))
    y_mat <- matrix(as.integer(y), n_sites, n_species)
  } else {
    stop("For JSDM, y must be a matrix (n_sites x n_species) or named list")
  }

  if (nrow(data) != n_sites) {
    stop(sprintf("y has %d sites but data has %d rows", n_sites, nrow(data)))
  }

  X_occ <- model.matrix(occ_formula, data)

  # Expand to site-species (same as community but no detection)
  N <- n_sites * n_species
  X_occ_expanded <- X_occ[rep(seq_len(n_sites), each = n_species), , drop = FALSE]

  # Flatten y to vector: site-species order
  y_flat <- integer(N)
  for (i in seq_len(n_sites)) {
    for (s in seq_len(n_species)) {
      obs <- (i - 1) * n_species + s
      y_flat[obs] <- y_mat[i, s]
    }
  }
  y_flat[is.na(y_flat)] <- 0L

  species_group <- rep(seq_len(n_species), times = n_sites)

  structure(list(
    model_type = "jsdm",
    y_jsdm = y_flat,
    X_processes = list(X_occ_expanded),
    formulas = list(occ = occ_formula),
    n_sites = n_sites,
    n_species = n_species,
    N = N,
    species_group = as.integer(species_group),
    species_names = species_names,
    process_info = list(
      list(name = "psi", p = ncol(X_occ), coef_names = colnames(X_occ))
    )
  ), class = "TulpaObs")
}
