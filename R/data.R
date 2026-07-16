#' Format occupancy data from matrices/lists
#'
#' @param y Detection history matrix (n_sites x max_visits) or 3D array.
#' @param occ.covs Data.frame of site-level covariates, or named list.
#' @param det.covs Named list of detection covariates. Each element is a vector
#'   (length n_sites, constant across visits) or matrix (n_sites x max_visits).
#' @param coords Optional n_sites x 2 coordinate matrix.
#' @param species Optional character or integer species identifier.
#' @return An `tobs_data` object.
#' @export
tobs_format <- function(y, occ.covs = NULL, det.covs = NULL,
                        coords = NULL, species = NULL) {
  if (is.list(occ.covs) && !is.data.frame(occ.covs)) {
    occ.covs <- as.data.frame(occ.covs)
  }

  if (!is.null(det.covs) && !is.list(det.covs)) {
    stop("det.covs must be a named list")
  }

  structure(list(
    y = y,
    occ.covs = occ.covs,
    det.covs = det.covs,
    coords = coords,
    species = species
  ), class = "tobs_data")
}

#' Convert long-format data to a site x visit observation object
#'
#' @param df Data.frame in long format (one row per site-visit).
#' @param y Character, name of the response column. Its meaning is set by
#'   `type`: a 0/1 detection (`"occurrence"`), an integer count
#'   (`"abundance"`), or a continuous cover proportion in `[0, 1]` (`"cover"`).
#' @param site Character, name of site identifier column.
#' @param visit Character, name of visit/replicate column.
#' @param type Response kind. `"occurrence"` (default) and `"abundance"` build
#'   an integer site x visit matrix; `"cover"` and `"positive"` build a double
#'   matrix and preserve continuous values (no coercion to integer). `"cover"`
#'   is a proportion in `[0, 1]` (the beta cover arm); `"positive"` is a
#'   positive real `(0, Inf)` (the lognormal / gamma cover arm), validated as
#'   non-negative with no upper bound. Both apply the `cover.floor` absence
#'   policy below.
#' @param cover.floor For `type = "cover"` / `"positive"`, the threshold at or
#'   below which a value is stored as `NA` rather than as a positive observation
#'   (default `0`). The positive arm of a hurdle is positive-only, so a `0` is
#'   an absence handled by the occurrence arm, not a positive observation;
#'   sending it to `NA` keeps a sampled-absent or unsampled cell from entering
#'   the positive arm as a fabricated zero (which, padded across a grid,
#'   flattens the spatial field). Set `cover.floor = -Inf` to keep every value
#'   verbatim.
#' @param occ.covs Character vector of site-level covariate names.
#' @param det.covs Character vector of visit-level covariate names. A named
#'   column that is a factor or character is preserved as categorical: the
#'   downstream detection / positive design expands it to k - 1 dummies for a
#'   k-level factor, with the factor's first level (for a character column, the
#'   first level after sorting the unique values) as the reference. Numeric
#'   columns are kept as a continuous covariate.
#' @param coords Character vector of length 2 for coordinate columns.
#' @param sites Optional site level set. When supplied, the site x visit grid
#'   uses these site identifiers in this order (rather than `df`'s
#'   first-appearance order), so a subset of `df` pivots onto a fixed,
#'   externally-defined site grid; site values in `df` outside the set error.
#' @param visits Optional visit level set, analogous to `sites` for the visit
#'   axis (default the sorted unique visits in `df`).
#' @param compact Logical (default `FALSE`). When `TRUE`, return a compact
#'   (ragged) `tobs_data`: the response is stored as one row per valid
#'   site-visit (a `tobs_ragged` carrier) rather than as a padded
#'   `[n_sites x max_visits]` matrix, and each detection covariate as a
#'   length-V vector in the same canonical `order(site, visit)`. The compact
#'   layout has no per-site visit cap (its memory is the number of observations,
#'   not the padded grid) and is consumed by the joint nested-Laplace
#'   `occu_cover()` engine, which works one valid visit at a time. Two compact
#'   calls on the same `df` / `site` / `visit` align row-for-row, so an
#'   occurrence response and a cover response can be paired directly.
#' @return An `tobs_data` object. With `compact = TRUE` its `$y` is a
#'   `tobs_ragged` carrier and its `$det.covs` are length-V vectors.
#' @export
tobs_data <- function(df, y, site, visit,
                      type = c("occurrence", "abundance", "cover", "positive"),
                      occ.covs = NULL, det.covs = NULL,
                      coords = NULL,
                      sites = NULL, visits = NULL,
                      cover.floor = 0,
                      compact = FALSE) {
  type <- match.arg(type)
  if (!is.data.frame(df)) stop("df must be a data.frame")
  for (col in c(y, site, visit)) {
    if (!col %in% names(df)) stop(sprintf("column '%s' not found in df", col))
  }

  # Site / visit level sets. By default they are read from `df` (first-appearance
  # order for sites, sorted for visits). An external set (`sites` / `visits`) is
  # supplied by the per-species batch builder so every species pivots onto ONE
  # canonical site x visit grid (a species absent at a site must still occupy
  # that site's row, aligned with the shared cell-level design); rows / columns
  # the species never visits stay NA and are validated against the level set.
  if (is.null(sites)) sites <- unique(df[[site]])
  if (is.null(visits)) visits <- sort(unique(df[[visit]]))
  n_sites <- length(sites)
  max_visits <- length(visits)

  si <- match(df[[site]], sites)
  vi <- match(df[[visit]], visits)
  if (anyNA(si)) {
    stop("tobs_data(): site value(s) in df are not in the supplied `sites` set.",
         call. = FALSE)
  }
  if (anyNA(vi)) {
    stop("tobs_data(): visit value(s) in df are not in the supplied `visits` set.",
         call. = FALSE)
  }

  # Compact (ragged) layout. The dense pivot below stores a [n_sites x max_visits]
  # matrix per response / detection covariate, so its memory is the PADDED grid
  # (NA at every unsampled site x visit). When one site holds tens of thousands of
  # visits, that grid forces a per-site visit cap before the data can be built.
  # The joint nested-Laplace occu_cover engine consumes one row per VALID visit
  # (it compacts the dense grid right back to `which(valid)` rows), so the padding
  # is pure waste there. `compact = TRUE` keeps the valid rows directly: the
  # response is a `tobs_ragged` carrier (value + site + visit, in canonical
  # order(site, visit)) and each detection covariate a length-V vector in the same
  # order. Two compact calls on the same `df` / `site` / `visit` share the order,
  # so an occurrence response and a cover response align row-for-row.
  if (compact) {
    ord  <- order(si, vi)
    si_o <- si[ord]; vi_o <- vi[ord]
    yv   <- df[[y]][ord]
    if (type %in% c("cover", "positive")) {
      values <- .tobs_floor_continuous(yv, type, cover.floor)
    } else {
      yi <- as.integer(yv)
      if (type == "occurrence" && !all(is.na(yi) | yi %in% c(0L, 1L)))
        stop("tobs_data(type = 'occurrence'): y must be 0/1")
      if (type == "abundance" && any(yi < 0L, na.rm = TRUE))
        stop("tobs_data(type = 'abundance'): y must be non-negative integer counts")
      values <- yi
    }

    occ_df <- NULL
    if (!is.null(occ.covs)) {
      site_rows <- match(sites, df[[site]])
      occ_df <- df[site_rows, occ.covs, drop = FALSE]
      rownames(occ_df) <- NULL
    }

    det_list <- NULL
    if (!is.null(det.covs)) {
      det_list <- lapply(det.covs, function(dc) {
        col <- df[[dc]][ord]
        if (is.factor(col) || is.character(col)) {
          levs <- if (is.factor(df[[dc]])) levels(df[[dc]]) else sort(unique(df[[dc]]))
          v <- as.character(col)
          attr(v, "tobs_factor") <- TRUE
          attr(v, "tobs_levels") <- levs
          v
        } else {
          as.numeric(col)
        }
      })
      names(det_list) <- det.covs
    }

    coord_mat <- NULL
    if (!is.null(coords)) {
      site_rows <- match(sites, df[[site]])
      coord_mat <- as.matrix(df[site_rows, coords])
    }

    y_ragged <- structure(
      list(values = values, site = si_o, visit = vi_o,
           n_sites = n_sites, max_visits = max_visits,
           n_visits = length(values), sites = sites, visits = visits,
           type = type),
      class = "tobs_ragged")

    return(structure(
      list(y = y_ragged, occ.covs = occ_df, det.covs = det_list,
           coords = coord_mat, compact = TRUE),
      class = "tobs_data"))
  }

  y_mat <- .tobs_long_response_matrix(df[[y]], si, vi, n_sites, max_visits, type,
                                      cover.floor = cover.floor)
  rownames(y_mat) <- sites

  # Extract site-level covariates
  occ_df <- NULL
  if (!is.null(occ.covs)) {
    # Take first occurrence per site
    site_rows <- match(sites, df[[site]])
    occ_df <- df[site_rows, occ.covs, drop = FALSE]
    rownames(occ_df) <- NULL
  }

  # Extract visit-level detection covariates (same site x visit layout as y_mat).
  # Numeric columns become a double matrix; factor / character columns become a
  # character matrix tagged with their level set so the downstream design keeps
  # them categorical (k - 1 dummies, first level the reference). Both share the
  # 2D-index fill (column-major-safe: `cbind(si, vi)`, never a linear slot).
  det_list <- NULL
  if (!is.null(det.covs)) {
    det_list <- lapply(det.covs, function(dc) {
      col <- df[[dc]]
      if (is.factor(col) || is.character(col)) {
        levs <- if (is.factor(col)) levels(col) else sort(unique(col))
        mat <- matrix(NA_character_, n_sites, max_visits)
        mat[cbind(si, vi)] <- as.character(col)
        attr(mat, "tobs_factor") <- TRUE
        attr(mat, "tobs_levels") <- levs
        mat
      } else {
        mat <- matrix(NA_real_, n_sites, max_visits)
        mat[cbind(si, vi)] <- as.numeric(col)
        mat
      }
    })
    names(det_list) <- det.covs
  }

  # Extract coordinates
  coord_mat <- NULL
  if (!is.null(coords)) {
    site_rows <- match(sites, df[[site]])
    coord_mat <- as.matrix(df[site_rows, coords])
  }

  tobs_format(y = y_mat, occ.covs = occ_df, det.covs = det_list,
              coords = coord_mat)
}

# Fill a site x visit response matrix from a long response vector at 2D indices
# (`si`, `vi`). The single source of truth for the long -> matrix pivot used by
# `tobs_data()` and the per-species `by=` batch builder. The response `type`
# drives both storage and validation: occurrence / abundance are integer
# (counts), cover is a continuous proportion kept as double -- coercing it to
# integer would truncate every value < 1 to zero. Cells not addressed by
# (`si`, `vi`) stay NA (column-major-safe: 2D indexing, never a linear slot).
.tobs_long_response_matrix <- function(yv, si, vi, n_sites, max_visits, type,
                                       cover.floor = 0) {
  if (type %in% c("cover", "positive")) {
    y_mat <- matrix(NA_real_, n_sites, max_visits)
    y_mat[cbind(si, vi)] <- .tobs_floor_continuous(yv, type, cover.floor)
  } else {
    yi <- as.integer(yv)
    if (type == "occurrence" && !all(is.na(yi) | yi %in% c(0L, 1L)))
      stop("tobs_data(type = 'occurrence'): y must be 0/1")
    if (type == "abundance" && any(yi < 0L, na.rm = TRUE))
      stop("tobs_data(type = 'abundance'): y must be non-negative integer counts")
    y_mat <- matrix(NA_integer_, n_sites, max_visits)
    y_mat[cbind(si, vi)] <- yi
  }
  y_mat
}

# Validate + floor a continuous positive-arm response, the single source for the
# dense and compact `tobs_data()` paths. `type = "cover"` is a proportion in
# `[0, 1]` (beta); `type = "positive"` is a positive real `(0, Inf)` (lognormal
# / gamma), validated non-negative with no upper bound. A value at or below
# `cover.floor` (default 0) is an absence handled by the occurrence arm, so it
# becomes NA rather than a fabricated zero (which, padded across a grid, flattens
# the spatial field). Returns the floored numeric vector.
.tobs_floor_continuous <- function(yv, type, cover.floor) {
  if (type == "cover") {
    if (any(yv < 0 | yv > 1, na.rm = TRUE))
      stop("tobs_data(type = 'cover'): y must lie in [0, 1]", call. = FALSE)
  } else {
    if (any(yv < 0, na.rm = TRUE))
      stop("tobs_data(type = 'positive'): y must be non-negative (a 0 is an ",
           "absence floored to NA; a negative value is invalid).",
           call. = FALSE)
  }
  noun <- if (type == "cover") "cover value" else "positive value"
  floored <- !is.na(yv) & yv <= cover.floor
  n_floor <- sum(floored)
  if (n_floor > 0L) {
    yv[floored] <- NA_real_
    message(sprintf(paste0(
      "tobs_data(type = '%s'): %d %s%s <= %g treated as absent (NA); the ",
      "positive arm sees only positive values, so a 0 is an absence for the ",
      "occurrence arm, not a fabricated zero. Set cover.floor = -Inf to keep ",
      "them."), type, n_floor, noun, if (n_floor == 1L) "" else "s",
      cover.floor))
  }
  as.numeric(yv)
}

#' Format multi-species occupancy data
#'
#' @param y 3D array (n_sites x max_visits x n_species) or named list of matrices.
#' @param occ.covs Data.frame of site-level covariates.
#' @param det.covs Named list of detection covariates.
#' @param coords Optional n_sites x 2 coordinate matrix.
#' @param species_names Optional character vector of species names.
#' @return An `tobs_data` object with multi-species structure.
#' @export
tobs_format_ms <- function(y, occ.covs = NULL, det.covs = NULL,
                           coords = NULL, species_names = NULL) {
  if (is.list(y) && !is.array(y)) {
    n_species <- length(y)
    n_sites <- nrow(y[[1]])
    max_visits <- ncol(y[[1]])
    if (is.null(species_names)) species_names <- names(y)
    y_array <- array(NA_integer_, dim = c(n_sites, max_visits, n_species))
    for (s in seq_len(n_species)) {
      y_array[, , s] <- as.integer(y[[s]])
    }
    y <- y_array
  }
  if (is.null(species_names)) {
    species_names <- paste0("sp", seq_len(dim(y)[3]))
  }
  if (is.list(occ.covs) && !is.data.frame(occ.covs)) {
    occ.covs <- as.data.frame(occ.covs)
  }

  structure(list(
    y = y,
    occ.covs = occ.covs,
    det.covs = det.covs,
    coords = coords,
    species_names = species_names,
    n_species = dim(y)[3]
  ), class = "tobs_data")
}

#' @export
print.tobs_ragged <- function(x, ...) {
  cat(sprintf("tobs_ragged (%s): %d sites, %d valid visits (max %d / site)\n",
              x$type, x$n_sites, x$n_visits, x$max_visits))
  invisible(x)
}

#' @export
print.tobs_data <- function(x, ...) {
  if (inherits(x$y, "tobs_ragged")) {
    cat(sprintf("tobs_data (compact): %d sites, %d valid visits (max %d / site)\n",
                x$y$n_sites, x$y$n_visits, x$y$max_visits))
    if (!is.null(x$occ.covs))
      cat(sprintf("  Occupancy covariates: %s\n", paste(names(x$occ.covs), collapse = ", ")))
    if (!is.null(x$det.covs))
      cat(sprintf("  Detection covariates: %s\n", paste(names(x$det.covs), collapse = ", ")))
    if (!is.null(x$coords)) cat("  Coordinates: yes\n")
    return(invisible(x))
  }
  n_sites <- nrow(x$y)
  max_visits <- if (is.matrix(x$y)) ncol(x$y) else dim(x$y)[2]
  cat(sprintf("tobs_data: %d sites, %d visits\n", n_sites, max_visits))
  if (!is.null(x$occ.covs)) {
    cat(sprintf("  Occupancy covariates: %s\n", paste(names(x$occ.covs), collapse = ", ")))
  }
  if (!is.null(x$det.covs)) {
    cat(sprintf("  Detection covariates: %s\n", paste(names(x$det.covs), collapse = ", ")))
  }
  if (!is.null(x$coords)) cat("  Coordinates: yes\n")
  invisible(x)
}

#' Summarise occupancy data
#'
#' Returns detection statistics including naive occupancy, naive detection,
#' per-visit rates, and detection frequency table.
#'
#' @param object An `tobs_data` object.
#' @param ... Ignored.
#' @return An `tobs_data_summary` object (printed automatically).
#' @export
summary.tobs_data <- function(object, ...) {
  y <- if (is.matrix(object$y)) object$y else object$y[, , 1]
  N <- nrow(y)
  J <- ncol(y)

  not_na <- y >= 0 | !is.na(y)
  # Treat -1 and NA as missing
  valid <- !is.na(y) & y >= 0
  n_obs <- sum(valid)
  n_missing <- N * J - n_obs

  det_count <- rowSums(y == 1 & valid)
  n_visits <- rowSums(valid)
  detected <- det_count > 0

  naive_psi <- mean(detected)
  naive_p <- if (sum(n_visits[detected]) > 0) {
    sum(det_count[detected]) / sum(n_visits[detected])
  } else NA_real_

  det_freq <- table(factor(det_count, levels = 0:J))
  det_per_visit <- vapply(seq_len(J), function(j) {
    v <- y[, j]
    ok <- !is.na(v) & v >= 0
    if (sum(ok) > 0) mean(v[ok]) else NA_real_
  }, numeric(1))

  out <- list(
    N = N, J = J, n_obs = n_obs, n_missing = n_missing,
    naive_psi = naive_psi, naive_p = naive_p,
    det_freq = det_freq, det_count = det_count,
    n_visits = n_visits, det_per_visit = det_per_visit,
    has_coords = !is.null(object$coords)
  )
  class(out) <- "tobs_data_summary"
  out
}

#' @export
print.tobs_data_summary <- function(x, ...) {
  cat("Occupancy data summary\n")
  cat(sprintf("  Sites: %d | Max visits: %d\n", x$N, x$J))
  cat(sprintf("  Observations: %d | Missing: %d (%.1f%%)\n",
              x$n_obs, x$n_missing,
              100 * x$n_missing / (x$n_obs + x$n_missing)))
  cat(sprintf("  Naive occupancy: %.3f (%d / %d sites)\n",
              x$naive_psi, sum(x$det_freq[-1]), x$N))
  if (!is.na(x$naive_p)) {
    cat(sprintf("  Naive detection: %.3f\n", x$naive_p))
  }
  cat("\n  Detection frequency (detections per site):\n")
  df <- as.data.frame(x$det_freq)
  names(df) <- c("detections", "sites")
  print(df, row.names = FALSE)
  cat(sprintf("\n  Per-visit detection rate: %s\n",
              paste(sprintf("V%d=%.2f", seq_along(x$det_per_visit),
                            x$det_per_visit), collapse = "  ")))
  if (x$has_coords) cat("  Coordinates: available\n")
  invisible(x)
}

#' Plot detection history patterns
#'
#' 2x2 panel: detection frequency histogram, per-visit detection rates,
#' visit completeness, and spatial detection map (if coordinates available).
#'
#' @param x An `tobs_data` object.
#' @param ... Ignored.
#' @return Invisible `NULL`.
#' @importFrom graphics hist barplot par polygon legend
#' @importFrom grDevices rgb
#' @export
plot.tobs_data <- function(x, ...) {
  y <- if (is.matrix(x$y)) x$y else x$y[, , 1]
  N <- nrow(y)
  J <- ncol(y)

  valid <- !is.na(y) & y >= 0
  det_count <- rowSums(y == 1 & valid)
  detected <- det_count > 0

  has_coords <- !is.null(x$coords)
  old_par <- par(mfrow = if (has_coords) c(2, 2) else c(1, 3),
                 mar = c(4, 4, 2.5, 1))
  on.exit(par(old_par))

  # Panel 1: detection frequency
  hist(det_count, breaks = seq(-0.5, max(det_count) + 0.5, by = 1),
       main = "Detections per site", xlab = "Number of detections",
       col = "grey80", border = "grey50")

  # Panel 2: per-visit detection rate
  det_per_visit <- vapply(seq_len(J), function(j) {
    v <- y[, j]; ok <- !is.na(v) & v >= 0
    if (sum(ok) > 0) mean(v[ok]) else 0
  }, numeric(1))
  barplot(det_per_visit, names.arg = paste0("V", seq_len(J)),
          main = "Detection rate by visit",
          ylab = "P(detect)", col = "steelblue",
          ylim = c(0, max(det_per_visit) * 1.2 + 0.01))

  # Panel 3: visit completeness
  completeness <- vapply(seq_len(J), function(j) {
    mean(!is.na(y[, j]) & y[, j] >= 0)
  }, numeric(1))
  barplot(completeness, names.arg = paste0("V", seq_len(J)),
          main = "Visit completeness",
          ylab = "Proportion surveyed", col = "darkseagreen",
          ylim = c(0, 1))

  # Panel 4: spatial map
  if (has_coords) {
    cols <- ifelse(detected, "tomato", "grey70")
    plot(x$coords[, 1], x$coords[, 2],
         col = cols, pch = 19,
         cex = 0.5 + det_count / max(max(det_count), 1),
         xlab = "X", ylab = "Y", main = "Detection map")
    legend("topright",
           legend = c("Detected", "Not detected"),
           col = c("tomato", "grey70"), pch = 19, cex = 0.8)
  }
  invisible(NULL)
}

# ============================================================================
# Simulation functions
# ============================================================================

#' Simulate single-species occupancy data
#'
#' @param N Number of sites (default 100).
#' @param J Number of visits (default 4).
#' @param n_occ_covs Number of occupancy covariates (default 2).
#' @param n_det_covs Number of detection covariates (default 1).
#' @param beta_occ Occupancy coefficients (auto-generated if NULL).
#' @param beta_det Detection coefficients (auto-generated if NULL).
#' @param seed Random seed.
#' @return A list with `y`, `data`, and `truth`.
#' @export
simulate_occu <- function(N = 100, J = 4,
                          n_occ_covs = 2, n_det_covs = 1,
                          beta_occ = NULL, beta_det = NULL,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(beta_occ)) beta_occ <- c(0, runif(n_occ_covs, -1, 1))
  if (is.null(beta_det)) beta_det <- c(0, runif(n_det_covs, -1, 1))

  # Covariates
  occ_covs <- data.frame(matrix(rnorm(N * n_occ_covs), N, n_occ_covs))
  names(occ_covs) <- paste0("occ_cov", seq_len(n_occ_covs))
  det_covs <- data.frame(matrix(rnorm(N * n_det_covs), N, n_det_covs))
  names(det_covs) <- paste0("det_cov", seq_len(n_det_covs))
  data <- cbind(occ_covs, det_covs)

  X_occ <- model.matrix(~ ., occ_covs)
  X_det <- model.matrix(~ ., det_covs)

  psi <- plogis(as.vector(X_occ %*% beta_occ))
  p <- plogis(as.vector(X_det %*% beta_det))
  z <- rbinom(N, 1, psi)

  y <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) {
    y[i, ] <- rbinom(J, 1, z[i] * p[i])
  }

  list(
    y = y,
    data = data,
    truth = list(beta_occ = beta_occ, beta_det = beta_det, psi = psi, p = p, z = z)
  )
}

#' Simulate multi-species occupancy data
#'
#' @param N Number of sites (default 100).
#' @param J Number of visits (default 4).
#' @param n_species Number of species (default 10).
#' @param beta_comm_mean Community mean for occupancy (default c(0, 0.5)).
#' @param beta_comm_sd Community SD for occupancy (default c(0.5, 0.3)).
#' @param alpha_comm_mean Community mean for detection (default c(0)).
#' @param alpha_comm_sd Community SD for detection (default c(0.5)).
#' @param seed Random seed.
#' @return A list with `y` (3D array), `data`, and `truth`.
#' @export
simulate_ms_occu <- function(N = 100, J = 4, n_species = 10,
                     beta_comm_mean = c(0, 0.5),
                     beta_comm_sd = c(0.5, 0.3),
                     alpha_comm_mean = c(0),
                     alpha_comm_sd = c(0.5),
                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_occ_covs <- length(beta_comm_mean) - 1
  n_det_covs <- length(alpha_comm_mean) - 1

  data <- data.frame(x = rnorm(N))
  if (n_occ_covs > 1) {
    for (k in 2:n_occ_covs) data[[paste0("occ_cov", k)]] <- rnorm(N)
  }

  X_occ <- model.matrix(~ ., data[, seq_len(n_occ_covs + (n_occ_covs == 0)), drop = FALSE])

  # Species-specific coefficients
  beta_species <- matrix(NA_real_, n_species, length(beta_comm_mean))
  alpha_species <- matrix(NA_real_, n_species, length(alpha_comm_mean))
  for (j in seq_along(beta_comm_mean)) {
    beta_species[, j] <- rnorm(n_species, beta_comm_mean[j], beta_comm_sd[j])
  }
  for (j in seq_along(alpha_comm_mean)) {
    alpha_species[, j] <- rnorm(n_species, alpha_comm_mean[j], alpha_comm_sd[j])
  }

  y <- array(NA_integer_, dim = c(N, J, n_species))
  z <- matrix(NA_integer_, N, n_species)
  for (s in seq_len(n_species)) {
    psi_s <- plogis(as.vector(X_occ %*% beta_species[s, ]))
    p_s <- plogis(alpha_species[s, 1])
    z[, s] <- rbinom(N, 1, psi_s)
    for (i in seq_len(N)) {
      y[i, , s] <- rbinom(J, 1, z[i, s] * p_s)
    }
  }

  list(
    y = y,
    data = data,
    truth = list(
      beta_species = beta_species,
      alpha_species = alpha_species,
      beta_comm_mean = beta_comm_mean,
      beta_comm_sd = beta_comm_sd,
      z = z
    )
  )
}

#' Simulate community relative-abundance (count / continuous) data
#'
#' Per-species GLMM with Gaussian community hyperpriors on the coefficients (the
#' `ms_count()` / spAbundance `msAbund` model): `g(mu_{s,i}) = X_i beta_s`,
#' `beta_s ~ N(beta_comm_mean, diag(beta_comm_sd^2))`, `y_{s,i}` drawn Poisson /
#' negative-binomial (log link) or Gaussian (identity). No detection.
#'
#' @param N Number of sites.
#' @param n_species Number of species.
#' @param beta_comm_mean,beta_comm_sd Community mean and SD of the per-species
#'   coefficients (intercept first). Length sets the number of covariates.
#' @param response One of `"poisson"`, `"negbin"`, `"gaussian"`, `"binomial"`.
#' @param size Negative-binomial community size (mean of the per-species
#'   `log_r`); `size.log.sd` is its across-species SD.
#' @param size.log.sd Across-species SD of `log(size)` (negbin only).
#' @param sd Gaussian residual SD.
#' @param trials Binomial trial count (`response = "binomial"`): a scalar
#'   (shared across sites and species) or a length-`N` per-site vector. Default
#'   10.
#' @param seed Optional RNG seed.
#' @return A list with `y` (an `N x n_species` matrix), `data`, and `truth`.
#' @export
simulate_ms_count <- function(N = 120, n_species = 10,
                              beta_comm_mean = c(1, 0.5),
                              beta_comm_sd = c(0.4, 0.3),
                              response = c("poisson", "negbin", "gaussian",
                                           "binomial"),
                              size = 2, size.log.sd = 0.3, sd = 1,
                              trials = 10, seed = NULL) {
  response <- match.arg(response)
  if (!is.null(seed)) set.seed(seed)
  n_cov <- length(beta_comm_mean) - 1L

  data <- data.frame(x = stats::rnorm(N))
  if (n_cov > 1L) for (k in 2:n_cov) data[[paste0("cov", k)]] <- stats::rnorm(N)
  X <- stats::model.matrix(~ ., data[, seq_len(max(n_cov, 1L)), drop = FALSE])

  beta_species <- matrix(NA_real_, n_species, length(beta_comm_mean))
  for (j in seq_along(beta_comm_mean)) {
    beta_species[, j] <- stats::rnorm(n_species, beta_comm_mean[j],
                                      beta_comm_sd[j])
  }
  r_s <- if (identical(response, "negbin"))
    exp(stats::rnorm(n_species, log(size), size.log.sd)) else rep(NA_real_, n_species)

  is_binom <- identical(response, "binomial")
  n_trials <- NULL
  if (is_binom) {
    nt <- if (length(trials) == 1L) rep(as.integer(trials), N)
          else as.integer(trials)
    if (length(nt) != N) stop("simulate_ms_count(): `trials` must be a scalar ",
                              "or length N.", call. = FALSE)
    n_trials <- matrix(nt, N, n_species)
  }

  is_gauss <- identical(response, "gaussian")
  y <- matrix(NA_real_, N, n_species,
              dimnames = list(NULL, paste0("sp", seq_len(n_species))))
  for (s in seq_len(n_species)) {
    eta <- as.numeric(X %*% beta_species[s, ])
    mu  <- switch(response, gaussian = eta, binomial = stats::plogis(eta),
                  exp(eta))
    y[, s] <- switch(response,
      poisson  = stats::rpois(N, mu),
      negbin   = stats::rnbinom(N, size = r_s[s], mu = mu),
      gaussian = stats::rnorm(N, mu, sd),
      binomial = stats::rbinom(N, size = n_trials[, s], prob = mu))
  }

  list(y = y, data = data,
       truth = list(beta_species = beta_species,
                    beta_comm_mean = beta_comm_mean,
                    beta_comm_sd = beta_comm_sd,
                    response = response, r_s = r_s, sd = sd,
                    trials = n_trials))
}

#' Simulate joint species distribution (presence/absence) data
#'
#' Presence `y_{i,s} ~ Bernoulli(psi_{i,s})` with `logit psi = X beta_s`, the
#' per-species occupancy coefficients drawn from Gaussian community hyperpriors.
#' Matches the [jsdm()] family, which observes presence directly (no detection
#' process), so the response is an `N x n_species` presence matrix.
#'
#' @param N Number of sites (default 100).
#' @param n_species Number of species (default 10).
#' @param beta_comm_mean Community-mean occupancy coefficients, length
#'   `1 + n_occ_covs` (intercept first).
#' @param beta_comm_sd Between-species SD of each occupancy coefficient (same
#'   length as `beta_comm_mean`).
#' @param seed Random seed.
#' @return A list with `y` (an `N x n_species` 0/1 presence matrix), `data`, and
#'   `truth` (per-species coefficients and the community hyperparameters).
#' @examples
#' sim <- simulate_jsdm(N = 60, n_species = 5, seed = 1)
#' dim(sim$y)   # 60 sites x 5 species presence matrix
#' @export
simulate_jsdm <- function(N = 100, n_species = 10,
                          beta_comm_mean = c(0, 0.5),
                          beta_comm_sd = c(0.5, 0.3),
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_occ_covs <- length(beta_comm_mean) - 1

  data <- data.frame(x = rnorm(N))
  if (n_occ_covs > 1) {
    for (k in 2:n_occ_covs) data[[paste0("occ_cov", k)]] <- rnorm(N)
  }
  X_occ <- model.matrix(~ .,
    data[, seq_len(n_occ_covs + (n_occ_covs == 0)), drop = FALSE])

  beta_species <- matrix(NA_real_, n_species, length(beta_comm_mean))
  for (j in seq_along(beta_comm_mean)) {
    beta_species[, j] <- rnorm(n_species, beta_comm_mean[j], beta_comm_sd[j])
  }

  y <- matrix(NA_integer_, N, n_species)
  for (s in seq_len(n_species)) {
    y[, s] <- rbinom(N, 1, plogis(as.vector(X_occ %*% beta_species[s, ])))
  }
  colnames(y) <- paste0("sp", seq_len(n_species))

  list(y = y, data = data,
       truth = list(beta_species   = beta_species,
                    beta_comm_mean = beta_comm_mean,
                    beta_comm_sd   = beta_comm_sd))
}

#' Simulate count / relative-abundance GLMM data
#'
#' One observed value per site from a GLMM with no detection process, matching
#' the [count()] family: `log`-link Poisson / negative-binomial counts, an
#' identity-link Gaussian response, or a logit-link binomial response (`k`
#' successes out of `n` trials per site). The mean predictor is `X beta` with
#' `beta` the fixed-effect coefficients (intercept first).
#'
#' @param N Number of sites (default 200).
#' @param beta Fixed-effect coefficients, length `1 + n_covs` (intercept first).
#' @param response The response distribution: `"poisson"`, `"negbin"`,
#'   `"gaussian"`, or `"binomial"`.
#' @param size Negative-binomial size (dispersion); larger is closer to Poisson.
#'   Used only for `response = "negbin"`.
#' @param sd Gaussian residual SD. Used only for `response = "gaussian"`.
#' @param trials Per-site trial count for `response = "binomial"`: a scalar
#'   (recycled) or a length-`N` vector. Default 10 (a scalar 1 gives Bernoulli
#'   data, the `svcPGBinom` `trials = 1` setting).
#' @param seed Random seed.
#' @return A list with `y` (a length-`N` numeric response vector; success counts
#'   for the binomial response), `data` (the site covariates), and `truth` (the
#'   coefficients, response, dispersion, and the binomial `trials`).
#' @examples
#' sim <- simulate_count(N = 100, beta = c(1, 0.5), seed = 1)
#' length(sim$y)
#' @export
simulate_count <- function(N = 200, beta = c(1, 0.5),
                           response = c("poisson", "negbin", "gaussian",
                                        "binomial"),
                           size = 2, sd = 1, trials = 10, seed = NULL) {
  response <- match.arg(response)
  if (!is.null(seed)) set.seed(seed)
  n_cov <- length(beta) - 1L

  data <- data.frame(x = rnorm(N))
  if (n_cov > 1L) {
    for (k in 2:n_cov) data[[paste0("x", k)]] <- rnorm(N)
  }
  X   <- model.matrix(~ ., data)
  eta <- as.vector(X %*% beta)
  mu  <- switch(response,
    gaussian = eta,
    binomial = stats::plogis(eta),
    exp(eta))

  n_trials <- NULL
  if (identical(response, "binomial")) {
    n_trials <- if (length(trials) == 1L) rep(as.integer(trials), N)
                else as.integer(trials)
    if (length(n_trials) != N) {
      stop("simulate_count(): `trials` must be a scalar or length N.",
           call. = FALSE)
    }
  }

  y <- switch(response,
    poisson  = stats::rpois(N, mu),
    negbin   = stats::rnbinom(N, size = size, mu = mu),
    gaussian = stats::rnorm(N, mu, sd),
    binomial = stats::rbinom(N, size = n_trials, prob = mu))

  list(y = y, data = data,
       truth = list(beta = beta, response = response,
                    size = size, sd = sd, trials = n_trials,
                    link = switch(response, gaussian = "identity",
                                  binomial = "logit", "log")))
}

#' Simulate Royle-Nichols occupancy data
#'
#' Latent abundance `N_i ~ Poisson(lambda_i)`, `log lambda = X beta_lambda`, and
#' per-visit detection `y_ij ~ Bernoulli(1 - (1 - r_i)^{N_i})` with per-individual
#' detection `logit r_i = beta_r` (site-level). Matches the [royle_nichols()]
#' family; the response is an `N x J` 0/1 detection-history matrix.
#'
#' @param N Number of sites (default 200).
#' @param J Number of visits per site (default 5).
#' @param beta_lambda Log-abundance coefficients `c(intercept, slope_on_x)`.
#' @param beta_r Per-individual detection logit (a scalar intercept).
#' @param seed Random seed.
#' @return A list with `y` (`N x J` 0/1 matrix), `data`, and `truth`
#'   (`beta_lambda`, `beta_r`, realised abundance `N`, per-site `lambda` / `r`).
#' @examples
#' sim <- simulate_royle_nichols(N = 100, J = 4, seed = 1)
#' dim(sim$y)
#' @export
simulate_royle_nichols <- function(N = 200, J = 5,
                                   beta_lambda = c(0.3, 0.5),
                                   beta_r = -0.8,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data   <- data.frame(x = rnorm(N))
  X_l    <- model.matrix(~ x, data)
  lambda <- exp(as.vector(X_l %*% beta_lambda))
  r      <- plogis(rep(beta_r[1], N))

  Ni <- rpois(N, lambda)
  y  <- matrix(0L, N, J)
  for (i in seq_len(N)) {
    p_i <- 1 - (1 - r[i])^Ni[i]
    y[i, ] <- rbinom(J, 1, p_i)
  }
  list(y = y, data = data,
       truth = list(beta_lambda = beta_lambda, beta_r = beta_r,
                    N = Ni, lambda = lambda, r = r))
}

#' Simulate temporal (multi-season) occupancy data
#'
#' Colonization and extinction are constant across the `n_seasons - 1` transition
#' intervals by default. Supplying `beta_gamma` and/or `beta_epsilon` makes the
#' corresponding rate SEASON-VARYING: a per-`(site, interval)` covariate is drawn
#' into an `[N x (n_seasons - 1)]` matrix column (`gamma_cov` / `eps_cov`) of the
#' returned `data`, and the rate is `plogis(beta[1] + beta[2] * cov)`. Fit these
#' with `colonization = ~ gamma_cov` / `extinction = ~ eps_cov` (the matrix
#' column drives the interval-indexed rate, gcol33/tulpaObs#124).
#'
#' @param N Number of sites (default 100).
#' @param J Number of visits per season (default 4).
#' @param n_seasons Number of seasons (default 5).
#' @param beta_occ Initial occupancy coefficients.
#' @param beta_det Detection coefficients.
#' @param gamma Colonization probability (default 0.2); ignored when
#'   `beta_gamma` is given.
#' @param epsilon Extinction probability (default 0.1); ignored when
#'   `beta_epsilon` is given.
#' @param beta_gamma Optional `c(intercept, slope)` for a season-varying
#'   colonization logit driven by a drawn per-`(site, interval)` covariate.
#' @param beta_epsilon Optional `c(intercept, slope)` for a season-varying
#'   extinction logit driven by a drawn per-`(site, interval)` covariate.
#' @param seed Random seed.
#' @return A list with `y` (3D array), `data`, and `truth`.
#' @export
simulate_dyn_occu <- function(N = 100, J = 4, n_seasons = 5,
                    beta_occ = c(0.5), beta_det = c(0),
                    gamma = 0.2, epsilon = 0.1,
                    beta_gamma = NULL, beta_epsilon = NULL,
                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  data <- data.frame(x = rnorm(N))
  psi1 <- plogis(beta_occ[1])
  p <- plogis(beta_det[1])
  n_int <- n_seasons - 1L

  # Season-varying rate matrices [N x (T-1)] when covariate coefficients are
  # given; otherwise the scalar rate broadcast over intervals. The covariate
  # draws happen only in the season-varying branch, so the constant-rate RNG
  # stream (and its output) is unchanged.
  gamma_cov <- eps_cov <- NULL
  gam_it <- matrix(gamma,   N, n_int)
  eps_it <- matrix(epsilon, N, n_int)
  if (!is.null(beta_gamma)) {
    gamma_cov <- matrix(rnorm(N * n_int), N, n_int)
    gam_it <- plogis(beta_gamma[1] + beta_gamma[2] * gamma_cov)
    data$gamma_cov <- gamma_cov
  }
  if (!is.null(beta_epsilon)) {
    eps_cov <- matrix(rnorm(N * n_int), N, n_int)
    eps_it <- plogis(beta_epsilon[1] + beta_epsilon[2] * eps_cov)
    data$eps_cov <- eps_cov
  }

  z <- matrix(NA_integer_, N, n_seasons)
  z[, 1] <- rbinom(N, 1, psi1)
  for (t in 2:n_seasons) {
    iv <- t - 1L
    z[, t] <- z[, t-1] * (1 - rbinom(N, 1, eps_it[, iv])) +
              (1 - z[, t-1]) * rbinom(N, 1, gam_it[, iv])
  }

  y <- array(NA_integer_, dim = c(N, J, n_seasons))
  for (i in seq_len(N)) {
    for (t in seq_len(n_seasons)) {
      y[i, , t] <- rbinom(J, 1, z[i, t] * p)
    }
  }

  list(
    y = y,
    data = data,
    truth = list(
      psi1 = psi1, p = p, gamma = gamma, epsilon = epsilon,
      beta_gamma = beta_gamma, beta_epsilon = beta_epsilon,
      z = z, beta_occ = beta_occ, beta_det = beta_det
    )
  )
}

#' Simulate integrated (multi-source) occupancy data
#'
#' @param N_total Total number of unique sites (default 150).
#' @param n_data Number of data sources (default 2).
#' @param J Vector of visits per source (default c(4, 3)).
#' @param n_shared Number of sites shared across sources (default 20).
#' @param beta_occ Occupancy coefficients (default c(0.5, 0.3)).
#' @param beta_det List of detection coefficient vectors per source.
#' @param seed Random seed.
#' @return A list with `y` (list of matrices), `data`, `site_maps`, and `truth`.
#' @export
simulate_int_occu <- function(N_total = 150, n_data = 2, J = c(4, 3),
                      n_shared = 20,
                      beta_occ = c(0.5, 0.3),
                      beta_det = list(c(0.2, -0.4), c(-0.1, 0.3)),
                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(J) != n_data) J <- rep(J[1], n_data)
  if (length(beta_det) != n_data) {
    stop("beta_det must be a list of length n_data")
  }

  data <- data.frame(x = rnorm(N_total))
  X_occ <- model.matrix(~ x, data)
  psi <- plogis(as.vector(X_occ %*% beta_occ))
  z <- rbinom(N_total, 1, psi)

  # Assign sites to sources
  shared <- seq_len(n_shared)
  remaining <- setdiff(seq_len(N_total), shared)
  n_per_source <- (N_total - n_shared) %/% n_data
  site_maps <- vector("list", n_data)
  y_list <- vector("list", n_data)

  for (s in seq_len(n_data)) {
    start <- (s - 1) * n_per_source + 1
    end <- min(s * n_per_source, length(remaining))
    source_sites <- sort(c(shared, remaining[start:end]))
    site_maps[[s]] <- source_sites

    ns <- length(source_sites)
    p_s <- plogis(beta_det[[s]][1])
    y_s <- matrix(NA_integer_, ns, J[s])
    for (i in seq_len(ns)) {
      y_s[i, ] <- rbinom(J[s], 1, z[source_sites[i]] * p_s)
    }
    y_list[[s]] <- y_s
  }

  list(
    y = y_list,
    data = data,
    site_maps = site_maps,
    truth = list(beta_occ = beta_occ, beta_det = beta_det, psi = psi, z = z)
  )
}

#' Simulate temporal multi-species occupancy data
#'
#' @param N Number of sites (default 50).
#' @param J Visits per season (default 3).
#' @param n_species Number of species (default 5).
#' @param n_seasons Number of seasons (default 4).
#' @param beta_comm_mean Community mean for occupancy (default c(0)).
#' @param beta_comm_sd Community SD for occupancy (default c(0.5)).
#' @param gamma Colonization probability (default 0.15).
#' @param epsilon Extinction probability (default 0.1).
#' @param field Optional per-site shared areal field (length `N`) added to the
#'   first-season occupancy logit of every species -- the shared field of the
#'   community dynamic-spatial model (stMsPGOcc). Default `NULL` (no field).
#' @param seed Random seed.
#' @return A list with `y` (4D array), `data`, and `truth`.
#' @seealso [ms_dyn_occu()], the family this simulates for.
#' @export
simulate_ms_dyn_occu <- function(N = 50, J = 3, n_species = 5, n_seasons = 4,
                      beta_comm_mean = c(0), beta_comm_sd = c(0.5),
                      gamma = 0.15, epsilon = 0.1,
                      field = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  data <- data.frame(x = rnorm(N))
  # Draw the per-species first-season occupancy logits (same RNG draw as before;
  # plogis of an rnorm), then apply the optional shared field on the logit scale.
  logit_psi1_species <- rnorm(n_species, beta_comm_mean[1], beta_comm_sd[1])
  psi1_species <- plogis(logit_psi1_species)
  p_species <- plogis(rnorm(n_species, 0, 0.5))
  if (!is.null(field) && length(field) != N) {
    stop("simulate_ms_dyn_occu(): `field` must have length N.", call. = FALSE)
  }

  z <- array(NA_integer_, dim = c(N, n_seasons, n_species))
  y <- array(NA_integer_, dim = c(N, J, n_seasons, n_species))

  for (sp in seq_len(n_species)) {
    psi1_i <- if (is.null(field)) rep(psi1_species[sp], N)
              else plogis(logit_psi1_species[sp] + field)
    z[, 1, sp] <- rbinom(N, 1, psi1_i)
    for (t in 2:n_seasons) {
      z[, t, sp] <- z[, t-1, sp] * (1 - rbinom(N, 1, epsilon)) +
                    (1 - z[, t-1, sp]) * rbinom(N, 1, gamma)
    }
    for (i in seq_len(N)) {
      for (t in seq_len(n_seasons)) {
        y[i, , t, sp] <- rbinom(J, 1, z[i, t, sp] * p_species[sp])
      }
    }
  }

  list(
    y = y,
    data = data,
    truth = list(
      psi1_species = psi1_species, p_species = p_species,
      gamma = gamma, epsilon = epsilon, z = z, field = field
    )
  )
}

#' Simulate integrated multi-species occupancy data
#'
#' @param N Number of sites (default 100).
#' @param J Vector of visits per source (default c(3, 4)).
#' @param n_species Number of species (default 5).
#' @param n_data Number of data sources (default 2).
#' @param seed Random seed.
#' @return A list with `y` (list of 3D arrays), `data`, and `truth`.
#' @seealso [ms_int_occu()], the family this simulates for.
#' @export
simulate_ms_int_occu <- function(N = 100, J = c(3, 4), n_species = 5,
                        n_data = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(J) != n_data) J <- rep(J[1], n_data)

  data <- data.frame(x = rnorm(N))
  psi_species <- plogis(rnorm(n_species, 0, 0.5))
  z <- matrix(NA_integer_, N, n_species)
  for (sp in seq_len(n_species)) z[, sp] <- rbinom(N, 1, psi_species[sp])

  y_list <- vector("list", n_data)
  p_det <- vector("list", n_data)
  for (s in seq_len(n_data)) {
    p_s <- plogis(rnorm(n_species, 0, 0.3))
    p_det[[s]] <- p_s
    y_s <- array(NA_integer_, dim = c(N, J[s], n_species))
    for (sp in seq_len(n_species)) {
      for (i in seq_len(N)) {
        y_s[i, , sp] <- rbinom(J[s], 1, z[i, sp] * p_s[sp])
      }
    }
    y_list[[s]] <- y_s
  }

  list(
    y = y_list,
    data = data,
    truth = list(psi_species = psi_species, p_det = p_det, z = z)
  )
}
