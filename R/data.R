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

#' Convert long-format data to occupancy format
#'
#' @param df Data.frame in long format (one row per site-visit).
#' @param y Character, name of detection column (0/1/NA).
#' @param site Character, name of site identifier column.
#' @param visit Character, name of visit/replicate column.
#' @param occ.covs Character vector of site-level covariate names.
#' @param det.covs Character vector of visit-level covariate names.
#' @param coords Character vector of length 2 for coordinate columns.
#' @return An `tobs_data` object.
#' @export
tobs_data <- function(df, y, site, visit,
                      occ.covs = NULL, det.covs = NULL,
                      coords = NULL) {
  if (!is.data.frame(df)) stop("df must be a data.frame")
  for (col in c(y, site, visit)) {
    if (!col %in% names(df)) stop(sprintf("column '%s' not found in df", col))
  }

  sites <- unique(df[[site]])
  n_sites <- length(sites)
  visits <- sort(unique(df[[visit]]))
  max_visits <- length(visits)

  # Build detection history matrix
  y_mat <- matrix(NA_integer_, n_sites, max_visits)
  rownames(y_mat) <- sites
  for (i in seq_len(nrow(df))) {
    si <- match(df[[site]][i], sites)
    vi <- match(df[[visit]][i], visits)
    y_mat[si, vi] <- as.integer(df[[y]][i])
  }

  # Extract site-level covariates
  occ_df <- NULL
  if (!is.null(occ.covs)) {
    # Take first occurrence per site
    site_rows <- match(sites, df[[site]])
    occ_df <- df[site_rows, occ.covs, drop = FALSE]
    rownames(occ_df) <- NULL
  }

  # Extract visit-level detection covariates
  det_list <- NULL
  if (!is.null(det.covs)) {
    det_list <- list()
    for (dc in det.covs) {
      mat <- matrix(NA_real_, n_sites, max_visits)
      for (i in seq_len(nrow(df))) {
        si <- match(df[[site]][i], sites)
        vi <- match(df[[visit]][i], visits)
        mat[si, vi] <- as.numeric(df[[dc]][i])
      }
      det_list[[dc]] <- mat
    }
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
print.tobs_data <- function(x, ...) {
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
#' @return A list with `y`, `data`, `coords`, and `truth`.
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

#' Simulate temporal (multi-season) occupancy data
#'
#' @param N Number of sites (default 100).
#' @param J Number of visits per season (default 4).
#' @param n_seasons Number of seasons (default 5).
#' @param beta_occ Initial occupancy coefficients.
#' @param beta_det Detection coefficients.
#' @param gamma Colonization probability (default 0.2).
#' @param epsilon Extinction probability (default 0.1).
#' @param seed Random seed.
#' @return A list with `y` (3D array), `data`, and `truth`.
#' @export
simulate_dyn_occu <- function(N = 100, J = 4, n_seasons = 5,
                    beta_occ = c(0.5), beta_det = c(0),
                    gamma = 0.2, epsilon = 0.1,
                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  data <- data.frame(x = rnorm(N))
  psi1 <- plogis(beta_occ[1])
  p <- plogis(beta_det[1])

  z <- matrix(NA_integer_, N, n_seasons)
  z[, 1] <- rbinom(N, 1, psi1)
  for (t in 2:n_seasons) {
    z[, t] <- z[, t-1] * (1 - rbinom(N, 1, epsilon)) +
              (1 - z[, t-1]) * rbinom(N, 1, gamma)
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
#' @param seed Random seed.
#' @return A list with `y` (4D array), `data`, and `truth`.
#' @export
simulate_dyn_ms_occu <- function(N = 50, J = 3, n_species = 5, n_seasons = 4,
                      beta_comm_mean = c(0), beta_comm_sd = c(0.5),
                      gamma = 0.15, epsilon = 0.1,
                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  data <- data.frame(x = rnorm(N))
  psi1_species <- plogis(rnorm(n_species, beta_comm_mean[1], beta_comm_sd[1]))
  p_species <- plogis(rnorm(n_species, 0, 0.5))

  z <- array(NA_integer_, dim = c(N, n_seasons, n_species))
  y <- array(NA_integer_, dim = c(N, J, n_seasons, n_species))

  for (sp in seq_len(n_species)) {
    z[, 1, sp] <- rbinom(N, 1, psi1_species[sp])
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
      gamma = gamma, epsilon = epsilon, z = z
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
#' @export
simulate_int_ms_occu <- function(N = 100, J = c(3, 4), n_species = 5,
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
