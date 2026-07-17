# =============================================================================
# occu.R — Internal model builders for tobs()
#
# `.tobs_build_model()` is the data-binding constructor: it consumes formulas
# + data + response and returns a `tobs_model` structure with the design
# matrices, response, and process metadata that the engine needs. Dispatched
# from `tobs()` via the per-family `.dispatch_*` helpers.
# =============================================================================


#' Build a tobs model object
#'
#' Inferred model type from arguments:
#' - **Single-season**: no `col_formula`, no `species`
#' - **Dynamic**: `col_formula` and/or `ext_formula` provided
#' - **Integrated**: `integrated = TRUE`, `y` a list of matrices
#'
#' Community families (ms_occu / ms_dyn_occu / ms_int_occu / ms_occu_cover /
#' ms_count / jsdm) have their own in-tree binders + Laplace-EM fitter and do not
#' pass through here.
#'
#' @keywords internal
.tobs_build_model <- function(occ_formula, det_formula = NULL, data, y,
                              col_formula = NULL, ext_formula = NULL,
                              species = NULL, integrated = FALSE,
                              abundance = FALSE, count = FALSE,
                              count_response = "poisson", count_trials = NULL,
                              det_visit_formula = NULL, det_visit_data = NULL) {

  is_dynamic   <- !is.null(col_formula) || !is.null(ext_formula)
  is_integrated <- isTRUE(integrated)
  is_abundance <- isTRUE(abundance)
  is_count     <- isTRUE(count)

  if (is_abundance) {
    if (is.null(det_formula)) stop("det_formula required for N-mixture models")
    return(.tobs_build_abun(occ_formula, det_formula, data, y,
                            det_visit_formula, det_visit_data))
  }
  if (is_count)      return(.tobs_build_count(occ_formula, data, y, count_response,
                                              trials = count_trials))
  if (is_integrated) return(.tobs_build_integrated(occ_formula, det_formula, data, y))
  if (is_dynamic)    return(.tobs_build_dynamic(occ_formula, det_formula, data, y,
                                                col_formula, ext_formula))

  if (is.null(det_formula)) stop("det_formula required for occupancy models")
  .tobs_build_single(occ_formula, det_formula, data, y,
                     det_visit_formula, det_visit_data)
}


# ============================================================================
# Per-model-type builders
# ============================================================================

.tobs_build_single <- function(occ_formula, det_formula, data, y,
                               det_visit_formula = NULL, det_visit_data = NULL) {
  if (!is.matrix(y)) {
    stop("y must be a matrix (n_sites x max_visits)")
  }
  .tobs_check_site_count(nrow(y), nrow(data), "rows")

  bind  <- .tobs_bind_formulas(list(psi = occ_formula, p = det_formula), data)
  X_occ <- model.matrix(bind$fe$psi, data)
  X_det <- model.matrix(bind$fe$p, data)

  X_det_visit <- .tobs_build_visit_X(det_visit_formula, det_visit_data,
                                     nrow(y), ncol(y), arm = "detection")

  y_int <- matrix(as.integer(y), nrow = nrow(y), ncol = ncol(y))
  y_int[is.na(y_int)] <- -1L

  n_detected <- sum(apply(y_int, 1, function(row) any(row[row >= 0] == 1)))

  structure(list(
    model_type = "single",
    y = y_int,
    X_processes = list(X_occ, X_det),
    X_det_visit = X_det_visit,
    formulas = list(occ = bind$fe$psi, det = bind$fe$p),
    structured_terms = bind$terms,
    data = data,
    n_sites = nrow(y),
    max_visits = ncol(y),
    process_info = list(
      list(name = "psi", p = ncol(X_occ), coef_names = colnames(X_occ)),
      list(name = "p",   p = ncol(X_det), coef_names = colnames(X_det))
    ),
    det_visit_names = if (!is.null(X_det_visit)) colnames(X_det_visit) else character(0),
    naive_occ = n_detected / nrow(y)
  ), class = "tobs_model")
}

.tobs_build_dynamic <- function(occ_formula, det_formula, data, y,
                                col_formula, ext_formula) {
  if (is.null(col_formula)) col_formula <- ~ 1
  if (is.null(ext_formula)) ext_formula <- ~ 1

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

  .tobs_check_site_count(n_sites, nrow(data), "sites")

  bind  <- .tobs_bind_formulas(
    list(psi1 = occ_formula, p = det_formula,
         gamma = col_formula, epsilon = ext_formula), data)
  X_occ <- model.matrix(bind$fe$psi1, data)
  # Detection may vary by primary season (gcol33/tulpaObs#124): a covariate
  # supplied as a [n_sites x T] matrix column of `data` drives per-season
  # detection, unrolled long-form over (site, season); a plain per-site covariate
  # keeps the site-level design (byte-identical to the constant-detection path).
  det_ad <- .tobs_season_arm_design(bind$fe$p, data, n_sites, n_seasons,
                                    "detection", fam = "dyn_occu")
  X_det <- det_ad$X

  # Colonization (gamma) and extinction (epsilon) span the T-1 transition
  # intervals. A rate that varies by interval is supplied as a [n_sites x (T-1)]
  # matrix column of `data`; the arm design is then long-form over
  # (site, interval), otherwise it is the site-level design and the fit is
  # byte-identical to the constant-rate path (gcol33/tulpaObs#124, the dyn_abun
  # #80 recipe). One transition interval per pair of adjacent seasons.
  n_intervals <- n_seasons - 1L
  col_ad <- .tobs_interval_arm_design(bind$fe$gamma, data, n_sites, n_intervals,
                                      "colonization", fam = "dyn_occu")
  ext_ad <- .tobs_interval_arm_design(bind$fe$epsilon, data, n_sites, n_intervals,
                                      "extinction", fam = "dyn_occu")
  X_col <- col_ad$X
  X_ext <- ext_ad$X

  # y_flat layout is site-major: y_flat[i*T*K + t*K + j] (0-indexed)
  # = y_flat[(i-1)*T*K + (t-1)*K + j] (1-indexed). This matches what
  # src/dyn_occ_likelihood.h and R/laplace.R::build_dynamic_callbacks
  # expect. Achieved by permuting the 3D y[site, visit, season] array to
  # [visit, season, site] so that column-major flattening makes site the
  # slowest-varying dimension.
  y_int <- as.integer(aperm(y, c(2, 3, 1)))
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
    # Site-level frame, kept for the same reason the single-season binder keeps
    # it: a structured term resolves its node index / group_var against it (an
    # areal `group_var`, a spatial bar's grouping factor), and those columns
    # need not appear in any arm's fixed-effect formula.
    data = data,
    y = y,
    y_flat = y_int,
    n_visits = n_visits,
    any_detected = any_detected,
    X_processes = list(X_occ, X_det, X_col, X_ext),
    formulas = list(occ = bind$fe$psi1, det = bind$fe$p,
                    col = bind$fe$gamma, ext = bind$fe$epsilon),
    structured_terms = bind$terms,
    n_sites = n_sites,
    n_seasons = n_seasons,
    n_intervals = n_intervals,
    col_season_varying = col_ad$season_varying,
    ext_season_varying = ext_ad$season_varying,
    det_season_varying = det_ad$season_varying,
    max_visits = max_visits,
    process_info = list(
      list(name = "psi1",    p = ncol(X_occ), coef_names = colnames(X_occ)),
      list(name = "p",       p = ncol(X_det), coef_names = colnames(X_det)),
      list(name = "gamma",   p = ncol(X_col), coef_names = colnames(X_col)),
      list(name = "epsilon", p = ncol(X_ext), coef_names = colnames(X_ext))
    )
  ), class = "tobs_model")
}

.tobs_build_integrated <- function(occ_formula, det_formula, data, y) {
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
  # Structured terms (spatial / re / temporal) are supported on the shared
  # occupancy field only; the per-source detection formulas are fixed-effects.
  occ_bind <- .tobs_bind_formulas(list(psi = occ_formula), data)
  X_occ    <- model.matrix(occ_bind$fe$psi, data)

  y_sources <- vector("list", n_sources)
  X_det_list <- vector("list", n_sources)
  det_fe     <- vector("list", n_sources)
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

    det_parsed <- .tobs_parse_formula(det_formulas[[s]], data = data)
    if (length(det_parsed$terms)) {
      stop("Structured terms on integrated detection sources are not ",
           "supported; place spatial / re / temporal terms on the ",
           "occupancy formula.", call. = FALSE)
    }
    det_fe[[s]] <- det_parsed$fe_formula
    src_rows <- site_maps[[s]] + 1L
    X_det_list[[s]] <- model.matrix(det_fe[[s]], data[src_rows, , drop = FALSE])
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
    # Site-level frame: a structured term resolves its node index / group_var
    # against it, and those columns need not appear in any arm's formula.
    data = data,
    y_sources = y_sources,
    site_maps = site_maps,
    X_processes = X_processes,
    formulas = c(list(occ = occ_bind$fe$psi), setNames(det_fe, names(y))),
    structured_terms = occ_bind$terms,
    n_sites = n_sites,
    n_sources = n_sources,
    process_info = process_info
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# count() -- a GLMM on the observed count / continuous response directly (no
# detection, no latent state; the abundance analogue of the JSDM). One value per
# site; a log-link Poisson / negative-binomial or an identity-link Gaussian. The
# response family travels on the model as `response`; the negbin size / Gaussian
# residual variance is estimated by an outer dispersion loop in .dispatch_count
# (tulpa_laplace takes a fixed phi per fit).
# ---------------------------------------------------------------------------
.tobs_build_count <- function(occ_formula, data, y, response = "poisson",
                              trials = NULL) {
  response <- match.arg(response, c("poisson", "negbin", "gaussian", "binomial"))
  if (is.matrix(y)) {
    if (ncol(y) != 1L) {
      stop("count(): y must be a vector or a one-column matrix (one value ",
           "per site).", call. = FALSE)
    }
    y <- as.vector(y)
  }
  if (!is.numeric(y)) {
    stop("count(): y must be a numeric vector (one value per site).",
         call. = FALSE)
  }
  n_data <- if (is.data.frame(data)) nrow(data) else length(y)
  .tobs_check_site_count(length(y), n_data, "sites")

  is_count_fam <- response %in% c("poisson", "negbin")
  is_binom     <- identical(response, "binomial")

  # Per-site trial count for the binomial response. A scalar recycles; a vector
  # is one trial count per site (before NA dropping). trials is meaningless for
  # the other responses -- ignore it there rather than error, so a stray default
  # never blocks a Poisson fit.
  if (is_binom) {
    n_trials <- if (is.null(trials)) rep(1L, length(y)) else as.numeric(trials)
    if (length(n_trials) == 1L) n_trials <- rep(n_trials, length(y))
    if (length(n_trials) != length(y)) {
      stop("count(response = \"binomial\"): `trials` must be a scalar or one ",
           "value per site (length ", length(y), ").", call. = FALSE)
    }
    if (any(!is.na(n_trials) & (n_trials < 1 | abs(n_trials - round(n_trials)) >
                                1e-8))) {
      stop("count(response = \"binomial\"): `trials` must be positive integers.",
           call. = FALSE)
    }
  }

  # Complete-case: drop sites with a missing response (and their design rows).
  # For a count response NA is genuinely missing (0 is a real count), so it is
  # dropped rather than coerced.
  bind  <- .tobs_bind_formulas(list(psi = occ_formula), data)
  X_occ <- model.matrix(bind$fe$psi, data)
  keep  <- !is.na(y)
  if (is_binom) keep <- keep & !is.na(n_trials)
  if (!all(keep)) {
    y     <- y[keep]
    X_occ <- X_occ[keep, , drop = FALSE]
    if (is_binom) n_trials <- n_trials[keep]
  }
  n_sites <- length(y)

  if (is_count_fam &&
      (any(y < 0) || any(abs(y - round(y)) > 1e-8))) {
    stop(sprintf(paste0("count(response = \"%s\"): y must be non-negative ",
                        "integer counts."), response), call. = FALSE)
  }
  if (is_binom) {
    if (any(y < 0) || any(abs(y - round(y)) > 1e-8)) {
      stop("count(response = \"binomial\"): y must be non-negative integer ",
           "success counts.", call. = FALSE)
    }
    if (any(y > n_trials)) {
      stop("count(response = \"binomial\"): every success count y must be <= ",
           "its trial count (0 <= k <= n).", call. = FALSE)
    }
  }

  link    <- if (identical(response, "gaussian")) "identity"
             else if (is_binom) "logit" else "log"
  y_store <- if (identical(response, "gaussian")) as.numeric(y)
             else as.integer(round(y))

  structure(list(
    model_type = "count",
    y_count = y_store,
    n_trials = if (is_binom) as.integer(round(n_trials)) else NULL,
    response = response,
    link = link,
    X_processes = list(X_occ),
    X_occ = X_occ,
    formulas = list(occ = bind$fe$psi),
    structured_terms = bind$terms,
    data = if (is.data.frame(data)) data[keep, , drop = FALSE] else data,
    n_sites = n_sites,
    N = n_sites,
    process_info = list(
      list(name = "mu", p = ncol(X_occ), coef_names = colnames(X_occ),
           link = link)
    )
  ), class = "tobs_model")
}


# ============================================================================
# Print method
# ============================================================================

#' @export
print.tobs_model <- function(x, ...) {
  type_label <- switch(x$model_type,
    single = "Single-season occupancy model",
    dynamic = "Multi-season dynamic occupancy model",
    integrated = sprintf("Integrated occupancy model (%d sources)", x$n_sources)
  )
  cat(type_label, "\n")

  if (x$model_type == "single") {
    cat(sprintf("  Sites: %d, Max visits: %d\n", x$n_sites, x$max_visits))
  } else if (x$model_type == "dynamic") {
    cat(sprintf("  Sites: %d, Seasons: %d, Max visits: %d\n",
                x$n_sites, x$n_seasons, x$max_visits))
  } else if (x$model_type == "integrated") {
    cat(sprintf("  Sites: %d, Sources: %d\n", x$n_sites, x$n_sources))
  }

  for (pi in x$process_info) {
    cat(sprintf("  %s covariates (%d): %s\n",
                pi$name, pi$p, paste(pi$coef_names, collapse = ", ")))
  }

  if (x$model_type == "single" && !is.null(x$naive_occ)) {
    cat(sprintf("  Naive occupancy: %.1f%%\n", 100 * x$naive_occ))
  }

  invisible(x)
}
