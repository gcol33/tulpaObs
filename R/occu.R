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
#' A structured term enters the process whose formula it is written in. For an
#' integrated model this makes `detection = ~ x + spde(lon, lat, ...)` a
#' **detection-arm** field: the occupancy state stays whatever the occupancy
#' formula says, and the continuous Matern field sits on the detection logit.
#' The detection arm there is one arm observed through S sources, so the field
#' is one structure fit once per source block -- source `s` fits its own
#' realization at its own sites, and two sources covering the same location each
#' carry their own value there, the same way each source carries its own
#' detection coefficients. `fit$spatial_field_det` is the per-source named list
#' of mesh fields (source names as given on `y`); project a source's field to
#' its sites with `fit$spatial$tulpa_spec$A`.
#'
#' The arm accepts a continuous `spde()` field under `method = "laplace"`. The
#' areal kinds and the temporal / re / svc / latent classes are grid-integrated
#' as latent blocks on the state arm, so they are written on the occupancy
#' formula; on a detection formula they error with a pointer rather than being
#' fit against the wrong arm.
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

# Gate the structured terms an integrated model's DETECTION arm may carry.
#
# The state arm is the wide one: every spatial / temporal / re / svc / latent
# term the occupancy formula can carry reaches its own fitter through
# `.tobs_structures_from_model()`, and is untouched here. The detection arm is
# reached by exactly one route -- the per-source detection blocks of the
# single-Laplace EM (`build_integrated_callbacks`), which attaches a continuous
# Matern (SPDE) mesh field to each source's binomial block. Every other term
# class, and the areal spatial kinds, would be built into the multi-block latent
# prior that the nested-Laplace path attaches to the STATE block, so a
# detection-arm term there would be fit against the wrong arm rather than the
# one it was written on. Those error here with a pointer instead.
#
# A single realization shared across both arms (a term written on the occupancy
# formula and `copy()`d onto detection, or a second field on the other arm) is
# the same limit the single-season path has: the block fitter fits one field
# realization per submodel block, so one shared realization is not expressible.
.tobs_validate_integrated_terms <- function(terms) {
  if (!length(terms)) return(invisible())
  on_arm <- function(t, k) k %in% t$processes
  sp_occ <- sp_det <- 0L
  for (t in terms) {
    spec <- t$spec
    is_sp <- inherits(spec, "tobs_spatial")
    if (on_arm(t, 1L)) sp_occ <- sp_occ + is_sp
    if (!on_arm(t, 2L)) next
    sp_det <- sp_det + is_sp
    label <- spec$label %||% class(spec)[1L]
    if (on_arm(t, 1L)) {
      stop("A structured term shared across the occupancy and detection arms ",
           "of an integrated model is not supported (each submodel block fits ",
           "its own field realization, so one shared realization is not ",
           "expressible). Write the term on one formula.", call. = FALSE)
    }
    if (!is_sp) {
      stop(sprintf(paste0(
        "`%s` on an integrated detection formula is not supported; the ",
        "detection arm carries a continuous spde() Matern field only. Place ",
        "re() / temporal() / svc() / latent() terms on the occupancy formula, ",
        "or use method = \"nuts\" for a family that samples them."), label),
        call. = FALSE)
    }
    if (!identical(spec$type, "spde")) {
      stop(sprintf(paste0(
        "A '%s' field on an integrated detection formula is not supported; ",
        "the detection arm is fit by the single-Laplace EM (method = ",
        "\"laplace\"), which carries a continuous spde() Matern field. The ",
        "areal kinds are grid-integrated on the state arm under method = ",
        "\"nested_laplace\" -- write the areal term on the occupancy formula."),
        spec$type), call. = FALSE)
    }
  }
  if (sp_det > 1L) {
    stop("Only one spatial term is supported on the integrated detection arm.",
         call. = FALSE)
  }
  if (sp_det > 0L && sp_occ > 0L) {
    stop("An integrated model carries one spatial field: a state field (on the ",
         "occupancy formula) or a detection field (on the detection formula), ",
         "not both. Fit the arms' fields separately, or use one arm.",
         call. = FALSE)
  }
  invisible()
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
  # A structured term's process is the formula it is written in, here as
  # everywhere else: a term on the occupancy formula is the shared state field,
  # a term on the detection formula is a detection-arm field. The detection arm
  # of an integrated model is one arm spread over S source blocks -- every
  # source measures the same latent occupancy state through its own detection
  # process -- so the arm is bound as a single process and each source block
  # carries its own realization of that field at its own sites. That needs ONE
  # detection formula, which is what "the same detection model at every source"
  # means; per-source formulas that differ carry fixed effects only. Sameness is
  # a property of the formula, not of the object: two sources written the same
  # way are one arm even though each `~` call carries its own environment, so
  # the comparison is on the deparsed formula.
  det_key <- vapply(det_formulas,
                    function(f) paste(deparse(f), collapse = " "), character(1))
  det_shared <- all(det_key == det_key[1L])
  bind <- .tobs_bind_formulas(
    if (det_shared) list(psi = occ_formula, p = det_formulas[[1L]])
    else list(psi = occ_formula),
    data)
  .tobs_validate_integrated_terms(bind$terms)
  X_occ <- model.matrix(bind$fe$psi, data)

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

    if (det_shared) {
      det_fe[[s]] <- bind$fe$p
    } else {
      det_parsed <- .tobs_parse_formula(det_formulas[[s]], data = data)
      if (length(det_parsed$terms)) {
        stop("A structured term on a per-source detection formula needs the ",
             "same detection formula at every source (the detection arm is ",
             "one arm, fit as one field realization per source block). Pass a ",
             "single `detection` formula carrying the term, or drop the term ",
             "from the per-source formulas.", call. = FALSE)
      }
      det_fe[[s]] <- det_parsed$fe_formula
    }
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
    formulas = c(list(occ = bind$fe$psi), setNames(det_fe, names(y))),
    structured_terms = bind$terms,
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
