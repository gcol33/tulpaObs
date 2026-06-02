# =============================================================================
# aggregation_scan.R - Suggest spatial cell size and yearly clustering so a
# single-season occupancy model is identifiable.
#
# Single-visit plot data (one record per plot per year) carries no within-unit
# replication, so psi and p are confounded until records are pooled. Pooling
# happens on two axes: plots into a spatial cell (shared-psi assumption) and
# years into a contiguous block (closure assumption). An occasion is then any
# record landing in a (cell, block) bucket, so the same bucket-counting code
# covers single-visit plots, repeat-visited plots, and the no-pooling limit.
#
# `occu_aggregation_scan()` scans candidate (cell size x year block) pairs and
# scores each by identifiability: structural replication ("count") or the
# curvature of the constant-model occupancy likelihood ("info", the smallest
# eigenvalue / posterior SEs of the 2x2 (logit psi, logit p) information).
# =============================================================================


# ---------------------------------------------------------------------------
# Spatial gridding
# ---------------------------------------------------------------------------

# Assign each record to a square cell of edge `size`, keyed by integer grid
# coordinates relative to `origin`. Returns a character cell id per record.
.scan_cell_assign <- function(x, y, size, origin) {
  ix <- floor((x - origin[1]) / size)
  iy <- floor((y - origin[2]) / size)
  paste(ix, iy, sep = "_")
}

# Median nearest-neighbour distance among unique plot locations. Caps the
# pairwise computation at `n_sample` sampled locations so the O(n^2) distance
# matrix stays bounded on large plot networks.
.scan_nn_spacing <- function(x, y, n_sample = 1500L) {
  loc <- unique(cbind(x, y))
  n <- nrow(loc)
  if (n < 2L) return(NA_real_)
  if (n > n_sample) {
    idx <- seq_len(n)
    keep <- idx[round(seq(1, n, length.out = n_sample))]
    loc <- loc[keep, , drop = FALSE]
  }
  d <- as.matrix(stats::dist(loc))
  diag(d) <- Inf
  stats::median(apply(d, 1, min))
}

# Auto sequence of candidate cell sizes: a geometric ladder from a size that
# merges nearest neighbours up to one spanning a fraction of the extent.
.scan_auto_cell_sizes <- function(x, y, n, nn_sample) {
  ext <- max(diff(range(x)), diff(range(y)))
  if (!is.finite(ext) || ext <= 0) return(numeric(0))
  nn <- .scan_nn_spacing(x, y, nn_sample)
  lo <- if (is.finite(nn) && nn > 0) nn else ext / 100
  hi <- ext / 2
  if (hi <= lo) hi <- lo * 4
  signif(exp(seq(log(lo), log(hi), length.out = n)), 3)
}


# ---------------------------------------------------------------------------
# Temporal segmentation
# ---------------------------------------------------------------------------

# Fixed contiguous blocks of `L` consecutive calendar years (by value, so gaps
# in the year axis are respected). Returns an integer block id per record.
.scan_fixed_blocks <- function(year, L, origin_year) {
  as.integer((year - origin_year) %/% L)
}

# Mean-shift segmentation of the per-year naive occupancy series by exact
# optimal partitioning (Jackson et al. 2005): minimize within-segment SSE plus
# `penalty` per added changepoint, subject to a minimum segment length. Returns
# the 1-based segment id for each entry of the sorted `years` vector.
.scan_changepoints <- function(rate, min_seg, penalty) {
  n <- length(rate)
  if (n <= 1L) return(rep(1L, n))
  cs  <- c(0, cumsum(rate))
  cs2 <- c(0, cumsum(rate^2))
  seg_cost <- function(i, j) {
    # SSE of rate[i..j]
    m <- j - i + 1L
    s <- cs[j + 1L] - cs[i]
    s2 <- cs2[j + 1L] - cs2[i]
    s2 - s * s / m
  }
  F  <- rep(Inf, n + 1L)
  bk <- integer(n + 1L)
  F[1L] <- -penalty
  for (t in seq_len(n)) {
    lo <- 1L
    hi <- t - min_seg + 1L
    if (hi < lo) next
    for (s in lo:hi) {
      cand <- F[s] + seg_cost(s, t) + penalty
      if (cand < F[t + 1L]) {
        F[t + 1L] <- cand
        bk[t + 1L] <- s - 1L
      }
    }
  }
  # backtrack
  seg <- integer(n)
  e <- n
  k <- 0L
  while (e > 0L) {
    s <- bk[e + 1L] + 1L
    seg[s:e] <- k
    k <- k + 1L
    e <- s - 1L
  }
  # contiguous runs in chronological order, relabelled 1, 2, ...
  as.integer(cumsum(c(TRUE, diff(seg) != 0)))
}

# Default segmentation penalty: a modified-BIC scale set from a robust noise
# estimate (MAD of first differences of the yearly rate), so the number of
# changepoints adapts to the series' own variability.
.scan_cp_penalty <- function(rate, mult) {
  n <- length(rate)
  if (n <= 2L) return(Inf)
  dd <- diff(rate)
  sigma2 <- (stats::mad(dd) / sqrt(2))^2
  if (!is.finite(sigma2) || sigma2 <= 0) sigma2 <- stats::var(rate)
  if (!is.finite(sigma2) || sigma2 <= 0) sigma2 <- 1e-6
  mult * sigma2 * log(n)
}


# ---------------------------------------------------------------------------
# Constant-model occupancy identifiability (the "info" score)
# ---------------------------------------------------------------------------

# Negative log-likelihood of the constant (intercept-only) single-season
# occupancy model on the logit scale. Each unit contributes K occasions with d
# detections: a detected unit gives psi * p^d (1-p)^(K-d); an undetected unit
# gives psi (1-p)^K + (1-psi). Units with K == 0 contribute nothing.
.scan_occu_negll <- function(par, K, d) {
  psi <- stats::plogis(par[1L])
  p   <- stats::plogis(par[2L])
  detected <- d > 0L
  ll <- 0
  if (any(detected)) {
    Kd <- K[detected]; dd <- d[detected]
    ll <- ll + sum(log(psi) + dd * log(p) + (Kd - dd) * log1p(-p))
  }
  und <- !detected & K > 0L
  if (any(und)) {
    Ku <- K[und]
    ll <- ll + sum(log(psi * (1 - p)^Ku + (1 - psi)))
  }
  -ll
}

# Fit the constant model and return identifiability diagnostics: posterior SEs
# of psi and p (delta method off the logit-scale Hessian), the smallest
# eigenvalue and condition number of the 2x2 observed information, and whether
# the information is well-conditioned.
.scan_fit_const <- function(K, d, max_condition, optim_control) {
  use <- K > 0L
  K <- K[use]; d <- d[use]
  out <- list(se_psi = NA_real_, se_p = NA_real_, min_eig = NA_real_,
              cond = NA_real_, psi_hat = NA_real_, p_hat = NA_real_,
              converged = FALSE, identifiable = FALSE)
  if (length(K) == 0L || !any(d > 0L) || !any(K >= 2L)) return(out)

  fit <- tryCatch(
    stats::optim(c(0, 0), .scan_occu_negll, K = K, d = d,
                 method = "BFGS", hessian = TRUE, control = optim_control),
    error = function(e) NULL)
  if (is.null(fit)) return(out)

  out$converged <- isTRUE(fit$convergence == 0L)
  psi <- stats::plogis(fit$par[1L]); p <- stats::plogis(fit$par[2L])
  out$psi_hat <- psi; out$p_hat <- p

  H <- fit$hessian
  ev <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  out$min_eig <- min(ev)
  out$cond <- if (min(ev) > 0) max(ev) / min(ev) else Inf

  V <- tryCatch(solve(H), error = function(e) NULL)
  if (!is.null(V) && all(diag(V) > 0)) {
    se_logit <- sqrt(diag(V))
    out$se_psi <- se_logit[1L] * psi * (1 - psi)
    out$se_p   <- se_logit[2L] * p * (1 - p)
  }
  out$identifiable <- out$converged && is.finite(out$cond) &&
    out$cond < max_condition && min(ev) > 0 &&
    is.finite(out$se_p) && is.finite(out$se_psi)
  out
}


# ---------------------------------------------------------------------------
# Per-candidate scoring
# ---------------------------------------------------------------------------

# Build (cell, block) units for one candidate and compute its statistics.
# `cell` / `block` are per-record ids; `resp` is the 0/1 response; `plot` is an
# optional per-record plot id used for the within-cell homogeneity proxy.
.scan_score_candidate <- function(cell, block, year, resp, plot,
                                   score, max_condition, optim_control) {
  unit <- paste(cell, block, sep = "@")
  uf <- factor(unit)
  K <- as.integer(table(uf))
  d <- as.integer(tapply(resp, uf, sum))
  detected <- d > 0L

  n_units <- length(K)
  n_rep <- sum(K >= 2L)
  n_det <- sum(detected)
  naive_psi <- mean(detected)
  naive_p <- if (sum(K[detected]) > 0) sum(d[detected]) / sum(K[detected]) else NA_real_

  # closure-violation proxy: within-block spread of the per-year detection rate,
  # averaged over blocks (high = occupancy moving within a supposed closed block)
  bf <- factor(block)
  turnover <- mean(tapply(seq_along(resp), bf, function(ix) {
    yr <- factor(year[ix])
    rate_by_year <- tapply(resp[ix], yr, mean)
    if (length(rate_by_year) <= 1L) 0 else stats::sd(rate_by_year)
  }), na.rm = TRUE)

  # within-cell homogeneity proxy: spread of per-plot detection rate inside a
  # cell, averaged over multi-plot cells (high = heterogeneous sites pooled)
  heterogeneity <- NA_real_
  if (!is.null(plot)) {
    cf <- factor(cell)
    hv <- tapply(seq_along(resp), cf, function(ix) {
      pf <- factor(plot[ix])
      if (nlevels(pf) <= 1L) return(NA_real_)
      stats::sd(tapply(resp[ix], pf, mean))
    })
    heterogeneity <- mean(hv, na.rm = TRUE)
  }

  row <- list(
    n_units = n_units, n_units_rep = n_rep, frac_rep = n_rep / n_units,
    mean_K = mean(K), n_detected = n_det, naive_psi = naive_psi,
    naive_p = naive_p, turnover = turnover, heterogeneity = heterogeneity,
    se_psi = NA_real_, se_p = NA_real_, min_eig = NA_real_, cond = NA_real_,
    converged = NA, identifiable = NA)

  if (identical(score, "info")) {
    f <- .scan_fit_const(K, d, max_condition, optim_control)
    row$se_psi <- f$se_psi; row$se_p <- f$se_p
    row$min_eig <- f$min_eig; row$cond <- f$cond
    row$converged <- f$converged; row$identifiable <- f$identifiable
  } else {
    # structural necessary conditions only
    row$identifiable <- n_rep > 0L && n_det > 0L &&
      isTRUE(naive_p > 0 && naive_p < 1)
  }
  row
}


# ---------------------------------------------------------------------------
# Front door
# ---------------------------------------------------------------------------

#' Suggest spatial and temporal aggregation for an identifiable occupancy model
#'
#' Scans candidate spatial cell sizes and yearly clusterings of long-format,
#' typically single-visit plot data and scores each by how well a single-season
#' occupancy model separates occupancy (`psi`) from detection (`p`). Replication
#' is manufactured by pooling plots into a spatial cell (assumes plots in a cell
#' share `psi`) and years into a contiguous block (assumes occupancy is closed
#' across the block); an occasion is any record in a (cell, block) bucket, so
#' single-visit plots, repeat-visited plots, and the no-pooling limit all flow
#' through one code path.
#'
#' Two scoring modes:
#' - `"info"` fits the constant (intercept-only) occupancy model per candidate
#'   and reports the posterior SEs of `psi`/`p` and the smallest eigenvalue and
#'   condition number of the 2x2 `(logit psi, logit p)` observed information. A
#'   confounding ridge shows up directly as a near-zero eigenvalue.
#' - `"count"` checks only the structural necessary conditions (units with two
#'   or more occasions, detected units, a non-degenerate naive detection rate).
#'
#' The recommended candidate is the **least-pooling** identifiable one (smallest
#' mean occasions per unit, ties broken by largest information eigenvalue): pool
#' only as much as identifiability requires. When no candidate is identifiable
#' the data cannot support a standard occupancy model and a Royle-Nichols or
#' count model is the alternative.
#'
#' @param data Long data.frame, one row per plot-year record (or plot-year-visit).
#' @param response Name of the 0/1 detection column.
#' @param coords Length-2 character vector naming the x and y coordinate columns.
#' @param year Name of the integer year column.
#' @param plot Optional name of the plot-id column; enables the within-cell
#'   homogeneity proxy.
#' @param cell_sizes Numeric vector of candidate cell edge lengths. `NULL`
#'   auto-proposes a geometric ladder from the nearest-neighbour spacing to half
#'   the extent.
#' @param block_lengths Integer vector of candidate contiguous block lengths (in
#'   years). `NULL` runs mean-shift changepoint segmentation of the per-year
#'   naive occupancy series and uses the resulting blocks as the single temporal
#'   candidate.
#' @param score `"info"` (curvature-based, default) or `"count"` (structural).
#' @param family Observation family. Only `"occupancy"` is supported; the
#'   per-family constant-model scorer is the extension point for others.
#' @param control List of tuning knobs: `n_cell_sizes` (auto ladder length, 6),
#'   `nn_sample` (NN-spacing sample cap, 1500), `min_block` (min years per
#'   segment, 1), `cp_penalty` (explicit changepoint penalty, else auto),
#'   `cp_penalty_mult` (auto-penalty scale, 1), `max_condition` (conditioning
#'   cap for the identifiable flag, 1e6), `optim_control` (passed to `optim`).
#' @return A `tobs_aggregation_scan` object: `candidates` (scored data.frame),
#'   `recommended` (chosen row or `NULL`), `segmentation` (auto changepoint
#'   blocks, if used), and scan metadata.
#' @export
occu_aggregation_scan <- function(data, response, coords, year, plot = NULL,
                                  cell_sizes = NULL, block_lengths = NULL,
                                  score = c("info", "count"),
                                  family = c("occupancy"),
                                  control = list()) {
  score <- match.arg(score)
  family <- match.arg(family)
  if (!is.data.frame(data)) stop("data must be a data.frame")
  if (length(coords) != 2L) stop("coords must name two columns (x, y)")
  for (col in c(response, coords, year, plot)) {
    if (!is.null(col) && !col %in% names(data))
      stop(sprintf("column '%s' not found in data", col))
  }

  ctrl <- modifyList(list(
    n_cell_sizes = 6L, nn_sample = 1500L, min_block = 1L,
    cp_penalty = NULL, cp_penalty_mult = 1, max_condition = 1e6,
    optim_control = list(maxit = 200L)), control)

  resp <- as.integer(data[[response]])
  if (!all(is.na(resp) | resp %in% c(0L, 1L)))
    stop("response must be 0/1 (occurrence)")
  x  <- as.numeric(data[[coords[1L]]])
  y  <- as.numeric(data[[coords[2L]]])
  yr <- as.integer(data[[year]])
  pl <- if (!is.null(plot)) as.character(data[[plot]]) else NULL

  ok <- !(is.na(resp) | is.na(x) | is.na(y) | is.na(yr))
  resp <- resp[ok]; x <- x[ok]; y <- y[ok]; yr <- yr[ok]
  if (!is.null(pl)) pl <- pl[ok]
  if (length(resp) == 0L) stop("no complete records after dropping NA")

  origin <- c(min(x), min(y))
  origin_year <- min(yr)

  if (is.null(cell_sizes)) {
    cell_sizes <- .scan_auto_cell_sizes(x, y, ctrl$n_cell_sizes, ctrl$nn_sample)
    if (length(cell_sizes) == 0L)
      stop("could not auto-propose cell sizes; supply cell_sizes")
  }

  # Temporal candidates: a named list of per-record block-id vectors plus a
  # human-readable definition string.
  blocks <- list()
  segmentation <- NULL
  if (is.null(block_lengths)) {
    yrs <- sort(unique(yr))
    rate <- as.numeric(tapply(resp, factor(yr, levels = yrs),
                              function(v) mean(v > 0L)))
    penalty <- if (!is.null(ctrl$cp_penalty)) ctrl$cp_penalty
               else .scan_cp_penalty(rate, ctrl$cp_penalty_mult)
    seg <- .scan_changepoints(rate, ctrl$min_block, penalty)
    year2seg <- stats::setNames(seg, yrs)
    blk <- as.integer(year2seg[as.character(yr)])
    segmentation <- data.frame(year = yrs, segment = seg, naive_psi = rate)
    blocks[["changepoint"]] <- list(id = blk, def = "changepoint")
  } else {
    for (L in block_lengths) {
      blk <- .scan_fixed_blocks(yr, as.integer(L), origin_year)
      blocks[[paste0("L", L)]] <- list(id = blk, def = paste0(L, "yr"))
    }
  }

  rows <- list()
  for (s in cell_sizes) {
    cell <- .scan_cell_assign(x, y, s, origin)
    for (bn in names(blocks)) {
      blk <- blocks[[bn]]
      sc <- .scan_score_candidate(cell, blk$id, yr, resp, pl,
                                  score, ctrl$max_condition, ctrl$optim_control)
      rows[[length(rows) + 1L]] <- c(
        list(cell_size = s, block = blk$def), sc)
    }
  }

  cand <- do.call(rbind, lapply(rows, function(r)
    as.data.frame(r, stringsAsFactors = FALSE)))
  rownames(cand) <- NULL

  # Recommendation: least pooling among identifiable candidates.
  rec <- NULL
  idok <- which(.scan_is_true(cand$identifiable))
  if (length(idok) > 0L) {
    sub <- cand[idok, , drop = FALSE]
    ord <- order(sub$mean_K, -.na0(sub$min_eig))
    rec <- sub[ord[1L], , drop = FALSE]
    rownames(rec) <- NULL
  }

  structure(list(
    candidates = cand,
    recommended = rec,
    segmentation = segmentation,
    score = score, family = family,
    n_records = length(resp),
    n_years = length(unique(yr)),
    year_range = range(yr),
    coord_extent = c(diff(range(x)), diff(range(y))),
    cell_sizes = cell_sizes
  ), class = "tobs_aggregation_scan")
}

# TRUE only where x is logical TRUE (NA -> FALSE); robust for the count mode
# where identifiable may be NA on degenerate candidates.
.scan_is_true <- function(x) !is.na(x) & x

.na0 <- function(x) ifelse(is.finite(x), x, 0)


# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

#' @export
print.tobs_aggregation_scan <- function(x, ...) {
  cat("Occupancy aggregation scan\n")
  cat(sprintf("  Records: %d | Years: %d (%d-%d) | Score: %s\n",
              x$n_records, x$n_years, x$year_range[1L], x$year_range[2L], x$score))
  cat(sprintf("  Candidates: %d (%d cell sizes x %d temporal blockings)\n",
              nrow(x$candidates), length(x$cell_sizes),
              nrow(x$candidates) / max(length(x$cell_sizes), 1L)))
  n_id <- sum(.scan_is_true(x$candidates$identifiable))
  cat(sprintf("  Identifiable candidates: %d / %d\n", n_id, nrow(x$candidates)))

  if (!is.null(x$segmentation)) {
    nseg <- length(unique(x$segmentation$segment))
    cat(sprintf("  Auto year segmentation: %d block(s)\n", nseg))
  }

  if (!is.null(x$recommended)) {
    r <- x$recommended
    cat("\n  Recommended (least pooling that identifies):\n")
    cat(sprintf("    cell size = %s | years = %s | occasions/unit = %.2f | units = %d\n",
                format(r$cell_size), r$block, r$mean_K, r$n_units))
    if (x$score == "info" && is.finite(r$se_p)) {
      cat(sprintf("    SE(psi) = %.3f | SE(p) = %.3f | info min-eig = %.3g\n",
                  r$se_psi, r$se_p, r$min_eig))
    }
  } else {
    cat("\n  No identifiable candidate found.\n")
    cat("    The data cannot support a standard occupancy model at any scanned\n")
    cat("    aggregation; consider a Royle-Nichols or count (N-mixture) model.\n")
  }
  invisible(x)
}

#' Plot an occupancy aggregation scan
#'
#' Heatmap of the identifiability score over the scanned (cell size x temporal
#' block) grid: information smallest-eigenvalue for `score = "info"`, fraction
#' of replicated units for `score = "count"`. Identifiable cells are outlined.
#'
#' @param x A `tobs_aggregation_scan` object.
#' @param ... Ignored.
#' @return Invisible `NULL`.
#' @importFrom graphics image axis box rect
#' @export
plot.tobs_aggregation_scan <- function(x, ...) {
  cand <- x$candidates
  cs <- sort(unique(cand$cell_size))
  bl <- unique(cand$block)
  val <- matrix(NA_real_, length(cs), length(bl),
                dimnames = list(format(cs), bl))
  idm <- matrix(FALSE, length(cs), length(bl))
  metric <- if (x$score == "info") cand$min_eig else cand$frac_rep
  for (i in seq_len(nrow(cand))) {
    ri <- match(cand$cell_size[i], cs)
    ci <- match(cand$block[i], bl)
    val[ri, ci] <- metric[i]
    idm[ri, ci] <- .scan_is_true(cand$identifiable[i])
  }
  graphics::image(seq_along(cs), seq_along(bl), val,
                  axes = FALSE, xlab = "cell size", ylab = "year block",
                  main = sprintf("Identifiability (%s)", x$score),
                  col = grDevices::hcl.colors(20, "YlGnBu", rev = TRUE))
  graphics::axis(1, at = seq_along(cs), labels = format(cs), las = 2)
  graphics::axis(2, at = seq_along(bl), labels = bl, las = 1)
  graphics::box()
  for (i in seq_along(cs)) for (j in seq_along(bl)) {
    if (idm[i, j])
      graphics::rect(i - 0.5, j - 0.5, i + 0.5, j + 0.5,
                     border = "firebrick", lwd = 2)
  }
  invisible(NULL)
}
