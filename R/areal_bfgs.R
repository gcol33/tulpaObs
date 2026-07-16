# =============================================================================
# areal_bfgs.R - shared nested-Laplace driver for areal fields on a gradient-only
# family marginal (gcol33/tulpaObs#51).
#
# Families whose marginal exposes an analytic gradient but NO analytic per-site
# Hessian (the open N-mixture dyn_abun, the false-positive occupancy fp_occu)
# carry an areal field on one arm by: an outer grid over the field hyperparameters;
# per cell, BFGS over (fixed params, field params) with the family's analytic
# gradient + the field log-prior; the Laplace marginal from a finite-difference
# Hessian of that gradient at the mode (the observed-information route the
# families' non-spatial fits use).
#
# The driver owns ALL field handling through a `field` spec, so a new field
# structure (ICAR / proper-CAR single block, or BYM2 two-block v/w) is a new spec
# rather than a copied driver. The family supplies only
#   eval(theta_fixed, offset) -> list(log_lik, grad_fixed, grad_eta)
# where `offset` is the per-site eta offset on the field arm and `grad_eta` the
# per-site d log L / d eta_field. The field spec maps its parameters to that
# per-site offset and scatters grad_eta back, owns the prior, the sum-to-zero
# constraint, and the reported field. Single source of truth for the areal-BFGS
# families.
# =============================================================================

.areal_Q <- function(adj, rho) { deg <- rowSums(adj != 0); diag(deg) - rho * adj }

# Inverse Cholesky of the FIXED field precision tau Q(rho) for the non-centered
# areal-field NUTS path (gcol33/tulpaObs#72): z = Linv %*% raw, raw ~ N(0, I) has
# covariance (tau Q(rho))^{-1}. A small ridge proper-ises an intrinsic ICAR
# (proper-CAR is already full rank). Returns the n x n inverse-Cholesky matrix
# (the field block passes it to the C++ sampler as field_Linv); errors if the
# precision is not positive definite. Single source of truth shared by every
# count / occupancy family's areal NUTS fitter (abun, removal, distance, fp_occu,
# dyn_abun).
.tobs_field_linv <- function(adj, tau, rho, n, ridge = 1e-4) {
  Q  <- .areal_Q(adj, rho)
  Qr <- tau * Q + diag(ridge * tau, n)
  L  <- tryCatch(chol(Qr), error = function(e) NULL)   # upper: L'L = Qr
  if (is.null(L)) stop("areal NUTS: fixed field precision not positive definite.",
                       call. = FALSE)
  backsolve(L, diag(n))                                # (L)^{-1}; z = Linv %*% raw
}

# Whitened-field loading L for the non-centered areal-field NUTS path
# (gcol33/tulpaObs#71): z = L %*% raw, raw ~ N(0, I), Cov(z) = (tau Q(rho))^{+}.
# `type` selects the field. A proper-CAR field is full rank -> the square inverse
# Cholesky (n x n; z covers the whole space). An intrinsic icar / bym2 field has
# the constant vector in the precision null space, so a square whitening leaves a
# flat field-mean direction that maxes the NUTS tree depth; the sum-to-zero
# reparameterisation drops that direction by keeping only the non-null eigenpairs
# of tau Q -> L is n x (n - 1), z is automatically centred (sum z = 0). The
# eigen-loading L = U_+ diag(1 / sqrt(tau lambda_+)) satisfies L L' = (tau Q)^{+}
# restricted to the sum-to-zero subspace. bym2 (Riebler 2016) scales its
# structured ICAR block by sigma sqrt(rho / scale_factor) and adds an
# unstructured iid block sigma sqrt(1 - rho) (full-rank, square); the structured
# block is the same eigen-loading.
.tobs_field_load <- function(adj, type, tau, rho, n, ridge = 1e-4, tol = 1e-8) {
  if (identical(type, "car_proper"))
    return(.tobs_field_linv(adj, tau, rho, n, ridge))
  # Intrinsic ICAR precision tau Q (rho = 1); reduced eigen-loading on the
  # non-null subspace (drop the constant direction).
  Q  <- .areal_Q(adj, 1.0)
  ev <- eigen(tau * Q, symmetric = TRUE)
  keep <- ev$values > tol * max(ev$values)             # non-null eigenpairs
  U <- ev$vectors[, keep, drop = FALSE]
  d <- ev$values[keep]
  U %*% diag(1 / sqrt(d), nrow = length(d))             # n x (n - 1)
}

# Whitened-field loading + fixed hyper for a non-centered areal-field NUTS path
# (gcol33/tulpaObs#71, #113). Single source of truth shared by every observation-
# family spatial NUTS fitter (abun / ms_abun / removal / distance / fp_occu /
# dyn_abun), so the intrinsic icar / bym2 sum-to-zero reparameterisation is
# derived once. The fixed precision (tau, rho) and, for bym2, the marginal SD
# sigma are supplied by the caller from its own nested-Laplace warm-start (the
# families carry them under different names). Returns:
#   * field_load: the n x n_raw loading L (z = L %*% raw, raw ~ N(0, I_{n_raw}))
#   * tau, rho:   the fixed precision hyperparameters (tau = NA for bym2, whose
#                 amplitude rides sigma inside the two-block loading)
#   * n_raw:      ncol(field_load) -- n for car_proper, n - 1 for icar,
#                 2n - 1 for bym2.
# car_proper is the square inverse Cholesky; icar is the sum-to-zero eigen-loading
# (drops the constant null direction); bym2 stacks the structured (centred ICAR,
# scaled by sigma sqrt(rho / scale_factor)) and unstructured (iid sigma
# sqrt(1 - rho)) blocks columnwise (Riebler 2016).
.tobs_nuts_field_loading <- function(adj, type, n, tau = NA_real_, rho = NA_real_,
                                     sigma = NA_real_, scale_factor = NULL) {
  if (identical(type, "bym2")) {
    sigma <- max(if (is.finite(sigma)) sigma else sqrt(max(tau, 1e-3, na.rm = TRUE)),
                 1e-3)
    rho   <- min(max(if (is.finite(rho)) rho else 0.5, 0.01), 0.99)
    sf    <- scale_factor %||% compute_bym2_scale(adj)
    Lstr  <- .tobs_field_load(adj, "icar", 1, 1, n)      # centred ICAR basis
    a <- sigma * sqrt(rho / sf); b <- sigma * sqrt(1 - rho)
    L <- cbind(a * Lstr, b * diag(n))                    # [structured | iid]
    return(list(field_load = L, tau = NA_real_, rho = rho, n_raw = ncol(L)))
  }
  tau <- max(tau, 1e-3)
  rho <- if (identical(type, "car_proper")) min(max(rho, 0.01), 0.99) else 1.0
  L   <- .tobs_field_load(adj, type, tau, rho, n)
  list(field_load = L, tau = tau, rho = rho, n_raw = ncol(L))
}

# Whitened-field loading + fixed hyper for a non-centered TEMPORAL-field NUTS path
# (gcol33/tulpaObs#114). The temporal analogue of `.tobs_nuts_field_loading`: the
# fixed precision is tau Q(rho) with Q the ar1 / rw1 / rw2 / iid structure matrix
# (`.tobs_temporal_Q`), and the whitened loading is the reduced eigen-loading over
# the non-null eigenpairs of tau Q (z = L %*% raw, raw ~ N(0, I), Cov(z) =
# (tau Q)^{+}). Full-rank ar1 / iid keep all T eigenpairs (square L, n_raw = T);
# rank-deficient rw1 / rw2 drop the 1 / 2 null directions (sum-to-zero, n_raw =
# T - 1 / T - 2) exactly like the intrinsic-icar areal case. Uniform across all
# four types, so the C++ field block (which only sees an n_field_units x n_raw
# loading + a per-site field_map) needs no temporal-specific branch.
.tobs_nuts_temporal_loading <- function(type, T, tau = NA_real_, rho = NA_real_,
                                        tol = 1e-8) {
  tau <- max(tau, 1e-3)
  Q  <- .tobs_temporal_Q(type, T, rho = if (identical(type, "ar1")) rho else NULL)
  ev <- eigen(tau * Q, symmetric = TRUE)
  keep <- ev$values > tol * max(ev$values)              # non-null eigenpairs
  U <- ev$vectors[, keep, drop = FALSE]
  d <- ev$values[keep]
  L <- U %*% diag(1 / sqrt(d), nrow = length(d))         # T x n_raw
  list(field_load = L, tau = tau, rho = rho, n_raw = ncol(L))
}

# ICAR / proper-CAR single-block field (parameter z, length n_sp; eta += z[map]).
.areal_field_car <- function(adj, kind, map, n_sp) {
  tau_grid <- exp(seq(log(0.3), log(30), length.out = 9L))
  rho_grid <- if (kind == "car_proper") seq(0.1, 0.95, length.out = 6L) else 1.0
  # One cell from physical hyperparameters: ICAR -> (tau); proper CAR -> (tau, rho).
  make_cell <- function(theta) {
    tau <- theta[1L]
    rho <- if (kind == "car_proper") theta[2L] else 1.0
    ldQ <- if (kind == "icar") 0 else {
      ch <- tryCatch(chol(.areal_Q(adj, rho)), error = function(e) NULL)
      if (is.null(ch)) -Inf else 2 * sum(log(diag(ch)))
    }
    list(tau = tau, rho = rho, ldQ = ldQ, Q = .areal_Q(adj, rho))
  }
  cells <- list()
  for (rho in rho_grid) for (tau in tau_grid)
    cells[[length(cells) + 1L]] <- make_cell(if (kind == "car_proper") c(tau, rho) else tau)
  axes <- c(list(.tobs_ccd_axis("tau", "log", lower = 0.3, upper = 30, start = 3)),
            if (kind == "car_proper")
              list(.tobs_ccd_axis("rho", "identity", lower = 0.1, upper = 0.95, start = 0.5)))
  list(
    n_field = n_sp, n_sp = n_sp, cells = cells, axes = axes, make_cell = make_cell,
    valid = function(cell) is.finite(cell$ldQ),
    offset = function(fp, cell) fp[map],
    scatter = function(grad_eta) {
      g <- numeric(n_sp); for (s in seq_along(map)) g[map[s]] <- g[map[s]] + grad_eta[s]; g
    },
    prior_logp = function(fp, cell) {
      quad <- as.numeric(t(fp) %*% cell$Q %*% fp)
      if (kind == "icar") -0.5 * cell$tau * quad + 0.5 * (n_sp - 1) * log(cell$tau)
      else 0.5 * cell$ldQ + 0.5 * n_sp * log(cell$tau) - 0.5 * cell$tau * quad
    },
    prior_grad = function(fp, cell) cell$tau * as.numeric(cell$Q %*% fp),  # d(-logp)/dfp
    center = function(fp) if (kind == "icar") fp - mean(fp) else fp,
    constrain = if (kind == "icar") rep(TRUE, n_sp) else rep(FALSE, n_sp),
    to_phi = function(fp, cell) fp,
    # Physical field hyperparameters per cell, for fixing the precision tau Q(rho)
    # on the NUTS path (gcol33/tulpaObs#72). ICAR pins rho = 1.
    type = if (kind == "icar") "icar" else "car_proper",
    to_hyper = function(cell) c(tau = cell$tau, rho = cell$rho)
  )
}

# BYM2 two-block field (v = ICAR, w = iid; eta += a v[map] + b w[map],
# a = sigma sqrt(rho/scale), b = sigma sqrt(1-rho)). The (v, w) priors are
# independent of (sigma, rho) (Riebler 2016), so the prior is constant across
# cells; the hyperparameters enter only the likelihood through (a, b).
.areal_field_bym2 <- function(adj, scale_factor, map, n_sp) {
  sigma_grid <- exp(seq(log(0.2), log(3), length.out = 5L))
  rho_grid   <- c(0.05, 0.3, 0.5, 0.7, 0.95)
  Q <- .areal_Q(adj, 1.0)                       # ICAR precision for v
  make_cell <- function(theta) {
    sg <- theta[1L]; rho <- theta[2L]
    list(sigma = sg, rho = rho,
         a = sg * sqrt(rho / scale_factor), b = sg * sqrt(1 - rho))
  }
  cells <- list()
  for (sg in sigma_grid) for (rho in rho_grid)
    cells[[length(cells) + 1L]] <- make_cell(c(sg, rho))
  axes <- list(
    .tobs_ccd_axis("sigma", "log",      lower = 0.2,  upper = 3,    start = 0.77),
    .tobs_ccd_axis("rho",   "identity", lower = 0.05, upper = 0.95, start = 0.5))
  scat1 <- function(grad_eta) {
    g <- numeric(n_sp); for (s in seq_along(map)) g[map[s]] <- g[map[s]] + grad_eta[s]; g
  }
  list(
    n_field = 2L * n_sp, n_sp = n_sp, cells = cells, axes = axes, make_cell = make_cell,
    valid = function(cell) cell$sigma > 0 && cell$rho >= 0 && cell$rho <= 1,
    offset = function(fp, cell) {
      v <- fp[seq_len(n_sp)]; w <- fp[n_sp + seq_len(n_sp)]
      cell$a * v[map] + cell$b * w[map]
    },
    scatter = function(grad_eta, cell) {
      s <- scat1(grad_eta); c(cell$a * s, cell$b * s)
    },
    prior_logp = function(fp, cell) {
      v <- fp[seq_len(n_sp)]; w <- fp[n_sp + seq_len(n_sp)]
      -0.5 * as.numeric(t(v) %*% Q %*% v) - 0.5 * sum(w^2)
    },
    prior_grad = function(fp, cell) {
      v <- fp[seq_len(n_sp)]; w <- fp[n_sp + seq_len(n_sp)]
      c(as.numeric(Q %*% v), w)
    },
    center = function(fp) { fp[seq_len(n_sp)] <- fp[seq_len(n_sp)] - mean(fp[seq_len(n_sp)]); fp },
    constrain = c(rep(TRUE, n_sp), rep(FALSE, n_sp)),  # sum-to-zero on v only
    to_phi = function(fp, cell) cell$a * fp[seq_len(n_sp)] + cell$b * fp[n_sp + seq_len(n_sp)],
    type = "bym2",
    to_hyper = function(cell) c(sigma = cell$sigma, rho = cell$rho),
    bym2 = TRUE
  )
}

# Areal-BFGS nested-Laplace fit over one OR several latent field blocks (#78).
#
# `field` is a single field spec (the historical single-block call) or a LIST of
# field specs (e.g. spatial + temporal). Each block owns a contiguous slice of the
# concatenated field-parameter vector and supplies the same closure interface
# (offset / scatter / prior_logp / prior_grad / center / constrain / to_phi /
# cells / to_hyper). The per-cell BFGS optimises (fixed params, field_1, ...,
# field_K); a cell is a tuple of per-block cells; the outer grid is the Cartesian
# product of the blocks' grids. A length-1 list is numerically identical to the
# single-block kernel. CCD outer integration applies only to the single-block
# case; multi-block uses the product grid.
# Gate a temporal() term on a count family (removal / distance / fp_occu /
# dyn_abun). A temporal AR1/RW1/RW2/iid field composes with the areal field on the
# arm under method = "laplace" / "nested_laplace" via the shared areal-BFGS driver
# (gcol33/tulpaObs#78). A temporal term WITHOUT a spatial field, or under NUTS, is
# not wired; those raise a clear error here so the family dispatch can call the
# spatial fitter unconditionally once the gate passes.
.tobs_check_count_temporal <- function(temporal, spatial, method, family, arm,
                                       allow_temporal_only = FALSE,
                                       allow_nuts_temporal = FALSE) {
  # NUTS + temporal is wired only where a fixed-hyper non-centered temporal field
  # rides the family's NUTS field block (dyn_abun; gcol33/tulpaObs#114). It runs
  # temporal-only (no simultaneous areal field) on that path.
  if (identical(method, "nuts") && isTRUE(allow_nuts_temporal)) {
    if (!is.null(spatial))
      stop(sprintf(paste0("%s() NUTS supports a temporal() field on its own (no ",
                          "simultaneous areal field); combine areal + temporal ",
                          "under method = \"nested_laplace\". (tulpaObs#114)"),
                   family), call. = FALSE)
    return(invisible(TRUE))
  }
  if (identical(method, "nuts"))
    stop(sprintf(paste0("%s() does not support a temporal() term under method = ",
                        "\"nuts\"; the temporal field composes with the areal ",
                        "field on the %s arm under method = \"nested_laplace\". ",
                        "(tulpaObs#78)"), family, arm), call. = FALSE)
  # A temporal-only field (no areal term) is wired on families whose spatial
  # fitter builds the areal-BFGS block list from either term (gcol33/tulpaObs#114).
  if (is.null(spatial) && !isTRUE(allow_temporal_only))
    stop(sprintf(paste0("%s() supports a temporal() term composed WITH an areal ",
                        "field on the %s arm (e.g. icar()/car_proper()/bym2() + ",
                        "temporal()) under method = \"nested_laplace\"; a temporal ",
                        "term on its own is not yet wired. (tulpaObs#78)"),
                 family, arm), call. = FALSE)
  invisible(TRUE)
}

# Structure precision matrix Q (tau = 1, marginal-precision scaled) of a temporal
# block over T points: ar1(rho) tridiagonal full-rank; rw1 tridiagonal (rank T-1);
# rw2 pentadiagonal (rank T-2); iid identity. Single source of truth for the
# temporal-block prior used alongside the areal field on the count families (#78).
.tobs_temporal_Q <- function(type, T, rho = NULL, cyclic = FALSE) {
  if (type == "iid") return(diag(1, T))
  if (type == "ar1") {
    if (is.null(rho) || abs(rho) >= 1) stop("ar1 rho must be in (-1, 1).", call. = FALSE)
    Q <- matrix(0, T, T)
    Q[1, 1] <- 1; Q[T, T] <- 1
    if (T > 2) for (t in 2:(T - 1)) Q[t, t] <- 1 + rho^2
    for (t in seq_len(T - 1)) { Q[t, t + 1] <- -rho; Q[t + 1, t] <- -rho }
    return(Q / (1 - rho^2))
  }
  if (type == "rw1") {
    Q <- matrix(0, T, T)
    if (cyclic) {
      for (t in seq_len(T)) {
        Q[t, t] <- 2
        nt <- if (t == T) 1L else t + 1L; pt <- if (t == 1L) T else t - 1L
        Q[t, nt] <- Q[t, nt] - 1; Q[t, pt] <- Q[t, pt] - 1
      }
    } else {
      Q[1, 1] <- 1; Q[1, 2] <- -1; Q[T, T] <- 1; Q[T, T - 1] <- -1
      if (T > 2) for (t in 2:(T - 1)) { Q[t, t] <- 2; Q[t, t - 1] <- -1; Q[t, t + 1] <- -1 }
    }
    return(Q)
  }
  if (type == "rw2") {
    if (T < 4) stop("rw2 needs at least 4 time points.", call. = FALSE)
    # Second-difference operator D (T-2 x T); Q = D' D, rank T - 2.
    D <- matrix(0, T - 2L, T)
    for (t in seq_len(T - 2L)) { D[t, t] <- 1; D[t, t + 1L] <- -2; D[t, t + 2L] <- 1 }
    return(crossprod(D))
  }
  stop(sprintf("Unsupported temporal type '%s'.", type), call. = FALSE)
}

# Temporal latent-field block for the areal-BFGS driver (#78). Mirrors
# `.areal_field_car`: a single block over `n_t` time points, eta += z[map] with
# `map = time_idx` (the per-site time index), the prior `0.5 (rank) log tau -
# 0.5 tau z' Q z` (rank-deficient rw1/rw2 use rank = T - deficiency and a
# sum-to-zero constraint; ar1/iid are full rank, no constraint). The grid is
# (tau) for rw1/rw2/iid and (tau, rho) for ar1.
.tobs_temporal_field <- function(type, map, n_t, cyclic = FALSE,
                                  tau_grid = NULL, rho_grid = NULL) {
  type <- match.arg(type, c("ar1", "rw1", "rw2", "iid"))
  # Narrow per-block grids keep the product with the spatial block's grid small
  # (the count-family areal field already sweeps ~9 tau, or tau x rho for proper-
  # CAR / BYM2). A 3-point tau (x 2 rho for AR1) mirrors the occupancy multi-block
  # path's narrowed defaults; users can override via the temporal() term.
  if (is.null(tau_grid)) tau_grid <- exp(seq(log(0.5), log(20), length.out = 3L))
  is_ar1 <- identical(type, "ar1")
  if (is_ar1 && is.null(rho_grid)) rho_grid <- c(0.4, 0.8)
  deficiency <- switch(type, rw1 = 1L, rw2 = 2L, 0L)
  rank_eff   <- n_t - deficiency
  constrain  <- rep(deficiency > 0L, n_t)         # sum-to-zero on rw1 / rw2 only

  make_cell <- function(theta) {
    tau <- theta[1L]
    rho <- if (is_ar1) theta[2L] else NA_real_
    Q   <- .tobs_temporal_Q(type, n_t, rho = if (is_ar1) rho else NULL, cyclic = cyclic)
    ldQ <- if (deficiency > 0L) 0 else {
      ch <- tryCatch(chol(Q), error = function(e) NULL)
      if (is.null(ch)) -Inf else 2 * sum(log(diag(ch)))
    }
    list(tau = tau, rho = rho, Q = Q, ldQ = ldQ)
  }
  cells <- list()
  if (is_ar1) { for (rho in rho_grid) for (tau in tau_grid)
                  cells[[length(cells) + 1L]] <- make_cell(c(tau, rho)) }
  else        { for (tau in tau_grid)
                  cells[[length(cells) + 1L]] <- make_cell(tau) }

  scat1 <- function(grad_eta) {
    g <- numeric(n_t)
    for (s in seq_along(map)) g[map[s]] <- g[map[s]] + grad_eta[s]
    g
  }
  list(
    n_field = n_t, n_sp = n_t, cells = cells, axes = NULL, make_cell = make_cell,
    valid = function(cell) is.finite(cell$ldQ) &&
                           (!is_ar1 || (cell$rho > -1 && cell$rho < 1)),
    offset = function(fp, cell) fp[map],
    scatter = scat1,
    prior_logp = function(fp, cell) {
      quad <- as.numeric(t(fp) %*% cell$Q %*% fp)
      0.5 * rank_eff * log(cell$tau) +
        (if (deficiency > 0L) 0 else 0.5 * cell$ldQ) - 0.5 * cell$tau * quad
    },
    prior_grad = function(fp, cell) cell$tau * as.numeric(cell$Q %*% fp),
    center = function(fp) if (deficiency > 0L) fp - mean(fp) else fp,
    constrain = constrain,
    to_phi = function(fp, cell) fp,
    type = type,
    to_hyper = function(cell) if (is_ar1) c(tau = cell$tau, rho = cell$rho)
                              else c(tau = cell$tau)
  )
}

# Resolve a temporal() term carried on a count family's arm into a temporal field
# spec for the areal-BFGS driver (#78). `temporal$time_idx` is the per-site time
# index (one entry per data row / site); it becomes the block's `map`.
.tobs_temporal_field_spec <- function(temporal, n_sites, family) {
  if (!inherits(temporal, "tobs_temporal"))
    stop("`temporal` must be a tobs_temporal object.", call. = FALSE)
  if (!is.null(temporal$group_idx))
    stop(sprintf(paste0("%s() temporal + areal spatial supports a single temporal ",
                        "field; grouped temporal( , group = ) is not wired on this ",
                        "path. (tulpaObs#78)"), family), call. = FALSE)
  ti <- as.integer(temporal$time_idx)
  if (length(ti) != n_sites)
    stop(sprintf(paste0("temporal term has %d time indices but the model has %d ",
                        "sites; one time index per site is required for %s."),
                 length(ti), n_sites, family), call. = FALSE)
  n_t <- if (!is.null(temporal$n_times)) as.integer(temporal$n_times)
         else max(ti, na.rm = TRUE)
  .tobs_temporal_field(temporal$type, map = ti, n_t = n_t,
                       cyclic = isTRUE(temporal$cyclic))
}

# Areal-BFGS nested-Laplace fit over one OR several latent field blocks (#78).
.tobs_areal_bfgs_fit <- function(eval, n_fixed, field, theta0_fix,
                                 max_iter = 300L, tol = 1e-8, label = "areal-bfgs",
                                 integration = c("auto", "ccd", "grid")) {
  integration <- match.arg(integration)

  blocks <- if (!is.null(field$n_field)) list(field) else field
  n_blk  <- length(blocks)
  nf_vec <- vapply(blocks, function(b) as.integer(b$n_field), 0L)
  nf     <- sum(nf_vec)
  blk_off <- cumsum(c(0L, nf_vec))            # field-slice starts within fi
  fi     <- n_fixed + seq_len(nf)
  # Per-block field index slice (into the global theta vector).
  blk_fi <- lapply(seq_len(n_blk),
                   function(b) n_fixed + blk_off[b] + seq_len(nf_vec[b]))
  is_bym2 <- vapply(blocks, function(b) isTRUE(b$bym2), TRUE)

  # One Laplace fit at a fixed field-hyperparameter cell (one per-block cell each).
  # Returns the joint marginal of (fixed params, field) plus the fixed-parameter
  # mode / cov block / per-block reported fields, or NULL on an invalid /
  # numerically failed cell.
  solve_cell <- function(cells_b) {
    if (!all(vapply(seq_len(n_blk),
                    function(b) blocks[[b]]$valid(cells_b[[b]]), TRUE)))
      return(NULL)
    block_offset <- function(theta) {
      off <- numeric(0L)
      for (b in seq_len(n_blk)) {
        fpb <- theta[blk_fi[[b]]]
        ob  <- blocks[[b]]$offset(fpb, cells_b[[b]])
        off <- if (length(off)) off + ob else ob
      }
      off
    }
    grad_wp <- function(theta) {              # grad of (log_lik + log_prior) over (fixed, field)
      e  <- eval(theta[seq_len(n_fixed)], block_offset(theta))
      gf <- numeric(nf)
      for (b in seq_len(n_blk)) {
        fpb <- theta[blk_fi[[b]]]
        sc  <- if (is_bym2[b]) blocks[[b]]$scatter(e$grad_eta, cells_b[[b]])
               else blocks[[b]]$scatter(e$grad_eta)
        gf[blk_off[b] + seq_len(nf_vec[b])] <- sc - blocks[[b]]$prior_grad(fpb, cells_b[[b]])
      }
      c(e$grad_fixed, gf)
    }
    ll_fn <- function(theta) {
      e  <- eval(theta[seq_len(n_fixed)], block_offset(theta))
      lp <- 0
      for (b in seq_len(n_blk))
        lp <- lp + blocks[[b]]$prior_logp(theta[blk_fi[[b]]], cells_b[[b]])
      e$log_lik + lp
    }
    th0 <- c(theta0_fix, numeric(nf))
    opt <- tryCatch(stats::optim(th0, function(t) -ll_fn(t), function(t) -grad_wp(t),
                   method = "BFGS", control = list(maxit = as.integer(max_iter), reltol = tol)),
                   error = function(e) NULL)
    if (is.null(opt)) return(NULL)
    th <- opt$par
    for (b in seq_len(n_blk))
      th[blk_fi[[b]]] <- blocks[[b]]$center(th[blk_fi[[b]]])
    H <- tryCatch(-.fp_fd_jacobian(grad_wp, th), error = function(e) NULL)
    if (is.null(H)) return(NULL)
    H <- 0.5 * (H + t(H)); ridge <- max(1e-8 * mean(abs(diag(H))), 1e-10); diag(H) <- diag(H) + ridge
    ch <- tryCatch(chol(H), error = function(e) NULL); if (is.null(ch)) return(NULL)
    ldH <- 2 * sum(log(diag(ch)))
    Hc <- H
    cc <- unlist(lapply(seq_len(n_blk),
                        function(b) blk_off[b] + which(blocks[[b]]$constrain)))
    if (length(cc)) {                     # sum-to-zero penalty on the constrained field block(s)
      pen <- 1e6 * mean(abs(diag(H))); idx <- n_fixed + cc
      Hc[idx, idx] <- Hc[idx, idx] + pen
    }
    cov_full <- tryCatch(solve(Hc), error = function(e) NULL)
    if (is.null(cov_full)) return(NULL)
    phi_b   <- lapply(seq_len(n_blk),
                      function(b) blocks[[b]]$to_phi(th[blk_fi[[b]]], cells_b[[b]]))
    hyper_b <- lapply(seq_len(n_blk),
                      function(b) if (is.function(blocks[[b]]$to_hyper))
                                    blocks[[b]]$to_hyper(cells_b[[b]]) else NULL)
    list(logm = ll_fn(th) - 0.5 * ldH,
         mode = th[seq_len(n_fixed)],
         phi  = phi_b, hyper = hyper_b,
         cov  = cov_full[seq_len(n_fixed), seq_len(n_fixed), drop = FALSE])
  }

  # Weighted (marginalised) summaries from a set of evaluated cells + weights.
  summarise <- function(res, w, method, pareto_k = NA_real_) {
    ok  <- vapply(res, Negate(is.null), TRUE) & is.finite(w) & w > 0
    if (!any(ok)) return(list(ok = FALSE))
    w[!ok] <- 0; w <- w / sum(w)
    wk    <- w[ok]; ik <- which(ok)
    modes <- t(vapply(ik, function(k) res[[k]]$mode, numeric(n_fixed)))
    beta_mean  <- as.numeric(crossprod(wk, modes))
    V <- matrix(0, n_fixed, n_fixed)
    for (j in seq_along(wk)) {
      dk <- modes[j, ] - beta_mean
      V <- V + wk[j] * (res[[ik[j]]]$cov + outer(dk, dk))
    }
    logm <- vapply(ik, function(k) res[[k]]$logm, numeric(1))
    # Posterior-mean field + field hyperparameters, per block (tau / rho or
    # sigma / rho), for reporting and for fixing the field precision on the NUTS
    # path (gcol33/tulpaObs#72). Block 1 is the spatial field (kept on the legacy
    # scalar slots `field_mean` / `hyper`); a temporal block 2 is reported under
    # `temporal_field` / `temporal_hyper`.
    n_sp_b <- vapply(blocks, function(b) as.integer(b$n_sp), 0L)
    field_means <- lapply(seq_len(n_blk), function(b) {
      phis <- t(vapply(ik, function(k) res[[k]]$phi[[b]], numeric(n_sp_b[b])))
      as.numeric(crossprod(wk, phis))
    })
    hyper_means <- lapply(seq_len(n_blk), function(b) {
      h1 <- res[[ik[1L]]]$hyper[[b]]
      if (is.null(h1)) return(NULL)
      Hm <- t(vapply(ik, function(k) res[[k]]$hyper[[b]], numeric(length(h1))))
      hm <- as.numeric(crossprod(wk, Hm)); names(hm) <- names(h1); hm
    })
    out <- list(ok = TRUE, beta_mean = beta_mean,
                field_mean = field_means[[1L]], hyper = hyper_means[[1L]],
                vcov = V, log_lik = sum(wk * logm),
                integration = method, pareto_k = pareto_k)
    if (n_blk >= 2L) {
      out$temporal_field <- field_means[[2L]]
      out$temporal_hyper <- hyper_means[[2L]]
    }
    out
  }

  # ---- outer integration: opt-in mode-centred CCD (single-block only,
  # gcol33/tulpaObs#60), silently declining to the fixed tensor grid when the
  # outer curvature is ill-conditioned. A multi-block fit uses the product grid.
  if (n_blk == 1L && identical(integration, "ccd") && !is.null(field$axes)) {
    eval_logm <- function(theta_phys) {
      r <- solve_cell(list(field$make_cell(theta_phys)))
      if (is.null(r) || !is.finite(r$logm)) NA_real_ else r$logm
    }
    cc <- tryCatch(.tobs_ccd_outer_grid(eval_logm, field$axes),
                   error = function(e) NULL)
    if (!is.null(cc)) {
      nn  <- nrow(cc$nodes)
      prog <- tulpa:::.tulpa_iter_progress(label, nn, unit = "cells")
      res <- vector("list", nn); logm <- rep(-Inf, nn)
      for (k in seq_len(nn)) {
        r <- solve_cell(list(field$make_cell(cc$nodes[k, ])))
        if (!is.null(r) && is.finite(r$logm)) { res[[k]] <- r; logm[k] <- r$logm }
        prog$tick()
      }
      prog$finish()
      if (any(is.finite(logm))) {
        w <- cc$dnode * exp(logm - max(logm[is.finite(logm)]))
        out <- summarise(res, w, "ccd", cc$pareto_k)
        if (isTRUE(out$ok)) return(out)
      }
    }
  }

  # Product grid over the blocks' per-cell grids.
  cell_grid <- .tobs_block_cell_product(blocks)
  n_grid <- length(cell_grid)
  prog <- tulpa:::.tulpa_iter_progress(label, n_grid, unit = "cells")
  res <- vector("list", n_grid); logm <- rep(-Inf, n_grid)
  for (k in seq_len(n_grid)) {
    r <- solve_cell(cell_grid[[k]])
    if (!is.null(r) && is.finite(r$logm)) { res[[k]] <- r; logm[k] <- r$logm }
    prog$tick()
  }
  prog$finish()
  if (!any(is.finite(logm))) return(list(ok = FALSE))
  w <- tulpa:::.nl_normalise_weights_safe(logm, "tau_grid / data")
  summarise(res, w, "grid")
}

# Cartesian product of each block's `cells` list -> a list of per-block cell
# tuples (each tuple a length-`n_blk` list, one cell per block). A single block
# returns one cell per element (identical to iterating `field$cells`).
.tobs_block_cell_product <- function(blocks) {
  grids <- lapply(blocks, function(b) b$cells)
  n_per <- vapply(grids, length, 0L)
  idx <- expand.grid(lapply(n_per, seq_len), KEEP.OUT.ATTRS = FALSE)
  lapply(seq_len(nrow(idx)), function(r)
    lapply(seq_along(grids), function(b) grids[[b]][[idx[r, b]]]))
}

# Resolve the field spec for an areal-BFGS family from the spatial term.
.tobs_areal_field_spec <- function(spatial, n_sites, family, map) {
  if (spatial$type %in% c("spde", "gp", "multiscale_gp"))
    stop(sprintf(paste0("%s() areal spatial supports icar() / car_proper() / bym2() ",
                        "under method = \"nested_laplace\"; the '%s' field is not yet ",
                        "wired for %s. (tulpaObs#51)"), family, spatial$type, family),
         call. = FALSE)
  if (!spatial$type %in% c("icar", "car_proper", "bym2"))
    stop(sprintf("%s() areal spatial supports icar() / car_proper() / bym2(); got '%s'.",
                 family, spatial$type), call. = FALSE)
  if (spatial$n_units != n_sites)
    stop(sprintf(paste0("spatial term has %d units but the model has %d sites; one ",
                        "spatial unit per site is required for %s."),
                 spatial$n_units, n_sites, family), call. = FALSE)
  adj <- if (!is.null(spatial$graph)) as.matrix(spatial$graph) else
    stop(sprintf("%s() spatial term must carry an adjacency graph.", family), call. = FALSE)
  if (identical(spatial$type, "bym2")) {
    sf <- spatial$scale_factor %||% compute_bym2_scale(spatial$graph)
    .areal_field_bym2(adj, sf, map, spatial$n_units)
  } else {
    .areal_field_car(adj, if (identical(spatial$type, "icar")) "icar" else "car_proper",
                     map, spatial$n_units)
  }
}
