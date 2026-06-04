# ms_occu_cover_spatial.R - reduced-rank spatial-factor community occu_cover
# (gcol33/tulpa#67, Stage 1: K = 1).
#
# The community / joint-SDM form of occu_cover: one model over many species that
# shares a small set of latent spatial factors and borrows strength via
# community priors on the per-species coefficients, instead of a separate fit
# per species. Stage 1 is the K = 1 case -- a single shared ICAR field w on the
# cell graph, with per-species loadings L_s on the OCCUPANCY (state) predictor:
#
#   logit psi_{s,c} = X_occ_c . (mu_occ + b_occ_s) + L_s * w_c
#   logit p_{s,i,j} = X_p_{ij} . (mu_p  + b_p_s)
#   g(cover)_{s,i,j} = X_pos_{ij} . (mu_pos + b_pos_s)
#   b_.s ~ N(0, Sigma_.)   (community RE)   w ~ ICAR (unit marginal scale)
#
# A single shared field with only a species intercept can shift each species'
# map level but not its SHAPE, so it forces every species onto one map; the
# per-species loading L_s gives each species its own spatial shape as a scaling
# of the shared factor -- the reduced-rank (HMSC / spatial-gllvm) structure that
# borrows strength and lets rare taxa get a calibrated map. Stage 1 places the
# factor on the state predictor only (a cover-arm factor is a later stage).
#
# This file currently provides the Stage-1 SIMULATOR (ground truth for the
# parameter-recovery harness the fitter is built against). The community-Newton
# + Sigma M-step fitter extends .tobs_fit_ms_occu_cover with the shared w + L_s
# block and lands in subsequent increments.

# Identifiability (K = 1). The bilinear term L_s * w_c is invariant under
# (L_s -> c L_s, w -> w / c) for any c != 0, and under the joint sign flip
# (L -> -L, w -> -w). Stage 1 resolves the scale by drawing w at unit marginal
# scale (Sorbye-Rue), so L_s carries the amplitude, and resolves the sign with a
# reference-species anchor: the loading of the first species is made positive
# (negating w and every L_s together if needed). The fitter recovers the truth
# in this same canonical form.

#' Simulate a reduced-rank spatial-factor community occu_cover data set (K = 1)
#'
#' Generates occupancy / detection / cover data for `n_species` species sharing
#' one latent ICAR spatial factor `w` on the cell graph, with per-species
#' loadings `L_s` on the occupancy state predictor and Gaussian community priors
#' on the per-species arm coefficients. This is the ground-truth generator for
#' the Stage-1 reduced-rank spatial JSDM (gcol33/tulpa#67).
#'
#' @param adj N x N 0/1 adjacency matrix of the cell graph (required); `N` cells.
#' @param n_species Number of species.
#' @param J Number of detection visits per cell.
#' @param n_occ_covs,n_det_covs,n_pos_covs Number of (Gaussian) covariates on the
#'   occupancy, detection, and cover arms; each arm also has an intercept.
#' @param mu_occ,mu_p,mu_pos Community mean coefficient vectors (intercept first).
#'   `NULL` picks sensible defaults of the right length.
#' @param sd_occ,sd_p,sd_pos Community RE SDs (diagonal `Sigma_.`); scalar
#'   (recycled) or per-coefficient.
#' @param mean_load,sd_load Mean and SD of the per-species loadings `L_s`.
#' @param sigma_pos Lognormal cover residual SD (on the log scale).
#' @param positive Cover family; only `"lognormal"` in Stage 1.
#' @param seed Optional RNG seed.
#'
#' @return A list with `y` (N x J x n_species 0/1 detections), `y_pos`
#'   (N x J x n_species cover, `NA` off detected visits), `data` (cell-level
#'   covariate frame), `species`, and `truth` (the canonical-form generating
#'   parameters: `mu_*`, `sd_*`, per-species `b_*`, loadings `L`, field `w`,
#'   `psi`, `z`).
#' @export
simulate_ms_occu_cover_spatial <- function(adj,
                                           n_species  = 8L,
                                           J          = 4L,
                                           n_occ_covs = 1L,
                                           n_det_covs = 1L,
                                           n_pos_covs = 1L,
                                           mu_occ     = NULL,
                                           mu_p       = NULL,
                                           mu_pos     = NULL,
                                           sd_occ     = 0.4,
                                           sd_p       = 0.4,
                                           sd_pos     = 0.3,
                                           mean_load  = 0,
                                           sd_load    = 1.0,
                                           sigma_pos  = 0.4,
                                           positive   = c("lognormal", "beta"),
                                           seed       = NULL) {
  positive <- match.arg(positive)
  if (identical(positive, "beta")) {
    stop("Stage 1 of the spatial community occu_cover supports ",
         "positive = 'lognormal' only.", call. = FALSE)
  }
  if (missing(adj) || is.null(adj) || !is.matrix(adj) ||
      nrow(adj) != ncol(adj)) {
    stop("adj must be a square N x N adjacency matrix.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  N <- nrow(adj)

  p_occ <- n_occ_covs + 1L
  p_p   <- n_det_covs + 1L
  p_pos <- n_pos_covs + 1L
  if (is.null(mu_occ)) mu_occ <- c(0.0, rep(0.6, n_occ_covs))
  if (is.null(mu_p))   mu_p   <- c(0.2, rep(-0.4, n_det_covs))
  if (is.null(mu_pos)) mu_pos <- c(log(5), rep(0.3, n_pos_covs))
  stopifnot(length(mu_occ) == p_occ, length(mu_p) == p_p,
            length(mu_pos) == p_pos)
  rec <- function(sd, p) if (length(sd) == 1L) rep(sd, p) else sd
  sd_occ <- rec(sd_occ, p_occ); sd_p <- rec(sd_p, p_p); sd_pos <- rec(sd_pos, p_pos)

  # Cell-level (occupancy) covariate frame.
  occ_covs <- if (n_occ_covs > 0L) {
    m <- matrix(stats::rnorm(N * n_occ_covs), N, n_occ_covs)
    df <- as.data.frame(m); names(df) <- paste0("occ_cov", seq_len(n_occ_covs)); df
  } else data.frame(row.names = seq_len(N))
  X_occ <- if (ncol(occ_covs)) stats::model.matrix(~ ., occ_covs)
           else stats::model.matrix(~ 1, data.frame(row.names = seq_len(N)))

  # Cell-level detection / cover designs (community covariates, shared across
  # species and visits in Stage 1).
  det_covs <- if (n_det_covs > 0L) {
    m <- matrix(stats::rnorm(N * n_det_covs), N, n_det_covs)
    df <- as.data.frame(m); names(df) <- paste0("det_cov", seq_len(n_det_covs)); df
  } else data.frame(row.names = seq_len(N))
  pos_covs <- if (n_pos_covs > 0L) {
    m <- matrix(stats::rnorm(N * n_pos_covs), N, n_pos_covs)
    df <- as.data.frame(m); names(df) <- paste0("pos_cov", seq_len(n_pos_covs)); df
  } else data.frame(row.names = seq_len(N))
  X_p   <- if (ncol(det_covs)) stats::model.matrix(~ ., det_covs)
           else stats::model.matrix(~ 1, data.frame(row.names = seq_len(N)))
  X_pos <- if (ncol(pos_covs)) stats::model.matrix(~ ., pos_covs)
           else stats::model.matrix(~ 1, data.frame(row.names = seq_len(N)))

  # Shared latent factor w: unit-marginal-scale ICAR draw (Sorbye-Rue), the same
  # convention simulate_occu_cover uses for its single-species field.
  Q       <- .occu_cover_icar_Q(adj)
  scale_q <- .occu_cover_icar_scale(adj)
  eig  <- eigen(Q, symmetric = TRUE)
  keep <- eig$values > 1e-8
  z_white <- stats::rnorm(sum(keep))
  w <- as.numeric(eig$vectors[, keep, drop = FALSE] %*%
                    (z_white / sqrt(eig$values[keep])))
  w <- w - mean(w)
  w <- w / sqrt(scale_q)

  # Per-species loadings and community RE coefficients.
  L     <- stats::rnorm(n_species, mean_load, sd_load)
  b_occ <- matrix(stats::rnorm(n_species * p_occ, 0, rep(sd_occ, each = n_species)),
                  n_species, p_occ)
  b_p   <- matrix(stats::rnorm(n_species * p_p, 0, rep(sd_p, each = n_species)),
                  n_species, p_p)
  b_pos <- matrix(stats::rnorm(n_species * p_pos, 0, rep(sd_pos, each = n_species)),
                  n_species, p_pos)

  # Sign anchor: make the reference species' (sp1) loading positive, flipping
  # (w, L) together so the truth is in the fitter's canonical form.
  if (L[1L] < 0) { L <- -L; w <- -w }

  species_names <- paste0("sp", seq_len(n_species))
  y     <- array(NA_integer_, dim = c(N, J, n_species),
                 dimnames = list(NULL, NULL, species_names))
  y_pos <- array(NA_real_,    dim = c(N, J, n_species),
                 dimnames = list(NULL, NULL, species_names))
  psi   <- matrix(NA_real_, N, n_species, dimnames = list(NULL, species_names))
  zmat  <- matrix(NA_integer_, N, n_species, dimnames = list(NULL, species_names))

  for (s in seq_len(n_species)) {
    eta_psi <- as.vector(X_occ %*% (mu_occ + b_occ[s, ])) + L[s] * w
    ps      <- stats::plogis(eta_psi)
    z       <- stats::rbinom(N, 1L, ps)
    pp      <- stats::plogis(as.vector(X_p %*% (mu_p + b_p[s, ])))
    eta_pos <- as.vector(X_pos %*% (mu_pos + b_pos[s, ]))
    for (i in seq_len(N)) {
      det_ij <- stats::rbinom(J, 1L, z[i] * pp[i])
      y[i, , s] <- det_ij
      hit <- det_ij == 1L
      if (any(hit)) {
        y_pos[i, hit, s] <- exp(eta_pos[i] + stats::rnorm(sum(hit), 0, sigma_pos))
      }
    }
    psi[, s] <- ps; zmat[, s] <- z
  }

  data <- data.frame(row.names = seq_len(N))
  for (df in list(occ_covs, det_covs, pos_covs)) if (ncol(df)) data <- cbind(data, df)

  list(
    y = y, y_pos = y_pos, data = data, species = species_names, adj = adj,
    truth = list(
      K = 1L, positive = positive,
      mu_occ = mu_occ, mu_p = mu_p, mu_pos = mu_pos,
      sd_occ = sd_occ, sd_p = sd_p, sd_pos = sd_pos,
      sigma_pos = sigma_pos,
      b_occ = b_occ, b_p = b_p, b_pos = b_pos,
      L = L, w = w, psi = psi, z = zmat
    )
  )
}
