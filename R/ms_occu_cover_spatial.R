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
#' @param K Number of shared latent spatial factors (`1 <= K <= n_species`).
#'   `K = 1` is the Stage-1 single-field case (loading vector, field vector);
#'   `K > 1` draws `K` ICAR fields with lower-triangular, positive-diagonal
#'   loadings and returns the `S x K` loading matrix / `N x K` field matrix.
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
                                           K          = 1L,
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
  K <- as.integer(K)
  if (K < 1L || K > n_species) {
    stop("K must satisfy 1 <= K <= n_species.", call. = FALSE)
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

  # K shared latent factors W (N x K): each a unit-marginal-scale ICAR draw
  # (Sorbye-Rue), the same convention simulate_occu_cover uses for its
  # single-species field. K = 1 reproduces the single-field random stream.
  Q       <- .occu_cover_icar_Q(adj)
  scale_q <- .occu_cover_icar_scale(adj)
  eig  <- eigen(Q, symmetric = TRUE)
  keep <- eig$values > 1e-8
  W <- matrix(0, N, K)
  for (k in seq_len(K)) {
    z_white <- stats::rnorm(sum(keep))
    wk <- as.numeric(eig$vectors[, keep, drop = FALSE] %*%
                       (z_white / sqrt(eig$values[keep])))
    wk <- wk - mean(wk)
    W[, k] <- wk / sqrt(scale_q)
  }

  # Per-species loadings L (S x K), lower-triangular with positive diagonal --
  # the rotation/sign/ordering-identified canonical form (gllvm/HMSC): factor k
  # loads on species k..S only, and L[k, k] > 0. K = 1 is the length-S loading
  # vector with the reference-species (sp1) sign anchor.
  L <- matrix(0, n_species, K)
  for (k in seq_len(K)) {
    rows <- k:n_species
    L[rows, k] <- stats::rnorm(length(rows), mean_load, sd_load)
  }
  b_occ <- matrix(stats::rnorm(n_species * p_occ, 0, rep(sd_occ, each = n_species)),
                  n_species, p_occ)
  b_p   <- matrix(stats::rnorm(n_species * p_p, 0, rep(sd_p, each = n_species)),
                  n_species, p_p)
  b_pos <- matrix(stats::rnorm(n_species * p_pos, 0, rep(sd_pos, each = n_species)),
                  n_species, p_pos)

  # Canonical sign anchor: make each factor's diagonal loading positive, flipping
  # that factor's (W[, k], L[, k]) pair together.
  for (k in seq_len(K)) {
    if (L[k, k] < 0) { L[, k] <- -L[, k]; W[, k] <- -W[, k] }
  }

  species_names <- paste0("sp", seq_len(n_species))
  y     <- array(NA_integer_, dim = c(N, J, n_species),
                 dimnames = list(NULL, NULL, species_names))
  y_pos <- array(NA_real_,    dim = c(N, J, n_species),
                 dimnames = list(NULL, NULL, species_names))
  psi   <- matrix(NA_real_, N, n_species, dimnames = list(NULL, species_names))
  zmat  <- matrix(NA_integer_, N, n_species, dimnames = list(NULL, species_names))

  for (s in seq_len(n_species)) {
    eta_psi <- as.vector(X_occ %*% (mu_occ + b_occ[s, ])) +
               as.numeric(W %*% L[s, ])
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
      K = K, positive = positive,
      mu_occ = mu_occ, mu_p = mu_p, mu_pos = mu_pos,
      sd_occ = sd_occ, sd_p = sd_p, sd_pos = sd_pos,
      sigma_pos = sigma_pos,
      b_occ = b_occ, b_p = b_p, b_pos = b_pos,
      # K = 1 keeps the Stage-1 shapes (loading vector, field vector); K > 1
      # returns the S x K loading matrix and the N x K field matrix.
      L = if (K == 1L) L[, 1L] else L,
      w = if (K == 1L) W[, 1L] else W,
      psi = psi, z = zmat
    )
  )
}


# ---------------------------------------------------------------------------
# Model binding
# ---------------------------------------------------------------------------

# Bind a K-factor spatial community occu_cover model. Wraps
# .tobs_build_ms_occu_cover for the (community) designs + per-species masked
# detection/cover matrices, and attaches the cell-graph ICAR structure (Q +
# Sorbye-Rue scale) shared by all K factors. `adj` is the N x N cell adjacency;
# N must equal the number of occupancy cells (rows of `data`). `K` is the number
# of shared latent spatial factors (1 <= K <= n_species).
.tobs_build_ms_occu_cover_spatial <- function(occ_formula, det_formula,
                                              pos_formula, data, y, y_pos,
                                              positive, species, adj, K = 1L,
                                              det_visit_formula = NULL,
                                              det_visit_data    = NULL,
                                              pos_visit_formula = NULL,
                                              pos_visit_data    = NULL) {
  if (identical(positive, "beta")) {
    stop("Stage 1 of the spatial community occu_cover supports ",
         "positive = 'lognormal' only.", call. = FALSE)
  }
  model <- .tobs_build_ms_occu_cover(
    occ_formula, det_formula, pos_formula, data, y, y_pos, positive, species,
    det_visit_formula = det_visit_formula, det_visit_data = det_visit_data,
    pos_visit_formula = pos_visit_formula, pos_visit_data = pos_visit_data)
  N <- model$n_sites
  if (is.null(adj) || !is.matrix(adj) || nrow(adj) != N || ncol(adj) != N) {
    stop(sprintf("adj must be an N x N adjacency with N = %d cells.", N),
         call. = FALSE)
  }
  K <- as.integer(K)
  if (K < 1L || K > model$n_species) {
    stop("K must satisfy 1 <= K <= n_species.", call. = FALSE)
  }
  model$model_type <- "ms_occu_cover_spatial"
  model$K          <- K
  model$adj        <- adj
  model$icar_Q     <- .occu_cover_icar_Q(adj)
  model$icar_scale <- .occu_cover_icar_scale(adj)
  model
}


# ---------------------------------------------------------------------------
# Penalised joint log-likelihood + gradient (inner mode-find target)
# ---------------------------------------------------------------------------

# Packed inner-latent parameter for the K-factor spatial community model, at
# fixed community covariances `Sigma`, loading prior SD `sd_L`, and per-factor
# field precisions `tau_w` (length K):
#   par = c(mu[P], b[S*P], vec(L)[S*K], vec(W)[N*K], log_disp)
# where P = P_occ + P_p + P_pos, b is species-major (species s occupies
# b[(s-1)*P + 1:P]), L (S x K, column-major) are the per-species loadings on the
# K shared fields W (N x K, column-major). K = 1 recovers the Stage-1 layout
# (L a length-S vector, W a length-N vector) bit for bit.
.ms_ocs_dims <- function(model) {
  pil <- model$process_info
  P_occ <- pil[[1L]]$p; P_p <- pil[[2L]]$p; P_pos <- pil[[3L]]$p
  list(P_occ = P_occ, P_p = P_p, P_pos = P_pos, P = P_occ + P_p + P_pos,
       S = model$n_species, N = model$n_sites, K = model$K %||% 1L,
       occ_idx = seq_len(P_occ), p_idx = P_occ + seq_len(P_p),
       pos_idx = P_occ + P_p + seq_len(P_pos))
}

.ms_ocs_unpack <- function(par, d) {
  off <- 0L
  mu  <- par[off + seq_len(d$P)]; off <- off + d$P
  b   <- lapply(seq_len(d$S), function(s) par[off + (s - 1L) * d$P + seq_len(d$P)])
  off <- off + d$S * d$P
  L   <- matrix(par[off + seq_len(d$S * d$K)], d$S, d$K); off <- off + d$S * d$K
  W   <- matrix(par[off + seq_len(d$N * d$K)], d$N, d$K); off <- off + d$N * d$K
  ld  <- par[off + 1L]
  list(mu = mu, b = b, L = L, W = W, ld = ld)
}

# Penalised joint log-likelihood and its gradient w.r.t. par. `Sinv` is the
# block-diagonal inverse community covariance over (occ, p, pos); `Pmu` the weak
# Gaussian precision on the community means; `inv_sdL2 = 1/sd_L^2` the loading
# ridge; `tau_w` the per-factor field precisions (scalar recycled, or length K;
# prior sum_k 0.5 * tau_w[k] * W[, k]' Q W[, k]). The shared-factor offset on the
# occupancy predictor is sum_k L[s, k] * W[, k] = W %*% L[s, ].
.ms_ocs_penll_grad <- function(model, par, Sinv, Pmu, inv_sdL2, tau_w,
                               grad = TRUE) {
  d <- .ms_ocs_dims(model); up <- .ms_ocs_unpack(par, d)
  mu <- up$mu; b <- up$b; L <- up$L; W <- up$W; ld <- up$ld
  Q  <- model$icar_Q
  tau_w <- rep_len(tau_w, d$K)
  cl <- function(e) pmin(pmax(e, -30), 30)

  ll <- 0
  g_mu <- numeric(d$P); g_b <- vector("list", d$S)
  g_L  <- matrix(0, d$S, d$K); g_W <- matrix(0, d$N, d$K); g_ld <- 0
  for (s in seq_len(d$S)) {
    v   <- .ms_occu_cover_species_view(model, s)
    th  <- mu + b[[s]]
    eta <- .occu_cover_eta_from_par(v, th[d$occ_idx], th[d$p_idx], th[d$pos_idx])
    # Inject the shared-factor offset on the occupancy predictor.
    eta$psi <- stats::plogis(cl(as.numeric(v$X_occ %*% th[d$occ_idx]) +
                                  as.numeric(W %*% L[s, ])))
    ll <- ll + sum(.occu_cover_site_ll(v, eta$psi, eta$p_mat, eta$ep_mat, ld))
    if (grad) {
      eg <- .occu_cover_eta_grad(v, eta$psi, eta$p_mat, eta$ep_mat, ld)
      cg <- .occu_cover_coef_grad(v, eg)              # c(g_occ, g_det, g_pos, g_ld)
      cg_coef <- cg[seq_len(d$P)]
      g_mu    <- g_mu + cg_coef
      g_b[[s]] <- cg_coef                              # RE prior added below
      g_L[s, ] <- as.numeric(crossprod(W, eg$g_psi))  # t(W) %*% g_psi_s
      g_W      <- g_W + outer(eg$g_psi, L[s, ])
      g_ld    <- g_ld + cg[d$P + 1L]
    }
  }

  # ---- priors (penalty) ----
  ll <- ll - 0.5 * sum(vapply(b, function(bb) as.numeric(bb %*% Sinv %*% bb), 0))
  ll <- ll - 0.5 * as.numeric(mu %*% Pmu %*% mu)
  ll <- ll - 0.5 * inv_sdL2 * sum(L^2)
  QW <- Q %*% W                                        # N x K
  ll <- ll - 0.5 * sum(tau_w * colSums(W * QW))

  if (!grad) return(list(ll = ll))

  g_mu <- g_mu - as.numeric(Pmu %*% mu)
  for (s in seq_len(d$S)) g_b[[s]] <- g_b[[s]] - as.numeric(Sinv %*% b[[s]])
  g_L  <- g_L - inv_sdL2 * L
  g_W  <- g_W - sweep(as.matrix(QW), 2L, tau_w, "*")

  grad_vec <- c(g_mu, unlist(g_b), as.numeric(g_L), as.numeric(g_W), g_ld)
  list(ll = ll, grad = grad_vec)
}


# ---------------------------------------------------------------------------
# Inner mode-find (conditional on the hyperparameters)
# ---------------------------------------------------------------------------

# Find the joint posterior mode of the inner latent
# par = c(mu, {b_s}, vec(L), vec(W), log_disp) at fixed community covariances
# `Sigma`, loading SD `sd_L`, and per-factor field precisions `tau_w`, by
# quasi-Newton ascent on the penalised joint log-likelihood with the analytic
# gradient (.ms_ocs_penll_grad). Correctness first: the arrowhead / Schur-folded
# Newton of the design doc is a later speed increment; here a BFGS solve on the
# verified gradient locates the same mode. The per-factor sign symmetry
# (L[, k], W[, k]) -> (-L[, k], -W[, k]) is resolved by making each factor's
# diagonal loading L[k, k] positive (matching the simulator's canonical form);
# for K > 1 the residual rotational freedom is resolved post hoc (psi is
# rotation-invariant, so the fit is unaffected). Returns the unpacked mode (w / L
# as a vector at K = 1, a matrix at K > 1) plus the achieved penalised
# log-likelihood.
.ms_ocs_inner_mode <- function(model, Sigma, sd_L = 1.0, tau_w = 1.0,
                               par_init = NULL, sigma.beta = 5,
                               maxit = 400L, hessian = FALSE) {
  d <- .ms_ocs_dims(model); K <- d$K
  Sinv <- matrix(0, d$P, d$P)
  Sinv[d$occ_idx, d$occ_idx] <- solve(Sigma$occ)
  Sinv[d$p_idx,   d$p_idx]   <- solve(Sigma$p)
  Sinv[d$pos_idx, d$pos_idx] <- solve(Sigma$pos)
  Pmu      <- diag(1 / sigma.beta^2, d$P)
  inv_sdL2 <- 1 / sd_L^2

  if (is.null(par_init)) {
    views <- lapply(seq_len(d$S), function(s) .ms_occu_cover_species_view(model, s))
    any_det <- mean(vapply(views, function(v) mean(rowSums(v$y * v$valid) > 0),
                           numeric(1)))
    pos_vals <- unlist(lapply(views, function(v) v$y_pos[v$valid & v$y == 1L]))
    mu0 <- numeric(d$P)
    mu0[d$occ_idx][1L] <- stats::qlogis(min(max(any_det, 1e-3), 1 - 1e-3))
    ld0 <- log(0.4)
    if (length(pos_vals)) {
      mu0[d$pos_idx][1L] <- mean(log(pos_vals))
      ld0 <- log(stats::sd(log(pos_vals)) + 0.1)
    }
    # The bilinear sum_k L_sk w_kc term has a saddle at (L, W) = (0, 0): there
    # both gradients vanish, so a zero start never leaves it. Warm-start from the
    # leading K EOFs of the per-species empirical-occupancy field (the design-doc
    # recipe) to land in the right basin. D[c, s] = 1 if species s was ever
    # detected at cell c; centring per species removes prevalence so the top
    # singular vectors capture the shared spatial gradients + their loadings.
    D <- vapply(views, function(v) as.numeric(rowSums(v$y * v$valid) > 0),
                numeric(d$N))                          # N x S
    Dc <- sweep(D, 2L, colMeans(D))
    W0 <- matrix(0, d$N, K); L0 <- matrix(0, d$S, K)
    sv <- tryCatch(svd(Dc, nu = K, nv = K), error = function(e) NULL)
    if (!is.null(sv)) {
      for (k in seq_len(K)) {
        if (sv$d[k] <= 1e-8) next
        u <- as.numeric(sv$u[, k]); u <- u - mean(u)
        sdu <- stats::sd(u)
        if (sdu > 1e-8) {
          W0[, k] <- u / sdu                           # unit-scale field start
          L0[, k] <- as.numeric(sv$v[, k]) * sv$d[k] * sdu / sqrt(d$N)
        }
      }
    }
    par_init <- c(mu0, numeric(d$S * d$P), as.numeric(L0), as.numeric(W0), ld0)
  }

  fn <- function(p) -.ms_ocs_penll_grad(model, p, Sinv, Pmu, inv_sdL2, tau_w,
                                        grad = FALSE)$ll
  gr <- function(p) -.ms_ocs_penll_grad(model, p, Sinv, Pmu, inv_sdL2, tau_w,
                                        grad = TRUE)$grad
  opt <- stats::optim(par_init, fn, gr, method = "BFGS",
                      control = list(maxit = maxit, reltol = 1e-10),
                      hessian = hessian)
  up <- .ms_ocs_unpack(opt$par, d)
  # Canonical sign anchor: make each factor's diagonal loading L[k, k] positive,
  # flipping (L[, k], W[, k]) together. The Hessian over the packed par is
  # sign-flip invariant in those blocks, so it stays valid for the M-step
  # covariances after the anchor.
  for (k in seq_len(K)) {
    if (up$L[k, k] < 0) { up$L[, k] <- -up$L[, k]; up$W[, k] <- -up$W[, k] }
  }
  list(mu = up$mu, b = up$b, ld = up$ld,
       L = if (K == 1L) as.numeric(up$L) else up$L,
       w = if (K == 1L) as.numeric(up$W) else up$W,
       par = opt$par, logpen = -opt$value, convergence = opt$convergence,
       hessian = if (hessian) opt$hessian else NULL, d = d)
}


# ---------------------------------------------------------------------------
# Laplace-EM fitter (Stage 1)
# ---------------------------------------------------------------------------

# Fit the K-factor spatial community occu_cover model by Laplace-EM: the inner
# mode-find above is the E-step; the M-step is the closed-form community
# covariance update per arm (Sigma_arm = mean_s[b_s b_s' + Cov(b_s)], the
# posterior second moment from the mode + the Hessian-block covariance, Louis
# 1982) and the per-factor rank-deficient ICAR field-precision update
# tau_w[k] = (N - 1) / (W[, k]' Q W[, k] + tr(Q Cov_wk)). The loading SD `sd_L`
# fixes the amplitude split and is held at its supplied value (a hyperparameter).
.tobs_fit_ms_occu_cover_spatial <- function(model, sd_L = 1.0,
                                            max.em = 30L, tol = 1e-3,
                                            sigma.beta = 5, verbose = FALSE) {
  d <- .ms_ocs_dims(model); K <- d$K
  Q <- model$icar_Q
  rank_w <- d$N - 1L                       # ICAR rank (one null direction)

  Sigma <- list(occ = diag(0.3^2, d$P_occ), p = diag(0.3^2, d$P_p),
                pos = diag(0.3^2, d$P_pos))
  tau_w <- rep(1.0, K)
  par_init <- NULL
  prev <- -Inf

  # Packed offsets of each species' b block and the K shared field blocks.
  b_off <- function(s) d$P + (s - 1L) * d$P
  w_off <- d$P + d$S * d$P + d$S * K

  ridge_inv <- function(H) {
    Hs <- (H + t(H)) / 2
    for (eps in c(0, 1e-8, 1e-6, 1e-4, 1e-2)) {
      ch <- tryCatch(chol(Hs + diag(eps, nrow(Hs))), error = function(e) NULL)
      if (!is.null(ch)) return(chol2inv(ch))
    }
    # SVD pseudo-inverse fallback (no MASS dependency).
    sv <- svd(Hs)
    pos <- sv$d > max(sv$d) * 1e-10
    sv$v[, pos, drop = FALSE] %*% ((1 / sv$d[pos]) *
        t(sv$u[, pos, drop = FALSE]))
  }

  for (em in seq_len(max.em)) {
    fit <- .ms_ocs_inner_mode(model, Sigma, sd_L = sd_L, tau_w = tau_w,
                              par_init = par_init, sigma.beta = sigma.beta,
                              hessian = TRUE)
    par_init <- fit$par
    Cov <- ridge_inv(fit$hessian)          # posterior covariance at the mode

    # M-step: community covariances (second moment = mode outer product + block
    # posterior covariance), per arm, averaged over species.
    acc <- list(occ = matrix(0, d$P_occ, d$P_occ),
                p   = matrix(0, d$P_p,   d$P_p),
                pos = matrix(0, d$P_pos, d$P_pos))
    for (s in seq_len(d$S)) {
      bs  <- fit$b[[s]]
      off <- b_off(s)
      for (arm in c("occ", "p", "pos")) {
        ai  <- d[[paste0(arm, "_idx")]]
        idx <- off + ai
        bb  <- bs[ai]
        acc[[arm]] <- acc[[arm]] + outer(bb, bb) + Cov[idx, idx, drop = FALSE]
      }
    }
    Sigma <- lapply(acc, function(A) {
      A <- A / d$S
      (A + t(A)) / 2
    })
    names(Sigma) <- c("occ", "p", "pos")

    # M-step: per-factor ICAR field precision (rank-deficient GMRF update).
    Wm <- matrix(fit$w, d$N, K)
    for (k in seq_len(K)) {
      widx   <- w_off + (k - 1L) * d$N + seq_len(d$N)
      Cov_wk <- Cov[widx, widx, drop = FALSE]
      quad   <- as.numeric(Wm[, k] %*% (Q %*% Wm[, k])) + sum(Q * Cov_wk)
      tau_w[k] <- rank_w / max(quad, 1e-8)
    }

    if (verbose) {
      cat(sprintf("EM %2d  logpen=%.3f  tau_w[1]=%.3f  sd_occ1=%.3f\n",
                  em, fit$logpen, tau_w[1L], sqrt(Sigma$occ[1L, 1L])))
    }
    if (is.finite(prev) && abs(fit$logpen - prev) < tol * (abs(prev) + tol)) {
      prev <- fit$logpen; break
    }
    prev <- fit$logpen
  }

  c(fit, list(Sigma = Sigma, tau_w = tau_w, sd_L = sd_L, em_logpen = prev,
              cov = Cov))
}


# ---------------------------------------------------------------------------
# Front-door wrapper (tobs_fit) + spatial-term detector
# ---------------------------------------------------------------------------

# Wrap the Laplace-EM output (.tobs_fit_ms_occu_cover_spatial) into a tobs_fit,
# mirroring build_ms_occu_cover_fit (the non-spatial community wrapper) and
# adding the shared-factor block (field w, per-species loadings L, field
# precision tau_w). The community-mean + dispersion posterior covariance comes
# from the packed-par Hessian block (par = c(mu[P], b[S*P], L[S], w[N], ld)).
build_ms_occu_cover_spatial_fit <- function(model, fit) {
  d   <- fit$d
  pil <- model$process_info
  P   <- d$P

  beta_names <- c(
    paste0("psi_", pil[[1L]]$coef_names),
    paste0("p_",   pil[[2L]]$coef_names),
    paste0("pos_", pil[[3L]]$coef_names)
  )
  disp_name <- "log_sigma_pos"
  par_names <- c(beta_names, disp_name)

  mu <- fit$mu; ld <- fit$ld
  means <- c(mu, ld); names(means) <- par_names

  Cov  <- fit$cov
  npar <- length(fit$par)
  sel  <- c(seq_len(P), npar)                 # community means + log-dispersion
  V <- Cov[sel, sel, drop = FALSE]; V <- (V + t(V)) / 2
  dimnames(V) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  # Per-species community structure (mu + BLUP deviations) per arm.
  B <- do.call(rbind, fit$b)                  # S x P
  arm_idx <- list(occ = d$occ_idx, p = d$p_idx, pos = d$pos_idx)
  arm_block <- function(arm) {
    idx  <- arm_idx[[arm]]
    blup <- B[, idx, drop = FALSE]
    coef <- sweep(blup, 2L, mu[idx], "+")
    rownames(blup) <- rownames(coef) <- model$species_names
    list(blup = blup, coef = coef)
  }
  occ_b <- arm_block("occ"); p_b <- arm_block("p"); pos_b <- arm_block("pos")
  colnames(occ_b$blup) <- colnames(occ_b$coef) <- pil[[1L]]$coef_names
  colnames(p_b$blup)   <- colnames(p_b$coef)   <- pil[[2L]]$coef_names
  colnames(pos_b$blup) <- colnames(pos_b$coef) <- pil[[3L]]$coef_names

  Sigma_occ <- fit$Sigma$occ; Sigma_p <- fit$Sigma$p; Sigma_pos <- fit$Sigma$pos
  dimnames(Sigma_occ) <- list(pil[[1L]]$coef_names, pil[[1L]]$coef_names)
  dimnames(Sigma_p)   <- list(pil[[2L]]$coef_names, pil[[2L]]$coef_names)
  dimnames(Sigma_pos) <- list(pil[[3L]]$coef_names, pil[[3L]]$coef_names)

  # Shared-factor block: field posterior mean + per-cell marginal SD (from the
  # K field blocks of the joint posterior covariance) and the per-species
  # loadings. K = 1 keeps the Stage-1 vector shapes; K > 1 returns N x K / S x K.
  K      <- d$K
  w_off  <- P + d$S * P + d$S * K
  widx   <- w_off + seq_len(d$N * K)
  field_sd <- matrix(sqrt(pmax(diag(Cov)[widx], 0)), d$N, K)
  L <- fit$L
  if (K == 1L) {
    field_sd <- as.numeric(field_sd)
    names(L) <- model$species_names
  } else {
    rownames(L) <- model$species_names
    colnames(L) <- paste0("factor", seq_len(K))
    colnames(field_sd) <- paste0("factor", seq_len(K))
  }

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(fit$logpen, n_draws),
    log_lik      = fit$logpen,
    N            = sum(model$valid)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    n_fixed      = length(means),
    fixed_names  = par_names,
    process_info = pil,
    model        = model,
    method       = "laplace-em",
    positive     = model$positive,
    spatial      = list(
      type     = "icar",
      K        = K,
      field    = fit$w,
      field_sd = field_sd,
      loadings = L,
      tau_w    = fit$tau_w,
      sd_L     = fit$sd_L
    ),
    ms_community = list(
      Sigma_occ = Sigma_occ, Sigma_p = Sigma_p, Sigma_pos = Sigma_pos,
      sd_occ = sqrt(pmax(diag(Sigma_occ), 0)),
      sd_p   = sqrt(pmax(diag(Sigma_p),   0)),
      sd_pos = sqrt(pmax(diag(Sigma_pos), 0)),
      coef_occ = occ_b$coef, coef_p = p_b$coef, coef_pos = pos_b$coef,
      blup_occ = occ_b$blup, blup_p = p_b$blup, blup_pos = pos_b$blup
    ),
    convergence  = list(converged = identical(fit$convergence, 0L),
                        n_iter = NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}


# Posterior draws of the per-species per-cell occupancy probability psi for the
# K-factor spatial community fit. Samples the packed inner latent from its
# Gaussian Laplace posterior N(mode, Cov) -- via an eigen PSD square root, since
# the joint Hessian is only PSD in the confounded ICAR-level / intercept (and,
# for K > 1, rotational) directions -- and maps each draw through the occupancy
# predictor eta_s = X_occ (mu_occ + b_s_occ) + sum_k L_sk w_kc. psi depends on
# (L, W) only through the product W L[s, ]', so it is invariant to the per-factor
# sign and (K > 1) rotational symmetries and needs no anchoring. Returns an
# [n_draws x N x S] array.
.ms_ocs_psi_posterior <- function(model, fit, n_draws = 300L) {
  d <- fit$d
  Sg <- (fit$cov + t(fit$cov)) / 2
  eg <- eigen(Sg, symmetric = TRUE)
  rt <- eg$vectors %*% (sqrt(pmax(eg$values, 0)) * t(eg$vectors))   # PSD sqrt
  npar <- length(fit$par)
  Z <- matrix(stats::rnorm(n_draws * npar), n_draws, npar)
  draws_par <- sweep(Z %*% rt, 2L, fit$par, "+")

  X_occ <- lapply(seq_len(d$S),
                  function(s) .ms_occu_cover_species_view(model, s)$X_occ)
  cl  <- function(e) pmin(pmax(e, -30), 30)
  psi <- array(0, dim = c(n_draws, d$N, d$S))
  for (i in seq_len(n_draws)) {
    up <- .ms_ocs_unpack(draws_par[i, ], d)
    for (s in seq_len(d$S)) {
      th_occ <- (up$mu + up$b[[s]])[d$occ_idx]
      eta <- as.numeric(X_occ[[s]] %*% th_occ) + as.numeric(up$W %*% up$L[s, ])
      psi[i, , s] <- stats::plogis(cl(eta))
    }
  }
  psi
}


# Detect a Stage-1 spatial request on the three occu_cover arms. Returns the
# shared-field adjacency + the fixed-effects occupancy formula when the
# occupancy arm carries a single icar() term (and detection / cover are plain),
# NULL when no arm carries a structured term (the non-spatial path), and errors
# on any other structured term (the supported surface is icar() on occupancy).
.tobs_ms_ocs_spatial_request <- function(occ_formula, det_formula, pos_formula,
                                         data) {
  parse_terms <- function(f) {
    if (is.null(f)) return(list())
    .tobs_parse_formula(f, data = data)$terms
  }
  occ_terms <- parse_terms(occ_formula)
  det_terms <- parse_terms(det_formula)
  pos_terms <- parse_terms(pos_formula)

  if (length(det_terms) || length(pos_terms)) {
    stop("ms_occu_cover(): structured terms are supported on the occupancy arm ",
         "only (Stage 1: a single icar() shared field). Use a plain formula on ",
         "detection / cover.", call. = FALSE)
  }
  if (length(occ_terms) == 0L) return(NULL)             # non-spatial path
  if (length(occ_terms) > 1L) {
    stop("ms_occu_cover(): Stage 1 supports a single icar() term on the ",
         "occupancy arm.", call. = FALSE)
  }
  spec <- occ_terms[[1L]]
  if (!inherits(spec, "tobs_spatial") || !identical(spec$label, "icar")) {
    stop(sprintf("ms_occu_cover(): Stage 1 spatial supports icar() only; got %s().",
                 spec$label %||% class(spec)[1L]), call. = FALSE)
  }
  list(graph  = spec$graph,
       fe_occ = .tobs_parse_formula(occ_formula, data = data)$fe_formula)
}
