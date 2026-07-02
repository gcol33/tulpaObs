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
#   b_.s ~ N(0, Sigma_.)   (community RE)   w ~ GMRF field (unit marginal scale)
#
# A single shared field with only a species intercept can shift each species'
# map level but not its SHAPE, so it forces every species onto one map; the
# per-species loading L_s gives each species its own spatial shape as a scaling
# of the shared factor -- the reduced-rank (HMSC / spatial-gllvm) structure that
# borrows strength and lets rare taxa get a calibrated map. The factor sits on
# the occupancy (state) predictor by default; a cover-arm factor lets the SAME
# field also load on the cover predictor through a free loading matrix L_pos, so
# the latent spatial structure is shared across the two processes.
#
# The shared field is an areal GMRF whose structure is set by the front-door term:
# icar() (improper intrinsic CAR, the default), car_proper() (proper CAR with a
# per-factor correlation rho), or bym2() (the Riebler 2016 convolution with a
# variance fraction phi); the hyperparameter is estimated by EM. The fitter is a
# Laplace-EM over the joint mode of c(mu, {b_s}, L, W, log_disp) with closed-form
# community-cov and field (tau, hyperparameter) M-steps; rank selection and a
# varimax post-rotation sit on top. The simulator below is the ground truth for
# the parameter-recovery harness.

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
#' the Stage-1 reduced-rank spatial JSDM.
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
#' @param cover_factor Logical; when `TRUE` the same shared fields `W` also load
#'   on the cover (positive) predictor through a free `S x K` loading matrix
#'   `L_pos` (the cover-arm factor). The cover-factor
#'   draws are gated, so `FALSE` (the default) reproduces the no-factor RNG stream
#'   exactly. `truth$L_pos` carries the generating cover loadings.
#' @param mean_load_pos,sd_load_pos Mean and SD of the cover-arm loadings `L_pos`
#'   (used only when `cover_factor = TRUE`).
#' @param field Areal structure of the shared latent factors: `"icar"` (improper
#'   intrinsic CAR, the default), `"car_proper"` (proper CAR with a correlation
#'   `rho`, field precision `tau (D - rho A)`), or `"bym2"` (the Riebler 2016
#'   reparameterised convolution with a spatial-variance fraction `phi`). `"icar"`
#'   reproduces the Stage 1-2 RNG stream byte for byte; `truth$field` / `truth$rho`
#'   / `truth$phi` carry the choice.
#' @param rho Proper-CAR correlation in `[0, 1)`, used only when
#'   `field = "car_proper"` (the `rho -> 1` limit is the ICAR field).
#' @param phi BYM2 spatial-variance fraction in `[0, 1]`, used only when
#'   `field = "bym2"` (`phi = 1` is the pure ICAR field, `phi = 0` pure iid).
#' @param sigma_pos Lognormal cover residual SD (on the log scale).
#' @param response Cover family; only `"lognormal"` in Stage 1.
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
                                           cover_factor  = FALSE,
                                           mean_load_pos = 0,
                                           sd_load_pos   = 1.0,
                                           field      = c("icar", "car_proper",
                                                          "bym2"),
                                           rho        = 0.9,
                                           phi        = 0.7,
                                           sigma_pos  = 0.4,
                                           positive   = c("lognormal", "beta"),
                                           seed       = NULL) {
  positive <- match.arg(positive)
  field    <- match.arg(field)
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

  # K shared latent factors W (N x K), each drawn at unit-marginal scale (the geo
  # mean of the field's marginal variances is 1, the Sorbye-Rue convention
  # simulate_occu_cover uses) so the loadings carry the amplitude. An icar field
  # uses the improper eigen draw (K = 1 reproduces the single-field stream); a
  # proper-CAR field draws from R(rho)^{-1}; a bym2 field is the convolution
  # sqrt(1-phi) v + sqrt(phi) u*. The car / bym2 branches are gated, so the
  # default icar stream stays byte-identical.
  W <- matrix(0, N, K)
  if (identical(field, "car_proper")) {
    spec    <- .ms_ocs_field_spec(adj, "car_proper")
    R       <- .ms_ocs_field_R(spec, rho)
    scale_c <- exp(mean(log(diag(solve(R)))))    # geo-mean marginal variance
    Rchol   <- chol(R)                           # R = U'U
    for (k in seq_len(K)) {
      wk <- backsolve(Rchol, stats::rnorm(N))    # ~ N(0, R^{-1})
      W[, k] <- wk / sqrt(scale_c)
    }
  } else if (identical(field, "bym2")) {
    # u* = unit-marginal scaled ICAR; v = iid N(0,1); both unit marginal, so the
    # convolution b = sqrt(1-phi) v + sqrt(phi) u* is unit marginal too.
    Q       <- .occu_cover_icar_Q(adj)
    scale_q <- .occu_cover_icar_scale(adj)
    eig  <- eigen(Q, symmetric = TRUE)
    keep <- eig$values > 1e-8
    for (k in seq_len(K)) {
      z_white <- stats::rnorm(sum(keep))
      u <- as.numeric(eig$vectors[, keep, drop = FALSE] %*%
                        (z_white / sqrt(eig$values[keep])))
      u <- (u - mean(u)) / sqrt(scale_q)         # scaled ICAR component u*
      v <- stats::rnorm(N)                       # iid component
      W[, k] <- sqrt(1 - phi) * v + sqrt(phi) * u
    }
  } else {
    Q       <- .occu_cover_icar_Q(adj)
    scale_q <- .occu_cover_icar_scale(adj)
    eig  <- eigen(Q, symmetric = TRUE)
    keep <- eig$values > 1e-8
    for (k in seq_len(K)) {
      z_white <- stats::rnorm(sum(keep))
      wk <- as.numeric(eig$vectors[, keep, drop = FALSE] %*%
                         (z_white / sqrt(eig$values[keep])))
      wk <- wk - mean(wk)
      W[, k] <- wk / sqrt(scale_q)
    }
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

  # Cover-arm shared-factor loadings (Stage 3): the SAME fields W also load on the
  # cover predictor through a free S x K loading matrix L_pos. Drawn only when
  # requested, after every Stage 1-2 draw, so the no-cover-factor RNG stream (and
  # every existing fixture) stays byte-identical.
  Lpos <- if (isTRUE(cover_factor)) {
    matrix(stats::rnorm(n_species * K, mean_load_pos, sd_load_pos), n_species, K)
  } else NULL

  # Canonical sign anchor: make each factor's diagonal loading positive, flipping
  # that factor's (W[, k], L[, k]) pair together -- and the cover loading L_pos[, k]
  # with them, since it shares the field (F_pos = W L_pos' is sign-invariant).
  for (k in seq_len(K)) {
    if (L[k, k] < 0) {
      L[, k] <- -L[, k]; W[, k] <- -W[, k]
      if (!is.null(Lpos)) Lpos[, k] <- -Lpos[, k]
    }
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
    if (!is.null(Lpos)) eta_pos <- eta_pos + as.numeric(W %*% Lpos[s, ])
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
      field = field,
      rho = if (identical(field, "car_proper")) rho else NULL,
      phi = if (identical(field, "bym2")) phi else NULL,
      mu_occ = mu_occ, mu_p = mu_p, mu_pos = mu_pos,
      sd_occ = sd_occ, sd_p = sd_p, sd_pos = sd_pos,
      sigma_pos = sigma_pos,
      b_occ = b_occ, b_p = b_p, b_pos = b_pos,
      # K = 1 keeps the Stage-1 shapes (loading vector, field vector); K > 1
      # returns the S x K loading matrix and the N x K field matrix.
      L = if (K == 1L) L[, 1L] else L,
      w = if (K == 1L) W[, 1L] else W,
      cover_factor = isTRUE(cover_factor),
      L_pos = if (is.null(Lpos)) NULL else if (K == 1L) Lpos[, 1L] else Lpos,
      psi = psi, z = zmat
    )
  )
}


# ---------------------------------------------------------------------------
# Areal field structure (icar / proper-CAR / BYM2)
# ---------------------------------------------------------------------------
#
# The K shared latent fields all share ONE areal structure, abstracted here so
# the reduced-rank machinery is not hardwired to a single GMRF. A field type
# fixes the precision (up to a per-factor scale tau) as tau * R(h), where h is an
# optional per-factor hyperparameter estimated by EM (a profile M-step) alongside
# tau; the Laplace evidence then treats (tau, h) as the type-II MLEs of an
# empirical-Bayes field prior. The three types:
#
#   icar       improper, fixed R = Q = D - A (rank N - 1, constant null space).
#   car_proper proper, R(rho) = D - rho A, PD for rho in [0, 1) (rho -> 1 = ICAR);
#              hyperparameter rho (correlation), rank N.
#   bym2       Riebler 2016 reparameterised convolution. The combined effect
#              b = (1/sqrt(tau)) (sqrt(1-phi) v + sqrt(phi) u*) (v iid, u* the
#              SCALED ICAR) has marginal precision tau * R(phi) with
#              R(phi) = [(1-phi) I + phi Sigma_u]^{-1}, Sigma_u = Q*^+ the scaled-
#              ICAR generalised inverse. The likelihood sees only the combined b,
#              so the structured / unstructured split need not be represented:
#              the marginal precision IS the prior. Hyperparameter phi (the
#              spatial variance fraction), rank N (PD for phi < 1).
#
# Eigendecomposing Sigma_u once makes R(phi) = V diag(1/((1-phi)+phi s)) V' and
# log|R(phi)| = -sum log((1-phi)+phi s) cheap to evaluate over the phi profile.

# log pseudo-determinant of a GMRF precision (sum of log nonzero eigenvalues).
.ms_ocs_logpdet_from_Q <- function(Q) {
  ev <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  sum(log(ev[ev > max(ev) * 1e-8]))
}

# Field-structure spec: everything the prior / M-step / evidence need to treat a
# field generically. `Q` = D - A (ICAR structure, also the simulator base); `A`,
# `deg` the adjacency / degree; `rank` the GMRF rank (N - 1 improper icar, N
# otherwise); `has_hyper` / `hyper_name` / `hyper_bounds` the per-factor
# hyperparameter. car_proper precomputes `gamma` (eigenvalues of the normalised
# adjacency D^{-1/2} A D^{-1/2}, so log|R(rho)| = sum(log deg) + sum(log(1 -
# rho gamma))); bym2 precomputes the eigenbasis `V` and eigenvalues `s` of the
# scaled-ICAR covariance Sigma_u.
.ms_ocs_field_spec <- function(adj, type = c("icar", "car_proper", "bym2")) {
  type <- match.arg(type)
  N   <- nrow(adj)
  deg <- rowSums(adj)
  if (any(deg == 0)) {
    iso <- which(deg == 0)
    stop(sprintf("areal graph has %d isolated node(s): %s. Connect or drop them.",
                 length(iso), paste(utils::head(iso, 5L), collapse = ", ")),
         call. = FALSE)
  }
  A <- unname(as.matrix(adj)); Q <- .occu_cover_icar_Q(adj)
  spec <- list(type = type, N = N, deg = deg, A = A, Q = Q,
               has_hyper = !identical(type, "icar"))
  if (identical(type, "car_proper")) {
    ds <- 1 / sqrt(deg)
    An <- (ds * A) * rep(ds, each = N)            # D^{-1/2} A D^{-1/2}
    spec$gamma       <- eigen((An + t(An)) / 2, symmetric = TRUE,
                              only.values = TRUE)$values
    spec$logdetD     <- sum(log(deg))
    spec$hyper_name  <- "rho"
    spec$hyper_bounds <- c(0, 1 - 1e-4)
    spec$rank        <- N
  } else if (identical(type, "bym2")) {
    # Scaled-ICAR covariance Sigma_u = Q^+ / scale_q (Sorbye-Rue): eigenvectors of
    # Q, eigenvalues s_j = (1/lambda_j)/scale_q on the structured space and 0 on
    # the constant null direction.
    scale_q <- .occu_cover_icar_scale(adj)
    eig <- eigen(Q, symmetric = TRUE)
    pos <- eig$values > max(eig$values) * 1e-8
    s   <- numeric(N); s[pos] <- (1 / eig$values[pos]) / scale_q
    spec$V           <- eig$vectors
    spec$s           <- s
    spec$hyper_name  <- "phi"
    spec$hyper_bounds <- c(0, 1 - 1e-4)
    spec$rank        <- N
  } else {
    spec$hyper_name <- NA_character_
    spec$rank       <- N - 1L
    spec$logpdetQ   <- .ms_ocs_logpdet_from_Q(Q)
  }
  spec
}

# Structure matrix R(h) for a field (h ignored for icar, whose structure is the
# fixed Q). car_proper: D - rho A (diag(A) = 0 -> diagonal = deg). bym2:
# V diag(1/((1-phi)+phi s)) V'.
.ms_ocs_field_R <- function(spec, h = NULL) {
  if (identical(spec$type, "car_proper")) {
    R <- -h * spec$A; diag(R) <- spec$deg
    return(R)
  }
  if (identical(spec$type, "bym2")) {
    m <- (1 - h) + h * spec$s
    return(spec$V %*% (t(spec$V) * (1 / m)))
  }
  spec$Q
}

# log (pseudo) determinant of R(h): fixed pseudo-det (icar), eigenvalue-sum forms
# for car_proper / bym2.
.ms_ocs_field_logdetR <- function(spec, h = NULL) {
  if (identical(spec$type, "car_proper")) {
    return(spec$logdetD + sum(log1p(-h * spec$gamma)))
  }
  if (identical(spec$type, "bym2")) {
    return(-sum(log((1 - h) + h * spec$s)))
  }
  spec$logpdetQ
}

# Per-field-type closed-form / profile M-step at the field second moment
# Sk = W_k W_k' + Cov(W_k). Returns the updated hyperparameter `h` (NULL for
# icar), the precision scale `tau = rank / tr(R(h) Sk)`, and `quad = tr(R(h) Sk)`.
# The hyperparameter profiles 0.5 log|R(h)| - 0.5 rank log tr(R(h) Sk) (tau
# profiled out) over hyper_bounds.
.ms_ocs_field_mstep <- function(spec, Wk, Cov_wk) {
  N <- spec$N
  if (identical(spec$type, "car_proper")) {
    cD <- sum(spec$deg * (Wk^2 + diag(Cov_wk)))                  # tr(D Sk)
    cA <- as.numeric(Wk %*% (spec$A %*% Wk)) + sum(spec$A * Cov_wk)  # tr(A Sk)
    f  <- function(h) 0.5 * .ms_ocs_field_logdetR(spec, h) -
                      0.5 * N * log(max(cD - h * cA, 1e-8))
    h    <- stats::optimize(f, spec$hyper_bounds, maximum = TRUE)$maximum
    quad <- max(cD - h * cA, 1e-8)
    return(list(h = h, tau = N / quad, quad = quad))
  }
  if (identical(spec$type, "bym2")) {
    # q_j = v_j' Sk v_j = (v_j' Wk)^2 + v_j' Cov v_j, so tr(R(phi) Sk) =
    # sum_j q_j / ((1-phi) + phi s_j).
    Vw <- as.numeric(crossprod(spec$V, Wk))
    q  <- Vw^2 + colSums(spec$V * (Cov_wk %*% spec$V))
    f  <- function(h) {
      m <- (1 - h) + h * spec$s
      -0.5 * sum(log(m)) - 0.5 * N * log(max(sum(q / m), 1e-8))
    }
    h    <- stats::optimize(f, spec$hyper_bounds, maximum = TRUE)$maximum
    quad <- max(sum(q / ((1 - h) + h * spec$s)), 1e-8)
    return(list(h = h, tau = N / quad, quad = quad))
  }
  # icar: rank-deficient (N - 1) precision update, no hyperparameter.
  quad <- as.numeric(Wk %*% (spec$Q %*% Wk)) + sum(spec$Q * Cov_wk)
  list(h = NULL, tau = (N - 1L) / max(quad, 1e-8), quad = quad)
}

# Per-factor structure list (length K) for the objective's field prior: the same
# Q for every icar factor, or the per-factor R(h_k) for a hyperparameterised type.
.ms_ocs_build_field_R <- function(model, hyper_w = NULL) {
  spec <- model$field_spec
  K    <- model$K %||% 1L
  if (is.null(spec) || !isTRUE(spec$has_hyper)) {
    return(rep(list(model$icar_Q), K))
  }
  hyper_w <- rep_len(hyper_w %||% 0.5, K)
  lapply(hyper_w, function(h) .ms_ocs_field_R(spec, h))
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
                                              cover_factor = FALSE,
                                              field_type = "icar",
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
  model$model_type   <- "ms_occu_cover_spatial"
  model$K            <- K
  model$cover_factor <- isTRUE(cover_factor)
  model$adj          <- adj
  model$field_type   <- field_type
  model$field_spec   <- .ms_ocs_field_spec(adj, field_type)
  model$icar_Q       <- model$field_spec$Q     # ICAR base (objective fallback,
  model$icar_scale   <- .occu_cover_icar_scale(adj)  # simulator scale)
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
  K  <- model$K %||% 1L
  S  <- model$n_species
  cf <- isTRUE(model$cover_factor)
  list(P_occ = P_occ, P_p = P_p, P_pos = P_pos, P = P_occ + P_p + P_pos,
       S = S, N = model$n_sites, K = K,
       cover_factor = cf, Lpos_w = if (cf) S * K else 0L,
       occ_idx = seq_len(P_occ), p_idx = P_occ + seq_len(P_p),
       pos_idx = P_occ + P_p + seq_len(P_pos))
}

.ms_ocs_unpack <- function(par, d) {
  off <- 0L
  mu  <- par[off + seq_len(d$P)]; off <- off + d$P
  b   <- lapply(seq_len(d$S), function(s) par[off + (s - 1L) * d$P + seq_len(d$P)])
  off <- off + d$S * d$P
  L   <- matrix(par[off + seq_len(d$S * d$K)], d$S, d$K); off <- off + d$S * d$K
  Lpos <- if (d$cover_factor) {
    m <- matrix(par[off + seq_len(d$Lpos_w)], d$S, d$K); off <- off + d$Lpos_w; m
  } else NULL
  W   <- matrix(par[off + seq_len(d$N * d$K)], d$N, d$K); off <- off + d$N * d$K
  ld  <- par[off + 1L]
  list(mu = mu, b = b, L = L, Lpos = Lpos, W = W, ld = ld)
}

# Penalised joint log-likelihood and its gradient w.r.t. par. `Sinv` is the
# block-diagonal inverse community covariance over (occ, p, pos); `Pmu` the weak
# Gaussian precision on the community means; `inv_sdL2 = 1/sd_L^2` the loading
# ridge; `tau_w` the per-factor field precisions (scalar recycled, or length K;
# prior sum_k 0.5 * tau_w[k] * W[, k]' Q W[, k]). The shared-factor offset on the
# occupancy predictor is sum_k L[s, k] * W[, k] = W %*% L[s, ]; when the model
# carries a cover-arm factor (model$cover_factor), the SAME fields W also load on
# the cover (positive) predictor through a free loading matrix Lpos (S x K),
# adding the per-cell offset sum_k Lpos[s, k] * W[, k] to every visit of the cell.
# Both loading matrices share the weakly-informative ridge inv_sdL2; the field's
# scale / rotation / sign is anchored by the occupancy loadings (triangular in the
# constrained parameterisation), so Lpos rides free.
.ms_ocs_penll_grad <- function(model, par, Sinv, Pmu, inv_sdL2, tau_w,
                               grad = TRUE) {
  d <- .ms_ocs_dims(model); up <- .ms_ocs_unpack(par, d)
  mu <- up$mu; b <- up$b; L <- up$L; Lpos <- up$Lpos; W <- up$W; ld <- up$ld
  # Per-factor field structure R_k (D - rho_k A for car_proper; the fixed ICAR
  # Q for every factor on the icar / fallback path -- byte-identical in value).
  Rk <- model$field_R %||% rep(list(model$icar_Q), d$K)
  tau_w <- rep_len(tau_w, d$K)
  cl <- .tobs_clamp_eta

  ll <- 0
  g_mu <- numeric(d$P); g_b <- vector("list", d$S)
  g_L  <- matrix(0, d$S, d$K); g_W <- matrix(0, d$N, d$K); g_ld <- 0
  g_Lpos <- if (d$cover_factor) matrix(0, d$S, d$K) else NULL
  for (s in seq_len(d$S)) {
    v   <- .ms_occu_cover_species_view(model, s)
    th  <- mu + b[[s]]
    eta <- .occu_cover_eta_from_par(v, th[d$occ_idx], th[d$p_idx], th[d$pos_idx])
    # Inject the shared-factor offset on the occupancy predictor.
    eta$psi <- stats::plogis(cl(as.numeric(v$X_occ %*% th[d$occ_idx]) +
                                  as.numeric(W %*% L[s, ])))
    # ...and, with a cover-arm factor, on the cover predictor: one per-cell offset
    # broadcast across the cell's visits (ep_mat is n_sites x max_visits, so a
    # length-N vector recycles down the rows = cells).
    if (d$cover_factor) {
      eta$ep_mat <- eta$ep_mat + as.numeric(W %*% Lpos[s, ])
    }
    ll <- ll + sum(.occu_cover_site_ll(v, eta$psi, eta$p_mat, eta$ep_mat, ld))
    if (grad) {
      eg <- .occu_cover_eta_grad(v, eta$psi, eta$p_mat, eta$ep_mat, ld)
      cg <- .occu_cover_coef_grad(v, eg)              # c(g_occ, g_det, g_pos, g_ld)
      cg_coef <- cg[seq_len(d$P)]
      g_mu    <- g_mu + cg_coef
      g_b[[s]] <- cg_coef                              # RE prior added below
      g_L[s, ] <- as.numeric(crossprod(W, eg$g_psi))  # t(W) %*% g_psi_s
      g_W      <- g_W + outer(eg$g_psi, L[s, ])
      if (d$cover_factor) {
        # The cover-offset gradient is the per-cell sum of the per-visit cover
        # eta-gradient (the offset enters every visit of the cell identically).
        gpos_cell    <- rowSums(eg$g_pos)
        g_Lpos[s, ]  <- as.numeric(crossprod(W, gpos_cell))
        g_W          <- g_W + outer(gpos_cell, Lpos[s, ])
      }
      g_ld    <- g_ld + cg[d$P + 1L]
    }
  }

  # ---- priors (penalty) ----
  ll <- ll - 0.5 * sum(vapply(b, function(bb) as.numeric(bb %*% Sinv %*% bb), 0))
  ll <- ll - 0.5 * as.numeric(mu %*% Pmu %*% mu)
  ll <- ll - 0.5 * inv_sdL2 * sum(L^2)
  if (d$cover_factor) ll <- ll - 0.5 * inv_sdL2 * sum(Lpos^2)
  QW <- matrix(0, d$N, d$K)                            # R_k W[, k] per factor
  for (k in seq_len(d$K)) QW[, k] <- as.numeric(Rk[[k]] %*% W[, k])
  ll <- ll - 0.5 * sum(tau_w * colSums(W * QW))

  if (!grad) return(list(ll = ll))

  g_mu <- g_mu - as.numeric(Pmu %*% mu)
  for (s in seq_len(d$S)) g_b[[s]] <- g_b[[s]] - as.numeric(Sinv %*% b[[s]])
  g_L  <- g_L - inv_sdL2 * L
  if (d$cover_factor) g_Lpos <- g_Lpos - inv_sdL2 * Lpos
  g_W  <- g_W - sweep(QW, 2L, tau_w, "*")

  grad_vec <- c(g_mu, unlist(g_b), as.numeric(g_L),
                if (d$cover_factor) as.numeric(g_Lpos) else numeric(0),
                as.numeric(g_W), g_ld)
  list(ll = ll, grad = grad_vec)
}


# ---------------------------------------------------------------------------
# Identifiability-constrained loading parameterisation (gllvm/HMSC)
# ---------------------------------------------------------------------------
#
# For K > 1 the unconstrained loadings are identified only up to an orthogonal
# rotation, so the unconstrained posterior is improper along the rotation
# manifold -- fine for the rotation-invariant point quantities (F = W L', psi),
# but it corrupts per-factor uncertainty and any posterior-based model
# comparison (WAIC explodes). The standard fix is the lower-triangular,
# positive-diagonal constraint: factor k loads on species k..S only, and
# L[k, k] = exp(l_kk) > 0. That uniquely fixes the rotation, giving a full-rank
# Hessian. The free loading vector packs, column by column,
#   [l_kk (log-diagonal), L[k+1, k], ..., L[S, k]].

.ms_ocs_lfree_dim <- function(S, K) as.integer(K * S - K * (K - 1L) / 2L)

# Free triangular loading vector -> full S x K loading matrix (zeros above the
# diagonal, exp() on the diagonal).
.ms_ocs_lfree_to_L <- function(lfree, S, K) {
  L <- matrix(0, S, K); pos <- 0L
  for (k in seq_len(K)) {
    nfree <- S - k + 1L
    blk <- lfree[pos + seq_len(nfree)]; pos <- pos + nfree
    L[k, k] <- exp(blk[1L])
    if (nfree > 1L) L[(k + 1L):S, k] <- blk[-1L]
  }
  L
}

# Full S x K loading matrix -> free triangular vector (log the diagonal).
.ms_ocs_L_to_lfree <- function(L, S, K) {
  out <- numeric(0)
  for (k in seq_len(K)) {
    diag_l <- log(max(L[k, k], 1e-6))
    off    <- if (k < S) L[(k + 1L):S, k] else numeric(0)
    out <- c(out, diag_l, off)
  }
  out
}

# Gradient w.r.t. full L (S x K) -> gradient w.r.t. the free triangular vector.
# The structural zeros (s < k) are dropped; the diagonal carries the log-link
# chain factor dL[k,k]/dl_kk = exp(l_kk) = L[k, k].
.ms_ocs_gL_to_glfree <- function(gL, L, S, K) {
  out <- numeric(0)
  for (k in seq_len(K)) {
    diag_g <- gL[k, k] * L[k, k]
    off_g  <- if (k < S) gL[(k + 1L):S, k] else numeric(0)
    out <- c(out, diag_g, off_g)
  }
  out
}

# Constrained penalised log-lik + gradient in the triangular parameterisation.
# par_c = c(mu[P], b[S*P], lfree[nL], vec(W)[N*K], log_disp); it expands lfree to
# the full L, defers to the (validated) unconstrained objective, then maps the
# L-block gradient back to lfree by the chain rule. No likelihood math is
# duplicated -- this is purely the reparameterisation adapter.
.ms_ocs_penll_grad_c <- function(model, par_c, Sinv, Pmu, inv_sdL2, tau_w,
                                 grad = TRUE) {
  d <- .ms_ocs_dims(model); S <- d$S; K <- d$K; P <- d$P
  nL <- .ms_ocs_lfree_dim(S, K)
  head_n <- P + S * P                       # mu + b prefix (shared layout)
  tail_n <- d$Lpos_w + d$N * K + 1L         # vec(Lpos) + vec(W) + log_disp suffix

  pre   <- par_c[seq_len(head_n)]
  lfree <- par_c[head_n + seq_len(nL)]
  suf   <- par_c[head_n + nL + seq_len(tail_n)]
  L     <- .ms_ocs_lfree_to_L(lfree, S, K)
  par_full <- c(pre, as.numeric(L), suf)

  res <- .ms_ocs_penll_grad(model, par_full, Sinv, Pmu, inv_sdL2, tau_w,
                            grad = grad)
  if (!grad) return(list(ll = res$ll))

  g <- res$grad
  g_pre <- g[seq_len(head_n)]
  g_L   <- matrix(g[head_n + seq_len(S * K)], S, K)
  g_suf <- g[head_n + S * K + seq_len(tail_n)]
  g_lfree <- .ms_ocs_gL_to_glfree(g_L, L, S, K)
  list(ll = res$ll, grad = c(g_pre, g_lfree, g_suf))
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
# Unpack a constrained packed par_c = c(mu, b, lfree, vec(W), log_disp) into the
# same fields as .ms_ocs_unpack, expanding lfree -> the lower-triangular L.
.ms_ocs_unpack_c <- function(par_c, d) {
  S <- d$S; K <- d$K; P <- d$P; N <- d$N; nL <- .ms_ocs_lfree_dim(S, K)
  off <- 0L
  mu  <- par_c[off + seq_len(P)]; off <- off + P
  b   <- lapply(seq_len(S), function(s) par_c[off + (s - 1L) * P + seq_len(P)])
  off <- off + S * P
  lfree <- par_c[off + seq_len(nL)]; off <- off + nL
  Lpos <- if (d$cover_factor) {
    m <- matrix(par_c[off + seq_len(d$Lpos_w)], S, K); off <- off + d$Lpos_w; m
  } else NULL
  W   <- matrix(par_c[off + seq_len(N * K)], N, K); off <- off + N * K
  ld  <- par_c[off + 1L]
  list(mu = mu, b = b, L = .ms_ocs_lfree_to_L(lfree, S, K), Lpos = Lpos,
       W = W, ld = ld)
}

.ms_ocs_inner_mode <- function(model, Sigma, sd_L = 1.0, tau_w = 1.0,
                               par_init = NULL, sigma.beta = 5,
                               maxit = 400L, hessian = FALSE,
                               constrain = FALSE) {
  d <- .ms_ocs_dims(model); K <- d$K
  Sinv <- matrix(0, d$P, d$P)
  Sinv[d$occ_idx, d$occ_idx] <- solve(Sigma$occ)
  Sinv[d$p_idx,   d$p_idx]   <- solve(Sigma$p)
  Sinv[d$pos_idx, d$pos_idx] <- solve(Sigma$pos)
  Pmu      <- diag(1 / sigma.beta^2, d$P)
  inv_sdL2 <- 1 / sd_L^2
  objective <- if (constrain) .ms_ocs_penll_grad_c else .ms_ocs_penll_grad

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
    if (constrain) {
      # Project the EOF loadings to the lower-triangular, positive-diagonal
      # canonical form for the constrained start: zero the upper triangle, make
      # each diagonal positive (flipping its field column), and pack to lfree.
      L0[upper.tri(L0)] <- 0
      for (k in seq_len(K)) {
        if (L0[k, k] < 0) { L0[, k] <- -L0[, k]; W0[, k] <- -W0[, k] }
        if (abs(L0[k, k]) < 1e-3) L0[k, k] <- 0.1
      }
      par_init <- c(mu0, numeric(d$S * d$P),
                    .ms_ocs_L_to_lfree(L0, d$S, K),
                    if (d$cover_factor) numeric(d$Lpos_w) else numeric(0),
                    as.numeric(W0), ld0)
    } else {
      par_init <- c(mu0, numeric(d$S * d$P), as.numeric(L0),
                    if (d$cover_factor) numeric(d$Lpos_w) else numeric(0),
                    as.numeric(W0), ld0)
    }
  }

  fn <- function(p) -objective(model, p, Sinv, Pmu, inv_sdL2, tau_w,
                               grad = FALSE)$ll
  gr <- function(p) -objective(model, p, Sinv, Pmu, inv_sdL2, tau_w,
                               grad = TRUE)$grad
  opt <- stats::optim(par_init, fn, gr, method = "BFGS",
                      control = list(maxit = maxit, reltol = 1e-10),
                      hessian = hessian)
  if (constrain) {
    # The triangular parameterisation is already identified (positive diagonal,
    # zeros above) -- no post-hoc sign / rotation step needed.
    up <- .ms_ocs_unpack_c(opt$par, d)
  } else {
    up <- .ms_ocs_unpack(opt$par, d)
    # Canonical sign anchor: make each factor's diagonal loading L[k, k]
    # positive, flipping (L[, k], W[, k]) -- and the cover loading Lpos[, k], which
    # shares the field -- together. The penalised objective is even in this joint
    # flip, so the Hessian over the packed par stays valid for the M-step
    # covariances (the reported field / loadings just adopt one canonical sign).
    for (k in seq_len(K)) {
      if (up$L[k, k] < 0) {
        up$L[, k] <- -up$L[, k]; up$W[, k] <- -up$W[, k]
        if (!is.null(up$Lpos)) up$Lpos[, k] <- -up$Lpos[, k]
      }
    }
  }
  Lpos_out <- if (is.null(up$Lpos)) NULL else if (K == 1L) as.numeric(up$Lpos) else up$Lpos
  list(mu = up$mu, b = up$b, ld = up$ld, constrained = constrain,
       L = if (K == 1L) as.numeric(up$L) else up$L,
       w = if (K == 1L) as.numeric(up$W) else up$W,
       Lpos = Lpos_out,
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
# 1982) and the per-factor field update (.ms_ocs_field_mstep): the precision
# scale tau_w[k] = rank / tr(R(h_k) Sk) and, for a hyperparameterised field
# (car_proper rho, bym2 phi), the profile update of h_k. The loading SD `sd_L`
# fixes the amplitude split and is held at its supplied value (a hyperparameter).
.tobs_fit_ms_occu_cover_spatial <- function(model, sd_L = 1.0,
                                            max.em = 30L, tol = 1e-3,
                                            sigma.beta = 5, verbose = FALSE,
                                            constrain = FALSE) {
  d <- .ms_ocs_dims(model); K <- d$K
  spec <- model$field_spec %||%
    .ms_ocs_field_spec(model$adj, model$field_type %||% "icar")
  model$field_spec <- spec
  hyper_w <- if (isTRUE(spec$has_hyper)) rep(0.5, K) else NULL
  # Per-factor structure fed to the E-step objective; rebuilt after each M-step.
  model$field_R <- .ms_ocs_build_field_R(model, hyper_w)

  Sigma <- list(occ = diag(0.3^2, d$P_occ), p = diag(0.3^2, d$P_p),
                pos = diag(0.3^2, d$P_pos))
  tau_w <- rep(1.0, K)
  par_init <- NULL
  prev <- -Inf

  # Packed offsets of each species' b block and the K shared field blocks. The
  # field blocks sit after the loadings, whose width depends on the
  # parameterisation: S*K unconstrained, or the triangular free count when
  # constrained.
  L_width <- if (constrain) .ms_ocs_lfree_dim(d$S, K) else d$S * K
  b_off <- function(s) d$P + (s - 1L) * d$P
  w_off <- d$P + d$S * d$P + L_width + d$Lpos_w

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
                              hessian = TRUE, constrain = constrain)
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

    # M-step: per-factor field hyperparameters at the field second moment
    # Sk = W_k W_k' + Cov(W_k). .ms_ocs_field_mstep returns the precision scale
    # tau_k = rank / tr(R(h_k) Sk) and, for a hyperparameterised field, the
    # profiled h_k (car_proper rho / bym2 phi); an icar field keeps the
    # rank-deficient (N - 1) update with no hyperparameter.
    Wm <- matrix(fit$w, d$N, K)
    for (k in seq_len(K)) {
      widx   <- w_off + (k - 1L) * d$N + seq_len(d$N)
      Cov_wk <- Cov[widx, widx, drop = FALSE]
      up_k   <- .ms_ocs_field_mstep(spec, Wm[, k], Cov_wk)
      tau_w[k] <- up_k$tau
      if (!is.null(up_k$h)) hyper_w[k] <- up_k$h
    }
    model$field_R <- .ms_ocs_build_field_R(model, hyper_w)

    if (verbose) {
      cat(sprintf("EM %2d  logpen=%.3f  tau_w[1]=%.3f%s  sd_occ1=%.3f\n",
                  em, fit$logpen, tau_w[1L],
                  if (isTRUE(spec$has_hyper))
                    sprintf("  %s[1]=%.3f", spec$hyper_name, hyper_w[1L]) else "",
                  sqrt(Sigma$occ[1L, 1L])))
    }
    if (is.finite(prev) && abs(fit$logpen - prev) < tol * (abs(prev) + tol)) {
      prev <- fit$logpen; break
    }
    prev <- fit$logpen
  }

  c(fit, list(Sigma = Sigma, tau_w = tau_w, field_hyper = hyper_w,
              field_type = spec$type, hyper_name = spec$hyper_name,
              sd_L = sd_L, em_logpen = prev, cov = Cov))
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
  L_width <- if (isTRUE(fit$constrained)) .ms_ocs_lfree_dim(d$S, K) else d$S * K
  w_off  <- P + d$S * P + L_width + d$Lpos_w
  widx   <- w_off + seq_len(d$N * K)
  field_sd <- matrix(sqrt(pmax(diag(Cov)[widx], 0)), d$N, K)
  L <- fit$L

  # Cover-arm loadings (the same shared fields W loading on the cover predictor),
  # carried from the inner mode in the SAME canonical sign as the reported field
  # (so the cover spatial contribution F_pos = W Lpos' has the right sign). K = 1
  # keeps the length-S vector shape; K > 1 returns S x K.
  Lpos <- fit$Lpos
  if (!is.null(Lpos)) {
    if (K == 1L) {
      Lpos <- as.numeric(Lpos); names(Lpos) <- model$species_names
    } else {
      Lpos <- matrix(Lpos, d$S, K)
      rownames(Lpos) <- model$species_names
      colnames(Lpos) <- paste0("factor", seq_len(K))
    }
  }
  rot <- NULL
  if (K == 1L) {
    field_sd <- as.numeric(field_sd)
    names(L) <- model$species_names
  } else {
    rownames(L) <- model$species_names
    colnames(L) <- paste0("factor", seq_len(K))
    colnames(field_sd) <- paste0("factor", seq_len(K))
    # Post-hoc varimax rotation for interpretable factors. The fit is invariant
    # to an orthogonal rotation R of the factors (F = W L' = (W R)(L R)'), so
    # rotate the loadings to a simple structure and the fields by the same R --
    # the predictor, the psi posterior, and every recovery quantity are
    # unchanged; only the per-factor labelling becomes interpretable.
    vm <- tryCatch(stats::varimax(L, normalize = FALSE), error = function(e) NULL)
    if (!is.null(vm)) {
      R <- matrix(as.numeric(vm$rotmat), K, K)
      L_rot <- L %*% R; W_rot <- fit$w %*% R
      dimnames(L_rot) <- dimnames(L)
      colnames(W_rot) <- colnames(L)
      # The cover loadings share W, so rotate them by the same R to keep the
      # cover spatial contribution F_pos = W Lpos' invariant under the relabelling.
      Lpos_rot <- if (!is.null(Lpos)) {
        lr <- Lpos %*% R; dimnames(lr) <- dimnames(Lpos); lr
      } else NULL
      rot <- list(rotmat = R, loadings = L_rot, field = W_rot,
                  loadings_cover = Lpos_rot)
    }
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
    spatial      = {
      ft   <- fit$field_type %||% model$field_type %||% "icar"
      hypr <- fit$field_hyper
      list(
        type     = if (isTRUE(d$cover_factor)) paste0(ft, "+cover") else ft,
        field_type = ft,
        K        = K,
        field    = fit$w,
        field_sd = field_sd,
        loadings = L,
        loadings_cover = Lpos,
        cover_factor   = isTRUE(d$cover_factor),
        rotation = rot,
        tau_w    = fit$tau_w,
        # The field hyperparameter is reported under its conventional name: the
        # CAR correlation rho or the BYM2 spatial-variance fraction phi.
        rho_w    = if (identical(ft, "car_proper")) hypr else NULL,
        phi_w    = if (identical(ft, "bym2"))       hypr else NULL,
        sd_L     = fit$sd_L,
        associations = .ms_ocs_associations(fit, d, model$species_names),
        maps     = .ms_ocs_map_summary(model, fit)
      )
    },
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


# Joint per-species per-cell posterior of the three map quantities for the
# K-factor spatial community fit. Samples the packed inner latent from its
# Gaussian Laplace posterior N(mode, Cov) -- via an eigen PSD square root, since
# the joint Hessian is only PSD in the confounded ICAR-level / intercept (and,
# for K > 1, rotational) directions -- and maps each draw through the per-arm
# predictors: occupancy eta_occ = X_occ (mu_occ + b_s_occ) + sum_k L_sk w_kc, and
# cover eta_pos = X_pos (mu_pos + b_s_pos) [+ sum_k Lpos_sk w_kc, with a cover-arm
# factor]. Both depend on (L, W) / (Lpos, W) only through W L[s, ]', so they are
# invariant to the per-factor sign and (K > 1) rotational symmetries and need no
# anchoring. Returns [n_draws x N x S] arrays for psi (occupancy), cover_cond (the
# hurdle conditional cover mean E[cover | present] -- lognormal exp(eta + sigma^2
# / 2), beta plogis(eta)), and cover_exp (the unconditional expected cover psi *
# cover_cond). The derived covers are formed per draw (the nonlinear cover mean +
# the psi product), then summarised -- the marginalize-derived-quantities rule,
# not a plug-in of the posterior-mean eta.
# Draw the packed inner latent from its Gaussian Laplace posterior N(par, cov)
# via an eigen PSD square root (the joint Hessian is only PSD along the confounded
# ICAR-level / intercept / rotational directions). Returns an [n_draws x npar]
# matrix. Shared by the map posterior and the pointwise-log-likelihood kernel.
.ms_ocs_draw_par <- function(par, cov, n_draws) {
  Sg <- (cov + t(cov)) / 2
  eg <- eigen(Sg, symmetric = TRUE)
  rt <- eg$vectors %*% (sqrt(pmax(eg$values, 0)) * t(eg$vectors))   # PSD sqrt
  npar <- length(par)
  Z <- matrix(stats::rnorm(n_draws * npar), n_draws, npar)
  sweep(Z %*% rt, 2L, par, "+")
}

# Posterior draws of the packed inner latent (optionally a column subset `cols`)
# from a NUTS fit: the actual draws (the exact posterior), thinned to ~n_draws.
# Returns NULL for a non-NUTS fit so the caller falls back to the Gaussian
# Laplace draw. Keeps the NUTS associations / maps on the same exact-posterior
# draws the WAIC uses, rather than a moment-matched N(par, cov).
.ms_ocs_nuts_inner <- function(fit, n_draws, cols = NULL) {
  nd <- fit$nuts
  if (is.null(nd) || is.null(nd$draws)) return(NULL)
  inner <- nd$draws[, nd$layout$inner, drop = FALSE]
  if (!is.null(cols)) inner <- inner[, cols, drop = FALSE]
  M <- nrow(inner)
  if (!is.null(n_draws) && as.integer(n_draws) < M) {
    idx <- unique(round(seq(1, M, length.out = as.integer(n_draws))))
    inner <- inner[idx, , drop = FALSE]
  }
  inner
}

.ms_ocs_joint_posterior <- function(model, fit, n_draws = 300L) {
  d <- fit$d
  draws_par <- .ms_ocs_nuts_inner(fit, n_draws) %||%
    .ms_ocs_draw_par(fit$par, fit$cov, n_draws)
  n_draws <- nrow(draws_par)

  unpack    <- if (isTRUE(fit$constrained)) .ms_ocs_unpack_c else .ms_ocs_unpack
  lognormal <- !identical(model$positive, "beta")
  X_occ <- model$X_occ; X_pos <- model$X_pos_site
  cl <- .tobs_clamp_eta
  psi <- cc <- ce <- array(0, dim = c(n_draws, d$N, d$S))
  for (i in seq_len(n_draws)) {
    up <- unpack(draws_par[i, ], d)
    sigma2 <- exp(up$ld)^2
    for (s in seq_len(d$S)) {
      th     <- up$mu + up$b[[s]]
      p_s    <- stats::plogis(cl(as.numeric(X_occ %*% th[d$occ_idx]) +
                                   as.numeric(up$W %*% up$L[s, ])))
      eta_pos <- as.numeric(X_pos %*% th[d$pos_idx])
      if (d$cover_factor) eta_pos <- eta_pos + as.numeric(up$W %*% up$Lpos[s, ])
      cond <- if (lognormal) exp(cl(eta_pos) + sigma2 / 2) else stats::plogis(cl(eta_pos))
      psi[i, , s] <- p_s
      cc[i, , s]  <- cond
      ce[i, , s]  <- p_s * cond
    }
  }
  list(psi = psi, cover_cond = cc, cover_exp = ce)
}


# Per-cell per-species posterior summary (the spatial-JSDM map output): the mean,
# median, and a central interval of each map quantity at every cell c and species
# s, from the joint-posterior draws above. The shared latent fields let a
# species' map borrow strength across the community, so a rare species gets a
# calibrated map (mean + interval) rather than the ragged empirical rate. Returns
# a list with `psi`, `cover_cond`, and `cover_exp`, each a list of N x S matrices.
.ms_ocs_map_summary <- function(model, fit, n_draws = 300L,
                                probs = c(0.025, 0.975)) {
  post <- .ms_ocs_joint_posterior(model, fit, n_draws)
  sp   <- model$species_names
  summ_one <- function(arr) {
    pull <- function(f) { m <- apply(arr, c(2L, 3L), f); dimnames(m) <- list(NULL, sp); m }
    list(mean   = pull(mean),
         median = pull(stats::median),
         lower  = pull(function(x) stats::quantile(x, probs[1L], names = FALSE)),
         upper  = pull(function(x) stats::quantile(x, probs[2L], names = FALSE)))
  }
  lapply(post, summ_one)
}

# Tidy long form (one row per cell x species) of a requested map quantity,
# returned by predict() for a spatial-factor ms_occu_cover() fit: `type`
# "occupancy" -> psi, "cover_cond" -> conditional cover mean, "cover_exp" ->
# unconditional expected cover.
.tobs_ms_ocs_predict_state <- function(object, type = "occupancy") {
  maps <- object$spatial$maps
  if (is.null(maps))
    stop("This fit carries no per-cell map posterior (refit with a current ",
         "tulpaObs).", call. = FALSE)
  key <- switch(type,
                occupancy = , state = , occurrence = "psi",
                cover_cond = "cover_cond", cover_exp = "cover_exp",
                stop(sprintf(paste0("predict() type '%s' is not available for a ",
                             "spatial-factor ms_occu_cover() fit; use ",
                             "\"occupancy\", \"cover_cond\", or \"cover_exp\"."),
                             type), call. = FALSE))
  m  <- maps[[key]]
  N  <- nrow(m$mean); S <- ncol(m$mean)
  sp <- colnames(m$mean) %||% paste0("sp", seq_len(S))
  nm <- if (identical(key, "psi")) "psi" else key
  out <- data.frame(cell = rep(seq_len(N), times = S),
                    species = rep(sp, each = N), stringsAsFactors = FALSE)
  out[[nm]]                  <- as.numeric(m$mean)
  out[[paste0(nm, "_median")]] <- as.numeric(m$median)
  out[[paste0(nm, "_lower")]]  <- as.numeric(m$lower)
  out[[paste0(nm, "_upper")]]  <- as.numeric(m$upper)
  out
}

# Fitted per-species surfaces (posterior means): the field-augmented occupancy
# psi and conditional cover (from the map posterior, so they carry the latent
# field) and the per-species detection probability p (no field on detection).
# Returns list(psi, p, cover), each N x S -- the shape the non-spatial community
# fitted() uses.
.tobs_fitted_ms_occu_cover_spatial <- function(object) {
  model <- object$model
  maps  <- object$spatial$maps
  if (is.null(maps))
    stop("This fit carries no fitted surfaces (refit with a current tulpaObs).",
         call. = FALSE)
  cm <- object$ms_community
  X_det_site <- model$X_det_site
  p <- stats::plogis(X_det_site %*%
                     t(cm$coef_p[, seq_len(ncol(X_det_site)), drop = FALSE]))
  dimnames(p) <- list(NULL, model$species_names)
  list(psi = maps$psi$mean, p = p, cover = maps$cover_cond$mean)
}

# Plug-in (posterior-mean) data simulation for the spatial-factor community fit:
# per species, draw z ~ Bernoulli(psi) with the field-augmented psi, the visit
# detections, and the cover hurdle at the observed visit pattern. Mirrors the
# non-spatial community simulate, adding the shared field on psi (and, with a
# cover-arm factor, on the cover predictor). Returns 3D arrays matching y / y_pos.
.tobs_simulate_ms_occu_cover_spatial <- function(object, nsim = 1) {
  model <- object$model; cm <- object$ms_community; sp <- object$spatial
  asK <- function(M) if (is.null(dim(M))) matrix(M, ncol = 1L) else M
  W    <- asK(sp$field); L <- asK(sp$loadings)
  Lpos <- if (is.null(sp$loadings_cover)) NULL else asK(sp$loadings_cover)
  n_sites <- model$n_sites; max_visits <- model$max_visits
  n_species <- model$n_species
  is_beta <- identical(model$positive, "beta")
  disp <- exp(object$means[[length(object$means)]])
  cl <- .tobs_clamp_eta

  # Per-species predictors with the shared-factor field offset (W L[s,] on psi,
  # W Lpos[s,] on cover) computed here; the z + detection + cover draws run in
  # cpp_simulate_ms_occu_cover from R's RNG stream in the former order.
  psi <- matrix(0, n_sites, n_species)
  p_mat <- array(0, c(n_sites, max_visits, n_species))
  ep_mat <- array(0, c(n_sites, max_visits, n_species))
  for (s in seq_len(n_species)) {
    eta <- .occu_cover_eta_from_par(model, cm$coef_occ[s, ], cm$coef_p[s, ],
                                    cm$coef_pos[s, ])
    psi[, s] <- stats::plogis(cl(as.numeric(model$X_occ %*% cm$coef_occ[s, ]) +
                                   as.numeric(W %*% L[s, ])))
    if (!is.null(Lpos)) eta$ep_mat <- eta$ep_mat + as.numeric(W %*% Lpos[s, ])
    p_mat[, , s] <- eta$p_mat; ep_mat[, , s] <- eta$ep_mat
  }
  res <- cpp_simulate_ms_occu_cover(psi, as.numeric(p_mat), as.numeric(ep_mat),
    as.integer(model$valid), as.numeric(disp), is_beta,
    n_sites, max_visits, n_species, as.integer(nsim))
  res <- lapply(res, function(r) {
    dn <- list(NULL, NULL, model$species_names)
    dimnames(r$y) <- dn; dimnames(r$y_pos) <- dn; r
  })
  if (nsim == 1L) res[[1]] else res
}


# ---------------------------------------------------------------------------
# Residual species-association matrices (spatial-JSDM / HMSC output)
# ---------------------------------------------------------------------------
#
# The K shared latent fields induce residual co-occurrence among species through
# the loadings: with the fields at unit marginal variance (the simulator's
# Sorbye-Rue convention, amplitude carried by the loadings), the species'
# occupancy factor contribution sum_k L_sk w_kc has cross-species covariance
# Omega_occ = L_occ L_occ' at a unit-variance cell, and the reported association
# is its correlation form R[s, s'] = Omega[s, s'] / sqrt(Omega[s, s] Omega[s', s']).
# This is invariant to the factor rotation / sign / column-scale symmetry
# (L Q Q' L' = L L' for orthogonal Q), so it is identified without anchoring. With
# a cover-arm factor the cover association is corr(L_pos L_pos') and the joint
# cross-arm association is the standardized L_occ L_pos' (species-s occupancy vs
# species-s' cover, sharing the field) -- the genuinely joint-model quantity.
#
# Per the marginalize-derived-quantities rule the matrices are summarised from
# draws of the loading posterior N(mode, cov), not the plug-in mode: each draw's
# correlation matrix is formed, then per-element median + central interval. The
# loading block is contiguous in the (constrained or unconstrained) packed
# coordinates at head_n = P + S*P, so its marginal posterior is the matching
# sub-block of the joint Laplace covariance. The rotation symmetry lives in the
# flat directions of that block, and corr(L L') ignores them, so the interval
# reflects genuine loading uncertainty rather than the gauge freedom.
.ms_ocs_corr_self <- function(L) {
  Om  <- tcrossprod(L)
  s   <- sqrt(pmax(diag(Om), 0))
  den <- outer(s, s)
  R   <- Om / den
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  R
}

.ms_ocs_corr_cross <- function(Locc, Lpos) {
  so  <- sqrt(pmax(rowSums(Locc^2), 0))
  sp  <- sqrt(pmax(rowSums(Lpos^2), 0))
  den <- outer(so, sp)
  R   <- (Locc %*% t(Lpos)) / den
  R[!is.finite(R)] <- 0
  R
}

.ms_ocs_associations <- function(fit, d, species_names,
                                  n_draws = 500L, prob = 0.95) {
  S <- d$S; K <- d$K
  constrained <- isTRUE(fit$constrained)
  head_n  <- d$P + d$S * d$P
  L_width <- if (constrained) .ms_ocs_lfree_dim(S, K) else S * K
  idx  <- head_n + seq_len(L_width + d$Lpos_w)
  mode <- fit$par[idx]
  # Exact NUTS draws of the loading block when present, else the Gaussian Laplace
  # draw -- the loading posterior is non-Gaussian under the triangular constraint,
  # so the NUTS draws give faithful association intervals (cf. the WAIC path).
  draws <- .ms_ocs_nuts_inner(fit, n_draws, idx)
  if (is.null(draws)) {
    V <- fit$cov[idx, idx, drop = FALSE]; V <- (V + t(V)) / 2
    draws <- .occu_cover_rmvn(n_draws, mode, V)
  }
  n_draws <- nrow(draws)

  to_L    <- function(blk) if (constrained) .ms_ocs_lfree_to_L(blk[seq_len(L_width)], S, K)
                           else matrix(blk[seq_len(L_width)], S, K)
  to_Lpos <- function(blk) matrix(blk[L_width + seq_len(d$Lpos_w)], S, K)

  arms <- c("occupancy", if (d$cover_factor) c("cover", "cross"))
  acc  <- stats::setNames(lapply(arms, function(a) array(0, c(n_draws, S, S))), arms)
  for (i in seq_len(n_draws)) {
    Lo <- to_L(draws[i, ])
    acc$occupancy[i, , ] <- .ms_ocs_corr_self(Lo)
    if (d$cover_factor) {
      Lp <- to_Lpos(draws[i, ])
      acc$cover[i, , ] <- .ms_ocs_corr_self(Lp)
      acc$cross[i, , ] <- .ms_ocs_corr_cross(Lo, Lp)
    }
  }

  qlo <- (1 - prob) / 2; qhi <- 1 - qlo
  dn  <- list(species_names, species_names)
  summ <- function(arr, est) {
    med <- apply(arr, c(2L, 3L), stats::median)
    lo  <- apply(arr, c(2L, 3L), stats::quantile, probs = qlo, names = FALSE)
    hi  <- apply(arr, c(2L, 3L), stats::quantile, probs = qhi, names = FALSE)
    dimnames(est) <- dimnames(med) <- dimnames(lo) <- dimnames(hi) <- dn
    list(estimate = est, median = med, lower = lo, upper = hi)
  }
  L0  <- to_L(mode)
  out <- list(occupancy = summ(acc$occupancy, .ms_ocs_corr_self(L0)))
  if (d$cover_factor) {
    Lp0 <- to_Lpos(mode)
    out$cover <- summ(acc$cover, .ms_ocs_corr_self(Lp0))
    out$cross <- summ(acc$cross, .ms_ocs_corr_cross(L0, Lp0))
  }
  out$prob <- prob; out$n_draws <- n_draws
  out
}


#' Residual species-association matrices from a spatial-factor community fit
#'
#' For a reduced-rank spatial-factor community occupancy-cover fit (a
#' [ms_occu_cover()] model carrying an `icar()`, `car_proper()`, or `bym2()`
#' field on the occupancy arm), return the residual species associations the
#' shared latent fields imply -- the spatial-JSDM / HMSC output. The K unit-scale
#' fields induce a residual occupancy covariance `L L'` across species; the
#' reported matrix is its correlation form, identified (invariant to the factor
#' rotation, sign, and column scale). When the fit also carries a cover-arm
#' factor, the cover association `corr(L_pos L_pos')` and the joint cross-arm
#' association (standardized `L_occ L_pos'`, species occupancy vs species cover)
#' are available too.
#'
#' Each matrix is summarised from draws of the loading posterior (not the plug-in
#' mode), so a central interval accompanies the point estimate.
#'
#' @param object A fitted `tobs_fit` from a spatial-factor `ms_occu_cover()`.
#' @param type Which association to return: `"occupancy"` (default), `"cover"`,
#'   or `"cross"` (occupancy vs cover). `"cover"` and `"cross"` require a
#'   cover-arm factor (an `icar()`/`car_proper()`/`bym2()` term on the cover
#'   formula).
#' @param summary One of `"median"`, `"estimate"` (rotation-invariant point
#'   estimate at the posterior mode), `"lower"`, `"upper"`; or `NULL` (default)
#'   to return the full list of all four `S x S` matrices plus the interval
#'   `prob` and draw count.
#' @return An `S x S` correlation matrix (when `summary` is given) or a list of
#'   the point estimate, posterior median, and interval bounds.
#' @export
tobs_associations <- function(object,
                              type = c("occupancy", "cover", "cross"),
                              summary = NULL) {
  assoc <- object$spatial$associations
  if (is.null(assoc))
    stop("No spatial-factor associations: 'object' is not a reduced-rank ",
         "spatial-factor community fit. Fit ms_occu_cover() with an ",
         "icar()/car_proper()/bym2() field on the occupancy arm.", call. = FALSE)
  type <- match.arg(type)
  arm  <- assoc[[type]]
  if (is.null(arm))
    stop(sprintf("No '%s' association in this fit%s.", type,
                 if (type %in% c("cover", "cross"))
                   " -- it carries no cover-arm factor" else ""),
         call. = FALSE)
  if (is.null(summary)) return(arm)
  arm[[match.arg(summary, c("median", "estimate", "lower", "upper"))]]
}


# ---------------------------------------------------------------------------
# Laplace marginal likelihood + rank (K) selection
# ---------------------------------------------------------------------------
#
# Selecting the number of latent factors K needs a criterion that INTEGRATES the
# field out, so the field prior supplies the Occam penalty. Latent-level
# pointwise criteria (held-out cells, WAIC / DIC) fail here: each extra ICAR
# field adds ~N effective latent parameters, so p_waic rises by ~N whether the
# rank-K signal is real or not -- they measure the field's effective dimension,
# not the rank.
#
# The right tool is the empirical-Bayes Laplace marginal likelihood log Z(K):
# the EM hyperparameters (Sigma, tau_w, and the field hyperparameter h for a
# car_proper / bym2 field) are the type-II MLEs, and at the converged theta the
# latent par = c(mu, {b_s}, L, W,
# log_disp) is integrated out by a Laplace approximation around the joint mode:
#
#   log Z(K) ~= logpen(mode)                                 [data LL + prior kernels]
#            + 0.5*npar*log(2pi) - 0.5*log|H|                 [Laplace volume, H = precision]
#            + 0.5*r*sum_k log tau_w[k] + 0.5*sum_k log|R(rho_k)|
#                - 0.5*K*r*log(2pi)                           [field prior NC, rank r]
#            - 0.5*nL*log(2pi*sd_L^2)                         [loading prior NC]
#            - 0.5*S*(P*log(2pi) + log|Sigma_occ|+log|Sigma_p|+log|Sigma_pos|)  [RE prior NC]
#            - 0.5*P*log(2pi) - P*log(sigma.beta)             [community-mean prior NC]
#            + sum_k log L_kk(mode)                           [log-diagonal Jacobian]
#
# The prior normalisers are not optional: a likelihood-flat direction (e.g. the
# soft field-level / intercept confound, pinned only by the RE prior) contributes
# a large posterior volume that the matching prior normaliser cancels, so only
# genuine rank-K signal moves log Z. This requires the IDENTIFIED (constrained,
# lower-triangular L) parameterisation -- the unconstrained Hessian is
# rank-deficient along the rotation manifold and |H| is then ill-defined.
#
# The field normaliser uses |R(rho_k)| at the field rank r: a full determinant at
# rank r = N for a proper-CAR field, or the generalized (pseudo) determinant of Q
# (product of its N - 1 nonzero eigenvalues) at rank r = N - 1 for an improper
# ICAR field, whose single constant null direction is counted out.

# Symmetric positive-definite log-determinant via a Cholesky ridge ladder (the
# constrained joint precision is PD up to the soft ICAR-level directions); eigen
# fallback if every ridge fails.
.ms_ocs_logdet_pd <- function(H) {
  Hs <- (H + t(H)) / 2
  for (eps in c(0, 1e-8, 1e-6, 1e-4, 1e-2)) {
    ch <- tryCatch(chol(Hs + diag(eps, nrow(Hs))), error = function(e) NULL)
    if (!is.null(ch)) return(2 * sum(log(diag(ch))))
  }
  ev <- eigen(Hs, symmetric = TRUE, only.values = TRUE)$values
  sum(log(pmax(ev, 1e-12)))
}

# Joint precision H = -Hessian(logpen) at the mode, by central finite differences
# of the ANALYTIC gradient (more accurate than differencing the objective). gfun
# returns the penalised-log-lik gradient, which vanishes at the mode, so
# H[, j] = -(grad(par + h e_j) - grad(par - h e_j)) / (2h) is the (PD) precision.
.ms_ocs_hess_fd <- function(gfun, par, h = 1e-4) {
  n <- length(par); H <- matrix(0, n, n)
  for (j in seq_len(n)) {
    pp <- par; pp[j] <- pp[j] + h
    pm <- par; pm[j] <- pm[j] - h
    H[, j] <- -(gfun(pp) - gfun(pm)) / (2 * h)
  }
  (H + t(H)) / 2
}

# Empirical-Bayes Laplace marginal log-likelihood log Z(K) for a CONVERGED,
# identifiability-constrained spatial community fit (.tobs_fit_ms_occu_cover_spatial
# with constrain = TRUE). `fit` carries the converged mode (`par`, triangular `L`),
# the type-II MLE hyperparameters (`Sigma`, `tau_w`, `sd_L`), and dims (`d`).
# logpen and H are (re)evaluated at the SAME converged theta so the mode, the
# curvature, and the prior normalisers are mutually consistent (the EM returns
# the E-step mode at theta_t alongside the M-step theta_{t+1}; at convergence they
# coincide, but recomputing removes the half-step). Returns log Z and a per-term
# breakdown for diagnostics.
.ms_ocs_log_evidence <- function(model, fit, sigma.beta = 5) {
  if (!isTRUE(fit$constrained)) {
    stop("log Z(K) requires the identified (constrained) fit; refit with ",
         "constrain = TRUE.", call. = FALSE)
  }
  d <- .ms_ocs_dims(model); K <- d$K; S <- d$S; P <- d$P; N <- d$N
  tau_w <- rep_len(fit$tau_w, K); sd_L <- fit$sd_L
  Sigma <- fit$Sigma
  spec    <- model$field_spec %||%
    .ms_ocs_field_spec(model$adj, model$field_type %||% "icar")
  hyper_w <- fit$field_hyper
  # Recompute logpen / H at the SAME converged field structure R(h).
  model$field_R <- .ms_ocs_build_field_R(model, hyper_w)

  Sinv <- matrix(0, P, P)
  Sinv[d$occ_idx, d$occ_idx] <- solve(Sigma$occ)
  Sinv[d$p_idx,   d$p_idx]   <- solve(Sigma$p)
  Sinv[d$pos_idx, d$pos_idx] <- solve(Sigma$pos)
  Pmu      <- diag(1 / sigma.beta^2, P)
  inv_sdL2 <- 1 / sd_L^2

  obj  <- function(p, grad) .ms_ocs_penll_grad_c(model, p, Sinv, Pmu, inv_sdL2,
                                                 tau_w, grad = grad)
  logpen <- obj(fit$par, FALSE)$ll
  H      <- .ms_ocs_hess_fd(function(p) obj(p, TRUE)$grad, fit$par)
  npar   <- length(fit$par)

  logdetH  <- .ms_ocs_logdet_pd(H)
  nL       <- .ms_ocs_lfree_dim(S, K)
  # The cover-arm loadings (if present) share the occupancy loadings' ridge, so
  # they add S*K free parameters to the loading-prior normaliser.
  n_load   <- nL + if (isTRUE(d$cover_factor)) d$Lpos_w else 0L
  ld_occ <- as.numeric(determinant(Sigma$occ, logarithm = TRUE)$modulus)
  ld_p   <- as.numeric(determinant(Sigma$p,   logarithm = TRUE)$modulus)
  ld_pos <- as.numeric(determinant(Sigma$pos, logarithm = TRUE)$modulus)

  Lmat <- if (K == 1L) matrix(fit$L, S, 1L) else fit$L
  diagL <- diag(Lmat)[seq_len(K)]

  vol      <- 0.5 * npar * log(2 * pi) - 0.5 * logdetH
  # Field prior normaliser: per factor 0.5 r log(tau) + 0.5 log|R| - 0.5 r log2pi,
  # with rank r = N (proper car) or N - 1 (improper icar) and the matching (full
  # or pseudo) determinant of R(rho_k). The improper null direction is counted
  # out (r = N - 1), so a flat field level contributes no spurious volume.
  r        <- spec$rank
  logdetR  <- if (isTRUE(spec$has_hyper)) {
    vapply(rep_len(hyper_w, K), function(hh) .ms_ocs_field_logdetR(spec, hh), 0)
  } else rep(.ms_ocs_field_logdetR(spec), K)
  nc_field <- 0.5 * r * sum(log(tau_w)) + 0.5 * sum(logdetR) -
              0.5 * K * r * log(2 * pi)
  nc_load  <- -0.5 * n_load * log(2 * pi * sd_L^2)
  nc_b     <- -0.5 * S * (P * log(2 * pi) + ld_occ + ld_p + ld_pos)
  nc_mu    <- -0.5 * P * log(2 * pi) - P * log(sigma.beta)
  jac      <- sum(log(pmax(diagL, 1e-12)))

  logZ <- logpen + vol + nc_field + nc_load + nc_b + nc_mu + jac
  list(logZ = as.numeric(logZ), K = K, npar = npar,
       logpen = logpen, vol = vol, logdetH = logdetH,
       nc_field = nc_field, nc_load = nc_load, nc_b = nc_b, nc_mu = nc_mu,
       jac = jac)
}

# Fit a ladder of K and pick the rank by the Laplace marginal likelihood. Each K
# is fit with the identifiability-constrained Laplace-EM (so log Z is well posed);
# the designs are K-invariant, so only model$K changes between fits. Returns a
# per-K evidence table (logZ + the breakdown terms), the selected K (argmax
# logZ), and the converged fit at the selected K (ready for the front-door
# wrapper). `K.max` defaults to a small ladder; selection stops early once log Z
# has decreased for two consecutive K (the evidence is unimodal in K once the
# signal is captured).
.ms_ocs_select_K <- function(model, K.max = 4L, sd_L = 1.0, sigma.beta = 5,
                             max.em = 40L, tol = 1e-4, verbose = FALSE) {
  K.max <- min(as.integer(K.max), model$n_species)
  rows <- list(); fits <- list(); best_drop <- 0L
  for (K in seq_len(K.max)) {
    m2 <- model; m2$K <- K
    f  <- .tobs_fit_ms_occu_cover_spatial(m2, sd_L = sd_L, max.em = max.em,
                                          tol = tol, sigma.beta = sigma.beta,
                                          constrain = TRUE)
    ev <- .ms_ocs_log_evidence(m2, f, sigma.beta = sigma.beta)
    fits[[K]] <- f
    rows[[K]] <- data.frame(K = K, logZ = ev$logZ, logpen = ev$logpen,
                            vol = ev$vol, nc_field = ev$nc_field,
                            nc_load = ev$nc_load, npar = ev$npar)
    if (verbose) {
      cat(sprintf("K=%d  logZ=%.2f  logpen=%.2f  vol=%.2f  nc_field=%.2f\n",
                  K, ev$logZ, ev$logpen, ev$vol, ev$nc_field))
    }
    if (K > 1L && ev$logZ < rows[[K - 1L]]$logZ) {
      best_drop <- best_drop + 1L
      if (best_drop >= 2L) break
    } else {
      best_drop <- 0L
    }
  }
  tab <- do.call(rbind, rows)
  K_sel <- tab$K[which.max(tab$logZ)]
  tab$best <- tab$K == K_sel
  list(table = tab, K = K_sel, fit = fits[[K_sel]])
}


# Detect a spatial request on the three occu_cover arms. The supported surface is
# a single areal field term -- icar() (improper), car_proper() (proper CAR with a
# correlation rho), or bym2() (the convolution with a variance fraction phi) -- on
# the occupancy arm and, optionally, the SAME field on the cover (positive) arm
# (a cover-arm factor, gcol33/tulpa#67 Stage 3). Returns the shared-field
# adjacency, the `field_type` ("icar" / "car_proper" / "bym2"), the
# fixed-effects occupancy / cover formulas (with the field term stripped), and
# `cover_factor` (TRUE when the cover arm also carries the field); NULL when no
# arm carries a structured term (the non-spatial path). Detection terms,
# unsupported terms, a multi-term arm, a cover-arm field without a matching
# occupancy field, or a cover field of a different type / graph all error (the
# field is shared, so the two arms must name one term on one graph).
.tobs_ms_ocs_spatial_request <- function(occ_formula, det_formula, pos_formula,
                                         data) {
  parse_terms <- function(f) {
    if (is.null(f)) return(list())
    .tobs_parse_formula(f, data = data)$terms
  }
  occ_terms <- parse_terms(occ_formula)
  det_terms <- parse_terms(det_formula)
  pos_terms <- parse_terms(pos_formula)

  if (length(det_terms)) {
    stop("ms_occu_cover(): structured terms are supported on the occupancy and ",
         "cover arms only; the detection arm must use a plain formula.",
         call. = FALSE)
  }
  if (length(occ_terms) == 0L) {
    if (length(pos_terms)) {
      stop("ms_occu_cover(): a cover-arm spatial factor shares the occupancy ",
           "field, so it requires a matching icar() term on the occupancy arm.",
           call. = FALSE)
    }
    return(NULL)                                        # non-spatial path
  }

  # Supported field terms: icar() (improper), car_proper() (proper CAR with a
  # correlation rho), and bym2() (the convolution with a variance fraction phi).
  # Validate one such term per spatial arm.
  supported <- c("icar", "car_proper", "bym2")
  one_field <- function(terms, arm) {
    if (length(terms) > 1L) {
      stop(sprintf("ms_occu_cover(): the %s arm supports a single areal field term.",
                   arm), call. = FALSE)
    }
    spec <- terms[[1L]]
    if (!inherits(spec, "tobs_spatial") || !(spec$label %in% supported)) {
      stop(sprintf(
        "ms_occu_cover(): spatial supports icar() / car_proper() / bym2() only; got %s() on the %s arm.",
        spec$label %||% class(spec)[1L], arm), call. = FALSE)
    }
    spec
  }
  occ_spec   <- one_field(occ_terms, "occupancy")
  field_type <- occ_spec$label
  graph      <- occ_spec$graph

  cover_factor <- length(pos_terms) > 0L
  if (cover_factor) {
    pos_spec <- one_field(pos_terms, "cover")
    if (!identical(pos_spec$label, field_type)) {
      stop(sprintf(paste0("ms_occu_cover(): the cover-arm factor shares the ",
                          "occupancy field, so it must use the same term (%s()) ",
                          "as the occupancy arm."), field_type), call. = FALSE)
    }
    if (!isTRUE(all.equal(unname(as.matrix(pos_spec$graph)),
                          unname(as.matrix(graph))))) {
      stop("ms_occu_cover(): the cover-arm factor shares the occupancy field, so ",
           "the areal term must name the same graph on both arms.", call. = FALSE)
    }
  }

  list(graph        = graph,
       field_type   = field_type,
       cover_factor = cover_factor,
       fe_occ = .tobs_parse_formula(occ_formula, data = data)$fe_formula,
       fe_pos = if (cover_factor) {
         .tobs_parse_formula(pos_formula, data = data)$fe_formula
       } else NULL)
}
