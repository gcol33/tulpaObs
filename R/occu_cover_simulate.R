# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate joint occupancy-detection + cover data
#'
#' Per-cell mixture: latent z_i ~ Bernoulli(psi_i), per-visit detection
#' y_ij | z_i = 1 ~ Bernoulli(p_ij), per-visit cover y_pos_ij | y_ij = 1
#' drawn from `positive` (`"beta"` or `"lognormal"`) on the cover-arm linear
#' predictor. Used by the recovery test and as the generator for the
#' joint occupancy + cover hurdle family (see [occu_cover()]).
#'
#' When `adj` is supplied (a square adjacency matrix), an ICAR field
#' `f[1..N]` is drawn from `MVN(0, Q^-)` (with sum-to-zero constraint),
#' and the linear predictors become
#'
#'     eta_psi_i = X_psi[i, ] %*% beta_occ + sigma * f[i]
#'     eta_pos_ij = X_pos[i, ] %*% beta_pos + alpha * sigma * f[i]
#'
#' This is the shared-field parameterisation (field SD `sigma`, cover-arm copy
#' amplitude `alpha`) that the nested-Laplace [occu_cover()] fit estimates.
#'
#' @param N Number of sites (cells).
#' @param J Number of visits per site.
#' @param n.occ.covs,n.det.covs,n.pos.covs Number of covariates on each arm
#'   (drawn IID standard normal).
#' @param beta.occ,beta.p,beta.pos Coefficient vectors c(intercept, slopes).
#'   Defaults pick weakly-informative values: psi intercept at logit(0.4),
#'   p intercept at logit(0.5), cover intercept on the appropriate link.
#' @param positive `"beta"` or `"lognormal"`.
#' @param phi Beta precision when `positive = "beta"` (default 30).
#' @param sigma.pos Lognormal residual SD when `positive = "lognormal"`
#'   (default 0.4).
#' @param adj Optional N x N adjacency matrix. When supplied, generates the
#'   shared ICAR field; when NULL, the simulator is non-spatial (matches v1).
#' @param sigma Spatial field amplitude (used only when `adj` is supplied), in
#'   the geo-mean-marginal-SD (Sorbye-Rue) convention: the field is drawn with
#'   geometric-mean marginal variance 1 before scaling by `sigma`, so `sigma`
#'   is directly comparable to a fitted `field_sd` (`fit$spatial$field_sd_mean`
#'   on the `nested_laplace` path, `fit$hyper_draws[, "field_sd"]` on `nuts`) --
#'   never to the raw `sigma` a `nested_laplace` fit reports, which is the
#'   amplitude against the unscaled intrinsic precision and differs by
#'   `sqrt(scale_q)`.
#' @param alpha Cover-arm scaling on the shared field (used only when `adj`
#'   is supplied). 1.0 = arms see the field identically; positive = same sign,
#'   negative = opposite.
#' @param trend Logical; when `TRUE` (and `adj` is supplied) a SECOND shared
#'   ICAR field `f2` (a spatially-varying temporal trend) is generated on the
#'   same graph, weighted by a per-cell covariate `time` drawn IID standard
#'   normal. The trend enters the occupancy and cover predictors as
#'   `sigma.trend * time_i * f2[i]` (occupancy) and
#'   `alpha.trend * sigma.trend * time_i * f2[i]` (cover); the detection
#'   predictor is unaffected. The `time` covariate is per-cell and broadcast
#'   to every visit of that cell.
#' @param sigma.trend Trend-field amplitude (used only when `trend = TRUE`).
#' @param alpha.trend Cover-arm scaling on the trend field (used only when
#'   `trend = TRUE`).
#' @param pos.field Logical; when `TRUE` (and `adj` is supplied) draw an
#'   INDEPENDENT areal field on the cover (positive) arm only -- an intercept
#'   field plus a time-weighted trend field, each unrelated to the occupancy
#'   field and with no alpha copy. Adds a `time` column to
#'   the returned `data` and reports `g0` / `g1` (the two fields) and their SDs in
#'   `truth`. Fit by placing `spatial(~ 1 + time || cell, graph = adj)` in the
#'   `positive` formula.
#' @param sigma.pos.int Cover-arm intercept-field SD (used only when
#'   `pos.field = TRUE`).
#' @param sigma.pos.trend Cover-arm trend-field SD (used only when
#'   `pos.field = TRUE`).
#' @param det.field Logical; when `TRUE` (and `adj` is supplied) draw an
#'   INDEPENDENT areal field on the detection arm only -- an intercept field
#'   plus a time-weighted trend field, unrelated to the occupancy and cover
#'   fields and with no alpha copy, so detection varies spatially on its own.
#'   Fit by placing `spatial(~ 0 + time || cell, graph = adj, to = "detection")`
#'   in the `detection` formula.
#' @param sigma.p.int Detection-arm intercept-field SD (used only when
#'   `det.field = TRUE`).
#' @param sigma.p.trend Detection-arm trend-field SD (used only when
#'   `det.field = TRUE`).
#' @param re.det.groups Optional integer `>= 2`: the number of levels of a
#'   per-visit detection random intercept (a `habitat` factor on `visit_data`,
#'   levels `hab1..K`), drawn `b_g ~ N(0, sigma.re.p^2)` and centred. `NULL`
#'   (default) adds no detection RE. Recover it with
#'   `detection = ~ ... + (1 | habitat)`.
#' @param sigma.re.p SD of the `re.det.groups` random intercept (default 0.7).
#' @param re.pos.groups Optional integer `>= 2`: the number of levels of a
#'   per-visit COVER-arm random intercept (a `habitat` factor on `visit_data`,
#'   levels `hab1..K`), drawn `b_g ~ N(0, sigma.re.pos^2)` and centred, added to
#'   the positive-cover linear predictor. `NULL` (default) adds no cover RE.
#'   Recover it with `positive = ~ ... + (1 | habitat)` under
#'   `cover.aggregate = "none"`. Truth in `truth$b_pos_re` / `truth$sigma_re_pos`.
#' @param sigma.re.pos SD of the `re.pos.groups` random intercept (default 0.7).
#' @param re.det Optional named list of FURTHER per-visit detection random
#'   effects, for crossed / nested / slope designs. Each element
#'   `list(K =, sigma =, prefix =, nested_in =, slope_cov =, sigma_slope =, rho =)`
#'   adds a factor column (levels `<prefix>1..K`). Without `slope_cov` it is a
#'   centred `N(0, sigma^2)` random intercept; `nested_in = "<name>"` nests its
#'   codes within a previously listed grouping (matching `(1 | parent/child)`),
#'   otherwise crossed. With `slope_cov = "<column>"` (a per-visit covariate,
#'   generated `N(0, slope_sd)` if absent, `slope_sd` default 1) it is a random
#'   slope: a slope-only uncorrelated block when `rho` is unset, or a correlated
#'   intercept + slope block (covariance from `sigma` / `sigma_slope` / `rho`)
#'   when `rho` is given.
#'   Truth is returned in `truth$re_det[[name]]` (named by the level label a fit
#'   reconstructs): `b` / `b_slope` BLUP vectors, or the `B` BLUP matrix plus
#'   `s0` / `s1` / `rho` for a correlated slope.
#' @param seed Optional integer seed.
#' @return A list with `y` (N x J detection matrix), `y_pos` (N x J cover
#'   matrix, NA where not detected), `data` (per-site covariate frame, gaining
#'   a `time` column when `trend = TRUE`), `visit_data` (per-visit covariate
#'   frame, N*J rows in site-major order), and `truth` (the coefficients,
#'   dispersion, and field(s) if generated; `f2`, `sigma_trend`, `alpha_trend`,
#'   and `time` when `trend = TRUE`).
#' @examples
#' sim <- simulate_occu_cover(N = 50, J = 3, seed = 1)
#' dim(sim$y)
#' @export
simulate_occu_cover <- function(N             = 200L,
                                 J             = 4L,
                                 n.occ.covs    = 1L,
                                 n.det.covs    = 1L,
                                 n.pos.covs    = 1L,
                                 beta.occ      = NULL,
                                 beta.p        = NULL,
                                 beta.pos      = NULL,
                                 positive      = c("lognormal", "beta", "gaussian"),
                                 phi           = 30,
                                 sigma.pos     = 0.4,
                                 adj           = NULL,
                                 sigma         = 0.6,
                                 alpha         = 1.0,
                                 trend         = FALSE,
                                 sigma.trend   = 0.6,
                                 alpha.trend   = 1.0,
                                 pos.field       = FALSE,
                                 sigma.pos.int   = 0.5,
                                 sigma.pos.trend = 0.6,
                                 det.field       = FALSE,
                                 sigma.p.int     = 0.5,
                                 sigma.p.trend   = 0.6,
                                 re.det.groups = NULL,
                                 sigma.re.p    = 0.7,
                                 re.pos.groups = NULL,
                                 sigma.re.pos  = 0.7,
                                 re.det        = NULL,
                                 seed          = NULL) {
  positive <- match.arg(positive)
  if (!is.null(seed)) set.seed(seed)
  N <- as.integer(N); J <- as.integer(J)

  if (is.null(beta.occ)) beta.occ <- c(stats::qlogis(0.4), stats::runif(n.occ.covs, -0.5, 0.5))
  if (is.null(beta.p))   beta.p   <- c(0, stats::runif(n.det.covs, -0.5, 0.5))
  if (is.null(beta.pos)) {
    pos_int <- switch(positive,
                      beta     = stats::qlogis(0.3),
                      gaussian = 2.0,
                      log(0.1))
    beta.pos <- c(pos_int, stats::runif(n.pos.covs, -0.5, 0.5))
  }

  # Optional shared ICAR field(s). Draw each f as N(0, Q^-) via the
  # eigendecomposition of Q; the constant (null) component is dropped, giving a
  # zero-mean draw on the constrained space, then divide by sqrt(scale_q) so the
  # field has geo-mean marginal variance 1 (the Sorbye-Rue convention; `sigma *
  # f` then has geo-mean marginal SD sigma, matching INLA's `scale.model = TRUE`
  # and the `field_sd` the sampled-hyper NUTS route reports).
  #
  # THAT IS NOT THE UNIT THE JOINT NESTED-LAPLACE FIT REPORTS `sigma` IN. That
  # engine hands its ICAR block the RAW graph precision Q = D - W at tau = 1 and
  # carries the amplitude in the arm scale, so its `sigma` is the amplitude on a
  # field of geo-mean marginal SD sqrt(scale_q) -- a factor of 2.1 on a 30-node
  # chain, 3.0 on a 60-node one. Compare a spatial occu_cover fit against this
  # simulator on the FIELD (sd / cor of the recovered surface), not by reading
  # `sigma` off both sides.
  f  <- numeric(N)
  f2 <- numeric(N)
  g0 <- numeric(N)   # arm-specific cover intercept field
  g1 <- numeric(N)   # arm-specific cover trend field
  h0 <- numeric(N)   # arm-specific detection intercept field
  h1 <- numeric(N)   # arm-specific detection trend field
  time_cov <- numeric(N)
  if (!is.null(adj)) {
    if (!is.matrix(adj) || nrow(adj) != N || ncol(adj) != N) {
      stop("adj must be an N x N adjacency matrix.", call. = FALSE)
    }
    draw_field <- function() as.numeric(.occu_cover_draw_icar_field(adj, 1L))
    f <- draw_field()
    # A time covariate is needed by either the shared trend field or the
    # arm-specific cover trend field.
    if (isTRUE(trend) || isTRUE(pos.field) || isTRUE(det.field)) {
      time_cov <- as.numeric(scale(stats::rnorm(N)))
    }
    if (isTRUE(trend)) f2 <- draw_field()
    # Arm-specific cover field(s): an INDEPENDENT cover-arm intercept field g0
    # and time-weighted trend field g1, each unrelated to the occupancy field f.
    # They enter the cover linear predictor only (no psi contribution, no alpha
    # copy), so delta_cover_cond carries a spatial structure the occupancy
    # field's alpha copy cannot express.
    if (isTRUE(pos.field)) {
      g0 <- draw_field()
      g1 <- draw_field()
    }
    # Arm-specific detection field(s): an INDEPENDENT detection-arm intercept field
    # h0 and time-weighted trend field h1, unrelated to the occupancy and cover
    # fields. They enter the detection linear predictor only (no psi / cover
    # contribution, no copy), so p varies spatially on its own.
    if (isTRUE(det.field)) {
      h0 <- draw_field()
      h1 <- draw_field()
    }
  }

  # Site-level covariates (psi predictor).
  occ_covs <- data.frame(matrix(stats::rnorm(N * n.occ.covs), N, n.occ.covs))
  names(occ_covs) <- paste0("occ_cov", seq_len(n.occ.covs))
  X_occ <- stats::model.matrix(~ ., occ_covs)
  eta_psi <- as.vector(X_occ %*% beta.occ) + sigma * f
  if (!is.null(adj) && isTRUE(trend)) {
    eta_psi <- eta_psi + sigma.trend * time_cov * f2
  }
  psi <- stats::plogis(eta_psi)
  z_state <- stats::rbinom(N, 1L, psi)

  # Visit-level covariates (p and cover predictors). Same draw used for both
  # arms, mirroring how `tobs_data()`'s `det.covs` matrices feed both formulas.
  det_covs <- data.frame(matrix(stats::rnorm(N * J * n.det.covs), N * J, n.det.covs))
  names(det_covs) <- paste0("det_cov", seq_len(n.det.covs))
  pos_covs <- data.frame(matrix(stats::rnorm(N * J * n.pos.covs), N * J, n.pos.covs))
  names(pos_covs) <- paste0("pos_cov", seq_len(n.pos.covs))
  visit_data <- cbind(det_covs, pos_covs)

  X_p   <- stats::model.matrix(~ ., det_covs)
  X_pos <- stats::model.matrix(~ ., pos_covs)
  eta_p   <- as.vector(X_p   %*% beta.p)
  eta_pos_base <- as.vector(X_pos %*% beta.pos)

  # Optional per-visit detection random effects. Each grouping is a categorical
  # visit-level factor (e.g. an EUNIS habitat class) with a random intercept b_g ~
  # N(0, sigma^2) on the detection linear predictor; the factor rides
  # `visit_data`, so a fit reads it via `visits` and `detection = ~ ... + (1 |
  # <factor>)`. `re_det_groups` / `sigma_re_p` set the first ("habitat") grouping
  # (back-compat); `re_det` is a named list of further groupings, each `list(K =,
  # sigma =, prefix =, nested_in =)`, for CROSSED (`nested_in = NULL`) or NESTED
  # (`nested_in = "<parent>"`, sub-codes nested within the parent factor's codes
  # -- matching `(1 | parent/child)`) designs. BLUPs are centred so the detection
  # intercept stays identified. Truth BLUPs are stored NAMED by the level label a
  # fit reconstructs (the interaction label for a nested grouping), so recovery
  # checks align by name, not factor sort order.
  grp_specs <- list()
  if (!is.null(re.det.groups)) {
    grp_specs[["habitat"]] <- list(K = as.integer(re.det.groups),
                                   sigma = sigma.re.p, prefix = "hab",
                                   nested_in = NULL)
  }
  if (!is.null(re.det)) {
    for (nm in names(re.det)) {
      s <- re.det[[nm]]
      grp_specs[[nm]] <- list(K = as.integer(s$K), sigma = s$sigma %||% 0.7,
                              prefix = s$prefix %||% nm, nested_in = s$nested_in,
                              # Random-slope fields: `slope_cov` names a per-visit
                              # covariate column (generated N(0, 1) if absent). With
                              # `rho` set it is a correlated intercept + slope block
                              # (Sigma from sigma / sigma_slope / rho); without, a
                              # slope-only uncorrelated block.
                              slope_cov = s$slope_cov,
                              sigma_slope = s$sigma_slope,
                              rho = s$rho, slope_sd = s$slope_sd)
    }
  }
  re_codes  <- list()    # within-grouping per-visit code (for nesting)
  re_truth  <- list()
  b_p_re <- NULL; re_det_levels <- NULL
  for (nm in names(grp_specs)) {
    s <- grp_specs[[nm]]
    if (s$K < 2L)
      stop(sprintf("re_det grouping '%s' needs K >= 2.", nm), call. = FALSE)
    sub <- sample.int(s$K, N * J, replace = TRUE)
    re_codes[[nm]] <- sub
    visit_data[[nm]] <- factor(paste0(s$prefix, sub),
                               levels = paste0(s$prefix, seq_len(s$K)))
    if (is.null(s$nested_in)) {
      code   <- sub
      labels <- paste0(s$prefix, seq_len(s$K))
    } else {
      parent <- grp_specs[[s$nested_in]]
      pcode  <- re_codes[[s$nested_in]]
      code   <- (pcode - 1L) * s$K + sub                 # unique (parent, sub)
      pc     <- rep(seq_len(parent$K), each = s$K)
      sc     <- rep(seq_len(s$K),      times = parent$K)
      labels <- paste0(parent$prefix, pc, ".", s$prefix, sc)  # interaction label
    }
    if (is.null(s$slope_cov)) {
      # Random intercept.
      b <- stats::rnorm(length(labels), 0, s$sigma); b <- b - mean(b)
      names(b) <- labels
      eta_p <- eta_p + b[code]
      re_truth[[nm]] <- list(b = b, sigma = s$sigma, levels = labels,
                             kind = "intercept")
      if (identical(nm, "habitat")) {                    # back-compat truth slots
        b_p_re <- unname(b); re_det_levels <- labels
      }
    } else {
      # Random slope on a per-visit covariate (generated N(0, slope_sd) if
      # absent; `slope_sd` != 1 exercises the slope-covariate standardization).
      if (is.null(visit_data[[s$slope_cov]]))
        visit_data[[s$slope_cov]] <- stats::rnorm(N * J, 0, s$slope_sd %||% 1)
      xv <- as.numeric(visit_data[[s$slope_cov]])
      if (is.null(s$rho)) {
        # Slope-only uncorrelated block (0 + x | g).
        b1 <- stats::rnorm(length(labels), 0, s$sigma); b1 <- b1 - mean(b1)
        names(b1) <- labels
        eta_p <- eta_p + xv * b1[code]
        re_truth[[nm]] <- list(b_slope = b1, sigma = s$sigma, levels = labels,
                               cov = s$slope_cov, kind = "slope")
      } else {
        # Correlated intercept + slope block (1 + x | g): (b0, b1) ~ N(0, Sigma).
        s0 <- s$sigma; s1 <- s$sigma_slope %||% s$sigma; rho <- s$rho
        Sig <- matrix(c(s0^2, rho * s0 * s1, rho * s0 * s1, s1^2), 2L, 2L)
        L   <- t(chol(Sig))
        B2  <- matrix(stats::rnorm(2L * length(labels)), length(labels), 2L) %*% t(L)
        B2  <- sweep(B2, 2L, colMeans(B2))               # centre each coef
        rownames(B2) <- labels; colnames(B2) <- c("(Intercept)", s$slope_cov)
        eta_p <- eta_p + B2[code, 1L] + xv * B2[code, 2L]
        re_truth[[nm]] <- list(B = B2, s0 = s0, s1 = s1, rho = rho,
                               levels = labels, cov = s$slope_cov, kind = "corr")
      }
    }
  }

  # Optional per-visit COVER-arm random intercept. A categorical visit-level
  # grouping carries a random intercept on the positive- cover linear predictor,
  # mirroring the detection-arm RE above; it rides `visit_data` so a fit reads
  # it via `positive = ~ ... + (1 | habitat)` under cover_aggregate = "none".
  # BLUPs are centred so the cover intercept stays identified. Kept separate
  # from the detection block: a pos-arm recovery test sets re_pos_groups alone,
  # so the `habitat` factor carries a genuine cover offset (the detection RE
  # would otherwise put the offset on eta_p).
  b_pos_re <- NULL; re_pos_levels <- NULL
  if (!is.null(re.pos.groups)) {
    Kp <- as.integer(re.pos.groups)
    if (Kp < 2L) stop("re_pos_groups needs K >= 2.", call. = FALSE)
    sub_pos    <- sample.int(Kp, N * J, replace = TRUE)
    re_pos_levels <- paste0("hab", seq_len(Kp))
    visit_data[["habitat"]] <- factor(paste0("hab", sub_pos),
                                      levels = re_pos_levels)
    b_pos <- stats::rnorm(Kp, 0, sigma.re.pos); b_pos <- b_pos - mean(b_pos)
    names(b_pos)  <- re_pos_levels
    eta_pos_base  <- eta_pos_base + b_pos[sub_pos]
    b_pos_re      <- b_pos
  }

  y     <- matrix(0L, N, J)
  y_pos <- matrix(NA_real_, N, J)

  for (i in seq_len(N)) {
    for (j in seq_len(J)) {
      idx <- (i - 1L) * J + j
      if (z_state[i] == 1L) {
        eta_p_ij <- eta_p[idx]
        if (!is.null(adj) && isTRUE(det.field)) {
          eta_p_ij <- eta_p_ij + sigma.p.int * h0[i] +
                      sigma.p.trend * time_cov[i] * h1[i]
        }
        p_ij <- stats::plogis(eta_p_ij)
        d <- stats::rbinom(1L, 1L, p_ij)
        y[i, j] <- d
        if (d == 1L) {
          eta_pos_ij <- eta_pos_base[idx] + alpha * sigma * f[i]
          if (!is.null(adj) && isTRUE(trend)) {
            eta_pos_ij <- eta_pos_ij + alpha.trend * sigma.trend * time_cov[i] * f2[i]
          }
          if (!is.null(adj) && isTRUE(pos.field)) {
            eta_pos_ij <- eta_pos_ij +
              sigma.pos.int * g0[i] + sigma.pos.trend * time_cov[i] * g1[i]
          }
          if (positive == "beta") {
            mu <- stats::plogis(eta_pos_ij)
            y_pos[i, j] <- stats::rbeta(1L, mu * phi, (1 - mu) * phi)
          } else if (positive == "gaussian") {
            # Identity-Gaussian arm: the cover magnitude is a plain Gaussian on
            # the raw response (no log transform).
            y_pos[i, j] <- stats::rnorm(1L, eta_pos_ij, sigma.pos)
          } else {
            y_pos[i, j] <- exp(stats::rnorm(1L, eta_pos_ij, sigma.pos))
          }
        }
      }
    }
  }

  occ_out <- occ_covs
  has_trend    <- !is.null(adj) && isTRUE(trend)
  has_posfield <- !is.null(adj) && isTRUE(pos.field)
  has_detfield <- !is.null(adj) && isTRUE(det.field)
  if (has_trend || has_posfield || has_detfield) occ_out$time <- time_cov
  # The arm-specific field bars index their graph node by a `cell` column.
  if (has_posfield || has_detfield) occ_out$cell <- seq_len(N)

  list(
    y          = y,
    y_pos      = y_pos,
    data       = occ_out,
    visit_data = visit_data,
    adj        = adj,
    truth      = list(
      beta_occ    = beta.occ,
      beta_p      = beta.p,
      beta_pos    = beta.pos,
      psi         = psi,
      z           = z_state,
      positive    = positive,
      phi         = if (positive == "beta")      phi       else NA_real_,
      sigma_pos   = if (positive == "lognormal") sigma.pos else NA_real_,
      f           = f,
      sigma       = if (!is.null(adj)) sigma else NA_real_,
      alpha       = if (!is.null(adj)) alpha else NA_real_,
      f2          = if (has_trend) f2          else NULL,
      time        = if (has_trend || has_posfield) time_cov else NULL,
      sigma_trend = if (has_trend) sigma.trend else NA_real_,
      alpha_trend = if (has_trend) alpha.trend else NA_real_,
      g0              = if (has_posfield) g0              else NULL,
      g1              = if (has_posfield) g1              else NULL,
      sigma_pos_int   = if (has_posfield) sigma.pos.int   else NA_real_,
      sigma_pos_trend = if (has_posfield) sigma.pos.trend else NA_real_,
      h0              = if (has_detfield) h0              else NULL,
      h1              = if (has_detfield) h1              else NULL,
      sigma_p_int     = if (has_detfield) sigma.p.int     else NA_real_,
      sigma_p_trend   = if (has_detfield) sigma.p.trend   else NA_real_,
      sigma_re_p  = if (!is.null(re.det.groups)) sigma.re.p else NA_real_,
      b_p_re      = b_p_re,
      re_det_levels = re_det_levels,
      sigma_re_pos = if (!is.null(re.pos.groups)) sigma.re.pos else NA_real_,
      b_pos_re     = b_pos_re,
      re_pos_levels = re_pos_levels,
      re_det      = if (length(re_truth)) re_truth else NULL
    )
  )
}
