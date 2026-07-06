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
