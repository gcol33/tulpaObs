# ============================================================================
# Shared state-arm encoding
#
# The occupancy state arm is encoded for the M-step as a pseudo-binomial:
# y = round(M * w), n_trials = M, where w = E[z_i | y] from the E-step. M is the
# pseudo-trial count that encoding pretends each site carries, and all three
# families (single / dynamic / integrated) pick it by the same rule:
#
#   latent_prior  -> M = 1. A nested-Laplace latent block. A state row carries
#     exactly ONE binary occupancy observation, so M = 1 is the site's real
#     information content. Anything larger overstates it M-fold and swamps the
#     field prior: between-cell binomial noise reads as a real field, the field
#     inflates, and the state slope inflates with it through the logistic
#     conditional-vs-marginal factor sqrt(1 + 0.346 sigma^2).
#   spatial_occ   -> M = 4. A single-Laplace SPDE mesh field, which needs the
#     same protection; M = 4 keeps the per-site effective sample size O(1) while
#     leaving some fractional resolution on the weights.
#   neither       -> M = 1000. The inflation makes the M-step a sharp binomial
#     whose mode equals the weighted mean, and there is no prior for it to swamp.
#
# The two arguments are mutually exclusive by construction --
# .tobs_laplace_nested() always passes spatial = NULL -- so the order of the
# first two branches is not reachable-input-sensitive.
#
# Measured on a 40-cell chain, 6 sites/cell, 12-20 seeds (dev_notes/_run_m_final.R,
# _run_m_single.R), truth slope 0.5 / f0 sd 1.0 / f1 sd 0.8:
#
#                      slope    coverage   f0 sd    f1 sd
#   dyn + icar  M=4    0.6232   0.83       1.3750   -
#   dyn + icar  M=1    0.5253   1.00       0.8618   -
#   dyn + SVC   M=4    0.7413   0.58       1.8212   1.8877
#   dyn + SVC   M=1    0.4872   0.92       0.9231   0.8601
#   single      M=4    0.5669   0.95       0.9295   -
#   single      M=1    0.5207   0.95       0.9294   -
#
# Monotone in M in every arm, and M = 1 is uniformly best: no arm regresses. The
# apparent downside -- round(w * 1) is 0/1, so the fractional resolution the
# inflation exists for is lost -- does not materialise: the slope IMPROVES,
# because the information-content error dominates the rounding loss.
#
# An areal icar/bym2/car_proper block is strongly informative at the grid scale,
# which does NOT let it tolerate a sharp encoding: the rows it is being fit
# against still carry one bit each, and that is what M has to match.
# ============================================================================
.tobs_state_M <- function(spatial_occ, latent_prior) {
  if (!is.null(latent_prior)) 1L
  else if (!is.null(spatial_occ)) 4L
  else 1000L
}

# Pseudo-binomial state block at the M above. `any_det` sites are pinned to w = 1
# (a detection proves occupancy). The SPDE mesh projection attaches only on the
# single-Laplace spatial path; a nested latent block's prior is attached to
# occ$prior upstream by .tobs_laplace_nested().
.tobs_encode_state_block <- function(weights, any_det, X_occ, spatial_occ,
                                     latent_prior) {
  M <- .tobs_state_M(spatial_occ, latent_prior)
  y_occ <- ifelse(any_det, M, as.integer(round(weights * M)))
  y_occ <- pmin(pmax(y_occ, 0L), M)
  occ_block <- list(y = y_occ, n_trials = rep(M, length(y_occ)), X = X_occ,
                    family = "binomial")
  # A no-op unless spatial_occ is an SPDE term; a nested latent block never
  # reaches here with a non-NULL spatial_occ.
  .attach_spatial_spde(occ_block, spatial_occ)
}


# ============================================================================
# Single-season callbacks
# ============================================================================
build_single_callbacks <- function(model, spatial = NULL, latent_prior = NULL) {
  y <- model$y
  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  X_det_visit <- model$X_det_visit  # NULL when no visit-level covariates
  max_visits <- ncol(y)
  n_sites <- model$n_sites
  p_occ <- ncol(X_occ)
  p_det <- ncol(X_det)
  p_det_visit <- if (is.null(X_det_visit)) 0L else ncol(X_det_visit)
  p_det_total <- p_det + p_det_visit

  # An SPDE term may enter the state arm, the detection arm, or both (its
  # `$shared = c(occ, det)` membership). Resolve the per-arm field once: the
  # state field attaches to the occ block, the detection field to the det
  # block. A detection field with visit-level detection covariates is not yet
  # wired (the field is site-indexed; the det block is per (site, visit) and
  # would need a row-expanded mesh projection).
  spatial_occ <- .spatial_for_arm(spatial, 1L)
  spatial_det <- .spatial_for_arm(spatial, 2L)
  if (!is.null(spatial_det) && p_det_visit > 0L) {
    stop("SPDE on the detection process with visit-level detection covariates ",
         "is not yet plumbed in .tobs_laplace; use shared occupancy-arm SPDE ",
         "or method = 'nuts'.", call. = FALSE)
  }

  n_valid <- integer(n_sites)
  n_det <- integer(n_sites)
  any_det <- logical(n_sites)
  valid_mat <- matrix(FALSE, n_sites, max_visits)
  for (i in seq_len(n_sites)) {
    v <- y[i, ] >= 0
    valid_mat[i, ] <- v
    n_valid[i] <- sum(v)
    n_det[i] <- sum(y[i, v] == 1)
    any_det[i] <- n_det[i] > 0
  }
  keep <- n_valid > 0

  # Per-visit indexing for the X_det_visit path. Row r of X_det_visit
  # corresponds to (site = (r-1) %/% max_visits + 1, visit = (r-1) %% max_visits + 1).
  if (p_det_visit > 0L) {
    site_idx_all <- rep(seq_len(n_sites), each = max_visits)
    visit_idx_all <- rep(seq_len(max_visits), times = n_sites)
  }

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    eta_occ <- as.vector(X_occ %*% beta_occ)
    sp_off <- .spatial_eta_offset(spatial_occ, fits$occ, p_occ)
    if (length(sp_off) == n_sites) eta_occ <- eta_occ + sp_off
    # Nested-Laplace: make the E-step weight P(z_i = 1 | y_i) field-aware.
    eta_occ <- .nested_state_eta(eta_occ, latent_prior, fits$occ, p_occ, n_sites)
    psi <- plogis(eta_occ)

    if (p_det_visit == 0L) {
      beta_det <- extract_beta(fits$det, p_det)
      eta_det <- as.vector(X_det %*% beta_det)
      det_off <- .spatial_eta_offset(spatial_det, fits$det, p_det)
      if (length(det_off) == n_sites) eta_det <- eta_det + det_off
      p_site <- plogis(eta_det)
      return(list(weights = occ_weights(psi, p_site, n_sites,
                                        n_valid, n_det, any_det)))
    }

    # Visit-level path: logit(p_ij) = X_det[i,] beta_site + X_det_visit[(i-1)*J + j,] beta_visit
    beta_det <- extract_beta(fits$det, p_det_total)
    eta_site <- as.vector(X_det %*% beta_det[seq_len(p_det)])
    eta_visit_long <- as.vector(X_det_visit %*%
                                  beta_det[(p_det + 1L):p_det_total])
    # eta_visit_long is in site-major order: reshape so [i, j] = visit (i, j)
    eta_visit_mat <- matrix(eta_visit_long, n_sites, max_visits, byrow = TRUE)
    logit_p_ij <- matrix(eta_site, n_sites, max_visits) + eta_visit_mat
    logit_p_ij <- .tobs_clamp_eta(logit_p_ij)
    # log(1 - plogis(eta)) = -log1pexp(eta) computed stably as -pmax(eta,0) - log1p(exp(-|eta|))
    log_1mp <- -(pmax(logit_p_ij, 0) + log1p(exp(-abs(logit_p_ij))))
    log_1mp[!valid_mat] <- 0
    log_prod_1mp <- rowSums(log_1mp)
    weights <- numeric(n_sites)
    for (i in seq_len(n_sites)) {
      if (any_det[i]) {
        weights[i] <- 1
      } else if (n_valid[i] == 0L) {
        weights[i] <- psi[i]
      } else {
        num <- psi[i] * exp(log_prod_1mp[i])
        weights[i] <- num / (num + (1 - psi[i]))
      }
    }
    list(weights = weights)
  }

  m_step_encode <- function(weights, ...) {
    occ_block <- .tobs_encode_state_block(weights, any_det, X_occ, spatial_occ,
                                          latent_prior)
    # Detection block: weight by w_i = P(z_i = 1 | y_i, theta). Sites that
    # the E-step thinks are likely empty (w_i ~ 0) must drop out of the
    # detection fit, otherwise they bias p_hat downward by feeding their
    # all-zero detection history as evidence about (1 - p)^J. Sites with
    # any detection have w_i = 1 (the E-step sets this).
    w_det <- weights
    w_det[any_det] <- 1

    if (!is.null(spatial_det)) {
      # SPDE detection field: the single-Laplace spatial solver consumes no
      # per-observation `weights`, so the occupancy weight is folded into the
      # binomial response by scaling both successes and trials by w_i
      # (y = round(w_i n_det_i), n = round(w_i n_valid_i)). This is the
      # frequency-weight-as-counts identity for a binomial mode/Hessian, and it
      # keeps ALL n_sites rows so the per-site rows stay aligned with the full
      # mesh projection A (n_sites x n_mesh). A near-empty site (w_i ~ 0)
      # collapses to a (0, 0) row that contributes nothing to the likelihood,
      # score, or Hessian -- the analogue of dropping it under the explicit
      # weight on the non-spatial path.
      y_det_w <- as.integer(round(w_det * n_det))
      n_det_w <- as.integer(round(w_det * n_valid))
      y_det_w <- pmin(pmax(y_det_w, 0L), n_det_w)
      det_block <- list(y = y_det_w, n_trials = n_det_w, X = X_det,
                        family = "binomial")
      det_block <- .attach_spatial_spde(det_block, spatial_det)
    } else if (p_det_visit == 0L) {
      keep_det <- keep & (w_det > 1e-6)
      det_block <- list(y = n_det[keep_det], n_trials = n_valid[keep_det],
                        X = X_det[keep_det, , drop = FALSE],
                        weights = w_det[keep_det], family = "binomial")
    } else {
      # Per-visit detection block: one Bernoulli row per (site, valid visit)
      # whose weight is the site's posterior occupancy w_i. Combined design
      # matrix stacks site-level X_det (replicated across visits) and
      # visit-level X_det_visit (already in site-major order).
      keep_visit <- valid_mat & (w_det >= 1e-6)
      site_kept <- site_idx_all[as.vector(t(keep_visit))]
      visit_kept <- visit_idx_all[as.vector(t(keep_visit))]
      n_kept <- length(site_kept)
      if (n_kept > 0L) {
        visit_row_idx <- (site_kept - 1L) * max_visits + visit_kept
        X_combined <- cbind(
          X_det[site_kept, , drop = FALSE],
          X_det_visit[visit_row_idx, , drop = FALSE]
        )
        y_kept <- y[cbind(site_kept, visit_kept)]
        det_block <- list(y = as.integer(y_kept),
                          n_trials = rep(1L, n_kept),
                          X = X_combined,
                          weights = w_det[site_kept],
                          family = "binomial")
      } else {
        det_block <- list(y = integer(0),
                          n_trials = integer(0),
                          X = matrix(0, 0, p_det_total),
                          weights = numeric(0),
                          family = "binomial")
      }
    }
    list(occ = occ_block, det = det_block)
  }

  z_draw <- function(weights, ...) {
    z <- as.integer(any_det)
    z[!any_det] <- rbinom(sum(!any_det), 1, clamp_w(weights[!any_det]))
    z
  }

  hard_encode <- function(z, ...) {
    occ_sites <- which(z == 1L)
    det_keep <- occ_sites[n_valid[occ_sites] > 0]
    occ_block <- list(y = z, n_trials = rep(1L, n_sites), X = X_occ,
                      family = "binomial")
    occ_block <- .attach_spatial_spde(occ_block, spatial_occ)
    if (!is.null(spatial_det)) {
      # Hard-z detection field: keep ALL n_sites rows aligned with the mesh
      # projection by zeroing the trials of sites that contribute no detection
      # evidence (z = 0 or no valid visits); a (0, 0) row drops out cleanly.
      keep_det_mask <- (z == 1L) & (n_valid > 0L)
      n_det_h <- ifelse(keep_det_mask, n_valid, 0L)
      y_det_h <- ifelse(keep_det_mask, n_det, 0L)
      det_block <- list(y = as.integer(y_det_h), n_trials = as.integer(n_det_h),
                        X = X_det, family = "binomial")
      det_block <- .attach_spatial_spde(det_block, spatial_det)
    } else if (p_det_visit == 0L) {
      det_block <- if (length(det_keep) > 0)
        list(y = n_det[det_keep], n_trials = n_valid[det_keep],
             X = X_det[det_keep, , drop = FALSE], family = "binomial")
      else NULL
    } else {
      keep_visit <- valid_mat & matrix(z == 1L, n_sites, max_visits)
      site_kept <- site_idx_all[as.vector(t(keep_visit))]
      visit_kept <- visit_idx_all[as.vector(t(keep_visit))]
      n_kept <- length(site_kept)
      det_block <- if (n_kept > 0L) {
        visit_row_idx <- (site_kept - 1L) * max_visits + visit_kept
        X_combined <- cbind(
          X_det[site_kept, , drop = FALSE],
          X_det_visit[visit_row_idx, , drop = FALSE]
        )
        y_kept <- y[cbind(site_kept, visit_kept)]
        list(y = as.integer(y_kept), n_trials = rep(1L, n_kept),
             X = X_combined, family = "binomial")
      } else NULL
    }
    list(occ = occ_block, det = det_block)
  }

  init <- glm_init(X_occ, X_det, any_det, n_det, n_valid, keep, p_occ, p_det)
  if (p_det_visit > 0L) {
    # Pad det init with zeros for visit-level cols so warm-start shapes
    # match the combined design when the penalized driver pulls beta_init
    # from the previous EM iteration.
    init$det$beta <- c(init$det$beta, rep(0, p_det_visit))
    init$det$se   <- c(init$det$se,   rep(1, p_det_visit))
  }

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init,
       p_per_submodel = c(occ = p_occ, det = p_det_total))
}

# ============================================================================
# Dynamic occupancy callbacks
# ============================================================================
build_dynamic_callbacks <- function(model, spatial = NULL, latent_prior = NULL) {
  y_flat <- model$y_flat
  n_sites <- model$n_sites
  n_seasons <- model$n_seasons
  max_visits <- model$max_visits
  X_occ <- model$X_processes[[1]]  # psi1
  X_det <- model$X_processes[[2]]  # p
  X_col <- model$X_processes[[3]]  # gamma  (site-level, OR long-form when SV)
  X_ext <- model$X_processes[[4]]  # epsilon
  p_occ <- ncol(X_occ); p_det <- ncol(X_det)
  p_col <- ncol(X_col); p_ext <- ncol(X_ext)

  # Season-varying colonization / extinction (gcol33/tulpaObs#124): a rate that
  # varies by interval carries a long-form [(n_sites x (T-1)) x p] design (built
  # by .tobs_interval_arm_design, site-major interval-minor). The E-step then
  # uses a per-interval transition matrix and the M-step encodes one logistic
  # row per (site, interval); constant rates keep the site-level design, the
  # per-interval transition collapses to a constant one, and the M-step
  # aggregates over a site's intervals -- byte-identical to the pre-#124 path.
  n_int  <- model$n_intervals %||% (n_seasons - 1L)
  col_sv <- isTRUE(model$col_season_varying)
  ext_sv <- isTRUE(model$ext_season_varying)
  # Season-varying detection (gcol33/tulpaObs#124): X_det is then the long-form
  # [(n_sites x n_seasons) x p] design (site-major season-minor); the per-season
  # detection probability is a [n_sites x n_seasons] matrix. A constant-detection
  # arm keeps the site-level design and the matrix broadcasts the site's single
  # probability across seasons -- byte-identical to the pre-#124 path.
  det_sv <- isTRUE(model$det_season_varying)

  # The state field enters season-1 occupancy psi1 only (one psi1 row per
  # site, the identity map). The colonization (gamma) and extinction
  # (epsilon) transition predictors are separate latent processes whose own
  # mesh fields are not wired here; `.validate_spatial_laplace` only routes a
  # state-arm (shared[1]) SPDE term to the dynamic path, and the term
  # constructor maps it to the psi1 predictor.
  spatial_occ <- .spatial_for_arm(spatial, 1L)

  # Precompute per site-season
  nv <- model$n_visits
  ad <- model$any_detected

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    beta_det <- extract_beta(fits$det, p_det)
    beta_col <- extract_beta(fits$col, p_col)
    beta_ext <- extract_beta(fits$ext, p_ext)

    eta_psi1 <- as.vector(X_occ %*% beta_occ)
    sp_off <- .spatial_eta_offset(spatial_occ, fits$occ, p_occ)
    if (length(sp_off) == n_sites) eta_psi1 <- eta_psi1 + sp_off
    psi1 <- plogis(eta_psi1)
    # Per-season detection probability [n_sites x n_seasons]: season-varying reads
    # the long-form eta site-major season-minor; constant broadcasts the per-site p.
    p_mat <- if (det_sv)
      matrix(plogis(as.vector(X_det %*% beta_det)), n_sites, n_seasons, byrow = TRUE)
    else
      matrix(plogis(as.vector(X_det %*% beta_det)), n_sites, n_seasons)
    # gamma / epsilon as [n_sites x n_int] per-interval rate matrices. Constant
    # rates recycle a length-n_sites vector down every column (each row constant),
    # so the per-interval transition below is byte-identical to the scalar path;
    # a season-varying arm reads its long-form eta site-major interval-minor.
    if (col_sv) {
      gam_mat <- matrix(plogis(as.vector(X_col %*% beta_col)), n_sites, n_int,
                        byrow = TRUE)
    } else {
      gam_mat <- matrix(plogis(as.vector(X_col %*% beta_col)), n_sites, n_int)
    }
    if (ext_sv) {
      eps_mat <- matrix(plogis(as.vector(X_ext %*% beta_ext)), n_sites, n_int,
                        byrow = TRUE)
    } else {
      eps_mat <- matrix(plogis(as.vector(X_ext %*% beta_ext)), n_sites, n_int)
    }

    # Exact forward-backward (Baum-Welch) E-step over the 2-state occupancy chain
    # of each site. The smoothed marginal gamma_t(z) =
    # P(z_it = z | y_{1:T}) feeds the psi1 / detection sufficient statistics, and
    # the smoothed pairwise joint xi_t(z, z') = P(z_it = z, z_{i,t+1} = z' | y_{1:T})
    # feeds colonization (xi_t(0, 1)) and extinction (xi_t(1, 0)) -- the exact
    # transition sufficient statistics, replacing the earlier forward-FILTERED
    # marginal-PRODUCT approximation (1 - w_{t-1}) w_t, which converged to a biased
    # non-ML fixed point. The reduced occupancy emission (a detected season forces
    # z = 1) drops only a z-independent detection-likelihood factor that cancels in
    # every smoothed ratio, so gamma and xi are exact. The detection-rate factors
    # are fit separately from the per-visit counts in the M-step.
    A00 <- 1 - gam_mat; A01 <- gam_mat                       # transition rows,
    A10 <- eps_mat;     A11 <- 1 - eps_mat                   #   [n_sites x n_int]
    w <- matrix(NA_real_, n_sites, n_seasons)                # smoothed P(z = 1 | y)
    col_y <- numeric(n_sites); col_n <- numeric(n_sites)
    ext_y <- numeric(n_sites); ext_n <- numeric(n_sites)
    # Per-interval expected transition counts, kept only when an arm is
    # season-varying (the M-step encodes one logistic row per (site, interval)).
    col_y_mat <- if (col_sv) matrix(0, n_sites, n_int) else NULL
    col_n_mat <- if (col_sv) matrix(0, n_sites, n_int) else NULL
    ext_y_mat <- if (ext_sv) matrix(0, n_sites, n_int) else NULL
    ext_n_mat <- if (ext_sv) matrix(0, n_sites, n_int) else NULL
    a <- matrix(0, n_seasons, 2)                             # scaled forward ahat_t
    b0v <- numeric(n_seasons); b1v <- numeric(n_seasons)     # emission b_t(0), b_t(1)
    cs <- numeric(n_seasons)                                 # per-season normaliser
    for (i in seq_len(n_sites)) {
      for (t in seq_len(n_seasons)) {
        idx <- (i - 1) * n_seasons + t
        nv_it <- nv[idx]; det_it <- ad[idx]
        if (nv_it == 0L)      { b0v[t] <- 1;  b1v[t] <- 1 }
        else if (det_it)      { b0v[t] <- 0;  b1v[t] <- 1 }
        else                  { b0v[t] <- 1;  b1v[t] <- (1 - p_mat[i, t])^nv_it }
      }
      # forward (scaled): ahat_t(z) = b_t(z) sum_z' ahat_{t-1}(z') A(z',z) / c_t.
      # The step t-1 -> t uses interval (t - 1)'s transition rates.
      u0 <- (1 - psi1[i]) * b0v[1]; u1 <- psi1[i] * b1v[1]
      c1 <- u0 + u1; cs[1] <- c1; a[1, 1] <- u0 / c1; a[1, 2] <- u1 / c1
      for (t in 2:n_seasons) {
        iv <- t - 1L
        pr0 <- a[t - 1, 1] * A00[i, iv] + a[t - 1, 2] * A10[i, iv]
        pr1 <- a[t - 1, 1] * A01[i, iv] + a[t - 1, 2] * A11[i, iv]
        v0 <- b0v[t] * pr0; v1 <- b1v[t] * pr1
        ct <- v0 + v1; cs[t] <- ct; a[t, 1] <- v0 / ct; a[t, 2] <- v1 / ct
      }
      # backward (scaled) + smoothed marginals / pairwise joints, T-1 .. 1. The
      # joint at backward step t is over seasons (t, t+1) = interval t.
      bw0 <- 1; bw1 <- 1                                     # beta_T(z) = 1
      w[i, n_seasons] <- a[n_seasons, 2]                     # gamma_T(1) = ahat_T(1)
      for (t in (n_seasons - 1):1) {
        bb0 <- b0v[t + 1] * bw0; bb1 <- b1v[t + 1] * bw1
        inv_c <- 1 / cs[t + 1]
        xi01 <- a[t, 1] * A01[i, t] * bb1 * inv_c            # colonization event
        xi00 <- a[t, 1] * A00[i, t] * bb0 * inv_c
        xi10 <- a[t, 2] * A10[i, t] * bb0 * inv_c            # extinction event
        xi11 <- a[t, 2] * A11[i, t] * bb1 * inv_c
        col_y[i] <- col_y[i] + xi01; col_n[i] <- col_n[i] + (xi00 + xi01)
        ext_y[i] <- ext_y[i] + xi10; ext_n[i] <- ext_n[i] + (xi10 + xi11)
        if (col_sv) { col_y_mat[i, t] <- xi01; col_n_mat[i, t] <- xi00 + xi01 }
        if (ext_sv) { ext_y_mat[i, t] <- xi10; ext_n_mat[i, t] <- xi10 + xi11 }
        bw0 <- (A00[i, t] * bb0 + A01[i, t] * bb1) * inv_c   # beta_t(0)
        bw1 <- (A10[i, t] * bb0 + A11[i, t] * bb1) * inv_c   # beta_t(1)
        w[i, t] <- a[t, 2] * bw1                             # gamma_t(1)
      }
    }
    attr(w, "col_y") <- col_y; attr(w, "col_n") <- col_n
    attr(w, "ext_y") <- ext_y; attr(w, "ext_n") <- ext_n
    attr(w, "col_y_mat") <- col_y_mat; attr(w, "col_n_mat") <- col_n_mat
    attr(w, "ext_y_mat") <- ext_y_mat; attr(w, "ext_n_mat") <- ext_n_mat
    list(weights = w)
  }

  m_step_encode <- function(weights, ...) {
    w <- weights  # n_sites x n_seasons matrix
    # Occupancy: psi1 from the season-1 weights, at the shared state-arm
    # encoding. The colonization / extinction arms carry no field and keep the
    # plain M = 1000 inflation.
    M <- 1000L
    occ_block <- .tobs_encode_state_block(
      weights      = w[, 1],
      any_det      = ad[seq(1, by = n_seasons, length.out = n_sites)],
      X_occ        = X_occ,
      spatial_occ  = spatial_occ,
      latent_prior = latent_prior
    )

    # Colonization (z_{t-1}=0 -> z_t=1) and extinction (z_{t-1}=1 -> z_t=0) from
    # the EXACT smoothed pairwise joints the E-step accumulated. For a
    # constant-rate arm col_y / ext_y are the expected colonization / extinction
    # events summed over a site's T-1 intervals, and the same M pseudo-binomial
    # inflation encodes ONE logistic observation per site (as for the occupancy
    # arm). For a SEASON-VARYING arm (gcol33/tulpaObs#124) the per-interval joints
    # are kept and encode ONE logistic row per (site, interval), site-major
    # interval-minor to match the long-form design; a covariate on that interval
    # then drives the rate. `t(mat)` flattens site-major interval-minor.
    trans_block <- function(sv, y_agg, n_agg, y_mat, n_mat, Xarm) {
      if (sv) {
        # Per-interval M-step = a WEIGHTED logistic regression: for interval
        # (i, t) the response is the transition probability given the origin
        # state, r = E[event] / P(origin state) = y_mat / n_mat, and the weight is
        # P(origin state | y) = n_mat. This is the exact conditional-expectation
        # M-step for a per-interval covariate. It replaces the aggregated arm's
        # M-pseudo-binomial inflation, which -- applied per interval -- makes each
        # row individually near-separable (a confident interval has r ~ 0 or 1),
        # so the inner Newton hits numerical separation and returns ~the initial
        # zero slope. A fractional response with a fractional weight is exactly a
        # soft-label logistic and is well conditioned.
        nn <- as.vector(t(n_mat))
        r  <- as.vector(t(y_mat)) / pmax(nn, 1e-12)
        r  <- pmin(pmax(r, 0), 1)
        list(y = r, n_trials = rep(1, length(r)), weights = nn,
             X = Xarm, family = "binomial")
      } else {
        # Aggregated per-site: guard against an all-other-origin site (no rows).
        yv2 <- as.integer(round(y_agg * M))
        nv2 <- pmax(as.integer(round(n_agg * M)), 1L)
        yv2 <- pmin(pmax(yv2, 0L), nv2)
        list(y = yv2, n_trials = nv2, X = Xarm, family = "binomial")
      }
    }
    col_block <- trans_block(col_sv, attr(w, "col_y"), attr(w, "col_n"),
                             attr(w, "col_y_mat"), attr(w, "col_n_mat"), X_col)
    ext_block <- trans_block(ext_sv, attr(w, "ext_y"), attr(w, "ext_n"),
                             attr(w, "ext_y_mat"), attr(w, "ext_n_mat"), X_ext)

    # Detection: per-(site, season) rows weighted by w[i, t] = P(z_it = 1 | y).
    # Replaces the legacy hard threshold (w > 0.5) which silently dropped
    # site-seasons in the boundary regime and double-counted detection
    # evidence for site-seasons in the high-confidence regime. A constant-detection
    # arm's X_det is site-indexed, so per-season rows read the site's covariates
    # (design_row = i); a season-varying arm's X_det is the long-form
    # [(site x season) x p] design, so the row is the (site, season) index
    # (i - 1) * n_seasons + t (gcol33/tulpaObs#124).
    rows_i <- integer(n_sites * n_seasons)
    det_count <- integer(n_sites * n_seasons)
    vis_count <- integer(n_sites * n_seasons)
    w_it <- numeric(n_sites * n_seasons)
    n_rows <- 0L
    for (i in seq_len(n_sites)) {
      for (t in seq_len(n_seasons)) {
        idx <- (i - 1) * n_seasons + t
        if (nv[idx] <= 0) next
        base <- (i - 1) * n_seasons * max_visits + (t - 1) * max_visits
        dc <- 0L; vc <- 0L
        for (j in seq_len(nv[idx])) {
          v <- y_flat[base + j]
          if (v >= 0) { vc <- vc + 1L; if (v == 1) dc <- dc + 1L }
        }
        if (vc == 0L) next
        w_eff <- if (dc > 0L) 1 else w[i, t]
        if (w_eff <= 1e-6) next
        n_rows <- n_rows + 1L
        rows_i[n_rows] <- if (det_sv) (i - 1L) * n_seasons + t else i
        det_count[n_rows] <- dc
        vis_count[n_rows] <- vc
        w_it[n_rows] <- w_eff
      }
    }
    if (n_rows > 0L) {
      rows_i <- rows_i[seq_len(n_rows)]
      det_count <- det_count[seq_len(n_rows)]
      vis_count <- vis_count[seq_len(n_rows)]
      w_it <- w_it[seq_len(n_rows)]
      det_block <- list(y = det_count, n_trials = vis_count,
                        X = X_det[rows_i, , drop = FALSE],
                        weights = w_it, family = "binomial")
    } else {
      det_block <- list(y = integer(0), n_trials = integer(0),
                        X = X_det[integer(0), , drop = FALSE],
                        weights = numeric(0), family = "binomial")
    }

    list(
      occ = occ_block,
      det = det_block,
      col = col_block,
      ext = ext_block
    )
  }

  z_draw <- function(weights, ...) {
    w <- weights
    z <- matrix(0L, n_sites, n_seasons)
    for (i in seq_len(n_sites)) {
      for (t in seq_len(n_seasons)) {
        idx <- (i - 1) * n_seasons + t
        if (ad[idx]) z[i, t] <- 1L
        else z[i, t] <- rbinom(1, 1, clamp_w(w[i, t]))
      }
    }
    z
  }

  hard_encode <- function(z, ...) {
    z1 <- z[, 1]
    # Colonization/extinction from hard transitions. For a constant-rate arm the
    # per-site transition counts encode one logistic row per site; a
    # season-varying arm (gcol33/tulpaObs#124) keeps each interval separate (one
    # row per (site, interval), n_trials = 1 at the matching origin state and 0
    # elsewhere so the engine drops the non-matching rows), site-major
    # interval-minor to match the long-form design.
    col_y <- integer(n_sites); col_n <- integer(n_sites)
    ext_y <- integer(n_sites); ext_n <- integer(n_sites)
    cym <- if (col_sv) matrix(0L, n_sites, n_int) else NULL
    cnm <- if (col_sv) matrix(0L, n_sites, n_int) else NULL
    eym <- if (ext_sv) matrix(0L, n_sites, n_int) else NULL
    enm <- if (ext_sv) matrix(0L, n_sites, n_int) else NULL
    for (i in seq_len(n_sites)) {
      for (t in 2:n_seasons) {
        iv <- t - 1L
        if (z[i, t - 1] == 0) {
          col_n[i] <- col_n[i] + 1L; if (z[i, t] == 1) col_y[i] <- col_y[i] + 1L
          if (col_sv) { cnm[i, iv] <- 1L; if (z[i, t] == 1) cym[i, iv] <- 1L }
        }
        if (z[i, t - 1] == 1) {
          ext_n[i] <- ext_n[i] + 1L; if (z[i, t] == 0) ext_y[i] <- ext_y[i] + 1L
          if (ext_sv) { enm[i, iv] <- 1L; if (z[i, t] == 0) eym[i, iv] <- 1L }
        }
      }
    }
    col_n <- pmax(col_n, 1L); ext_n <- pmax(ext_n, 1L)
    col_hard <- if (col_sv)
      list(y = as.integer(as.vector(t(cym))),
           n_trials = as.integer(as.vector(t(cnm))), X = X_col,
           family = "binomial")
      else list(y = col_y, n_trials = col_n, X = X_col, family = "binomial")
    ext_hard <- if (ext_sv)
      list(y = as.integer(as.vector(t(eym))),
           n_trials = as.integer(as.vector(t(enm))), X = X_ext,
           family = "binomial")
      else list(y = ext_y, n_trials = ext_n, X = X_ext, family = "binomial")

    # Detection counts among occupied (z = 1) site-seasons. A constant-detection
    # arm aggregates a site's detections over its seasons onto one per-site row;
    # a season-varying arm (gcol33/tulpaObs#124) keeps one row per (site, season)
    # so each reads the season's own detection covariate (the long-form X_det row
    # (i - 1) * n_seasons + t).
    if (det_sv) {
      row_map <- integer(n_sites * n_seasons)
      dcv <- integer(n_sites * n_seasons); vcv <- integer(n_sites * n_seasons)
      nr <- 0L
      for (i in seq_len(n_sites)) {
        for (t in seq_len(n_seasons)) {
          if (z[i, t] != 1 || nv[(i-1)*n_seasons+t] <= 0) next
          base <- (i-1)*n_seasons*max_visits + (t-1)*max_visits
          dc <- 0L; vc <- 0L
          for (j in seq_len(nv[(i-1)*n_seasons+t])) {
            v <- y_flat[base + j]
            if (v >= 0) { vc <- vc + 1L; if (v == 1) dc <- dc + 1L }
          }
          if (vc == 0L) next
          nr <- nr + 1L
          row_map[nr] <- (i - 1L) * n_seasons + t; dcv[nr] <- dc; vcv[nr] <- vc
        }
      }
      det_hard <- if (nr > 0L)
        list(y = dcv[seq_len(nr)], n_trials = vcv[seq_len(nr)],
             X = X_det[row_map[seq_len(nr)], , drop = FALSE], family = "binomial")
      else NULL
    } else {
      total_det <- integer(n_sites); total_vis <- integer(n_sites)
      for (i in seq_len(n_sites)) {
        for (t in seq_len(n_seasons)) {
          if (z[i, t] == 1 && nv[(i-1)*n_seasons+t] > 0) {
            base <- (i-1)*n_seasons*max_visits + (t-1)*max_visits
            for (j in seq_len(nv[(i-1)*n_seasons+t])) {
              v <- y_flat[base + j]
              if (v >= 0) { total_vis[i] <- total_vis[i]+1L; if (v==1) total_det[i] <- total_det[i]+1L }
            }
          }
        }
      }
      dk <- total_vis > 0
      det_hard <- if (sum(dk) > 0)
        list(y = total_det[dk], n_trials = total_vis[dk],
             X = X_det[dk,,drop=FALSE], family = "binomial")
      else NULL
    }

    occ_block <- list(y = z1, n_trials = rep(1L, n_sites), X = X_occ,
                      family = "binomial")
    occ_block <- .attach_spatial_spde(occ_block, spatial_occ)
    list(
      occ = occ_block,
      det = det_hard,
      col = col_hard,
      ext = ext_hard
    )
  }

  init <- list(
    occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)),
    det = list(beta = rep(0, p_det), se = rep(1, p_det)),
    col = list(beta = rep(0, p_col), se = rep(1, p_col)),
    ext = list(beta = rep(0, p_ext), se = rep(1, p_ext))
  )

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init,
       p_per_submodel = c(occ = p_occ, det = p_det, col = p_col, ext = p_ext))
}

# ============================================================================
# Integrated occupancy callbacks
# ============================================================================
build_integrated_callbacks <- function(model, spatial = NULL,
                                       latent_prior = NULL) {
  y_sources <- model$y_sources
  site_maps <- model$site_maps
  X_occ <- model$X_processes[[1]]
  n_sites <- model$n_sites
  n_sources <- model$n_sources
  p_occ <- ncol(X_occ)

  # The shared psi field enters the state arm (one state row per site, the
  # identity map -- the proven single-season path). A field on the detection
  # arm enters every source's per-source detection block, broadcast onto that
  # source's sites (`src_rows`) and folded into the response by count-scaling,
  # exactly as the single-season detection arm does. The field is fit
  # independently per source block (one realization per submodel block; a
  # genuinely shared realization across sources needs the copy() path).
  spatial_occ <- .spatial_for_arm(spatial, 1L)
  spatial_det <- .spatial_for_arm(spatial, 2L)

  # Per-source detection info
  src_info <- lapply(seq_len(n_sources), function(s) {
    ys <- y_sources[[s]]; ns <- nrow(ys); mv <- ncol(ys)
    nv <- integer(ns); nd <- integer(ns); ad <- logical(ns)
    for (i in seq_len(ns)) {
      valid <- ys[i, ] >= 0; nv[i] <- sum(valid)
      nd[i] <- sum(ys[i, valid] == 1); ad[i] <- nd[i] > 0
    }
    X_det <- model$X_processes[[1 + s]]
    src_rows <- site_maps[[s]] + 1L
    # Detection-arm field projection broadcast onto this source's sites.
    spatial_det_s <- if (!is.null(spatial_det))
      .tobs_spde_broadcast_spec(spatial_det, src_rows) else NULL
    list(nv = nv, nd = nd, ad = ad, X_det = X_det[src_rows, , drop = FALSE],
         p_det = ncol(X_det), keep = nv > 0, src_rows = src_rows,
         spatial_det = spatial_det_s)
  })

  # Global detection status per site
  any_det_global <- logical(n_sites)
  for (s in seq_len(n_sources)) {
    for (j in seq_along(src_info[[s]]$src_rows)) {
      if (src_info[[s]]$ad[j]) any_det_global[src_info[[s]]$src_rows[j]] <- TRUE
    }
  }

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    eta_occ <- as.vector(X_occ %*% beta_occ)
    sp_off <- .spatial_eta_offset(spatial_occ, fits$occ, p_occ)
    if (length(sp_off) == n_sites) eta_occ <- eta_occ + sp_off
    # Nested-Laplace: the shared psi field informs which undetected sites are
    # occupied, so the E-step weight must see it (as the single-season arm does).
    eta_occ <- .nested_state_eta(eta_occ, latent_prior, fits$occ, p_occ, n_sites)
    psi <- plogis(eta_occ)
    weights <- psi  # Prior occupancy
    for (s in seq_len(n_sources)) {
      si <- src_info[[s]]
      beta_det <- extract_beta(fits[[paste0("det", s)]], si$p_det)
      eta_det <- as.vector(si$X_det %*% beta_det)
      det_off <- .spatial_eta_offset(si$spatial_det, fits[[paste0("det", s)]],
                                     si$p_det)
      if (length(det_off) == length(si$src_rows)) eta_det <- eta_det + det_off
      p_s <- plogis(eta_det)
      for (j in seq_along(si$src_rows)) {
        i <- si$src_rows[j]
        if (si$ad[j]) { weights[i] <- 1 }
        else if (si$nv[j] > 0) {
          prod_1mp <- (1 - p_s[j])^si$nv[j]
          num <- weights[i] * prod_1mp
          weights[i] <- num / (num + (1 - weights[i]) + 1e-10)
        }
      }
    }
    list(weights = weights)
  }

  m_step_encode <- function(weights, ...) {
    occ_block <- .tobs_encode_state_block(weights, any_det_global, X_occ,
                                          spatial_occ, latent_prior)
    specs <- list(occ = occ_block)
    # Per-source detection blocks: weight each row by w_i at the global
    # site mapped through src_rows. Sites where the E-step says "almost
    # certainly empty" drop out of every source's detection fit.
    for (s in seq_len(n_sources)) {
      si <- src_info[[s]]
      w_src <- weights[si$src_rows]
      w_src[si$ad] <- 1
      if (!is.null(si$spatial_det)) {
        # SPDE detection field on this source: fold the occupancy weight into
        # the binomial response by count-scaling (y = round(w nd), n =
        # round(w nv)) so every source row stays aligned with the broadcast
        # mesh projection A; a near-empty site collapses to a (0, 0) row.
        y_det_w <- as.integer(round(w_src * si$nd))
        n_det_w <- as.integer(round(w_src * si$nv))
        y_det_w <- pmin(pmax(y_det_w, 0L), n_det_w)
        det_block <- list(y = y_det_w, n_trials = n_det_w, X = si$X_det,
                          family = "binomial")
        specs[[paste0("det", s)]] <- .attach_spatial_spde(det_block,
                                                          si$spatial_det)
      } else {
        dk <- si$keep & (w_src > 1e-6)
        specs[[paste0("det", s)]] <- list(y = si$nd[dk], n_trials = si$nv[dk],
                                          X = si$X_det[dk,,drop=FALSE],
                                          weights = w_src[dk],
                                          family = "binomial")
      }
    }
    specs
  }

  z_draw <- function(weights, ...) {
    z <- as.integer(any_det_global)
    undet <- !any_det_global
    z[undet] <- rbinom(sum(undet), 1, clamp_w(weights[undet]))
    z
  }

  hard_encode <- function(z, ...) {
    occ_block <- list(y = z, n_trials = rep(1L, n_sites), X = X_occ,
                      family = "binomial")
    occ_block <- .attach_spatial_spde(occ_block, spatial_occ)
    specs <- list(occ = occ_block)
    for (s in seq_len(n_sources)) {
      si <- src_info[[s]]
      if (!is.null(si$spatial_det)) {
        keep_det_mask <- (z[si$src_rows] == 1L) & (si$nv > 0L)
        n_det_h <- ifelse(keep_det_mask, si$nv, 0L)
        y_det_h <- ifelse(keep_det_mask, si$nd, 0L)
        det_block <- list(y = as.integer(y_det_h),
                          n_trials = as.integer(n_det_h),
                          X = si$X_det, family = "binomial")
        specs[[paste0("det", s)]] <- .attach_spatial_spde(det_block,
                                                          si$spatial_det)
      } else {
        occ_local <- z[si$src_rows] == 1L & si$nv > 0
        if (any(occ_local)) {
          specs[[paste0("det", s)]] <- list(y = si$nd[occ_local],
                                            n_trials = si$nv[occ_local],
                                            X = si$X_det[occ_local,,drop=FALSE],
                                            family = "binomial")
        }
      }
    }
    specs
  }

  init <- list(occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)))
  p_sub <- c(occ = p_occ)
  for (s in seq_len(n_sources)) {
    nm <- paste0("det", s)
    init[[nm]] <- list(beta = rep(0, src_info[[s]]$p_det), se = rep(1, src_info[[s]]$p_det))
    p_sub[nm] <- src_info[[s]]$p_det
  }

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init, p_per_submodel = p_sub)
}


# ============================================================================
# Count / relative-abundance GLMM callbacks (no detection, no latent state)
# ============================================================================
# The observed count is the response of a single GLMM block. There is no latent
# variable, so the "E-step" just returns the observed y and the M-step fits the
# block; the EM loop converges immediately (each M-step block fit IS the MLE for
# the given dispersion). The negbin size / Gaussian residual variance is a fixed
# `phi` supplied on the model (`model$count_phi`); .dispatch_count updates it in
# an outer loop and refits.
build_count_callbacks <- function(model, spatial = NULL) {
  X        <- model$X_processes[[1]]
  N        <- model$N
  p        <- ncol(X)
  response <- model$response %||% "poisson"
  fam <- switch(response,
    poisson  = "poisson",
    negbin   = "neg_binomial_2",
    gaussian = "gaussian",
    binomial = "binomial",
    stop(sprintf("count(): unsupported response '%s'.", response),
         call. = FALSE))
  phi    <- model$count_phi %||% 1.0
  is_int <- response %in% c("poisson", "negbin", "binomial")
  yv     <- if (is_int) as.integer(model$y_count) else as.numeric(model$y_count)
  # Binomial carries a per-site trial count; every other response is a plain
  # single-value GLMM block.
  n_trials <- if (identical(response, "binomial"))
                as.integer(model$n_trials %||% rep(1L, N)) else NULL

  mk_block <- function() {
    blk <- list(y = yv, X = X, family = fam)
    if (identical(fam, "binomial")) blk$n_trials <- n_trials
    else if (fam != "poisson") blk$phi <- phi
    blk
  }

  e_step        <- function(fits, ...) list(weights = as.numeric(yv))
  m_step_encode <- function(weights, ...) list(occ = mk_block())
  z_draw        <- function(weights, ...) yv
  hard_encode   <- function(z, ...) list(occ = mk_block())
  init          <- list(occ = list(beta = rep(0, p), se = rep(1, p)))

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init, p_per_submodel = c(occ = p))
}

