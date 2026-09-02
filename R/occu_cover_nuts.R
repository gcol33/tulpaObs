# occu_cover_nuts.R - NUTS target for the non-spatial joint occupancy +
# cover-hurdle family (occu_cover()).
#
# The Laplace fit (.tobs_fit_occu_cover) sums the latent occupancy state z out in
# closed form (two states per cell) and returns a Gaussian observed-Fisher
# posterior over the packed coefficient vector
#   theta = c(beta_psi, beta_p, beta_pos, log_dispersion).
# NUTS instead samples the exact marginal posterior of that vector, giving
# calibrated (non-Gaussian) intervals and the per-draw pointwise likelihood that
# WAIC / LOO need. The bare target carries no latent structure, so its parameter
# vector is just the flat three-arm coefficient block plus one log-dispersion
# scalar; two optional blocks trail it, each a no-op when absent -- a coupled
# areal field on the latent state, and one observation-arm random intercept per
# grouping factor on the detection / positive-cover arms.
#
# .tobs_occu_cover_nuts_logpost is the R oracle: it recomputes the joint
# log-posterior and gradient exactly as the C++ FullGradFn (src/occu_cover_nuts.cpp,
# cpp_occu_cover_nuts) does, and a byte-exact test cross-checks the two before the
# sampler is trusted. The per-cell marginal mirrors .occu_cover_site_ll; the cover
# arm and the no-detection mixture reuse the same closed forms the Laplace path uses.


# ---------------------------------------------------------------------------
# Coupled areal field with sampled hyperparameters
#
# R mirror of src/nuts_field_hyper.h. The field is
#   z = sigma * ( B1 %*% (s1(rho) * raw1) + s2(rho) * raw2 )
# over a FIXED basis B1, with the field SD sigma, the mixing / spatial-
# correlation rho and the cross-arm copy amplitude alpha either sampled as
# bounded coordinates or pinned. A sampled hyper rides an unconstrained u,
#   t = t_lo + (t_hi - t_lo) * plogis(u),   value = inv_link(t),
# with `t` the coordinate the nested-Laplace outer grid spaces its nodes in
# (log for sigma / alpha, logit for rho) and [t_lo, t_hi] that grid's own span.
#
# Each coordinate carries a DECLARED prior, and the sampler's log-density is
# that prior carried to u,
#   log p(u) = log p_t(t) + log(dt/du),   dt/du = (t_hi - t_lo) e (1 - e),
# with e = plogis(u). Flat in t over [t_lo, t_hi] is the measure a grid axis
# with no declared density integrates against, and leaves the normalised
# log(e) + log(1 - e). An Exponential(rate) declared on the NATURAL scale of a
# log-link coordinate carries the log-link Jacobian too, p_t(t) = rate
# exp(-rate v) v, so its log-density is
#   log(e) + log(1 - e) + log(rate) - rate v + log(v) + log(t_hi - t_lo)
# and its natural-scale score 1/v - rate joins the data score in the chain rule.
# That is the copy scale's penalized-complexity slab (tulpa's
# .hyper_copy_slab_density): the outer grid weighs the same Exponential against
# a point mass at v = 0, which a gradient sampler cannot visit, so a sampled
# copy scale carries the continuum alone and is the posterior conditional on a
# coupled field.
#
# A block may carry a per-site design WEIGHT, which is what makes it a
# spatially-varying coefficient rather than a second intercept field: site i
# loads w_i * z[unit(i)] on psi and alpha * w_i * z[unit(i)] on the cover arm.
# Several blocks stack, each with its own basis, map, weight and (sigma, rho,
# alpha) coordinates, laid out back to back.
# ---------------------------------------------------------------------------

.OCHF_CONST <- 0L; .OCHF_BYM2_STR <- 1L; .OCHF_BYM2_IID <- 2L; .OCHF_CAR <- 3L

# Declared prior codes, mirroring HyperPrior in src/nuts_field_hyper.h.
.OCHF_PRIOR_FLAT <- 0L; .OCHF_PRIOR_EXP <- 1L

# Rate of the copy scale's penalized-complexity slab for a declared alpha grid
# whose largest node is `upper`. Read off tulpa's own declared density rather
# than restated, so the sampler's slab cannot drift from the measure the outer
# grid integrates: that density is log p(x) = log(rate) - rate x, so
# p(0) - p(1) on the log scale IS the rate. NULL when the grid declares no
# positive node to anchor it.
.ochf_copy_slab_rate <- function(upper) {
  if (length(upper) != 1L || !is.finite(upper) || upper <= 0) return(NULL)
  d <- tulpa:::.hyper_copy_slab_density(upper)
  if (is.null(d)) return(NULL)
  as.numeric(d(0) - d(1))
}

# `control$copy.slab` for the sampled copy amplitude. Unset, the sampler takes
# the measure the outer grid declares for that axis, so both engines target the
# same model unless the caller says otherwise; "flat" asks for the alternative
# the grid also accepts, flat in log alpha over the span its nodes tile.
.occu_cover_nuts_copy_slab <- function(x) {
  if (is.null(x)) return(tulpa:::.hyper_check_copy_slab(NULL))
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !x %in% c("flat", "exponential")) {
    stop("control$copy.slab must be \"flat\" or \"exponential\".", call. = FALSE)
  }
  x
}

# Does this field block carry a copy onto the cover arm? A resolved amplitude of
# 0 (or a non-finite one) IS the decoupled model -- there is no coupling, so
# there is no amplitude to report. The sampled-path spelling of
# `.tobs_alpha_axis_decoupled()`.
.ochf_has_copy <- function(alpha) isTRUE(is.finite(alpha) && alpha > 0)

# A block's reported copy amplitude, or 0 when it carries no copy and therefore
# no `alpha` entry. `[[` errors on a missing name, so every reader of a
# possibly-absent amplitude goes through this.
.ochf_block_alpha <- function(hyper, suffix) {
  nm <- paste0("alpha", suffix)
  if (nm %in% names(hyper)) hyper[[nm]] else 0
}

.ochf_inv_link <- function(t, link) if (link == 1L) stats::plogis(t) else exp(t)
.ochf_link     <- function(v, link) if (link == 1L) stats::qlogis(v) else log(v)

# One hyper's value at `theta` plus the chain-rule pieces (see HyperValue).
.ochf_value <- function(h, theta) {
  if (is.null(h$coord))
    return(list(value = h$fixed, dvalue_dt = 0, dt_du = 0, e = 0))
  e <- stats::plogis(theta[h$coord])
  t <- h$t_lo + (h$t_hi - h$t_lo) * e
  v <- .ochf_inv_link(t, h$link)
  list(value = v, dvalue_dt = if (h$link == 1L) v * (1 - v) else v,
       dt_du = (h$t_hi - h$t_lo) * e * (1 - e), e = e)
}

# Canonical view of a field description. Accepts the sampled-hyper spec entries
# the C++ block reads AND the legacy fixed-hyper form (`Linv` / `field_load` +
# `alpha`), which resolves to sigma pinned at 1 over a constant scaling -- the
# fixed loading already carries sigma and rho in its columns, so that
# configuration reproduces it exactly. `total` is the 1-based index of the
# log-dispersion coordinate; the whitened field follows it, then each sampled
# hyper in the order (sigma, rho, alpha).
.ochf_view <- function(field, total) {
  B1      <- field$field_load %||% field$Linv
  n_units <- as.integer(field$n_field_units %||% nrow(B1))
  m1      <- ncol(B1)
  has_iid <- isTRUE(as.logical(field$field_has_iid %||% FALSE))
  n_raw   <- m1 + if (has_iid) n_units else 0L
  k       <- total + n_raw
  slot <- function(nm, link, dflt) {
    key <- function(suffix) field[[paste0("field_", nm, "_", suffix)]]
    h <- list(link = link, fixed = key("fixed") %||% dflt, coord = NULL,
              prior = .OCHF_PRIOR_FLAT, rate = 0)
    lo <- key("lo"); hi <- key("hi")
    if (!is.null(lo) && !is.null(hi)) {
      k <<- k + 1L
      h$t_lo <- lo; h$t_hi <- hi; h$coord <- k
      pr <- key("prior")
      if (!is.null(pr) && as.integer(pr) == .OCHF_PRIOR_EXP) {
        if (link != 0L)
          stop(sprintf(paste0(
            "field_%s_prior: an Exponential prior is declared on a hyper's ",
            "natural scale, which needs a log-link coordinate."), nm),
            call. = FALSE)
        rate <- key("rate")
        if (is.null(rate) || !is.finite(rate) || rate <= 0)
          stop(sprintf("field_%s_rate must be a positive rate.", nm),
               call. = FALSE)
        h$prior <- .OCHF_PRIOR_EXP; h$rate <- as.numeric(rate)
      }
    }
    h
  }
  sigma <- slot("sigma", 0L, 1)
  rho   <- slot("rho",   1L, 1)
  alpha <- slot("alpha", 0L, 0)
  if (is.null(alpha$coord) && is.null(field$field_alpha_fixed))
    alpha$fixed <- field$field_alpha %||% field$alpha %||% 0
  list(n_units = n_units, m1 = m1, B1 = B1,
       scale1 = as.integer(field$field_scale1 %||% .OCHF_CONST),
       has_iid = has_iid, sf = as.numeric(field$field_sf %||% 1),
       lambda = as.numeric(field$field_lambda %||% numeric(0)),
       n_raw = n_raw, o_raw = total, n_coord = k - total,
       field_map = as.integer(field$field_map),
       weight = if (is.null(field$field_weight)) NULL
                else as.numeric(field$field_weight),
       sigma = sigma, rho = rho, alpha = alpha)
}

# Every field block a spec (or a fit's field description) declares, each viewed
# at its own offset. Accepts the plural spelling (`field_blocks`, or a bare list
# of descriptions) and the single-block one, which is the same list of length
# one -- so the oracle reads exactly what the C++ builder reads.
.ochf_views <- function(field, total) {
  if (is.null(field)) return(list())
  blocks <- if (!is.null(field$field_blocks)) field$field_blocks
            else if (is.null(field$field_map)) field    # bare list of blocks
            else list(field)
  out <- vector("list", length(blocks))
  k <- total
  for (b in seq_along(blocks)) {
    out[[b]] <- .ochf_view(blocks[[b]], k)
    k <- k + out[[b]]$n_coord
  }
  out
}

# Per-site loading of one block on the state arm: w_i * z[unit(i)].
.ochf_site_value <- function(fv, fs) {
  z <- fs$z[fv$field_map]
  if (is.null(fv$weight)) z else fv$weight * z
}

# Forward pass: hyper values, per-column block scalings, and the field z.
.ochf_forward <- function(fv, theta) {
  sg <- .ochf_value(fv$sigma, theta)
  rh <- .ochf_value(fv$rho,   theta)
  al <- .ochf_value(fv$alpha, theta)
  rho <- rh$value
  s1 <- rep(1, fv$m1); ds1 <- rep(0, fv$m1)
  if (fv$scale1 == .OCHF_BYM2_STR) {
    s1[]  <- sqrt(rho / fv$sf)
    ds1[] <- 1 / (2 * sqrt(rho * fv$sf))
  } else if (fv$scale1 == .OCHF_CAR) {
    a   <- 1 - rho * fv$lambda
    s1  <- 1 / sqrt(a)
    ds1 <- 0.5 * fv$lambda * s1 / a
  }
  s2 <- if (fv$has_iid) sqrt(1 - rho) else 0
  ds2 <- if (fv$has_iid) -1 / (2 * sqrt(1 - rho)) else 0
  raw1 <- theta[fv$o_raw + seq_len(fv$m1)]
  raw2 <- if (fv$has_iid) theta[fv$o_raw + fv$m1 + seq_len(fv$n_units)]
          else numeric(0)
  z <- as.numeric(fv$B1 %*% (s1 * raw1))
  if (fv$has_iid) z <- z + s2 * raw2
  z <- sg$value * z
  list(sigma = sg, rho = rh, alpha = al, s1 = s1, ds1 = ds1, s2 = s2, ds2 = ds2,
       raw1 = raw1, raw2 = raw2, z = z)
}

# Backward pass. `g_z` is d log L / d z per unit (the caller folds in alpha times
# the copied arm's score); `g_alpha_data` is d log L / d alpha from the data.
# Returns the field coordinates' gradient (raw block then sampled hypers) and
# the whitened-field prior plus the sampled hypers' log-densities.
.ochf_backward <- function(fv, fs, g_z, g_alpha_data) {
  sigma <- fs$sigma$value
  BtG   <- as.numeric(crossprod(fv$B1, g_z))
  lp    <- -0.5 * sum(fs$raw1^2)
  g_raw <- sigma * fs$s1 * BtG - fs$raw1
  g_rho <- sum(sigma * fs$ds1 * fs$raw1 * BtG)
  if (fv$has_iid) {
    lp    <- lp - 0.5 * sum(fs$raw2^2)
    g_raw <- c(g_raw, sigma * fs$s2 * g_z - fs$raw2)
    g_rho <- g_rho + sum(sigma * fs$ds2 * fs$raw2 * g_z)
  }
  g_hyper <- numeric(0)
  # A declared natural-scale prior enters through the SAME chain rule as the
  # data score -- its natural-scale derivative is added to `dlp_dvalue` -- so
  # the transform's own (1 - 2e) term is written once.
  add <- function(h, hv, dlp_dvalue) {
    if (is.null(h$coord)) return(invisible(NULL))
    lp <<- lp + log(hv$e) + log(1 - hv$e)
    if (h$prior == .OCHF_PRIOR_EXP) {
      v <- hv$value
      dlp_dvalue <- dlp_dvalue + 1 / v - h$rate
      lp <<- lp + log(h$rate) - h$rate * v + log(v) + log(h$t_hi - h$t_lo)
    }
    g_hyper <<- c(g_hyper, dlp_dvalue * hv$dvalue_dt * hv$dt_du + (1 - 2 * hv$e))
  }
  add(fv$sigma, fs$sigma, if (sigma > 0) sum(g_z * fs$z) / sigma else 0)
  add(fv$rho,   fs$rho,   g_rho)
  add(fv$alpha, fs$alpha, g_alpha_data)
  list(grad = c(g_raw, g_hyper), lp = lp)
}


# ---------------------------------------------------------------------------
# Observation-arm random intercepts
#
# R mirror of src/nuts_re_block.h. Each grouping factor on the detection or
# positive-cover formula is ONE non-centered block: n_groups whitened
# coordinates z_g ~ N(0, 1) plus the SAMPLED log_sigma_re under a
# N(0, sigma_lsd^2) prior, adding sigma_re * z[code(row)] to that arm's per-row
# linear predictor. The prior is written on the sampled coordinate itself
# (log sigma_re), so the target needs no change-of-variables term: at z = 0 the
# data term does not depend on log_sigma_re and the target reduces to exactly
# that Gaussian. A code of 0 (padded visit / unseen level) carries no effect,
# matching the deterministic engine's 0 scatter sentinel.
#
# `blocks` is a list of list(arm = 1 detection / 2 positive, group = per-row
# codes, n_groups, sigma_lsd); `base` is the last coordinate before the first
# block (the log-dispersion index plus the field block's own coordinates).
# ---------------------------------------------------------------------------

.ocre_view <- function(blocks, base) {
  if (is.null(blocks) || length(blocks) == 0L) return(list())
  k <- base
  lapply(blocks, function(b) {
    G <- as.integer(b$n_groups)
    # `o_z` is a 0-based OFFSET (the whitened block is theta[o_z + 1:G]);
    # `o_logsig` is the 1-based INDEX of the log-SD coordinate that follows it.
    out <- list(arm = as.integer(b$arm), codes = as.integer(b$group), n_groups = G,
                sigma_lsd = as.numeric(b$sigma_lsd %||% 1.5),
                o_z = k, o_logsig = k + G + 1L)
    k <<- k + G + 1L
    out
  })
}

# Per-row offset sigma_re * z[code], 0 where the code is 0.
.ocre_offset <- function(rb, theta) {
  z   <- theta[rb$o_z + seq_len(rb$n_groups)]
  sig <- exp(theta[rb$o_logsig])
  ifelse(rb$codes > 0L, sig * z[pmax(rb$codes, 1L)], 0)
}

# Gradient + log-prior of one block. `g_row` is d log L / d eta_arm on the same
# per-row layout as `codes`.
.ocre_backward <- function(rb, theta, g_row) {
  z   <- theta[rb$o_z + seq_len(rb$n_groups)]
  ls  <- theta[rb$o_logsig]
  sig <- exp(ls)
  gz  <- vapply(seq_len(rb$n_groups),
                function(g) sum(g_row[rb$codes == g]), numeric(1))
  ils2 <- 1 / rb$sigma_lsd^2
  list(grad = c(sig * gz - z,
                sum(g_row * .ocre_offset(rb, theta)) - ils2 * ls),
       lp = -0.5 * sum(z^2) - 0.5 * ils2 * ls^2)
}


# Joint log-posterior + gradient of the non-spatial occu_cover coefficient vector
# theta = c(beta_psi, beta_p, beta_pos, log_disp) (beta_p / beta_pos each packed as
# the site-level block then the optional visit-level block, exactly as the fitter
# stacks them). Weak Gaussian priors N(0, sigma.beta^2) on every coefficient and a
# broad N(0, sigma.logdisp^2) on log_disp keep the ridge / dispersion proper
# without materially shifting the data-dominated optimum. Returns list(lp, grad)
# over the packed coordinates. This is the oracle the C++ FullGradFn mirrors.
.tobs_occu_cover_nuts_logpost <- function(theta, model, sigma.beta = 5,
                                          sigma.logdisp = 5, field = NULL,
                                          re = NULL) {
  pin         <- model$process_info
  p_occ       <- pin[[1L]]$p
  p_det_site  <- ncol(model$X_det_site)
  p_det_visit <- if (!is.null(model$X_det_visit)) ncol(model$X_det_visit) else 0L
  p_pos_site  <- ncol(model$X_pos_site)
  p_pos_visit <- if (!is.null(model$X_pos_visit)) ncol(model$X_pos_visit) else 0L
  p_p   <- p_det_site + p_det_visit
  p_pos <- p_pos_site + p_pos_visit
  n_coef <- p_occ + p_p + p_pos
  total  <- n_coef + 1L

  # Optional coupled field(s): each block's z over its own trailing slice of
  # theta, entering psi weighted by its per-site design column and the cover arm
  # scaled by its own copy amplitude alpha. sigma / rho / alpha are sampled
  # coordinates or pinned, as each block declares. The R oracle mirrors the C++
  # FullGradFn byte-for-byte; absent, the non-spatial target.
  fvs       <- .ochf_views(field, total)
  has_field <- length(fvs) > 0L
  # Per-block state, the per-site loading on psi, and the coordinates the field
  # blocks occupy in total (so any trailing block -- the RE blocks -- starts
  # where the C++ layout starts it).
  fss     <- lapply(fvs, .ochf_forward, theta = theta)
  f_site  <- lapply(seq_along(fvs), function(b) .ochf_site_value(fvs[[b]], fss[[b]]))
  n_field <- sum(vapply(fvs, `[[`, numeric(1), "n_coord"))
  psi_off <- if (has_field) Reduce(`+`, f_site) else numeric(model$n_sites)
  pos_off <- if (has_field)
    Reduce(`+`, lapply(seq_along(fvs),
                       function(b) fss[[b]]$alpha$value * f_site[[b]]))
    else numeric(model$n_sites)

  bo         <- theta[seq_len(p_occ)]
  bp_site    <- theta[p_occ + seq_len(p_det_site)]
  bp_visit   <- if (p_det_visit > 0L) theta[p_occ + p_det_site + seq_len(p_det_visit)] else numeric(0)
  bpos_site  <- theta[p_occ + p_p + seq_len(p_pos_site)]
  bpos_visit <- if (p_pos_visit > 0L) theta[p_occ + p_p + p_pos_site + seq_len(p_pos_visit)] else numeric(0)
  log_disp   <- theta[total]
  disp       <- exp(log_disp)
  pos_code   <- .occu_cover_pos_code(model$positive)

  N <- model$n_sites; J <- model$max_visits
  sgm <- function(e) 1 / (1 + exp(-e))

  eta_psi <- as.numeric(model$X_occ %*% bo) + psi_off
  psi     <- sgm(eta_psi)

  eta_p <- matrix(as.numeric(model$X_det_site %*% bp_site), N, J)   # site block broadcast
  if (p_det_visit > 0L) {
    eta_p <- eta_p + matrix(as.numeric(model$X_det_visit %*% bp_visit), N, J, byrow = TRUE)
  }
  eta_pos <- matrix(as.numeric(model$X_pos_site %*% bpos_site) + pos_off, N, J)
  if (p_pos_visit > 0L) {
    eta_pos <- eta_pos + matrix(as.numeric(model$X_pos_visit %*% bpos_visit), N, J, byrow = TRUE)
  }

  # Observation-arm random intercepts. Each block owns n_groups whitened z
  # coordinates and its own SAMPLED log_sigma_re, laid out back to back AFTER
  # the field block, and adds the per-row offset sigma_re * z[code] to its arm.
  # `codes` is site-major over the padded grid (row (i - 1) * J + v), so
  # `matrix(..., byrow = TRUE)` puts each code on its own (site, visit); a 0
  # code carries no effect.
  re_view <- .ocre_view(re, total + n_field)
  for (rb in re_view) {
    off <- matrix(.ocre_offset(rb, theta), N, J, byrow = TRUE)
    if (rb$arm == 1L) eta_p   <- eta_p   + off
    else              eta_pos <- eta_pos + off
  }

  valid <- model$valid; y <- model$y; y_pos <- model$y_pos

  g_eta_psi <- numeric(N)
  g_eta_p   <- matrix(0, N, J)
  g_eta_pos <- matrix(0, N, J)
  g_logdisp <- 0
  lp_data   <- 0

  for (i in seq_len(N)) {
    vv <- which(valid[i, ])
    if (length(vv) == 0L) next
    any_det <- any(y[i, vv] == 1L)
    if (any_det) {
      lp_data <- lp_data + log(psi[i])
      g_eta_psi[i] <- 1 - psi[i]
      for (v in vv) {
        pv <- sgm(eta_p[i, v])
        if (y[i, v] == 1L) { lp_data <- lp_data + log(pv);     g_eta_p[i, v] <- 1 - pv }
        else               { lp_data <- lp_data + log(1 - pv); g_eta_p[i, v] <- -pv }
      }
      # Cover factor at detected visits with an observed cover; a missing (NA)
      # cover drops out (missing-at-random cover), the detection loop above still
      # counts the visit.
      for (v in vv[y[i, vv] == 1L & is.finite(y_pos[i, vv])]) {
        ev <- eta_pos[i, v]; yy <- y_pos[i, v]
        if (pos_code == 3L) {            # beta
          mu <- sgm(ev); a <- mu * disp; b <- (1 - mu) * disp
          # `.tobs_log_safe`, matching src/tobs_math.h::log_safe, which the
          # C++ target this oracle is asserted byte-for-byte against uses. A
          # raw log() made the two differ by -Inf against a finite value at a
          # cover of exactly 0 or 1, so an FD gradient check on boundary data
          # compared a different function than the sampler ran.
          ly <- .tobs_log_safe(yy); l1my <- .tobs_log_safe(1 - yy)
          lp_data <- lp_data + lgamma(disp) - lgamma(a) - lgamma(b) +
                     (a - 1) * ly + (b - 1) * l1my
          g_eta_pos[i, v] <- disp * mu * (1 - mu) *
                             (-digamma(a) + digamma(b) + ly - l1my)
          g_logdisp <- g_logdisp + disp * (digamma(disp) - mu * digamma(a) -
                       (1 - mu) * digamma(b) + mu * ly + (1 - mu) * l1my)
        } else if (pos_code == 4L) {     # identity-Gaussian (#112): raw response
          sig <- disp; r <- (yy - ev) / sig
          lp_data <- lp_data - .tobs_log_safe(sig) -
                     0.5 * log(2 * pi) - 0.5 * r * r
          g_eta_pos[i, v] <- r / sig
          g_logdisp <- g_logdisp + (r * r - 1)
        } else {                         # lognormal: Gaussian on log(cover)
          sig <- disp; lyy <- .tobs_log_safe(yy); r <- (lyy - ev) / sig
          lp_data <- lp_data - lyy - .tobs_log_safe(sig) -
                     0.5 * log(2 * pi) - 0.5 * r * r
          g_eta_pos[i, v] <- r / sig
          g_logdisp <- g_logdisp + (r * r - 1)
        }
      }
    } else {
      P0 <- exp(sum(log(1 - sgm(eta_p[i, vv]))))
      L  <- psi[i] * P0 + (1 - psi[i])
      lp_data <- lp_data + log(L)
      g_eta_psi[i] <- -psi[i] * (1 - psi[i]) * (1 - P0) / L
      for (v in vv) g_eta_p[i, v] <- -psi[i] * P0 * sgm(eta_p[i, v]) / L
    }
  }

  # Design-sandwich the eta-gradients onto the coefficient blocks. The visit-level
  # long vector follows site-major order (row i visit v -> (i-1)*J + v), matching
  # the byrow = TRUE design broadcast and the C++ row map i * J + v.
  grad_bo  <- as.numeric(crossprod(model$X_occ, g_eta_psi))
  grad_bp  <- as.numeric(crossprod(model$X_det_site, rowSums(g_eta_p)))
  if (p_det_visit > 0L)
    grad_bp <- c(grad_bp, as.numeric(crossprod(model$X_det_visit, as.numeric(t(g_eta_p)))))
  grad_bpos <- as.numeric(crossprod(model$X_pos_site, rowSums(g_eta_pos)))
  if (p_pos_visit > 0L)
    grad_bpos <- c(grad_bpos, as.numeric(crossprod(model$X_pos_visit, as.numeric(t(g_eta_pos)))))
  grad <- c(grad_bo, grad_bp, grad_bpos, g_logdisp)

  # Field blocks: per block accumulate the per-cell field score -- the psi-arm
  # score plus alpha times the cover-arm score, both weighted by that block's own
  # per-site design column -- and its copy amplitude's own data-score, then hand
  # both to the field backward pass. Each block's gradient (whitened field, then
  # any sampled hyper) appends in block order, as the C++ layout lays them out.
  lp_field <- 0
  if (has_field) {
    g_pos_site <- rowSums(g_eta_pos)
    for (b in seq_along(fvs)) {
      fv <- fvs[[b]]
      w  <- fv$weight %||% rep(1, N)
      g_f_site <- w * (g_eta_psi + fss[[b]]$alpha$value * g_pos_site)
      g_field  <- numeric(fv$n_units)
      for (s in seq_len(N)) {
        cell <- fv$field_map[s]
        g_field[cell] <- g_field[cell] + g_f_site[s]
      }
      bk       <- .ochf_backward(fv, fss[[b]], g_field,
                                 sum(f_site[[b]] * g_pos_site))
      grad     <- c(grad, bk$grad)
      lp_field <- lp_field + bk$lp
    }
  }

  # RE blocks trail the field block, each scoring against its own arm's per-row
  # eta-gradient (site-major, matching `codes`).
  lp_re <- 0
  for (rb in re_view) {
    g_row <- if (rb$arm == 1L) as.numeric(t(g_eta_p)) else as.numeric(t(g_eta_pos))
    bk    <- .ocre_backward(rb, theta, g_row)
    grad  <- c(grad, bk$grad)
    lp_re <- lp_re + bk$lp
  }

  lp  <- lp_data + lp_field + lp_re
  ib2 <- 1 / sigma.beta^2
  nb  <- n_coef
  bv  <- theta[seq_len(nb)]
  lp  <- lp - 0.5 * ib2 * sum(bv^2)
  grad[seq_len(nb)] <- grad[seq_len(nb)] - ib2 * bv
  ild2 <- 1 / sigma.logdisp^2
  lp   <- lp - 0.5 * ild2 * log_disp^2
  grad[total] <- grad[total] - ild2 * log_disp
  list(lp = lp, grad = grad)
}

# Resolve the model's observation-arm RE term designs
# (.occu_cover_obs_re_design, shared with the nested-Laplace path) into the
# sampler's block descriptions: arm code 1 = detection, 2 = positive cover, one
# block per grouping factor in formula order (so crossed / nested groupings are
# simply several blocks). `sigma_lsd` is the prior SD of log sigma_re.
#
# The sampler's block is a scalar per group, so a random SLOPE -- which needs a
# per-row design weight (uncorrelated) or a multivariate free-Sigma block
# (correlated) -- is rejected here rather than silently fitted as an intercept.
# Returns NULL when the model carries no observation-arm RE.
.occu_cover_nuts_re_blocks <- function(model, sigma_lsd = 1.5) {
  arms <- list(list(tag = "p",   code = 1L, terms = model$re_det),
               list(tag = "pos", code = 2L, terms = model$re_pos))
  out <- list()
  for (a in arms) {
    for (d in a$terms %||% list()) {
      if (!identical(d$type, "intercept") || !identical(d$n_coefs, 1L)) {
        stop(sprintf(paste0(
          "occu_cover() method = \"nuts\" samples random INTERCEPT blocks on the ",
          "detection / positive-cover arms: one scalar deviation per level of a ",
          "grouping factor. The %s term `%s` is a random SLOPE, which needs the ",
          "grid-integrated method = \"nested_laplace\" path (weighted iid blocks ",
          "for `||`, a free-Sigma `miid` block for `|`)."),
          a$tag, d$term_label), call. = FALSE)
      }
      out[[length(out) + 1L]] <- list(
        arm = a$code, arm_tag = a$tag, group = as.integer(d$codes_flat),
        n_groups = as.integer(d$n_groups), sigma_lsd = as.numeric(sigma_lsd),
        var = d$var, levels = d$levels)
    }
  }
  if (!length(out)) NULL else out
}

# Public names for the RE coordinates, derived from the SAME arm/term naming the
# nested-Laplace postprocess uses (.occu_cover_re_sigma_names): a lone term on an
# arm keeps the bare `sigma_re_p` / `sigma_re_pos`, crossed / nested terms
# sharing an arm are disambiguated by the grouping variable. The sampled
# coordinate is the log SD, and the whitened deviations reuse the same stem.
.occu_cover_nuts_re_names <- function(blocks) {
  sig  <- .occu_cover_re_sigma_names(
    lapply(blocks, function(b) list(arm = b$arm_tag, var = b$var)))
  stem <- sub("^sigma_re", "re", sig)
  unlist(lapply(seq_along(blocks), function(i)
    c(paste0(stem[i], "_z", seq_len(blocks[[i]]$n_groups)),
      paste0("log_", sig[i]))))
}

# Per-(site, visit) RE offset matrices at a given set of block deviations
# `b_list` (one numeric vector of per-group offsets per block), on the padded
# [n_sites x max_visits] grid the arm predictors live on. Returns one matrix per
# observation arm, so the fit's log-likelihood is evaluated at the predictor the
# sampler ran on.
.occu_cover_nuts_re_offsets <- function(blocks, b_list, N, J) {
  z <- matrix(0, N, J)
  out <- list(p = z, pos = z)
  for (i in seq_along(blocks)) {
    rb <- blocks[[i]]
    off <- ifelse(rb$group > 0L, b_list[[i]][pmax(rb$group, 1L)], 0)
    tag <- rb$arm_tag
    out[[tag]] <- out[[tag]] + matrix(off, N, J, byrow = TRUE)
  }
  out
}

# Build the C++ NUTS spec list from a bound non-spatial occu_cover model. NULL
# visit-level designs become explicit 0-column matrices so the C++ side reads a
# uniform layout. The matrices are passed straight through (raw design scale; the
# occu_cover path is not autoscaled), so the draws land on the natural scale.
.tobs_occu_cover_nuts_spec <- function(model) {
  N <- model$n_sites; J <- model$max_visits
  empty_visit <- function(X) if (is.null(X)) matrix(0, N * J, 0L) else X
  list(
    n_sites     = as.integer(N),
    max_visits  = as.integer(J),
    pos_code    = .occu_cover_pos_code(model$positive),
    y           = matrix(as.integer(model$y), N, J),
    y_pos       = matrix(as.numeric(model$y_pos), N, J),
    valid       = matrix(as.integer(model$valid), N, J),
    X_occ       = model$X_occ,
    X_det_site  = model$X_det_site,
    X_det_visit = empty_visit(model$X_det_visit),
    X_pos_site  = model$X_pos_site,
    X_pos_visit = empty_visit(model$X_pos_visit)
  )
}


# ---------------------------------------------------------------------------
# Front-door NUTS fitter for the non-spatial joint occupancy + cover hurdle
# ---------------------------------------------------------------------------

# Sample the exact non-spatial occu_cover coefficient posterior via tulpa's NUTS
# engine and the in-tree C++ FullGradFn (cpp_occu_cover_nuts), warm-started at the
# Laplace mode with a diagonal Laplace metric, then package the draws into the
# same tobs_fit shape .tobs_fit_occu_cover returns so coef / vcov / confint /
# predict / WAIC read the NUTS posterior. The occu_cover path is not autoscaled,
# so the returned draws / means / vcov are already on the natural coefficient
# scale.
#
# An observation-arm random intercept rides the same target as one non-centered
# block per grouping factor, with its group SD SAMPLED rather than pinned or
# grid-integrated. `sigma.logdisp` and `sigma.logre` are internal weak-prior
# widths (no control knob, like abun's sigma.logr): the log-dispersion and each
# log group SD. `...` absorbs unused sampler controls (progress.*).
.tobs_fit_occu_cover_nuts <- function(model, priors = NULL,
                                      sigma.beta = NULL, sigma.logdisp = 5,
                                      sigma.logre = 1.5,
                                      n.iter = NULL, n.warmup = NULL,
                                      n.chains = NULL, n.thin = NULL,
                                      n.threads = NULL, max.treedepth = NULL,
                                      adapt.delta = NULL, seed = NULL,
                                      verbose = FALSE, ...) {
  # Sampler defaults come from the one engine table.
  .tobs_fill_sampler(environment(), "nuts")

  pin   <- model$process_info
  p_occ <- pin[[1L]]$p; p_p <- pin[[2L]]$p; p_pos <- pin[[3L]]$p
  n_par <- p_occ + p_p + p_pos + 1L
  N <- model$n_sites; J <- model$max_visits

  par_names <- c(
    paste0("psi_", pin[[1L]]$coef_names),
    paste0("p_",   pin[[2L]]$coef_names),
    paste0("pos_", pin[[3L]]$coef_names),
    if (identical(model$positive, "beta")) "log_phi" else "log_sigma_pos"
  )

  # Observation-arm random intercepts: one non-centered block per grouping
  # factor on the detection / positive-cover arm, each with its own SAMPLED log
  # sigma_re. Resolved before the warm fit so an unsupported term errors
  # immediately.
  re_blocks <- .occu_cover_nuts_re_blocks(model, sigma_lsd = sigma.logre)
  re_names  <- if (is.null(re_blocks)) NULL else .occu_cover_nuts_re_names(re_blocks)
  n_re_coef <- if (is.null(re_blocks)) 0L
               else sum(vapply(re_blocks, function(b) b$n_groups + 1L, numeric(1)))

  # Warm start at the Laplace mode + a diagonal Laplace metric from its vcov.
  # The RE coordinates start whitened at 0 with a moderate sigma_re, the
  # convention the observation-family samplers share (nuts_chains.R).
  warm <- .tobs_fit_occu_cover(model, method = "laplace", priors = priors,
                               sigma.beta = sigma.beta, verbose = FALSE)
  theta0 <- as.numeric(warm$means)
  V <- as.matrix(warm$vcov)
  inv_metric <- if (!is.null(V) && all(dim(V) == n_par) && all(is.finite(diag(V))))
                  pmax(diag(V), 1e-6) else rep(1, n_par)
  for (b in re_blocks %||% list()) {
    theta0     <- c(theta0, rep(0, b$n_groups), log(0.5))
    inv_metric <- c(inv_metric, rep(1, b$n_groups), 0.25)
  }

  spec <- .tobs_occu_cover_nuts_spec(model)
  if (!is.null(re_blocks)) {
    spec$re_blocks <- lapply(re_blocks, function(b)
      list(re_arm = b$arm, re_group = b$group, n_re_groups = b$n_groups,
           sigma_re_lsd = b$sigma_lsd))
  }

  run_chain <- function(ch) {
    cpp_occu_cover_nuts(
      spec, theta0 = theta0, sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp,
      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
      n_warmup = as.integer(n.warmup), max_treedepth = as.integer(max.treedepth),
      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
      verbose = isTRUE(verbose) && ch == 1L)
  }
  n_chains <- max(1L, as.integer(n.chains))
  chains   <- lapply(.tobs_nuts_run_parallel(run_chain, n_chains, n.threads),
                     .tobs_nuts_thin_chain, n.thin = n.thin)
  per_chain_draws <- lapply(chains, `[[`, "draws")
  all_draws <- do.call(rbind, per_chain_draws)
  colnames(all_draws) <- c(par_names, re_names)
  # The reported draw matrix is the coefficient block: `coef` / `vcov` / the WAIC
  # substrate read it positionally, with the trailing column the log-dispersion.
  # The RE coordinates are summarised on fit$re instead (fit$re_draws keeps them
  # per draw), the same split the spatial sampler makes for its field.
  draws    <- all_draws[, seq_len(n_par), drop = FALSE]
  n_draws  <- nrow(draws)

  means  <- colMeans(draws); names(means) <- par_names
  V_post <- stats::cov(draws); dimnames(V_post) <- list(par_names, par_names)
  sds    <- sqrt(pmax(diag(V_post), 0)); names(sds) <- par_names

  # Per-block RE posterior: sigma_re = exp(log_sigma_re) and the per-group
  # deviations b_g = sigma_re * z_g, on the natural scale. A positive variance
  # component at a few dozen groups is right-skewed, so the median is reported
  # alongside the mean as the summary to quote against a known truth.
  re_out <- NULL
  if (!is.null(re_blocks)) {
    k <- n_par
    re_out <- lapply(seq_along(re_blocks), function(i) {
      rb  <- re_blocks[[i]]
      idx <- k + seq_len(rb$n_groups); k <<- k + rb$n_groups + 1L
      sig <- exp(all_draws[, k])
      blup <- sig * all_draws[, idx, drop = FALSE]
      list(arm = rb$arm_tag, var = rb$var, levels = rb$levels,
           n_groups = rb$n_groups, n_coefs = 1L, coef_names = "(Intercept)",
           correlated = FALSE,
           sigma = mean(sig), sigma_median = stats::median(sig),
           sigma_sd = stats::sd(sig),
           sigma_draws = as.numeric(sig),
           blup = colMeans(blup), blup_sd = apply(blup, 2L, stats::sd))
    })
    names(re_out) <- .occu_cover_re_keys(
      vapply(re_blocks, `[[`, character(1), "arm_tag"),
      vapply(re_blocks, function(b) as.character(b$var), character(1)))
  }

  # Data log-likelihood at the posterior mean (scale-invariant), so logLik() on
  # the NUTS fit matches the laplace-path convention. Reuses the shared marginal,
  # evaluated at the predictor the sampler ran on: the RE offsets enter the
  # detection arm on its logit scale and the cover arm on its link scale.
  bo   <- means[seq_len(p_occ)]
  bp   <- means[p_occ + seq_len(p_p)]
  bpos <- means[p_occ + p_p + seq_len(p_pos)]
  eta  <- .occu_cover_eta_from_par(model, bo, bp, bpos)
  p_mat <- eta$p_mat; ep_mat <- eta$ep_mat
  if (!is.null(re_blocks)) {
    off <- .occu_cover_nuts_re_offsets(
      re_blocks, lapply(re_out, `[[`, "blup"), N, J)
    p_mat  <- stats::plogis(.tobs_clamp_eta(stats::qlogis(p_mat) + off$p))
    ep_mat <- ep_mat + off$pos
  }
  ll_mean <- sum(.occu_cover_site_ll(model, eta$psi, p_mat, ep_mat, means[n_par]))

  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- as.integer(unlist(lapply(chains, `[[`, "divergent")))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- mean(vapply(chains, function(ch) ch$epsilon %||% NA_real_, numeric(1)),
                    na.rm = TRUE)

  nuts <- list(accept_prob = accept, divergent = divergent, treedepth = treedepth,
               epsilon = epsilon, n_chains = n_chains, divergent_total = sum(divergent),
               sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp)

  # Per-group SD convergence over the natural-scale sigma_re draws, chain by
  # chain (exp is monotone, so this is the split-Rhat of the sampled coordinate).
  if (!is.null(re_blocks)) {
    ends   <- cumsum(vapply(per_chain_draws, nrow, integer(1)))
    starts <- c(1L, ends[-length(ends)] + 1L)
    sig_ch <- lapply(seq_len(n_chains), function(ch)
      do.call(cbind, lapply(re_out, function(r)
        r$sigma_draws[starts[ch]:ends[ch]])))
    dg <- .tobs_nuts_rhat_ess(sig_ch)
    names(dg$rhat) <- names(dg$ess) <- names(re_out)
    nuts$sigma_logre  <- sigma.logre
    nuts$re_sigma     <- vapply(re_out, `[[`, numeric(1), "sigma")
    nuts$re_sigma_rhat <- dg$rhat
    nuts$re_sigma_ess  <- dg$ess
  }

  fit <- structure(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V_post,
    n_samples    = n_draws,
    n_params     = n_par,
    log_prob     = rep(ll_mean, n_draws),
    log_lik      = ll_mean,
    N            = sum(model$valid),
    accept_prob  = accept,
    divergent    = divergent,
    treedepth    = treedepth,
    epsilon      = epsilon,
    col_names    = par_names,
    param_names  = par_names,
    process_info = pin,
    model        = model,
    spatial      = NULL,
    re           = re_out,
    re_draws     = if (is.null(re_blocks)) NULL
                   else all_draws[, n_par + seq_len(n_re_coef), drop = FALSE],
    method       = "nuts",
    positive     = model$positive,
    nuts         = nuts,
    convergence  = list(converged = NA, n_iter = as.integer(n.iter))
  ), class = c("tobs_fit", "tulpa_fit"))

  # Per-parameter split-R-hat / bulk + tail ESS, through the writer every sampled
  # path shares, so summary.tobs_fit surfaces them per parameter; `fit$nuts`
  # carries the same two vectors alongside the sampler diagnostics. Reported over
  # the coefficient block, the coordinates `summary()` puts on its rows.
  fit <- .tobs_nuts_attach_convergence(fit, per_chain_draws,
                                       par_names = par_names,
                                       cols = seq_len(n_par))
  fit$nuts$rhat <- fit$convergence$rhat
  fit$nuts$ess  <- fit$convergence$ess_bulk
  fit
}


# ---------------------------------------------------------------------------
# Spatial occu_cover NUTS: fixed-hyper non-centered coupled proper-CAR field
# ---------------------------------------------------------------------------

# Desugar a varying-coefficient bar (`spatial(~ 1 || node, graph = adj)`) on the
# psi formula into the plain areal field spec the NUTS sampler takes, through the
# SAME expansion the nested-Laplace path uses (.tobs_expand_spatial_bar): one
# unweighted intercept field plus one weight-scaled field per bar covariate
# column, each identical to what `icar(graph = adj, group_var = node)` /
# `icar(graph = adj, weight = col, group_var = node)` builds. A single-column bar
# therefore IS the plain areal term, and routes unchanged.
#
# The sampler carries one block PER FIELD -- each with its own loading, site ->
# node map, per-site design weight and copy amplitude -- so a bar declaring an
# intercept field plus varying-coefficient field(s) expands to that many blocks. A
# correlated bar (`|`) is a different structure (one free-Sigma MCAR block across
# the fields) and is not sampled here.
.occu_cover_nuts_bar_fields <- function(spec, data) {
  if (!is.null(spec$by_var)) {
    stop(paste0(
      "occu_cover() NUTS + areal spatial samples one field over one graph; a ",
      "replicated field (spatial(<bar>, by = \"", spec$by_var, "\")) needs ",
      "method = \"nested_laplace\"."), call. = FALSE)
  }
  if (isTRUE(spec$correlated)) {
    stop(paste0(
      "occu_cover() NUTS + areal spatial samples INDEPENDENT areal fields, one ",
      "block each (`||`); a correlated bar (`|`) is one free-Sigma MCAR block ",
      "across the fields and needs method = \"nested_laplace\"."), call. = FALSE)
  }
  .tobs_expand_spatial_bar(spec, data)
}

# Resolve the single spatial term + fixed-effect psi formula for the NUTS spatial
# path. Unlike .occu_cover_spatial_fields (which gates to icar/bym2 for the
# grid-integrated nested-Laplace engine), the NUTS path also accepts car_proper()
# (the full-rank precision the fixed-hyper non-centered field is best conditioned
# on) and rejects temporal / RE terms with a pointer to the nested-Laplace route.
# The bar form is desugared first, so `spatial(~ 1 || cell, graph = adj)` and
# `icar(graph = adj, group_var = "cell")` reach the sampler as one field
# description, and `spatial(~ 1 + w || cell, graph = adj)` as the intercept field
# plus one varying-coefficient field per covariate column. Returns NULL when the
# psi formula carries no spatial term (the non-spatial NUTS sampler), or list(fe,
# spatial, group_var, trends) with `spatial` the intercept field and `trends` the
# weighted ones.
.occu_cover_nuts_spatial_term <- function(formula, data) {
  bind <- .tobs_bind_formulas(list(psi = formula), data)
  if (length(bind$terms) == 0L) return(NULL)
  spatial <- Filter(function(t) inherits(t$spec, "tobs_spatial"), bind$terms)
  re_terms <- Filter(function(t) inherits(t$spec, "tobs_re"), bind$terms)
  other   <- Filter(function(t) !inherits(t$spec, "tobs_spatial") &&
                                !inherits(t$spec, "tobs_re"), bind$terms)
  if (length(spatial) == 0L) return(NULL)
  if (length(re_terms) > 0L || length(other) > 0L) {
    stop(paste0(
      "occu_cover() NUTS + areal spatial samples coupled areal field(s) on the ",
      "psi formula; a per-group RE / temporal term composes only on the ",
      "grid-integrated method = \"nested_laplace\" path."), call. = FALSE)
  }
  # Every areal term on the formula, bars desugared, in declaration order.
  specs <- unlist(lapply(spatial, function(t) {
    s <- t$spec
    if (isTRUE(s$is_bar)) .occu_cover_nuts_bar_fields(s, data) else list(s)
  }), recursive = FALSE)
  if (any(vapply(specs, function(s) isTRUE(s$is_multifield), logical(1)))) {
    stop(paste0(
      "occu_cover() NUTS + areal spatial takes each coupled field as its own ",
      "term or bar column; a pre-combined multi-field container needs ",
      "method = \"nested_laplace\"."), call. = FALSE)
  }
  weighted <- vapply(specs, function(s) !is.null(s$weight), logical(1))
  base <- specs[!weighted]
  if (length(base) != 1L) {
    stop(sprintf(paste0(
      "occu_cover() NUTS + areal spatial needs exactly one unweighted intercept ",
      "field beside any varying-coefficient field(s); got %d."), length(base)),
      call. = FALSE)
  }
  trends <- specs[weighted]
  # Every coupled field shares one node set: the sampler holds one graph and one
  # site -> node map, and the varying-coefficient fields differ only in weight.
  base_graph <- base[[1L]]$graph
  for (s in trends) {
    if (!identical(dim(s$graph), dim(base_graph)) || !all(s$graph == base_graph))
      stop("occu_cover() NUTS coupled fields must share the same areal graph ",
           "as the intercept field (same nodes / adjacency).", call. = FALSE)
    if (!identical(s$type, base[[1L]]$type))
      stop("occu_cover() NUTS coupled fields must share one areal kind; got ",
           base[[1L]]$type, "() and ", s$type, "().", call. = FALSE)
  }
  gvs <- unique(Filter(Negate(is.null), lapply(specs, `[[`, "group_var")))
  if (length(gvs) > 1L)
    stop("occu_cover() NUTS coupled fields must share a single group_var ",
         "(or none).", call. = FALSE)
  list(fe = bind$fe$psi, spatial = base[[1L]],
       group_var = if (length(gvs) == 1L) gvs[[1L]] else NULL,
       trends = lapply(trends, function(s) list(
         weight = as.numeric(s$weight),
         label  = s$weight_label %||% s$component %||% "trend")))
}


# Run the nested-Laplace joint engine once with a proper-CAR copy block
# to obtain the FIXED field hyperparameters (sigma, rho_car) and the copy
# amplitude (alpha), plus the grid-weighted posterior-mean coupled field. This is
# the occu_cover analogue of the warm fit .tobs_fit_abun_nuts_spatial reads from
# nmix_laplace_car_proper: the same proper-CAR precision tau Q(rho) the NUTS field
# block then fixes, here estimated through the shared two-state cell-coupling
# marginal (the occu_cover_{lognormal,beta} spec) so the field hyper sits at the
# nested-Laplace estimate of THIS model (the marginalize-then-fix recipe). Returns
# the betas, log-dispersion, (sigma, rho, alpha), and field f at its grid-weighted
# posterior mean.
#
# `trends` are the varying-coefficient fields beside the intercept field, each a
# per-site weight column. With any of them the warm fit is the MULTI-BLOCK
# coupled path -- one areal block per field, each carrying its per-arm design
# weight and its own copy amplitude axis -- which is the same structure the
# nested-Laplace SVC route fits (R/occu_cover_joint.R), so the sampler's prior
# support and its starting values come from the grid the deterministic backend
# would have integrated.
.tobs_occu_cover_nuts_carproper_warm <- function(model, adj, priors,
                                                 max.iter = 200L, tol = 1e-6,
                                                 type = "car_proper",
                                                 sigma.grid = NULL,
                                                 rho.car.grid = NULL,
                                                 alpha.grid = NULL,
                                                 alpha.n = NULL,
                                                 trends = list(),
                                                 alpha.grid.trend = NULL,
                                                 alpha.n.trend = NULL,
                                                 copy.slab = NULL,
                                                 copy.atom.mass = NULL) {
  is_beta   <- identical(model$positive, "beta")
  spec_name <- switch(model$positive,
                      beta     = "occu_cover_beta",
                      gaussian = "occu_cover_gaussian",
                      "occu_cover_lognormal")
  n_sites   <- model$n_sites
  site_cell <- model$site_cell %||% seq_len(n_sites)
  n_cells   <- nrow(adj)

  # Pre-fit the pos-arm dispersion at the empirical cover spread (matching the
  # joint non-latent path); it rides the spec's phi slot, fixed here.
  # Observed covers only (a detected visit may carry a missing cover).
  pos_vals <- model$y_pos[model$valid & model$y == 1L]
  pos_vals <- pos_vals[is.finite(pos_vals)]
  sigma_pos_init <- if (is_beta) {
    if (length(pos_vals) >= 2L) {
      mu_hat <- mean(pos_vals); var_hat <- max(stats::var(pos_vals), 1e-6)
      max((mu_hat * (1 - mu_hat)) / var_hat - 1, 1)
    } else 10
  } else if (identical(model$positive, "gaussian")) {
    # Residual SD on the raw response scale (no log; the response may be negative).
    if (length(pos_vals) >= 2L) max(stats::sd(pos_vals), 0.05) else 1
  } else {
    if (length(pos_vals) > 0L) max(stats::sd(log(pos_vals)), 0.05) + 0.05 else 0.4
  }

  # The sampler reads the axis's NODES -- their span is the support of alpha's
  # flat prior -- so the resolution knob is resolved to nodes here rather than
  # handed to the engine as `alpha_n`, the way the grid-integrated routes do.
  alpha_axis <- .tobs_alpha_axis(alpha.grid, alpha.n)
  alpha_grid <- .tobs_alpha_nodes(alpha_axis)
  sigma_grid <- sigma.grid %||% .tobs_default_sigma_grid()
  rho_car_grid <- rho.car.grid %||% .tobs_default_rho_car_grid()

  n_trend <- length(trends)
  multi   <- n_trend > 0L
  if (multi) {
    # Each coupled field adds its own (sigma, alpha[, rho]) axes, so the outer
    # tensor is the PRODUCT over blocks, and the engine caps it. The sampler
    # reads each axis's SPAN as the support of that hyper's flat prior, and
    # nothing else off the grid, so a DEFAULTED axis is thinned to three nodes
    # over the SAME span: the prior is unchanged and the warm fit stays
    # affordable. Provenance is the auto-grid mark, not whether the argument
    # arrived -- `share(spatial())` with no amplitude hands this path the default
    # alpha axis through `control$alpha.grid`, and that is a defaulted axis. An
    # axis the caller chose keeps its nodes.
    thin <- function(g) {
      if (!tulpa::is_auto_grid(g)) return(g)
      v   <- sort(unique(as.numeric(g)))
      pos <- v[v != 0]
      if (length(pos) > 3L)
        pos <- pos[unique(round(seq(1, length(pos), length.out = 3L)))]
      .tobs_mark_auto(c(v[v == 0], pos), TRUE)
    }
    alpha_grid   <- thin(alpha_grid)
    sigma_grid   <- thin(sigma_grid)
    rho_car_grid <- thin(rho_car_grid)
    # bym2's mixing axis is defaulted by the engine rather than passed, so it is
    # named here to be thinned on the same terms -- over its own span, so the
    # sampled rho's prior support is the one the engine would have integrated.
    bym2_rho_grid <- if (identical(type, "bym2"))
      thin(tulpa::auto_grid(tulpa:::.nl_grid_axis("bym2_rho"))) else NULL
  }
  # Single-block (multi = FALSE): the pos arm carries
  # field_coef = list(name = "alpha", grid = alpha_grid), so the copy alpha axis
  # rides the one shared field. The single-block joint path takes the copy
  # coefficient on the arm, not a top-level `copy` block. With a varying-
  # coefficient field the fit takes the multi-block driver instead, where each
  # block's amplitude is an explicit copy spec.
  arms_out <- .occu_cover_build_joint_arms(
    model = model, sigma_pos_init = sigma_pos_init,
    alpha_axis = .tobs_alpha_axis(grid = alpha_grid),
    positive = model$positive, multi = multi, n_cells = n_cells,
    site_cell = site_cell, cover_aggregate = "none")
  responses <- arms_out$responses

  arm_priors <- .occu_cover_coupled_arm_priors(priors, responses)
  for (nm in c("psi", "p", "pos")) {
    ap <- arm_priors[[nm]]
    if (!is.null(ap)) {
      responses[[nm]]$beta_prior_mean <- ap$mean
      responses[[nm]]$beta_prior_prec <- ap$prec
    }
  }

  csr <- .occu_cover_adj_to_csr(adj)
  block_template <- function(sigma_g, extra = list()) {
    p <- c(list(
      type            = type,
      n_spatial_units = csr$n_spatial_units,
      adj_row_ptr     = csr$adj_row_ptr,
      adj_col_idx     = csr$adj_col_idx,
      n_neighbors     = csr$n_neighbors,
      sigma_grid      = sigma_g), extra)
    # Only the proper-CAR field carries a spatial-correlation rho grid; the
    # intrinsic icar / bym2 fields fix rho (icar rho = 1, bym2 mixing gridded by
    # the engine's own bym2 axis, thinned above when several blocks multiply it).
    if (identical(type, "car_proper")) p$rho_car_grid <- rho_car_grid
    if (multi && identical(type, "bym2")) p$rho_grid <- bym2_rho_grid
    p
  }

  copy_arg <- NULL
  if (!multi) {
    prior_arg <- block_template(sigma_grid, list(
      spatial_idx = lapply(responses, function(a) as.integer(a$spatial_idx))))
  } else {
    # One areal block per coupled field, all on the same graph. The per-arm field
    # node is the site's cell (psi rows are sites, pos rows are visits behind a
    # site); the detection arm is excluded by its field_coef = 0. `svc_weight`
    # carries the per-row design weight, which is what makes a block a varying
    # coefficient; each block's cross-arm amplitude is its own copy spec.
    pos_site     <- arms_out$pos_site
    n_v          <- arms_out$n_visits_valid
    spatial_idx  <- list(as.integer(site_cell), arms_out$cell_of_visit,
                         as.integer(site_cell[pos_site]))
    field_block  <- function(w_site, sigma_g) {
      w <- if (is.null(w_site)) rep(1.0, n_sites) else as.numeric(w_site)
      block_template(sigma_g, list(
        spatial_idx = spatial_idx,
        svc_weight  = list(w, rep(1.0, n_v), w[pos_site])))
    }
    alpha_axis_trend <- .tobs_alpha_axis_trend(
      list(alpha.grid.trend = alpha.grid.trend, alpha.n.trend = alpha.n.trend),
      alpha_axis)
    alpha_grid_trend <- thin(.tobs_alpha_nodes(alpha_axis_trend))
    prior_arg <- c(list(field_block(NULL, sigma_grid)),
                   lapply(trends, function(tf)
                     field_block(tf$weight, sigma_grid)))
    copy_arg <- c(
      list(list(arm = "pos", block = 1L, alpha_grid = alpha_grid)),
      lapply(seq_len(n_trend), function(j)
        list(arm = "pos", block = j + 1L, alpha_grid = alpha_grid_trend)))
  }

  # Largest positive node of the alpha axis each block hands the engine. The
  # copy scale's declared slab is anchored on it (tulpa's
  # .hyper_copy_slab_density), so the sampler reads its rate off the same node
  # set the outer grid's density is declared over.
  alpha_upper <- function(b) {
    g <- as.numeric(if (multi && b > 1L) alpha_grid_trend else alpha_grid)
    pos <- g[is.finite(g) & g > 0]
    if (length(pos) == 0L) NA_real_ else max(pos)
  }

  fit_call <- list(
    responses = responses, prior = prior_arg, cell_coupling = spec_name,
    control = list(max_iter = as.integer(max.iter), tol = as.numeric(tol),
                   n_threads = 1L, store_Q = FALSE, adaptive_grid = FALSE,
                   var_of_means_consistency = FALSE, diagnose_k = FALSE,
                   progress = FALSE))
  # The sampler takes each axis's realised span as the support of that hyper's
  # flat prior, so the design has to BE a tensor of the axis nodes. A second
  # field takes the axis count past three, where the engine would otherwise
  # switch to a mode-centred CCD star whose column range is a design radius
  # rather than an integrated span -- and can sit outside the nodes entirely.
  # One field never reaches that threshold, so the request is made only where it
  # changes something.
  if (multi) fit_call$control$integration <- "grid"
  # A declared copy-scale slab is the measure BOTH backends integrate, so the
  # warm fit's outer grid carries the shape the sampler is given. The point mass
  # at alpha = 0 is the warm fit's alone -- a gradient sampler cannot visit it --
  # so its declared share rides here and nowhere else.
  if (!is.null(copy.slab)) fit_call$control$copy_slab <- copy.slab
  if (!is.null(copy.atom.mass))
    fit_call$control$copy_atom_mass <- copy.atom.mass
  if (!is.null(copy_arg)) fit_call$copy <- copy_arg
  fit <- do.call(tulpa::tulpa_nested_laplace_joint, fit_call)

  # Through the shared reconciliation so the warm start is read against the same
  # outer-grid measure the fit integrated, declared slab and node spacing
  # included, rather than a flat weight per cell.
  oc  <- .tobs_joint_ok_cells(fit, "occu_cover NUTS spatial: warm car_proper fit")
  ok  <- oc$ok_cells
  w   <- oc$w
  fit <- oc$fit

  layout <- fit$arm_layout
  p_psi  <- layout$p[1L]; p_p <- layout$p[2L]; p_pos <- layout$p[3L]
  bpsi_idx <- layout$beta_start[1L] + seq_len(p_psi)
  bp_idx   <- layout$beta_start[2L] + seq_len(p_p)
  bpos_idx <- layout$beta_start[3L] + seq_len(p_pos)
  f0       <- (layout$field_starts %||% layout$phi_start)
  modes <- fit$modes[ok, , drop = FALSE]

  beta_psi <- as.numeric(crossprod(w, modes[, bpsi_idx, drop = FALSE]))
  beta_p   <- as.numeric(crossprod(w, modes[, bp_idx,   drop = FALSE]))
  beta_pos <- as.numeric(crossprod(w, modes[, bpos_idx, drop = FALSE]))

  tg     <- fit$theta_grid[ok, , drop = FALSE]
  # Block b's axis is named bare on the single-block grid and `b<k>.` prefixed on
  # the multi-block one.
  pick   <- function(nm, b) {
    j <- match(sprintf("b%d.%s", b, nm), colnames(fit$theta_grid))
    if (is.na(j) && b == 1L) j <- match(nm, colnames(fit$theta_grid))
    if (is.na(j)) return(NA_real_)
    sum(w * as.numeric(tg[, j]))
  }
  blocks <- lapply(seq_len(1L + n_trend), function(b) {
    field_idx <- f0[[b]] + seq_len(n_cells)
    fld <- as.numeric(crossprod(w, modes[, field_idx, drop = FALSE]))
    list(field = fld - mean(fld),                 # sum-to-zero convention
         sigma = pick("sigma", b), alpha = pick("alpha", b),
         rho   = if (identical(type, "car_proper")) pick("rho_car", b)
                 else if (identical(type, "bym2")) pick("rho", b) else 1.0,
         alpha_upper = alpha_upper(b),
         weight = if (b == 1L) NULL else trends[[b - 1L]]$weight,
         label  = if (b == 1L) "intercept" else trends[[b - 1L]]$label)
  })

  b1 <- blocks[[1L]]
  list(beta_psi = beta_psi, beta_p = beta_p, beta_pos = beta_pos,
       log_disp = log(sigma_pos_init), field = b1$field, type = type,
       sigma = b1$sigma, rho = b1$rho, alpha = b1$alpha,
       alpha_upper = b1$alpha_upper,
       blocks = blocks, joint_fit = fit)
}


# Fixed basis B1 of a coupled areal field, together with the rho-scaling its
# columns carry. Every areal kind factors into a basis that does NOT depend on
# the sampled hypers:
#   icar / bym2 structured  the sum-to-zero eigen-loading of the intrinsic
#                           precision Q at unit precision, so sigma is a
#                           scalar multiply and bym2's rho a scalar
#                           re-weight against the unstructured block.
#   car_proper              Q(rho) = D - rho W = D^{1/2}(I - rho Lambda)D^{1/2}
#                           in the eigenbasis of the symmetrically normalised
#                           adjacency D^{-1/2} W D^{-1/2} = U Lambda U'. Hence
#                           B1 = D^{-1/2} U is fixed and rho only rescales the
#                           column weights (1 - rho lambda_j)^{-1/2}: no
#                           per-leapfrog Cholesky of an irregular graph.
# D and W follow .areal_Q: W = adj, D = diag(number of non-zero neighbours).
.occu_cover_nuts_field_basis <- function(adj, type, n, scale_factor = NULL) {
  if (identical(type, "car_proper")) {
    deg <- rowSums(adj != 0)
    if (any(deg <= 0))
      stop("occu_cover NUTS spatial: the graph has an isolated node, so the ",
           "proper-CAR precision is singular.", call. = FALSE)
    dm12 <- 1 / sqrt(deg)
    M    <- adj * outer(dm12, dm12)
    ev   <- eigen((M + t(M)) / 2, symmetric = TRUE)
    return(list(B1 = dm12 * ev$vectors, scale1 = .OCHF_CAR, has_iid = FALSE,
                sf = 1, lambda = ev$values))
  }
  Lstr <- .tobs_field_load(adj, "icar", 1, 1, n)
  if (identical(type, "bym2"))
    return(list(B1 = Lstr, scale1 = .OCHF_BYM2_STR, has_iid = TRUE,
                sf = scale_factor %||% .bym2_scale(adj),
                lambda = numeric(0)))
  list(B1 = Lstr, scale1 = .OCHF_CONST, has_iid = FALSE, sf = 1,
       lambda = numeric(0))
}

# Bounds for each sampled hyper, read off the warm nested-Laplace fit's OWN
# outer grid: the span its nodes reach in the coordinate the axis is spaced in
# (log for sigma / alpha, logit for rho). The nodes are where the deterministic
# backend evaluates that axis, so bounding the sampled coordinate to their span
# holds both backends over the same region and makes the same
# `control$sigma.grid` / `alpha.grid` / `rho.car.grid` knob move both. The
# measure ON that span is the axis's declared prior, which the block builder
# attaches. An axis the grid pinned to one node (or that the fit does not carry)
# returns NULL and the hyper stays pinned.
.occu_cover_nuts_hyper_bounds <- function(warm, type, block = 1L) {
  tg <- warm$joint_fit$theta_grid
  # The single-block warm grid names its axes bare (`sigma`); the multi-block one
  # prefixes each block (`b2.sigma`), so a block reads its own axis under either
  # spelling.
  sup <- warm$joint_fit$axis_support
  ax <- function(nm, positive_only = FALSE) {
    if (is.null(tg) || is.null(colnames(tg))) return(NULL)
    a <- sprintf("b%d.%s", block, nm)
    j <- match(a, colnames(tg))
    if (is.na(j) && block == 1L) { a <- nm; j <- match(a, colnames(tg)) }
    if (is.na(j)) return(NULL)
    v <- as.numeric(tg[, j]); v <- v[is.finite(v)]
    if (positive_only) v <- v[v > 0]
    if (length(v) < 2L) return(NULL)
    # The span the outer quadrature integrates over, which reaches half a node
    # step past the outermost node rather than stopping on it. The fit records
    # it; recompute from the settled grid through the same helper when an older
    # fit object does not carry one.
    r <- sup[[a]]
    if (is.null(r)) {
      specs <- tulpa:::.joint_axis_specs_from_grid(tg)
      k <- match(a, vapply(specs, function(z) z$name, character(1)))
      if (is.na(k)) return(NULL)
      r <- tulpa:::.hyper_axis_support(v, specs[[k]])
    }
    if (is.null(r)) return(NULL)
    r <- sort(as.numeric(r))
    if (!all(is.finite(r)) || r[1L] <= 0 && positive_only) return(NULL)
    if (r[2L] <= r[1L] * (1 + 1e-8)) return(NULL)
    # The node range travels with the span so a bounded axis can fall back to it
    # (`.ochf_rho_support()`); it is an attribute rather than a second element so
    # the span itself stays the declared support this function returns.
    attr(r, "nodes") <- range(v)
    r
  }
  list(sigma = ax("sigma"),
       alpha = ax("alpha", positive_only = TRUE),
       rho   = switch(type, bym2 = ax("rho"), car_proper = ax("rho_car"), NULL))
}

# The span the sampled mixing parameter is given: what the outer grid declared
# for its axis, intersected with the interval the field is defined on.
#
# A rho axis is BOUNDED -- the bym2 mixing weight lives in (0, 1), the
# proper-CAR correlation in (0, 1 / lambda_max) for the largest eigenvalue of
# the symmetrically normalised adjacency -- and the declared span can leave it.
# tulpa classifies every `rho*` axis as evenly spaced, so its support reaches
# half a node step past the outermost node on the NATURAL scale, while the
# proper-CAR default nodes (0.5, 0.8, 0.95, 0.99) are laid out logit-spaced:
# that step lands at 1.01, which is not a correlation (gcol33/tulpa#657).
#
# Where the span leaves the interval it falls back to the outermost NODE, the
# last value the outer quadrature evaluated. Taking the interval's own open edge
# instead is what a clamp to 1 - 1e-4 does, and the proper-CAR loading scales
# its j-th column by (1 - rho lambda_j)^(-1/2): at rho lambda_max = 1 - 1e-4 the
# field's leading direction carries a scale of 100, so the sampler has a
# near-intrinsic funnel inside its own support. Measured on the coupled
# fixture, the difference is in NOTES_measurements.md.
.ochf_rho_support <- function(bounds, lambda) {
  if (is.null(bounds)) return(NULL)
  r  <- sort(as.numeric(bounds))
  nd <- sort(as.numeric(attr(bounds, "nodes") %||% r))
  lam_max <- if (length(lambda)) max(lambda, 0) else 0
  hi <- if (lam_max > 0) 1 / lam_max else 1
  if (!(r[2L] < hi)) r[2L] <- min(nd[2L], r[2L])
  if (!(r[1L] > 0))  r[1L] <- max(nd[1L], r[1L])
  if (!(r[1L] > 0 && r[2L] < hi && r[2L] > r[1L])) return(NULL)
  r
}

# Assemble the field spec entries the C++ block (and the R oracle) read, plus the
# warm-start values of the sampled hyper coordinates. `sample_hyper = FALSE`
# reproduces the fixed-hyper block: the loading built at the warm estimate, with
# sigma and rho baked into its columns and alpha a constant.
#
# `block` names which of the warm fit's coupled fields this is, so its hypers and
# its outer-grid axes are read off that block; `weight` is the per-site design
# column of a varying-coefficient field (NULL = the unweighted intercept field).
#
# `copy.slab` declares the measure the sampled copy amplitude carries over its
# span: "flat" in log alpha, or "exponential" for the penalized-complexity slab
# the outer grid weighs against its point mass at alpha = 0. The sampler has no
# way to visit that point mass, so under "exponential" it targets the coupled
# component of the engine's mixture prior, not the mixture.
.occu_cover_nuts_field_block <- function(adj, type, n_cells, site_cell, warm,
                                         scale_factor = NULL,
                                         sample_hyper = TRUE,
                                         block = 1L, weight = NULL,
                                         copy.slab = NULL) {
  copy.slab <- .occu_cover_nuts_copy_slab(copy.slab)
  wh <- if (length(warm$blocks) >= block) warm$blocks[[block]] else warm
  wt <- if (is.null(weight)) NULL else as.numeric(weight)
  if (!sample_hyper) {
    fl <- .tobs_nuts_field_loading(adj, type, n_cells,
                                   tau = 1 / max(wh$sigma, 1e-3)^2,
                                   rho = wh$rho, sigma = wh$sigma,
                                   scale_factor = scale_factor)
    entries <- list(n_field_units = as.integer(n_cells),
                    field_map = as.integer(site_cell),
                    field_load = fl$field_load,
                    field_alpha = as.numeric(wh$alpha))
    if (!is.null(wt)) entries$field_weight <- wt
    return(list(entries = entries, n_raw = fl$n_raw, theta0_hyper = numeric(0),
                sampled = character(0),
                pinned = c(sigma = wh$sigma, rho = fl$rho,
                           if (.ochf_has_copy(wh$alpha)) c(alpha = wh$alpha)),
                tau = fl$tau, rho = fl$rho))
  }

  bas <- .occu_cover_nuts_field_basis(adj, type, n_cells, scale_factor)
  bnd <- .occu_cover_nuts_hyper_bounds(warm, type, block = block)
  bnd$rho <- .ochf_rho_support(
    bnd$rho, if (bas$scale1 == .OCHF_CAR) bas$lambda else numeric(0))
  entries <- list(n_field_units = as.integer(n_cells),
                  field_map = as.integer(site_cell),
                  field_load = bas$B1,
                  field_scale1 = as.integer(bas$scale1),
                  field_has_iid = as.integer(bas$has_iid),
                  field_sf = as.numeric(bas$sf),
                  field_lambda = as.numeric(bas$lambda))
  if (!is.null(wt)) entries$field_weight <- wt
  n_raw <- ncol(bas$B1) + if (bas$has_iid) n_cells else 0L

  theta0_hyper <- numeric(0); sampled <- character(0); pinned <- numeric(0)
  # Warm-start a sampled coordinate at the grid-integrated estimate, held off the
  # transform's flat tails so the first leapfrog step is informative.
  start_u <- function(v, lo, hi, link) {
    fr <- (.ochf_link(v, link) - .ochf_link(lo, link)) /
          (.ochf_link(hi, link) - .ochf_link(lo, link))
    stats::qlogis(min(max(fr, 0.05), 0.95))
  }
  add <- function(nm, value, bounds, link, dflt) {
    if (is.null(bounds) || !is.finite(value) || value <= 0) {
      entries[[paste0("field_", nm, "_fixed")]] <<-
        as.numeric(if (is.finite(value) && value > 0) value else dflt)
      pinned[[nm]] <<- entries[[paste0("field_", nm, "_fixed")]]
      return(invisible(NULL))
    }
    v <- min(max(value, bounds[1L]), bounds[2L])
    entries[[paste0("field_", nm, "_lo")]] <<- .ochf_link(bounds[1L], link)
    entries[[paste0("field_", nm, "_hi")]] <<- .ochf_link(bounds[2L], link)
    theta0_hyper <<- c(theta0_hyper, start_u(v, bounds[1L], bounds[2L], link))
    sampled <<- c(sampled, nm)
  }
  add("sigma", wh$sigma, bnd$sigma, 0L, 1)
  # icar carries no mixing parameter: rho = 1 is the intrinsic precision itself.
  if (identical(type, "icar")) { entries$field_rho_fixed <- 1; pinned[["rho"]] <- 1 }
  else add("rho", wh$rho, bnd$rho, 1L, 1)
  # `alpha` is the amplitude of a copy, so it is a hyperparameter of this model
  # only when there IS a copy. With none, `.occu_cover_apply_copy_coupling()`
  # resolves the amplitude to 0 -- the decoupled model -- and the loading must
  # still be defined for the target, but recording it as a PINNED hyper reports
  # an absent term as one the fit conditioned on. Set the entry and leave it out
  # of `sampled` / `pinned`, matching what the deterministic backend reports
  # (#293). A user-requested pin at a positive value (`share(spatial(),
  # alpha = 0.5)`) is a real conditioning choice and is still reported.
  #
  # NOT the same case as `rho = 1` under icar above: that is the intrinsic
  # precision, part of the icar model itself, not an absent term.
  if (.ochf_has_copy(wh$alpha)) add("alpha", wh$alpha, bnd$alpha, 0L, 0)
  else entries$field_alpha_fixed <- 0
  # The copy amplitude's declared slab, anchored on the same alpha node set the
  # engine's own density is declared over.
  if ("alpha" %in% sampled && identical(copy.slab, "exponential")) {
    rate <- .ochf_copy_slab_rate(wh$alpha_upper)
    if (is.null(rate))
      stop("occu_cover NUTS spatial: copy.slab = \"exponential\" needs a copy ",
           "grid with a positive node to anchor its rate.", call. = FALSE)
    entries$field_alpha_prior <- .OCHF_PRIOR_EXP
    entries$field_alpha_rate  <- rate
  }

  # The span each sampled hyper actually got, in natural units. A stated grid
  # names NODES, and the axis reaches half a node step past the outermost one at
  # each end, so the span is wider than the nodes and widest for the shortest
  # grid (#301). Reported so a fit says what region it worked over rather than
  # leaving it to be recomputed from the grid.
  support <- lapply(bnd[intersect(names(bnd), sampled)], as.numeric)
  list(entries = entries, n_raw = n_raw, theta0_hyper = theta0_hyper,
       sampled = sampled, pinned = pinned, tau = NA_real_, rho = wh$rho,
       support = support)
}

# Natural-scale hyper draws from the sampled field coordinates, and the per-draw
# field. `draws` is the full sampler matrix; `total` the log-dispersion index.
#
# `field_sd` is the fourth reported quantity: the geometric-mean marginal SD the
# block's own covariance implies at that draw's hypers,
#   Cov(z) = sigma^2 (B1 diag(s1^2) B1' + s2^2 I),
# so diag(Cov) needs only the row sums of (B1 s1)^2. It is the one field-scale
# summary whose meaning does not depend on the areal kind (icar's precision, the
# bym2 mixing weights and the proper-CAR eigen weights all normalise differently),
# which makes it the quantity a simulation truth can be stated in.
.ochf_hyper_draws <- function(entries, draws, total) {
  fvs <- .ochf_views(entries, total)
  n_draws <- nrow(draws)
  # Column suffix per block: the intercept field keeps the bare names, a lone
  # varying-coefficient field is `_trend`, several are indexed -- the same
  # spelling the grid-integrated route's hyper table uses.
  n_trend <- length(fvs) - 1L
  suffix  <- c("", if (n_trend == 1L) "_trend"
                   else if (n_trend > 1L) paste0("_trend", seq_len(n_trend)))
  hyp <- list()
  fields <- vector("list", length(fvs))
  for (b in seq_along(fvs)) {
    fv <- fvs[[b]]
    for (nm in c("sigma", "rho", "alpha")) {
      h <- fv[[nm]]
      # A block with no copy has no amplitude: emit no `alpha` column rather
      # than a constant zero one, so the draws table agrees with
      # `sampled_hyper` / `fixed_hyper` about what this model has (#293). The
      # loading reader defaults a missing column to 0, which is the same value
      # it would have read.
      if (identical(nm, "alpha") && is.null(h$coord) &&
          !.ochf_has_copy(h$fixed)) next
      hyp[[paste0(nm, suffix[b])]] <-
        if (is.null(h$coord)) rep(h$fixed, n_draws)
        else .ochf_inv_link(h$t_lo + (h$t_hi - h$t_lo) *
                              stats::plogis(draws[, h$coord]), h$link)
    }
    z   <- matrix(NA_real_, n_draws, fv$n_units)
    fsd <- numeric(n_draws)
    B2  <- fv$B1^2
    for (i in seq_len(n_draws)) {
      fs <- .ochf_forward(fv, draws[i, ])
      z[i, ] <- fs$z
      v <- as.numeric(B2 %*% (fs$s1^2)) + if (fv$has_iid) fs$s2^2 else 0
      fsd[i] <- fs$sigma$value * exp(mean(log(pmax(v, 1e-300))) / 2)
    }
    hyp[[paste0("field_sd", suffix[b])]] <- fsd
    fields[[b]] <- z
  }
  list(hyper = do.call(cbind, hyp), field = fields[[1L]], fields = fields,
       suffix = suffix)
}


# Sample the exact occu_cover coefficient posterior jointly with a coupled
# non-centered areal field on the latent state z: the psi-arm field z (one value
# per cell) enters psi linearly and is copied to the cover (positive) arm with
# the amplitude alpha. Parameter vector: c(beta_psi, beta_p, beta_pos, log_disp,
# raw_field, u_sigma?, u_rho?, u_alpha?).
#
# The field SD sigma, the mixing / spatial-correlation rho and the copy amplitude
# alpha are SAMPLED: every areal kind's loading factors as a fixed basis with
# hyper-dependent column weights, so the joint density costs a scalar (or
# per-column) rescale per leapfrog step and no re-decomposition. Their priors are
# flat in the coordinate the nested-Laplace outer grid spaces its nodes in, over
# that grid's own span, so the sampler integrates the same hyper measure the
# deterministic backend does -- and the sampler is then an INDEPENDENT reference
# for it rather than a fit conditioned on its point estimate. icar carries no
# mixing parameter (rho = 1 is the intrinsic precision), and an axis the grid
# pinned to a single node stays pinned; `fit$nuts$sampled_hyper` / `$fixed_hyper`
# report which is which per fit. `fixed.hyper = TRUE` restores the #74 / #113
# behaviour, conditioning on the warm fit's (sigma, rho, alpha).
#
# This is the occu_cover analogue of .tobs_fit_abun_nuts_spatial, with the field
# block in src/occu_cover_nuts.cpp byte-exact vs the R oracle's field branch.
# Intrinsic icar / bym2 fields sample through the sum-to-zero eigen-basis that
# drops the precision null-space (constant) direction (/#113).
#
# `trends` adds one varying-coefficient field per per-site weight column beside
# the intercept field. Each is its own block -- own whitened field, own (sigma,
# rho, alpha) coordinates -- so the parameter vector grows to c(beta_psi, beta_p,
# beta_pos, log_disp, [raw_b, u_sigma_b?, u_rho_b?, u_alpha_b?] per block). `...`
# absorbs unused sampler controls.
.tobs_fit_occu_cover_nuts_spatial <- function(model, spatial, priors = NULL,
                                              trends = list(),
                                              sigma.beta = NULL, sigma.logdisp = 5,
                                              n.iter = NULL, n.warmup = NULL,
                                              n.chains = NULL, n.thin = NULL,
                                              n.threads = NULL,
                                              max.treedepth = NULL,
                                              adapt.delta = NULL, seed = NULL,
                                              max.iter = 200L, tol = 1e-6,
                                              verbose = FALSE,
                                              sigma.grid = NULL, rho.car.grid = NULL,
                                              alpha.grid = NULL,
                                              alpha.n = NULL,
                                              alpha.grid.trend = NULL,
                                              alpha.n.trend = NULL,
                                              fixed.hyper = FALSE,
                                              copy.slab = NULL,
                                              copy.atom.mass = NULL, ...) {
  # Sampler defaults come from the one engine table. This path carries its own
  # adaptation target there (the sampled proper-CAR rho reaches the
  # near-intrinsic boundary; see .TOBS_FAMILY_DEFAULTS).
  .tobs_fill_sampler(environment(), "nuts", "occu_cover_spatial")

  slab <- .occu_cover_nuts_copy_slab(copy.slab)
  .tobs_reject_weighted_spatial(spatial, "occu_cover NUTS psi spatial")
  if (!spatial$type %in% c("icar", "car_proper", "bym2")) {
    stop(sprintf(paste0(
      "occu_cover() NUTS + areal spatial supports icar() / car_proper() / ",
      "bym2() on the psi arm (coupled to the cover arm); got '%s'. Use ",
      "method = \"nested_laplace\" for other field kinds."),
      spatial$type), call. = FALSE)
  }
  pin   <- model$process_info
  p_occ <- pin[[1L]]$p; p_p <- pin[[2L]]$p; p_pos <- pin[[3L]]$p
  n_coef <- p_occ + p_p + p_pos
  n_base <- n_coef + 1L                    # + log_dispersion

  n_sites   <- model$n_sites
  site_cell <- model$site_cell %||% seq_len(n_sites)
  adj       <- as.matrix(spatial$graph)
  n_cells   <- nrow(adj)
  if (length(site_cell) != n_sites || max(site_cell) > n_cells ||
      min(site_cell) < 1L) {
    stop(sprintf(paste0(
      "occu_cover NUTS spatial: site_cell must map %d sites into 1..%d ",
      "graph nodes."), n_sites, n_cells), call. = FALSE)
  }

  # Warm nested-Laplace joint fit: betas, field, and the grid-integrated
  # (sigma, rho, alpha). Under fixed.hyper it supplies the pinned values; when
  # the hypers are sampled it supplies their starting values and, through its own
  # outer grid, the support of their priors.
  warm <- .tobs_occu_cover_nuts_carproper_warm(
    model, adj, priors, type = spatial$type, max.iter = max.iter, tol = tol,
    sigma.grid = sigma.grid, rho.car.grid = rho.car.grid, alpha.grid = alpha.grid,
    alpha.n = alpha.n, trends = trends, alpha.grid.trend = alpha.grid.trend,
    alpha.n.trend = alpha.n.trend, copy.slab = copy.slab,
    copy.atom.mass = copy.atom.mass)

  # One block per coupled field: the unweighted intercept field, then one per
  # varying-coefficient weight column. Each reads its own warm hypers and its own
  # outer-grid axes, so a trend field's amplitude is its own coordinate rather
  # than the intercept field's.
  n_trend  <- length(trends)
  weights  <- c(list(NULL), lapply(trends, `[[`, "weight"))
  for (b in seq_along(trends)) {
    w <- weights[[b + 1L]]
    if (length(w) != n_sites || anyNA(w))
      stop(sprintf(paste0(
        "occu_cover NUTS spatial: the varying-coefficient field '%s' needs one ",
        "finite weight per site (%d); got %d."),
        trends[[b]]$label %||% "trend", n_sites, length(w)), call. = FALSE)
  }
  fbs <- lapply(seq_along(weights), function(b)
    .occu_cover_nuts_field_block(
      adj, spatial$type, n_cells, site_cell, warm,
      scale_factor = spatial$scale_factor,
      sample_hyper = !isTRUE(fixed.hyper), block = b, weight = weights[[b]],
      copy.slab = slab))
  n_raw   <- vapply(fbs, `[[`, numeric(1), "n_raw")
  n_hyper <- vapply(fbs, function(f) length(f$sampled), numeric(1))
  # Suffix per block, matching the hyper table's own column spelling.
  suffix <- c("", if (n_trend == 1L) "_trend"
                  else if (n_trend > 1L) paste0("_trend", seq_len(n_trend)))

  # car_proper warm-starts raw at the integrated field, projected onto the
  # block's own basis (B1 %*% (s1 * raw) = f / sigma); the intrinsic icar / bym2
  # loadings start raw at 0.
  raw_start <- function(b) {
    fb <- fbs[[b]]
    wb <- warm$blocks[[b]]
    if (!identical(spatial$type, "car_proper")) return(numeric(fb$n_raw))
    if (isTRUE(fixed.hyper)) {
      Q <- .areal_Q(adj, fb$rho)
      L <- tryCatch(chol(fb$tau * Q + diag(1e-4 * fb$tau, n_cells)),
                    error = function(e) NULL)
      if (is.null(L)) stop("occu_cover NUTS spatial: field precision not PD.",
                           call. = FALSE)
      return(as.numeric(L %*% wb$field))
    }
    # B1 = D^{-1/2} U with orthonormal U, so B1^{-1} = U' D^{1/2} and the column
    # scaling divides out: raw = (1 - rho lambda)^{1/2} U' D^{1/2} f / sigma.
    deg <- rowSums(adj != 0)
    lam <- fb$entries$field_lambda
    s1  <- 1 / sqrt(pmax(1 - wb$rho * lam, 1e-8))
    as.numeric(crossprod(fb$entries$field_load * deg,
                         wb$field)) / (s1 * max(wb$sigma, 1e-3))
  }
  # The flat vector lays each block out as [raw, sampled hypers], block by block,
  # exactly as hyper_field_build_list() reads it.
  block_start <- unlist(lapply(seq_along(fbs), function(b)
    c(raw_start(b), fbs[[b]]$theta0_hyper)))
  theta0 <- c(warm$beta_psi, warm$beta_p, warm$beta_pos, warm$log_disp,
              block_start)

  spec <- .tobs_occu_cover_nuts_spec(model)
  spec$field_blocks <- lapply(fbs, `[[`, "entries")

  # Sized from the ACTUAL coordinate count -- the engine takes the metric as a
  # bare pointer, so a short vector reads past its end.
  inv_metric <- c(rep(0.1, n_base), rep(1, sum(n_raw) + sum(n_hyper)))
  if (length(inv_metric) != length(theta0))
    stop("occu_cover NUTS spatial: internal layout mismatch (metric ",
         length(inv_metric), " vs theta ", length(theta0), ").", call. = FALSE)

  par_names <- c(
    paste0("psi_", pin[[1L]]$coef_names),
    paste0("p_",   pin[[2L]]$coef_names),
    paste0("pos_", pin[[3L]]$coef_names),
    if (identical(model$positive, "beta")) "log_phi" else "log_sigma_pos")
  # paste0() recycles a zero-length argument to "", so a block with every hyper
  # pinned would otherwise gain a phantom "u_" column name.
  block_names <- unlist(lapply(seq_along(fbs), function(b) c(
    paste0("raw", suffix[b], "_", seq_len(fbs[[b]]$n_raw)),
    if (length(fbs[[b]]$sampled) > 0L)
      paste0("u_", fbs[[b]]$sampled, suffix[b]) else character(0))))
  all_names <- c(par_names, block_names)

  run_chain <- function(ch) {
    cpp_occu_cover_nuts(
      spec, theta0 = theta0, sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp,
      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
      n_warmup = as.integer(n.warmup), max_treedepth = as.integer(max.treedepth),
      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
      verbose = isTRUE(verbose) && ch == 1L)
  }
  n_chains <- max(1L, as.integer(n.chains))
  chains   <- lapply(.tobs_nuts_run_parallel(run_chain, n_chains, n.threads),
                     .tobs_nuts_thin_chain, n.thin = n.thin)
  per_chain_draws <- lapply(chains, `[[`, "draws")
  draws    <- do.call(rbind, per_chain_draws)
  colnames(draws) <- all_names
  n_draws  <- nrow(draws)

  b_idx <- seq_len(n_base)
  means <- colMeans(draws)
  V_post <- stats::cov(draws[, b_idx, drop = FALSE])
  dimnames(V_post) <- list(par_names, par_names)
  sds <- sqrt(pmax(diag(V_post), 0)); names(sds) <- par_names
  par_means <- means[b_idx]; names(par_means) <- par_names

  # Per-draw field and natural-scale hypers. The field is a nonlinear function of
  # the sampled hypers, so its posterior mean is the mean of the per-draw field,
  # not the field at the mean coordinate.
  hd <- .ochf_hyper_draws(spec$field_blocks, draws, n_base)
  field_means <- lapply(hd$fields, colMeans)
  field_mean  <- field_means[[1L]]
  hyper_draws <- hd$hyper
  # A pinned block folds sigma and rho into the loading's columns rather than
  # carrying them as factors, so the block reports them as 1. Restate the values
  # the fit actually conditioned on -- `field_sd` is unaffected, since the scaled
  # loading already carries them.
  for (b in seq_along(fbs)) {
    pin_b <- fbs[[b]]$pinned
    for (nm in names(pin_b) %||% character(0)) {
      col <- paste0(nm, suffix[b])
      if (col %in% colnames(hyper_draws)) hyper_draws[, col] <- pin_b[[nm]]
    }
  }
  hyper_mean  <- colMeans(hyper_draws)
  hyper_sd    <- apply(hyper_draws, 2L, stats::sd)
  # A positive variance component at a few dozen binary sites has a right-skewed
  # posterior, where the mean sits above the bulk; the median is the summary to
  # quote against a known truth.
  hyper_median <- apply(hyper_draws, 2L, stats::median)
  alpha       <- .ochf_block_alpha(hyper_mean, "")
  # Per-block sampled / pinned hyper names, suffixed as their columns are.
  # paste0() recycles a zero-length argument to "", so a block that samples (or
  # pins) nothing would otherwise contribute a phantom bare-suffix name.
  suffixed <- function(nm, b)
    if (length(nm)) paste0(nm, suffix[b]) else character(0)
  sampled_nm <- unlist(lapply(seq_along(fbs), function(b)
    suffixed(fbs[[b]]$sampled, b)))
  pinned_nm  <- unlist(lapply(seq_along(fbs), function(b)
    suffixed(names(fbs[[b]]$pinned), b)))
  pinned_val <- unlist(lapply(seq_along(fbs), function(b) {
    p <- fbs[[b]]$pinned
    if (!length(p)) return(numeric(0))
    stats::setNames(as.numeric(p), paste0(names(p), suffix[b]))
  }))
  # unlist() over per-block character(0) returns NULL, and `fixed_hyper` is a
  # character vector by contract (empty when nothing is pinned).
  if (is.null(sampled_nm)) sampled_nm <- character(0)
  if (is.null(pinned_nm))  pinned_nm  <- character(0)
  if (is.null(pinned_val)) pinned_val <- numeric(0)

  # Data log-likelihood at the posterior mean (every field on psi, each scaled by
  # its own alpha on cover), so logLik() matches the laplace convention.
  bo   <- par_means[seq_len(p_occ)]
  bp   <- par_means[p_occ + seq_len(p_p)]
  bpos <- par_means[p_occ + p_p + seq_len(p_pos)]
  eta  <- .occu_cover_eta_from_par(model, bo, bp, bpos)
  site_load <- lapply(seq_along(fbs), function(b) {
    w <- weights[[b]] %||% rep(1, n_sites)
    as.numeric(w) * field_means[[b]][site_cell]
  })
  f_site <- Reduce(`+`, site_load)
  # A block with no copy loads nothing onto the cover arm, and reports no
  # `alpha` to read; `.ochf_block_alpha()` supplies the 0 it would have carried.
  pos_site_load <- Reduce(`+`, lapply(seq_along(fbs), function(b)
    .ochf_block_alpha(hyper_mean, suffix[b]) * site_load[[b]]))
  psi_f  <- stats::plogis(.tobs_clamp_eta(stats::qlogis(eta$psi) + f_site))
  ep_f   <- eta$ep_mat + pos_site_load
  ll_mean <- sum(.occu_cover_site_ll(model, psi_f, eta$p_mat, ep_f, par_means[n_base]))

  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- as.integer(unlist(lapply(chains, `[[`, "divergent")))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- mean(vapply(chains, function(ch) ch$epsilon %||% NA_real_, numeric(1)),
                    na.rm = TRUE)

  # Per-hyper convergence over the natural-scale draws, chain by chain (the
  # bounded transform is monotone, so this is the split-Rhat of the coordinate
  # itself). Reported for the sampled hypers only.
  hyper_diag <- NULL
  if (length(sampled_nm) > 0L) {
    ends <- cumsum(vapply(per_chain_draws, nrow, integer(1)))
    starts <- c(1L, ends[-length(ends)] + 1L)
    hyper_chain <- lapply(seq_len(n_chains), function(ch)
      hyper_draws[starts[ch]:ends[ch], sampled_nm, drop = FALSE])
    hyper_diag <- .tobs_nuts_rhat_ess(hyper_chain)
    names(hyper_diag$rhat) <- names(hyper_diag$ess) <- sampled_nm
  }

  # The honest hyper report: `sampled_hyper` names the hypers this fit integrated
  # over, `fixed_hyper` those it conditioned on, and `fixed_hyper_values` says at
  # what. `fixed_hyper` is a character vector -- empty when nothing is pinned --
  # so a fit can never claim to have sampled a hyper it did not. With a
  # varying-coefficient field beside the intercept one each block's hypers carry
  # that block's suffix, so the two are never conflated.
  nuts <- list(accept_prob = accept, divergent = divergent, treedepth = treedepth,
               epsilon = epsilon, n_chains = n_chains,
               divergent_total = sum(divergent),
               sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp,
               field_rho = hyper_mean[["rho"]], field_alpha = alpha,
               field_sigma = hyper_mean[["sigma"]],
               sampled_hyper = sampled_nm,
               fixed_hyper = pinned_nm,
               fixed_hyper_values = pinned_val,
               # The span each sampled hyper was worked over, natural units,
               # suffixed per block like `sampled_hyper` (#301). A stated
               # `alpha = grid(c(...))` names nodes; this is the axis those
               # nodes imply, which reaches half a node step past each end.
               hyper_support = unlist(lapply(seq_along(fbs), function(b) {
                 s <- fbs[[b]]$support
                 if (!length(s)) return(NULL)
                 stats::setNames(s, paste0(names(s), suffix[b]))
               }), recursive = FALSE),
               hyper_mean = hyper_mean, hyper_median = hyper_median,
               hyper_sd = hyper_sd,
               hyper_rhat = hyper_diag$rhat, hyper_ess = hyper_diag$ess,
               warm_hyper = c(sigma = warm$sigma, rho = warm$rho,
                              alpha = warm$alpha),
               n_fields = length(fbs),
               n_field_units = n_cells)

  fit <- structure(list(
    draws        = draws[, b_idx, drop = FALSE],
    means        = par_means,
    sds          = sds,
    vcov         = V_post,
    n_samples    = n_draws,
    n_params     = n_base,
    log_prob     = rep(ll_mean, n_draws),
    log_lik      = ll_mean,
    N            = sum(model$valid),
    accept_prob  = accept,
    divergent    = divergent,
    treedepth    = treedepth,
    epsilon      = epsilon,
    col_names    = par_names,
    param_names  = par_names,
    process_info = pin,
    model        = model,
    spatial      = list(type = spatial$type, graph = adj,
                        sigma_mean = hyper_mean[["sigma"]],
                        sigma_sd = hyper_sd[["sigma"]],
                        # NULL when the field is not copied onto the cover arm:
                        # there is no amplitude, so there is none to summarise.
                        alpha_mean = if ("alpha" %in% names(hyper_mean)) alpha
                                     else NULL,
                        alpha_sd = if ("alpha" %in% names(hyper_sd))
                                     hyper_sd[["alpha"]] else NULL,
                        rho_mean = hyper_mean[["rho"]],
                        rho_sd = hyper_sd[["rho"]],
                        n_fields = length(fbs),
                        field_labels = c("intercept",
                                         vapply(trends, function(tf)
                                           as.character(tf$label %||% "trend"),
                                           character(1))),
                        # Per-block suffix and per-site design weight, so the
                        # criteria read every block's loading off the fit rather
                        # than re-deriving it.
                        field_suffix = suffix, field_weights = weights,
                        sampled_hyper = sampled_nm, fixed_hyper = pinned_nm),
    spatial_field = field_mean,
    # Varying-coefficient surfaces, named by their weight column, beside the
    # intercept field.
    trend_field   = if (n_trend > 0L) field_means[[2L]] else NULL,
    trend_fields  = if (n_trend > 0L)
      stats::setNames(field_means[-1L],
                      vapply(trends, function(tf)
                        as.character(tf$label %||% "trend"), character(1)))
      else NULL,
    field_draws   = hd$field,
    trend_field_draws = if (n_trend > 0L) hd$fields[-1L] else NULL,
    hyper_draws   = hyper_draws,
    method       = "nuts",
    positive     = model$positive,
    nuts         = nuts,
    convergence  = list(converged = NA, n_iter = as.integer(n.iter))
  ), class = c("tobs_fit", "tulpa_fit"))

  # Diagnostics over the coefficient block (the coordinates the fit reports);
  # the whitened field `raw` coordinates carry no named parameter.
  fit <- .tobs_nuts_attach_convergence(fit, per_chain_draws, par_names = par_names,
                                       cols = b_idx)
  fit$nuts$rhat <- fit$convergence$rhat
  fit$nuts$ess  <- fit$convergence$ess_bulk
  fit
}
