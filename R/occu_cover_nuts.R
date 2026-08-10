# occu_cover_nuts.R - NUTS target for the non-spatial joint occupancy +
# cover-hurdle family (occu_cover()).
#
# The Laplace fit (.tobs_fit_occu_cover) sums the latent occupancy state z out in
# closed form (two states per cell) and returns a Gaussian observed-Fisher
# posterior over the packed coefficient vector
#   theta = c(beta_psi, beta_p, beta_pos, log_dispersion).
# NUTS instead samples the exact marginal posterior of that vector, giving
# calibrated (non-Gaussian) intervals and the per-draw pointwise likelihood that
# WAIC / LOO need. There is no latent field, random effect, or community
# covariance on the non-spatial path, so the parameter vector is just the flat
# three-arm coefficient block plus one log-dispersion scalar.
#
# .tobs_occu_cover_nuts_logpost is the R oracle: it recomputes the joint
# log-posterior and gradient exactly as the C++ FullGradFn (src/occu_cover_nuts.cpp,
# cpp_occu_cover_nuts) does, and a byte-exact test cross-checks the two before the
# sampler is trusted. The per-cell marginal mirrors .occu_cover_site_ll; the cover
# arm and the no-detection mixture reuse the same closed forms the Laplace path uses.


# ---------------------------------------------------------------------------
# Coupled areal field with sampled hyperparameters (gcol33/tulpaObs#204)
#
# R mirror of src/nuts_field_hyper.h. The field is
#   z = sigma * ( B1 %*% (s1(rho) * raw1) + s2(rho) * raw2 )
# over a FIXED basis B1, with the field SD sigma, the mixing / spatial-
# correlation rho and the cross-arm copy amplitude alpha either sampled as
# bounded coordinates or pinned. A sampled hyper rides an unconstrained u,
#   t = t_lo + (t_hi - t_lo) * plogis(u),   value = inv_link(t),
# with `t` the coordinate the nested-Laplace outer grid spaces its nodes in
# (log for sigma / alpha, logit for rho) and [t_lo, t_hi] that grid's own span.
# The prior is flat in t, which is the measure the grid integrates against
# (equal weight per cell); flat prior + change of variables leaves the
# normalised log-density log(e) + log(1 - e), e = plogis(u).
# ---------------------------------------------------------------------------

.OCHF_CONST <- 0L; .OCHF_BYM2_STR <- 1L; .OCHF_BYM2_IID <- 2L; .OCHF_CAR <- 3L

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
  slot <- function(fixed_key, lo_key, hi_key, link, dflt) {
    h <- list(link = link, fixed = field[[fixed_key]] %||% dflt, coord = NULL)
    lo <- field[[lo_key]]; hi <- field[[hi_key]]
    if (!is.null(lo) && !is.null(hi)) {
      k <<- k + 1L
      h$t_lo <- lo; h$t_hi <- hi; h$coord <- k
    }
    h
  }
  sigma <- slot("field_sigma_fixed", "field_sigma_lo", "field_sigma_hi", 0L, 1)
  rho   <- slot("field_rho_fixed",   "field_rho_lo",   "field_rho_hi",   1L, 1)
  alpha <- slot("field_alpha_fixed", "field_alpha_lo", "field_alpha_hi", 0L, 0)
  if (is.null(alpha$coord) && is.null(field$field_alpha_fixed))
    alpha$fixed <- field$field_alpha %||% field$alpha %||% 0
  list(n_units = n_units, m1 = m1, B1 = B1,
       scale1 = as.integer(field$field_scale1 %||% .OCHF_CONST),
       has_iid = has_iid, sf = as.numeric(field$field_sf %||% 1),
       lambda = as.numeric(field$field_lambda %||% numeric(0)),
       n_raw = n_raw, o_raw = total,
       field_map = as.integer(field$field_map),
       sigma = sigma, rho = rho, alpha = alpha)
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
  add <- function(h, hv, dlp_dvalue) {
    if (is.null(h$coord)) return(invisible(NULL))
    g_hyper <<- c(g_hyper, dlp_dvalue * hv$dvalue_dt * hv$dt_du + (1 - 2 * hv$e))
    lp <<- lp + log(hv$e) + log(1 - hv$e)
  }
  add(fv$sigma, fs$sigma, if (sigma > 0) sum(g_z * fs$z) / sigma else 0)
  add(fv$rho,   fs$rho,   g_rho)
  add(fv$alpha, fs$alpha, g_alpha_data)
  list(grad = c(g_raw, g_hyper), lp = lp)
}


# Joint log-posterior + gradient of the non-spatial occu_cover coefficient vector
# theta = c(beta_psi, beta_p, beta_pos, log_disp) (beta_p / beta_pos each packed as
# the site-level block then the optional visit-level block, exactly as the fitter
# stacks them). Weak Gaussian priors N(0, sigma.beta^2) on every coefficient and a
# broad N(0, sigma.logdisp^2) on log_disp keep the ridge / dispersion proper
# without materially shifting the data-dominated optimum. Returns list(lp, grad)
# over the packed coordinates. This is the oracle the C++ FullGradFn mirrors.
.tobs_occu_cover_nuts_logpost <- function(theta, model, sigma.beta = 5,
                                          sigma.logdisp = 5, field = NULL) {
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

  # Optional coupled field (gcol33/tulpaObs#74, #204): z over the trailing field
  # block of theta, entering psi additively and the cover arm scaled by the copy
  # amplitude alpha. sigma / rho / alpha are sampled coordinates or pinned, as
  # `field` declares. The R oracle mirrors the C++ FullGradFn byte-for-byte;
  # absent, the non-spatial target.
  has_field <- !is.null(field)
  if (has_field) {
    fv     <- .ochf_view(field, total)
    fs     <- .ochf_forward(fv, theta)
    f_site <- fs$z[fv$field_map]                  # per-site field value
    alpha  <- fs$alpha$value
  } else {
    f_site <- numeric(model$n_sites)
    alpha  <- 0
  }

  bo         <- theta[seq_len(p_occ)]
  bp_site    <- theta[p_occ + seq_len(p_det_site)]
  bp_visit   <- if (p_det_visit > 0L) theta[p_occ + p_det_site + seq_len(p_det_visit)] else numeric(0)
  bpos_site  <- theta[p_occ + p_p + seq_len(p_pos_site)]
  bpos_visit <- if (p_pos_visit > 0L) theta[p_occ + p_p + p_pos_site + seq_len(p_pos_visit)] else numeric(0)
  log_disp   <- theta[total]
  disp       <- exp(log_disp)
  pos_code   <- .tobs_cover_pos_code(model$positive)

  N <- model$n_sites; J <- model$max_visits
  sgm <- function(e) 1 / (1 + exp(-e))

  eta_psi <- as.numeric(model$X_occ %*% bo) + f_site
  psi     <- sgm(eta_psi)

  eta_p <- matrix(as.numeric(model$X_det_site %*% bp_site), N, J)   # site block broadcast
  if (p_det_visit > 0L) {
    eta_p <- eta_p + matrix(as.numeric(model$X_det_visit %*% bp_visit), N, J, byrow = TRUE)
  }
  eta_pos <- matrix(as.numeric(model$X_pos_site %*% bpos_site) + alpha * f_site, N, J)
  if (p_pos_visit > 0L) {
    eta_pos <- eta_pos + matrix(as.numeric(model$X_pos_visit %*% bpos_visit), N, J, byrow = TRUE)
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
          ly <- log(yy); l1my <- log(1 - yy)
          lp_data <- lp_data + lgamma(disp) - lgamma(a) - lgamma(b) +
                     (a - 1) * ly + (b - 1) * l1my
          g_eta_pos[i, v] <- disp * mu * (1 - mu) *
                             (-digamma(a) + digamma(b) + ly - l1my)
          g_logdisp <- g_logdisp + disp * (digamma(disp) - mu * digamma(a) -
                       (1 - mu) * digamma(b) + mu * ly + (1 - mu) * l1my)
        } else if (pos_code == 4L) {     # identity-Gaussian (#112): raw response
          sig <- disp; r <- (yy - ev) / sig
          lp_data <- lp_data - log(sig) - 0.5 * log(2 * pi) - 0.5 * r * r
          g_eta_pos[i, v] <- r / sig
          g_logdisp <- g_logdisp + (r * r - 1)
        } else {                         # lognormal: Gaussian on log(cover)
          sig <- disp; r <- (log(yy) - ev) / sig
          lp_data <- lp_data - log(yy) - log(sig) - 0.5 * log(2 * pi) - 0.5 * r * r
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

  # Field block: accumulate the per-cell field score (psi arm + alpha * cover
  # arm) and the copy amplitude's own data-score, then hand both to the field
  # backward pass. Its gradient (whitened field, then any sampled hyper) appends
  # to the coefficient block.
  lp_field <- 0
  if (has_field) {
    g_pos_site <- rowSums(g_eta_pos)
    g_f_site   <- g_eta_psi + alpha * g_pos_site
    g_field    <- numeric(fv$n_units)
    for (s in seq_len(N)) {
      cell <- fv$field_map[s]
      g_field[cell] <- g_field[cell] + g_f_site[s]
    }
    bk       <- .ochf_backward(fv, fs, g_field, sum(f_site * g_pos_site))
    grad     <- c(grad, bk$grad)
    lp_field <- bk$lp
  }

  lp  <- lp_data + lp_field
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
    pos_code    = .tobs_cover_pos_code(model$positive),
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
# scale. `sigma.logdisp` is an internal weak-prior width (no control knob, like
# abun's sigma.logr). `...` absorbs unused sampler controls (n.thin / n.threads /
# progress.*).
.tobs_fit_occu_cover_nuts <- function(model, priors = NULL,
                                      sigma.beta = NULL, sigma.logdisp = 5,
                                      n.iter = NULL, n.warmup = NULL,
                                      n.chains = NULL, max.treedepth = NULL,
                                      adapt.delta = NULL, seed = NULL,
                                      verbose = FALSE, ...) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188).
  .tobs_fill_sampler(environment(), "nuts")

  pin   <- model$process_info
  p_occ <- pin[[1L]]$p; p_p <- pin[[2L]]$p; p_pos <- pin[[3L]]$p
  n_par <- p_occ + p_p + p_pos + 1L

  par_names <- c(
    paste0("psi_", pin[[1L]]$coef_names),
    paste0("p_",   pin[[2L]]$coef_names),
    paste0("pos_", pin[[3L]]$coef_names),
    if (identical(model$positive, "beta")) "log_phi" else "log_sigma_pos"
  )

  # Warm start at the Laplace mode + a diagonal Laplace metric from its vcov.
  warm <- .tobs_fit_occu_cover(model, method = "laplace", priors = priors,
                               sigma.beta = sigma.beta, verbose = FALSE)
  theta0 <- as.numeric(warm$means)
  V <- as.matrix(warm$vcov)
  inv_metric <- if (!is.null(V) && all(dim(V) == n_par) && all(is.finite(diag(V))))
                  pmax(diag(V), 1e-6) else rep(1, n_par)

  spec <- .tobs_occu_cover_nuts_spec(model)

  run_chain <- function(ch) {
    cpp_occu_cover_nuts(
      spec, theta0 = theta0, sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp,
      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
      n_warmup = as.integer(n.warmup), max_treedepth = as.integer(max.treedepth),
      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
      verbose = isTRUE(verbose) && ch == 1L)
  }
  n_chains <- max(1L, as.integer(n.chains))
  chains   <- lapply(seq_len(n_chains), run_chain)
  per_chain_draws <- lapply(chains, `[[`, "draws")
  draws    <- do.call(rbind, per_chain_draws)
  colnames(draws) <- par_names
  n_draws  <- nrow(draws)

  means  <- colMeans(draws); names(means) <- par_names
  V_post <- stats::cov(draws); dimnames(V_post) <- list(par_names, par_names)
  sds    <- sqrt(pmax(diag(V_post), 0)); names(sds) <- par_names

  # Data log-likelihood at the posterior mean (scale-invariant), so logLik() on
  # the NUTS fit matches the laplace-path convention. Reuses the shared marginal.
  bo   <- means[seq_len(p_occ)]
  bp   <- means[p_occ + seq_len(p_p)]
  bpos <- means[p_occ + p_p + seq_len(p_pos)]
  eta  <- .occu_cover_eta_from_par(model, bo, bp, bpos)
  ll_mean <- sum(.occu_cover_site_ll(model, eta$psi, eta$p_mat, eta$ep_mat, means[n_par]))

  accept    <- unlist(lapply(chains, `[[`, "accept_prob"))
  divergent <- as.integer(unlist(lapply(chains, `[[`, "divergent")))
  treedepth <- as.integer(unlist(lapply(chains, `[[`, "treedepth")))
  epsilon   <- mean(vapply(chains, function(ch) ch$epsilon %||% NA_real_, numeric(1)),
                    na.rm = TRUE)

  nuts <- list(accept_prob = accept, divergent = divergent, treedepth = treedepth,
               epsilon = epsilon, n_chains = n_chains, divergent_total = sum(divergent),
               sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp)

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
    method       = "nuts",
    positive     = model$positive,
    nuts         = nuts,
    convergence  = list(converged = NA, n_iter = as.integer(n.iter))
  ), class = c("tobs_fit", "tulpa_fit"))

  # Per-parameter split-R-hat / bulk + tail ESS, through the writer every sampled
  # path shares, so summary.tobs_fit surfaces them per parameter; `fit$nuts`
  # carries the same two vectors alongside the sampler diagnostics.
  fit <- .tobs_nuts_attach_convergence(fit, per_chain_draws,
                                       par_names = par_names)
  fit$nuts$rhat <- fit$convergence$rhat
  fit$nuts$ess  <- fit$convergence$ess_bulk
  fit
}


# ---------------------------------------------------------------------------
# Spatial occu_cover NUTS: fixed-hyper non-centered coupled proper-CAR field
# (gcol33/tulpaObs#74)
# ---------------------------------------------------------------------------

# Desugar a varying-coefficient bar (`spatial(~ 1 || node, graph = adj)`) on the
# psi formula into the plain areal field spec the NUTS sampler takes, through the
# SAME expansion the nested-Laplace path uses (.tobs_expand_spatial_bar): one
# unweighted intercept field plus one weight-scaled field per bar covariate
# column, each identical to what `icar(graph = adj, group_var = node)` /
# `icar(graph = adj, weight = col, group_var = node)` builds. A single-column bar
# therefore IS the plain areal term, and routes unchanged (gcol33/tulpaObs#203).
#
# The sampler carries ONE field block -- one loading, one site -> node map and one
# copy amplitude (src/occu_cover_nuts.cpp, the shared FieldBlock layout in
# src/nuts_field_block.h) -- so a bar declaring a second field has nowhere to put
# it; that is the limitation named in the error, not the bar spelling.
.occu_cover_nuts_bar_field <- function(spec, data) {
  if (!is.null(spec$by_var)) {
    stop(paste0(
      "occu_cover() NUTS + areal spatial samples one field over one graph; a ",
      "replicated field (spatial(<bar>, by = \"", spec$by_var, "\")) needs ",
      "method = \"nested_laplace\"."), call. = FALSE)
  }
  fields <- .tobs_expand_spatial_bar(spec, data)
  if (length(fields) > 1L) {
    stop(sprintf(paste0(
      "occu_cover() NUTS + areal spatial samples a SINGLE areal field: the ",
      "sampler carries one field block (one loading, one site -> node map, one ",
      "copy amplitude onto the cover arm). This bar declares %d fields (the ",
      "intercept field plus %d varying-coefficient field(s)); the additional ",
      "field(s) need the grid-integrated method = \"nested_laplace\" path."),
      length(fields), length(fields) - 1L), call. = FALSE)
  }
  if (isTRUE(spec$correlated)) {
    stop(paste0(
      "occu_cover(): a correlated spatial bar (`|`) needs at least one ",
      "coefficient beyond the intercept (e.g. spatial(~ 1 + x | cell, ",
      "graph = adj)); a single field has no cross-covariance to estimate. ",
      "Use the independent spelling `||` for a single field."), call. = FALSE)
  }
  fields[[1L]]
}

# Resolve the single spatial term + fixed-effect psi formula for the NUTS
# spatial path. Unlike .occu_cover_spatial_fields (which gates to icar/bym2 for
# the grid-integrated nested-Laplace engine), the NUTS path also accepts
# car_proper() (the full-rank precision the fixed-hyper non-centered field is
# best conditioned on) and rejects SVC/trend/temporal/RE terms with a pointer to
# the nested-Laplace route. The bar
# form is desugared first, so `spatial(~ 1 || cell, graph = adj)` and
# `icar(graph = adj, group_var = "cell")` reach the sampler as one field
# description (gcol33/tulpaObs#203). Returns NULL when the psi formula carries no
# spatial term (the non-spatial NUTS sampler), or list(fe, spatial, group_var).
.occu_cover_nuts_spatial_term <- function(formula, data) {
  bind <- .tobs_bind_formulas(list(psi = formula), data)
  if (length(bind$terms) == 0L) return(NULL)
  spatial <- Filter(function(t) inherits(t$spec, "tobs_spatial"), bind$terms)
  re_terms <- Filter(function(t) inherits(t$spec, "tobs_re"), bind$terms)
  other   <- Filter(function(t) !inherits(t$spec, "tobs_spatial") &&
                                !inherits(t$spec, "tobs_re"), bind$terms)
  if (length(spatial) == 0L) return(NULL)
  if (length(spatial) > 1L || length(re_terms) > 0L || length(other) > 0L) {
    stop(paste0(
      "occu_cover() NUTS + areal spatial samples a SINGLE shared coupled field ",
      "(one areal term on the psi formula); a second field / per-group RE / ",
      "temporal term composes only on the grid-integrated ",
      "method = \"nested_laplace\" path."), call. = FALSE)
  }
  spec <- spatial[[1L]]$spec
  if (isTRUE(spec$is_bar)) spec <- .occu_cover_nuts_bar_field(spec, data)
  if (!is.null(spec$weight) || isTRUE(spec$is_multifield)) {
    stop(paste0(
      "occu_cover() NUTS + areal spatial samples a single unweighted shared ",
      "field; a weighted SVC field needs method = \"nested_laplace\"."),
      call. = FALSE)
  }
  list(fe = bind$fe$psi, spatial = spec, group_var = spec$group_var)
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
.tobs_occu_cover_nuts_carproper_warm <- function(model, adj, priors,
                                                 max.iter = 200L, tol = 1e-6,
                                                 type = "car_proper",
                                                 sigma.grid = NULL,
                                                 rho.car.grid = NULL,
                                                 alpha.grid = NULL) {
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

  alpha_grid <- alpha.grid %||% .tobs_default_alpha_grid()
  sigma_grid <- sigma.grid %||% .tobs_default_sigma_grid()
  rho_car_grid <- rho.car.grid %||% .tobs_default_rho_car_grid()

  # Single-block (multi = FALSE): the pos arm carries
  # field_coef = list(name = "alpha", grid = alpha_grid), so the copy alpha axis
  # rides the one shared field. The single-block joint path takes the copy
  # coefficient on the arm, not a top-level `copy` block.
  arms_out <- .occu_cover_build_joint_arms(
    model = model, sigma_pos_init = sigma_pos_init, alpha_grid = alpha_grid,
    positive = model$positive, multi = FALSE, n_cells = n_cells,
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
  prior_arg <- list(
    type            = type,
    n_spatial_units = csr$n_spatial_units,
    adj_row_ptr     = csr$adj_row_ptr,
    adj_col_idx     = csr$adj_col_idx,
    n_neighbors     = csr$n_neighbors,
    sigma_grid      = sigma_grid,
    spatial_idx     = lapply(responses, function(a) as.integer(a$spatial_idx)))
  # Only the proper-CAR field carries a spatial-correlation rho grid; the
  # intrinsic icar / bym2 fields fix rho (icar rho = 1, bym2 mixing gridded by the
  # engine's own bym2 axis).
  if (identical(type, "car_proper")) prior_arg$rho_car_grid <- rho_car_grid

  fit <- tulpa::tulpa_nested_laplace_joint(
    responses = responses, prior = prior_arg,
    cell_coupling = spec_name,
    control = list(max_iter = as.integer(max.iter), tol = as.numeric(tol),
                   n_threads = 1L, store_Q = FALSE, adaptive_grid = FALSE,
                   var_of_means_consistency = FALSE, diagnose_k = FALSE,
                   progress = FALSE))

  ok <- which(is.finite(fit$log_marginal))
  if (length(ok) == 0L)
    stop("occu_cover NUTS spatial: warm car_proper fit failed at every grid ",
         "cell. Bump control$max.iter or tighten control$tol.", call. = FALSE)
  w_raw <- exp(fit$log_marginal[ok] - max(fit$log_marginal[ok]))
  w     <- w_raw / sum(w_raw)

  layout <- fit$arm_layout
  p_psi  <- layout$p[1L]; p_p <- layout$p[2L]; p_pos <- layout$p[3L]
  bpsi_idx <- layout$beta_start[1L] + seq_len(p_psi)
  bp_idx   <- layout$beta_start[2L] + seq_len(p_p)
  bpos_idx <- layout$beta_start[3L] + seq_len(p_pos)
  f0       <- (layout$field_starts %||% layout$phi_start)
  field_idx <- f0[[1L]] + seq_len(n_cells)
  modes <- fit$modes[ok, , drop = FALSE]

  beta_psi <- as.numeric(crossprod(w, modes[, bpsi_idx, drop = FALSE]))
  beta_p   <- as.numeric(crossprod(w, modes[, bp_idx,   drop = FALSE]))
  beta_pos <- as.numeric(crossprod(w, modes[, bpos_idx, drop = FALSE]))
  field    <- as.numeric(crossprod(w, modes[, field_idx, drop = FALSE]))
  field    <- field - mean(field)            # sum-to-zero convention

  tg     <- fit$theta_grid[ok, , drop = FALSE]
  pick   <- function(nm) {
    j <- match(nm, colnames(fit$theta_grid)); if (is.na(j)) return(NA_real_)
    sum(w * as.numeric(tg[, j]))
  }
  sigma <- pick("sigma"); alpha <- pick("alpha")
  rho   <- if (identical(type, "car_proper")) pick("rho_car")
           else if (identical(type, "bym2")) pick("rho") else 1.0

  list(beta_psi = beta_psi, beta_p = beta_p, beta_pos = beta_pos,
       log_disp = log(sigma_pos_init), field = field, type = type,
       sigma = sigma, rho = rho, alpha = alpha, joint_fit = fit)
}


# Fixed basis B1 of a coupled areal field, together with the rho-scaling its
# columns carry (gcol33/tulpaObs#204). Every areal kind factors into a basis that
# does NOT depend on the sampled hypers:
#   icar / bym2 structured  the sum-to-zero eigen-loading of the intrinsic
#                           precision Q at unit precision (gcol33/tulpaObs#71),
#                           so sigma is a scalar multiply and bym2's rho a scalar
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
                sf = scale_factor %||% compute_bym2_scale(adj),
                lambda = numeric(0)))
  list(B1 = Lstr, scale1 = .OCHF_CONST, has_iid = FALSE, sf = 1,
       lambda = numeric(0))
}

# Bounds for each sampled hyper, read off the warm nested-Laplace fit's OWN
# outer grid. That grid IS the measure the deterministic backend integrates
# against -- equal weight per cell over log-spaced sigma / alpha nodes and
# logit-spaced rho nodes -- so taking its span as the support of a flat prior in
# the same coordinate makes the two backends integrate the same hyper prior, and
# makes the same `control$sigma.grid` / `alpha.grid` / `rho.car.grid` knob move
# both. An axis the grid pinned to one node (or that the fit does not carry)
# returns NULL and the hyper stays pinned.
.occu_cover_nuts_hyper_bounds <- function(warm, type) {
  tg <- warm$joint_fit$theta_grid
  ax <- function(nm, positive_only = FALSE) {
    if (is.null(tg) || is.null(colnames(tg))) return(NULL)
    j <- match(nm, colnames(tg))
    if (is.na(j)) return(NULL)
    v <- as.numeric(tg[, j]); v <- v[is.finite(v)]
    if (positive_only) v <- v[v > 0]
    if (length(v) < 2L) return(NULL)
    r <- range(v)
    if (r[2L] <= r[1L] * (1 + 1e-8)) return(NULL)
    r
  }
  list(sigma = ax("sigma"),
       alpha = ax("alpha", positive_only = TRUE),
       rho   = switch(type, bym2 = ax("rho"), car_proper = ax("rho_car"), NULL))
}

# Assemble the field spec entries the C++ block (and the R oracle) read, plus the
# warm-start values of the sampled hyper coordinates. `sample_hyper = FALSE`
# reproduces the fixed-hyper block: the loading built at the warm estimate, with
# sigma and rho baked into its columns and alpha a constant.
.occu_cover_nuts_field_block <- function(adj, type, n_cells, site_cell, warm,
                                         scale_factor = NULL,
                                         sample_hyper = TRUE) {
  if (!sample_hyper) {
    fl <- .tobs_nuts_field_loading(adj, type, n_cells,
                                   tau = 1 / max(warm$sigma, 1e-3)^2,
                                   rho = warm$rho, sigma = warm$sigma,
                                   scale_factor = scale_factor)
    entries <- list(n_field_units = as.integer(n_cells),
                    field_map = as.integer(site_cell),
                    field_load = fl$field_load,
                    field_alpha = as.numeric(warm$alpha))
    return(list(entries = entries, n_raw = fl$n_raw, theta0_hyper = numeric(0),
                sampled = character(0),
                pinned = c(sigma = warm$sigma, rho = fl$rho, alpha = warm$alpha),
                tau = fl$tau, rho = fl$rho))
  }

  bas <- .occu_cover_nuts_field_basis(adj, type, n_cells, scale_factor)
  bnd <- .occu_cover_nuts_hyper_bounds(warm, type)
  # rho rides a logit coordinate, so a grid node at 0 or 1 would put a bound at
  # infinity; both ends of the mixing / correlation range are degenerate anyway.
  if (!is.null(bnd$rho)) bnd$rho <- pmin(pmax(bnd$rho, 1e-4), 1 - 1e-4)
  # A proper-CAR precision stays positive definite only while rho lambda_j < 1.
  if (bas$scale1 == .OCHF_CAR && !is.null(bnd$rho)) {
    lam_max <- max(bas$lambda, 0)
    if (lam_max > 0)
      bnd$rho[2L] <- min(bnd$rho[2L], (1 - 1e-6) / lam_max)
    if (bnd$rho[2L] <= bnd$rho[1L]) bnd$rho <- NULL
  }
  entries <- list(n_field_units = as.integer(n_cells),
                  field_map = as.integer(site_cell),
                  field_load = bas$B1,
                  field_scale1 = as.integer(bas$scale1),
                  field_has_iid = as.integer(bas$has_iid),
                  field_sf = as.numeric(bas$sf),
                  field_lambda = as.numeric(bas$lambda))
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
  add("sigma", warm$sigma, bnd$sigma, 0L, 1)
  # icar carries no mixing parameter: rho = 1 is the intrinsic precision itself.
  if (identical(type, "icar")) { entries$field_rho_fixed <- 1; pinned[["rho"]] <- 1 }
  else add("rho", warm$rho, bnd$rho, 1L, 1)
  add("alpha", warm$alpha, bnd$alpha, 0L, 0)

  list(entries = entries, n_raw = n_raw, theta0_hyper = theta0_hyper,
       sampled = sampled, pinned = pinned, tau = NA_real_, rho = warm$rho)
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
  fv <- .ochf_view(entries, total)
  hyp <- list()
  for (nm in c("sigma", "rho", "alpha")) {
    h <- fv[[nm]]
    hyp[[nm]] <- if (is.null(h$coord)) rep(h$fixed, nrow(draws))
                 else .ochf_inv_link(h$t_lo + (h$t_hi - h$t_lo) *
                                       stats::plogis(draws[, h$coord]), h$link)
  }
  n_draws <- nrow(draws)
  z <- matrix(NA_real_, n_draws, fv$n_units)
  fsd <- numeric(n_draws)
  B2 <- fv$B1^2
  for (i in seq_len(n_draws)) {
    fs <- .ochf_forward(fv, draws[i, ])
    z[i, ] <- fs$z
    v <- as.numeric(B2 %*% (fs$s1^2)) + if (fv$has_iid) fs$s2^2 else 0
    fsd[i] <- fs$sigma$value * exp(mean(log(pmax(v, 1e-300))) / 2)
  }
  hyp$field_sd <- fsd
  list(hyper = do.call(cbind, hyp), field = z)
}


# Sample the exact occu_cover coefficient posterior jointly with a coupled
# non-centered areal field on the latent state z (gcol33/tulpaObs#74, #204):
# the psi-arm field z (one value per cell) enters psi linearly and is copied to
# the cover (positive) arm with the amplitude alpha. Parameter vector:
# c(beta_psi, beta_p, beta_pos, log_disp, raw_field, u_sigma?, u_rho?, u_alpha?).
#
# The field SD sigma, the mixing / spatial-correlation rho and the copy amplitude
# alpha are SAMPLED (gcol33/tulpaObs#204): every areal kind's loading factors as a
# fixed basis with hyper-dependent column weights, so the joint density costs a
# scalar (or per-column) rescale per leapfrog step and no re-decomposition. Their
# priors are flat in the coordinate the nested-Laplace outer grid spaces its nodes
# in, over that grid's own span, so the sampler integrates the same hyper measure
# the deterministic backend does -- and the sampler is then an INDEPENDENT
# reference for it rather than a fit conditioned on its point estimate. icar
# carries no mixing parameter (rho = 1 is the intrinsic precision), and an axis
# the grid pinned to a single node stays pinned; `fit$nuts$sampled_hyper` /
# `$fixed_hyper` report which is which per fit. `fixed.hyper = TRUE` restores the
# #74 / #113 behaviour, conditioning on the warm fit's (sigma, rho, alpha).
#
# This is the occu_cover analogue of .tobs_fit_abun_nuts_spatial (tulpa#87), with
# the field block in src/occu_cover_nuts.cpp byte-exact vs the R oracle's field
# branch. Intrinsic icar / bym2 fields sample through the sum-to-zero eigen-basis
# that drops the precision null-space (constant) direction
# (gcol33/tulpaObs#71/#113). `...` absorbs unused sampler controls.
.tobs_fit_occu_cover_nuts_spatial <- function(model, spatial, priors = NULL,
                                              sigma.beta = NULL, sigma.logdisp = 5,
                                              n.iter = NULL, n.warmup = NULL,
                                              n.chains = NULL, max.treedepth = NULL,
                                              adapt.delta = NULL, seed = NULL,
                                              max.iter = 200L, tol = 1e-6,
                                              verbose = FALSE,
                                              sigma.grid = NULL, rho.car.grid = NULL,
                                              alpha.grid = NULL,
                                              fixed.hyper = FALSE, ...) {
  # Sampler defaults come from the one engine table (gcol33/tulpaObs#188). This
  # path carries its own adaptation target there (the sampled proper-CAR rho
  # reaches the near-intrinsic boundary; see .TOBS_FAMILY_DEFAULTS).
  .tobs_fill_sampler(environment(), "nuts", "occu_cover_spatial")

  .tobs_reject_weighted_spatial(spatial, "occu_cover NUTS psi spatial")
  if (!spatial$type %in% c("icar", "car_proper", "bym2")) {
    stop(sprintf(paste0(
      "occu_cover() NUTS + areal spatial supports icar() / car_proper() / ",
      "bym2() on the psi arm (coupled to the cover arm); got '%s'. Use ",
      "method = \"nested_laplace\" for other field kinds. (gcol33/tulpaObs#74, #113)"),
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
    sigma.grid = sigma.grid, rho.car.grid = rho.car.grid, alpha.grid = alpha.grid)

  fb <- .occu_cover_nuts_field_block(
    adj, spatial$type, n_cells, site_cell, warm,
    scale_factor = spatial$scale_factor, sample_hyper = !isTRUE(fixed.hyper))
  n_raw   <- fb$n_raw
  n_hyper <- length(fb$sampled)
  alpha   <- warm$alpha

  # car_proper warm-starts raw at the integrated field, projected onto the
  # block's own basis (B1 %*% (s1 * raw) = f / sigma); the intrinsic icar / bym2
  # loadings start raw at 0.
  raw0 <- if (identical(spatial$type, "car_proper")) {
    if (isTRUE(fixed.hyper)) {
      Q <- .areal_Q(adj, fb$rho)
      L <- tryCatch(chol(fb$tau * Q + diag(1e-4 * fb$tau, n_cells)),
                    error = function(e) NULL)
      if (is.null(L)) stop("occu_cover NUTS spatial: field precision not PD.",
                           call. = FALSE)
      as.numeric(L %*% warm$field)
    } else {
      # B1 = D^{-1/2} U with orthonormal U, so B1^{-1} = U' D^{1/2} and the
      # column scaling divides out: raw = (1 - rho lambda)^{1/2} U' D^{1/2} f / sigma.
      deg  <- rowSums(adj != 0)
      lam  <- fb$entries$field_lambda
      s1   <- 1 / sqrt(pmax(1 - warm$rho * lam, 1e-8))
      as.numeric(crossprod(fb$entries$field_load * deg,
                           warm$field)) / (s1 * max(warm$sigma, 1e-3))
    }
  } else numeric(n_raw)
  theta0 <- c(warm$beta_psi, warm$beta_p, warm$beta_pos, warm$log_disp,
              raw0, fb$theta0_hyper)

  spec <- .tobs_occu_cover_nuts_spec(model)
  spec[names(fb$entries)] <- fb$entries

  inv_metric <- c(rep(0.1, n_base), rep(1, n_raw), rep(1, n_hyper))

  par_names <- c(
    paste0("psi_", pin[[1L]]$coef_names),
    paste0("p_",   pin[[2L]]$coef_names),
    paste0("pos_", pin[[3L]]$coef_names),
    if (identical(model$positive, "beta")) "log_phi" else "log_sigma_pos")
  raw_names <- paste0("raw_", seq_len(n_raw))
  # paste0() recycles a zero-length argument to "", so a fit with every hyper
  # pinned would otherwise gain a phantom "u_" column name.
  hyper_names <- if (n_hyper > 0L) paste0("u_", fb$sampled) else character(0)
  all_names <- c(par_names, raw_names, hyper_names)

  run_chain <- function(ch) {
    cpp_occu_cover_nuts(
      spec, theta0 = theta0, sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp,
      inv_metric = inv_metric, n_iter = as.integer(n.iter + n.warmup),
      n_warmup = as.integer(n.warmup), max_treedepth = as.integer(max.treedepth),
      adapt_delta = adapt.delta, seed = as.integer(seed + ch - 1L),
      verbose = isTRUE(verbose) && ch == 1L)
  }
  n_chains <- max(1L, as.integer(n.chains))
  chains   <- lapply(seq_len(n_chains), run_chain)
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
  hd <- .ochf_hyper_draws(fb$entries, draws, n_base)
  field_mean  <- colMeans(hd$field)
  hyper_draws <- hd$hyper
  # A pinned block folds sigma and rho into the loading's columns rather than
  # carrying them as factors, so the block reports them as 1. Restate the values
  # the fit actually conditioned on -- `field_sd` is unaffected, since the scaled
  # loading already carries them.
  for (nm in intersect(names(fb$pinned), colnames(hyper_draws)))
    hyper_draws[, nm] <- fb$pinned[[nm]]
  hyper_mean  <- colMeans(hyper_draws)
  hyper_sd    <- apply(hyper_draws, 2L, stats::sd)
  # A positive variance component at a few dozen binary sites has a right-skewed
  # posterior, where the mean sits above the bulk; the median is the summary to
  # quote against a known truth.
  hyper_median <- apply(hyper_draws, 2L, stats::median)
  alpha       <- hyper_mean[["alpha"]]

  # Data log-likelihood at the posterior mean (field on psi + alpha*field on
  # cover), so logLik() matches the laplace convention.
  bo   <- par_means[seq_len(p_occ)]
  bp   <- par_means[p_occ + seq_len(p_p)]
  bpos <- par_means[p_occ + p_p + seq_len(p_pos)]
  eta  <- .occu_cover_eta_from_par(model, bo, bp, bpos)
  f_site <- field_mean[site_cell]
  psi_f  <- stats::plogis(.tobs_clamp_eta(stats::qlogis(eta$psi) + f_site))
  ep_f   <- eta$ep_mat + alpha * f_site
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
  if (n_hyper > 0L) {
    ends <- cumsum(vapply(per_chain_draws, nrow, integer(1)))
    starts <- c(1L, ends[-length(ends)] + 1L)
    hyper_chain <- lapply(seq_len(n_chains), function(ch)
      hyper_draws[starts[ch]:ends[ch], fb$sampled, drop = FALSE])
    hyper_diag <- .tobs_nuts_rhat_ess(hyper_chain)
    names(hyper_diag$rhat) <- names(hyper_diag$ess) <- fb$sampled
  }

  # The honest hyper report (gcol33/tulpaObs#204): `sampled_hyper` names the
  # hypers this fit integrated over, `fixed_hyper` those it conditioned on, and
  # `fixed_hyper_values` says at what. `fixed_hyper` is a character vector --
  # empty when nothing is pinned -- so a fit can never claim to have sampled a
  # hyper it did not.
  pinned_nm <- names(fb$pinned) %||% character(0)
  nuts <- list(accept_prob = accept, divergent = divergent, treedepth = treedepth,
               epsilon = epsilon, n_chains = n_chains,
               divergent_total = sum(divergent),
               sigma_beta = sigma.beta, sigma_logdisp = sigma.logdisp,
               field_rho = hyper_mean[["rho"]], field_alpha = alpha,
               field_sigma = hyper_mean[["sigma"]],
               sampled_hyper = fb$sampled,
               fixed_hyper = pinned_nm,
               fixed_hyper_values = fb$pinned,
               hyper_mean = hyper_mean, hyper_median = hyper_median,
               hyper_sd = hyper_sd,
               hyper_rhat = hyper_diag$rhat, hyper_ess = hyper_diag$ess,
               warm_hyper = c(sigma = warm$sigma, rho = warm$rho,
                              alpha = warm$alpha),
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
                        alpha_mean = alpha, alpha_sd = hyper_sd[["alpha"]],
                        rho_mean = hyper_mean[["rho"]],
                        rho_sd = hyper_sd[["rho"]],
                        sampled_hyper = fb$sampled, fixed_hyper = pinned_nm),
    spatial_field = field_mean,
    field_draws   = hd$field,
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
