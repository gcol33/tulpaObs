# distsamp_open.R - open-population distance sampling (unmarked distsampOpen).
# A Dail-Madsen open N-mixture (as dyn_abun(), #37) with a distance-bin
# MULTINOMIAL emission at each primary period instead of the binomial. This is the
# open-population counterpart of the single-season gdistremoval().
#
#   N_{i,1}   ~ Poisson(lambda_i)                                   (initial abundance)
#   N_{i,t}   = Binomial(N_{i,t-1}, omega_i) + Poisson(gamma_i)     (Dail-Madsen "constant")
#   n_{i,t}   ~ Binomial(N_{i,t}, pdist_i)                          (detected total, distance)
#   yDist_{i,t,.} | n_{i,t} ~ Multinomial(cpd_i / pdist_i)         (distance-band allocation)
#
# The band allocation is CONDITIONAL on the period total, independent of the latent
# N, so it factors out of the HMM sum (the gdistremoval trick):
#
#   site_ll = sum_t dmultinom(yDist[i,t,.] ; cpd_i/pdist_i)                 (band allocations)
#           + HMM_forward( y = n_{i,.}, p_t = pdist_i, lambda, omega, gamma )
#
# and the HMM-forward over the latent N sequence with a per-period binomial
# emission dbinom(n_{i,t}; N_{i,t}, pdist_i) is EXACTLY the dyn_abun() marginal, so
# the existing (validated) cpp_dyn_abun_total_log_lik kernel is reused: feed it
# eta_p = logit(pdist_i) (the overall distance detection enters as the detection
# logit) and y = the per-period detected totals. No new HMM kernel; the only new
# pieces are pdist(sigma), the per-period band multinomials, and their assembly.
#
#   log lambda_i = X_lambda_i . beta_lambda    (initial abundance)
#   log sigma_i  = X_sigma_i  . beta_sigma     (distance detection scale)
#   logit omega_i= X_omega_i  . beta_omega     (apparent survival)
#   log gamma_i  = X_gamma_i  . beta_gamma     (recruitment)
#
#   .tobs_build_distsamp_open()   data binder -> model_type = "distsamp_open"
#   .tobs_fit_distsamp_open()     optim over the composed marginal
#   .dispatch_distsamp_open()     tobs() entry
#
# Scope (v1): half-normal key, line / point transect, Poisson initial abundance,
# constant Dail-Madsen dynamics, site-level arms. NB / ZIP initial abundance, other
# dynamics (autoreg / ricker / gompertz), and season-varying sigma are follow-ups.

# ---------------------------------------------------------------------------
# Composed marginal
# ---------------------------------------------------------------------------

# Per-site log-likelihood of the open-population distance model. `theta` packs
# c(beta_lambda, beta_sigma, beta_omega, beta_gamma); returns the summed
# log-likelihood (or a large penalty on an infeasible detection).
# Is this a negative-binomial initial-abundance fit (a trailing `log_r` coord)?
# The base fitter handles Poisson / negbin; zip / zinb take .tobs_fit_distsamp_open_zip.
.dso_use_nb <- function(model) identical(model$mixture %||% "poisson", "negbin")

# ---------------------------------------------------------------------------
# Alternative population dynamics (unmarked distsampOpen tp2..tp5)
# ---------------------------------------------------------------------------

# Per-dynamics active-arm layout + kernel code. The abundance / distance-scale
# arms (lambda / sigma) are always present; the third / fourth arms and how they
# feed the shared dyn_abun HMM kernel vary by dynamics (see src/dyn_abun_kernel.h
# compute_dyn_abun_site_dyn). `code` 0 uses the analytic constant kernel with a
# TIED gamma for notrend; codes 2..5 use the value-only density-dependent kernel.
#   constant : lambda sigma omega(logit surv) gamma(log)                 -> analytic
#   notrend  : lambda sigma omega(logit surv); gamma = (1-omega)*lambda  -> code 0, numeric
#   trend    : lambda sigma gamma(log); no survival                      -> code 3, numeric
#   autoreg  : lambda sigma omega(logit surv) gamma(log)                 -> code 2, numeric
#   ricker   : lambda sigma K(log) r(identity growth)                    -> code 4, numeric
#   gompertz : lambda sigma K(log) r(identity growth)                    -> code 5, numeric
.dso_dyn_meta <- function(dynamics) {
  base <- list(list(name = "lambda", link = "log", src = "lambda"),
               list(name = "sigma",  link = "log", src = "sigma"))
  arms <- switch(dynamics,
    constant = c(base, list(list(name = "omega", link = "logit",    src = "omega"),
                            list(name = "gamma", link = "log",      src = "gamma"))),
    notrend  = c(base, list(list(name = "omega", link = "logit",    src = "omega"))),
    trend    = c(base, list(list(name = "gamma", link = "log",      src = "gamma"))),
    autoreg  = c(base, list(list(name = "omega", link = "logit",    src = "omega"),
                            list(name = "gamma", link = "log",      src = "gamma"))),
    ricker   = c(base, list(list(name = "K",     link = "log",      src = "omega"),
                            list(name = "r",     link = "identity", src = "gamma"))),
    gompertz = c(base, list(list(name = "K",     link = "log",      src = "omega"),
                            list(name = "r",     link = "identity", src = "gamma"))),
    stop("distsamp_open(): unknown dynamics '", dynamics, "'.", call. = FALSE))
  list(arms = arms,
       code = switch(dynamics, autoreg = 2L, trend = 3L, ricker = 4L, gompertz = 5L, 0L),
       analytic   = identical(dynamics, "constant"),
       tied_gamma = identical(dynamics, "notrend"),
       has_omega  = !identical(dynamics, "trend"))
}

# Per-site band-allocation log-likelihood (the distance multinomial, a function of
# sigma only), shared by every dynamics. `cpd` = [n x Jbin] cell probs, pdist =
# rowSums. Returns a length-n_sites vector.
.dso_band_ll_vec <- function(cpd, pdist, model) {
  pid  <- cpd / pdist
  band <- numeric(model$n_sites)
  for (t in seq_len(model$n_seasons)) {
    yb   <- model$y[, , t]
    band <- band + lgamma(rowSums(yb) + 1) - rowSums(lgamma(yb + 1)) +
            rowSums(yb * log(pmax(pid, 1e-300)))
  }
  band
}

# Map a coefficient vector to the kernel eta inputs for an alternative-dynamics
# fit. `model$dyn_arms` is the active-arm list; each arm's design is X_processes
# in the same order. Returns eta_lambda, sigma, the kernel omega / gamma slots
# (raw eta -- the kernel applies each dynamics' link), and pdist / cpd.
.dso_dyn_eta <- function(theta, model) {
  meta <- model$dyn_meta; arms <- meta$arms
  po <- 0L; eta <- vector("list", length(arms)); names(eta) <- vapply(arms, `[[`, "", "name")
  for (a in seq_along(arms)) {
    X <- model$X_processes[[a]]; p <- ncol(X)
    eta[[a]] <- as.vector(X %*% theta[po + seq_len(p)]); po <- po + p
  }
  eta_lambda <- eta[["lambda"]]
  sigma <- exp(eta[["sigma"]])
  cpd   <- .gdr_dist_cp(sigma, model$cutpoints, model$transect)
  pdist <- rowSums(cpd)
  # Kernel omega slot: survival-eta (constant/notrend/autoreg), log-K eta
  # (ricker/gompertz), or 0 (trend, no survival term).
  eta_omega <- if (!meta$has_omega) rep(0, model$n_sites)
               else if (!is.null(eta[["omega"]])) eta[["omega"]] else eta[["K"]]
  # Kernel gamma slot: tied for notrend, else the free gamma / r arm.
  eta_gamma <- if (meta$tied_gamma)
                 log(pmax((1 - stats::plogis(eta[["omega"]])) * exp(eta_lambda), 1e-300))
               else if (!is.null(eta[["gamma"]])) eta[["gamma"]] else eta[["r"]]
  list(eta_lambda = eta_lambda, sigma = sigma, cpd = cpd, pdist = pdist,
       eta_omega = eta_omega, eta_gamma = eta_gamma)
}

# Negative log-likelihood for an alternative-dynamics fit (numeric-gradient BFGS).
.dso_dyn_negll <- function(theta, model) {
  e <- .dso_dyn_eta(theta, model)
  if (any(!is.finite(e$pdist)) || any(e$pdist <= 1e-8) || any(e$pdist >= 1 - 1e-10))
    return(1e10)
  pd  <- pmin(pmax(e$pdist, 1e-10), 1 - 1e-10)
  code <- model$dyn_meta$code
  hmm <- if (code == 0L)
    cpp_dyn_abun_total_log_lik(model$y_flat, model$n_sites, model$n_seasons, 1L,
      model$K_max, e$eta_lambda, stats::qlogis(pd), e$eta_omega, e$eta_gamma,
      use_nb = FALSE, eta_logr = 0)$log_lik
  else
    cpp_dyn_abun_dynamics_log_lik(model$y_flat, model$n_sites, model$n_seasons, 1L,
      model$K_max, e$eta_lambda, stats::qlogis(pd), e$eta_omega, e$eta_gamma,
      dynamics = code, use_nb = FALSE, eta_logr = 0)$log_lik
  val <- -(hmm + sum(.dso_band_ll_vec(e$cpd, pd, model)))
  if (is.finite(val)) val else 1e10
}

# Fitter for the alternative-dynamics path: numeric-gradient BFGS over the exact
# forward-HMM + band marginal, with a numeric-Hessian observed-information vcov.
# `model` carries the active-arm designs + `dyn_meta` (built by the binder).
.tobs_fit_distsamp_open_dyn <- function(model, verbose = TRUE, ...) {
  meta <- model$dyn_meta
  ps   <- vapply(model$process_info, function(pp) pp$p, integer(1))
  # Moment init: detected-total scale for lambda / sigma; neutral for the rest.
  ntot <- model$ntot; sig0 <- stats::median(model$cutpoints[-1])
  pd0  <- sum(.gdr_dist_cp(sig0, model$cutpoints, model$transect))
  lam0 <- mean(ntot[, 1L]) / max(pd0, 0.05)
  init <- numeric(sum(ps)); off <- cumsum(c(0L, ps))
  for (a in seq_along(meta$arms)) {
    nm <- meta$arms[[a]]$name
    init[off[a] + 1L] <- switch(nm,
      lambda = log(max(lam0, 1e-2)), sigma = log(max(sig0, 1e-2)),
      omega  = stats::qlogis(0.6), gamma = log(max(mean(ntot), 1)),
      K = log(max(2 * lam0, 5)), r = 0.1, 0)
  }

  par_names <- unlist(lapply(model$process_info, function(pp)
    paste0(pp$name, "_", pp$coef_names)))

  .tobs_bfgs_marginal_fit(
    function(theta) .dso_dyn_negll(theta, model),
    init, par_names, model, N = model$n_sites,
    control = list(maxit = 500L, reltol = 1e-9),
    extra = function(means) list(mixture = "poisson",
                                 dynamics = model$dynamics,
                                 zero_inflated = FALSE))
}

.dso_negll <- function(theta, model) {
  Xl <- model$X_processes[[1L]]; Xs <- model$X_processes[[2L]]
  Xo <- model$X_processes[[3L]]; Xg <- model$X_processes[[4L]]
  pl <- ncol(Xl); ps <- ncol(Xs); po <- ncol(Xo); pg <- ncol(Xg)
  use_nb <- .dso_use_nb(model); ir <- pl + ps + po + pg + 1L
  bl <- theta[seq_len(pl)]
  bs <- theta[pl + seq_len(ps)]
  bo <- theta[pl + ps + seq_len(po)]
  bg <- theta[pl + ps + po + seq_len(pg)]
  elr <- if (use_nb) theta[ir] else 0

  sigma <- exp(as.vector(Xs %*% bs))
  cpd   <- .gdr_dist_cp(sigma, model$cutpoints, model$transect)   # [n x Jbin]
  pdist <- rowSums(cpd)
  if (any(!is.finite(pdist)) || any(pdist <= 1e-8) || any(pdist >= 1 - 1e-10))
    return(1e10)

  ev <- cpp_dyn_abun_total_log_lik(
    model$y_flat, model$n_sites, model$n_seasons, 1L, model$K_max,
    as.vector(Xl %*% bl),                 # eta_lambda (log link)
    stats::qlogis(pdist),                 # eta_p = logit(pdist)  (the reuse)
    as.vector(Xo %*% bo),                 # eta_omega (logit link)
    as.vector(Xg %*% bg),                 # eta_gamma (log link)
    use_nb = use_nb, eta_logr = elr)
  hmm <- sum(ev$log_lik)

  # Per-period distance-band multinomials (a function of sigma only).
  pid  <- cpd / pdist                                             # [n x Jbin]
  band <- 0
  for (t in seq_len(model$n_seasons)) {
    yb   <- model$y[, , t]                                        # [n x Jbin]
    band <- band + sum(lgamma(rowSums(yb) + 1) - rowSums(lgamma(yb + 1)) +
                       rowSums(yb * log(pmax(pid, 1e-300))))
  }
  val <- -(hmm + band)
  if (is.finite(val)) val else 1e10
}

# Analytic gradient of the LOG-likelihood wrt theta = c(beta_lambda, beta_sigma,
# beta_omega, beta_gamma). The dyn_abun kernel returns grad_eta_{lambda,p,omega,
# gamma} in the same call as the value, so lambda / omega / gamma chain straight
# through their designs. sigma enters only via pdist (the detection logit) and the
# band multinomials; the distance-integral derivative dcpd/dsigma is a cheap
# central finite difference on .gdr_dist_cp (NO extra HMM evaluation), so the whole
# gradient costs one kernel call. Returns a length-n_par numeric vector.
.dso_grad <- function(theta, model) {
  Xl <- model$X_processes[[1L]]; Xs <- model$X_processes[[2L]]
  Xo <- model$X_processes[[3L]]; Xg <- model$X_processes[[4L]]
  pl <- ncol(Xl); ps <- ncol(Xs); po <- ncol(Xo); pg <- ncol(Xg)
  use_nb <- .dso_use_nb(model); ir <- pl + ps + po + pg + 1L
  bl <- theta[seq_len(pl)]
  bs <- theta[pl + seq_len(ps)]
  bo <- theta[pl + ps + seq_len(po)]
  bg <- theta[pl + ps + po + seq_len(pg)]
  elr <- if (use_nb) theta[ir] else 0

  sigma <- exp(as.vector(Xs %*% bs))
  cpd   <- .gdr_dist_cp(sigma, model$cutpoints, model$transect)
  pdist <- rowSums(cpd)
  pd    <- pmin(pmax(pdist, 1e-10), 1 - 1e-10)

  ev <- cpp_dyn_abun_total_log_lik(
    model$y_flat, model$n_sites, model$n_seasons, 1L, model$K_max,
    as.vector(Xl %*% bl), stats::qlogis(pd),
    as.vector(Xo %*% bo), as.vector(Xg %*% bg), use_nb = use_nb, eta_logr = elr)

  # Cheap central FD of the distance cell probs wrt sigma (distance integral only).
  h    <- 1e-5 * pmax(sigma, 1)
  cpp1 <- .gdr_dist_cp(sigma + h, model$cutpoints, model$transect)
  cpm1 <- .gdr_dist_cp(sigma - h, model$cutpoints, model$transect)
  dcpd <- (cpp1 - cpm1) / (2 * h)                    # [n x Jbin]
  dpd  <- rowSums(dcpd)                              # dpdist/dsigma

  # HMM contribution to dL/dsigma via eta_p = logit(pdist).
  dL_dsig <- ev$grad_eta_p / (pd * (1 - pd)) * dpd
  # Band contribution: sum_b Y_b d(log cpd_b)/dsigma - Ntot d(log pdist)/dsigma,
  # with Y_b the band totals across periods and Ntot the total detected per site.
  Yb   <- apply(model$y, c(1L, 2L), sum)             # [n x Jbin] band totals
  Ntot <- rowSums(Yb)
  dL_dsig <- dL_dsig +
    rowSums(Yb / pmax(cpd, 1e-300) * dcpd) - Ntot / pd * dpd

  g <- numeric(pl + ps + po + pg + if (use_nb) 1L else 0L)
  g[seq_len(pl)]                <- as.numeric(crossprod(Xl, ev$grad_eta_lambda))
  g[pl + seq_len(ps)]           <- as.numeric(crossprod(Xs, dL_dsig * sigma))
  g[pl + ps + seq_len(po)]      <- as.numeric(crossprod(Xo, ev$grad_eta_omega))
  g[pl + ps + po + seq_len(pg)] <- as.numeric(crossprod(Xg, ev$grad_eta_gamma))
  # log_r score is returned summed across sites (the one dispersion coordinate).
  if (use_nb) g[ir] <- as.numeric(ev$grad_eta_logr)
  g
}

# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# `y` is a 3D array [n_sites x n_bins x n_seasons] of per-distance-band counts at
# each primary period. abund_formula = log lambda, det_formula = log sigma,
# omega_formula = logit survival, gamma_formula = log recruitment.
.tobs_build_distsamp_open <- function(abund_formula, det_formula, omega_formula,
                                      gamma_formula, data, y, cutpoints, transect,
                                      K_max = NULL, mixture = "poisson",
                                      dynamics = "constant") {
  if (length(dim(y)) != 3L) {
    stop("distsamp_open() y must be a 3D array [n_sites x n_bins x n_seasons] of ",
         "per-distance-band counts.", call. = FALSE)
  }
  if (any(y < 0 | y != round(y), na.rm = TRUE)) {
    stop("distsamp_open() y must be non-negative integer counts.", call. = FALSE)
  }
  n_sites   <- dim(y)[1L]; n_bins <- dim(y)[2L]; n_seasons <- dim(y)[3L]
  if (n_seasons < 2L) {
    stop("distsamp_open() needs >= 2 primary periods (n_seasons = dim(y)[3]); ",
         "for a single period use distance().", call. = FALSE)
  }
  if (length(cutpoints) != n_bins + 1L) {
    stop(sprintf(paste0("distsamp_open() cutpoints must have length dim(y)[2] + 1 ",
         "= %d (the distance-bin edges)."), n_bins + 1L), call. = FALSE)
  }
  cutpoints <- as.numeric(cutpoints)
  if (any(diff(cutpoints) <= 0) || cutpoints[1] < 0) {
    stop("distsamp_open() cutpoints must be strictly increasing and start >= 0.",
         call. = FALSE)
  }
  .tobs_check_site_count(n_sites, nrow(data), "sites")
  storage.mode(y) <- "integer"

  # Per-period detected totals feed the HMM kernel (secondary occasions absorbed,
  # J = 1). Layout matches cpp_dyn_abun_total_log_lik: aperm(y[,1,,], c(2,3,1)).
  ntot   <- apply(y, c(1L, 3L), sum)                     # [n_sites x n_seasons]
  y_kern <- array(as.integer(ntot), c(n_sites, 1L, n_seasons))
  y_flat <- as.integer(aperm(y_kern, c(2L, 3L, 1L)))
  # The abundance-HMM forward is cubic in K_max, so the truncation is chosen from
  # the DETECTION-CORRECTED abundance scale, not a blunt multiple of the detected
  # total: max_i N_i ~ max(ntot) / pdist, plus a few Poisson SDs of headroom. A
  # rough pdist at the median-cutpoint sigma sets the scale (the fit re-integrates
  # the true pdist per iteration; the truncation only needs to bound N's tail).
  if (is.null(K_max)) {
    pd0  <- sum(.gdr_dist_cp(stats::median(cutpoints[-1]), cutpoints, transect))
    maxN <- max(ntot) / max(pd0, 0.1)
    # The detection-corrected maxN already reflects the realised (heavier) upper
    # tail of an NB / zero-inflated initial abundance through max(ntot); a modest
    # extra headroom bounds the marginal tail without inflating the cubic-in-K cost.
    infl  <- if (mixture %in% c("negbin", "zinb", "zip")) 5 else 4
    K_max <- as.integer(ceiling(maxN + infl * sqrt(maxN) + 10))
  }
  K_max <- as.integer(K_max)

  bind <- .tobs_bind_formulas(
    list(lambda = abund_formula, sigma = det_formula,
         omega = omega_formula, gamma = gamma_formula), data)

  base <- list(
    model_type  = "distsamp_open",
    y           = y,
    y_flat      = y_flat,
    ntot        = ntot,
    cutpoints   = cutpoints,
    transect    = transect,
    n_bins      = n_bins,
    n_seasons   = n_seasons,
    K_max       = K_max,
    mixture     = mixture,
    dynamics    = dynamics,
    structured_terms = bind$terms,
    data        = data,
    n_sites     = n_sites)

  if (identical(dynamics, "constant")) {
    # Constant Dail-Madsen: the four standard arms + the analytic-gradient path
    # (kept byte-identical; the NB / ZIP mixtures layer over this).
    X_lambda <- stats::model.matrix(bind$fe$lambda, data)
    X_sigma  <- stats::model.matrix(bind$fe$sigma, data)
    X_omega  <- stats::model.matrix(bind$fe$omega, data)
    X_gamma  <- stats::model.matrix(bind$fe$gamma, data)
    model <- c(base, list(
      X_processes = list(X_lambda, X_sigma, X_omega, X_gamma),
      formulas    = list(lambda = bind$fe$lambda, sigma = bind$fe$sigma,
                         omega = bind$fe$omega, gamma = bind$fe$gamma),
      process_info = list(
        list(name = "lambda", p = ncol(X_lambda),
             coef_names = colnames(X_lambda), link = "log"),
        list(name = "sigma", p = ncol(X_sigma),
             coef_names = colnames(X_sigma), link = "log"),
        list(name = "omega", p = ncol(X_omega),
             coef_names = colnames(X_omega), link = "logit"),
        list(name = "gamma", p = ncol(X_gamma),
             coef_names = colnames(X_gamma), link = "log"))))
    return(structure(model, class = "tobs_model"))
  }

  # Alternative dynamics: build only the ACTIVE arms per the dynamics layout, each
  # sourced from the matching formula (lambda<-formula, sigma<-detection,
  # omega/K<-omega=~, gamma/r<-gamma=~). The kernel-slot links are applied inside
  # the dyn kernel; here process_info records the reported name / link.
  meta <- .dso_dyn_meta(dynamics)
  fe_by_src <- list(lambda = bind$fe$lambda, sigma = bind$fe$sigma,
                    omega = bind$fe$omega, gamma = bind$fe$gamma)
  X_list <- vector("list", length(meta$arms))
  pinfo  <- vector("list", length(meta$arms))
  formulas <- list()
  for (a in seq_along(meta$arms)) {
    arm <- meta$arms[[a]]
    f   <- fe_by_src[[arm$src]]
    X   <- stats::model.matrix(f, data)
    X_list[[a]] <- X
    pinfo[[a]]  <- list(name = arm$name, p = ncol(X),
                        coef_names = colnames(X), link = arm$link)
    formulas[[arm$name]] <- f
  }
  model <- c(base, list(
    X_processes  = X_list,
    formulas     = formulas,
    process_info = pinfo,
    dyn_meta     = meta))
  structure(model, class = "tobs_model")
}

.dso_unpack <- function(theta, model) {
  Xl <- model$X_processes[[1L]]; Xs <- model$X_processes[[2L]]
  Xo <- model$X_processes[[3L]]; Xg <- model$X_processes[[4L]]
  pl <- ncol(Xl); ps <- ncol(Xs); po <- ncol(Xo); pg <- ncol(Xg)
  mix    <- model$mixture %||% "poisson"
  use_nb <- mix %in% c("negbin", "zinb")
  is_zi  <- mix %in% c("zip", "zinb")
  out <- list(
    lambda = exp(as.vector(Xl %*% theta[seq_len(pl)])),
    sigma  = exp(as.vector(Xs %*% theta[pl + seq_len(ps)])),
    omega  = stats::plogis(as.vector(Xo %*% theta[pl + ps + seq_len(po)])),
    gamma  = exp(as.vector(Xg %*% theta[pl + ps + po + seq_len(pg)]))
  )
  base <- pl + ps + po + pg
  if (is_zi) { out$zi_logit <- theta[base + 1L]; out$zi_omega <- stats::plogis(out$zi_logit) }
  if (use_nb) { out$log_r <- theta[base + if (is_zi) 2L else 1L]; out$r <- exp(out$log_r) }
  out
}

# ---------------------------------------------------------------------------
# Fitter
# ---------------------------------------------------------------------------

.tobs_fit_distsamp_open <- function(model, verbose = TRUE, ...) {
  use_nb <- .dso_use_nb(model)
  ps <- vapply(model$process_info, function(pp) pp$p, integer(1))
  # Moment init: mean detected total seeds lambda given a rough detection; the
  # distance scale seeds at the median cutpoint; moderate survival / recruitment.
  ntot   <- model$ntot
  sig0   <- stats::median(model$cutpoints[-1])
  pdist0 <- sum(.gdr_dist_cp(sig0, model$cutpoints, model$transect))
  lam0   <- mean(ntot[, 1L]) / max(pdist0, 0.05)
  init <- c(log(max(lam0, 1e-2)), rep(0, ps[1] - 1L),
            log(max(sig0, 1e-2)),  rep(0, ps[2] - 1L),
            stats::qlogis(0.6),     rep(0, ps[3] - 1L),
            log(max(mean(ntot), 1)), rep(0, ps[4] - 1L))
  if (use_nb) init <- c(init, log(2))   # trailing log_r seed (moderate overdispersion)

  # Analytic gradient (one HMM kernel call per evaluation) drives BFGS; the
  # observed information is the FD-Jacobian of the negative gradient at the mode
  # (2*n_par cheap gradient calls), as in fp_occu() -- far cheaper than a numeric
  # Hessian of the value over the cubic-in-K forward recursion.
  par_names <- unlist(lapply(model$process_info, function(pp)
    paste0(pp$name, "_", pp$coef_names)))
  if (use_nb) par_names <- c(par_names, "log_r")

  .tobs_bfgs_marginal_fit(
    function(theta) .dso_negll(theta, model), init, par_names, model,
    N = model$n_sites,
    gr = function(theta) -.dso_grad(theta, model),
    extra = function(means) list(
      mixture       = model$mixture,
      log_r         = if (use_nb) unname(means[["log_r"]]) else NULL,
      r             = if (use_nb) exp(unname(means[["log_r"]])) else NULL,
      zero_inflated = FALSE))
}

# ---------------------------------------------------------------------------
# Zero-inflated open-population distance sampling (ZIP / ZINB)
# ---------------------------------------------------------------------------

# A structural-zero site is never occupied across any primary period, so all its
# distance-band counts are zero. The observed per-site marginal is the same
# two-component mixture the static / dyn_abun ZIP layers use, over the composed
# open-population distance marginal (forward-HMM + band multinomials):
#   L_i = omega * 1{all y_i = 0} + (1 - omega) * L_open_i,
# with L_open_i = exp(hmm_i + band_i) the exact marginal. On an all-zero site the
# band multinomial contributes 0 (ysum = 0), so L_open_i reduces to the HMM part
# there. This is a PURE-R layer over the C++ per-site marginal; the Poisson /
# negbin paths are untouched. omega is an intercept-only structural-zero
# probability (logit), named `zi_logit` (distinct from `omega`, the survival arm).
#
# Scope (v1): non-spatial laplace only, intercept-only zi. The additive marginal
# and its per-site gradient are shared verbatim with the .dso_negll / .dso_grad
# base fitter (weighted by the structural-zero posterior w_i).
.tobs_fit_distsamp_open_zip <- function(model, verbose = TRUE, ...) {
  is_nb <- identical(model$mixture, "zinb")
  Xl <- model$X_processes[[1L]]; Xs <- model$X_processes[[2L]]
  Xo <- model$X_processes[[3L]]; Xg <- model$X_processes[[4L]]
  pl <- ncol(Xl); ps <- ncol(Xs); po <- ncol(Xo); pg <- ncol(Xg)
  cutpoints <- model$cutpoints; transect <- model$transect
  N <- model$n_sites; T <- model$n_seasons; K <- model$K_max

  idx <- list(lambda = seq_len(pl), sigma = pl + seq_len(ps),
              omega = pl + ps + seq_len(po), gamma = pl + ps + po + seq_len(pg))
  izi <- pl + ps + po + pg + 1L
  ir  <- if (is_nb) izi + 1L else NA_integer_

  # Per-site all-band-zero indicator over the whole [bin x season] block.
  Yb  <- apply(model$y, c(1L, 2L), sum)          # [n x Jbin] band totals
  az  <- rowSums(Yb) == 0

  # One composed-marginal evaluation: the open-population per-site log-lik (HMM +
  # band) AND the per-site eta / sigma gradient pieces, so BFGS runs on analytic
  # gradients weighted by the structural-zero posterior.
  eval_open <- function(theta) {
    sigma <- exp(as.vector(Xs %*% theta[idx$sigma]))
    cpd   <- .gdr_dist_cp(sigma, cutpoints, transect)
    pdist <- rowSums(cpd)
    if (any(!is.finite(pdist)) || any(pdist <= 1e-8) || any(pdist >= 1 - 1e-10))
      return(NULL)
    pd  <- pmin(pmax(pdist, 1e-10), 1 - 1e-10)
    ev  <- tryCatch(cpp_dyn_abun_total_log_lik(
      model$y_flat, N, T, 1L, K,
      as.vector(Xl %*% theta[idx$lambda]), stats::qlogis(pd),
      as.vector(Xo %*% theta[idx$omega]),  as.vector(Xg %*% theta[idx$gamma]),
      use_nb = is_nb, eta_logr = if (is_nb) theta[ir] else 0),
      error = function(e) NULL)
    if (is.null(ev) || any(!is.finite(ev$log_lik_site))) return(NULL)
    # Per-period band multinomials (a function of sigma only), per site.
    pid  <- cpd / pdist
    band <- numeric(N)
    for (t in seq_len(T)) {
      yb   <- model$y[, , t]
      band <- band + lgamma(rowSums(yb) + 1) - rowSums(lgamma(yb + 1)) +
              rowSums(yb * log(pmax(pid, 1e-300)))
    }
    list(ev = ev, cpd = cpd, pdist = pd, sigma = sigma, band = band,
         llr = as.numeric(ev$log_lik_site) + band)
  }

  # ZIP marginal log-lik + posterior "not-a-structural-zero" weight w_i.
  zip_pieces <- function(theta, op) {
    om <- stats::plogis(theta[izi]); log1m <- log1p(-om)
    llr <- op$llr; ll <- numeric(N); w <- rep(1, N)
    ll[!az] <- log1m + llr[!az]
    a <- log1m + llr[az]; b <- log(om); mx <- pmax(a, b)
    Li <- mx + log(exp(a - mx) + exp(b - mx))
    ll[az] <- Li; w[az] <- exp(a - Li)
    list(log_lik = sum(ll), w = w, om = om)
  }

  neg_ll <- function(theta) {
    op <- eval_open(theta); if (is.null(op)) return(1e10)
    val <- -zip_pieces(theta, op)$log_lik
    if (is.finite(val)) val else 1e10
  }
  neg_grad <- function(theta) {
    op <- eval_open(theta); if (is.null(op)) return(rep(0, length(theta)))
    zp <- zip_pieces(theta, op); w <- zp$w; om <- zp$om; ev <- op$ev
    # sigma enters via eta_p (pdist) AND the band multinomials; the per-site
    # dL/dsigma mirrors .dso_grad, weighted by the structural-zero posterior w_i.
    sig <- op$sigma; cpd <- op$cpd; pd <- op$pdist
    h   <- 1e-5 * pmax(sig, 1)
    dcpd <- (.gdr_dist_cp(sig + h, cutpoints, transect) -
             .gdr_dist_cp(sig - h, cutpoints, transect)) / (2 * h)
    dpd  <- rowSums(dcpd)
    dL_dsig <- ev$grad_eta_p / (pd * (1 - pd)) * dpd +
      rowSums(Yb / pmax(cpd, 1e-300) * dcpd) - rowSums(Yb) / pd * dpd
    g <- numeric(length(theta))
    g[idx$lambda] <- as.numeric(crossprod(Xl, w * ev$grad_eta_lambda))
    g[idx$sigma]  <- as.numeric(crossprod(Xs, w * dL_dsig * sig))
    g[idx$omega]  <- as.numeric(crossprod(Xo, w * as.numeric(ev$grad_eta_omega)))
    g[idx$gamma]  <- as.numeric(crossprod(Xg, w * as.numeric(ev$grad_eta_gamma)))
    g[izi] <- sum((1 - om) - w)
    # The NB log_r score is returned summed (not per-site), so it cannot be
    # ZIP-weighted analytically; central-difference just this coordinate.
    if (is_nb) {
      th <- theta; hh <- 1e-4
      th[ir] <- theta[ir] + hh; fp <- -neg_ll(th)
      th[ir] <- theta[ir] - hh; fm <- -neg_ll(th)
      g[ir] <- (fp - fm) / (2 * hh)
    }
    -g
  }

  # Warm start: lambda from the detected-site scale, moderate dynamics, the
  # structural-zero logit from a modest share of the all-zero sites.
  theta0 <- numeric(if (is_nb) ir else izi)
  sig0   <- stats::median(cutpoints[-1])
  pd0    <- sum(.gdr_dist_cp(sig0, cutpoints, transect))
  nz1    <- model$ntot[!az, 1L]
  theta0[idx$lambda[1]] <- log(max(mean(if (length(nz1)) nz1 else model$ntot[, 1L]) /
                                     max(pd0, 0.05), 0.5))
  theta0[idx$sigma[1]]  <- log(max(sig0, 1e-2))
  theta0[idx$omega[1]]  <- stats::qlogis(0.6)
  theta0[idx$gamma[1]]  <- log(max(mean(model$ntot), 1))
  theta0[izi]           <- stats::qlogis(min(max(mean(az) * 0.5, 0.05), 0.7))
  if (is_nb) theta0[ir] <- log(2)

  par_names <- c(paste0("lambda_", colnames(Xl)), paste0("sigma_", colnames(Xs)),
                 paste0("omega_", colnames(Xo)), paste0("gamma_", colnames(Xg)),
                 "zi_logit")
  if (is_nb) par_names <- c(par_names, "log_r")

  .tobs_bfgs_marginal_fit(
    neg_ll, theta0, par_names, model, N = model$n_sites, gr = neg_grad,
    control = list(maxit = 500L, reltol = 1e-8),
    extra = function(means) list(
      mixture       = model$mixture,
      zero_inflated = TRUE,
      zi_logit      = unname(means[["zi_logit"]]),
      zi_omega      = stats::plogis(unname(means[["zi_logit"]])),
      log_r         = if (is_nb) unname(means[["log_r"]]) else NULL,
      r             = if (is_nb) exp(unname(means[["log_r"]])) else NULL))
}

# ---------------------------------------------------------------------------
# tobs() dispatcher
# ---------------------------------------------------------------------------

.dispatch_distsamp_open <- function(formula, data, family, detection, y, visits,
                                    engine, priors, control,
                                    approx = "gaussian_laplace",
                                    correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection))
    stop("distsamp_open() requires a `detection` formula (the site-level ",
         "log-sigma distance-scale model).", call. = FALSE)
  if (is.null(y))
    stop("distsamp_open() requires `y` (a 3D array [n_sites x n_bins x ",
         "n_seasons] of per-distance-band counts).", call. = FALSE)
  if (!is.null(visits))
    stop("distsamp_open() detection is site-level; visit-level covariates ",
         "(`visits`) are not yet supported.", call. = FALSE)
  cutpoints <- family$params$cutpoints
  if (is.null(cutpoints))
    stop("distsamp_open() requires `cutpoints` (the distance-bin edges); pass ",
         "distsamp_open(cutpoints = ...).", call. = FALSE)
  mixture  <- family$params$mixture %||% "poisson"
  dynamics <- family$params$dynamics %||% "constant"
  if (!identical(dynamics, "constant") && !identical(mixture, "poisson")) {
    stop(sprintf(paste0("distsamp_open(dynamics = \"%s\") is Poisson-only for now; ",
         "the negbin / zero-inflated initial abundance is layered on the constant ",
         "Dail-Madsen marginal. Use dynamics = \"constant\" with mixture = \"%s\", ",
         "or mixture = \"poisson\" with the alternative dynamics."),
         dynamics, mixture), call. = FALSE)
  }
  model <- .tobs_build_distsamp_open(
    abund_formula = formula, det_formula = detection,
    omega_formula = dots$omega %||% ~1, gamma_formula = dots$gamma %||% ~1,
    data = data, y = y, cutpoints = cutpoints,
    transect = family$params$transect, K_max = family$params$K_max,
    mixture = mixture, dynamics = dynamics)
  .tobs_reject_unwired_structs(
    model, "distsamp_open()",
    hint = paste0("the open-population distance marginal is fitted on fixed ",
                  "effects only, so drop the term"))
  if (!identical(dynamics, "constant"))
    .tobs_fit_distsamp_open_dyn(model, verbose = isTRUE(control$verbose))
  else if (mixture %in% c("zip", "zinb"))
    .tobs_fit_distsamp_open_zip(model, verbose = isTRUE(control$verbose))
  else
    .tobs_fit_distsamp_open(model, verbose = isTRUE(control$verbose))
}

# ---------------------------------------------------------------------------
# fitted / predict / residuals / pointwise log-likelihood / simulate
# ---------------------------------------------------------------------------

# fitted(): per-site initial abundance lambda, distance scale sigma, survival
# omega, recruitment gamma, and the overall distance detection pdist.
.tobs_fitted_distsamp_open <- function(object) {
  model <- object$model
  if (!identical(model$dynamics %||% "constant", "constant")) {
    e <- .dso_dyn_eta(object$means, model)
    out <- list(lambda = exp(e$eta_lambda), sigma = e$sigma, pdist = e$pdist)
    # Report the dynamics arm(s) on their natural scale (from the raw kernel slots):
    # omega -> survival prob, K -> carrying capacity, gamma -> recruitment rate,
    # r -> growth rate (identity).
    for (arm in model$dyn_meta$arms) {
      if (arm$name %in% c("lambda", "sigma")) next
      out[[arm$name]] <- switch(arm$name,
        omega  = stats::plogis(e$eta_omega),
        K      = exp(e$eta_omega),
        gamma  = exp(e$eta_gamma),
        r      = e$eta_gamma)
    }
    return(out)
  }
  up <- .dso_unpack(object$means, model)
  pdist <- rowSums(.gdr_dist_cp(up$sigma, model$cutpoints, model$transect))
  list(lambda = up$lambda, sigma = up$sigma, omega = up$omega, gamma = up$gamma,
       pdist = pdist, r = up$r %||% NULL, zi_omega = up$zi_omega %||% NULL)
}

# predict(): abundance (default), distance scale, survival, or recruitment. For an
# alternative-dynamics fit the requested type maps to the active arm by name
# (survival -> omega, recruitment -> gamma / r, and "K" for the carrying capacity);
# a type whose arm is absent for that dynamics errors with the available set.
.tobs_predict_distsamp_open <- function(object, newdata = NULL,
                                        type = c("abundance", "distance",
                                                 "survival", "recruitment")) {
  type  <- match.arg(type)
  model <- object$model; b <- object$means
  p <- vapply(model$process_info, function(pp) pp$p, integer(1))
  off <- cumsum(c(0L, p))
  arm_names <- vapply(model$process_info, function(pp) pp$name, character(1))
  want <- switch(type, abundance = "lambda", distance = "sigma",
                 survival = "omega", recruitment = "gamma")
  arm <- match(want, arm_names)
  # For ricker / gompertz the third / fourth arms are K / r; map the generic
  # "survival" / "recruitment" requests onto them.
  if (is.na(arm) && identical(type, "survival"))    arm <- match("K", arm_names)
  if (is.na(arm) && identical(type, "recruitment")) arm <- match("r", arm_names)
  if (is.na(arm)) {
    stop(sprintf(paste0("distsamp_open(dynamics = \"%s\"): predict type '%s' has no arm; ",
                        "available: %s."),
         model$dynamics %||% "constant", type,
         paste(arm_names, collapse = ", ")), call. = FALSE)
  }
  X <- if (is.null(newdata)) model$X_processes[[arm]]
       else stats::model.matrix(model$formulas[[arm_names[arm]]], newdata)
  eta <- as.vector(X %*% b[off[arm] + seq_len(p[arm])])
  link <- model$process_info[[arm]]$link
  if (identical(link, "logit")) stats::plogis(eta)
  else if (identical(link, "identity")) eta
  else exp(eta)
}

# residuals(): per-site Pearson / deviance on the first-period detected total
# against its expected value lambda * pdist.
.tobs_residuals_distsamp_open <- function(object, type) {
  fv   <- .tobs_fitted_distsamp_open(object)
  m1   <- fv$lambda * fv$pdist
  obs  <- object$model$ntot[, 1L]
  eps  <- 1e-10
  res <- switch(type,
    response = obs - m1,
    pearson  = (obs - m1) / sqrt(m1 + eps),
    deviance = sign(obs - m1) * sqrt(2 * abs(
      ifelse(obs > 0, obs * log(obs / pmax(m1, eps)), 0) - (obs - m1))))
  list(occ = res, det = NULL)
}

# Pointwise (per-site) log-likelihood [n_draws x n_sites] over posterior draws.
.tobs_ploglik_distsamp_open <- function(object, n.draws = 1000L, n.threads = 1L) {
  model <- object$model
  draws <- object$draws
  if (!is.null(n.draws) && as.integer(n.draws) < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  t(vapply(seq_len(nrow(draws)), function(d)
    .dso_site_loglik(draws[d, ], model), numeric(model$n_sites)))
}

# Per-site log-likelihood vector at one coefficient draw (for WAIC / ploglik).
# The open-population marginal per site is the exact forward-HMM log-lik
# (`log_lik_site`, NOT the summed `log_lik`) plus the per-period band multinomials;
# a zero-inflated fit mixes the structural-zero point mass over the all-zero sites.
.dso_site_loglik <- function(theta, model) {
  if (!identical(model$dynamics %||% "constant", "constant")) {
    e  <- .dso_dyn_eta(theta, model)
    pd <- pmin(pmax(e$pdist, 1e-10), 1 - 1e-10)
    code <- model$dyn_meta$code
    ev <- if (code == 0L)
      cpp_dyn_abun_total_log_lik(model$y_flat, model$n_sites, model$n_seasons, 1L,
        model$K_max, e$eta_lambda, stats::qlogis(pd), e$eta_omega, e$eta_gamma,
        use_nb = FALSE, eta_logr = 0)
    else
      cpp_dyn_abun_dynamics_log_lik(model$y_flat, model$n_sites, model$n_seasons, 1L,
        model$K_max, e$eta_lambda, stats::qlogis(pd), e$eta_omega, e$eta_gamma,
        dynamics = code, use_nb = FALSE, eta_logr = 0)
    return(as.numeric(ev$log_lik_site) + .dso_band_ll_vec(e$cpd, pd, model))
  }
  up    <- .dso_unpack(theta, model)
  mix   <- model$mixture %||% "poisson"
  use_nb <- mix %in% c("negbin", "zinb")
  cpd   <- .gdr_dist_cp(up$sigma, model$cutpoints, model$transect)
  pdist <- rowSums(cpd)
  pdist <- pmin(pmax(pdist, 1e-10), 1 - 1e-10)
  ev <- cpp_dyn_abun_total_log_lik(
    model$y_flat, model$n_sites, model$n_seasons, 1L, model$K_max,
    log(up$lambda), stats::qlogis(pdist), stats::qlogis(up$omega),
    log(up$gamma), use_nb = use_nb, eta_logr = if (use_nb) up$log_r else 0)
  pid  <- cpd / pdist
  band <- numeric(model$n_sites)
  for (t in seq_len(model$n_seasons)) {
    yb   <- model$y[, , t]
    band <- band + lgamma(rowSums(yb) + 1) - rowSums(lgamma(yb + 1)) +
            rowSums(yb * log(pmax(pid, 1e-300)))
  }
  llr <- as.numeric(ev$log_lik_site) + band
  if (mix %in% c("zip", "zinb")) {
    om <- up$zi_omega; log1m <- log1p(-om)
    az <- rowSums(apply(model$y, c(1L, 2L), sum)) == 0
    ll <- log1m + llr
    a  <- log1m + llr[az]; b <- log(om); mx <- pmax(a, b)
    ll[az] <- mx + log(exp(a - mx) + exp(b - mx))
    return(ll)
  }
  llr
}

# Posterior replicate distance-bin arrays: draw a coefficient vector, then per
# site draw the open-population N sequence and the per-period distance obs.
.tobs_simulate_distsamp_open <- function(object, nsim = 1) {
  model <- object$model
  dyn   <- model$dynamics %||% "constant"
  draw_one <- function() {
    idx <- sample.int(nrow(object$draws), 1L)
    if (!identical(dyn, "constant")) {
      fv <- .tobs_fitted_distsamp_open(
        utils::modifyList(object, list(means = object$draws[idx, ])))
      return(.dso_draw_dyn(fv, dyn, model$cutpoints, model$transect,
                           model$n_seasons))
    }
    up  <- .dso_unpack(object$draws[idx, ], model)
    .dso_draw(up$lambda, up$sigma, up$omega, up$gamma, model$cutpoints,
              model$transect, model$n_seasons, r = up$r %||% NULL,
              zi = up$zi_omega %||% 0)
  }
  if (nsim == 1L) return(draw_one())
  lapply(seq_len(nsim), function(s) draw_one())
}

# Core draw for the alternative dynamics: per site the N sequence under the
# requested transition (unmarked tp2..tp5), then the per-period distance obs. `fv`
# is the per-site parameter list from .tobs_fitted_distsamp_open (lambda, sigma,
# pdist via cpd, plus omega / K / gamma / r as the dynamics uses).
.dso_draw_dyn <- function(fv, dynamics, cutpoints, transect, T) {
  lambda <- fv$lambda; sigma <- fv$sigma
  n   <- length(lambda)
  cpd <- .gdr_dist_cp(sigma, cutpoints, transect); pdist <- rowSums(cpd)
  Jb  <- ncol(cpd)
  om  <- fv$omega; gam <- fv$gamma; Kc <- fv$K; rr <- fv$r
  N <- matrix(0L, n, T); N[, 1L] <- stats::rpois(n, lambda)
  for (t in 2:T) {
    prev <- N[, t - 1L]
    N[, t] <- switch(dynamics,
      notrend  = stats::rbinom(n, prev, om) + stats::rpois(n, (1 - om) * lambda),
      autoreg  = stats::rbinom(n, prev, om) + stats::rpois(n, gam * prev),
      trend    = stats::rpois(n, prev * gam),
      ricker   = stats::rpois(n, prev * exp(rr * (1 - prev / Kc))),
      gompertz = stats::rpois(n, prev * exp(rr * (1 - log(prev + 1) / log(Kc + 1)))))
  }
  y <- array(0L, c(n, Jb, T))
  for (i in seq_len(n)) for (t in seq_len(T)) {
    det <- stats::rbinom(1L, N[i, t], pdist[i])
    if (det > 0L)
      y[i, , t] <- as.integer(stats::rmultinom(1L, det, cpd[i, ] / pdist[i]))
  }
  y
}

# Core draw: per site the (Dail-Madsen) N sequence + per-period distance obs. The
# initial abundance N_1 is Poisson (r = NULL) or negative-binomial (size r); a
# structural-zero fraction `zi` sets N_. = 0 for the whole site (never occupied).
.dso_draw <- function(lambda, sigma, omega, gamma, cutpoints, transect, T,
                      r = NULL, zi = 0) {
  n   <- length(lambda)
  cpd <- .gdr_dist_cp(sigma, cutpoints, transect); pdist <- rowSums(cpd)
  Jb  <- ncol(cpd)
  N1  <- if (is.null(r)) stats::rpois(n, lambda)
         else stats::rnbinom(n, size = r, mu = lambda)
  if (zi > 0) N1[stats::runif(n) < zi] <- 0L        # structural zeros
  N   <- matrix(0L, n, T); N[, 1L] <- N1
  for (t in 2:T)
    N[, t] <- stats::rbinom(n, N[, t - 1L], omega) + stats::rpois(n, gamma)
  # A structural-zero site stays empty across all periods (the ZI model's meaning).
  if (zi > 0) N[N1 == 0L, ] <- 0L
  y <- array(0L, c(n, Jb, T))
  for (i in seq_len(n)) for (t in seq_len(T)) {
    det <- stats::rbinom(1L, N[i, t], pdist[i])
    if (det > 0L)
      y[i, , t] <- as.integer(stats::rmultinom(1L, det, cpd[i, ] / pdist[i]))
  }
  y
}

# ---------------------------------------------------------------------------
# Simulator for recovery tests
# ---------------------------------------------------------------------------

#' Simulate an open-population distance-sampling data set
#'
#' Draws from the [distsamp_open()] model: an open metapopulation
#' (`N_1 ~ Poisson / NB(lambda)`, `N_t = Binomial(N_{t-1}, omega) +
#' Poisson(gamma)`) observed by distance sampling (half-normal key, scale
#' `sigma`) at each primary period.
#'
#' @param N Number of sites (default 200).
#' @param cutpoints Distance-bin edges `0 = c_0 < ... < c_B`.
#' @param n_seasons Number of primary periods (default 4).
#' @param transect `"line"` (default) or `"point"`.
#' @param n_abund_covs,n_det_covs Number of abundance / distance covariates.
#' @param beta_lambda,beta_sigma Coefficients on the log-abundance and log-scale
#'   arms. Defaults give moderate abundance / detection.
#' @param omega,gamma Apparent survival probability and recruitment rate
#'   (intercept-only defaults 0.7 / 2.5).
#' @param mixture Initial-abundance distribution: `"poisson"` (default),
#'   `"negbin"` / `"zinb"` (draws `N_1` negative-binomial with size `size`), or
#'   `"zip"` / `"zinb"` (a `zi` fraction of sites are structural zeros).
#' @param size Negative-binomial size (overdispersion `r`) when `mixture` draws
#'   an NB initial abundance (default 3).
#' @param zi Structural-zero probability when `mixture` is `"zip"` / `"zinb"`
#'   (default 0.3).
#' @param dynamics Population-dynamics form (see [distsamp_open()]): `"constant"`
#'   (default), `"notrend"`, `"trend"`, `"autoreg"`, `"ricker"`, or `"gompertz"`.
#'   The alternative dynamics draw the abundance sequence from the matching
#'   transition (Poisson initial abundance only).
#' @param K Carrying capacity for `"ricker"` / `"gompertz"` (default `4 * lambda0`);
#'   `r` the intrinsic growth rate for those forms (default 0.3). For `"trend"` /
#'   `"autoreg"` the recruitment multiplier is `gamma`.
#' @param r Intrinsic growth rate for `"ricker"` / `"gompertz"` (default 0.3).
#' @param seed Optional random seed.
#' @return A list with `y` (`[n_sites x n_bins x n_seasons]` distance-band
#'   counts), `data`, and `truth`.
#' @export
simulate_distsamp_open <- function(N = 200, cutpoints = c(0, 10, 20, 30, 40),
                                   n_seasons = 4L, transect = "line",
                                   n_abund_covs = 1, n_det_covs = 1,
                                   beta_lambda = NULL, beta_sigma = NULL,
                                   omega = 0.7, gamma = 2.5,
                                   mixture = c("poisson", "negbin", "zip", "zinb"),
                                   size = 3, zi = 0.3,
                                   dynamics = c("constant", "notrend", "trend",
                                                "autoreg", "ricker", "gompertz"),
                                   K = NULL, r = 0.3, seed = NULL) {
  mixture  <- match.arg(mixture)
  dynamics <- match.arg(dynamics)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(beta_lambda))
    beta_lambda <- c(log(15), stats::runif(n_abund_covs, -0.3, 0.3))
  if (is.null(beta_sigma))
    beta_sigma <- c(log(stats::median(cutpoints[-1])),
                    stats::runif(n_det_covs, -0.2, 0.2))

  mk <- function(k, tag) {
    d <- data.frame(matrix(stats::rnorm(N * k), N, k))
    names(d) <- paste0(tag, seq_len(k)); d
  }
  ac <- mk(n_abund_covs, "abund_cov"); dc <- mk(n_det_covs, "det_cov")
  data <- cbind(ac, dc)
  lambda <- exp(as.vector(stats::model.matrix(~ ., ac) %*% beta_lambda))
  sigma  <- exp(as.vector(stats::model.matrix(~ ., dc) %*% beta_sigma))

  if (!identical(dynamics, "constant")) {
    if (is.null(K)) K <- 4 * exp(beta_lambda[1])
    fv <- list(lambda = lambda, sigma = sigma, omega = rep(omega, N),
               gamma = rep(gamma, N), K = rep(K, N), r = rep(r, N))
    y <- .dso_draw_dyn(fv, dynamics, as.numeric(cutpoints), transect,
                       as.integer(n_seasons))
    dimnames(y) <- list(NULL, paste0("band", seq_len(dim(y)[2])),
                        paste0("period", seq_len(n_seasons)))
    return(list(y = y, data = data,
      truth = list(beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                   omega = omega, gamma = gamma, K = K, r = r, lambda = lambda,
                   sigma = sigma, mixture = "poisson", dynamics = dynamics)))
  }

  use_nb <- mixture %in% c("negbin", "zinb")
  is_zi  <- mixture %in% c("zip", "zinb")
  y <- .dso_draw(lambda, sigma, rep(omega, N), rep(gamma, N),
                 as.numeric(cutpoints), transect, as.integer(n_seasons),
                 r = if (use_nb) size else NULL, zi = if (is_zi) zi else 0)
  dimnames(y) <- list(NULL, paste0("band", seq_len(dim(y)[2])),
                      paste0("period", seq_len(n_seasons)))
  list(y = y, data = data,
       truth = list(beta_lambda = beta_lambda, beta_sigma = beta_sigma,
                    omega = omega, gamma = gamma, lambda = lambda, sigma = sigma,
                    mixture = mixture, dynamics = dynamics,
                    r = if (use_nb) size else NULL, zi = if (is_zi) zi else NULL))
}
