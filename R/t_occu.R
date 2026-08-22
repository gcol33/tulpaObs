# t_occu.R - multi-season occupancy with a temporal (AR1) year random effect on
# the occupancy state (spOccupancy tPGOcc). Distinct from dyn_occu() (colext):
# there is NO colonization / extinction transition -- the per-(site, season)
# occupancy state is a Bernoulli GLMM with a shared AR1 year effect on the
# logit,
#
#   z_{i,t}          ~ Bernoulli(psi_{i,t}),  logit psi_{i,t} = X_occ_i . beta + eta_t
#   eta_t            = rho eta_{t-1} + w_t,    w_t ~ N(0, sigma^2)   (AR1 year effect)
#   y_{i,t,j}|z=1    ~ Bernoulli(p_{i,t}),    logit p_{i,t}   = X_det_i . beta_p
#
# Conditional on the year effects eta the seasons factorise, so a Polya-Gamma
# Gibbs sampler is exact: draw z, then the joint (beta_occ, eta) as ONE Gaussian
# Markov random field update (the AR1 precision on eta is the GMRF prior, exactly
# as spPGOcc's ICAR precision is for a spatial field), then beta_p, then the AR1
# hyperparameters (sigma^2 conjugate Inverse-Gamma, rho on a grid). The year
# effect is centred sum-to-zero each sweep with its level moved to the intercept.
# PG draws use tulpa's Polson-Scott-Windle sampler (tulpa:::cpp_rpg).
#
# v1: site-level occupancy + detection covariates (broadcast across seasons),
# method = "pg_gibbs" only (this IS spOccupancy's engine for the family).

#' Multi-season occupancy with a temporal (AR1) year effect (tPGOcc)
#'
#' The spOccupancy `tPGOcc` model: per-`(site, season)` Bernoulli occupancy with a
#' shared AR1 year random effect on the occupancy logit and NO colonization /
#' extinction dynamics (that is [dyn_occu()]). Use it for an occupancy trend over
#' years with a temporal random effect. Fit with `method = "pg_gibbs"` (the
#' Polya-Gamma Gibbs sampler spOccupancy uses); `y` is a 3D array
#' `[n_sites x n_seasons x max_visits]` (or a list of per-season
#' `[n_sites x max_visits]` matrices), the occupancy `formula` and `detection`
#' formula are site-level.
#'
#' @return A `tobs_family` object.
#' @export
t_occu <- function() {
  obs_family(
    name           = "t_occu",
    class_long     = "multi-season occupancy (AR1 year effect)",
    latent         = "bernoulli_panel",
    observation    = "binomial_detection",
    replicates     = "required",
    default_engine = "pg_gibbs",
    status         = "working"
  )
}

.tobs_build_t_occu <- function(occ_formula, det_formula, data, y) {
  if (is.list(y) && !is.array(y)) {
    n_seasons <- length(y)
    n_sites <- nrow(y[[1L]]); max_visits <- ncol(y[[1L]])
    ya <- array(NA_integer_, dim = c(n_sites, n_seasons, max_visits))
    for (t in seq_len(n_seasons)) ya[, t, ] <- as.matrix(y[[t]])
    y <- ya
  }
  if (length(dim(y)) != 3L)
    stop("t_occu() y must be a 3D array [n_sites x n_seasons x max_visits] or a ",
         "list of per-season matrices.", call. = FALSE)
  n_sites <- dim(y)[1L]; n_seasons <- dim(y)[2L]; max_visits <- dim(y)[3L]
  .tobs_check_site_count(n_sites, nrow(data), "sites")
  if (n_seasons < 2L)
    stop("t_occu() needs >= 2 seasons; for one season use occu().", call. = FALSE)

  X_occ <- stats::model.matrix(occ_formula, data)
  X_det <- stats::model.matrix(det_formula, data)

  # Per-(site, season) detection sufficient statistics.
  nvis <- kdet <- matrix(0L, n_sites, n_seasons)
  for (i in seq_len(n_sites)) for (t in seq_len(n_seasons)) {
    raw <- y[i, t, ]; raw <- raw[!is.na(raw) & raw >= 0]
    nvis[i, t] <- length(raw); kdet[i, t] <- sum(raw == 1L)
  }
  anydet <- kdet > 0L

  structure(list(
    model_type  = "t_occu",
    y           = y,
    n_sites     = n_sites, n_seasons = n_seasons, max_visits = max_visits,
    X_occ       = X_occ, X_det = X_det,
    nvis        = nvis, kdet = kdet, anydet = anydet,
    formulas    = list(occ = occ_formula, det = det_formula),
    data        = data,
    process_info = list(
      list(name = "psi", p = ncol(X_occ), coef_names = colnames(X_occ), link = "logit"),
      list(name = "p",   p = ncol(X_det), coef_names = colnames(X_det), link = "logit"))
  ), class = "tobs_model")
}

# AR1 precision matrix on the T year effects (unit innovation variance): the
# standard stationary AR1 GMRF precision, tridiagonal.
.t_occu_ar1_Q <- function(T_s, rho) {
  Q <- matrix(0, T_s, T_s)
  d <- c(1, rep(1 + rho^2, T_s - 2L), 1)
  diag(Q) <- d
  for (t in seq_len(T_s - 1L)) { Q[t, t + 1L] <- -rho; Q[t + 1L, t] <- -rho }
  Q
}

.tobs_fit_t_occu_pg_gibbs <- function(model, priors = NULL, sigma.beta = NULL,
                                      n.iter = NULL, n.warmup = NULL,
                                      n.chains = NULL, n.thin = NULL, seed = NULL,
                                      verbose = FALSE, ...) {
  # Sampler defaults come from the one engine table.
  .tobs_fill_sampler(environment(), "pg_gibbs")

  rpg   <- get("cpp_rpg", envir = asNamespace("tulpa"))
  X_occ <- model$X_occ; X_det <- model$X_det
  n <- model$n_sites; T_s <- model$n_seasons
  p_psi <- ncol(X_occ); p_p <- ncol(X_det)
  nvis <- model$nvis; kdet <- model$kdet; anydet <- model$anydet
  B0inv_psi <- diag(1 / sigma.beta^2, p_psi)        # joins the (beta_psi, eta) GMRF
  prec      <- 1 / sigma.beta^2                      # detection-arm coefficient prior
  ig_a <- 0.1; ig_b <- 0.1                          # near-Jeffreys IG(sigma^2)
  rho_grid <- seq(-0.95, 0.95, by = 0.05)           # AR1 correlation grid

  # Long-form occupancy design over (site, season), site-major season-minor, plus
  # the season index per row (for the year-effect indicator Z).
  Xf  <- X_occ[rep(seq_len(n), each = T_s), , drop = FALSE]   # [(n*T) x p_psi]
  srow <- rep(seq_len(T_s), times = n)                        # season of each row
  N_rows <- n * T_s

  n_keep <- length(seq.int(n.warmup + 1L, n.iter, by = n.thin))
  par_names <- c(paste0("psi_", model$process_info[[1L]]$coef_names),
                 paste0("p_",   model$process_info[[2L]]$coef_names),
                 "log_sigma_ar1", "rho_ar1")

  run_chain <- function(chain_id) {
    set.seed(seed + chain_id)
    bpsi <- stats::rnorm(p_psi, 0, 0.2); bp <- stats::rnorm(p_p, 0, 0.2)
    eta <- rep(0, T_s); sigma2 <- 1; rho <- 0.3
    out <- matrix(NA_real_, n_keep, p_psi + p_p + 2L); ki <- 0L
    eta_sum <- numeric(T_s); nsum <- 0L
    for (it in seq_len(n.iter)) {
      psi_row <- stats::plogis(as.vector(Xf %*% bpsi) + eta[srow])
      p_it    <- stats::plogis(as.vector(X_det %*% bp))         # per-site detection
      # 1. latent occupancy z_{i,t}
      psi_m <- matrix(psi_row, n, T_s, byrow = TRUE)
      pmat  <- matrix(p_it, n, T_s)
      z <- .tobs_pg_draw_z(psi_m, (1 - pmat)^nvis, anydet)
      zf <- as.vector(t(z))                                    # row order (i,t)
      # 2. joint (beta_psi, eta) GMRF update given omega_psi, AR1 precision on eta
      eta_lin <- as.vector(Xf %*% bpsi) + eta[srow]
      om <- rpg(rep(1, N_rows), eta_lin); kap <- zf - 0.5
      Q <- .t_occu_ar1_Q(T_s, rho) / sigma2
      Pbb <- crossprod(Xf, Xf * om) + B0inv_psi                # p x p
      # Z' Omega Z is diagonal (per-season omega sum); Z' Omega Xf via rowsum.
      om_season <- as.numeric(tapply(om, srow, sum))
      Pee <- diag(om_season, T_s) + Q                          # T x T
      Pbe <- t(rowsum(Xf * om, srow))                          # p x T (X' Omega Z)
      Prec <- rbind(cbind(Pbb, Pbe), cbind(t(Pbe), Pee))
      rhs  <- c(crossprod(Xf, kap), as.numeric(tapply(kap, srow, sum)))
      L <- chol(Prec + diag(1e-8, p_psi + T_s))
      mean_v <- backsolve(L, forwardsolve(t(L), rhs))
      samp <- mean_v + backsolve(L, stats::rnorm(p_psi + T_s))
      bpsi <- samp[seq_len(p_psi)]; eta_raw <- samp[p_psi + seq_len(T_s)]
      # 3. AR1 hyperparameters on the RAW year effect (its actual AR1 structure).
      # sigma^2 conjugate Inverse-Gamma; rho on a grid over the AR1 log-density.
      ss <- (1 - rho^2) * eta_raw[1L]^2 +
            sum((eta_raw[-1L] - rho * eta_raw[-T_s])^2)
      sigma2 <- 1 / stats::rgamma(1, ig_a + T_s / 2, ig_b + 0.5 * ss)
      lp_rho <- vapply(rho_grid, function(rg) {
        s <- (1 - rg^2) * eta_raw[1L]^2 +
             sum((eta_raw[-1L] - rg * eta_raw[-T_s])^2)
        0.5 * log(1 - rg^2) - 0.5 * s / sigma2
      }, numeric(1))
      w <- exp(lp_rho - max(lp_rho)); rho <- sample(rho_grid, 1L, prob = w)
      # Now sum-to-zero the year effect for reporting / the next psi (its level is
      # confounded with the occupancy intercept), moving the mean to the intercept.
      # The AR1 hyperparameters above used the raw series, so rho is not whitened.
      me <- mean(eta_raw); eta <- eta_raw - me; bpsi[1L] <- bpsi[1L] + me
      # 4. detection coefficients at occupied (site, season) with visits
      occ <- which(z == 1L & nvis > 0L, arr.ind = TRUE)
      if (nrow(occ) >= p_p) {
        si <- occ[, 1L]; ti <- occ[, 2L]
        Xo <- X_det[si, , drop = FALSE]
        nv_o <- nvis[cbind(si, ti)]; kd_o <- kdet[cbind(si, ti)]
        om_p <- rpg(nv_o, as.vector(Xo %*% bp))
        bp <- .tobs_pg_draw_beta(Xo, om_p, kd_o - nv_o / 2, prec)
      }
      if (it > n.warmup && ((it - n.warmup - 1L) %% n.thin == 0L)) {
        ki <- ki + 1L
        out[ki, ] <- c(bpsi, bp, 0.5 * log(sigma2), rho)
        eta_sum <- eta_sum + eta; nsum <- nsum + 1L
      }
    }
    list(draws = out, eta = eta_sum / nsum)
  }

  chains <- lapply(seq_len(n.chains), run_chain)
  summ <- .tobs_pg_summarize(lapply(chains, `[[`, "draws"), par_names)
  eta_mean <- Reduce(`+`, lapply(chains, `[[`, "eta")) / n.chains

  .tobs_pg_finalize_fit(
    summ, par_names, model, model$process_info, N = sum(model$nvis),
    n.iter = n.iter, n.chains = n.chains,
    extra = list(temporal_field = eta_mean))
}

.dispatch_t_occu <- function(formula, data, family, detection, y, visits,
                             engine, priors, control,
                             approx = "gaussian_laplace",
                             correction = "none", ...) {
  if (is.null(detection))
    stop("t_occu() requires a `detection` formula.", call. = FALSE)
  if (is.null(y))
    stop("t_occu() requires `y` (a 3D array [n_sites x n_seasons x max_visits] ",
         "or a list of per-season matrices).", call. = FALSE)
  if (!identical(engine, "pg_gibbs"))
    stop("t_occu() supports method = \"pg_gibbs\" only (the Polya-Gamma Gibbs ",
         "sampler spOccupancy's tPGOcc uses).", call. = FALSE)
  model <- .tobs_build_t_occu(occ_formula = formula, det_formula = detection,
                              data = data, y = y)
  control <- .tobs_control_defaults(control, "pg_gibbs", "t_occu")
  .tobs_fit_t_occu_pg_gibbs(
    model, priors = priors,
    sigma.beta = control[["sigma.beta"]],
    n.iter   = as.integer(control[["n.iter"]]),
    n.warmup = as.integer(control[["n.warmup"]]),
    n.chains = max(as.integer(control[["n.chains"]]), 2L),
    n.thin   = as.integer(control[["n.thin"]]),
    seed     = as.integer(control[["seed"]]),
    verbose  = isTRUE(control[["verbose"]]))
}

#' Simulate multi-season occupancy with an AR1 year effect (tPGOcc)
#'
#' @param N Number of sites (default 150).
#' @param T_seasons Number of seasons / years (default 8).
#' @param J Visits per season (default 3).
#' @param beta_occ Occupancy coefficients (intercept first; default `c(0.2)`).
#' @param p Detection probability (default 0.4).
#' @param rho AR1 correlation of the year effect (default 0.6).
#' @param sigma AR1 innovation SD of the year effect (default 0.7).
#' @param seed Optional random seed.
#' @return A list with `y` (3D array `[N x T x J]`), `data`, and `truth`.
#' @seealso [t_occu()], the family this simulates for.
#' @export
simulate_t_occu <- function(N = 150, T_seasons = 8, J = 3, beta_occ = c(0.2),
                            p = 0.4, rho = 0.6, sigma = 0.7, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  data <- data.frame(x = stats::rnorm(N))
  X_occ <- if (length(beta_occ) > 1L) stats::model.matrix(~ x, data)
           else stats::model.matrix(~ 1, data)
  # AR1 year effect (stationary start).
  eta <- numeric(T_seasons)
  eta[1L] <- stats::rnorm(1, 0, sigma / sqrt(1 - rho^2))
  for (t in 2:T_seasons) eta[t] <- rho * eta[t - 1L] + stats::rnorm(1, 0, sigma)
  eta <- eta - mean(eta)
  lin <- as.vector(X_occ %*% beta_occ)
  z <- matrix(0L, N, T_seasons)
  for (t in seq_len(T_seasons))
    z[, t] <- stats::rbinom(N, 1L, stats::plogis(lin + eta[t]))
  y <- array(0L, dim = c(N, T_seasons, J))
  for (t in seq_len(T_seasons)) for (j in seq_len(J))
    y[, t, j] <- ifelse(z[, t] == 1L, stats::rbinom(N, 1L, p), 0L)
  list(y = y, data = data,
       truth = list(beta_occ = beta_occ, p = p, rho = rho, sigma = sigma,
                    eta = eta, z = z))
}
