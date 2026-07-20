# =============================================================================
# ms_abun_latent.R - community N-mixture with a shared latent structure: latent
# factors (the spAbundance lfMsNMix analogue), a shared field, or both (the
# spatial-factor case). gcol33/tulpaObs#117. Poisson.
#
#   N_{s,i}       ~ Poisson(lambda_{s,i})
#   y_{s,i,j} | N ~ Binomial(N_{s,i}, p_{s,i,j})
#   log lambda_{s,i} = X_i . (mu_lambda + b_lambda_s) + sum_k W[i,k] F[u(i),k]
#                                                     + sum_q lambda_{s,q} zeta_{q,i}
#   logit p_{s,i,j}  = X_p_{ij} . (mu_p + b_p_s)
#
# The field F is SHARED across species (it enters every species' predictor with
# loading 1); the per-species structure is carried by the factors zeta ~ N(0, I)
# through their per-species loadings. With both present the loadings are centred
# across species, so the field owns the shared spatial mean and the factors own
# the between-species residual co-occurrence. This is the same decomposition
# ms_count() / ms_occu() / jsdm() use, and it differs from spAbundance's sfMsNMix,
# which instead gives each FACTOR its own spatial process; the shared-field +
# iid-factor split here is the house convention across every community family.
#
# The latent N_{s,i} still integrates out in closed form per species-site, so the
# whole latent structure sits on eta_lambda and the family reduces to ONE
# callback for the shared driver (R/community_latent.R):
#
#   working(eta_lambda) -> (score, curv)
#
# the score and curvature of the Royle marginal with respect to an additive
# offset on the abundance log-predictor. R/nmix_site_marginal.R already
# differentiates that marginal through both arms, so
#
#   score = grad_eta_lambda
#   curv  = info_eta_lambda - var_N * score_wt_lambda^2
#
# is the (1,1) entry of the per-site marginal observed-information block
# B_i = diag(I_lambda, I_p) - Var(N_i|y) v v', v = (-w_i, p_i1..p_iJ) (Louis
# 1982) -- the abundance curvature with the detection arm profiled out. The block
# coordinate ascent, the field Newton, the factor update and the hyperparameter
# grids are all the driver's; this file supplies the oracle and the wiring.
#
# The plain shared field with NO factors keeps its dedicated C++ path
# (nmix_community_laplace_*, tulpaObs#12) -- that route is faster and already
# recovery-tested, and this driver does not replace it.
# =============================================================================


# Per-species Royle marginals over the shared abundance design. `lf` is the
# long-form (y, site_idx, species_idx, X_p); each species gets its own marginal
# closure over its own visit rows, sharing the site-level X_lambda.
.tobs_ms_abun_marginals <- function(model, K_max = NULL) {
  lf <- .tobs_ms_nmix_longform(model)
  X_lambda <- model$X_processes[[1L]]
  if (is.null(K_max)) K_max <- max(lf$y) + 100L
  K_max <- as.integer(K_max)
  marg <- lapply(seq_len(model$n_species), function(s) {
    k <- lf$species_idx == s
    nmix_site_marginal(y = lf$y[k], site_idx = lf$site_idx[k],
                       X_lambda = X_lambda, X_p = lf$X_p[k, , drop = FALSE],
                       mixture = "P", K_max = K_max)
  })
  list(lf = lf, marg = marg, K_max = K_max,
       X_p = lapply(seq_len(model$n_species),
                    function(s) lf$X_p[lf$species_idx == s, , drop = FALSE]))
}


# The working oracle at fixed detection predictors. `eta_p_list[[s]]` is species
# s's per-visit logit predictor; the Royle marginal is then a function of that
# species' abundance log-predictor alone, which is what the field / factor
# updates perturb.
.tobs_ms_abun_oracle <- function(marg, eta_p_list, Ns, S) {
  eval_all <- function(eta) {
    lapply(seq_len(S), function(s)
      marg[[s]]$eval(eta[, s], eta_p_list[[s]]))
  }
  list(
    n_sites = Ns, n_species = S,
    working = function(eta) {
      evs   <- eval_all(eta)
      score <- vapply(evs, function(ev) as.numeric(ev$grad_eta_lambda), numeric(Ns))
      curv  <- vapply(evs, function(ev)
        pmax(as.numeric(ev$info_eta_lambda -
                        ev$var_N * ev$score_wt_lambda^2), 1e-8), numeric(Ns))
      list(score = score, curv = curv)
    },
    ll_cell = function(eta) vapply(eval_all(eta),
      function(ev) as.numeric(ev$log_lik_site), numeric(Ns)),
    data_ll = function(eta) sum(vapply(eval_all(eta), function(ev) ev$log_lik, 0)))
}


# Fit a community N-mixture with latent factors (`latent`), a per-species-loaded
# shared field (`spatial`), or both, by the shared block coordinate ascent.
# Poisson only: a negative-binomial size is a second per-site dispersion, which is
# not identified against a per-site latent structure.
.tobs_fit_ms_abun_latent <- function(model, spatial = NULL, latent = NULL,
                                     mixture = "poisson", K_max = NULL,
                                     max.iter = 100L, tol = 1e-4,
                                     sigma.beta = 5, priors = NULL,
                                     max.outer = 25L, verbose = FALSE, ...) {
  if (!identical(mixture, "poisson")) {
    stop("ms_abun() with latent() factors is Poisson-only: a negative-binomial ",
         "size is a second per-site dispersion and is not identified against a ",
         "per-site latent structure (gcol33/tulpaObs#117). Use ",
         "mixture = \"poisson\", or drop latent() for the negbin community fit.",
         call. = FALSE)
  }
  S  <- model$n_species
  Ns <- model$n_sites
  ms <- .tobs_ms_abun_marginals(model, K_max = K_max)
  lf <- ms$lf
  X_lambda <- model$X_processes[[1L]]
  p_lam <- ncol(X_lambda)
  p_p   <- ncol(lf$X_p)
  P     <- p_lam + p_p
  lam_idx <- seq_len(p_lam)
  p_idx   <- p_lam + seq_len(p_p)
  arm_idx <- list(lambda = lam_idx, p = p_idx)

  # Warm-start the community means at the no-latent community fit: the N-mixture
  # marginal is poorly behaved from a naive start (lambda and p trade off), and
  # this is the same C++ EM the plain community path runs.
  warm <- nmix_laplace_re(
    y = lf$y, site_idx = lf$site_idx, species_idx = lf$species_idx,
    X_lambda = X_lambda, X_p = lf$X_p, n_sites = Ns, n_species = S,
    K_max = ms$K_max, max_iter = as.integer(max.iter), mixture = "P",
    optimizer = "em", n_quad = 1L, lkj_eta = 1, verbose = FALSE)
  mu0 <- c(as.numeric(warm$mu_lambda), as.numeric(warm$mu_p))

  em_fit <- function(site_off, fac_off, em_prev) {
    sp_ll <- function(s, theta, global) {
      e_lam <- as.numeric(X_lambda %*% theta[lam_idx]) + site_off + fac_off[, s]
      ms$marg[[s]]$eval(e_lam, as.numeric(ms$X_p[[s]] %*% theta[p_idx]))$log_lik
    }
    sp_grad <- function(s, theta, global) {
      e_lam <- as.numeric(X_lambda %*% theta[lam_idx]) + site_off + fac_off[, s]
      ev <- ms$marg[[s]]$eval(e_lam, as.numeric(ms$X_p[[s]] %*% theta[p_idx]))
      c(as.numeric(crossprod(X_lambda, ev$grad_eta_lambda)),
        as.numeric(crossprod(ms$X_p[[s]], ev$grad_eta_p)))
    }
    # Analytic per-species marginal observed information: the design-sandwiched
    # per-site Louis block
    #   B_i = diag(I_lambda, I_p) - Var(N_i|y) v v',  v = (-w_i, p_i1..p_iJ)
    # which nmix_site_marginal() already forms. Site i's coordinates are
    # (eta_lambda_i, eta_p over its visits), so D_i is the block design row
    # [X_lambda[i, ] | 0 ; 0 | X_p[visits, ]] and the information is
    # sum_i D_i' B_i D_i. This is what the EM's finite-difference Hessian spends
    # 2P marginal sweeps per species per Newton step approximating.
    sp_info <- function(s, theta, global) {
      e_lam <- as.numeric(X_lambda %*% theta[lam_idx]) + site_off + fac_off[, s]
      mg <- ms$marg[[s]]
      ev <- mg$eval(e_lam, as.numeric(ms$X_p[[s]] %*% theta[p_idx]))
      I  <- matrix(0, P, P)
      for (i in seq_len(Ns)) {
        obs <- mg$obs_by_site[[i]]
        D <- matrix(0, 1L + length(obs), P)
        D[1L, lam_idx] <- X_lambda[i, ]
        if (length(obs)) {
          D[-1L, p_idx] <- ms$X_p[[s]][obs, , drop = FALSE]
        }
        I <- I + crossprod(D, mg$obs_info_block(i, ev) %*% D)
      }
      (I + t(I)) / 2
    }
    # Each Royle-marginal evaluation sums over the latent N, so re-running the
    # community EM from a cold start on every outer pass dominates the fit. Warm
    # start the whole state (means, deviations, community covariances) from the
    # previous pass: the latent offset moves little once the ascent settles, so
    # the resumed EM converges in a few iterations.
    .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      sp_info = sp_info,
      init_mu = if (is.null(em_prev)) mu0 else em_prev$mu,
      init_b     = em_prev$b_list,
      init_Sigma = em_prev$Sigma,
      init_global = numeric(0), penalize_global = FALSE,
      sigma_beta = sigma.beta, priors = priors, sigma_init = 0.3,
      max_iter = min(as.integer(max.iter), 40L), tol = as.numeric(tol),
      newton_max = 20L, verbose = FALSE)
  }
  # The structured arm is the abundance arm: the driver perturbs eta_lambda.
  offset_of <- function(em) {
    vapply(seq_len(S), function(s)
      as.numeric(X_lambda %*% (em$mu + em$b_list[[s]])[lam_idx]), numeric(Ns))
  }
  # The detection predictors move with the community fit, so the oracle is
  # rebuilt from the current EM state each outer pass.
  make_oracle <- function(em) {
    eta_p_list <- lapply(seq_len(S), function(s)
      as.numeric(ms$X_p[[s]] %*% (em$mu + em$b_list[[s]])[p_idx]))
    .tobs_ms_abun_oracle(ms$marg, eta_p_list, Ns, S)
  }

  res <- .tobs_community_latent_ascent(
    spatial = spatial, latent = latent, model = model, what = "ms_abun()",
    make_oracle = make_oracle, em_fit = em_fit, offset_of = offset_of,
    allow = c("icar", "car_proper", "bym2", "spde"),
    tol = tol, max.outer = max.outer, verbose = verbose)

  em  <- res$em
  raw <- list(
    mu_lambda = em$mu[lam_idx], mu_p = em$mu[p_idx],
    vcov = em$Vf,
    Sigma_lambda = em$Sigma$lambda, Sigma_p = em$Sigma$p,
    b_lambda = do.call(rbind, lapply(em$b_list, function(b) b[lam_idx])),
    b_p      = do.call(rbind, lapply(em$b_list, function(b) b[p_idx])),
    log_lik = em$logML, converged = isTRUE(em$converged), n_iter = em$n_iter,
    optimizer = "block_coordinate", n_quad = 1L, lkj_eta = 1)
  fit <- build_ms_nmix_fit(raw, model, mixture = "poisson", spatial = NULL)
  fit$method <- "laplace"
  fit <- .tobs_latent_attach_field(fit, res, spatial, "nmix_field_offset")
  fit <- .tobs_latent_attach_factor(fit, res, latent, model,
                                    "nmix_factor_offset")
  fit
}
