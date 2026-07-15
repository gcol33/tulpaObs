# =============================================================================
# ms_count_spatial.R - community relative-abundance GLMM with a shared latent
# structure: a shared areal field (the spAbundance sfMsAbund / svcMsAbund
# analogues), latent factors (lfMsAbund), or BOTH (the spatial-factor case).
# gcol33/tulpaObs#117, #118. Poisson.
#
#   log mu_{s,i} = X_i . (mu_beta + b_s) + sum_k W[i,k] F[u(i),k]
#                                        + sum_q lambda_{s,q} zeta_{q,i}
#   b_s ~ N(0, Sigma_beta),  F ~ ICAR/CAR/BYM2(tau),  zeta_{q,i} ~ N(0, 1)
#
# The block coordinate ascent, the areal Newton, the factor update, and the field
# hyperparameter grids live in R/community_latent.R and are shared with every
# other community family. This file supplies only the Poisson working oracle and
# the model wiring.
# =============================================================================


# Community GLMM working oracle: the score and curvature of the per-(site,
# species) log-likelihood with respect to an additive offset on the linear
# predictor. Poisson (log link) gives (y - mu, mu); Bernoulli (logit link, the
# jsdm() response) gives (y - psi, psi (1 - psi)). `data_ll` may drop terms that
# are constant in eta (the Poisson lgamma(y + 1) normaliser) -- it is only ever
# compared across a field hyperparameter grid, over which they cancel.
.tobs_ms_count_oracle <- function(y_mat, link = "log") {
  dims <- list(n_sites = nrow(y_mat), n_species = ncol(y_mat))
  if (identical(link, "logit")) {
    return(c(dims, list(
      working = function(eta) {
        psi <- stats::plogis(eta)
        list(score = y_mat - psi, curv = psi * (1 - psi))
      },
      data_ll = function(eta) {
        sum(ifelse(y_mat > 0, stats::plogis(eta,  log.p = TRUE),
                              stats::plogis(-eta, log.p = TRUE)))
      })))
  }
  c(dims, list(
    working = function(eta) {
      mu <- exp(pmin(eta, 700))
      list(score = y_mat - mu, curv = mu)
    },
    data_ll = function(eta) {
      e <- pmin(eta, 700)
      sum(y_mat * e - exp(e))
    }))
}


# Fit a community count model with a shared areal field (`spatial`), latent
# factors (`latent`), or both. Single source of truth for every latent-count
# route; the field-only / factor-only models are the special cases.
.tobs_fit_ms_count_latent <- function(model, spatial = NULL, latent = NULL,
                                      max.iter = 200L, tol = 1e-4,
                                      sigma.beta = 5, priors = NULL,
                                      max.outer = 25L, verbose = FALSE, ...) {
  response <- model$response %||% "poisson"
  # Poisson (ms_count) and Bernoulli (jsdm) carry no dispersion parameter, so the
  # per-site latent structure is identified against them. A negbin size /
  # Gaussian residual variance would be confounded with a per-site field, exactly
  # as in the single-species areal count gate.
  if (!response %in% c("poisson", "bernoulli")) {
    stop(sprintf(paste0(
      "A shared field / latent() factors on a community GLMM support the ",
      "Poisson (ms_count) and Bernoulli (jsdm) responses; got \"%s\". The ",
      "negbin size / Gaussian residual variance is not identified against a ",
      "per-site latent structure (gcol33/tulpaObs#117)."), response),
      call. = FALSE)
  }
  what <- if (identical(response, "bernoulli")) "jsdm()" else "ms_count()"
  X <- model$X; P <- ncol(X); S <- model$n_species; Ns <- model$n_sites
  su <- model$summaries
  if (any(!model$valid)) {
    stop(what, " shared field / latent factors need a complete y (no NA ",
         "species-site cells) (gcol33/tulpaObs#117).", call. = FALSE)
  }
  y_mat  <- matrix(as.numeric(model$y), Ns, S)
  link   <- model$link %||% "log"
  oracle <- .tobs_ms_count_oracle(y_mat, link = link)
  is_bern <- identical(response, "bernoulli")

  arm_idx <- list(mu = seq_len(P))
  mbar <- mean(y_mat)
  mu0 <- numeric(P)
  mu0[1L] <- if (is_bern) stats::qlogis(min(max(mbar, 1e-3), 1 - 1e-3))
             else log(max(mean(rowSums(y_mat) / S), 0.1))

  em_fit <- function(site_off, fac_off, em_prev) {
    sp_ll <- function(s, theta, global) {
      e <- as.numeric(su[[s]]$X %*% theta) + site_off + fac_off[, s]
      if (is_bern) {
        sum(ifelse(su[[s]]$y > 0, stats::plogis(e,  log.p = TRUE),
                                  stats::plogis(-e, log.p = TRUE)))
      } else {
        sum(stats::dpois(su[[s]]$y, exp(pmin(e, 700)), log = TRUE))
      }
    }
    sp_grad <- function(s, theta, global) {
      e <- as.numeric(su[[s]]$X %*% theta) + site_off + fac_off[, s]
      mu_s <- if (is_bern) stats::plogis(e) else exp(pmin(e, 700))
      as.numeric(crossprod(su[[s]]$X, su[[s]]$y - mu_s))
    }
    .tobs_community_em(
      S = S, P = P, arm_idx = arm_idx, sp_ll = sp_ll, sp_grad = sp_grad,
      init_mu = if (is.null(em_prev)) mu0 else em_prev$mu,
      init_global = numeric(0), penalize_global = FALSE,
      sigma_beta = sigma.beta, priors = priors, sigma_init = 0.3,
      max_iter = min(as.integer(max.iter), 60L), tol = as.numeric(tol),
      newton_max = 30L, verbose = FALSE)
  }
  offset_of <- function(em) {
    vapply(seq_len(S), function(s) as.numeric(X %*% (em$mu + em$b_list[[s]])),
           numeric(Ns))
  }

  res <- .tobs_community_latent_ascent(
    spatial = spatial, latent = latent, model = model, what = what,
    make_oracle = function(em) oracle, em_fit = em_fit, offset_of = offset_of,
    allow = c("icar", "car_proper", "bym2", "spde"),
    tol = tol, max.outer = max.outer, verbose = verbose)

  fit <- build_ms_count_fit(model, res$em, arm_idx, disp = NULL)
  fit$method <- "laplace"
  fit <- .tobs_latent_attach_field(fit, res, spatial, "count_field_offset")
  fit <- .tobs_latent_attach_factor(fit, res, latent, model,
                                    "count_factor_offset")
  fit
}
