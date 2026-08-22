# ms_count_pg_gibbs.R - Polya-Gamma Gibbs for the community Bernoulli / binomial
# GLMM (jsdm() and count(response="binomial") community, ms_count(); spOccupancy
# msPGOcc-family / svcMsPGBinom). The community relative-abundance GLMM on a
# LOGISTIC response has no latent state and no detection sub-model -- y is the
# observed k-of-n (n = 1 for jsdm's presence/absence) -- so this is the simplest
# of the PG engines: per-species coefficients with Gaussian community
# hyperpriors, each an exactly conjugate Gaussian update conditional on the
# Polya-Gamma auxiliaries. Only the logistic responses (bernoulli / binomial)
# admit PG augmentation; the Poisson / negbin / gaussian responses of ms_count()
# are not routed here.
#
#   logit(p_{s,i}) = X_i . beta_s,   beta_s ~ N(mu, diag(tau^2))
#   y_{s,i} ~ Binomial(n_{s,i}, p_{s,i})           (n = 1 = jsdm bernoulli)
#
# Per sweep: per species omega ~ PG(n, eta) and the conjugate beta_s, then the
# conjugate community mean + near-Jeffreys Inverse-Gamma community variance. PG
# draws via tulpa's Polson-Scott-Windle sampler. Gives a calibrated
# community-variance posterior (the Laplace-EM attenuates it). v1: non-spatial,
# no latent factors (the sfMsPGBinom / lfJSDM spatial-factor PG variants are
# follow-ups).

.tobs_fit_ms_count_pg_gibbs <- function(model, priors = NULL, sigma.beta = NULL,
                                        n.iter = NULL, n.warmup = NULL,
                                        n.chains = NULL, n.thin = NULL, seed = NULL,
                                        verbose = FALSE, ...) {
  # Sampler defaults come from the one engine table.
  .tobs_fill_sampler(environment(), "pg_gibbs")

  if (!(model$response %in% c("bernoulli", "binomial"))) {
    stop("method = \"pg_gibbs\" applies to the logistic responses only ",
         "(jsdm() / count(response = \"binomial\")); for Poisson / negbin / ",
         "gaussian use method = \"laplace\" or \"nuts\".", call. = FALSE)
  }
  rpg <- get("cpp_rpg", envir = asNamespace("tulpa"))
  X   <- model$X; y <- model$y; nt <- model$n_trials; valid <- model$valid
  S   <- model$n_species; p <- ncol(X)

  # Per-species valid-row index + k (successes) / n (trials) on those rows.
  vi <- lapply(seq_len(S), function(s) which(valid[, s]))
  ks <- lapply(seq_len(S), function(s) as.numeric(y[vi[[s]], s]))
  # Bernoulli (jsdm) has no trials matrix -> one trial per observation.
  ns <- lapply(seq_len(S), function(s)
    if (is.null(nt)) rep(1, length(vi[[s]])) else as.numeric(nt[vi[[s]], s]))
  Xs <- lapply(seq_len(S), function(s) X[vi[[s]], , drop = FALSE])

  n_keep <- length(seq.int(n.warmup + 1L, n.iter, by = n.thin))
  cn <- model$process_info[[1L]]$coef_names
  run_chain <- function(chain_id) {
    set.seed(seed + chain_id)
    b  <- matrix(stats::rnorm(S * p, 0, 0.2), S, p)
    mu <- rep(0, p); tau2 <- rep(1, p)
    mu_draws <- matrix(NA_real_, n_keep, p); tau_draws <- matrix(NA_real_, n_keep, p)
    b_sum <- matrix(0, S, p); nsum <- 0L; ki <- 0L
    for (it in seq_len(n.iter)) {
      for (s in seq_len(S)) {
        if (!length(vi[[s]])) next
        eta <- as.vector(Xs[[s]] %*% b[s, ])
        om  <- rpg(ns[[s]], eta)
        b[s, ] <- .tobs_pg_draw_beta(Xs[[s]], om, ks[[s]] - ns[[s]] / 2,
                                     1 / tau2, mu / tau2)
      }
      cu <- .tobs_pg_community_update(b, mu, tau2, S); mu <- cu$mu; tau2 <- cu$tau2
      if (it > n.warmup && ((it - n.warmup - 1L) %% n.thin == 0L)) {
        ki <- ki + 1L; mu_draws[ki, ] <- mu; tau_draws[ki, ] <- sqrt(tau2)
        b_sum <- b_sum + b; nsum <- nsum + 1L
      }
    }
    list(mu = mu_draws, tau = tau_draws, b = b_sum / nsum)
  }

  chains <- lapply(seq_len(n.chains), run_chain)
  summ <- .tobs_pg_summarize(lapply(chains, `[[`, "mu"), cn)

  tau_all <- do.call(rbind, lapply(chains, `[[`, "tau"))
  sd_mu <- apply(tau_all, 2L, stats::median); names(sd_mu) <- cn

  coef_mu <- Reduce(`+`, lapply(chains, `[[`, "b")) / n.chains
  rownames(coef_mu) <- model$species_names; colnames(coef_mu) <- cn
  blup_mu <- sweep(coef_mu, 2L, summ$means, "-")

  .tobs_pg_finalize_fit(
    summ, cn, model, model$process_info, N = sum(model$valid),
    n.iter = n.iter, n.chains = n.chains,
    extra = list(ms_community = list(
      Sigma_mu = diag(sd_mu^2, p), sd_mu = sd_mu,
      coef_mu = coef_mu, blup_mu = blup_mu)))
}
