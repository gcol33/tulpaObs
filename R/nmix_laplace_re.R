#' Community / multispecies N-mixture by Laplace
#'
#' @description
#' Fits the community (spAbundance `msNMix`) N-mixture model: a per-species
#' Royle (2004) N-mixture with Gaussian community hyperpriors on the per-species
#' abundance and detection coefficients,
#' \deqn{N_{s,i} \sim \mathrm{Poisson}(\lambda_{s,i}), \quad
#'       y_{s,i,j} | N \sim \mathrm{Binomial}(N_{s,i}, p_{s,i,j}),}
#' \deqn{\log \lambda_{s,i} = X_\lambda^{(i)} (\mu_\lambda + b^\lambda_s), \quad
#'       \mathrm{logit}\, p_{s,i,j} = X_p^{(ij)} (\mu_p + b^p_s),}
#' \deqn{b^\lambda_s \sim N(0, \Sigma_\lambda), \quad b^p_s \sim N(0, \Sigma_p).}
#'
#' The latent abundances integrate out per species-site in closed form (the
#' shared N-mixture kernel exposed by [nmix_site_marginal()]); the
#' per-species coefficient deviations \eqn{b_s = (b^\lambda_s, b^p_s)} are the
#' random effects. Both solvers assemble the per-species marginal in compiled
#' code via a native oracle (no per-group round trip into R). At `n_quad = 1`
#' (the default, the joint Laplace / glmer `nAGQ = 1`) the fit is a Laplace-EM
#' (block-coordinate Newton mode + closed-form covariance M-step + Schur-
#' complement SEs) -- the fast production path. A higher `n_quad` routes the same
#' native oracle through the shared compiled AGHQ engine ([tulpa_re_aghq()]),
#' replacing each species' Laplace integral with adaptive Gauss-Hermite
#' quadrature to reduce the small-cluster (few species) downward bias of the
#' community covariances `Sigma_lambda` / `Sigma_p` -- at a
#' `n_quad^(p_lambda + p_p)` per-species grid cost. Each species'
#' marginal -- value, abundance/detection score, and the per-site observed-
#' information block carrying the \eqn{\mathrm{Var}(N_i\mid y_i)} abundance/
#' detection coupling -- is supplied as the per-group oracle, so there is one
#' marginal/quadrature/covariance implementation across the package. The
#' Gaussian community priors pin \eqn{(\mu_\lambda, \mu_p)} as fixed effects, so
#' no sum-to-zero constraint is needed; their standard errors come from the
#' marginal observed-information Hessian.
#'
#' The abundance mixing distribution is Poisson (`mixture = "P"`) or negative
#' binomial (`mixture = "NB"`). Under NB the dispersion is itself a per-species
#' random effect, `log_r_s ~ N(mu_log_r, sigma_log_r)` (i.e.
#' `r_s ~ LogNormal`): the per-species RE vector widens to
#' `b_s = (b_lambda_s, b_p_s, b_logr_s)`, the community log-dispersion `mu_log_r`
#' joins `(mu_lambda, mu_p)` as a fixed effect, and `sigma_log_r` joins the
#' community covariances as a third (scalar) block, all integrated jointly by the
#' same AGHQ engine. The per-species size is
#' \eqn{r_s = \exp(\mu_{\log r} + b^{\log r}_s)}. Poisson is the
#' \eqn{r \to \infty} limit (no dispersion coordinate).
#'
#' @param y Integer vector of counts, one entry per observed visit (long form,
#'   all species stacked).
#' @param site_idx Integer vector (same length as `y`), 1-based site index.
#' @param species_idx Integer vector (same length as `y`), 1-based species index.
#' @param X_lambda Numeric matrix `[n_sites x p_lambda]` of abundance covariates
#'   (shared across species; one row per site).
#' @param X_p Numeric matrix `[n_obs x p_p]` of detection covariates (long form,
#'   row order matching `y`).
#' @param n_sites,n_species Integer counts.
#' @param mu_lambda_init,mu_p_init Optional warm starts for the community means.
#'   Default: the column means of independent per-species [nmix_laplace()]
#'   fits.
#' @param Sigma_lambda_init,Sigma_p_init Optional warm starts for the community
#'   covariances. Default: the (ridge-regularized) sample covariance of the
#'   per-species coefficient estimates.
#' @param K_max Marginal-sum truncation (default `max(y) + 100`), applied per
#'   site (see `headroom`).
#' @param headroom Latent-N states summed above each site's own `max(y_i)`.
#'   `NULL` (default) derives it from `K_max`: an unset `K_max` caps each site at
#'   `max(y_i) + 100`, an explicit one truncates globally and uncapped. A caller
#'   that resolved the ceiling itself passes both.
#' @param max_iter Optimizer iteration cap (default 200).
#' @param optimizer Outer optimize driver over the shared native oracle:
#'   `"em"` (default) is the fast Laplace-EM (block-coordinate Newton mode +
#'   closed-form covariance M-step + Schur SE), the exact `n_quad = 1` solver.
#'   `"joint_grad"` and `"joint_fd"` are the joint `(theta, log-Cholesky Sigma)`
#'   optimizers ([tulpa_re_aghq()]) and both do the `n_quad > 1` AGHQ debias of
#'   the community covariances. `"joint_grad"` (the fast debias path) supplies
#'   the analytic Fisher-identity gradient -- one group sweep per step, no
#'   per-coordinate objective re-solve -- and requires `n_quad > 1` (at
#'   `n_quad = 1` use the EM). `"joint_fd"` finite-differences the objective
#'   (slower; the FD sweep re-solves every per-species mode per coordinate) and
#'   is kept for correctness / architecture validation and as the `n_quad = 1`
#'   joint reference.
#' @param mixture Abundance mixing distribution: `"P"` (Poisson, default),
#'   `"NB"` (negative binomial with a per-species dispersion random effect
#'   `log_r_s ~ N(mu_log_r, sigma_log_r)`), or their zero-inflated counterparts
#'   `"ZIP"` / `"ZINB"` (a per-species structural-zero random effect
#'   `logit_omega_s ~ N(mu_omega, sigma_omega)`; a share `omega_s` of a species'
#'   sites is structurally empty). None of `"NB"` / `"ZIP"` / `"ZINB"` has a
#'   closed-form EM, so each defaults `optimizer` to `"joint_grad"` and errors on
#'   `optimizer = "em"`. The default `n_quad` when unsupplied is `5` for `"NB"`
#'   and `3` for the zero-inflated `"ZIP"` / `"ZINB"` (which carry an extra
#'   per-species RE coordinate, so a coarser grid keeps the tensor tractable).
#' @param r_init Initial community-mean negative-binomial size for the joint
#'   optimizer (`mixture = "NB"` only; default `10`, a moderate overdispersion
#'   start). The optimizer carries `mu_log_r = log(r_init)` as the
#'   `(p_lambda + p_p + 1)`-th fixed effect.
#' @param sigma_logr_init Initial standard deviation of the per-species
#'   log-dispersion random effect `log_r_s` (`mixture = "NB"` / `"ZINB"` only;
#'   default `0.5`). Seeds the scalar log_r covariance block.
#' @param omega_init Initial community-mean structural-zero probability for the
#'   joint optimizer (`mixture = "ZIP"` / `"ZINB"` only; default `0.2`). The
#'   optimizer carries `mu_omega = qlogis(omega_init)` as a trailing fixed effect.
#' @param sigma_omega_init Initial standard deviation of the per-species
#'   structural-zero-logit random effect `logit_omega_s` (`mixture = "ZIP"` /
#'   `"ZINB"` only; default `0.5`). Seeds the scalar omega covariance block.
#' @param n_quad Quadrature points per random-effect dimension passed to
#'   [tulpa_re_aghq()] (default 1, the joint Laplace). A higher `n_quad`
#'   debiases the community covariances at a `n_quad^(p_lambda + p_p)`
#'   per-species grid cost, so keep it modest when the coefficient dimension is
#'   large. This sets the order of the correlated coefficient blocks
#'   (`lambda`, `p`); the scalar nuisance blocks use `n_quad_scalar`.
#' @param n_quad_scalar Quadrature points for the scalar nuisance random-effect
#'   blocks -- the negative-binomial dispersion `log_r` and the zero-inflation
#'   `logit_omega` (default 3). These integrate a 1-D posterior, so
#'   [tulpa_re_aghq()]'s per-block quadrature runs them at this coarser order
#'   rather than the full `n_quad`, trimming the `n_quad^(NB + ZI)` factor the
#'   nuisance axes would otherwise multiply into the tensor grid. Floored at 3
#'   (and the floor wins over a smaller `n_quad`): the rule is converged there,
#'   and below it the reported community-mean SE on that block can collapse by an
#'   order of magnitude with no other symptom.
#' @param lkj_eta LKJ shape regularizing each *correlated* community covariance
#'   block's correlation off the boundary (default 1, no penalty); passed
#'   through to [tulpa_re_aghq()]. Does not touch the marginal SDs.
#' @param sigma_beta Weak Gaussian ridge SD on the community means (default 100,
#'   i.e. `tau = 1e-4`, matching the other Laplace paths); stabilizes a
#'   weakly-identified community mean without materially shifting it.
#' @param omega_sigma_prior Penalized-Complexity prior `c(U, alpha)`
#'   (`P(sigma_omega > U) = alpha`) on the per-species structural-zero random
#'   effect SD (`mixture = "ZIP"` / `"ZINB"` only; default `c(1, 0.05)`, interior
#'   mode ~ 0.33). `sigma_omega` is the softest AGHQ direction and at few species
#'   can collapse to the boundary, flattening the marginal Hessian and attenuating
#'   the recovered SD; the weak prior adds curvature there (passed to
#'   [tulpa_re_aghq()]'s `sigma_prior`) without biasing an identified fit. `NULL`
#'   restores pure ML on the structural-zero variance. Ignored for Poisson / NB.
#' @param logr_sigma_prior Penalized-Complexity prior `c(U, alpha)`
#'   (`P(sigma_log_r > U) = alpha`) on the per-species log-dispersion random
#'   effect SD (`mixture = "NB"` / `"ZINB"` only; default `NULL`, pure ML).
#'   `sigma_log_r` is the same shape of parameter as `sigma_omega` -- one scalar
#'   variance over species -- and settles near its lower boundary the same way at
#'   few species. At 8 and 36 species with a simulated `sigma_logr = 0.5`,
#'   `n_quad = 3` and `n_quad_scalar = 3`, fits recovering `sigma_log_r >= 0.30`
#'   covered `mu_log_r` 33/34 while those below covered 2/5.
#'
#'   That split is a property of those two group counts and does not carry to a
#'   third. Re-measured at 18 species on the same fixture and seeds (19 fits),
#'   one fit sits below 0.30, the point-estimate error is uncorrelated with the
#'   recovered SD (Spearman -0.08), and the reported SE is very nearly a
#'   deterministic multiple of it (Spearman +0.98, R^2 = 0.99). So what a
#'   threshold split picks up at 18 species is the SE side alone: a low
#'   recovered SD buys a proportionally narrower interval around an error that
#'   did not shrink with it. Both halves moving together is what the 8-and-36
#'   pool showed and what is absent there.
#'
#'   The SE itself is
#'   \eqn{\mathrm{se}^2 = (\hat\sigma_{\log r}^2 + c) / S}, with `c` the
#'   per-species inverse information for `log_r_s`. Measured over 97 fits at 8 /
#'   18 / 36 species, `c` is 0.130 / 0.133 / 0.141 -- one constant across a
#'   factor of 4.5 in `S` -- so the shape is right and `sigma_hat` is the input
#'   that is off. Its attenuation (0.418 / 0.448 / 0.487 against a simulated
#'   0.5) is monotone in `S` and passes straight into the interval.
#'
#'   A 1.28x-too-narrow `mu_log_r` interval at 18 species on this fixture is
#'   mostly the seed block, not the estimator. `mu_log_r` is a population mean
#'   and each seed draws `S` log-dispersions around it, so the across-seed
#'   error carries a
#'   `sigma_logr^2 / S` term that the SE includes and that is itself measured on
#'   ~19 seeds: it supplies about two thirds of the spread, and the 18-species
#'   blocks drew it 18-21% wide (chi-square p = 0.12 and 0.03). Putting that
#'   draw at its expectation and rebuilding the SE at the simulated sigma gives
#'   a scale of 0.990 / 1.035 / 0.963 at 8 / 18 / 36 species. Use
#'   `simulate_ms_abun()`'s `truth$mu_log_r_real` to score against the seed's
#'   own realized species mean and keep the two apart.
#'
#'   A fit now says whether its dispersion variance is distinguishable from zero
#'   at all: `fit$ms_dispersion$sigma_log_r_boundary` carries the boundary test
#'   and a fit that fails it warns.
#'
#'   What the penalty reaches on this block appears to be the boundary rather
#'   than the calibration. Measured on 20 seeds of that fixture at 8 species with
#'   `c(U, alpha) = c(1, 0.05)`: a `sigma_log_r` of 0.01-0.06 lifts by an order of
#'   magnitude and a seed that stopped at a singular marginal Hessian under pure
#'   ML converges and covers, while paired coverage holds at 17 of 19 with the
#'   same two seeds missing, both of them at `sigma_log_r` around 0.2 where the
#'   penalty is weak. The fits that were already calibrated take a small
#'   systematic shift (`mu_log_r` by -0.006, p = 0.001; SEs a median 4.2%
#'   narrower). Hence the `NULL` default: it is a lever for a fit whose dispersion
#'   variance came back near zero or that failed outright, rather than a
#'   correction expected to hold across fits.
#'
#'   When both this and `omega_sigma_prior` are set they must be equal: the engine
#'   applies one Penalized-Complexity prior across every block it regularizes.
#'   Ignored for Poisson / ZIP.
#' @param verbose Unused (kept for backward compatibility); the engine is silent.
#'
#' @return A list of class `nmix_re_fit`: `mu_lambda`, `mu_p` (community
#'   means), `vcov` (their joint covariance from the AGHQ marginal Hessian,
#'   `(p_lambda + p_p)` square; marginalizes the community-covariance
#'   uncertainty rather than plugging in `Sigma`), `Sigma_lambda`, `Sigma_p`
#'   (community covariances), `b_lambda`, `b_p` (per-species BLUP deviations,
#'   `n_species` rows), `log_lik` (AGHQ marginal), `converged`, `K_max`,
#'   `n_quad`, `lkj_eta`, and (when `mixture = "NB"`) the dispersion summaries:
#'   `mu_log_r` (community-mean log-dispersion; its SE is the trailing `vcov`
#'   diagonal), `sigma_log_r` (the per-species log-dispersion SD), `b_logr`
#'   (per-species deviations), `r_s` (per-species sizes
#'   \eqn{\exp(\mu_{\log r} + b^{\log r}_s)}), and `r` equal to
#'   \eqn{\exp(\mu_{\log r})} (the community-mean size, the LogNormal median).
#'
#' @references
#' Royle, J. A. (2004). N-mixture models for estimating population size from
#'   spatially replicated counts. *Biometrics* 60, 108-115.
#' Doser, J. et al. (2023). spAbundance. `msNMix()`.
#'
#' @seealso [nmix_laplace()] (single species), [nmix_site_marginal()]
#'   (the per-species marginal primitive), [tulpa_re_aghq()] (the shared
#'   random-effect integrator).
nmix_laplace_re <- function(y, site_idx, species_idx,
                                  X_lambda, X_p, n_sites, n_species,
                                  mu_lambda_init = NULL, mu_p_init = NULL,
                                  Sigma_lambda_init = NULL, Sigma_p_init = NULL,
                                  K_max = NULL, headroom = NULL,
                                  max_iter = 200L,
                                  optimizer = c("em", "joint_fd", "joint_grad"),
                                  mixture = c("P", "NB", "ZIP", "ZINB"),
                                  r_init = 10, sigma_logr_init = 0.5,
                                  omega_init = 0.2, sigma_omega_init = 0.5,
                                  n_quad = 1L,
                                  n_quad_scalar = .TOBS_MIN_SCALAR_NQUAD,
                                  lkj_eta = 1,
                                  sigma_beta = 100,
                                  omega_sigma_prior = c(1, 0.05),
                                  logr_sigma_prior = NULL,
                                  verbose = FALSE) {
  mixture <- match.arg(mixture)
  is_nb <- mixture %in% c("NB", "ZINB")   # negative-binomial abundance (log_r RE)
  is_zi <- mixture %in% c("ZIP", "ZINB")  # structural-zero share (logit_omega RE)
  # Neither the NB dispersion nor the ZI structural-zero share has a closed-form
  # EM (the closed-form Sigma M-step is the Poisson Laplace special case): each is
  # an extra fixed effect / per-species RE the joint optimizer carries. So both
  # default to the analytic-gradient debias path, which (like any joint_grad)
  # needs n_quad > 1. Explicit optimizer = "em" with either is an error.
  needs_joint <- is_nb || is_zi
  if (needs_joint && missing(optimizer)) optimizer <- "joint_grad"
  # Default quadrature order: pure NB keeps 5; the zero-inflated families add a
  # further per-species RE coordinate (logit_omega), so the tensor grid is
  # n_quad^(p_lambda + p_p + [NB] + [ZI]) -- a 5^d grid is punishing there, so ZI
  # defaults to 3 (the user can raise it).
  if (needs_joint && missing(n_quad))    n_quad <- if (is_zi) 3L else 5L
  optimizer <- match.arg(optimizer)
  if (needs_joint && optimizer == "em") {
    stop("mixture = \"", mixture, "\" needs a joint optimizer (\"joint_grad\" ",
         "or \"joint_fd\"); the EM is the Poisson Laplace solver.", call. = FALSE)
  }
  y        <- as.integer(y)
  site_idx <- as.integer(site_idx)
  species_idx <- as.integer(species_idx)
  if (!is.matrix(X_lambda)) stop("`X_lambda` must be a numeric matrix.", call. = FALSE)
  if (!is.matrix(X_p))      stop("`X_p` must be a numeric matrix.", call. = FALSE)
  n_obs <- length(y)
  if (length(site_idx) != n_obs || length(species_idx) != n_obs) {
    stop("`site_idx` and `species_idx` must have the same length as `y`.", call. = FALSE)
  }
  if (nrow(X_p) != n_obs) stop("nrow(X_p) must equal length(y).", call. = FALSE)
  if (nrow(X_lambda) != n_sites) stop("nrow(X_lambda) must equal n_sites.", call. = FALSE)
  p_lam <- ncol(X_lambda); p_p <- ncol(X_p)
  # Marginal-sum truncation, per site. The truncation must cover the latent-N
  # posterior, which the observed counts pull ABOVE the prior-lambda tail, so a
  # qpois(lambda) cap under-covers. What it has to cover is each site's OWN
  # posterior, so the ceiling is per site (see .nmix_truncation): the lgamma
  # cache removes the per-N gamma cost but the log-sum-exp and moment passes stay
  # linear in the state count, which a shared ceiling inflates for every site but
  # the one holding the largest count.
  trunc    <- .nmix_truncation(K_max, y)
  K_max    <- trunc$K_max
  headroom <- if (is.null(headroom)) trunc$headroom else as.integer(headroom)
  # The guard verifies the truncation through a Poisson per-site marginal, which
  # would under-state what a heavier-tailed abundance needs, so the cap is not
  # taken on the NB / zero-inflated paths until the check carries the fitted
  # dispersion. Those keep the shared ceiling.
  if (is_nb || is_zi) headroom <- -1L
  # Captured for the truncation guard below, which re-enters this fit with a
  # wider window; taken here so the re-entry reproduces the original call in the
  # original caller's frame rather than whatever frame the guard runs in.
  self_call   <- match.call()
  self_caller <- parent.frame()
  # A count above K_max has zero probability under the truncated marginal (the
  # latent N is summed to K_max), so a user-supplied K_max below the largest
  # observed count makes the per-(species,site) marginal structurally -Inf and the
  # joint optimum singular. Fail with an actionable message rather than the opaque
  # "singular marginal Hessian" the AGHQ solve would otherwise return. The default
  # (max(y) + 100) never trips this.
  y_max <- max(y)
  if (K_max < y_max) {
    stop(sprintf(paste0("K_max (%d) is below the largest observed count (%d). ",
                        "The N-mixture marginal sums the latent N only to K_max, ",
                        "so a count above K_max has zero probability. Raise K_max ",
                        "above max(y)."), K_max, y_max), call. = FALSE)
  }

  # ---- warm start: independent per-species fixed-effect fits ----
  if (is.null(mu_lambda_init) || is.null(Sigma_lambda_init) ||
      is.null(mu_p_init) || is.null(Sigma_p_init)) {
    B_lam <- matrix(NA_real_, n_species, p_lam)
    B_p   <- matrix(NA_real_, n_species, p_p)
    for (s in seq_len(n_species)) {
      sel <- species_idx == s
      # Warm start only seeds the community fit, so a boundary-weight warning on
      # a throwaway per-species init fit is not actionable -- suppress it (the
      # community fit's own K_max governs the final truncation).
      f <- tryCatch(
        suppressWarnings(
          nmix_laplace(y = y[sel], site_idx = site_idx[sel],
                             X_lambda = X_lambda, X_p = X_p[sel, , drop = FALSE],
                             mixture = "P", K_max = K_max, headroom = headroom,
                             verbose = FALSE)),
        error = function(e) NULL)
      if (!is.null(f)) {
        if (all(is.finite(f$beta_lambda))) B_lam[s, ] <- f$beta_lambda
        if (all(is.finite(f$beta_p)))      B_p[s, ]   <- f$beta_p
      }
    }
    ok_l <- stats::complete.cases(B_lam)
    ok_p <- stats::complete.cases(B_p)
    if (is.null(mu_lambda_init))
      mu_lambda_init <- if (any(ok_l)) colMeans(B_lam[ok_l, , drop = FALSE])
                        else c(log(mean(y) + 1), rep(0, p_lam - 1L))
    if (is.null(mu_p_init))
      mu_p_init <- if (any(ok_p)) colMeans(B_p[ok_p, , drop = FALSE]) else rep(0, p_p)
    if (is.null(Sigma_lambda_init)) Sigma_lambda_init <- .nmix_re_cov0(B_lam[ok_l, , drop = FALSE], p_lam)
    if (is.null(Sigma_p_init))      Sigma_p_init      <- .nmix_re_cov0(B_p[ok_p, , drop = FALSE], p_p)
  }

  # One native oracle (NMixCommunityOracle) is the shared backend: it assembles
  # the per-species marginal / score / observed-info / complete-data Fisher in
  # C++, and BOTH optimize drivers consume it -- there is no second marginal
  # source. The driver is selected by `optimizer`.
  n_quad <- as.integer(n_quad)

  # Guard the per-site truncation: the fitted coefficients have to carry the same
  # score under the shared ceiling as under the capped window, or the cap moved
  # the answer. Checked here rather than in the oracle, which exposes no such
  # diagnostic; a fit that fails is redone wider, escalating to the uncapped
  # ceiling, never reported truncated.
  guard <- function(out) {
    if (headroom < 0L) return(out)
    gap <- tryCatch(
      .nmix_community_score_gap(
        lf = list(y = y, site_idx = site_idx, species_idx = species_idx,
                  X_p = X_p),
        X_lambda = X_lambda,
        coef_lambda = sweep(as.matrix(out$b_lambda), 2L,
                            as.numeric(out$mu_lambda), "+"),
        coef_p = sweep(as.matrix(out$b_p), 2L, as.numeric(out$mu_p), "+"),
        K_max = K_max, headroom = headroom),
      error = function(e) NA_real_)
    if (!is.finite(gap) || gap <= .NMIX_SCORE_TOL) return(out)
    h_next <- .nmix_widen_headroom(headroom, K_max)
    if (is.null(h_next)) return(out)
    self_call$headroom <- h_next
    eval(self_call, self_caller)
  }
  orc <- cpp_nmix_community_oracle(y, site_idx, species_idx, X_lambda, X_p,
                                   n_sites, n_species, K_max,
                                   nb = is_nb, zi = is_zi, headroom = headroom)

  if (optimizer == "em") {
    # Default. EM is an outer driver over the shared oracle: block-coordinate
    # Newton mode (complete-data Fisher) + closed-form Sigma M-step + Schur SE.
    # It reaches the same n_quad = 1 (joint Laplace) stationary point as the
    # joint optimizer (agq_plan.md 4.3) without the FD-gradient objective sweep
    # that dominates joint_fd's runtime, so it is the production path. Quadrature
    # (n_quad > 1) is a joint-driver feature; the EM is Laplace (n_quad = 1).
    if (n_quad != 1L) {
      stop("optimizer = \"em\" is the n_quad = 1 (Laplace) solver; AGHQ ",
           "(n_quad > 1) needs optimizer = \"joint_fd\".", call. = FALSE)
    }
    # tol = 1e-4 on max|dSigma| between EM iters: the M-step Sigma update
    # oscillates with amplitude > 1e-6 on sparse fixtures (few visits / many
    # species) even after the estimate is statistically stable, so the
    # historical 1e-6 default never fires and the EM grinds through max_iter,
    # eating ~250s per fit at S=12 / N=60 / J=4 (12-arm covariance updates,
    # Newton inner loop at inner_tol=1e-8). 1e-4 is still ~3 orders of
    # magnitude tighter than any downstream coverage / recovery test gate
    # (which read Sigma at ~0.1 precision) and shrinks the EM to a clean
    # converge in 10-30 iters across the recovery / coverage fixtures.
    em <- cpp_nmix_community_em(
      orc, mu_init = c(as.numeric(mu_lambda_init), as.numeric(mu_p_init)),
      Sigma_lambda_init = as.matrix(Sigma_lambda_init),
      Sigma_p_init      = as.matrix(Sigma_p_init),
      max_iter = as.integer(max_iter), tol = 1e-4,
      inner_max = 50L, inner_tol = 1e-6, sigma_beta = sigma_beta,
      verbose = isTRUE(verbose))
    out <- list(
      mu_lambda = as.numeric(em$mu_lambda), mu_p = as.numeric(em$mu_p),
      vcov = em$vcov, Sigma_lambda = em$Sigma_lambda, Sigma_p = em$Sigma_p,
      b_lambda = em$b_lambda, b_p = em$b_p, log_lik = em$log_lik,
      converged = isTRUE(em$converged), n_iter = as.integer(em$n_iter),
      K_max = K_max,
      n_quad = 1L, lkj_eta = lkj_eta, optimizer = "em", mixture = "P")
    class(out) <- c("nmix_re_fit", "list")
    out <- guard(out)
    return(out)
  }

  # Joint optimizers (joint_fd / joint_grad) drive the SAME native oracle through
  # the shared compiled AGHQ engine (tulpa_re_aghq) by quadrature, and both do the
  # n_quad > 1 variance-component debias the EM cannot.
  #   joint_grad: the analytic Fisher-identity gradient (cpp_aghq_objective_grad)
  #     -- one group sweep per optim step, no per-coordinate objective re-solve.
  #     This is the fast production debias path. It needs n_quad > 1: the analytic
  #     gradient omits the Laplace curvature term, so at n_quad = 1 it disagrees
  #     with the objective (use the EM, which is the exact n_quad = 1 solver).
  #   joint_fd: optim finite-differences the objective -- the FD sweep re-solves
  #     every per-species mode for each perturbed coordinate (~100% of the
  #     residual runtime), so it is slower; kept for correctness / architecture
  #     validation and as the n_quad = 1 joint reference.
  if (optimizer == "joint_grad" && n_quad <= 1L) {
    stop("optimizer = \"joint_grad\" is the analytic-gradient AGHQ debias path ",
         "(n_quad > 1); at n_quad = 1 the analytic gradient omits the Laplace ",
         "curvature term. Use optimizer = \"em\" (the exact n_quad = 1 solver).",
         call. = FALSE)
  }
  grad_mode <- if (optimizer == "joint_grad") "analytic" else "fd"
  # NB / ZI each add a trailing per-species random-effect coordinate mirroring the
  # abundance / detection blocks: NB the log-dispersion log_r_s ~ N(mu_log_r,
  # sigma_log_r) (r_s = exp(mu_log_r + b_logr_s)), ZI the structural-zero logit
  # logit_omega_s ~ N(mu_omega, sigma_omega) (omega_s = plogis(mu_omega +
  # b_omega_s)). The oracle widens the per-species RE vector to d = p_lambda + p_p
  # + (NB) + (ZI); each community mean joins theta as a trailing fixed effect and
  # each gets its own scalar (diagonal) covariance block. Coordinate / block order
  # is [lambda | p | log_r? | omega?], matching the oracle's idx_logr / idx_omega.
  theta0 <- c(as.numeric(mu_lambda_init), as.numeric(mu_p_init))
  re_terms <- list(
    list(n_groups = n_species, n_coefs = p_lam, correlated = p_lam > 1L),
    list(n_groups = n_species, n_coefs = p_p,   correlated = p_p   > 1L))
  Sigma0 <- list(as.matrix(Sigma_lambda_init), as.matrix(Sigma_p_init))
  scalar_block <- list(n_groups = n_species, n_coefs = 1L, correlated = FALSE)
  logr_blk <- omega_blk <- NA_integer_
  if (is_nb) {
    theta0   <- c(theta0, log(r_init))
    re_terms <- c(re_terms, list(scalar_block))
    Sigma0   <- c(Sigma0, list(matrix(sigma_logr_init^2, 1L, 1L)))
    logr_blk <- length(re_terms)
  }
  if (is_zi) {
    theta0    <- c(theta0, stats::qlogis(min(max(omega_init, 1e-3), 1 - 1e-3)))
    re_terms  <- c(re_terms, list(scalar_block))
    Sigma0    <- c(Sigma0, list(matrix(sigma_omega_init^2, 1L, 1L)))
    omega_blk <- length(re_terms)
  }
  # Per-block quadrature order (tulpa_re_aghq accepts a per-RE-block n_quad
  # vector). The correlated coefficient blocks (lambda, p) use the requested
  # n_quad; the scalar nuisance blocks (log_r, omega) integrate a 1-D posterior
  # on which the Gauss-Hermite rule is converged at .TOBS_MIN_SCALAR_NQUAD, so
  # they run at the coarser n_quad_scalar. The tensor grid is
  # prod_b n_quad_b^dim_b, so trimming each scalar axis removes the extra n_quad
  # factor it would otherwise multiply in (e.g. ZINB at n_quad = 5,
  # n_quad_scalar = 3 gives 5^2 * 5^2 * 3 * 3 = 5625 nodes vs 5^6 = 15625).
  # Floored at .TOBS_MIN_SCALAR_NQUAD rather than clamped into [2, n_quad]:
  # below that order the rule does not represent the integrand, and what it
  # returns is a silently collapsed community-mean SE (#234). The floor overrides
  # a smaller n_quad for the same reason -- a coefficient block may be run at the
  # plain Laplace order, a scalar block whose marginal SE is reported may not.
  nq_vec <- rep(as.integer(n_quad), length(re_terms))
  nqs    <- .nmix_scalar_nquad(n_quad, n_quad_scalar)
  if (!is.na(logr_blk))  nq_vec[logr_blk]  <- nqs
  if (!is.na(omega_blk)) nq_vec[omega_blk] <- nqs
  # Regularize a weakly-identified scalar variance component. Both scalar blocks
  # are a single variance over species, both are among the softest AGHQ
  # directions, and at few species either can settle near its lower boundary,
  # flattening the marginal Hessian and attenuating the recovered SD. A weak
  # Penalized-Complexity prior adds curvature there (the +log sigma Jacobian
  # repels sigma -> 0, the -lambda sigma term caps inflation). Passing NULL for a
  # block restores pure ML on that variance.
  #
  # The collapse is not cosmetic on the log_r block: at 8 and 36 species with a
  # simulated sigma_logr = 0.5, the fits recovering sigma_log_r >= 0.30 cover
  # mu_log_r 33/34 while those below cover 2/5, with the point estimate 2.2x
  # further out and the interval 28% narrower (#235, NOTES_measurements.md).
  # At 18 species that joint movement is absent -- the error is uncorrelated
  # with the recovered SD there and only the interval tracks it -- so the pooled
  # split is not a property of the family at every group count (#280, #250).
  # The separate report of a 1.28x-narrow mu_log_r interval at 18 species is
  # mostly that seed block's own species draw rather than the SE: se^2 is
  # (sigma_hat^2 + c)/S with one c across 8/18/36 species, and putting the draw
  # at its expectation leaves a scale of 0.990/1.035/0.963 (#285).
  #
  # The penalty is graded by proximity to zero, and on that block it reaches the
  # boundary rather than the calibration: at c(1, 0.05) a sigma_log_r of 0.01-0.06
  # lifts by an order of magnitude and a fit that stopped at a singular marginal
  # Hessian comes back, while paired coverage holds at 17 of 19 with the same two
  # seeds missing at sigma_log_r around 0.2, and the already-calibrated fits take
  # a small systematic shift (mu_log_r by -0.006, SEs a median 4.2% narrower).
  # Hence NULL here, against c(1, 0.05) on omega, where the identified fits
  # measured unmoved (dev_notes/finding_ms_abun_zip_regularization.md).
  #
  # tulpa_re_aghq()'s `sigma_prior` carries ONE `prior_sigma` for every block its
  # `blocks` vector lists, so two blocks can be regularized together only at the
  # same (U, alpha); asking for different ones is an error rather than a silent
  # choice of one.
  pri_blk <- integer(0)
  pri_val <- NULL
  if (is_zi && !is.na(omega_blk) && !is.null(omega_sigma_prior)) {
    pri_blk <- c(pri_blk, omega_blk)
    pri_val <- omega_sigma_prior
  }
  if (is_nb && !is.na(logr_blk) && !is.null(logr_sigma_prior)) {
    if (!is.null(pri_val) && !isTRUE(all.equal(as.numeric(pri_val),
                                               as.numeric(logr_sigma_prior)))) {
      stop("omega_sigma_prior and logr_sigma_prior must be equal when both are ",
           "set: the engine applies one Penalized-Complexity prior across every ",
           "regularized block. Got c(", paste(pri_val, collapse = ", "),
           ") and c(", paste(logr_sigma_prior, collapse = ", "), ").",
           call. = FALSE)
    }
    pri_blk <- c(pri_blk, logr_blk)
    pri_val <- logr_sigma_prior
  }
  sigma_prior <- if (length(pri_blk))
    list(blocks = pri_blk, prior_sigma = pri_val) else NULL
  fit <- tulpa::tulpa_re_aghq(
    theta0  = theta0,
    re_terms = re_terms,
    Sigma0  = Sigma0,
    oracle  = orc, gradient = grad_mode,
    n_quad = nq_vec, lkj_eta = lkj_eta,
    theta_prior_sd = sigma_beta, sigma_prior = sigma_prior,
    max_iter = as.integer(max_iter))

  if (is.null(fit)) {
    # The engine declines for three reasons and names the groups behind two of
    # them in the warning it has just raised: a singular / non-finite optimum,
    # an objective already undefined at the starting parameters, and an optimum
    # whose value is the per-group failure sentinel rather than an attained
    # marginal likelihood. Naming only the first sent a reader after a Hessian
    # that was never the cause.
    stop("Community N-mixture optimization failed: the AGHQ engine returned no ",
         "fit (a singular or non-finite optimum, or a per-species posterior ",
         "solve that failed -- the engine's warning names which species). Try a ",
         "different warm start, K_max, or n_quad.", call. = FALSE)
  }

  # Per-species solve status. A species the engine could not solve has NA BLUPs,
  # so the community means, the covariance blocks and the dispersion block are
  # not estimates of what it contributes; the fit carries the status and is not
  # reported as converged.
  # (Groups are species; this fitter is indexed, not named -- build_ms_nmix_fit()
  # attaches the species names when it assembles the reported fit.)
  gstat <- .tobs_aghq_group_status(fit)

  mu <- fit$theta
  out <- list(
    mu_lambda    = mu[seq_len(p_lam)],
    mu_p         = mu[p_lam + seq_len(p_p)],
    vcov         = fit$theta_cov,
    Sigma_lambda = fit$Sigma_list[[1L]],
    Sigma_p      = fit$Sigma_list[[2L]],
    b_lambda     = fit$blup[[1L]],
    b_p          = fit$blup[[2L]],
    # Per-species FULL joint posterior covariance/cross-Hessian across the
    # lambda + p RE terms (tulpa::tulpa_re_aghq()'s blup_cov_g/blup_cross_g
    # pt. 2) -- needed by sbc()'s posterior tier to draw a species'
    # (b_lambda_s, b_p_s) jointly with the community mean instead of
    # independently ( one level deeper: the lambda/p identifiability ridge
    # means a species' abundance and detection deviations are themselves
    # correlated). NULL when the community fit ran via the n_quad = 1
    # Laplace-EM path (cpp_nmix_community_em(), a different engine that does
    # not expose this) rather than tulpa_re_aghq().
    blup_cov_g   = fit$blup_cov_g,
    blup_cross_g = fit$blup_cross_g,
    log_lik      = fit$log_marginal,
    converged    = .tobs_aghq_converged(fit, gstat),
    group_ok     = gstat$group_ok,
    groups_failed = gstat$failed,
    n_iter       = .tobs_aghq_n_iter(fit),
    K_max        = K_max,
    # Headline n_quad is the coefficient-block order; the scalar nuisance blocks
    # ran at n_quad_scalar (the full per-block vector is nq_vec).
    n_quad       = as.integer(n_quad),
    n_quad_grid  = nq_vec,
    lkj_eta      = fit$lkj_eta,
    optimizer    = optimizer,
    mixture      = mixture
  )
  # NB: a trailing theta entry is mu_log_r (community-mean log-dispersion); its
  # log-scale SE is the corresponding marginal-Hessian diagonal in `vcov`. Its
  # scalar covariance block is sigma_log_r^2 and the per-species BLUPs give
  # r_s = exp(mu_log_r + b_logr_s). The theta index / block position depend on
  # whether ZI also added a coordinate (order [lambda | p | log_r? | omega?]).
  th_i <- p_lam + p_p
  if (is_nb) {
    th_i            <- th_i + 1L
    out$mu_log_r    <- unname(mu[th_i])
    out$sigma_log_r <- sqrt(pmax(as.numeric(fit$Sigma_list[[logr_blk]]), 0))
    out$sigma_log_r_boundary <- .tobs_aghq_variance_boundary(fit, logr_blk)
    out$b_logr      <- as.numeric(fit$blup[[logr_blk]])
    out$r_s         <- exp(out$mu_log_r + out$b_logr)
    # Community-mean size (the LogNormal median; the report's headline r).
    out$r           <- exp(out$mu_log_r)
  }
  # ZI: a trailing theta entry is mu_omega (community-mean structural-zero logit);
  # its SE is the marginal-Hessian diagonal in `vcov`. Its scalar covariance block
  # is sigma_omega^2 and the per-species BLUPs give
  # omega_s = plogis(mu_omega + b_omega_s). out$omega is the community-mean
  # structural-zero probability plogis(mu_omega).
  if (is_zi) {
    th_i            <- th_i + 1L
    out$mu_omega    <- unname(mu[th_i])
    out$sigma_omega <- sqrt(pmax(as.numeric(fit$Sigma_list[[omega_blk]]), 0))
    out$sigma_omega_boundary <- .tobs_aghq_variance_boundary(fit, omega_blk)
    out$b_omega     <- as.numeric(fit$blup[[omega_blk]])
    out$omega_s     <- stats::plogis(out$mu_omega + out$b_omega)
    out$omega       <- stats::plogis(out$mu_omega)
  }
  # Both scalar blocks are one variance over species and either can settle at
  # its lower boundary, where the fit still converges and the point estimate is
  # still ordinary. Raised once, after both are attached, so a fit collapsing on
  # both says so in one place.
  .tobs_warn_variance_boundary(Filter(Negate(is.null), list(
    sigma_log_r = out$sigma_log_r_boundary,
    sigma_omega = out$sigma_omega_boundary)))
  class(out) <- c("nmix_re_fit", "list")
  guard(out)
}

# Per-species community N-mixture oracle for the shared RE integrator. Each
# species' marginal (nmix_site_marginal) is built once; the oracle
# evaluates it at the full coefficients coef = (mu + b_s), mapping the RE
# deviation through the abundance / detection designs. The b-space curvature is
# the design-sandwiched per-site observed-information block (with the latent-N
# coupling) -- the covariate generalization of the intercept-only assembly.
.nmix_re_oracle <- function(y, site_idx, species_idx, X_lambda, X_p,
                            n_sites, n_species, p_lam, p_p, K_max,
                            headroom = NULL) {
  cl <- .tobs_clamp_eta
  d  <- p_lam + p_p
  rows_by_sp <- split(seq_len(length(y)), species_idx)
  marg <- lapply(seq_len(n_species), function(s) {
    sel <- rows_by_sp[[as.character(s)]]
    if (is.null(sel)) sel <- integer(0)
    nmix_site_marginal(
      y = y[sel], site_idx = site_idx[sel],
      X_lambda = X_lambda, X_p = X_p[sel, , drop = FALSE],
      mixture = "P", K_max = K_max, headroom = headroom)
  })

  eta_of <- function(m, coef) {
    list(lambda = cl(as.numeric(m$X_lambda %*% coef[seq_len(p_lam)])),
         p      = cl(as.numeric(m$X_p %*% coef[p_lam + seq_len(p_p)])))
  }

  list(
    grad_hess = function(s, coef) {
      m  <- marg[[s]]
      e  <- eta_of(m, coef)
      ev <- m$eval(e$lambda, e$p)
      grad <- c(as.numeric(crossprod(m$X_lambda, ev$grad_eta_lambda)),
                as.numeric(crossprod(m$X_p, ev$grad_eta_p)))
      # negH: marginal observed info (with the Var[N|y] abundance/detection
      # coupling), used by the Laplace marginal. fisher: complete-data Fisher
      # (block-diagonal, PSD), supplied for the mode-find Newton -- the observed
      # info can be indefinite away from the mode (latent-N coupling).
      negH <- matrix(0, d, d); fisher <- matrix(0, d, d)
      for (i in seq_len(m$n_sites)) {
        obs <- m$obs_by_site[[i]]; Ji <- length(obs)
        Zi  <- matrix(0, 1L + Ji, d)
        Zi[1L, seq_len(p_lam)] <- m$X_lambda[i, ]
        if (Ji > 0L) Zi[-1L, p_lam + seq_len(p_p)] <- m$X_p[obs, , drop = FALSE]
        negH   <- negH + crossprod(Zi, m$obs_info_block(i, ev) %*% Zi)
        Fdiag  <- c(ev$info_eta_lambda[i], if (Ji > 0L) ev$info_eta_p[obs] else numeric(0))
        fisher <- fisher + crossprod(Zi, Fdiag * Zi)
      }
      list(logL = ev$log_lik, grad = grad, negH = negH, fisher = fisher)
    },
    node_ll = function(s, COEF) {
      m <- marg[[s]]
      vapply(seq_len(nrow(COEF)), function(k) {
        e <- eta_of(m, COEF[k, ])
        m$eval(e$lambda, e$p)$log_lik
      }, numeric(1))
    })
}

# Method-of-moments community covariance seed: the sample covariance of the
# per-species coefficient estimates, ridge-regularized to PD.
.nmix_re_cov0 <- function(B, p) {
  if (is.null(B) || nrow(B) < 2L) return(diag(0.25, p))
  V <- stats::cov(B)
  V[!is.finite(V)] <- 0
  V + diag(max(1e-3, 1e-3 * mean(abs(diag(V)))), p)
}
