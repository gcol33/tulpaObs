# Characterize the small-sample MLE bias of the beta precision under
# *known* mu (no field, no spatial, no joint engine). Originally written
# to confirm that the residual -40% bias of phi_pos at n_pos~46 in D7
# cell B was Ferrari & Cribari-Neto (2004, JAS §3.4) small-sample MLE
# bias.
#
# The probe REFUTED that hypothesis: with known mu the 1D Brent MLE is
# slightly *upward*-biased (+9% at n=25, +4% at n=50, +1% at n=100) and
# the parametric-bootstrap correction kills it cleanly. The residual
# -40% in D7 must then come from upstream — posterior shrinkage of
# (sigma, rho, alpha, w_s) collapsing the variance of the field-
# corrected linear predictor. Confirmed by the companion probe at
# `INLAabun/example/validation/.probe_d7_phi_origin.R`: refitting phi
# against the *truth* linear predictor recovers cleanly; against the
# posterior linear predictor it does not, and bootstrap correction on
# the posterior path makes it *worse* (-44%).
#
# Conclusion: phi_pos at small n_pos needs a phi-on-the-outer-grid
# refactor, not a small-sample bias correction. Tracked at tulpaObs#7.
#
# Truth phi = 30 to match D7 cell B.


# Beta MLE for phi conditional on mu (1D Brent, same code path as the
# in-package .refit_beta_phi_postfield).
beta_phi_mle <- function(y, mu, bounds = c(0.1, 1e4)) {
  eps <- 1e-6
  mu_c <- pmin(pmax(mu, eps), 1 - eps)
  y_c  <- pmin(pmax(y,  eps), 1 - eps)
  nloglik <- function(phi) {
    a <- mu_c * phi
    b <- (1 - mu_c) * phi
    -sum(lgamma(phi) - lgamma(a) - lgamma(b) +
         (a - 1) * log(y_c) + (b - 1) * log(1 - y_c))
  }
  stats::optimize(nloglik, interval = bounds, tol = 1e-4)$minimum
}

# Parametric bootstrap bias correction. Returns 2*phi_hat - mean(phi_boot).
# B = 100 keeps Monte Carlo noise on the correction below ~3% at n_pos=46
# with truth_phi=30.
beta_phi_bcorrect <- function(y, mu, B = 100L, bounds = c(0.1, 1e4)) {
  phi_hat <- beta_phi_mle(y, mu, bounds)
  n <- length(y)
  eps <- 1e-6
  mu_c <- pmin(pmax(mu, eps), 1 - eps)
  a <- mu_c * phi_hat
  b <- (1 - mu_c) * phi_hat
  phi_boot <- numeric(B)
  for (k in seq_len(B)) {
    y_b <- rbeta(n, a, b)
    y_b <- pmin(pmax(y_b, eps), 1 - eps)
    phi_boot[k] <- beta_phi_mle(y_b, mu_c, bounds)
  }
  list(phi_hat = phi_hat,
       phi_bc  = 2 * phi_hat - mean(phi_boot),
       boot_mean = mean(phi_boot),
       boot_sd = sd(phi_boot))
}

# Cox-Snell-style analytical bias correction for 1D phi MLE with known mu.
# Derivation: third-order Bartlett identity gives
#   bias(phi_hat) = (E[U''] - 2 E[U^3]) / (6 I^2) + O(1/n^2)
# where U_i = psi(phi) - mu_i psi(mu_i phi) - (1-mu_i) psi((1-mu_i) phi)
#             + mu_i log(y_i) + (1-mu_i) log(1-y_i).
# Implementation uses Bartlett's identity E[U U'] = -(E[U''] + E[U^3])/3,
# so equivalently bias = E[U U']/I^2 + E[U'']/(2 I^2).
# All expectations are taken at phi = phi_hat, mu = mu_i (plug-in).
beta_phi_coxsnell <- function(y, mu, bounds = c(0.1, 1e4)) {
  phi <- beta_phi_mle(y, mu, bounds)
  eps <- 1e-6
  mu  <- pmin(pmax(mu, eps), 1 - eps)
  a <- mu * phi
  b <- (1 - mu) * phi
  # Per-obs Fisher info (= -E[U_i']):
  I_i <- digamma_trigamma <- mu^2 * trigamma(a) + (1 - mu)^2 * trigamma(b) - trigamma(phi)
  I_total <- sum(I_i)
  # E[U_i'']: third derivative of log-density w.r.t. phi
  #   U_i' = -psi'(phi) + mu^2 psi'(mu phi) + (1-mu)^2 psi'((1-mu) phi)
  #   U_i'' = -psi''(phi) + mu^3 psi''(mu phi) + (1-mu)^3 psi''((1-mu) phi)
  # psi'' = -tetragamma; in R: psigamma(x, deriv=2)
  Upp_i <- -psigamma(phi, 2) + mu^3 * psigamma(a, 2) + (1 - mu)^3 * psigamma(b, 2)
  EUpp <- sum(Upp_i)
  # E[U_i^3]: third moment of score under correct model.
  # Under Beta(a, b): log(y) ~ generalized log-gamma; let
  #   d1 = log(y) has E = psi(a) - psi(a+b), Var = psi'(a) - psi'(a+b), ...
  # But here a+b = phi for every i, so log y_i conditional on mu_i has
  #   E[log y_i]      = psi(mu_i phi) - psi(phi)
  #   E[log(1-y_i)]   = psi((1-mu_i) phi) - psi(phi)
  # Score per obs: U_i = c_i + mu_i log y_i + (1-mu_i) log(1-y_i)
  # where c_i is the constant part. We need central third moment of U_i.
  # Closed forms involve the third cumulants of log y_i and log(1-y_i)
  # and their covariance. Use polygamma identities:
  #   cum_3(log y | Beta(a,b)) = psigamma(a, 2) - psigamma(a+b, 2)
  #   cum_3(log(1-y) | Beta(a,b)) = psigamma(b, 2) - psigamma(a+b, 2)
  # Cross-third-cumulant cum_{1,2}(log y, log(1-y)) = -psigamma(a+b, 2)? No,
  # this needs the joint MGF. For Beta(a,b):
  #   log y and log(1-y) are NOT independent. Joint cumulant generating
  #   function: K(s, t) = log B(a+s, b+t) - log B(a, b)
  # so kappa_{p,q} = d^{p+q} K / d s^p d t^q at (0,0). This gives:
  #   kappa_{3,0} = psigamma(a, 2) - psigamma(a+b, 2)
  #   kappa_{0,3} = psigamma(b, 2) - psigamma(a+b, 2)
  #   kappa_{2,1} = -psigamma(a+b, 2)
  #   kappa_{1,2} = -psigamma(a+b, 2)
  # Since a + b = phi, psigamma(a+b, 2) = psigamma(phi, 2).
  k30 <- psigamma(a, 2) - psigamma(phi, 2)
  k03 <- psigamma(b, 2) - psigamma(phi, 2)
  k21 <- -psigamma(phi, 2)
  k12 <- -psigamma(phi, 2)
  # Third cumulant of U_i = mu_i * X + (1-mu_i) * Y where X = log y, Y = log(1-y):
  #   kappa_3(mu X + (1-mu) Y) = mu^3 k30 + 3 mu^2 (1-mu) k21
  #                              + 3 mu (1-mu)^2 k12 + (1-mu)^3 k03
  EU3_i <- mu^3 * k30 + 3 * mu^2 * (1 - mu) * k21 +
           3 * mu * (1 - mu)^2 * k12 + (1 - mu)^3 * k03
  EU3 <- sum(EU3_i)
  # Cox-Snell bias:
  bias <- (EUpp - 2 * EU3) / (6 * I_total^2)
  list(phi_hat = phi, phi_cs = phi - bias, bias_est = bias)
}

# Simulate a sparse-positive panel matching D7 cell B's logistic-link mu
# distribution: eta_pos ~ N(eta_mean, eta_sd) so mu spans roughly the same
# range you see in practice. With eta = N(0.4, 0.7), mu lives in [0.3, 0.85]
# most of the time, similar to D7 cell B's positive-arm support.
simulate_beta_known_mu <- function(n, phi = 30, eta_mean = 0.4, eta_sd = 0.7,
                                   seed = 1L) {
  set.seed(seed)
  eta <- rnorm(n, eta_mean, eta_sd)
  mu  <- plogis(eta)
  y   <- rbeta(n, mu * phi, (1 - mu) * phi)
  list(y = y, mu = mu)
}

phi_true <- 30
n_grid <- c(25L, 50L, 100L, 200L, 400L)
n_seeds <- 200L

cat(sprintf("\n--- Beta phi MLE small-sample bias, phi_true = %g ---\n", phi_true))
cat(sprintf("%6s %8s %8s %8s | %8s %8s | %8s %8s\n",
            "n", "MLE_med", "MLE_mean", "MLE_relB",
            "BC_mean", "BC_relB", "CS_mean", "CS_relB"))
cat(paste(rep("-", 78), collapse = ""), "\n")

for (n in n_grid) {
  mle <- numeric(n_seeds)
  bc  <- numeric(n_seeds)
  cs  <- numeric(n_seeds)
  for (s in seq_len(n_seeds)) {
    sim <- simulate_beta_known_mu(n, phi = phi_true, seed = 10000L + s)
    fit_mle <- beta_phi_mle(sim$y, sim$mu)
    fit_bc  <- beta_phi_bcorrect(sim$y, sim$mu, B = 60L)
    fit_cs  <- beta_phi_coxsnell(sim$y, sim$mu)
    mle[s]  <- fit_mle
    bc[s]   <- fit_bc$phi_bc
    cs[s]   <- fit_cs$phi_cs
  }
  cat(sprintf("%6d %8.2f %8.2f %+8.1f%% | %8.2f %+8.1f%% | %8.2f %+8.1f%%\n",
              n,
              median(mle), mean(mle), 100 * (mean(mle) - phi_true) / phi_true,
              mean(bc),    100 * (mean(bc) - phi_true) / phi_true,
              mean(cs),    100 * (mean(cs) - phi_true) / phi_true))
}

cat("\n--- D7 cell B regime: n_pos ~ 46, n=300 panel-scale ---\n")
cat("(matches D7 sparse cell B truth phi=30)\n")
