# Upstream tulpa spec — utilities for simplified Laplace

**Status:** Phase 3.2 design doc. To be filed as a tulpa GitHub issue
once tulpa WIP (current uncommitted Laplace-routing work in
`R/em_laplace.R`) is committed and rebased.

**Companion:** `dev_notes/simplified_laplace_derivation.md`.

---

## 1. Re-scoped goal

After deriving the SLA formula (companion note §2.1) the boundary
between tulpa and tulpaObs becomes much narrower than initially planned:

- **All assembly of $\gamma_j$ lives in tulpaObs.** The formula
  $\gamma_j = \sigma_j^{-3} \sum_i \ell_i'''(\eta_i^\star)\, v_{i,j}^3$
  needs $\ell_i'''$ (family-specific, owned by tulpaObs),
  $v_{i,j} = [\mathbf{J}\boldsymbol\Sigma]_{i,j}$ (computable from
  $\mathbf{J} = \mathbf{X}$ for fixed-effects marginals and
  $\boldsymbol\Sigma = \mathbf{H}_\beta^{-1}$, both already returned by
  `tulpa_laplace()`), and $\sigma_j^2 = \boldsymbol\Sigma_{jj}$.
- **tulpa does NOT need a `d3_diag` callback.** Original Phase 3.2 spec
  assumed it; the derivation shows it's unnecessary for the M-step-Laplace
  path because tulpaObs computes $\ell_i'''$ at the *original* likelihood
  (not the pseudo-binomial encoding — see §3 below for why this matters).

What tulpa *does* need to grow:

1. **Skew-normal utilities** — pure-R helpers callable from tulpaObs:
   `sn_match()`, `sn_quantile()`, `sn_cdf()`. Lightweight; could even
   live in tulpaObs as `:::`-private. **Place in tulpa** because they
   are generic and any tulpa-derived package (tulpaglmm, future) may
   want them.
2. **Cumulant-triple pooling** in MI/Gibbs correction — `tulpa::em_correction.R`
   currently pools $(\mu, \sigma)$ across imputations via Rubin's rules;
   extend to pool $(\mu, \sigma, \gamma)$ via the generalised formula in
   companion §7.4.
3. **Nested-Laplace per-grid-point Hessian exposure** — for SLA in
   `tulpa_nested_laplace()` we need $\mathbf{H}_\beta$ at each
   $\boldsymbol\theta$-grid point, not just the integrated marginals. The
   nested-Laplace driver currently discards per-grid Hessians after
   computing the integrated log-marginal. Exposing them is a small
   plumbing addition — `nl_fit$grid_hessians` as a list.

That's it. No likelihood-spec interface change, no new callbacks.

---

## 2. Minimal proposed API additions to tulpa

### 2.1 Skew-normal utilities (`R/skew_normal.R`, new file)

```r
#' Match cumulants (mu, sigma, gamma) to skew-normal (xi, omega, alpha)
#'
#' Inverse of the skew-normal moment formulas. Returns NULL with a
#' warning when |gamma| exceeds the skew-normal ceiling (~0.995).
#'
#' @param mu Mean
#' @param sigma Standard deviation (positive)
#' @param gamma Skewness coefficient
#' @return List with elements xi, omega, alpha; or NULL if |gamma| > 0.995
#' @export
sn_match <- function(mu, sigma, gamma) {
  gmax <- 0.9952717  # (4-pi)/2 * (2/pi)^(3/2) / (1 - 2/pi)^(3/2)
  if (abs(gamma) >= gmax) {
    warning("|gamma| = ", round(abs(gamma), 3),
            " exceeds skew-normal ceiling (", round(gmax, 3),
            "); returning NULL — caller should use direct-quadrature ",
            "quantiles instead.", call. = FALSE)
    return(NULL)
  }
  c1 <- ((4 - pi) / 2)^(2 / 3)
  delta_sq <- (pi / 2) * abs(gamma)^(2 / 3) / (abs(gamma)^(2 / 3) + c1)
  delta <- sign(gamma) * sqrt(delta_sq)
  omega <- sigma / sqrt(1 - 2 * delta^2 / pi)
  xi    <- mu - omega * delta * sqrt(2 / pi)
  alpha <- delta / sqrt(1 - delta^2)
  list(xi = xi, omega = omega, alpha = alpha)
}

#' Skew-normal quantile
#' @param p Probabilities (in [0, 1])
#' @param sn Skew-normal parameter list from sn_match()
#' @return Quantile vector
#' @export
sn_quantile <- function(p, sn) {
  # Newton iteration on sn_cdf(); fast for vectorised p
  ...
}

#' Skew-normal CDF (Owen's T via sn::psn fallback or direct)
#' @export
sn_cdf <- function(q, sn) { ... }
```

Suggests dependency on `sn` (Azzalini's package) for the Owen's T
implementation. Alternative: hand-roll Owen's T (~30 lines, well-known
recipe), zero new dependency. **Recommendation:** hand-roll. Per CLAUDE.md
"No Dependency Shortcuts" — `sn` is one package for ~50 lines of code.

### 2.2 Cumulant pooling in MI/Gibbs (`R/em_correction.R` extension)

Current MI pooling (Rubin's rules) in `tulpa_em_correction()`:

```r
# Existing:
mu_pooled    <- mean(mu_k)
sigma2_pooled <- mean(sigma_k^2) + (1 + 1/K) * var(mu_k)
```

Extension when imputations carry `gamma_k`:

```r
# New (engaged only when all imputations have non-NA gamma_k):
kappa3_pooled <- mean(sigma_k^3 * gamma_k) +
                 3 * mean((mu_k - mu_pooled) * sigma_k^2) +
                 mean((mu_k - mu_pooled)^3)
gamma_pooled  <- kappa3_pooled / sqrt(sigma2_pooled)^3
```

Backward-compatible: if `gamma_k` not present (existing callers), skip
the new code path and return as before.

### 2.3 Nested-Laplace per-grid-point Hessian (`R/nested_laplace.R`)

Add a `keep_grid_hessians = FALSE` argument. When TRUE, the
`tulpa_nested_laplace()` return list gains:

```r
fit$grid_hessians  # list of length n_grid; element [[k]] = H_beta at theta_k
fit$grid_modes     # list of length n_grid; element [[k]] = beta_hat at theta_k
fit$grid_weights   # already returned
```

Allows tulpaObs to compute per-grid SLA marginals and pool via 2.2.
Memory cost: $O(n_{\mathrm{grid}} \cdot p_{\mathrm{fixed}}^2)$ — tiny.

---

## 3. Why `d3_diag` callback is NOT proposed

Original Phase 3.2 spec (in `dev_notes/simplified_laplace_derivation.md`
§7) suggested adding a `d3_diag(eta, params, data) -> n-vector` callback
to tulpa's likelihood interface. Re-deriving shows this is unnecessary
and would in fact be **incorrect** for the EM-Laplace path:

| Path                          | What "$\ell_i'''$" should be                                          | Where evaluated |
|-------------------------------|------------------------------------------------------------------------|-----------------|
| Direct nested-Laplace (jsdm, future direct-marginal occu) | Third derivative of the actual observation likelihood (binomial, etc.) | tulpa likelihood module |
| **EM-Laplace M-step** (occu, ms_occu, int_occu, abun, cover) | Third derivative of the *original* marginal likelihood, NOT the pseudo-binomial encoding | tulpaObs family code |

For the EM-Laplace M-step: the inner Laplace sees a pseudo-binomial with
`n_trials = M = 1000`. The pseudo-binomial third derivative scales with
$M$, but $\sigma_j$ scales as $M^{-1/2}$ — so the pseudo-skewness
$\gamma^{\mathrm{pseudo}}_j = \sigma_j^3 \cdot \ell_i^{\mathrm{pseudo},'''}$
scales as $M^{-1/2} \to 0$. The pseudo-binomial smooths skewness away,
which is wrong for the original occupancy posterior.

**Correct approach:** compute $\ell_i'''$ at the *original* occupancy
log-likelihood at the EM-converged $\boldsymbol\beta$, using the actual
data and the actual occupancy likelihood structure (per-site closed forms
in companion §3 / §3.4 of the derivation note). tulpaObs already knows
this — tulpa shouldn't.

So tulpa stays out of the family-specific $\ell_i'''$ business entirely.

---

## 4. Minimal reproducer (commit to tulpa/dev_notes/)

```r
# tulpa/dev_notes/sla_minimal_check.R
# Shows that tulpa_laplace already returns everything tulpaObs needs to
# compute SLA marginals for fixed-effects, modulo the small utility
# additions in §2.

set.seed(1)
n <- 100
X <- cbind(1, rnorm(n))
beta_true <- c(-1, 0.5)
eta <- X %*% beta_true
p <- plogis(eta)
y <- rbinom(n, 1, p)

fit <- tulpa::tulpa_laplace(y = y, n_trials = rep(1L, n), X = X,
                            family = "binomial", return_hessian = TRUE)

# tulpa already returns these:
stopifnot(!is.null(fit$mode))     # MAP beta
stopifnot(!is.null(fit$H_beta))   # marginal Hessian (post-Schur for spatial)

# tulpaObs can compute SLA from these directly:
beta_hat <- fit$mode[seq_len(ncol(X))]
Sigma <- solve(fit$H_beta)
sigma_j <- sqrt(diag(Sigma))

eta_hat <- as.numeric(X %*% beta_hat)
p_hat   <- plogis(eta_hat)
l3 <- -p_hat * (1 - p_hat) * (1 - 2 * p_hat)  # n_trials = 1

# v_{i,j} = (X Sigma)_{i,j}
XSig <- X %*% Sigma
gamma_j <- vapply(seq_len(ncol(X)), function(j) {
  sum(l3 * XSig[, j]^3) / sigma_j[j]^3
}, numeric(1))

stopifnot(all(is.finite(gamma_j)))
print(gamma_j)
```

Run this in tulpa to confirm `H_beta` exposure is sufficient. If it
errors (e.g. `H_beta` missing in some path), that's the bug to fix.

---

## 5. Proposed tulpa issue text

> ### feat(laplace): utility helpers for simplified Laplace marginals
>
> ### Summary
> Downstream tulpaObs (Phase 3, simplified Laplace) needs three small
> additions to tulpa. None require touching the Laplace engine or
> likelihood-spec interface. All are backwards-compatible.
>
> ### Required additions
> 1. `sn_match()`, `sn_quantile()`, `sn_cdf()` skew-normal utilities
>    (~80 LOC pure R, no new deps; Owen's T hand-rolled).
> 2. Cumulant-triple pooling in `tulpa_em_correction()` — engaged only
>    when imputations carry `gamma_k`; existing two-moment path
>    unchanged.
> 3. `keep_grid_hessians = FALSE` arg in `tulpa_nested_laplace()` that
>    exposes per-grid-point `H_beta` and `mode` for downstream SLA.
>
> ### Why not a `d3_diag` callback
> Originally planned, ruled out: the EM-Laplace path needs
> $\ell_i'''$ at the *original* family likelihood, not the M-step's
> pseudo-binomial encoding. Family-specific $\ell_i'''$ stays in
> tulpaObs / tulpaglmm where the per-family structure is known.
>
> ### Downstream impact
> Blocks tulpaObs Phase 3.3 (per-family $\ell_i'''$) and Phase 3.4
> (skew-aware `confint`/`predict`). See
> `tulpaObs/dev_notes/upstream_tulpa_sla_spec.md` for the full design
> and `tulpaObs/dev_notes/simplified_laplace_derivation.md` for the
> math.
>
> ### Reproducer
> `dev_notes/sla_minimal_check.R` (in tulpa, to be added): shows that
> existing `tulpa_laplace()` already exposes enough for downstream SLA
> assembly, modulo the three utilities above.

---

## 6. What I will NOT do without checking first

- **Modify tulpa code while WIP changes exist there.** Current `git status`
  in tulpa shows 11 modified files (unrelated nested-Laplace EM routing
  work). I will not stack SLA changes on top of unfinished work.
- **File the GitHub issue automatically.** Per CLAUDE.md, file via
  `gh issue create --repo gcol33/tulpa` *after* the user reviews this spec
  and confirms.
- **Add a `Suggests: sn` dependency to tulpa.** Hand-roll Owen's T per §2.1.
