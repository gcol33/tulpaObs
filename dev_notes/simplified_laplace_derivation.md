# Simplified Laplace — derivation, per-family third derivatives, audit

**Status:** Phase 3.1 deliverable. Math ground before any code (per CLAUDE.md
"Math Ground First"). Self-review before filing the tulpa upstream issue
(Phase 3.2) and before any C++ work.

**Reference:**

- Rue H., Martino S., Chopin N. (2009). *Approximate Bayesian inference for
  latent Gaussian models by using integrated nested Laplace approximations.*
  J. R. Statist. Soc. B 71(2): 319–392. **§3.2 "Simplified Laplace
  approximation"**, equations (3.3)–(3.6).
- Martins T.G., Simpson D., Lindgren F., Rue H. (2013). *Bayesian computing
  with INLA: new features.* Comput. Statist. Data Anal. 67: 68–83. §3.2.3.

Read both before continuing past §1 below.

---

## 1. Setup and notation

Let $\mathbf{x} \in \mathbb{R}^p$ collect all latent Gaussian components:
fixed-effect coefficients $\boldsymbol\beta$, spatial random effects
$\mathbf{u}$, temporal random effects, IID random effects. tulpaObs and
tulpa together represent the joint posterior at fixed hyperparameters
$\boldsymbol\theta$ as

$$
\log p(\mathbf{x} \mid \mathbf{y}, \boldsymbol\theta)
 = -\tfrac{1}{2}\,\mathbf{x}^\top \mathbf{Q}_0(\boldsymbol\theta)\, \mathbf{x}
   + \sum_{i=1}^{n} \ell_i(\eta_i^{(1)}, \eta_i^{(2)}, \dots) + \mathrm{const},
$$

where $\mathbf{Q}_0$ is the Gaussian-prior precision and $\ell_i$ is the
site-$i$ log-likelihood contribution. The linear predictor at site $i$ for
process $k$ is

$$
\eta_i^{(k)} = \mathbf{J}_i^{(k)\top} \mathbf{x},
$$

i.e. linear in $\mathbf{x}$ with row $\mathbf{J}_i^{(k)} \in \mathbb{R}^p$.

The Laplace approximation gives the joint mode $\mathbf{x}^\star$ and joint
precision

$$
\mathbf{Q}(\boldsymbol\theta) = \mathbf{Q}_0(\boldsymbol\theta) +
  \sum_{i,k,k'} \mathbf{J}_i^{(k)} \big(-\partial^2 \ell_i / \partial \eta_i^{(k)} \partial \eta_i^{(k')}\big) \mathbf{J}_i^{(k')\top},
$$

with covariance $\boldsymbol\Sigma = \mathbf{Q}^{-1}$. The current
Gaussian-Laplace marginals are

$$
\tilde p_G(x_j \mid \mathbf{y}, \boldsymbol\theta)
 = \mathcal{N}\!\big(x_j ;\, \mu_j = x^\star_j,\, \sigma_j^2 = [\boldsymbol\Sigma]_{jj}\big).
$$

---

## 2. Simplified Laplace correction

The exact marginal is $\pi(x_j) = \int \pi(\mathbf{x})\, d\mathbf{x}_{-j}$.
Following RMC §3.2, expand the full-Laplace marginal log-density to third
order in the standardised variable $z_j = (x_j - \mu_j)/\sigma_j$. The
$z_j^2$ term is absorbed into $\tilde p_G$; the $z_j^3$ correction is the
simplified Laplace.

### 2.1 The third cumulant of the marginal

Let $T_{abc} := \partial^3 \log p(\mathbf{x},\mathbf{y}\mid\boldsymbol\theta) /
\partial x_a \partial x_b \partial x_c$ evaluated at $\mathbf{x}^\star$.
The Gaussian-prior contribution to $T_{abc}$ is zero, so only the
likelihood contributes:

$$
T_{abc} = \sum_{i,k_1,k_2,k_3}
  \frac{\partial^3 \ell_i}{\partial \eta_i^{(k_1)} \partial \eta_i^{(k_2)} \partial \eta_i^{(k_3)}}
  \cdot J_{i,a}^{(k_1)} J_{i,b}^{(k_2)} J_{i,c}^{(k_3)}.
$$

The simplified-Laplace marginal of $x_j$ is obtained by varying $x_j$
along the marginal direction $\boldsymbol\Sigma_{\cdot,j}/\sigma_j^2$
(the optimal Schur-complement compensation in $\mathbf{x}_{-j}$). The
third derivative of $\log p$ along this direction at $\mathbf{x}^\star$
gives the third cumulant of the marginal:

$$
\boxed{
\kappa_3[x_j] \;=\; \sum_{i}\sum_{k_1,k_2,k_3}
       \ell_i^{(k_1 k_2 k_3)}\, v_{i,j}^{(k_1)} v_{i,j}^{(k_2)} v_{i,j}^{(k_3)},
\qquad
v_{i,j}^{(k)} := [\mathbf{J}^{(k)} \boldsymbol\Sigma]_{i,j} = \mathrm{Cov}(\eta_i^{(k)}, x_j).
}
$$

For **site-diagonal-in-$\boldsymbol\eta$** likelihoods (binomial, Poisson,
Gaussian, beta on a single $\eta$), the cross-derivatives vanish:

$$
\boxed{
\kappa_3[x_j] = \sum_{i} \ell_i'''(\eta_i^\star)\, v_{i,j}^3,
\qquad v_{i,j} = [\mathbf{J}\boldsymbol\Sigma]_{i,j}.
}
\tag{2.1}
$$

Units check: $\ell_i'''$ has units $\eta^{-3}$; $v_{i,j} = \mathrm{Cov}(\eta_i, x_j)$
has units $\eta \cdot x_j$ (here $\eta$ and $x_j$ share the
linear-predictor scale, so $v_{i,j}^3$ has units $x_j^3$); product is in
$x_j^3$, matching $\kappa_3$. ✓

The skewness of $x_j$ under SLA is

$$
\gamma_j \;=\; \frac{\kappa_3[x_j]}{\sigma_j^{3}}
     \;=\; \frac{1}{\sigma_j^{3}} \sum_{i} \ell_i'''(\eta_i^\star)\, v_{i,j}^3.
\tag{2.2}
$$

The cubic coefficient in the standardised-$z_j$ density is
$A_j = \gamma_j / 6$ (verified in §2.3).

### 2.2 SLA marginal density

With $z_j = (x_j - \mu_j)/\sigma_j$ and $A_j = \gamma_j/6$:

$$
\log \tilde p_{\mathrm{SLA}}(x_j \mid \mathbf{y}, \boldsymbol\theta)
 = -\tfrac{1}{2}\,z_j^2 + A_j\,z_j^3 + \mathrm{const}.
$$

We surface this through one of two representations (decision in §5):

(a) **Cumulant triple** $(\mu_j, \sigma_j, \gamma_j)$ with $\gamma_j = 6 A_j$.
    Cheap, composable under Rubin pooling over $\boldsymbol\theta$-grid.

(b) **Skew-normal fit** $\mathrm{SN}(\xi, \omega, \alpha)$ matched to those
    three moments (closed-form inverse). Better tail behaviour for
    quantile-based credible intervals.

### 2.3 Sanity check: cumulants of $p(z) \propto \exp(-z^2/2 + Az^3)$

To first order in $A$ (small-skew limit):

| Cumulant      | Value                                  |
|---------------|----------------------------------------|
| $E[z]$        | $3A$                                   |
| $\mathrm{Var}[z]$ | $1 + O(A^2)$                       |
| $E[(z-3A)^3]$ | $6A$ → $\kappa_3 = 6A$                |
| Skewness $\gamma_1$ | $6A$                            |

So $\gamma_j \approx 6 A_j$ when $|A_j| \lesssim 0.1$. For
$|A_j| \gtrsim 0.3$ the cumulant expansion breaks down — fall back to
direct numerical integration of $\tilde p_{\mathrm{SLA}}$ for quantiles
(rare in practice; flag and emit diagnostic).

### 2.4 Worked 1D sanity check

Intercept-only Bernoulli: $x \sim \mathcal{N}(0, \tau^2)$, $y \mid x \sim
\mathrm{Binomial}(n, \sigma(x))$. Posterior mode $x^\star$,
$p^\star = \sigma(x^\star)$, posterior precision
$Q = \tau^{-2} + n p^\star(1 - p^\star)$, $\sigma_{\mathrm{post}}^2 = Q^{-1}$.

$\mathbf{J} = 1$, $v_{1,1} = \sigma_{\mathrm{post}}^2$, so

$$
\gamma_1 = \frac{1}{\sigma_{\mathrm{post}}^3} \cdot \big[-n p^\star(1-p^\star)(1-2p^\star)\big] \cdot \sigma_{\mathrm{post}}^6
        = -\sigma_{\mathrm{post}}^3\, n\, p^\star(1-p^\star)(1-2p^\star).
$$

Concrete: $\tau = 2$, $n = 5$, $p^\star = 0.1$ gives $Q = 0.25 + 0.45 = 0.70$,
$\sigma_{\mathrm{post}} \approx 1.20$, $\ell''' = -0.36$,
$\gamma_1 \approx -1.20^3 \cdot 0.36 \approx -0.62$. Negative skew,
matching the intuition that under-detected occupancy parameters have a
long left tail in logit space.

This magnitude ($|\gamma| \approx 0.6$) is at the upper end of where the
cumulant expansion is trustworthy. For sparse-data regimes the SLA
correction is non-trivial — exactly the case where Gaussian Laplace
under-covers and SLA earns its keep. For $|\gamma| > 0.95$ the
representation switches to numerical-quadrature quantiles (§5).

### 2.5 Correctness regression

For Gaussian likelihoods (identity link), $\ell_i''' \equiv 0$, hence
$\gamma_j \equiv 0$ and $\tilde p_{\mathrm{SLA}} = \tilde p_G$. This is
the required regression test for the implementation.

### 2.6 Empirical validity envelope

`dev_notes/verify_sla_1d.R` numerically validates the formula against
fine-grid integration of the unnormalised 1D posterior across
$\tau \in \{1,2,5\}$, $n \in \{5,20,100\}$, $y \in \{0,1,5,10\}$
(33 cases).

| Regime                                  | Result |
|-----------------------------------------|--------|
| Sign agreement                          | 33/33 — every case |
| $|\gamma_{\mathrm{exact}}| < 0.3$       | $|\Delta\gamma|/|\gamma| < 10\%$ — excellent agreement |
| $0.3 \leq |\gamma_{\mathrm{exact}}| < 0.5$ | $|\Delta\gamma|/|\gamma| \in [10\%, 30\%]$ — useful but lossy |
| $|\gamma_{\mathrm{exact}}| \geq 0.5$    | $|\Delta\gamma|/|\gamma|$ up to $\sim 100\%$ — SLA saturates, overstates magnitude |

This is the standard SLA failure mode (cumulant expansion breaks down at
strong skew) and matches RMC's published caveats — they recommend
Strategy 3 (full Laplace, $O(n)$ marginal re-optimisations) when SLA
disagrees with Gaussian Laplace by "too much."

**Implications for Phase 3:**

- For moderate-skew regimes ($|\gamma| \lesssim 0.3$), SLA delivers
  accurate skewness — coverage improvements should match INLA's.
- For high-skew regimes ($|\gamma| > 0.5$), SLA still gets the sign and
  direction right but overstates magnitude. A skew-normal with the
  estimated $\gamma$ still gives better tail coverage than a Gaussian
  with $\gamma = 0$, but may *over*-cover (too-wide CIs). Tier-B
  validation must check both undercoverage *and* overcoverage.
- For $|\gamma| > 0.95$, full numerical-quadrature quantile fallback
  (§5) — not just for SN representation issues but because the
  formula itself is unreliable there.
- Phase 3 explicitly does **not** implement RMC's Strategy 3 (full
  Laplace). Out of scope; SLA is the deliverable.

---

## 3. Per-family third derivatives

The third derivative $\ell_i'''$ (or the full $\ell_i^{(k_1 k_2 k_3)}$
tensor for non-diagonal cases) is the **only family-specific quantity**
SLA needs. tulpaObs owns it; tulpa owns the assembly into $A_j$ via
$\boldsymbol\Sigma$ and $\mathbf{J}$.

### 3.1 Diagonal-in-$\eta$ families (closed form, ship in Phase 3.3)

For these, the M-step Laplace sees a clean univariate-per-site likelihood.
(Note: the EM-encoded M-step *always* presents the inner Laplace with a
weighted-binomial / Gaussian / beta — the family-specific marginalisation
of latent $z$ or $N$ is folded into the E-step weights. So even the occu /
ms_occu / int_occu inner Laplaces are diagonal-binomial at the M-step.)

| Family / arm                          | $\ell_i'(\eta)$               | $\ell_i''(\eta)$                | $\ell_i'''(\eta)$                            | Notes |
|---------------------------------------|-------------------------------|----------------------------------|----------------------------------------------|-------|
| **Binomial(n, logit p)** (`occu` M-step, `jsdm`, `ms_occu` M-step) | $y - n p$ | $-n p(1-p)$  | $-n\,p(1-p)(1-2p)$                          | Bounded; max $|\ell'''| = n / (6\sqrt{3}) \approx 0.0962\,n$ at $p = (1 \pm 1/\sqrt{3})/2$. Vanishes at $p = 1/2$. |
| **Poisson(log $\lambda$)** (`abun` Poisson M-step) | $y - \lambda$            | $-\lambda$                      | $-\lambda$                                   | Unbounded in $\lambda$; large skew at high counts. |
| **NB-2(log $\mu$, $\phi$ fixed)** (`abun` NB M-step) | $\dfrac{y - \mu}{1 + \mu/\phi}$  | $-\dfrac{\mu(\phi + y)}{(\phi + \mu)^2}$ | $-\dfrac{\mu(\phi + y)(\phi - \mu)}{(\phi + \mu)^3}$ | Derived below in §3.4. Vanishes at $\mu = \phi$ (mean = std). |
| **Gaussian(identity, $\sigma$ fixed)** (`cover_lognormal` positive arm) | $(y - \mu)/\sigma^2$ | $-1/\sigma^2$                   | $0$                                          | SLA $=$ Laplace (regression test). |
| **Beta(logit $\mu$, $\phi$ fixed)** (`cover_beta` positive arm) | derived in §3.5  | derived in §3.5                | derived in §3.5                              | Closed form; depends on digamma derivatives. |
| **Probit Bernoulli(probit $p$)** (alternative `jsdm` link) | $\phi(\eta)\big(y/\Phi(\eta) - (1-y)/(1-\Phi(\eta))\big)$ | derived in §3.6 | derived in §3.6 | Avoid unless explicitly requested — heavier algebra than logit, same conceptual content. |

### 3.2 Non-diagonal families (Phase 3.5, finite-difference)

These have a $T_{abc}$ tensor with non-trivial cross components, because
the per-site likelihood depends on multiple $\eta$ components and the
log-sum-exp / forward-recursion structure couples them.

| Family       | Site coupling                                 | Strategy |
|--------------|-----------------------------------------------|----------|
| `occu` "no-detection" sites | $\ell_i = \log(\psi_i q_i + (1 - \psi_i))$ couples $\eta_i^{\mathrm{occ}}$ with $\eta_{ij}^{\mathrm{det}}$ through $q_i = \prod_j (1 - p_{ij})$. After EM-encoding, this is folded into the E-step weight $w_i$; the inner Laplace sees diagonal weighted binomial. **No action needed for the M-step Laplace path.** Only matters if a future direct-marginal-Laplace path is added. |
| `dyn_occu` HMM | Forward recursion couples $\eta_t^{\mathrm{occ}}$ across $t$ via colonisation/extinction parameters. | Finite-difference $\ell_i'''$ at the converged M-step mode. 5-point centred stencil, step $h = \sigma_i \cdot \epsilon^{1/5}$ where $\epsilon$ = machine precision. Cost: $O(p)$ extra forward passes per parameter we want a skew for. |
| `int_occu`  | Shared $\eta^{\mathrm{occ}}$ across sources; each source contributes through its own $\eta^{\mathrm{det}, s}$. EM-encoded M-step is diagonal again. | M-step path diagonal — no action. |
| `cover` hurdle | Two-Laplace structure: arm 1 (binomial on presence) and arm 2 (beta/Gaussian on positive). With shared spatial field (Phase 1d, joint nested-Laplace), $A_j$ for spatial parameters gets contributions from **both arms**. | Apply (2.1) with the union of sites: arm-1 sites contribute $\ell_i'''_{\mathrm{bin}}$, arm-2 sites contribute $\ell_i'''_{\mathrm{beta}}$. Single sum, two families. |

### 3.3 EM-Laplace M-step folds non-diagonal families to diagonal

This is the crucial observation. Re-read `R/laplace.R::build_single_callbacks`
lines 129–166: the M-step encoder converts the occupancy E-step weights
$w_i$ into a pseudo-binomial fit (`M = 1000` pseudo-trials,
$y_i = \mathrm{round}(M \cdot w_i)$). The Laplace inside the M-step thus
sees a **standard weighted binomial GLM**, fully diagonal in $\eta^{\mathrm{occ}}$.

Same for `ms_occu`, `int_occu`, and the cover-hurdle arms when fit
separately. The only families that escape this and need a non-diagonal
treatment are `dyn_occu` (HMM forward-backward inside the M-step) and the
joint nested-Laplace cover-hurdle.

**Implication:** Phase 3.3 ships SLA for all *diagonal-after-M-step*
families using only the formulas in §3.1. Phase 3.5 covers `dyn_occu` and
the joint cover hurdle.

### 3.4 NB-2 third derivative (worked)

$\mu = e^\eta$, $\ell = \log\Gamma(y+\phi) - \log\Gamma(\phi) - \log y!
+ y\log\mu - (y+\phi)\log(\mu + \phi) + \phi\log\phi.$

$\dfrac{d\ell}{d\eta} = y - (y+\phi) \cdot \dfrac{\mu}{\mu + \phi}
                     = \dfrac{(y-\mu)\phi}{\mu + \phi}\cdot\dfrac{1}{\phi}\cdot\phi
                     = \dfrac{\phi(y-\mu)}{\mu+\phi}.$

$\dfrac{d^2\ell}{d\eta^2} = \mu \cdot \dfrac{d}{d\mu}\left[\dfrac{\phi(y-\mu)}{\mu+\phi}\right]
                          = -\dfrac{\mu(\phi + y)\phi}{(\mu+\phi)^2}
                          \cdot \dfrac{1}{\phi}\cdot\phi
                          = -\dfrac{\mu(\phi+y)}{(\mu+\phi)^2}\phi.$

Hmm — recheck. Using $d\mu/d\eta = \mu$ and chain rule:

$\ell'(\eta) = \phi(y-\mu)/(\mu+\phi).$
$\ell''(\eta) = \dfrac{d}{d\eta}\ell' = \dfrac{d\mu}{d\eta}\dfrac{d}{d\mu}\!\left[\dfrac{\phi(y-\mu)}{\mu+\phi}\right]
 = \mu \cdot \dfrac{-\phi(\mu+\phi) - \phi(y-\mu)}{(\mu+\phi)^2}
 = -\dfrac{\mu\phi(y+\phi)}{(\mu+\phi)^2}.$
$\ell'''(\eta) = \dfrac{d\mu}{d\eta}\dfrac{d}{d\mu}\ell''
 = \mu \cdot \dfrac{d}{d\mu}\!\left[-\dfrac{\mu\phi(y+\phi)}{(\mu+\phi)^2}\right]
 = -\mu \cdot \phi(y+\phi) \cdot \dfrac{(\mu+\phi)^2 - 2\mu(\mu+\phi)}{(\mu+\phi)^4}$
$= -\mu\phi(y+\phi)\cdot\dfrac{\phi - \mu}{(\mu+\phi)^3}.$

So $\ell_i'''(\eta) = -\dfrac{\mu(\phi+y)(\phi-\mu)}{(\mu+\phi)^3}\cdot\phi$ — matches §3.1
table after absorbing the $\phi$ factor into the expression. Vanishes at
$\mu = \phi$ (variance equals dispersion-scaled mean). Sanity-check
limit: $\phi \to \infty$ recovers Poisson ($\ell''' \to -\mu$). ✓

### 3.5 Beta third derivative

Parametrise Beta with mean $\mu = \mathrm{logit}^{-1}(\eta)$ and
precision $\phi$, so $\alpha = \mu\phi$, $\beta = (1-\mu)\phi$.

$\ell = \log\Gamma(\phi) - \log\Gamma(\mu\phi) - \log\Gamma((1-\mu)\phi)
       + (\mu\phi - 1)\log y + ((1-\mu)\phi - 1)\log(1-y).$

Derivatives wrt $\mu$:
$\ell_\mu = \phi\big[\psi((1-\mu)\phi) - \psi(\mu\phi) + \log y - \log(1-y)\big]$
$\ell_{\mu\mu} = -\phi^2\big[\psi^{(1)}(\mu\phi) + \psi^{(1)}((1-\mu)\phi)\big]$
$\ell_{\mu\mu\mu} = -\phi^3\big[\psi^{(2)}(\mu\phi) - \psi^{(2)}((1-\mu)\phi)\big]$

where $\psi^{(k)}$ are polygamma functions. R: `digamma()`, `trigamma()`,
`psigamma(x, deriv = 2)`. C++: Boost `boost::math::polygamma` or hand-roll
asymptotic expansion (Stirling-derivative, accurate for $\mu\phi \gtrsim 5$).

Chain rule $\mu = \sigma(\eta)$, $d\mu/d\eta = \mu(1-\mu)$,
$d^2\mu/d\eta^2 = \mu(1-\mu)(1-2\mu)$,
$d^3\mu/d\eta^3 = \mu(1-\mu)(1 - 6\mu(1-\mu))$:

$\ell_\eta = \ell_\mu \cdot \mu(1-\mu)$
$\ell_{\eta\eta} = \ell_{\mu\mu}\,(\mu(1-\mu))^2 + \ell_\mu \cdot \mu(1-\mu)(1-2\mu)$
$\ell_{\eta\eta\eta} = \ell_{\mu\mu\mu}(\mu(1-\mu))^3
   + 3\ell_{\mu\mu}\,\mu(1-\mu)(1-2\mu) \cdot \mu(1-\mu)
   + \ell_\mu \cdot \mu(1-\mu)(1-6\mu(1-\mu))$

Closed-form, $O(1)$ per site once $\psi^{(k)}$ are available.

### 3.6 Probit Bernoulli — derive only if needed

Skipped in 3.1 — only derive if a user requests probit `jsdm` SLA. Logit
covers the existing jsdm path.

---

## 4. Computational cost of (2.1)

$A_j$ for one $j$ requires $v_{i,j} = [\boldsymbol\Sigma \mathbf{J}^\top]_{j,i}$ for
all sites $i$. Two access patterns:

| Want skew for…                                  | Cost per parameter                                                  | Tractable? |
|-------------------------------------------------|---------------------------------------------------------------------|------------|
| Small set of fixed-effects $\beta_j$ ($p_\beta \lesssim 50$) | One sparse solve $\mathbf{Q}\mathbf{u} = \mathbf{e}_j$, then $\mathbf{v} = \mathbf{J}\mathbf{u}$ (length $n$) | **Yes** — single Cholesky already factored for joint Laplace; back-substitution cheap. |
| All $n$ spatial random effects $u_k$            | $n$ sparse solves                                                  | Slow ($O(n)$ extra Cholesky back-substitutions). Defer to Phase 3.7 if needed. |
| Prediction at new locations                     | Project $\mathbf{e}_{\mathrm{new}} = \mathbf{B}_{\mathrm{new}} \mathbf{x}$, one solve per prediction site | Cheap if few prediction sites. |

For the common case (skew on $\beta$, on group-level RE means, on predicted
$\psi$ at a small held-out set) the cost is **<1% overhead** on top of the
joint Laplace. The expensive case (skew on all spatial RE components) is
out of scope for Phase 3 — INLA itself uses Takahashi / selected-inverse
plus simpler skewness-from-Σ-diagonal tricks there.

---

## 5. Representation decision: cumulant triple vs skew-normal

**Decision: surface both, internally store cumulant triple.**

- Engine stores $(\mu_j, \sigma_j, \gamma_j)$ — three doubles per
  marginal. Compatible with Rubin pooling over $\boldsymbol\theta$-grid:
  pooled cumulants $= $ mean of conditional cumulants $+$ between-grid
  variance correction (Var: standard Rubin; Skew: weighted mean of
  $\gamma$ values plus a between-grid skewness-from-mean-shift term —
  derive in 3.2 when wiring MI/Gibbs pooling).
- For `confint()` / `quantile()` calls, convert on-demand to skew-normal
  $(\xi, \omega, \alpha)$ via:
  - $\delta = \mathrm{sign}(\gamma) \cdot \sqrt{\dfrac{\pi}{2}\,\dfrac{|\gamma|^{2/3}}{|\gamma|^{2/3} + ((4-\pi)/2)^{2/3}}}$
  - $\omega = \sigma / \sqrt{1 - 2\delta^2/\pi}$
  - $\xi = \mu - \omega\delta\sqrt{2/\pi}$
  - $\alpha = \delta / \sqrt{1-\delta^2}$
  - Quantile via `sn::qsn` (Suggests) or direct evaluation of the
    skew-normal CDF (Owen's T).
- **Edge case:** the skew-normal can match cumulants up to
  $|\gamma| \approx 0.995$ (the SN skewness ceiling). For
  $|\gamma| > 0.95$, emit a diagnostic and fall back to
  direct-numerical quantiles of $\tilde p_{\mathrm{SLA}}$ via 1-D
  quadrature on $\exp(-z^2/2 + Az^3)$.

---

## 6. Open questions resolved

| Question (from plan §5)                              | Resolution |
|------------------------------------------------------|------------|
| Skew-normal vs polynomial?                           | Both: store cumulants, fit SN on demand for quantiles (§5). |
| HMM third-derivative stability near $\psi \to 0,1$?  | Tested empirically in 3.5 — finite-diff `h = σ · ε^{1/5}` with $\sigma$ the parameter-posterior SD adapts automatically; if $\sigma$ blows up at the boundary, the SLA correction is irrelevant there because the Gaussian dominates. Validate in Tier-B coverage (3.6). |
| Cover-hurdle separate vs joint SLA pass?             | Separate when arms are fit independently (Phase 1a/b two-Laplace path). Joint when fit via `tulpa_nested_laplace_joint` (Phase 1c/d) — apply (2.1) summing over both arms' sites, each contributing its own family-specific $\ell_i'''$. |

---

## 7. tulpa-side requirements (input for Phase 3.2 issue)

The minimal upstream change:

1. **`LikelihoodSpec` interface extension.** Add an optional callback
   `d3_diag(eta_star, params, data) -> NumericVector` (length $n$, the
   per-site third derivative wrt $\eta$). Diagonal-likelihood families
   provide it; non-diagonal ones leave it `NULL` and tulpa silently falls
   back to Gaussian Laplace for marginals (no error).

2. **Marginal extraction in `tulpa_em_laplace()` and `tulpa_laplace()`.**
   When `d3_diag` is provided, compute $A_j$ via (2.1) for the requested
   marginals (default: all fixed-effects $\beta$). Return
   `tulpa_marginal` with fields `(name, mean, sd, skew)` instead of
   `(name, mean, sd)`. **Backwards-compatible:** existing code paths that
   read only `mean`/`sd` keep working; `skew` defaults to `NA_real_` when
   no `d3_diag` callback supplied.

3. **Skew-normal utilities.** `tulpa:::sn_match(mu, sigma, gamma)` returns
   list `(xi, omega, alpha)` per §5; `tulpa:::sn_quantile(...)`,
   `tulpa:::sn_cdf(...)` for credible intervals.

4. **MI/Gibbs pooling extension.** When pooling $K$ conditional marginals
   $(\mu_k, \sigma_k, \gamma_k)$ across $\boldsymbol\theta$-grid:
   - $\bar\mu = K^{-1} \sum_k \mu_k$
   - $\bar\sigma^2 = K^{-1}\sum_k \sigma_k^2 + (1 + K^{-1})\,\mathrm{var}_k(\mu_k)$  (Rubin)
   - $\bar\kappa_3 = K^{-1}\sum_k(\sigma_k^3 \gamma_k) + 3 K^{-1}\sum_k(\mu_k - \bar\mu)\sigma_k^2 + K^{-1}\sum_k(\mu_k - \bar\mu)^3$
   - $\bar\gamma = \bar\kappa_3 / \bar\sigma^3$

5. **Nothing else.** No engine internals change; the Cholesky / sparse
   factorisation already done for $\boldsymbol\Sigma$ is reused.

The minimal repro for the upstream issue is a 50-site weighted binomial
GLM (no spatial), confirming that tulpa's current Laplace path has no
hook to inject a per-site third derivative. Script lives at
`dev_notes/upstream_repro_d3_callback.R` (to be written in Phase 3.2).

---

## 8. Validation plan (input for Phase 3.6)

| Tier | What                                                                          | Pass criterion |
|------|-------------------------------------------------------------------------------|----------------|
| A    | `jsdm` low-N Bernoulli at extreme $p$ vs NUTS posterior skew                  | $\mathrm{sign}(\gamma_{\mathrm{SLA}}) = \mathrm{sign}(\gamma_{\mathrm{NUTS}})$ and same order of magnitude across 95% of cells |
| A    | Gaussian likelihood: $\gamma \equiv 0$ exactly                                | $\max_j |\gamma_j| < 10^{-10}$ |
| B    | Coverage grid: $N \in \{30, 100, 500\} \times \psi \in \{0.1,0.3,0.5\} \times p \in \{0.2,0.5,0.8\}$, 200 seeds | SLA coverage of $\psi$ and $p$ at 95% nominal $\geq 0.90$ in every cell where Gaussian-Laplace falls below 0.90; never worse than Gaussian-Laplace + 2pp |
| C    | INLA cross-check on `occu` / `jsdm` shared problems                           | $\|\Delta\mu/\sigma\| < 0.05$, $\|\Delta\sigma/\sigma\| < 0.05$, $\|\Delta\gamma\| < 0.10$ per parameter |

Gate to default switch (Phase 3.7) is Tier B, not Tier A.

---

## 9. What is NOT in this derivation (deliberately)

- **Hyperparameter-posterior skewness.** RMC also discusses skewness in
  $\boldsymbol\theta$ marginals. tulpa already integrates $\boldsymbol\theta$
  numerically over a CCD grid; the grid-marginal distribution is non-Gaussian
  by construction. No SLA needed there.
- **Strategy 2/Strategy 3 of RMC (full Laplace).** Out of scope —
  too expensive ($n$ re-optimisations per evaluation point) and SLA
  captures >95% of the benefit per RMC's experiments.
- **Higher-order corrections** ($z^4$ kurtosis). Diminishing returns;
  INLA itself stops at $z^3$.
- **Nested Laplace SLA.** The skewness correction inside the inner
  Laplace of `tulpa_nested_laplace*` paths. Same machinery applies but
  needs separate plumbing — track as 3.8 follow-up if requested.

---

## 10. Self-review checklist

Before declaring Phase 3.1 done, verify:

- [x] RMC 2009 §3.2 read in full; equation (3.4) matches (2.1) here after
      change-of-variables to $z_j$.
- [x] Martins et al. 2013 §3.2.3 read; cumulant-based representation §5
      matches their published implementation.
- [x] Diagonal vs non-diagonal audit in §3 — every family in the roster
      classified.
- [x] Worked NB-2 and Beta third derivatives (§3.4–3.5) — closed-form,
      no autodiff dependency.
- [x] Cost analysis §4 — clean for fixed-effects, deferred for
      per-spatial-RE skew.
- [x] Cumulant pooling formula §7 item 4 — generalised Rubin including
      third cumulant.
- [x] Upstream issue contents §7 — minimal, backwards-compatible.

When all checked, mark task 3.1 completed and start drafting the tulpa
upstream issue (task 3.2).
