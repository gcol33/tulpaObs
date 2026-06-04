# =============================================================================
# obs_families.R — Family-object constructors for tulpaObs
#
# Each family is a small list with class "tobs_family". It carries the latent-
# state type, the observation likelihood, replicate requirements, default
# engines, and any family-specific hyperparameters (K_max, beta vs lognormal
# positive part, etc.).
#
# Families do not fit models. tobs() reads the family object and dispatches to
# the appropriate engine path.
# =============================================================================


#' Construct a tobs family object
#'
#' Low-level constructor for `tobs_family` objects. End users should call the
#' specific family functions ([occu()], [abun()], [cover()], ...) rather than
#' this directly.
#'
#' @param name short slug, e.g. `"occu"`.
#' @param class_long human-readable name, e.g. `"single-season occupancy"`.
#' @param latent latent-state distribution, e.g. `"bernoulli"`, `"poisson"`,
#'   `"lognormal"`, `"hurdle"`.
#' @param observation observation likelihood, e.g. `"binomial_detection"`,
#'   `"binomial_N"`, `"beta"`, `"distance_binned"`.
#' @param replicates one of `"required"`, `"optional"`, `"single"`.
#' @param default_engine `"laplace"`, `"nested_laplace"`, or `"nuts"`.
#' @param status `"working"`, `"planned"`, or `"experimental"`.
#' @param params named list of family-specific parameters carried with the
#'   object (K_max, positive-part link, etc.).
#' @param control_keys character vector of extra `control` names this family's
#'   dispatcher accepts beyond the engine/route controls. These are added to
#'   the allowlist `tobs()` validates `control` against, so a family with a
#'   bespoke dispatcher (e.g. the cover hurdle's grid controls) is not rejected.
#'
#' @return A `tobs_family` object.
#' @keywords internal
#' @export
obs_family <- function(name,
                       class_long,
                       latent,
                       observation,
                       replicates    = c("required", "optional", "single"),
                       default_engine = c("laplace", "nested_laplace", "nuts"),
                       status         = c("working", "planned", "experimental"),
                       params         = list(),
                       control_keys   = character(0)) {
  replicates     <- match.arg(replicates)
  default_engine <- match.arg(default_engine)
  status         <- match.arg(status)

  structure(
    list(
      name           = name,
      class_long     = class_long,
      latent         = latent,
      observation    = observation,
      replicates     = replicates,
      default_engine = default_engine,
      status         = status,
      params         = params,
      control_keys   = control_keys
    ),
    class = "tobs_family"
  )
}


# ---------------------------------------------------------------------------
# Working families — dispatch to existing engine paths
# ---------------------------------------------------------------------------

#' Single-season occupancy family
#'
#' Latent Bernoulli occupancy state with binomial detection per visit. The
#' MacKenzie et al. (2002) single-season model.
#'
#' @return A `tobs_family` object.
#' @export
#' @examples
#' f <- occu()
#' f
occu <- function() {
  obs_family(
    name           = "occu",
    class_long     = "single-season occupancy",
    latent         = "bernoulli",
    observation    = "binomial_detection",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working"
  )
}


#' Dynamic (multi-season) occupancy family
#'
#' Bernoulli occupancy state with colonisation + extinction transitions across
#' seasons (the MacKenzie et al. dynamic model).
#'
#' @return A `tobs_family` object.
#' @export
dyn_occu <- function() {
  obs_family(
    name           = "dyn_occu",
    class_long     = "dynamic occupancy (HMM)",
    latent         = "bernoulli_hmm",
    observation    = "binomial_detection",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working"
  )
}


#' Integrated occupancy family
#'
#' Multiple data sources informing a shared latent occupancy state, each with
#' its own detection process.
#'
#' @return A `tobs_family` object.
#' @export
int_occu <- function() {
  obs_family(
    name           = "int_occu",
    class_long     = "integrated occupancy",
    latent         = "bernoulli",
    observation    = "multisource_detection",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working"
  )
}


#' Joint species distribution family (no detection process)
#'
#' Multivariate occurrence with shared latent factors. No observation
#' replication — treats observed presence/absence as the response.
#'
#' @return A `tobs_family` object.
#' @export
jsdm <- function() {
  obs_family(
    name           = "jsdm",
    class_long     = "joint species distribution",
    latent         = "latent_factor",
    observation    = "probit",
    replicates     = "single",
    default_engine = "nuts",
    status         = "working"
  )
}


#' Joint occupancy-detection + cover hurdle family
#'
#' Cell-level latent presence `psi`, per-visit binomial detection `p`, and
#' per-visit positive cover `f_pos` (beta or lognormal) on a third linear
#' predictor. The N-mixture analogue for vegetation cover: where [abun()]
#' couples occupancy to counts via Royle's binomial-N marginal, [occu_cover()]
#' couples it to a positive-cover observation at each detected visit. The
#' latent presence z marginalises out in closed form (two states), so the
#' marginal log-likelihood is exact.
#'
#' Per-cell likelihood:
#'
#'     any_det_i : L_i = psi_i * prod_j h_ij
#'     no_det_i  : L_i = psi_i * prod_j (1 - p_ij) + (1 - psi_i)
#'
#'     h_ij = (1 - p_ij) * 1{y_ij = 0}
#'          + p_ij       * f_pos(y_pos_ij; eta_pos_ij, ...) * 1{y_ij = 1}
#'
#' Reduces to [occu()] when the cover arm is degenerate, and to the
#' plot-level cover hurdle ([cover()]) when J = 1 and detection is perfect.
#'
#' @section v1 scope (status `"experimental"`):
#' Laplace only; non-spatial. A shared spatial field across the three arms
#' (the analogue of [cover()]'s nested-Laplace joint engine) is v2. Structured
#' terms (`bym2()`, `icar()`, `re()`, ...) error from the dispatcher with a
#' pointer to this note.
#'
#' @section Coupled fields and spatially-varying trends:
#' The spatial engine (`method = "nested_laplace"`, default
#' `control$engine = "joint_coupled"`) couples one shared areal field (the
#' cell intercept, an unweighted `icar()` / `bym2()` term) across the occupancy
#' and cover arms. ADDITIONAL coupled fields - spatially-varying coefficients,
#' e.g. a temporal trend - are added as *weighted* areal terms in the formula:
#'
#' ```r
#' ~ elev + icar(graph = adj) + icar(graph = adj, weight = year)
#' ```
#'
#' Each weighted term `icar(graph, weight = col)` is a second coupled field
#' whose contribution to a predictor row is `weight_i * z[cell_i]`. All coupled
#' fields share the same graph; each couples onto the cover arm with its own
#' scale (`alpha` for the intercept field, `alpha_trend` for a trend field),
#' integrated over the outer grid. The intercept field is reported in
#' `fit$spatial_field`, trend fields in `fit$trend_field` / `fit$trend_fields`.
#' The trend coupling grid defaults to `control$alpha.grid`; override it with
#' `control$alpha.grid.trend`.
#'
#' A single trend field may also be requested out-of-band with
#' `control = list(trend = list(weight = "<col>"))`, naming a numeric per-cell
#' covariate in the cell `data`; this is the equivalent of one weighted formula
#' term. Specify the trend field one way or the other, not both.
#'
#' @section Sites larger than cells (`group_var`):
#' By default each site (one row of `y` / `data`, one latent occupancy state) is
#' its own field node, so the graph must have one node per site. Passing
#' `group_var = "<col>"` to the `icar()` / `bym2()` term maps each site to a
#' field node named by that integer column, so several sites can share one node.
#' The field then stays length `n_cells` (the graph) while occupancy, detection,
#' and cover run over `n_sites`. The motivating layout is a site = cell x
#' time-period: plots in a cell-period are the detection replicates, occupancy is
#' per cell-period, and a per-site time weight on a coupled trend field
#' (`icar(graph, weight = time, group_var = "cell")`) gives a detection-corrected
#' occupancy trend on a shared cell field. Supported on the default
#' `joint_coupled` engine; the `v2_joint` / `v3_nested` escape hatches bind the
#' field 1:1 to sites and reject `group_var`.
#'
#' @section Checkpoint / resume:
#' A full-field `occu_cover()` fit integrates over a large outer hyperparameter
#' grid and can run for hours, so a reboot or OOM kill otherwise loses the whole
#' run. `control$checkpoint = list(path = "fit.ckpt", resume = TRUE)` makes the
#' joint engine append each completed grid cell to `path`; a `resume = TRUE` run
#' loads the finished cells and solves only the remaining ones, reproducing the
#' from-scratch fit. `resume = FALSE` starts a fresh file (any stale checkpoint at
#' `path` is removed first). A torn tail from a killed write is discarded and
#' re-solved, and a checkpoint written for different data or settings is rejected
#' rather than resumed onto. Forwarded to [tulpa::tulpa_nested_laplace_joint()].
#'
#' @section Cell-aggregated cover (`cover_aggregate`):
#' On the shared-field spatial path the cover arm has one observation per
#' detected visit, so a cell with many detected plots informs the shared field
#' far more than the single occupancy observation for that cell; the field is
#' then driven almost entirely by the cover arm and the detection-corrected
#' occupancy surface flattens. `cover_aggregate` collapses the cover arm to a
#' single response per occupancy unit (cell-period) so the two arms inform the
#' field with comparable weight:
#' \itemize{
#'   \item `"mean"` -- model the mean cover over the unit's detected visits;
#'   \item `"median"` -- model the median cover;
#'   \item `"latent"` -- model a per-unit cover random effect
#'     \eqn{u_i \sim N(0, \sigma_u^2)} shared across the unit's detected visits
#'     and integrated out per unit, so the cover arm contributes one marginal
#'     observation per unit while keeping every detected visit (the principled
#'     counterpart of mean / median aggregation). The lognormal arm integrates
#'     in closed form (compound-symmetry); the beta arm uses adaptive
#'     Gauss-Hermite quadrature (`control$n.quad`, default 15). The within-unit
#'     dispersion is pre-fit from the within-unit spread and held fixed;
#'     \eqn{\sigma_u} is integrated on the outer grid (`control$sigma.u.grid`).
#'   \item `"none"` -- the per-visit cover arm (one cover observation per
#'     detected visit).
#' }
#' `NULL` (the default) selects `"mean"` on the shared-field spatial
#' `nested_laplace` path and `"none"` on the non-spatial `laplace` path (which
#' has no shared field to over-weight). Aggregation needs a cell-level positive
#' design: a visit-level covariate in the `positive` formula cannot be collapsed
#' to one value per cell and errors. Aggregation is currently wired on the
#' spatial path only; requesting `"mean"` / `"median"` on the non-spatial
#' `laplace` fit errors rather than silently using per-visit cover.
#'
#' @param positive likelihood for the positive cover arm. `"beta"` (cover
#'   in (0, 1)) or `"lognormal"` (log-cover Gaussian).
#' @param cover_aggregate how the cover arm collapses per occupancy unit on the
#'   shared-field spatial path: `"mean"`, `"median"`, `"latent"` (a per-unit
#'   cover random effect integrated out), or `"none"` (per-visit). `NULL`
#'   (default) is `"mean"` on the spatial path and `"none"` on the non-spatial
#'   path. See the *Cell-aggregated cover* section.
#' @return A `tobs_family` object.
#' @seealso [occu()] (no cover), [cover()] (plot-level hurdle, no detection),
#'   [abun()] (counts not cover).
#' @export
occu_cover <- function(positive = c("beta", "lognormal"),
                       cover_aggregate = NULL) {
  positive <- match.arg(positive)
  if (!is.null(cover_aggregate)) {
    cover_aggregate <- match.arg(cover_aggregate,
                                 c("mean", "median", "latent", "none"))
  }
  obs_family(
    name           = "occu_cover",
    class_long     = "joint occupancy-detection + cover hurdle",
    latent         = "bernoulli",
    observation    = if (positive == "beta")
                       "detection_plus_beta"
                     else
                       "detection_plus_lognormal",
    replicates     = "required",
    default_engine = "laplace",
    status         = "experimental",
    params         = list(positive = positive,
                          cover_aggregate = cover_aggregate),
    control_keys   = c(
      "max.iter", "tol", "sigma.beta", "engine",
      "sigma.grid", "alpha.grid", "alpha.grid.trend", "trend",
      "phi.grid.pos", "sigma.u.grid", "n.quad",
      "n.threads", "inner.refresh", "hessian",
      "n.threads.outer", "force.sparse", "integration",
      "adaptive.grid", "adaptive.grid.edge.thresh", "adaptive.grid.max.passes",
      "diagnose.k", "k.samples",
      "checkpoint"
    )
  )
}


#' Community (multispecies) joint occupancy-detection + cover family
#'
#' The community version of [occu_cover()]: a per-species joint occupancy-cover
#' model with Gaussian community hyperpriors on the per-species coefficients of
#' all three arms (occupancy `psi`, detection `p`, and positive cover). Rare
#' species borrow strength from common ones through the shared community means
#' and covariances, the same pooling that stabilises [ms_occu()] and
#' [ms_abun()] but on the joint occupancy + vegetation-cover response.
#'
#' Per species `s`, cell `i`, visit `j`:
#'
#'     z_{s,i}        ~ Bernoulli(psi_{s,i})
#'     y_{s,i,j} | z  ~ Bernoulli(p_{s,i,j})
#'     c_{s,i,j} | y  ~ f_pos(eta_pos_{s,i,j}, disp)
#'     logit psi_{s,i}  = X_occ_i . (mu_occ + b_occ_s)
#'     logit p_{s,i,j}  = X_p_{ij} . (mu_p   + b_p_s)
#'     g(cover)         = X_pos_{ij} . (mu_pos + b_pos_s)
#'     b_occ_s ~ N(0, Sigma_occ), b_p_s ~ N(0, Sigma_p),
#'     b_pos_s ~ N(0, Sigma_pos)
#'
#' The latent presence `z` integrates out per species-cell in closed form (the
#' same two-state mixture as [occu_cover()]); the per-species coefficient
#' deviations are the random effects, integrated by a Laplace-EM. The
#' positive-arm dispersion is a shared community parameter.
#'
#' @section Inputs:
#' `y` and `y_pos` are 3D arrays `[n_sites x max_visits x n_species]` (or named
#' lists of `n_sites x max_visits` matrices, one per species); `species =` is
#' required. The occupancy `formula`, the `detection` formula, and the cover
#' `positive` formula carry community covariates shared across species. `coef()`
#' returns the community means; `ranef()` the per-species coefficient deviations.
#'
#' @section Reduced-rank spatial factors:
#' A single areal field term (`icar(graph = adj)`, `car_proper(graph = adj)`, or
#' `bym2(graph = adj)`) on the occupancy `formula` fits the reduced-rank (HMSC /
#' spatial-gllvm) community model: `K` shared latent fields with per-species
#' loadings on the occupancy predictor,
#' `logit psi_sc = X_c (mu_occ + b_occ_s) + sum_k L_sk w_kc`. A single shared
#' field with only a species intercept can shift each map's level but not its
#' shape; the loadings give each species its own spatial shape as a combination
#' of the shared factors, borrowing strength so rare taxa get a calibrated map.
#' The number of factors is set by `control = list(n.factors = K)`, or chosen
#' automatically by `control = list(n.factors = "auto", n.factors.max = M)`,
#' which fits the identified (lower-triangular loading) ladder and picks the rank
#' that maximises the empirical-Bayes Laplace marginal likelihood (the field is
#' integrated out, so its prior supplies the Occam penalty -- latent-level
#' criteria such as WAIC track the field's effective dimension, not the rank, and
#' under-select). The per-rank evidence is returned in `fit$spatial$K_selection`.
#' `sd.load` sets the loading prior scale.
#'
#' Adding the SAME `icar(graph = adj)` term to the cover `positive` formula
#' shares the latent fields across the two processes: the fields then also load on
#' the cover predictor through a free loading matrix,
#' `g(cover)_sc = X_c (mu_pos + b_pos_s) + sum_k Lpos_sk w_kc`, so a species'
#' spatial occupancy pattern and its spatial cover pattern are linked through one
#' set of factors. The cover loadings are returned in `fit$spatial$loadings_cover`.
#' The field is shared, so the cover-arm term must match the occupancy arm (same
#' type and graph), and a cover-arm field without a matching occupancy field errors.
#'
#' The field structure is set by the term: `icar()` (improper intrinsic CAR, the
#' default), `car_proper(graph = adj)` (a proper CAR whose per-factor correlation
#' `rho` is returned in `fit$spatial$rho_w`), or `bym2(graph = adj)` (the
#' Riebler 2016 convolution whose per-factor spatial-variance fraction `phi` is
#' returned in `fit$spatial$phi_w`); the field hyperparameter is estimated by EM.
#' All three share one engine; `fit$spatial$field_type` records the choice.
#' Structured terms on the detection arm, or an unsupported field term, error from
#' the dispatcher.
#'
#' @section Scope (status `"experimental"`):
#' The non-spatial fit and the occupancy-arm reduced-rank spatial fit are Laplace
#' / Laplace-EM. A shared coupled spatial field across all three arms with the
#' per-species RE block layered on top (the community analogue of
#' [occu_cover()]'s `nested_laplace` joint-coupled engine) needs upstream engine
#' support that does not yet combine a per-group RE block with a shared latent
#' field on one predictor; such a structured term on any arm errors from the
#' dispatcher rather than being silently dropped.
#'
#' @param positive likelihood for the positive cover arm. `"beta"` (cover in
#'   (0, 1)) or `"lognormal"` (log-cover Gaussian).
#' @return A `tobs_family` object.
#' @seealso [occu_cover()] (single species), [ms_occu()] (community occupancy,
#'   no cover), [ms_abun()] (community N-mixture).
#' @export
ms_occu_cover <- function(positive = c("beta", "lognormal")) {
  positive <- match.arg(positive)
  obs_family(
    name           = "ms_occu_cover",
    class_long     = "community joint occupancy-detection + cover hurdle",
    latent         = "bernoulli",
    observation    = if (positive == "beta")
                       "detection_plus_beta"
                     else
                       "detection_plus_lognormal",
    replicates     = "required",
    default_engine = "laplace",
    status         = "experimental",
    params         = list(positive = positive),
    control_keys   = c("max.iter", "tol", "sigma.beta", "newton.max", "sd.load",
                       "n.factors", "n.factors.max")
  )
}


#' Three-level (multiscale) occupancy + cover hurdle family
#'
#' A cell-level occupancy gate, a plot-level availability gate, per-visit
#' detection, and the cover hurdle, for data where a site's "visits" are
#' spatially distinct plots aggregated into a `(cell, period)` rather than
#' temporal revisits (the EVA / MOTIVATE vegetation layout; Nichols et al.
#' 2008; Mordecai et al. 2011). [occu_cover()] treats plots as detection
#' replicates of one occupancy state, which on spatial subunits conflates
#' within-cell prevalence into the detection arm (Kendall & White 2009); this
#' family separates them with an explicit middle level:
#'
#'     z_c        ~ Bernoulli(psi_c)                 # cell / range occupancy
#'     a_cj | z=1 ~ Bernoulli(theta_cj)              # plot availability / use
#'     y_cjv|a=1  ~ Bernoulli(p_cjv)                 # detection
#'     cover|y=1  ~ f_pos(.; eta_pos, disp)          # cover hurdle (beta / lognormal)
#'
#' Both `z` (over cells) and `a` (over plots) marginalize in closed form (two
#' states each), so the joint marginal log-likelihood is exact and reuses the
#' same nested-Laplace cell-coupling machinery as [occu_cover()].
#'
#' @section Inputs:
#' `y` / `y_pos` are `[n_plots x max_visits]` matrices (plots are the rows, the
#' availability units; visits the columns). The state-process `formula` is the
#' cell-level occupancy predictor and MUST carry an areal field naming the
#' per-plot cell column, `icar(graph = adj, group_var = "cell")`. `availability
#' = ~ ...` is the plot-level theta predictor (default `~ 1`); `detection` the
#' per-visit p predictor; `positive = ~ ...` the cover predictor (default the
#' detection formula). `y_pos` is read only where `y == 1`.
#'
#' @section Identifiability:
#' The availability (`theta`) and detection (`p`) levels separate only with
#' replication WITHIN a plot. Single releves supply none, so a plain fit
#' identifies `psi` (cell) and the product `theta * p` (plot) -- it reduces to
#' [occu_cover()] with `p := theta * p`. Within-plot temporal replication (e.g.
#' a resurvey of the same plot in a later period) makes the third level
#' estimable.
#'
#' @section Scope (status `"experimental"`):
#' Spatial joint nested-Laplace only (`method = "nested_laplace"`): a single
#' shared areal field coupled across the occupancy (`sigma`) and cover
#' (`alpha * sigma`) arms, integrated over the outer `(sigma, alpha)` grid.
#' Spatially varying trend fields and a non-spatial Laplace path are not yet
#' wired.
#'
#' @param positive likelihood for the positive cover arm. `"beta"` (cover in
#'   (0, 1)) or `"lognormal"` (log-cover Gaussian).
#' @return A `tobs_family` object.
#' @seealso [occu_cover()] (two-level), [cover()] (plot hurdle, no detection).
#' @export
occu_multiscale_cover <- function(positive = c("beta", "lognormal")) {
  positive <- match.arg(positive)
  obs_family(
    name           = "occu_multiscale_cover",
    class_long     = "three-level occupancy + cover hurdle",
    latent         = "bernoulli",
    observation    = if (positive == "beta")
                       "availability_detection_plus_beta"
                     else
                       "availability_detection_plus_lognormal",
    replicates     = "required",
    default_engine = "nested_laplace",
    status         = "experimental",
    params         = list(positive = positive),
    control_keys   = c(
      "max.iter", "tol", "sigma.beta",
      "sigma.grid", "alpha.grid", "phi.grid.pos", "n.threads",
      "inner.refresh", "hessian", "n.threads.outer", "force.sparse",
      "adaptive.grid", "adaptive.grid.edge.thresh", "adaptive.grid.max.passes",
      "diagnose.k", "k.samples", "checkpoint"
    )
  )
}


# ---------------------------------------------------------------------------
# Planned families — error informatively when used until implemented
# ---------------------------------------------------------------------------

#' N-mixture abundance family
#'
#' Latent Poisson (or NB) abundance with binomial detection per visit
#' (Royle 2004).
#'
#' @param K_max upper bound for the latent-abundance marginal sum (the exact
#'   integration over `N` is truncated at `K_max`). `NULL` (default) lets the
#'   engine pick `max(y) + 100` (matching `unmarked::pcount()`); raise it if a
#'   fit warns that the posterior over `N` puts mass on the boundary.
#' @param mixture latent-abundance distribution, both fitted via tulpa's
#'   closed-form marginal Laplace engine. `"poisson"` is the default;
#'   `"negbin"` (negative binomial, `Var(N) = lambda + lambda^2 / r`) adds an
#'   overdispersion parameter `r`. Non-spatially the log size `log_r` is
#'   estimated jointly with the coefficients and reported with an SE; on the
#'   areal-spatial path (`method = "nested_laplace"`) `r` is integrated over the
#'   outer hyperparameter grid and reported as a posterior mean / sd. Named
#'   `mixture` (after `unmarked::pcount()`) to avoid collision with the
#'   model-type `family` argument of [tobs()].
#' @return A `tobs_family` object.
#' @details
#' Royle's (2004) N-mixture model: latent abundance `N_i ~ Poisson(lambda_i)`
#' with `log lambda_i = X_lambda beta_lambda`, and counts
#' `y_ij | N_i ~ Binomial(N_i, p_ij)` with `logit p_ij = X_p beta_p`. The
#' abundance formula is the `tobs()` `formula`; the per-visit detection formula
#' is `detection`. The marginal likelihood integrates `N` out exactly, so the
#' fit is a direct Laplace approximation (no EM): `method = "laplace"`
#' (fixed effects) or `method = "nested_laplace"` (an areal `icar()` / `bym2()`
#' / `car()` offset on the abundance arm).
#' @export
abun <- function(K_max = NULL, mixture = c("poisson", "negbin")) {
  mixture <- match.arg(mixture)
  obs_family(
    name           = "abun",
    class_long     = "N-mixture abundance",
    latent         = mixture,
    observation    = "binomial_N",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(K_max = K_max, mixture = mixture)
  )
}


#' Multispecies N-mixture family
#'
#' Per-species N-mixture with shared community-level hyperparameters.
#'
#' @inheritParams abun
#' @return A `tobs_family` object.
#' @export
ms_abun <- function(K_max = NULL, mixture = c("poisson", "negbin")) {
  mixture <- match.arg(mixture)
  obs_family(
    name           = "ms_abun",
    class_long     = "multispecies N-mixture",
    latent         = mixture,
    observation    = "binomial_N",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(K_max = K_max, mixture = mixture)
  )
}


#' Open-population (Dail-Madsen) N-mixture family
#'
#' Latent N evolves across seasons via survival + recruitment.
#'
#' @inheritParams abun
#' @return A `tobs_family` object.
#' @export
dyn_abun <- function(K_max = NULL, mixture = c("poisson", "negbin")) {
  mixture <- match.arg(mixture)
  obs_family(
    name           = "dyn_abun",
    class_long     = "Dail-Madsen open N-mixture",
    latent         = "dail_madsen",
    observation    = "binomial_N",
    replicates     = "required",
    default_engine = "nuts",
    status         = "planned",
    params         = list(K_max = K_max, mixture = mixture)
  )
}


#' Distance-sampling family
#'
#' Latent density with hazard-rate or half-normal detection over distance bins.
#'
#' @param key detection-function key. `"halfnorm"` or `"hazard"`.
#' @return A `tobs_family` object.
#' @export
distance <- function(key = c("halfnorm", "hazard")) {
  key <- match.arg(key)
  obs_family(
    name           = "distance",
    class_long     = "distance sampling",
    latent         = "density",
    observation    = "distance_binned",
    replicates     = "optional",
    default_engine = "laplace",
    status         = "planned",
    params         = list(key = key)
  )
}


#' Removal-sampling family
#'
#' Latent N with sequential removal observations.
#'
#' @return A `tobs_family` object.
#' @export
removal <- function() {
  obs_family(
    name           = "removal",
    class_long     = "removal sampling",
    latent         = "poisson",
    observation    = "removal_sequence",
    replicates     = "required",
    default_engine = "laplace",
    status         = "planned"
  )
}


#' False-positive occupancy family
#'
#' Occupancy with both false negatives and false positives in the detection
#' process (e.g. acoustic classifiers, citizen-science misidentification).
#'
#' @return A `tobs_family` object.
#' @export
fp_occu <- function() {
  obs_family(
    name           = "fp_occu",
    class_long     = "false-positive occupancy",
    latent         = "bernoulli",
    observation    = "binomial_with_fp",
    replicates     = "required",
    default_engine = "nuts",
    status         = "planned"
  )
}


#' Cover hurdle family (vegetation cover, MOTIVATE pattern)
#'
#' Latent presence (Bernoulli) plus conditional positive cover (beta or
#' lognormal). Does not share the replicate-detection assumption of the other
#' families — see `vignette("families")` for the conceptual caveat.
#'
#' @section Joint nested-Laplace engine — spatial-prior parameterisation:
#' When fitted with `method = "nested_laplace"` and an areal spatial term in
#' the latent-presence formula (`bym2(graph = adj)`, or `car()` /
#' `car_proper()`), the engine identifies a single latent field `z` per region
#' and parameterises the two arms as `eta_occ = X beta_occ + sigma * z`
#' and `eta_pos = X beta_pos + alpha * sigma * z`. Both BYM2 sub-blocks
#' (`phi`, `theta`) are subject to hard sum-to-zero constraints (see
#' `.joint_inner_var()` for the math): the cover-arm intercept
#' `beta_pos[1]` is identified as the population mean of `eta_pos`
#' under `mean(z) = 0`. Any simulator generating data for this engine
#' must demean both `phi_f` and `theta_f` before scaling, otherwise the
#' estimator targets `beta_pos_0_truth + alpha * mean(w_s_sim)` and
#' coverage of the *population* truth collapses with alpha. See
#' [simulate_cover_joint()] for a ready-made demeaned simulator.
#'
#' @section Checkpoint / resume:
#' A spatial cover-hurdle fit integrates over a large outer hyperparameter grid
#' and can run for hours. `control$checkpoint = list(path = "fit.ckpt", resume =
#' TRUE)` makes the joint engine append each completed grid cell to `path`; a
#' `resume = TRUE` run loads the finished cells and solves only the rest,
#' reproducing the from-scratch fit, so a killed or rebooted fit resumes instead
#' of restarting. `resume = FALSE` starts a fresh file. Forwarded to
#' [tulpa::tulpa_nested_laplace_joint()].
#'
#' @param positive likelihood for the positive part. `"beta"` (cover in
#'   (0, 1)) or `"lognormal"` (log-cover Gaussian).
#' @return A `tobs_family` object.
#' @export
cover <- function(positive = c("beta", "lognormal")) {
  positive <- match.arg(positive)
  obs_family(
    name           = "cover",
    class_long     = "vegetation cover hurdle",
    latent         = "hurdle",
    observation    = if (positive == "beta") "binomial_plus_beta" else "binomial_plus_lognormal",
    replicates     = "single",
    default_engine = "laplace",
    status         = "working",
    params         = list(positive = positive),
    # The cover hurdle has its own (.dispatch_cover) grid-based control surface,
    # separate from the occupancy fitter and named with underscores. Declaring
    # the keys keeps tobs()'s control validation from rejecting them.
    control_keys   = c(
      "max.iter", "tol", "n.threads", "n.threads.outer", "prior.sigma", "prior.alpha",
      "phi.grid", "sigma.grid", "sigma.pos.grid", "rho.grid", "tau.grid",
      "rho.car.grid", "tau.temporal.grid", "rho.temporal.grid",
      "sigma.temporal.grid", "sigma.re.grid",
      "trend", "alpha.grid", "alpha.grid.trend", "integration",
      "adaptive.grid", "adaptive.grid.edge.thresh", "adaptive.grid.max.passes",
      "prune", "prune.tol", "hessian", "aggregate.occ",
      "progress", "progress.every", "progress.throttle", "progress.file",
      "checkpoint"
    )
  )
}


# ---------------------------------------------------------------------------
# Print method
# ---------------------------------------------------------------------------

#' Print method for tobs_family
#' @param x a `tobs_family` object.
#' @param ... ignored.
#' @return `x`, invisibly.
#' @export
print.tobs_family <- function(x, ...) {
  cat(sprintf("<tobs_family: %s>\n", x$class_long))
  cat(sprintf("  name           : %s\n", x$name))
  cat(sprintf("  latent state   : %s\n", x$latent))
  cat(sprintf("  observation    : %s\n", x$observation))
  cat(sprintf("  replicates     : %s\n", x$replicates))
  cat(sprintf("  default method : %s  (method = \"auto\")\n", x$default_engine))
  cat(sprintf("  status         : %s\n", x$status))
  if (length(x$params)) {
    cat("  params         :\n")
    for (nm in names(x$params)) {
      cat(sprintf("    %s = %s\n", nm, format(x$params[[nm]])))
    }
  }
  invisible(x)
}
