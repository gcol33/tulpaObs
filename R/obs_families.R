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
#' @param default_engine `"laplace"`, `"nested_laplace"`, `"nuts"`, or
#'   `"pg_gibbs"`. This is the route `method = "auto"` resolves to.
#' @param status `"working"`, `"planned"`, or `"experimental"`.
#' @param response shape of the family's response: `"vector"` for a plain
#'   length-N response vector (the cover hurdle), or `"matrix"` (the default)
#'   for a detection-history / count matrix, 3D array, or list. A
#'   single-vector-response family accepts the response on the top formula
#'   left-hand side (`response ~ predictors`) so `y =` may be omitted; a
#'   matrix-response family takes the response via `y =` only (a matrix does
#'   not sit on a formula LHS). Consulted by [tobs()].
#' @param params named list of family-specific parameters carried with the
#'   object (K_max, positive-part link, etc.).
#' @param control_keys character vector of extra `control` names this family's
#'   dispatcher accepts beyond the engine/route controls. These are added to
#'   the allowlist `tobs()` validates `control` against, so a family with a
#'   bespoke dispatcher (e.g. the cover hurdle's grid controls) is not rejected.
#'   Admitted on every route the family supports.
#' @param control_groups character vector of capability groups (names of
#'   `.tobs_control_groups`) this family participates in beyond the ones its
#'   route admits unconditionally. Unlike `control_keys` these stay route-gated,
#'   so a group tied to the Laplace engines is still rejected under `"nuts"`.
#'   `"block_coordinate"` is the case this exists for: `max.outer` /
#'   `factor.starts` mean something only to a family whose latent structure is fit
#'   by the block-coordinate driver.
#'
#' @return A `tobs_family` object.
#' @keywords internal
#' @export
obs_family <- function(name,
                       class_long,
                       latent,
                       observation,
                       replicates    = c("required", "optional", "single"),
                       default_engine = c("laplace", "nested_laplace", "nuts",
                                          "pg_gibbs"),
                       status         = c("working", "planned", "experimental"),
                       params         = list(),
                       control_keys   = character(0),
                       control_groups = character(0),
                       response       = c("matrix", "vector")) {
  replicates     <- match.arg(replicates)
  default_engine <- match.arg(default_engine)
  status         <- match.arg(status)
  response       <- match.arg(response)

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
      control_keys   = control_keys,
      control_groups = control_groups,
      response       = response
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
#' seasons (the MacKenzie et al. dynamic model). Colonisation / extinction are
#' constant across a site's seasons by default; a season-varying rate is supplied
#' as a `[n_sites x (T - 1)]` matrix column of `data` (one column per transition
#' interval) on `colonization` / `extinction`. Detection is site-level by default;
#' a `[n_sites x T]` matrix column of `data` (one column per primary season) on
#' the `detection` formula gives per-season detection. All season-varying paths
#' run under `method = "laplace"` (the exact HMM-forward marginal refine
#' calibrates the coefficients); they are gated under `"nuts"` with a pointer.
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
#' A continuous-mesh `spde()` term on the `detection` formula loads a
#' spatially-structured detection field, broadcast to every source, under
#' `method = "laplace"`. The per-source fields are reported as a named list in
#' `fit$spatial_field_det`. All sources must share one detection formula; areal
#' terms belong on the state arm.
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
    status         = "working",
    # Shares the ms_count binder and fitter, so latent() factors and a shared
    # field are fit by the same block-coordinate driver.
    control_groups = c("block_coordinate", "block_coordinate_factor")
  )
}


#' Count / relative-abundance GLMM family (no detection process)
#'
#' A generalized linear model on an observed count (or continuous) response
#' directly, with no detection sub-model and no latent abundance to marginalize.
#' The abundance analogue of [jsdm()]: where [abun()] fits Royle's N-mixture
#' (latent `N` plus imperfect detection), `count()` fits the observed response
#' as a plain GLMM. This is the "relative abundance" model of `spAbundance`
#' (`abund`): one value per site, `log`-link Poisson / negative-binomial counts
#' or an identity-link Gaussian response.
#'
#' The response is a numeric vector (one value per site), supplied via `y =` or
#' on a two-sided `formula` left-hand side (`count.value ~ predictors`).
#'
#' @param response The response distribution: `"poisson"` (log link),
#'   `"negbin"` (negative binomial, log link, an estimated size / dispersion),
#'   `"gaussian"` (identity link, an estimated residual variance), or
#'   `"binomial"` (logit link, `k` successes out of `n` trials per site). The
#'   binomial response is the detection-free binomial GLMM of `spOccupancy`
#'   (`svcPGBinom`): supply the per-site trial count as `trials =` on [tobs()]
#'   (default 1, i.e. a Bernoulli response). Unlike the count / Gaussian areal
#'   fit, a binomial areal field is identified (the variance is pinned by `n`),
#'   so an `icar()` / `car_proper()` term composes with any `trials`.
#' @return A `tobs_family` object.
#' @seealso [abun()] (N-mixture: latent abundance + detection), [jsdm()]
#'   (occurrence GLMM, no detection).
#' @export
#' @examples
#' f <- count("poisson")
#' f
count <- function(response = c("poisson", "negbin", "gaussian", "binomial")) {
  response <- match.arg(response)
  obs_family(
    name           = "count",
    class_long     = "count / relative-abundance GLMM",
    latent         = "none",
    observation    = response,
    replicates     = "single",
    default_engine = "laplace",
    status         = "working",
    params         = list(response = response),
    response       = "vector"
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
#' A detected visit (`y_ij = 1`) with a missing cover (`y_pos_ij = NA`) keeps its
#' detection term but drops the `f_pos` factor: cover is taken missing-at-random,
#' so the cover likelihood runs over the detected visits with an observed cover.
#'
#' Reduces to [occu()] when the cover arm is degenerate, and to the
#' plot-level cover hurdle ([cover()]) when J = 1 and detection is perfect.
#'
#' @section Engines:
#' The non-spatial fit is a direct Laplace approximation (`method = "laplace"`)
#' or a NUTS sampler over the same exact two-state marginal (`method = "nuts"`,
#' beta or lognormal cover), giving calibrated intervals and a per-draw
#' pointwise likelihood for WAIC / LOO. A shared areal field across the
#' occupancy and cover arms is the `method = "nested_laplace"` path (a structured
#' `icar()` / `bym2()` term on the psi formula); a structured term on a
#' non-spatial route errors from the dispatcher with a pointer to it.
#'
#' `method = "nuts"` also samples the coupled areal field(s)
#' (`icar()` / `bym2()` / `car_proper()`), written either as areal terms or as
#' the bar `spatial(~ 1 || cell, graph = adj)`, TOGETHER WITH their
#' hyperparameters: the field SD, the mixing (`bym2`) or spatial-correlation
#' (`car_proper`) parameter, and - where the formula asks for one - the cover-arm
#' copy amplitude. Each is bounded to the span the `nested_laplace` path's
#' outer quadrature integrates that axis over, in that grid's own coordinate,
#' and the same `control$sigma.grid` / `alpha.grid` / `rho.car.grid` knobs set
#' it on both routes - so the sampler is an independent reference for the
#' hyperparameter layer rather than a fit conditioned on that layer's point
#' estimate. The measure over that span is the one the outer grid declares for
#' the axis: flat in the grid's coordinate for the field SD and the mixing
#' parameter, and for the copy amplitude the penalized-complexity Exponential
#' slab, with `control$copy.slab = "flat"` asking for the flat alternative the
#' grid also accepts. The grid weighs that slab against a point mass at
#' `alpha = 0` ("no coupling"), which a gradient sampler cannot visit, so the
#' sampled copy amplitude is the posterior conditional on a coupled field; the
#' point mass's share of the prior is read on the `nested_laplace` route. `fit$nuts$sampled_hyper` and
#' `fit$nuts$fixed_hyper` name which is which per fit (`icar` pins `rho` at 1:
#' the intrinsic precision has no mixing parameter; a field with no `copy()` pins
#' the copy amplitude at 0), `fit$hyper_draws` carries
#' their posterior alongside `field_sd`, the geometric-mean marginal SD the
#' field implies. `control$fixed.hyper = TRUE` conditions on the warm
#' nested-Laplace estimate instead.
#'
#' The grid-integrated `nested_laplace` path reports `field_sd` too
#' (`fit$spatial$field_sd_mean` / `field_sd_sd`, alongside `sigma_mean`), in
#' the SAME geo-mean-marginal-SD convention as the NUTS path's `field_sd`.
#' It is a fixed multiple of `sigma` (not an independent grid axis), so it is
#' NOT in `fit$means` / `fit$vcov` -- folding it in there would make the
#' joint parameter vcov exactly singular. Do not compare fits, or a fit
#' against [simulate_occu_cover()]'s `sigma`, by reading the raw `sigma` on
#' this path: it is the field's amplitude against the unscaled intrinsic ICAR
#' precision `Q = D - W`, which differs from `field_sd`
#' by `sqrt(scale_q)`, a graph-size-dependent factor (about 2.1 for a 30-node
#' chain graph). `field_sd` is the number comparable across fits and to a
#' simulation truth; `sigma` is the raw amplitude the engine's grid axis is
#' spelled in.
#'
#' A spatially-varying coefficient - a *weighted* areal term beside the
#' unweighted intercept field, or a bar column
#' `spatial(~ 1 + year || cell, graph = adj)` - is sampled as a SECOND field
#' block, with its own whitened surface and its own (SD, mixing, copy amplitude)
#' coordinates: two fields share no hyperparameter. The surfaces are
#' `fit$spatial_field` and `fit$trend_field` / `fit$trend_fields` (named by the
#' weight column), and each block's hypers carry that block's suffix in
#' `fit$hyper_draws` and in `fit$nuts$sampled_hyper` (`sigma_trend`,
#' `alpha_trend`, indexed when there are several). Each field's amplitude is
#' addressed on its own through `copy(spatial(), terms = list(...))` -- with
#' stated nodes, or `grid(n = )` for the engine's axis read more finely -- or
#' through `control$alpha.grid[.trend]` / `control$alpha.n[.trend]`, the same
#' knobs the grid-integrated route reads. A second field
#' multiplies the warm fit's outer grid, so a defaulted axis is thinned to three
#' nodes over the same span - the sampled prior is unchanged, the warm fit stays
#' affordable. A correlated bar (`|`, one free-Sigma MCAR block across the
#' fields) stays on the grid-integrated `nested_laplace` path; the community
#' spatial-factor route is [ms_occu_cover()].
#'
#' On the shared-field `nested_laplace` path the field-coupled occupancy slope
#' Wald interval is mildly anti-conservative at small N (the grid-integrated
#' Laplace under-disperses that coefficient; the pooled coverage across the three
#' arms stays near nominal). The non-spatial `nuts` path gives fully calibrated
#' occupancy intervals.
#'
#' @section What a bare areal term means for the cover arm:
#' A term's process is the formula it sits in. An areal term written in the
#' occurrence formula therefore loads on the OCCURRENCE arm alone: the cover arm
#' sees no field, and the copy amplitude `alpha` is 0. This holds under
#' `method = "nested_laplace"` and `method = "nuts"` alike, so the same input
#' fits the same model on either engine.
#'
#' To carry the occurrence field onto the cover arm as well, write the copy in
#' the `positive` formula:
#'
#' ```r
#' positive = ~ pos_cov + copy(spatial())                      # amplitude estimated
#' positive = ~ pos_cov + copy(spatial(), alpha = grid(c(...))) # over given nodes
#' positive = ~ pos_cov + copy(spatial(), alpha = grid(n = 9))  # engine axis, finer
#' positive = ~ pos_cov + copy(spatial(), alpha = 0.5)          # fixed
#' positive = ~ pos_cov + copy(spatial(), prior = list("pc.prec", c(4, 0.01)))
#' ```
#'
#' The amplitude is integrated over that axis on the `nested_laplace` path and
#' sampled over its span on the `nuts` path; a scalar `alpha =` pins it on both.
#' `prior =` regularizes the coefficient itself, in the joint driver's
#' `list(<family>, <params>)` shape (see `prior_alpha` in
#' [tulpa::tulpa_nested_laplace_joint()]); one reaches the engine per fit, so a
#' fit copying both an intercept and a trend block is refused it rather than
#' given it on the first (gcol33/tulpa#655).
#'
#' `control$alpha.grid[.trend]`, `control$alpha.n[.trend]` and
#' `control$prior.alpha` are the lower-level spelling of the same three
#' requests, and the representation the formula compiles into; set each in one
#' place, not both.
#'
#' Without a `copy()` (and without `control$alpha.grid`) the occurrence field is
#' NOT carried onto the cover arm -- there is no implicit coupling -- so the
#' amplitude is not a parameter of the model, and neither engine reports one.
#' Under `nested_laplace` `alpha` is absent from `coef()`, `vcov()` and the
#' `n_params` count; under `nuts` it is absent from `fit$nuts$sampled_hyper`,
#' `fit$nuts$fixed_hyper` and the `fit$hyper_draws` columns. `predict()` returns
#' a cover arm the occurrence field does not enter. Add `copy(spatial())` to
#' couple the two arms.
#'
#' The axis the amplitude rides carries prior structure -- a point mass at
#' `alpha = 0` ("no coupling") and a log-spaced slab above it -- so it is set
#' in one of two ways. `copy(alpha = grid(...))` STATES its nodes, and with them
#' that structure. `copy(alpha = grid(n = 9))` states a RESOLUTION: the engine
#' re-reads its own axis with that many slab nodes, point mass and slab bounds
#' unchanged, so sharpening the axis never restates its structure. `terms =`
#' gives either, per block. The resolution is the only way to raise this axis,
#' because it does not densify when the donor `control$sigma.grid` does: on an
#' informative data set the outer grid's quadrature effective sample size
#' saturates on the copy amplitude while every other axis tracks the request
#' (measured engine-side, `NOTES_measurements.md`). One block takes one of the
#' two; giving both is an error. `control$alpha.grid[.trend]` and
#' `control$alpha.n[.trend]` are the lower-level spelling of the same pair.
#'
#' Set it when the fit's hyperparameter intervals are reported. On the coupled
#' SBC fixture the declared resolution leaves the field SD miscalibrated once
#' the data are informative -- at 10 visits per site its uniformity p-value is
#' 9.1e-05, against 0.17 at 3 visits -- and `control$alpha.n = 21` returns every
#' scored quantity to nominal at both. Thirteen slab nodes do not (the field SD
#' still reads 4.5e-03), and the price is a 2.3-3x longer fit, since the outer
#' grid is a tensor. Point estimates are not affected the way the intervals are;
#' the measurement is in `NOTES_measurements.md`.
#'
#' @section Coupled fields and spatially-varying trends:
#' The spatial engines (`method = "nested_laplace"`, default
#' `control$engine = "joint"`, and `method = "nuts"`) share one areal field (the
#' cell intercept, an unweighted `icar()` / `bym2()` term) across the occupancy
#' and cover arms, when the `positive` formula copies it. ADDITIONAL coupled
#' fields - spatially-varying coefficients,
#' e.g. a temporal trend - are added as *weighted* areal terms in the formula:
#'
#' ```r
#' ~ elev + icar(graph = adj) + icar(graph = adj, weight = year)
#' ```
#'
#' Each weighted term `icar(graph, weight = col)` is a second field whose
#' contribution to a predictor row is `weight_i * z[cell_i]`. All such fields
#' share the same graph; a `copy()` carries one onto the cover arm with its own
#' scale (`alpha` for the intercept field, `alpha_trend` for a trend field),
#' integrated over the outer grid, and a field the `positive` formula does not
#' copy stays on occurrence. `copy(spatial(), terms = list(...))` gives a
#' per-block amplitude and must address every block. The intercept field is
#' reported in `fit$spatial_field`, trend fields in `fit$trend_field` /
#' `fit$trend_fields`. The trend coupling grid defaults to
#' `control$alpha.grid`; override it with `control$alpha.grid.trend`. Under
#' `method = "nuts"` each field is its own sampled block, with the scales
#' sampled over those same spans rather than integrated over the grid.
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
#' `joint` engine; the `v2_joint` / `v3_nested` escape hatches bind the
#' field 1:1 to sites and reject `group_var`.
#'
#' @section Per-group random intercept on the shared-field path (`re()` / `(1 | g)`):
#' On the spatial `nested_laplace` path a single random INTERCEPT term on the psi
#' formula -- `re(g)` or the `(1 | g)` bar -- layers a per-group occupancy offset
#' on top of the shared field. It joins the joint fit as one `iid` prior block
#' whose latent rides the occupancy arm only; its variance integrates on the
#' outer grid alongside the field sigma / alpha and is reported as the `sigma_re`
#' hyperparameter, with the per-group BLUPs in `fit$re` and via [ranef()]. A
#' default `re.sigma.grid` (log-spaced) sets the grid; pass `control$re.sigma.grid`
#' to override. Scope: one random-intercept term (a scalar per group maps onto the
#' one iid block); a random slope, a correlated multi-coefficient block, or an RE
#' without a shared field is rejected (the joint engine integrates every variance
#' component on its grid, so multiple variance components do not scale -- richer
#' RE is the non-spatial cover-hurdle EM's route). This is the single-species
#' analogue of the community spatial occu_cover ([ms_occu_cover()])
#' and the tulpaObs consumer of tulpa's field + per-group RE engine composition.
#'
#' A bar is therefore always a random effect, never a spatial field, even when the
#' grouping factor names the areal graph nodes. The engine's inline-MCAR call
#' `tulpa::spatial(graph, ~ 1 + x | cell)` reads `| cell` as a separable spatial
#' field, but the same spelling in an `occu_cover()` formula is an IID random
#' effect on `cell`. For a spatial field on the cells use an areal term --
#' `icar(graph = adj, group_var = "cell")` (plus
#' `icar(graph = adj, weight = x, group_var = "cell")` for a spatially-varying
#' trend) -- not a bar. When a formula carries a bar whose grouping factor is also
#' an areal term's `group_var`, `occu_cover()` emits a one-time message noting the
#' bar is fitted as a random effect; suppress it with [base::suppressMessages()].
#'
#' @section Independent field on the cover arm (placement):
#' A `copy(spatial())` in the `positive` formula shares the occupancy field: the
#' cover arm sees it as `alpha * (occupancy field)`, the coregionalization copy
#' on the outer `(sigma, alpha)` grid. When the cover trend is spatially
#' structured but is not
#' a scalar multiple of the occupancy field, `alpha` collapses toward 0 and the
#' cover arm inherits no field, so per-cell conditional cover (and its change over
#' time, `delta_cover_cond`) comes out flat. A `spatial()` bar placed in the
#' `positive` formula adds an INDEPENDENT, non-copied areal field on the cover arm
#' alone:
#'
#' ```r
#' tobs(occurrence = ~ x + icar(graph = adj, group_var = "cell"),
#'      detection  = ~ 1,
#'      positive   = ~ time + spatial(~ 1 + time || cell, graph = adj),
#'      data = cell_dat, y = y, y_pos = y_pos,
#'      family = occu_cover("lognormal"), method = "nested_laplace")
#' ```
#'
#' The occupancy field still drives psi and, via the alpha copy, `delta_cover_exp`;
#' the independent cover field carries a cover-specific structure the alpha copy
#' cannot express. Its intercept + coefficient columns become separate ICAR blocks
#' on the cover arm, each with its own precision integrated on the outer grid and
#' reported as `sigma_pos_field` (intercept) / `sigma_pos_field_<col>` (a covariate
#' column). The field SD grid defaults to `control$sigma.grid`; set
#' `control$sigma.grid.pos.field` to integrate it over its own (usually coarser)
#' grid, which keeps the added axis from multiplying the outer-grid cost.
#' The per-cell field posterior is in `fit$pos_field` / `fit$pos_field_table`
#' (and `fit$pos_field_tables` per column). Scope: it composes with the shared
#' occupancy field but uses per-visit cover (`cover_aggregate = "none"`) and does
#' not combine with the correlated `|` MCAR field, the latent cover RE, or the
#' batched fused path; like the shared field it is fitted as ICAR (bym2/car read as
#' ICAR). A static intercept cover field is only weakly identified against the
#' occupancy field's alpha copy (both are per-cell cover offsets); a time-weighted
#' trend cover field is identified separately and is what makes `delta_cover_cond`
#' spatially varying.
#'
#' The same placement works on the detection arm: a spatial-field term in the
#' `detection` formula (`detection = ~ 1 + spatial(~ 0 + time || cell, graph =
#' adj)`) fits an independent field on the detection predictor, for a
#' spatially-structured detection probability. Its SD is reported as
#' `sigma_p_field` / `sigma_p_field_<col>`.
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
#' @section Random effects on the observation arms:
#' Under `method = "nested_laplace"` (shared-field spatial) and `method = "nuts"`
#' (non-spatial, which samples the group SD), random
#' intercepts may be written on the **detection** or **positive-cover** formula
#' with the usual `lme4` bar or `re()` spelling, e.g.
#' `detection = ~ effort + (1 | habitat)`. The grouping is per visit -- one code
#' per `(site, visit)` entry, distinct from the per-site occupancy-arm
#' grouping -- so a many-level categorical visit covariate (an EUNIS
#' habitat class, an observer) enters as a partially pooled random intercept.
#' **Crossed** (`(1 | habitat) + (1 | observer)`) and **nested**
#' (`(1 | region/site)`, which expands to `re(region) + re(region:site)`)
#' groupings are supported: each term joins the joint fit as its own `iid` latent
#' block. Each block's variance integrates on the outer grid (`sigma_re_p` /
#' `sigma_re_pos` for a lone term on the detection / cover arm, suffixed by the
#' grouping variable -- `sigma_re_p_habitat` -- when several terms share an arm;
#' tune with `control$re.sigma.grid.p` / `re.sigma.grid.pos`). The per-group BLUPs
#' are reported in `fit$re` (one entry per term, keyed by arm or `"<arm>:<var>"`)
#' and via [ranef()]. `predict()` sums every term's BLUP offset on the predicted
#' arm when `newdata` carries the grouping columns and `type = "detection"` (or
#' `"occurrence"` / `"cover_cond"`); an unseen / held-out level shrinks that term
#' to the population mean. Several crossed terms multiply the outer grid, so set
#' `control$integration = "ccd"` (and/or coarsen the RE grids) at scale.
#'
#' Random **slopes** are also supported (tulpa engine >= 0.0.39): an uncorrelated
#' slope (`(x || g)`, `(0 + x | g)`) is one per-row weighted `iid` block per
#' coefficient, and a correlated slope (`(1 + x | g)`) one multivariate
#' free-Sigma `miid` block. A slope term's `fit$re` entry carries an
#' `[n_groups x n_coefs]` BLUP matrix, a per-coefficient `sigma`, and (correlated)
#' a `cor` matrix; the hyperparameters are `sigma_re_p_<coef>` and
#' `cor_re_p_<ci>_<cj>`, marginalized over the grid. `predict()` weights each
#' coefficient by its covariate column in `newdata` (intercept = 1). Each slope
#' covariate is **standardized** to unit SD before fitting (mirroring the
#' fixed-effect design), so the variance grid is scale-invariant -- the reported
#' slope BLUPs and `sigma` are back-transformed to the covariate's natural units
#' (correlation is scale-free). A correlated slope adds `p(p+1)/2` Sigma axes to
#' the outer grid, so its free-Sigma grid uses a compact principled default
#' (symmetric correlation nodes including 0); widen it with
#' `control$re.logchol.grid.p` / `re.logchol.grid.pos`.
#'
#' The positive-cover RE needs per-visit cover (`cover_aggregate = "none"`). As
#' with the occupancy-arm RE, each grid-integrated variance carries the binary /
#' small-cluster inner-Laplace attenuation and is a lower bound on the truth; the
#' BLUPs recover the per-group structure.
#'
#' `method = "nuts"` samples an observation-arm random **intercept** instead of
#' integrating it on a grid: each grouping factor is a non-centered block
#' (`b_g = sigma_re z_g`, `z_g ~ N(0, 1)`) whose group SD is a coordinate of the
#' sampled vector under a `N(0, 1.5^2)` prior on `log sigma_re` -- the same
#' weakly informative width the other observation-family samplers use, and the
#' reason the sampled SD carries no small-cluster attenuation to correct for.
#' Crossed and nested groupings are simply several blocks, on either or both
#' observation arms. `fit$re` keys and [ranef()] match the `nested_laplace`
#' path; each term additionally reports `sigma_median` and the per-draw
#' `sigma_draws`, the summaries to quote for a right-skewed variance component,
#' with per-SD convergence in `fit$nuts$re_sigma_rhat` / `re_sigma_ess`. The
#' sampled coefficient surface (`coef()`, `vcov()`, `fit$draws`) is unchanged.
#' Still on `nested_laplace` only: a random slope, and an observation-arm RE
#' composed with the coupled areal field. `waic()` / `cpo()` on a
#' sampled fit read the coefficient draws, so they score the RE at zero -- the
#' same limitation the sampled coupled field has.
#'
#' @param response likelihood for the positive cover arm. `"beta"` (cover
#'   in (0, 1)), `"lognormal"` (log-cover Gaussian), or `"gaussian"` (an
#'   identity-Gaussian magnitude on a real, unbounded scale -- the delta-normal
#'   hurdle; not for cover fractions, which stay on `"beta"`/`"lognormal"`).
#' @param cover_aggregate how the cover arm collapses per occupancy unit on the
#'   shared-field spatial path: `"mean"`, `"median"`, `"latent"` (a per-unit
#'   cover random effect integrated out), or `"none"` (per-visit). `NULL`
#'   (default) is `"mean"` on the spatial path and `"none"` on the non-spatial
#'   path. See the *Cell-aggregated cover* section.
#' @return A `tobs_family` object.
#' @seealso [occu()] (no cover), [cover()] (plot-level hurdle, no detection),
#'   [abun()] (counts not cover).
#' @examples
#' \donttest{
#' N <- 120; J <- 4
#' sim <- simulate_occu_cover(N = N, J = J, n_occ_covs = 1, n_det_covs = 1,
#'                            n_pos_covs = 1, positive = "beta", seed = 1)
#' long <- data.frame(site_id = rep(seq_len(N), each = J),
#'                    visit    = rep(seq_len(J), times = N),
#'                    y        = as.vector(t(sim$y)),
#'                    det_cov1 = sim$visit_data$det_cov1,
#'                    pos_cov1 = sim$visit_data$pos_cov1)
#' od  <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
#'                  det.covs = c("det_cov1", "pos_cov1"))
#' cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
#' y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
#' fit <- tobs(~ occ_cov1, data = cell_dat, family = occu_cover("beta"),
#'             detection = ~ det_cov1, positive = ~ pos_cov1,
#'             y = od$y, y_pos = y_pos, visits = od$det.covs, method = "laplace")
#' summary(fit)
#' }
#' @export
occu_cover <- function(response = c("beta", "lognormal", "gaussian"),
                       cover_aggregate = NULL) {
  positive <- match.arg(response)
  if (!is.null(cover_aggregate)) {
    cover_aggregate <- match.arg(cover_aggregate,
                                 c("mean", "median", "latent", "none"))
  }
  obs_family(
    name           = "occu_cover",
    class_long     = "joint occupancy-detection + cover hurdle",
    latent         = "bernoulli",
    observation    = switch(positive,
                            beta      = "detection_plus_beta",
                            gaussian  = "detection_plus_gaussian",
                            "detection_plus_lognormal"),
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(positive = positive,
                          cover_aggregate = cover_aggregate),
    control_keys   = c(
      "max.iter", "tol", "sigma.beta", "engine",
      "sigma.grid", "alpha.grid", "alpha.grid.trend", "trend",
      # Resolution of the copy coefficient's own outer axis: the engine re-reads
      # its declared axis -- the atom at alpha = 0 plus the log-spaced slab --
      # with this many slab nodes. `alpha.grid` states nodes instead, so a block
      # takes one of the two.
      "alpha.n", "alpha.n.trend",
      "rho.car.grid",
      "phi.grid.pos", "sigma.grid.pos.field", "sigma.u.grid", "n.quad",
      "n.threads", "inner.refresh", "hessian",
      "n.threads.outer", "force.sparse", "integration",
      "adaptive.grid", "adaptive.grid.edge.thresh", "adaptive.grid.max.passes",
      # Tuning for `integration = "grid_adaptive"`, a DIFFERENT mechanism from
      # the three refinement knobs above: the keep / expand radius from the
      # peak, the coarse-seed subsample stride per axis, the kept-fraction
      # ceiling past which the builder declines back to the dense tensor, and
      # the smallest dense tensor worth the machinery at all.
      "adaptive.grid.cutoff", "adaptive.grid.stride",
      "adaptive.grid.max.frac", "adaptive.grid.min.cells",
      "var.of.means.consistency", "var.of.means.min.ess",
      # Cheap-pass outer-grid screen (see the engine's `prune` / `prune_tol`),
      # and the two placement knobs: whether a default axis is recentred on the
      # hyperparameter mode at all, and whether that placement is detected on a
      # thinned pilot grid instead of a full extra solve.
      "prune", "prune.tol", "auto.recenter", "recenter.pilot",
      # Shape of the prior on the cross-arm copy scale: the continuum measure
      # ("exponential" or "flat" in log alpha over the `alpha.grid` span) and
      # the prior probability of the no-coupling point mass at alpha = 0.
      # `copy.slab` declares the same continuum measure on the NUTS + areal
      # sampler's own copy coordinate, where it defaults to "flat"; that
      # sampler has no way to visit the point mass, so `copy.atom.mass` reaches
      # only its warm nested-Laplace fit.
      "copy.slab", "copy.atom.mass",
      "diagnose.k", "diagnose.draws", "k.samples", "k.bootstrap",
      "k.tail.points", "k.conf.bands",
      "re.sigma.grid", "re.sigma.grid.p", "re.sigma.grid.pos",
      "re.logchol.grid.p", "re.logchol.grid.pos",
      # Condition the NUTS + areal sampler on the warm nested-Laplace fit's
      # (sigma, rho, alpha) instead of sampling them.
      "fixed.hyper",
      "checkpoint",
      # Diagnostic-parallelism thread count for the outer-grid Pareto-k pass
      # (control$diagnose.k), separate from n.threads.outer's inner-Newton use.
      "k.threads"
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
#' The shared fields imply a residual species-association matrix (the
#' spatial-JSDM / HMSC output): [tobs_associations()] returns the occupancy
#' association `corr(L L')` and, with a cover-arm factor, the cover association
#' and the cross-arm occupancy-vs-cover association, each marginalised over the
#' loading posterior so an interval accompanies the estimate. `predict()` returns
#' the other JSDM output -- the per-species per-cell occupancy maps (`psi` with
#' `psi_lower` / `psi_upper`), marginalised over the loading + field posterior so
#' a rare species borrows strength across the shared factors for a calibrated map.
#'
#' @section Scope:
#' The non-spatial fit is Laplace-EM. A community spatial occu_cover -- a shared
#' latent field coupled across the occupancy and cover arms with per-species RE on
#' all three arms -- is the reduced-rank spatial-factor fit: an `icar()` /
#' `car_proper()` / `bym2()` term on the occupancy formula (and, with a matching
#' term on the cover formula, a cover-arm factor) fits per-species loadings on the
#' shared field via Laplace-EM and NUTS, with the per-species community covariances
#' `Sigma_occ` / `Sigma_p` / `Sigma_pos` on top. The free
#' per-species loading form generalises a single common field amplitude, so it
#' subsumes the common-amplitude coupling of [occu_cover()]'s joint-coupled engine.
#' A community model whose per-arm variance components are integrated on an outer
#' grid (the engine route of [occu_cover()]'s joint-coupled path) is not used: the
#' joint nested-Laplace engine integrates every variance component on its outer
#' grid, so per-arm community RE variances plus the field hyperparameters exceed
#' the engine's grid cap -- the closed-form covariance M-step of the Laplace-EM is
#' the scaling route for community variance components.
#' Structured terms beyond the shared field error from the dispatcher rather than
#' being silently dropped.
#'
#' @section Community variance debias:
#' The community-MEAN estimates (`coef()`, `vcov()`, `confint()`) are unbiased.
#' The community-VARIANCE components -- the per-arm covariance matrices in
#' `fit$ms_community$Sigma_occ` / `Sigma_p` / `Sigma_pos` and their `sd_*`, which
#' set the spread of the per-species deviations -- pick up Laplace small-cluster
#' attenuation at small per-species n under the EM M-step. They are debiased by
#' default with an adaptive Gauss-Hermite quadrature of the exact per-species RE
#' posterior (`control$re.aghq = TRUE`, the same correction the single-arm AGHQ
#' path applies; disable with `re.aghq = FALSE`, set nodes with `control$n.quad`).
#' The cover hurdle ties a species' psi / p / cover coefficients through the data,
#' so the per-species RE posterior is not separable across arms and the AGHQ is a
#' tensor product over the joint RE vector -- `n.quad^P` nodes per species in the
#' total RE dimension `P`. That cost is exponential in `P`, so the debias is a
#' hard scope limit: it runs only up to a small total RE dimension
#' (`control$re.aghq.maxdim`, default 4, e.g. intercept-only on all three arms),
#' and larger RE designs keep the EM variance. The EM (Laplace) variance is a
#' documented lower bound on the true component -- attenuated toward zero at small
#' per-species n, monotonically less attenuated as n grows -- not a bias of
#' unknown sign. `fit$ms_community$var_attenuation$debias` records `"aghq"` or
#' `"none"`, `$affects` names the lower-bounded components, and `print()` flags
#' the EM case.
#'
#' @param response likelihood for the positive cover arm. `"beta"` (cover in
#'   (0, 1)), `"lognormal"` (log-cover Gaussian), or `"gaussian"` (an
#'   identity-link Gaussian magnitude, the delta-normal hurdle; for a
#'   pre-transformed / unbounded positive response, not raw cover fractions).
#' @return A `tobs_family` object.
#' @seealso [occu_cover()] (single species), [ms_occu()] (community occupancy,
#'   no cover), [ms_abun()] (community N-mixture).
#' @export
ms_occu_cover <- function(response = c("beta", "lognormal", "gaussian")) {
  positive <- match.arg(response)
  obs_family(
    name           = "ms_occu_cover",
    class_long     = "community joint occupancy-detection + cover hurdle",
    latent         = "bernoulli",
    observation    = switch(positive,
                            beta      = "detection_plus_beta",
                            gaussian  = "detection_plus_gaussian",
                            "detection_plus_lognormal"),
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(positive = positive),
    control_keys   = c("max.iter", "tol", "sigma.beta", "newton.max", "sd.load",
                       "n.factors", "n.factors.max", "constrain",
                       "re.aghq.maxdim")
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
#' A detected visit (`y_cjv = 1`) with a missing cover (`y_pos_cjv = NA`) keeps
#' its detection term but drops the `f_pos` factor: cover is taken
#' missing-at-random, so the cover likelihood runs over the detected visits with
#' an observed cover. This is the rule [occu_cover()] applies, and it holds on
#' every engine below.
#'
#' @section Identifiability:
#' The availability (`theta`) and detection (`p`) levels separate only with
#' replication WITHIN a plot. Single releves supply none, so a plain fit
#' identifies `psi` (cell) and the product `theta * p` (plot) -- it reduces to
#' [occu_cover()] with `p := theta * p`. Within-plot temporal replication (e.g.
#' a resurvey of the same plot in a later period) makes the third level
#' estimable.
#'
#' @section Scope:
#' Three engines. `method = "nested_laplace"` carries a single shared areal
#' field coupled across the occupancy (`sigma`) and cover (`alpha * sigma`)
#' arms, integrated over the outer `(sigma, alpha)` grid. `method = "laplace"`
#' and `method = "nuts"` are the non-spatial path (iid cells, no field): the
#' exact three-level marginal (z over cells, a over plots both summed in closed
#' form) optimised directly (`"laplace"`, a Gaussian observed-Fisher posterior)
#' or sampled (`"nuts"`, the exact coefficient posterior with calibrated
#' intervals and WAIC / LOO). Cells are declared the same way on every path, via
#' an `icar(graph = adj, group_var = "<cell>")` term (the graph drives the field
#' under `"nested_laplace"` and supplies only the plot -> cell map under
#' `"laplace"` / `"nuts"`). On the `"nested_laplace"` path additional weighted
#' areal terms in the psi formula
#' (`icar(graph = adj, group_var = "<cell>", weight = <cell covariate>)`) add
#' spatially-varying-coefficient trend fields, each coupled onto the cover arm
#' with its own `alpha_trend`; the fitted fields are in `fit$trend_field` /
#' `fit$trend_fields`. The coupled / trend field is not sampled, so
#' `method = "nuts"` takes a single cell-declaring areal term only.
#'
#' The coupling is written as `copy(spatial())` in the `positive` formula, the
#' same spelling [occu_cover()] takes; it needs a field to copy, so it is
#' accepted under `method = "nested_laplace"` and refused on the two
#' non-spatial engines, which fix the field at 0. Omitting it leaves the
#' amplitude on the engine's default axis, estimated -- unlike [occu_cover()],
#' where a missing `copy()` pins `alpha = 0`.
#'
#' The copy amplitude's axis carries prior structure -- a point mass at
#' `alpha = 0` and a log-spaced slab above it -- so `copy(alpha = grid(...))`
#' STATES its nodes, while `copy(alpha = grid(n = ))` states a RESOLUTION: the
#' engine re-reads its own axis with that many slab nodes, point mass and bounds
#' unchanged. `terms =` gives either, per block, and
#' `control$alpha.grid[.trend]` / `control$alpha.n[.trend]` are the lower-level
#' spelling of the same pair. One block takes one of the two.
#'
#' @param response likelihood for the positive cover arm. `"beta"` (cover in
#'   (0, 1)), `"lognormal"` (log-cover Gaussian), or `"gaussian"` (an
#'   identity-link Gaussian magnitude, the delta-normal hurdle; for a
#'   pre-transformed / unbounded positive response, not raw cover fractions).
#' @return A `tobs_family` object.
#' @seealso [occu_cover()] (two-level), [cover()] (plot hurdle, no detection).
#' @export
occu_multiscale_cover <- function(response = c("beta", "lognormal", "gaussian")) {
  positive <- match.arg(response)
  obs_family(
    name           = "occu_multiscale_cover",
    class_long     = "three-level occupancy + cover hurdle",
    latent         = "bernoulli",
    observation    = switch(positive,
                            beta      = "availability_detection_plus_beta",
                            gaussian  = "availability_detection_plus_gaussian",
                            "availability_detection_plus_lognormal"),
    replicates     = "required",
    default_engine = "nested_laplace",
    status         = "working",
    params         = list(positive = positive),
    control_keys   = c(
      "max.iter", "tol", "sigma.beta",
      "sigma.grid", "alpha.grid", "alpha.grid.trend",
      "alpha.n", "alpha.n.trend", "phi.grid.pos", "n.threads",
      "inner.refresh", "hessian", "n.threads.outer", "force.sparse",
      "adaptive.grid", "adaptive.grid.edge.thresh", "adaptive.grid.max.passes",
      # Tuning for `integration = "grid_adaptive"`, a DIFFERENT mechanism from
      # the three refinement knobs above: the keep / expand radius from the
      # peak, the coarse-seed subsample stride per axis, the kept-fraction
      # ceiling past which the builder declines back to the dense tensor, and
      # the smallest dense tensor worth the machinery at all.
      "adaptive.grid.cutoff", "adaptive.grid.stride",
      "adaptive.grid.max.frac", "adaptive.grid.min.cells",
      "diagnose.k", "diagnose.draws", "k.samples", "k.bootstrap",
      "k.tail.points", "k.conf.bands",
      "checkpoint", "k.threads"
    )
  )
}


# ---------------------------------------------------------------------------
# Abundance, distance, removal, and specialised occupancy families
# ---------------------------------------------------------------------------

# Validate a latent-count truncation. `K_max` is the FIRST formal of abun() /
# ms_abun(), so abun("negbin") -- reaching for the mixing distribution -- binds
# the string to K_max instead. Left alone it coerces to NA and resurfaces much
# later as an unrelated comparison error inside a kernel, so catch it here where
# the fix is obvious.
.tobs_check_K_max <- function(K_max, family_name) {
  if (is.null(K_max)) return(invisible(NULL))
  if (!is.numeric(K_max) || length(K_max) != 1L || !is.finite(K_max) ||
      K_max < 1) {
    stop(sprintf(paste0("%s(): `K_max` must be a single positive number (the ",
                        "latent-count truncation), got %s. `K_max` is the first ",
                        "argument -- to set the mixing distribution write ",
                        "%s(mixture = \"negbin\")."),
                 family_name, deparse(K_max)[1L], family_name), call. = FALSE)
  }
  invisible(NULL)
}

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
#' @examples
#' \donttest{
#' sim <- simulate_abun(N = 120, J = 4, n_abund_covs = 1, n_det_covs = 1, seed = 1)
#' fit <- tobs(~ abund_cov1, data = sim$data, family = abun(),
#'             detection = ~ det_cov1, y = sim$y, method = "laplace")
#' summary(fit)
#' }
#' @export
abun <- function(K_max = NULL, mixture = c("poisson", "negbin", "zip", "zinb")) {
  mixture <- match.arg(mixture)
  .tobs_check_K_max(K_max, "abun")
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


#' Royle-Nichols occupancy family
#'
#' Occupancy with abundance-induced detection heterogeneity (Royle & Nichols
#' 2003; \pkg{unmarked} `occuRN`). Latent abundance `N_i ~ Poisson(lambda_i)`
#' drives the per-visit detection probability `1 - (1 - r_ij)^{N_i}`, where
#' `r_ij` is the per-individual detection probability. The state `formula` models
#' `log lambda` (abundance); `detection` models `logit r`. Detection is site-level
#' by default; passing `visits` makes it visit-varying (`logit r_ij` gains the
#' visit-level covariates), exactly as for the occupancy / N-mixture front doors.
#' The latent `N` marginalises in closed form (a Poisson sum to `K_max`), so the
#' fit maximises the exact marginal with an observed-information vcov.
#'
#' @param K_max Upper summation bound for the latent abundance (default: a
#'   data-driven Poisson-tail guess).
#' @return A `tobs_family` object for [tobs()].
#' @examples
#' \donttest{
#' sim <- simulate_royle_nichols(N = 150, J = 5, seed = 1)
#' fit <- tobs(~ x, data = sim$data, family = royle_nichols(),
#'             detection = ~ 1, y = sim$y, control = list(verbose = FALSE))
#' coef(fit)
#'
#' # Visit-varying detection via `visits`:
#' sv <- simulate_royle_nichols(N = 150, J = 5, beta_r_visit = 0.8, seed = 1)
#' fv <- tobs(~ x, data = sv$data, family = royle_nichols(),
#'            detection = ~ w, y = sv$y, visits = sv$visits,
#'            control = list(verbose = FALSE))
#' coef(fv)
#' }
#' @export
royle_nichols <- function(K_max = NULL) {
  obs_family(
    name           = "royle_nichols",
    class_long     = "Royle-Nichols occupancy",
    latent         = "poisson",
    observation    = "bernoulli_N",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(K_max = K_max)
  )
}


#' Time-to-detection occupancy family
#'
#' Occupancy where a survey records the TIME to first detection rather than a
#' 0/1 outcome (Garrard et al. 2008; \pkg{unmarked} `occuTTD`). At an occupied
#' site the time-to-detection is exponential with rate `lambda` (constant
#' hazard); a survey of length `surveyLength` that reaches its end without a
#' detection is censored. An unoccupied site never detects. The state `formula`
#' models `logit psi` (occupancy); `detection` models `log lambda` (the
#' site-level detection rate). The latent occupancy state integrates out in
#' closed form (two states), so the fit maximises the exact marginal with an
#' observed-information vcov.
#'
#' `y` is an `N x J` matrix of detection times: a value in `(0, surveyLength)`
#' is a detection; a value `>= surveyLength` is a non-detection (censored);
#' `NA` is a survey not conducted.
#'
#' @param surveyLength Survey length `Tmax` (the censoring time): a scalar, a
#'   length-`N` vector, or an `N x J` matrix. Default 1.
#' @return A `tobs_family` object for [tobs()].
#' @examples
#' \donttest{
#' sim <- simulate_occu_ttd(N = 200, J = 4, Tmax = 3, seed = 1)
#' fit <- tobs(~ psi_cov1, data = sim$data, family = occu_ttd(surveyLength = 3),
#'             detection = ~ rate_cov1, y = sim$y, control = list(verbose = FALSE))
#' coef(fit)
#' }
#' @export
occu_ttd <- function(surveyLength = 1) {
  obs_family(
    name           = "occu_ttd",
    class_long     = "time-to-detection occupancy",
    latent         = "bernoulli",
    observation    = "exponential_ttd",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(surveyLength = surveyLength)
  )
}


#' Multi-species co-occurrence occupancy family
#'
#' Joint occupancy of `S` species where the occupancy state `z in {0,1}^S`
#' follows a log-linear model with first-order (per species) and second-order
#' (per species pair) natural parameters (Rota et al. 2016; \pkg{unmarked}
#' `occuMulti`). The second-order parameters are the species interactions
#' (positive = co-occur more than independent; negative = avoid). The state
#' `formula` is the shared occupancy covariate design (each natural parameter
#' carries its own coefficients); `detection` the shared site-level per-species
#' detection design. The latent state integrates out by enumerating the `2^S`
#' states, so the exact marginal is maximised with an observed-information vcov.
#'
#' `y` is a length-`S` list of `N x J` 0/1/NA detection matrices (or a 3D
#' `[sites x visits x species]` array); `species` names the arms.
#'
#' @return A `tobs_family` object for [tobs()].
#' @examples
#' \donttest{
#' sim <- simulate_occu_multi(S = 2, N = 300, seed = 1)
#' fit <- tobs(~ scov1, data = sim$data, family = occu_multi(),
#'             detection = ~ 1, y = sim$y, species = sim$species,
#'             control = list(verbose = FALSE))
#' coef(fit)
#' }
#' @export
occu_multi <- function() {
  obs_family(
    name           = "occu_multi",
    class_long     = "multi-species co-occurrence occupancy",
    latent         = "bernoulli",
    observation    = "multi_state",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list()
  )
}


#' Double-observer abundance family
#'
#' Abundance from a double-observer protocol (\pkg{unmarked} `multinomPois` with a
#' double-observer pi-function). Site abundance `N ~ Poisson(lambda)` is surveyed
#' by two observers with detection `p1` / `p2`. The state `formula` models
#' `log lambda`; `detection` the shared site-level per-observer detection design
#' (observers 1 and 2 carry separate coefficients). By Poisson-multinomial
#' thinning the observable cell counts are independent Poissons, so the marginal
#' is closed form with no latent-abundance summation.
#'
#' @section Independent vs dependent protocol:
#' `type = "independent"` (the default): the two observers detect independently and
#' each individual is recorded by observer 1 only, observer 2 only, or both, so `y`
#' is an `N x 3` matrix of cell counts in that column order and the three cells are
#' `pi = (p1(1 - p2), (1 - p1)p2, p1 p2)`.
#'
#' `type = "dependent"`: a removal-style protocol where a *primary* observer
#' records what it detects and a *secondary* observer records only what the primary
#' missed, so `y` is an `N x 2` matrix `(primary-detected, secondary-only)` with
#' cells `pi = (p_pri, (1 - p_pri)p_sec)`. A single fixed primary observer gives
#' only two cells for three parameters (`lambda`, `p1`, `p2`), a ridge that does
#' not separate the two observer detections; observer **role-swapping** identifies
#' them -- pass `primary =` (a length-`N` vector in `{1, 2}` naming each site's
#' primary observer), and with observer 1 primary at some sites and observer 2 at
#' others `(p_pri, p_sec)` alternates between `(p1, p2)` and `(p2, p1)`, giving the
#' four cell means that recover `(lambda, p1, p2)`. With `p1 = p2` (an
#' interchangeable pair) the dependent protocol reduces to a two-pass [removal()]
#' model.
#'
#' @param type `"independent"` (default; `N x 3` cell counts) or `"dependent"`
#'   (removal-style, `N x 2` cell counts, needs `primary =` for identifiability).
#' @return A `tobs_family` object for [tobs()].
#' @examples
#' \donttest{
#' sim <- simulate_double_observer(N = 200, seed = 1)
#' fit <- tobs(~ abund_cov1, data = sim$data, family = double_observer(),
#'             detection = ~ det_cov1, y = sim$y, control = list(verbose = FALSE))
#' coef(fit)
#'
#' # Dependent (role-swapping) protocol:
#' sd <- simulate_double_observer(N = 300, type = "dependent", seed = 1)
#' fd <- tobs(~ abund_cov1, data = sd$data, family = double_observer("dependent"),
#'            detection = ~ det_cov1, y = sd$y, primary = sd$primary,
#'            control = list(verbose = FALSE))
#' coef(fd)
#' }
#' @export
double_observer <- function(type = c("independent", "dependent")) {
  type <- match.arg(type)
  obs_family(
    name           = "double_observer",
    class_long     = "double-observer abundance",
    latent         = "poisson",
    observation    = "multinomial_cells",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(type = type)
  )
}


#' Multi-season integrated occupancy family
#'
#' Dynamic (multi-season) occupancy observed by several detection sources
#' (spOccupancy `tIntPGOcc`): the product of a dynamic
#' occupancy HMM (season-1 occupancy `psi1`, colonization `gamma`, extinction
#' `eps`) and integrated occupancy (a per-season emission that pools multiple
#' detection sources). Pooling sources across seasons is how colonization /
#' extinction estimates come out of individually-sparse opportunistic data. The
#' latent occupancy sequence integrates out by the two-state HMM forward
#' recursion; the exact marginal is maximised with an observed-information vcov.
#'
#' `y` is a length-`S` list of `[sites x visits x seasons]` detection arrays
#' (0/1/NA). The state `formula` models `logit psi1`; `colonization = ~ ...` and
#' `extinction = ~ ...` the site-level transitions (required, as in [dyn_occu()]);
#' `detection` the shared per-source detection design (each source carries its own
#' coefficients). v1: every source covers all sites and the same season grid,
#' constant transitions, site-level detection.
#'
#' @return A `tobs_family` object for [tobs()].
#' @examples
#' \donttest{
#' sim <- simulate_dyn_int_occu(N = 200, T_seasons = 4, S = 2, seed = 1)
#' fit <- tobs(~ 1, data = sim$data, family = dyn_int_occu(),
#'             detection = ~ 1, colonization = ~ 1, extinction = ~ 1,
#'             y = sim$y, sources = sim$sources, control = list(verbose = FALSE))
#' coef(fit)
#' }
#' @export
dyn_int_occu <- function() {
  obs_family(
    name           = "dyn_int_occu",
    class_long     = "multi-season integrated occupancy",
    latent         = "bernoulli",
    observation    = "hmm_multisource",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list()
  )
}


#' Multispecies N-mixture family
#'
#' Per-species N-mixture with shared community-level hyperparameters.
#'
#' @inheritParams abun
#' @param mixture Abundance mixing distribution: `"poisson"` (default),
#'   `"negbin"`, or their zero-inflated counterparts `"zip"` / `"zinb"`. The
#'   negative-binomial dispersion and the zero-inflation structural-zero share
#'   are each a per-species random effect (`log_r_s ~ N(mu_log_r, sigma_log_r)`,
#'   `logit_omega_s ~ N(mu_omega, sigma_omega)`) integrated by the community AGHQ
#'   path alongside the abundance / detection coefficients. Zero-inflation is a
#'   non-spatial `laplace` fit with an intercept-only structural-zero logit; a
#'   shared field stays Poisson / negbin.
#' @section Reading the negbin community dispersion:
#' Under `mixture = "negbin"` the community mean `mu_log_r` reaches
#' `coef()` / `vcov()` / `confint()` with a marginal Wald SE, and the variance
#' component `sigma_log_r` is reported on `fit$ms_dispersion`. That interval is
#' calibrated **conditional on the variance component being recovered**, and the
#' two are not independent: `sigma_log_r` is a scalar variance over species and
#' at few species it can settle near its lower boundary, which shifts
#' `mu_log_r` and narrows its SE at the same time.
#'
#' Measured on simulated data (39 `laplace` fits at 8 and 36 species,
#' `sigma_logr = 0.5`): where `sigma_log_r` came back at
#' least 0.30, the nominal 95% interval covered 33 of 34 with a
#' `sqrt(mean(z^2))` of 0.88; where it came back below, it covered 2 of 5, the
#' point estimate was 2.2x further from the truth and the SE 28% narrower. So
#' read `fit$ms_dispersion$sigma_log_r` before quoting a `mu_log_r` interval: a
#' value near zero, well under what the data should support, means the interval
#' is not trustworthy even though the fit converged and reported no warning.
#' Those thresholds are the ones this fixture separated on and are not a general
#' rule; what transfers is the conditioning, not the number.
#'
#' Away from that boundary the interval is calibrated, and calibrated at every
#' group count measured. Over 97 fits at 8 / 18 / 36 species the scale
#' `sd(mu_log_r) / mean(se)` reads 1.066 / 1.265 / 0.981, and the middle figure
#' is the seed block rather than the estimator: `mu_log_r` is a POPULATION mean,
#' each seed draws `S` log-dispersions around it, and that draw supplies about
#' two thirds of the across-seed spread. The 18-species blocks drew theirs
#' 18-21% wider than `sigma_logr / sqrt(S)`. Put the draw at its expectation and
#' the scale is 1.077 / 1.101 / 0.977; rebuild the SE at the simulated sigma as
#' well and it is 0.990 / 1.035 / 0.963 (NOTES_measurements.md). The residual
#' there is the attenuation of
#' `sigma_log_r` itself (0.418 / 0.448 / 0.487 against 0.5), which shrinks as
#' species are added and which the SE inherits, since
#' `se^2 = (sigma_log_r^2 + c) / S`.
#'
#' A simulation study measuring this needs `simulate_ms_abun()`'s
#' `truth$mu_log_r_real` -- the mean of the log-dispersions that seed actually
#' drew -- beside `truth$mu_log_r`. Score coverage of a community mean against
#' the constant; score point recovery and interval SCALE against the realized
#' mean, or the draw is charged to the estimator.
#'
#' `control$logr.sigma.prior` puts a Penalized-Complexity prior on that variance,
#' which adds curvature at `sigma_log_r -> 0`. Measured on the same fixture at
#' `c(1, 0.05)`, it is worth reaching for when a fit reports a `sigma_log_r` near
#' zero or fails outright with a singular marginal Hessian -- one seed that
#' errored under pure maximum likelihood converged and covered under the prior,
#' and variance components near 0.01-0.06 lifted by an order of magnitude. It
#' does **not** repair coverage: on 19 paired seeds the nominal 95% interval
#' covered 17 either way, with the same two seeds missing, because the prior
#' bites near the boundary while those misses sit at `sigma_log_r` around 0.2.
#' It is off by default because the fits that were already calibrated pick up a
#' small systematic shift under it (`mu_log_r` by -0.006, p = 0.001).
#' @return A `tobs_family` object.
#' @export
ms_abun <- function(K_max = NULL,
                    mixture = c("poisson", "negbin", "zip", "zinb")) {
  mixture <- match.arg(mixture)
  .tobs_check_K_max(K_max, "ms_abun")
  obs_family(
    name           = "ms_abun",
    class_long     = "multispecies N-mixture",
    latent         = mixture,
    observation    = "binomial_N",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(K_max = K_max, mixture = mixture),
    # latent() factors and the spatial-factor route are fit by the
    # block-coordinate driver. A plain shared field with no factors keeps the
    # dedicated C++ path and reaches no outer loop.
    control_groups = c("block_coordinate", "block_coordinate_factor")
  )
}


#' Open-population (Dail-Madsen) N-mixture family
#'
#' Latent abundance evolves across primary seasons via apparent survival and
#' recruitment (Dail & Madsen 2011): `N_1 ~ Poisson(lambda)`; for `t >= 2`,
#' `N_t = S_t + G_t` with survivors `S_t ~ Binomial(N_{t-1}, omega)` and recruits
#' `G_t ~ Poisson(gamma)`; observed via `Binomial(N_t, p)` over secondary visits.
#' The latent abundance sequence is summed out by an exact HMM forward recursion
#' (it is not closed form, unlike the static [abun()]); analytic gradients come
#' from forward-mode differentiation of the scaled forward algorithm, so the fit
#' is a direct maximum-likelihood / Laplace fit with a NUTS path over the same
#' marginal.
#'
#' Four arms: initial abundance `lambda` (the `tobs()` `formula`) and detection
#' `p` (`detection`) are site-level; apparent survival `omega` (the `omega`
#' argument, default `~ 1`) and recruitment `gamma` (the `gamma` argument,
#' default `~ 1`) span the `T - 1` transition intervals. A constant `omega` /
#' `gamma` is shared across a site's seasons; supplying a season-varying
#' covariate (a `[n_sites x (T - 1)]` matrix column of `data`, one column per
#' transition interval) on `omega` / `gamma` gives interval-specific vital
#' rates. The response `y` is a 3D array `[n_sites x max_visits x n_seasons]` (or
#' a list of per-season count matrices); missing visits are `NA`.
#'
#' @param K_max abundance-state truncation for the forward recursion (states
#'   `0..K_max`). `NULL` (default) uses `max(count) + 40`; raise it if abundance
#'   may exceed that (the forward cost is roughly cubic in `K_max`).
#' @param mixture initial-abundance distribution: `"poisson"` (default),
#'   `"negbin"` (negative-binomial `N_1 ~ NB(mean = lambda, size = r)`), or their
#'   zero-inflated counterparts `"zip"` / `"zinb"` (a structural-zero share
#'   `omega` of sites is never occupied across any season; the remaining sites
#'   follow the Dail-Madsen open-population process). Zero-inflation is
#'   non-spatial `laplace` with an intercept-only `omega`; a field / RE / NUTS
#'   stay Poisson / negbin.
#' @return A `tobs_family` object.
#' @references Dail, D., Madsen, L. (2011). Models for estimating abundance from
#'   repeated counts of an open metapopulation. *Biometrics* 67, 577-587.
#' @export
dyn_abun <- function(K_max = NULL, mixture = c("poisson", "negbin", "zip", "zinb")) {
  mixture <- match.arg(mixture)
  obs_family(
    name           = "dyn_abun",
    class_long     = "Dail-Madsen open N-mixture",
    latent         = "dail_madsen",
    observation    = "binomial_N",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(K_max = K_max, mixture = mixture)
  )
}


#' Binned distance-sampling family
#'
#' Latent abundance `N_i ~ Poisson(lambda_i)` (or negative binomial) in a covered
#' region, observed through a half-normal or hazard-rate detection function over
#' distance bins. With `B` bins the detected counts are multinomial over
#' `(bin 1, ..., bin B, undetected)` with cell probabilities
#' `pi_b = integral_bin g(x; sigma) f(x) dx` and `1 - sum_b pi_b`, where `f(x)` is
#' the distance density (uniform for a line transect, proportional to distance for
#' a point transect). The latent `N` is summed out in closed form (truncation
#' `K_max`), so the fit is a direct Laplace approximation (no EM), with a NUTS
#' path over the same marginal; the per-bin detection integrals are evaluated by
#' Gauss-Legendre quadrature.
#'
#' The `tobs()` `formula` is the abundance (`log lambda`) model; `detection` is
#' the site-level detection-scale (`log sigma`) model. The response `y` is an
#' `n_sites x n_bins` integer matrix of per-bin detected counts. The bin edges and
#' transect geometry travel with the family: `distance(cutpoints = ...)`.
#'
#' @param key detection-function key. `"halfnorm"` (default) or `"hazard"`
#'   (the hazard-rate shape `b` is estimated as a scalar, reported as
#'   `log_shape`).
#' @param transect `"line"` (default; perpendicular distances uniform) or
#'   `"point"` (radial distances, density proportional to distance).
#' @param cutpoints numeric distance-bin edges, length `n_bins + 1`
#'   (`0 = c_0 < c_1 < ... < c_B`). Required.
#' @inheritParams abun
#' @return A `tobs_family` object.
#' @references
#' Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L., Borchers,
#'   D. L., Thomas, L. (2001). Introduction to Distance Sampling. Oxford.
#' Royle, J. A., Dawson, D. K., Bates, S. (2004). Modeling abundance effects in
#'   distance sampling. *Ecology* 85, 1591-1597.
#' @examples
#' \donttest{
#' sim <- simulate_distance(N = 200, key = "halfnorm", transect = "line",
#'                          n_abund_covs = 1, n_sigma_covs = 1, seed = 1)
#' fit <- tobs(~ abund_cov1, data = sim$data,
#'             family = distance(key = "halfnorm", transect = "line",
#'                               cutpoints = sim$cutpoints),
#'             detection = ~ sigma_cov1, y = sim$y, method = "laplace")
#' summary(fit)
#' }
#' @export
distance <- function(key = c("halfnorm", "hazard"),
                     transect = c("line", "point"),
                     cutpoints = NULL,
                     K_max = NULL, mixture = c("poisson", "negbin")) {
  key      <- match.arg(key)
  transect <- match.arg(transect)
  mixture  <- match.arg(mixture)
  obs_family(
    name           = "distance",
    class_long     = "binned distance sampling",
    latent         = mixture,
    observation    = "distance_binned",
    replicates     = "optional",
    default_engine = "laplace",
    status         = "working",
    params         = list(key = key, transect = transect,
                          cutpoints = cutpoints, K_max = K_max,
                          mixture = mixture)
  )
}


#' Community (multispecies) binned distance-sampling family
#'
#' The community analogue of [distance()]: per-species binned distance sampling
#' with Gaussian community hyperpriors on the per-species abundance and
#' detection-scale coefficients, so rare species borrow strength from common ones
#' through the shared community means and covariances -- the same pooling that
#' stabilises [ms_abun()], on the distance-sampling response.
#'
#' Per species `s`, site `i`, distance bin `b`:
#'
#'     N_{s,i}        ~ Poisson(lambda_{s,i})
#'     y_{s,i,.} | N  ~ Multinomial(N_{s,i}; pi_{s,i,1..B}, 1 - p_det)
#'     log lambda_{s,i} = X_i     . (mu_lambda + b_lambda_s)
#'     log sigma_{s,i}  = X_sig_i . (mu_sigma  + b_sigma_s)
#'     b_lambda_s ~ N(0, Sigma_lambda),  b_sigma_s ~ N(0, Sigma_sigma)
#'
#' with `pi_b` the integral of the detection function over bin `b`. The latent
#' `N_{s,i}` integrates out per species-site in closed form, exactly as for
#' [distance()]; the per-species coefficient deviations are the random effects,
#' integrated by the shared community Laplace-EM.
#'
#' Adding `latent(n)` to the abundance formula gives residual species
#' co-occurrence through `n` per-site latent factors with per-species loadings;
#' adding an areal term (`icar()` / `car_proper()` / `bym2()`) or `spde()`
#' alongside gives a shared spatial field. Both route through the same
#' block-coordinate engine as the other community families.
#'
#' @section Inputs:
#' `y` is a 3D array `[n_sites x n_bins x n_species]` or a named list of
#' `n_sites x n_bins` per-bin count matrices. `formula` is the abundance
#' (`log lambda`) predictor, `detection` the detection-scale (`log sigma`)
#' predictor, and `species` the species labels.
#'
#' @param key Detection function: `"halfnorm"` (default) or `"hazard"`. Under the
#'   hazard-rate key the scalar log-shape is shared across species.
#' @param transect Transect geometry: `"line"` (default) or `"point"`.
#' @param cutpoints Distance-bin edges, length `dim(y)[2] + 1`, strictly
#'   increasing and starting at `>= 0`.
#' @param K_max Truncation for the latent abundance sum. Defaults to
#'   `3 * max(rowSums(y)) + 100`.
#' @param mixture Abundance mixing distribution. `"poisson"` only; the
#'   negative-binomial size is not yet carried as a per-species random effect.
#' @return A `tobs_family` object.
#' @seealso [distance()] (single species), [ms_abun()] (community N-mixture).
#' @export
ms_distance <- function(key = c("halfnorm", "hazard"),
                        transect = c("line", "point"),
                        cutpoints = NULL,
                        K_max = NULL, mixture = c("poisson", "negbin")) {
  key      <- match.arg(key)
  transect <- match.arg(transect)
  mixture  <- match.arg(mixture)
  obs_family(
    name           = "ms_distance",
    class_long     = "multispecies binned distance sampling",
    latent         = mixture,
    observation    = "distance_binned",
    replicates     = "optional",
    default_engine = "laplace",
    status         = "working",
    params         = list(key = key, transect = transect,
                          cutpoints = cutpoints, K_max = K_max,
                          mixture = mixture),
    # max.iter / tol / sigma.beta come from the laplace_em group; quad.order (the
    # Gauss-Legendre order for the per-bin detection integrals) is this family's
    # own knob.
    control_keys   = "quad.order",
    # latent() factors and a shared field are fit by the block-coordinate driver
    control_groups = c("block_coordinate", "block_coordinate_factor")
  )
}


#' Removal-sampling family
#'
#' Sequential-depletion removal sampling: latent abundance
#' `N_i ~ Poisson(lambda_i)` (or negative binomial) observed through `K` ordered
#' removal passes, where pass `k` removes
#' `Binomial(N_i - sum_{l<k} y_{il}, p_{ik})` of the individuals still present.
#' The declining catch sequence identifies detection `p` and abundance `N`; the
#' latent `N` is summed out in closed form (truncation `K_max`), so the fit is a
#' direct Laplace approximation (no EM), with a NUTS path over the same marginal.
#'
#' The `tobs()` `formula` is the abundance (`log lambda`) model; `detection` is
#' the per-pass detection (`logit p`) model. The response `y` is an
#' `n_sites x K` integer matrix of per-pass removals with the passes in column
#' order; complete pass sequences are required (no `NA`).
#'
#' @inheritParams abun
#' @return A `tobs_family` object.
#' @references
#' Royle, J. A. (2004). N-mixture models for estimating population size from
#'   spatially replicated counts. *Biometrics* 60, 108-115.
#' Dorazio, R. M., Jelks, H. L., Jordan, F. (2005). Improving removal-based
#'   estimates of abundance. *Biometrics* 61, 1093-1101.
#' @export
removal <- function(K_max = NULL, mixture = c("poisson", "negbin")) {
  mixture <- match.arg(mixture)
  obs_family(
    name           = "removal",
    class_long     = "removal sampling",
    latent         = mixture,
    observation    = "removal_sequence",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(K_max = K_max, mixture = mixture)
  )
}


#' Multistate false-positive occupancy family
#'
#' Occupancy with both false negatives and false positives in the detection
#' process (e.g. acoustic classifiers, citizen-science misidentification), using
#' the Miller et al. (2011) confirmed-detection design that makes the model
#' robustly identifiable. Each visit yields a detection state `y in {0, 1, 2}`:
#' `0` = no detection, `1` = ambiguous detection (a true detection OR a false
#' positive), `2` = certain / confirmed detection (only possible when the site is
#' truly occupied). The latent occupancy `z` marginalises in closed form (two
#' states), so the fit maximises the exact marginal likelihood directly
#' (analytic-gradient BFGS, observed-information covariance) with a NUTS path over
#' the same marginal.
#'
#' Four site-level logit arms: occupancy `psi` (the `tobs()` `formula`), true
#' detection `p11` (`detection`), false-positive rate `p10`, and the probability a
#' true detection is certain `b`. The `p10` and `b` predictors default to
#' intercept-only and are set with the `p10 = ~ ...` and `certainty = ~ ...`
#' arguments to [tobs()] (`certainty` is the `b` arm). The response `y` is an
#' `n_sites x J` integer matrix in `{0, 1, 2}` (NA visits dropped).
#'
#' @return A `tobs_family` object.
#' @references
#' Miller, D. A. W., Nichols, J. D., McClintock, B. T., Grant, E. H. C., Bailey,
#'   L. L., Weir, L. A. (2011). Improving occupancy estimation when two types of
#'   observational error occur. *Ecology* 92, 1422-1428.
#' Royle, J. A., Link, W. A. (2006). Generalized site occupancy models allowing
#'   for false positive and false negative errors. *Ecology* 87, 835-841.
#' @export
fp_occu <- function() {
  obs_family(
    name           = "fp_occu",
    class_long     = "multistate false-positive occupancy",
    latent         = "bernoulli",
    observation    = "multistate_detection",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working"
  )
}


#' Joint distance + removal sampling family
#'
#' The `unmarked` `gdistremoval` model (Amundson et al. 2014): a **single-season**
#' point-count design in which the detected individuals are recorded two ways at
#' once -- their distance band (distance sampling) and the removal period of first
#' detection (removal sampling). This is NOT the open-population distance model
#' (`distsampOpen`); the abundance is static.
#'
#' Site abundance `N_i ~ Poisson(lambda_i)`; the detected birds are
#' cross-classified by a distance band and a removal period. Writing `pdist_i` for
#' the overall distance detection, `prem_i` for the overall removal detection, the
#' total detected is a binomial thinning of `N_i`, and Poisson is closed under
#' binomial thinning, so the marginal is closed-form:
#' \deqn{\mathrm{ysum}_i \sim \mathrm{Poisson}(\lambda_i\, p^{dist}_i\, p^{rem}_i),}
#' with the distance-band counts and removal-period counts as two conditional
#' multinomials -- the [double_observer()] Poisson-multinomial pattern, here with
#' a distance multinomial (half-normal key band integrals) and a depleting-removal
#' multinomial `pi_k = r (1 - r)^{k-1}`.
#'
#' Three site-level arms: log abundance `lambda` (the [tobs()] `formula`), log
#' distance scale `sigma` (`detection`), and logit per-period removal capture `r`
#' (`removal = ~ ...`, default intercept-only).
#'
#' @section Inputs:
#' `y` is an `n_sites x n_bins` integer matrix of per-distance-band counts;
#' `y_rem` an `n_sites x n_periods` integer matrix of per-removal-period counts.
#' The per-site row totals must match (the same detected birds cross-classified).
#'
#' @param transect Transect geometry: `"line"` (default) or `"point"`.
#' @param cutpoints Distance-bin edges, length `ncol(y) + 1`, strictly increasing
#'   and starting at `>= 0`.
#' @return A `tobs_family` object.
#' @references
#' Amundson, C. L., Royle, J. A., Handel, C. M. (2014). A hierarchical model
#'   combining distance sampling and time removal to estimate detection
#'   probability during avian point counts. *The Auk* 131, 476-494.
#' @seealso [distance()], [removal()], [double_observer()].
#' @export
gdistremoval <- function(transect = c("line", "point"), cutpoints = NULL) {
  transect <- match.arg(transect)
  obs_family(
    name           = "gdistremoval",
    class_long     = "joint distance + removal sampling",
    latent         = "poisson",
    observation    = "distance_removal",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(transect = transect, cutpoints = cutpoints)
  )
}


#' Open-population distance-sampling family
#'
#' The `unmarked` `distsampOpen` model: a Dail-Madsen open N-mixture (as
#' [dyn_abun()]) observed by distance sampling at each primary period -- the
#' open-population counterpart of the single-season [gdistremoval()]. Initial
#' abundance `N_1 ~ Poisson(lambda)`, then `N_t = Binomial(N_{t-1}, omega) +
#' Poisson(gamma)` (apparent survival `omega` + recruitment `gamma`); at each
#' primary period the detected birds are distance-sampled into bins.
#'
#' The distance-band allocation is conditional on the period total detected, so it
#' factors out of the abundance HMM (the [gdistremoval()] trick): the marginal is
#' the [dyn_abun()] forward recursion with the detection probability set to the
#' overall distance detection `pdist`, plus the per-period band multinomials.
#'
#' Four site-level arms: log abundance `lambda` (the [tobs()] `formula`), log
#' distance scale `sigma` (`detection`), logit survival `omega` (`omega = ~ ...`),
#' and log recruitment `gamma` (`gamma = ~ ...`), the last two intercept-only by
#' default.
#'
#' @section Inputs:
#' `y` is a 3D array `[n_sites x n_bins x n_seasons]` of per-distance-band counts
#' at each primary period (secondary occasions absorbed into the period total).
#'
#' @param transect Transect geometry: `"line"` (default) or `"point"`.
#' @param cutpoints Distance-bin edges, length `dim(y)[2] + 1`, strictly
#'   increasing and starting at `>= 0`.
#' @param K_max Truncation for the latent abundance HMM. Defaults to
#'   `3 * max(period total) + 40`.
#' @param mixture Initial-abundance mixing distribution. `"poisson"` (default),
#'   `"negbin"` (negative binomial, `Var(N_1) = lambda + lambda^2 / r`, with an
#'   overdispersion `r` estimated jointly and reported as `log_r`), `"zip"`
#'   (zero-inflated Poisson) or `"zinb"` (zero-inflated negative binomial). A
#'   structural-zero site is never occupied across any primary period, so all its
#'   band counts are zero; the observed per-site marginal is the two-component
#'   mixture `omega * 1{all zero} + (1 - omega) * L_open`, with `L_open` the exact
#'   open-population distance marginal and `omega` an intercept-only structural-zero
#'   probability reported as `zi_logit` (distinct from `omega`, the survival arm).
#'   The negative-binomial size / zero-inflation is layered over the same
#'   forward-HMM marginal; the Poisson path is unchanged.
#' @param dynamics Population-dynamics form for the transition
#'   `N_t | N_{t-1}` (following `unmarked::distsampOpen()`):
#'   `"constant"` (default) is `N_t = Binomial(N_{t-1}, omega) + Poisson(gamma)`;
#'   `"notrend"` ties recruitment so the expected abundance is stationary
#'   (`gamma = (1 - omega) * lambda`, no free `gamma` arm); `"trend"` is an
#'   exponential-growth process `N_t ~ Poisson(N_{t-1} * gamma)` with no survival
#'   term (`omega` unused); `"autoreg"` adds density-dependent recruitment
#'   `Poisson(N_{t-1} * gamma)` on top of survival; `"ricker"` and `"gompertz"`
#'   are the density-regulated recruitment forms with an estimated carrying
#'   capacity `K` (reported as `K`, on the log scale, modelled through the
#'   `omega = ~ ...` formula slot) and a growth rate reported as `r` (modelled
#'   through the `gamma = ~ ...` slot).
#'   `"constant"` / `"notrend"` fit with the exact analytic gradient; the
#'   density-dependent forms (`trend` / `autoreg` / `ricker` / `gompertz`) use the
#'   exact forward-HMM marginal with a numeric gradient.
#' @return A `tobs_family` object.
#' @references
#' Dail, D., Madsen, L. (2011). Models for estimating abundance from repeated
#'   counts of an open metapopulation. *Biometrics* 67, 577-587.
#' Sollmann, R., Gardner, B., Chandler, R. B., Royle, J. A., Sillett, T. S.
#'   (2015). An open-population hierarchical distance sampling model. *Ecology*
#'   96, 325-331.
#' @seealso [gdistremoval()] (single-season), [dyn_abun()] (open N-mixture),
#'   [distance()].
#' @export
distsamp_open <- function(transect = c("line", "point"), cutpoints = NULL,
                          K_max = NULL,
                          mixture = c("poisson", "negbin", "zip", "zinb"),
                          dynamics = c("constant", "notrend", "trend", "autoreg",
                                       "ricker", "gompertz")) {
  transect <- match.arg(transect)
  mixture  <- match.arg(mixture)
  dynamics <- match.arg(dynamics)
  obs_family(
    name           = "distsamp_open",
    class_long     = "open-population distance sampling",
    latent         = mixture,
    observation    = "distance_open",
    replicates     = "required",
    default_engine = "laplace",
    status         = "working",
    params         = list(transect = transect, cutpoints = cutpoints,
                          K_max = K_max, mixture = mixture, dynamics = dynamics)
  )
}


#' Presence + nominal class hurdle family
#'
#' A hurdle for a response that is either absent or, when present, one of `K`
#' **nominal (unordered)** classes -- a colour morph, microhabitat / substrate
#' use, or a cryptic-species / classifier label observed given the organism is
#' present. The presence/absence part is a Bernoulli arm; the class given present
#' is a baseline-category multinomial logit (the last class is the baseline). The
#' two arms factorise the likelihood exactly:
#' \deqn{P(y = 0) = 1 - \psi, \qquad P(y = k) = \psi \, p_k,}
#' with \eqn{p = \mathrm{softmax}(X \beta_{class})}.
#'
#' This is the categorical counterpart of [cover()] (presence + a magnitude):
#' here the positive part is an unordered class, so it uses a multinomial logit
#' rather than beta / lognormal. For an *ordered* class response (Braun-Blanquet
#' cover bands) use `cover(response = "ordinal")`, which exploits the ordering;
#' this family is for classes with no ordering.
#'
#' @section Response:
#' `y` is a length-N integer vector in `0..K`: `0` is absent, `k` is class `k`.
#' It may sit on the formula left-hand side (`y ~ ...`), so `y =` can be dropped.
#' The `formula` predictors are shared by both arms.
#'
#' @section Scope:
#' Non-spatial Laplace (`method = "laplace"`). The multinomial math is the
#' FD-validated tulpa kernel (`multinomial_logit.h`); the non-spatial fit is the
#' vectorised R Newton over the same closed forms. Spatial fields / NUTS (the
#' native multi-process likelihood) and the latent-class *misclassification*
#' variant (the K-class generalisation of [fp_occu()], a confusion matrix on the
#' observed label) are documented follow-ups.
#'
#' @param classes optional character vector of class labels (length `K`), used
#'   only to name the coefficient blocks; when `NULL`, `K` is taken from
#'   `max(y)` and the classes are labelled `1..K`.
#' @return A `tobs_family` object.
#' @seealso [cover()] (presence + magnitude), [fp_occu()] (two-state
#'   false-positive detection).
#' @export
occu_categorical <- function(classes = NULL) {
  if (!is.null(classes) &&
      (!is.character(classes) || length(classes) < 2L || anyNA(classes))) {
    stop("`classes` must be a character vector of at least two class labels, ",
         "or NULL to infer them from the data.", call. = FALSE)
  }
  obs_family(
    name           = "occu_categorical",
    class_long     = "presence + nominal class hurdle",
    latent         = "hurdle",
    observation    = "binomial_plus_multinomial",
    replicates     = "single",
    default_engine = "laplace",
    status         = "working",
    response       = "vector",
    params         = list(classes = classes),
    control_keys   = c("max.iter", "tol", "prior.prec", "sigma.beta")
  )
}


#' Cover hurdle family (vegetation cover, MOTIVATE pattern)
#'
#' Latent presence (Bernoulli) plus conditional positive cover (beta or
#' lognormal). Does not share the replicate-detection assumption of the other
#' families — see `vignette("families")` for the conceptual caveat.
#'
#' @section Response on the formula left-hand side:
#' The cover response is a single length-N vector, so it may sit on the top
#' formula left-hand side and `y =` is dropped. The two
#' calls are equivalent:
#'
#' ```r
#' # response on the LHS (y = omitted)
#' tobs(cover.flat ~ time.sc + habitat +
#'        spatial(~ 1 + time.sc || cell_idx, graph = adj),
#'      data = dat, family = cover(response = "beta"),
#'      method = "nested_laplace")
#'
#' # the same fit with a one-sided formula and an explicit y =
#' tobs(~ time.sc + habitat +
#'        spatial(~ 1 + time.sc || cell_idx, graph = adj),
#'      y = dat$cover.flat, data = dat, family = cover(response = "beta"),
#'      method = "nested_laplace")
#' ```
#'
#' Naming the response makes the per-arm spatial labels read naturally: `cover()`
#' splits `cover.flat` into a `presence` arm and a `positive` arm, the arm names
#' that per-arm formulas and copy() address. The LHS is evaluated against `data`
#' (then the calling environment), so it may be a bare column or an expression.
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
#' @section Spatially varying trend (a weighted areal term):
#' A spatially varying trend is model structure, so it is declared in the
#' formula as a SECOND, weighted areal term on the same graph as the intercept
#' field, not through `control`:
#'
#' ```r
#' ~ time.sc +
#'   icar(graph = adj, group_var = "cell_idx") +
#'   icar(graph = adj, weight = time.sc, group_var = "cell_idx")
#' ```
#'
#' The unweighted term is the shared intercept field; the weighted term
#' `icar(..., weight = col)` is the trend field, whose contribution to each
#' arm's predictor is `weight_i * z[cell_i]`, coupled onto the cover arm with
#' its own scale (`alpha_trend`, reported in `fit$alpha_trend` /
#' `fit$sigma_trend`) integrated over the outer grid. The umbrella spelling
#' `spatial(graph = adj, model = "icar", weight = col)` resolves identically.
#' The trend coupling grid defaults to `control$alpha.grid`; override it with
#' `control$alpha.grid.trend`. Requires `method = "nested_laplace"`. A coupled
#' trend cannot currently combine with `temporal()` / `re()` blocks in the same
#' fit.
#'
#' The axis the amplitude rides carries prior structure -- a point mass at
#' `alpha = 0` ("no coupling") and a log-spaced slab above it -- so it is set
#' in one of two ways. `copy(alpha = grid(...))` STATES its nodes, and with them
#' that structure. `copy(alpha = grid(n = 9))` states a RESOLUTION: the engine
#' re-reads its own axis with that many slab nodes, point mass and slab bounds
#' unchanged, so sharpening the axis never restates its structure. `terms =`
#' gives either, per block. The resolution is the only way to raise this axis,
#' because it does not densify when the donor `control$sigma.grid` does: on an
#' informative data set the outer grid's quadrature effective sample size
#' saturates on the copy amplitude while every other axis tracks the request
#' (measured engine-side, `NOTES_measurements.md`). One block takes one of the
#' two; giving both is an error. `control$alpha.grid[.trend]` and
#' `control$alpha.n[.trend]` are the lower-level spelling of the same pair.
#'
#' Set it when the fit's hyperparameter intervals are reported. On the coupled
#' SBC fixture the declared resolution leaves the field SD miscalibrated once
#' the data are informative -- at 10 visits per site its uniformity p-value is
#' 9.1e-05, against 0.17 at 3 visits -- and `control$alpha.n = 21` returns every
#' scored quantity to nominal at both. Thirteen slab nodes do not (the field SD
#' still reads 4.5e-03), and the price is a 2.3-3x longer fit, since the outer
#' grid is a tensor. Point estimates are not affected the way the intervals are;
#' the measurement is in `NOTES_measurements.md`.
#'
#' @section Varying-coefficient spatial bar (the compact single-term form):
#' The intercept field plus its weighted trend field can also be written as one
#' `spatial()` term carrying an lme4-style coefficient formula:
#'
#' ```r
#' ~ time.sc + habitat +
#'   spatial(~ 1 + time.sc || cell_idx, graph = adj)
#' ```
#'
#' The bar left-hand side spells the coefficient fields: the intercept column
#' (`1`) is the unweighted field; each covariate column (`time.sc`) is a
#' weight-scaled coefficient field (`weight_i * z[cell_i]`). The bar right-hand
#' side (`cell_idx`) is the graph node index (the areal `group_var`); `||`
#' requests independent intercept and slope fields, a single `|` makes them
#' correlated. This desugars to exactly the two-term weighted-areal form above,
#' so the two spellings give the same fit.
#'
#' @section Choosing a field's arm: placement and copy():
#' The cover hurdle's two arms are `presence` (the `y > 0` Bernoulli arm) and
#' `positive` (the `y | y > 0` arm); `summary()` and the coefficient output print
#' these same labels. A field is placed on an arm by writing it in that arm's
#' per-arm formula (`presence = ~ ...`, `positive = ~ ...`), and shared across
#' arms with copy(). A field in the single shared `formula` reaches both arms.
#'
#' \emph{Both arms (shared / copied).} A field in the shared `formula` is one
#' presence-anchored latent copied onto the positive arm with an estimated
#' coupling per coefficient field (`alpha` for the intercept field, `alpha_trend`
#' for the trend field), marginalized on the outer grid:
#'
#' ```r
#' # one latent field, shared by both arms
#' tobs(cover.flat ~ time.sc + habitat +
#'        spatial(~ 1 + time.sc || cell_idx, graph = adj),
#'      data = dat, family = cover(response = "beta"), method = "nested_laplace")
#' ```
#'
#' ```
#' presence: eta_presence = ... + u_cell + time.sc * s_cell
#' positive: eta_positive = ... + alpha * u_cell + alpha_trend * time.sc * s_cell
#' ```
#'
#' Per-arm formulas make the shared field explicit: place it on `presence` and
#' copy it onto `positive`. `copy(spatial())` estimates the coupling on the
#' default axis; `copy(spatial(), alpha = grid(c(...)))` integrates it over
#' supplied nodes, `alpha = grid(n = 9)` over the engine's own axis read at `n`
#' slab nodes, and `alpha = 0.5` fixes it. `copy(spatial(), prior = list(...))`
#' regularizes the coefficient itself; a fit copying both the intercept and a
#' weighted trend block is refused it, since one prior reaches the engine per
#' fit (gcol33/tulpa#655).
#'
#' ```r
#' tobs(presence = ~ time.sc + habitat +
#'                   spatial(~ 1 + time.sc || cell_idx, graph = adj),
#'      positive = ~ time.sc + habitat + copy(spatial()),
#'      data = dat, family = cover(response = "beta"), method = "nested_laplace")
#' ```
#'
#' \emph{One arm only (free / separate).} A field written in one arm's formula is
#' a separate latent on that arm alone, with its own precision and no cross-arm
#' coupling. A field in each arm's formula gives two
#' independent latents:
#'
#' ```r
#' tobs(presence = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
#'      positive = ~ time.sc + spatial(~ 1 + time.sc || cell_idx, graph = adj),
#'      data = dat, family = cover(response = "beta"), method = "nested_laplace")
#' ```
#'
#' ```
#' presence: eta_presence = ... + u_presence_cell + time.sc * s_presence_cell
#' positive: eta_positive = ... + u_positive_cell + time.sc * s_positive_cell
#' ```
#'
#' The free fit reports `sigma_armspecific`. Arm-specific fields are their own
#' spatial structure: they do not combine with a shared field, a weighted trend,
#' or `temporal()` / `re()` in the same formula, and at most one targets each arm.
#'
#' The `||` and `|` axis is separate from shared / free: `||` makes the intercept
#' and slope fields independent, while a single `|` makes them correlated (a free
#' cross-covariance, MCAR). A correlated `|` bar shared across both arms (in the
#' shared `formula`, or on `presence` with `copy()` on `positive`) is copied with
#' one estimated amplitude; placed on one arm's formula it is a free-Sigma
#' correlated field on that arm alone, no cross-arm copy.
#'
#' ```
#' field in the shared formula (or copy())   one shared / copied latent (presence anchor, coupling estimated)
#' a field in each arm's formula             separate / free latents, no coupling
#' ||                                        independent intercept and slope coefficient fields
#' |                                         correlated (MCAR) coefficient fields, copy-only
#' ```
#'
#' @section A formula bar is a random effect, not a spatial field:
#' A bare lme4 bar in the formula -- `(1 | cell)`, `(1 + x | cell)`,
#' `(x || cell)` -- is a grouped random effect, not a spatial field, even when the
#' grouping factor names the areal graph nodes. The engine's inline-MCAR call
#' `tulpa::spatial(graph, ~ 1 + x | cell)` reads `| cell` as a separable spatial
#' field, but the same spelling inside a `cover()` formula is parsed as an IID
#' random effect on `cell`. For a spatial field on the cells write either the
#' compact `spatial()` bar of the section above,
#' `spatial(~ 1 + x || cell, graph = adj)`, or the two-term weighted-areal form,
#' `icar(graph = adj, group_var = "cell") +`
#' `icar(graph = adj, weight = x, group_var = "cell")`. When a formula carries a
#' bar whose grouping factor is also an areal term's `group_var`, `cover()` emits
#' a one-time message noting the bar is fitted as a random effect; suppress it
#' with [base::suppressMessages()].
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
#' @param response likelihood for the positive cover part:
#'   * `"beta"` -- cover in (0, 1) with a logit link;
#'   * `"beta_oi"` -- a one-inflated Beta: a point mass at cover = 1 (plots
#'     recorded at full cover, a genuine boundary mass rather than a near-1
#'     continuous value) plus a Beta on the interior (0, 1). The inflation
#'     probability is a constant, estimated as the share of positive plots at
#'     the ceiling; the interior Beta is fit on the (0, 1) plots. Reported as
#'     `pi_one`. (With the zero hurdle this is the zero-one-inflated Beta.);
#'   * `"lognormal"` -- a Gaussian on log-cover (unbounded above);
#'   * `"lognormal_trunc"` -- a Gaussian on log-cover upper-truncated at
#'     `log(1) = 0`, the bounded-support form for cover in (0, 1] (it cannot
#'     place mass above cover = 1 the way `"lognormal"` can). Requires
#'     `method = "nested_laplace"`;
#'   * `"ordinal"` -- an interval-censored Gaussian on log-cover with KNOWN
#'     class thresholds (an ordered probit for Braun-Blanquet-style class data),
#'     set with `breaks`. Requires `method = "nested_laplace"`;
#'   * `"gaussian"` -- an identity-Gaussian magnitude (the delta-normal hurdle):
#'     Bernoulli presence times a Gaussian on the raw response. For a
#'     pre-transformed or otherwise unbounded response on a real scale (it
#'     permits negative fitted values), NOT for cover fractions in (0, 1); those
#'     stay on `"beta"` / `"lognormal"` / `"ordinal"`. Presence is the nonzero
#'     sentinel (`y != 0`) rather than `y > 0`.
#' @param breaks for `positive = "ordinal"` only, the interior cover-class
#'   boundaries on the (0, 1) cover-fraction scale (strictly ascending, all in
#'   (0, 1)); the open outer classes are added automatically. `NULL` for the
#'   other families.
#' @return A `tobs_family` object.
#' @examples
#' \donttest{
#' sim <- simulate_cover(N = 200, seed = 1)
#' fit <- tobs(~ x, data = sim$data, family = cover("beta"), y = sim$y)
#' summary(fit)
#' }
#' @export
cover <- function(response = c("beta", "beta_oi", "lognormal", "lognormal_trunc",
                               "ordinal", "gaussian"),
                  breaks = NULL) {
  positive <- match.arg(response)
  # The ordinal positive arm is an interval-censored Gaussian on log-cover with
  # KNOWN class thresholds (an ordered probit, not free cutpoints): the cover is
  # recorded only as the ordinal class it falls in, so the likelihood is the
  # class probability MASS -- a genuine PMF with no change-of-variable Jacobian,
  # the measure-invariant counterpart of beta / lognormal for Braun-Blanquet
  # data. `breaks` are the interior class boundaries on the (0, 1) cover-fraction
  # scale; the open outer classes are added automatically.
  if (positive == "ordinal") {
    if (is.null(breaks) || !is.numeric(breaks) || length(breaks) < 1L ||
        anyNA(breaks) || is.unsorted(breaks, strictly = TRUE) ||
        any(breaks <= 0) || any(breaks >= 1)) {
      stop("cover(response = \"ordinal\") requires `breaks`: the interior ",
           "cover-class boundaries on the (0, 1) cover-fraction scale, strictly ",
           "ascending and all in (0, 1). For the MOTIVATE Braun-Blanquet scheme ",
           "(myscale): c(0.002, 0.015, 0.03, 0.05, 0.25, 0.50, 0.75). The open ",
           "outer classes (0, breaks[1]] and (breaks[K], 1) are added ",
           "automatically.", call. = FALSE)
    }
  } else if (!is.null(breaks)) {
    stop("`breaks` is only used for cover(response = \"ordinal\").", call. = FALSE)
  }
  obs_family(
    name           = "cover",
    class_long     = "vegetation cover hurdle",
    latent         = "hurdle",
    observation    = switch(positive,
                            beta            = "binomial_plus_beta",
                            beta_oi         = "binomial_plus_beta_one_inflated",
                            lognormal       = "binomial_plus_lognormal",
                            lognormal_trunc = "binomial_plus_lognormal_trunc",
                            ordinal         = "binomial_plus_ordinal",
                            gaussian        = "binomial_plus_gaussian"),
    replicates     = "single",
    default_engine = "laplace",
    status         = "working",
    # The cover response is a plain length-N cover vector, so it may sit on the
    # top formula LHS (`cover.flat ~ ...`) and drop `y =`.
    response       = "vector",
    params         = list(positive = positive, breaks = breaks),
    # The cover hurdle has its own (.dispatch_cover) grid-based control surface,
    # separate from the occupancy fitter and named with underscores. Declaring
    # the keys keeps tobs()'s control validation from rejecting them. `trend`
    # is retained here only so a `control$trend` left over from the old API
    # reaches .dispatch_cover's removal error (a weighted formula term now
    # carries the trend), rather than a generic unknown-control-key rejection.
    # `sigma.pos.grid` is retained for the same reason: its removal error names
    # `alpha.grid`, the axis that replaced it.
    control_keys   = c(
      "max.iter", "tol", "n.threads", "n.threads.outer", "prior.sigma", "prior.alpha",
      "prior.phi",
      "phi.grid", "sigma.grid", "sigma.pos.grid", "rho.grid", "tau.grid",
      "rho.car.grid", "tau.temporal.grid", "rho.temporal.grid",
      "sigma.temporal.grid", "sigma.re.grid",
      "trend", "alpha.grid", "alpha.grid.trend", "alpha.n", "alpha.n.trend",
      "integration",
      "adaptive.grid", "adaptive.grid.edge.thresh", "adaptive.grid.max.passes",
      # Tuning for `integration = "grid_adaptive"`, a DIFFERENT mechanism from
      # the three refinement knobs above: the keep / expand radius from the
      # peak, the coarse-seed subsample stride per axis, the kept-fraction
      # ceiling past which the builder declines back to the dense tensor, and
      # the smallest dense tensor worth the machinery at all.
      "adaptive.grid.cutoff", "adaptive.grid.stride",
      "adaptive.grid.max.frac", "adaptive.grid.min.cells",
      "var.of.means.consistency", "var.of.means.min.ess",
      "prune", "prune.tol", "hessian", "aggregate.occ", "aggregate.pos",
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
  cat(sprintf("  response       : %s\n", x$response %||% "matrix"))
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
