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


#' Multispecies (community) occupancy family
#'
#' Per-species occupancy with shared community-level hyperparameters on the
#' species random effects.
#'
#' @return A `tobs_family` object.
#' @export
ms_occu <- function() {
  obs_family(
    name           = "ms_occu",
    class_long     = "multispecies occupancy",
    latent         = "bernoulli",
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


# ---------------------------------------------------------------------------
# Planned families — error informatively when used until implemented
# ---------------------------------------------------------------------------

#' N-mixture abundance family
#'
#' Latent Poisson (or NB) abundance with binomial detection per visit
#' (Royle 2004).
#'
#' @param K_max upper bound for latent N in the EM E-step summation.
#' @param mixture latent-abundance distribution: `"poisson"` or `"negbin"`.
#'   Named `mixture` (after `unmarked::pcount()`) to avoid collision with the
#'   model-type `family` argument of [tobs()].
#' @return A `tobs_family` object.
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
    status         = "planned",
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
    status         = "planned",
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
      "max.iter", "tol", "n.threads", "prior.sigma", "prior.alpha",
      "phi.grid", "sigma.grid", "sigma.pos.grid", "rho.grid", "tau.grid",
      "rho.car.grid", "tau.temporal.grid", "rho.temporal.grid",
      "sigma.temporal.grid", "sigma.re.grid",
      "adaptive.grid", "adaptive.grid.edge.thresh", "adaptive.grid.max.passes"
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
