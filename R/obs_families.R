# =============================================================================
# obs_families.R — Family-object constructors for tulpaObs
#
# Each family is a small list with class "tulpa_obs_family". It carries
# the latent-state type, the observation likelihood, replicate requirements,
# default engines, and any family-specific hyperparameters (K_max, beta vs
# lognormal positive part, etc.).
#
# Families do not fit models. tulpa_obs() reads the family object and
# dispatches to the appropriate engine path.
# =============================================================================


#' Construct a tulpaObs family object
#'
#' Low-level constructor for `tulpa_obs_family` objects. End users should
#' call the specific family functions ([occ()], [nmixture()], etc.) rather
#' than this directly.
#'
#' @param name short name, e.g. `"occ"`.
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
#'
#' @return A `tulpa_obs_family` object.
#' @keywords internal
#' @export
obs_family <- function(name,
                       class_long,
                       latent,
                       observation,
                       replicates    = c("required", "optional", "single"),
                       default_engine = c("laplace", "nested_laplace", "nuts"),
                       status         = c("working", "planned", "experimental"),
                       params         = list()) {
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
      params         = params
    ),
    class = "tulpa_obs_family"
  )
}


# ---------------------------------------------------------------------------
# Working families — dispatch to existing tulpaObs builders
# ---------------------------------------------------------------------------

#' Single-season occupancy family
#'
#' Latent Bernoulli occupancy state with binomial detection per visit.
#' This is the standard MacKenzie et al. (2002) single-season model.
#'
#' @return A `tulpa_obs_family` object.
#' @export
#' @examples
#' f <- occ()
#' f
occ <- function() {
  obs_family(
    name           = "occ",
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
#' Bernoulli occupancy state with colonisation + extinction transitions
#' across seasons (the MacKenzie et al. dynamic model).
#'
#' @return A `tulpa_obs_family` object.
#' @export
dynamic_occ <- function() {
  obs_family(
    name           = "dynamic_occ",
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
#' Per-species occupancy with shared community-level hyperparameters on
#' the species random effects.
#'
#' @return A `tulpa_obs_family` object.
#' @export
multispecies_occ <- function() {
  obs_family(
    name           = "multispecies_occ",
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
#' Multiple data sources informing a shared latent occupancy state, each
#' with its own detection process.
#'
#' @return A `tulpa_obs_family` object.
#' @export
integrated_occ <- function() {
  obs_family(
    name           = "integrated_occ",
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
#' @return A `tulpa_obs_family` object.
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
#' @param family `"poisson"` or `"negbin"`.
#' @return A `tulpa_obs_family` object.
#' @export
nmixture <- function(K_max = NULL, family = c("poisson", "negbin")) {
  family <- match.arg(family)
  obs_family(
    name           = "nmixture",
    class_long     = "N-mixture abundance",
    latent         = if (family == "poisson") "poisson" else "negbin",
    observation    = "binomial_N",
    replicates     = "required",
    default_engine = "laplace",
    status         = "planned",
    params         = list(K_max = K_max, family = family)
  )
}


#' Multispecies N-mixture family
#'
#' Per-species N-mixture with shared community-level hyperparameters.
#'
#' @inheritParams nmixture
#' @return A `tulpa_obs_family` object.
#' @export
multispecies_nmix <- function(K_max = NULL, family = c("poisson", "negbin")) {
  family <- match.arg(family)
  obs_family(
    name           = "multispecies_nmix",
    class_long     = "multispecies N-mixture",
    latent         = if (family == "poisson") "poisson" else "negbin",
    observation    = "binomial_N",
    replicates     = "required",
    default_engine = "laplace",
    status         = "planned",
    params         = list(K_max = K_max, family = family)
  )
}


#' Open-population (Dail-Madsen) N-mixture family
#'
#' Latent N evolves across seasons via survival + recruitment.
#'
#' @inheritParams nmixture
#' @return A `tulpa_obs_family` object.
#' @export
dynamic_nmix <- function(K_max = NULL, family = c("poisson", "negbin")) {
  family <- match.arg(family)
  obs_family(
    name           = "dynamic_nmix",
    class_long     = "Dail-Madsen open N-mixture",
    latent         = "dail_madsen",
    observation    = "binomial_N",
    replicates     = "required",
    default_engine = "nuts",
    status         = "planned",
    params         = list(K_max = K_max, family = family)
  )
}


#' Distance-sampling family
#'
#' Latent density with hazard-rate or half-normal detection over distance
#' bins.
#'
#' @param key detection-function key. `"halfnorm"` or `"hazard"`.
#' @return A `tulpa_obs_family` object.
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
#' @return A `tulpa_obs_family` object.
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
#' Occupancy with both false negatives and false positives in the
#' detection process (e.g. acoustic classifiers, citizen-science misidentification).
#'
#' @return A `tulpa_obs_family` object.
#' @export
false_positive <- function() {
  obs_family(
    name           = "false_positive",
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
#' lognormal). Does not share the replicate-detection assumption of the
#' other families — see `vignette("families")` for the conceptual caveat.
#'
#' @param positive likelihood for the positive part. `"beta"` (cover in
#'   (0, 1)) or `"lognormal"` (log-cover Gaussian).
#' @return A `tulpa_obs_family` object.
#' @export
cover_hurdle <- function(positive = c("beta", "lognormal")) {
  positive <- match.arg(positive)
  # Phase 1a: lognormal-positive variant flips to "working" with the
  # two-Laplace dispatcher; the beta-positive variant stays "planned"
  # until Phase 1c/1d (joint multi-likelihood in tulpa_nested_laplace).
  is_working    <- positive == "lognormal"
  default_eng   <- if (is_working) "laplace" else "nested_laplace"
  status_value  <- if (is_working) "working" else "planned"
  obs_family(
    name           = "cover_hurdle",
    class_long     = "vegetation cover hurdle",
    latent         = "hurdle",
    observation    = if (positive == "beta") "binomial_plus_beta" else "binomial_plus_lognormal",
    replicates     = "single",
    default_engine = default_eng,
    status         = status_value,
    params         = list(positive = positive)
  )
}


# ---------------------------------------------------------------------------
# Print method
# ---------------------------------------------------------------------------

#' Print method for tulpa_obs_family
#' @param x a `tulpa_obs_family` object.
#' @param ... ignored.
#' @return `x`, invisibly.
#' @export
print.tulpa_obs_family <- function(x, ...) {
  cat(sprintf("<tulpa_obs_family: %s>\n", x$class_long))
  cat(sprintf("  name           : %s\n", x$name))
  cat(sprintf("  latent state   : %s\n", x$latent))
  cat(sprintf("  observation    : %s\n", x$observation))
  cat(sprintf("  replicates     : %s\n", x$replicates))
  cat(sprintf("  default engine : %s\n", x$default_engine))
  cat(sprintf("  status         : %s\n", x$status))
  if (length(x$params)) {
    cat("  params         :\n")
    for (nm in names(x$params)) {
      cat(sprintf("    %s = %s\n", nm, format(x$params[[nm]])))
    }
  }
  invisible(x)
}
