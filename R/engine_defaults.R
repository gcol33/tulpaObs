# =============================================================================
# engine_defaults.R - the one table of per-engine control defaults.
#
# The sampler knobs (chain length, warmup, chain count, thinning, seed, the HMC
# adaptation knobs, and the coefficient / dispersion prior scales) have defaults
# that are a property of the FITTING ENGINE rather than of the family being fit,
# and they used to be written out at every dispatcher branch that reads them --
# six sites for the Polya-Gamma Gibbs profile alone. Changing one meant editing
# all six, and missing one produced a family whose chain was shorter than its
# siblings with no test failing (gcol33/tulpaObs#183).
#
# Scope. This table covers the SAMPLER knobs only. `max.iter` and `tol` are
# deliberately NOT here: they are Laplace-EM knobs, they are read inside the
# NUTS branches too (for the EM warm start that precedes sampling), and their
# value genuinely differs by route within one family -- ms_occu_cover() iterates
# its own EM 30 times at tol 1e-3 but warm-starts its sampler with 200 at 1e-4,
# and ms_occu()'s plain areal C++ EM uses 100 where its block-coordinate latent
# fitter uses 200. Those are per-route values, not a profile, so they stay at
# their call sites where the route that owns them is visible.
# =============================================================================


# Per-engine profile. A knob listed here has the same default for every family
# that does not carry an explicit row in .TOBS_FAMILY_DEFAULTS below.
.TOBS_ENGINE_DEFAULTS <- list(

  # Laplace-EM: the coefficient ridge. Only the prior scale is an engine-level
  # value; the iteration budget is per-route (see the scope note above).
  laplace = list(sigma.beta = 5),

  # NUTS. One chain by default (the fits it backs are expensive and warm-started
  # at the Laplace mode), equal warmup and sampling, and the Stan-conventional
  # adaptation knobs. `n.iter` is POST-WARMUP draws kept per chain, so a run is
  # n.iter + n.warmup iterations long.
  nuts = list(n.iter = 1000L, n.warmup = 1000L, n.chains = 1L,
              max.treedepth = 10L, adapt.delta = 0.9, seed = 1L,
              sigma.beta = 5, sigma.logr = 1.5),

  # Polya-Gamma Gibbs. A conjugate sweep costs far less than a NUTS trajectory,
  # so the chain is three times longer; two chains so split-Rhat is available
  # without a second call. The coefficient prior is tighter than the NUTS one --
  # a conjugate update has no step-size adaptation to absorb a wide prior. There
  # is deliberately no adapt.delta or max.treedepth: those are HMC knobs, and
  # their absence here is the answer to "what is the default adapt.delta for
  # pg_gibbs".
  pg_gibbs = list(n.iter = 3000L, n.warmup = 1500L, n.chains = 2L,
                  n.thin = 1L, seed = 1L, sigma.beta = 2.5)
)


# Deliberate per-(engine, family) departures from the profile above. A row here
# is a decision on the record; a bare literal at a call site was not
# distinguishable from a typo.
.TOBS_FAMILY_DEFAULTS <- list(

  nuts = list(
    # Count / abundance families sample their coefficients on a log link, where
    # a unit change is multiplicative, so they carry a wider coefficient prior
    # than the logit-link occupancy families. This mirrors the C++ model
    # defaults exactly (MsCountNutsModel / MsNmixNutsModel sigma_beta = 10 vs
    # MsOccuNutsModel / MsIntOccuNutsModel / MsOccuCoverNutsModel = 5).
    ms_count = list(sigma.beta = 10),
    jsdm     = list(sigma.beta = 10),
    ms_abun  = list(sigma.beta = 10),

    # The spatial-factor community occu_cover sampler carries a coupled field
    # and per-species loadings, a harder geometry than the non-spatial targets:
    # a tighter adaptation target and a shorter warmup than the shared profile.
    ms_occu_cover_spatial = list(n.warmup = 500L, adapt.delta = 0.95)
  )
)


# The resolved profile for one (engine, family). An engine with no profile
# returns an empty list, so a caller can resolve unconditionally.
.tobs_engine_defaults <- function(engine, family = NULL) {
  base <- .TOBS_ENGINE_DEFAULTS[[engine]]
  if (is.null(base)) return(list())
  ov <- if (is.null(family)) NULL else .TOBS_FAMILY_DEFAULTS[[engine]][[family]]
  if (is.null(ov)) base else utils::modifyList(base, ov)
}


# One knob of the (engine, family) profile. For a branch that reads a single
# knob inline, where resolving the whole profile would say less than the one
# lookup does. Errors on an unknown knob rather than returning NULL, so a
# renamed row cannot silently become "no default".
.tobs_default <- function(engine, knob, family = NULL) {
  d <- .tobs_engine_defaults(engine, family)
  if (!knob %in% names(d)) {
    stop(sprintf("no default for control$%s under engine \"%s\"%s.", knob, engine,
                 if (is.null(family)) "" else sprintf(" (family \"%s\")", family)),
         call. = FALSE)
  }
  d[[knob]]
}


# Fill every knob of the (engine, family) profile the caller did not set. A user
# value always wins, and an absent or NULL entry counts as unset -- the same
# rule the `control[["knob"]] %||% default` sites used to apply one at a time.
.tobs_control_defaults <- function(control, engine, family = NULL) {
  if (is.null(control)) control <- list()
  d <- .tobs_engine_defaults(engine, family)
  for (k in names(d)) if (is.null(control[[k]])) control[[k]] <- d[[k]]
  control
}
