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

  # Laplace-EM: the prior / regularization scales. Only these are engine-level
  # values; the iteration budget is per-route (see the scope note above).
  #
  #   sigma.beta  coefficient ridge SD.
  #   sd.load     prior SD on a spatial-factor loading. A prior width like
  #               sigma.beta, and its copies must move together: the auto-K
  #               ladder selects a rank by marginal evidence under sd.load, so a
  #               value that drifted between the selection fit and the final fit
  #               would select a rank the fit does not use (gcol33/tulpaObs#189).
  #   re.lkj      LKJ shape regularizing a CORRELATED random slope's correlation
  #               in the AGHQ refine; pulls a weakly-identified correlation off
  #               the +-1 boundary without touching the marginal SDs. 1 disables
  #               it (uniform). No effect on intercept / uncorrelated terms.
  laplace = list(sigma.beta = 5, sd.load = 1.0, re.lkj = 1.5),

  # NUTS. One chain by default (the fits it backs are expensive and warm-started
  # at the Laplace mode), equal warmup and sampling, and the Stan-conventional
  # adaptation knobs. `n.iter` is POST-WARMUP draws kept per chain, so a run is
  # n.iter + n.warmup iterations long.
  nuts = list(n.iter = 1000L, n.warmup = 1000L, n.chains = 1L, n.thin = 1L,
              max.treedepth = 10L, adapt.delta = 0.9, seed = 1L,
              sigma.beta = 5, sigma.logr = 1.5),

  # Polya-Gamma Gibbs. A conjugate sweep costs far less than a NUTS trajectory,
  # so the chain is three times longer; two chains so split-Rhat is available
  # without a second call. The coefficient prior is tighter than the NUTS one --
  # a conjugate update has no step-size adaptation to absorb a wide prior. There
  # is deliberately no adapt.delta or max.treedepth: those are HMC knobs, and
  # their absence here is the answer to "what is the default adapt.delta for
  # pg_gibbs".
  #
  # CONVENTION, and it is the OPPOSITE of the nuts block above: here `n.iter` is
  # the TOTAL sweep count and warmup comes out of it, so the chain keeps
  # `n.iter - n.warmup` = 1500 draws. Every pg_gibbs fitter computes its kept
  # count as `length(seq.int(n.warmup + 1L, n.iter, by = n.thin))`; they agree
  # with each other, and the flip relative to NUTS is why it is stated here
  # rather than left to the reader (gcol33/tulpaObs#188).
  pg_gibbs = list(n.iter = 3000L, n.warmup = 1500L, n.chains = 2L,
                  n.thin = 1L, seed = 1L, sigma.beta = 2.5)
)


# The single-species NUTS entry, `.tobs_fit_model()`, which serves occupancy
# (occu / dyn_occu / int_occu) and every observation family (abun, removal,
# distance, fp_occu, dyn_abun, count). It forwards explicit sampler arguments
# down to each family fitter, so the family fitter's own formals never apply and
# THIS is the live answer for those families (gcol33/tulpaObs#188).
#
# The three knobs below predate the table and are kept on the record rather than
# aligned, because each is a calibrated value on a recovery-tested path and
# nothing shows it is wrong: a wider coefficient prior, a looser adaptation
# target, and a different stream seed. `n.iter` is deliberately NOT here -- see
# the note in `.tobs_fit_model()`.
.TOBS_SINGLE_SPECIES_NUTS <- list(sigma.beta = 10, adapt.delta = 0.8, seed = 42)

# The same entry's Laplace / nested-Laplace coefficient ridge. Wider than the
# community profile's 5 for the same reason the NUTS row is: it is the value
# every single-species recovery and coverage test was calibrated against, and
# the occupancy Laplace path is the most-exercised route in the package.
# Narrowing it is a change to those numbers, not a cleanup, so it is recorded.
.TOBS_SINGLE_SPECIES_LAPLACE <- list(sigma.beta = 10)


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
    ms_occu_cover_spatial = list(n.warmup = 500L, adapt.delta = 0.95),

    # The single-species occu_cover sampler with a coupled areal field, whose
    # field hyperparameters are sampled (gcol33/tulpaObs#204). A proper-CAR
    # precision Q(rho) = D - rho W approaches the INTRINSIC (rank-deficient)
    # limit as rho -> 1, and the data on an ICAR-simulated field pushes rho
    # there, so the field's near-null direction stretches against the psi
    # intercept -- a funnel the profile's step size cannot walk. Measured over
    # four seeds at N = 64, the divergence count falls 4 / 1 / 16 / 13 (0.80) ->
    # 6 / 0 / 2 / 0 (0.95) -> 0 / 0 / 0 / 0 (0.99) while the posterior does not
    # move (rho mean identical to three decimals, field_sd within Monte Carlo
    # noise) -- a step-size artifact, not a region the chain was missing. The
    # cost is roughly 2x wall time, which a reference posterior is worth.
    occu_cover_spatial = list(adapt.delta = 0.99)
  )
)


# `n.quad` -- one control name, five routes, and deliberately NOT one number
# (gcol33/tulpaObs#189). Every route below integrates a different marginal over
# a different latent dimension, so the node count that suffices differs by an
# order of magnitude; what was wrong was not the spread but that no reader could
# find out which number applied, while `?tobs` stated a single default that most
# routes do not use.
#
# `.tobs_n_quad(route)` is the accessor; each entry names the marginal.

# Minimum Gauss-Hermite order on a scalar nuisance random-effect block (the
# negative-binomial dispersion `log_r`, the zero-inflation `logit_omega`). A
# 2-node rule places two nodes and has no freedom left to represent curvature:
# where the 1-D posterior is not near-Gaussian, the marginal it returns can come
# out arbitrarily sharp and the reported community-mean SE collapses with it
# (gcol33/tulpaObs#234: a 17x collapse on `mu_log_r`, with the fit converged and
# the point estimate ordinary). The rule is converged at 3 -- on the `log_r`
# block, 3 / 5 / 9 nodes agree to ten decimal places -- so 3 is both the order
# below which the integrand is not represented and the order above which nothing
# on that axis moves. `nmix_laplace_re()` applies it as a floor, not only as a
# default, because nothing downstream can detect the collapse.
.TOBS_MIN_SCALAR_NQUAD <- 3L

.TOBS_NQUAD_ROUTES <- list(
  # Formula random effect under method = "laplace": adaptive Gauss-Hermite over
  # the exact per-group marginal, debiasing the Laplace small-cluster
  # attenuation of a binary occupancy variance component. Binary data carries
  # little information per group, so the refine wants many nodes.
  re_aghq = 9L,

  # Community N-mixture (`ms_abun()`) variance-component debias. 1 = the plain
  # Laplace (nAGQ = 1) marginal, i.e. the EM default: each species' count
  # marginal is already informative, so the AGHQ refine barely moves the
  # community covariances and is opt-in (via optimizer = "joint_fd") for the
  # sparse / rare-species regime.
  ms_nmix = 1L,

  # The trailing scalar coordinates of that same community fit (a per-species
  # log-dispersion under NB, a per-species structural-zero logit under ZI),
  # integrated separately. Both the default and the floor, since the rule is
  # converged here and unusable below.
  ms_nmix_scalar = .TOBS_MIN_SCALAR_NQUAD,

  # Community joint occupancy-cover (`ms_occu_cover()`): tensor AGHQ over the
  # joint per-species RE vector, so the node count is raised to a power of the
  # RE dimension and stays small.
  ms_occu_cover = 5L,

  # Latent cover-per-unit (`cover_aggregate = "latent"`) on the joint engine.
  # The lognormal per-unit marginal is CLOSED FORM, so it needs no quadrature at
  # all; the beta one is not, and is integrated numerically.
  cover_latent_lognormal = 1L,
  cover_latent_beta      = 15L,

  # Community latent() factors: Gauss-Hermite nodes for the joint site marginal
  # the factor scores integrate on. Both the loading magnitude and the
  # score-matched offset are insensitive to this (argmax stable to < 0.4% vs 21
  # nodes), so it stays small.
  community_latent = 5L
)

.tobs_n_quad <- function(route) {
  if (!route %in% names(.TOBS_NQUAD_ROUTES)) {
    stop(sprintf("no n.quad default for route \"%s\".", route), call. = FALSE)
  }
  .TOBS_NQUAD_ROUTES[[route]]
}

# Resolve the AGHQ order of a scalar nuisance random-effect block from the
# requested coefficient-block order and the requested scalar order. The scalar
# blocks track the coarser of the two, then take .TOBS_MIN_SCALAR_NQUAD as a hard
# floor -- including over `n_quad` itself, which may legitimately sit at the
# plain-Laplace order while a block whose marginal SE is reported may not.
.nmix_scalar_nquad <- function(n_quad, n_quad_scalar) {
  max(.TOBS_MIN_SCALAR_NQUAD,
      min(as.integer(n_quad), as.integer(n_quad_scalar)))
}


# The resolved profile for the single-species entry `.tobs_fit_model()`. Its
# departures belong to the ENTRY, not to the nine families that pass through it,
# so they are one override applied by engine rather than nine identical family
# rows. `pg_gibbs` deliberately has none: `occu(method = "pg_gibbs")` used to
# keep 1000 draws where every `ms_*` sibling kept 1500, which is verbatim the
# failure this table exists to prevent, and no reason for the shorter chain was
# ever recorded (gcol33/tulpaObs#188). It now reads the shared profile.
.tobs_single_species_defaults <- function(engine) {
  base <- .TOBS_ENGINE_DEFAULTS[[engine]] %||% list()
  ov <- switch(engine,
               nuts           = .TOBS_SINGLE_SPECIES_NUTS,
               laplace        = .TOBS_SINGLE_SPECIES_LAPLACE,
               nested_laplace = .TOBS_SINGLE_SPECIES_LAPLACE,
               NULL)
  if (is.null(ov)) base else utils::modifyList(base, ov)
}


# The resolved profile for one (engine, family). An engine with no profile
# returns an empty list, so a caller can resolve unconditionally.
.tobs_engine_defaults <- function(engine, family = NULL) {
  base <- .TOBS_ENGINE_DEFAULTS[[engine]] %||% list()
  ov <- if (is.null(family)) NULL else .TOBS_FAMILY_DEFAULTS[[engine]][[family]]
  # An engine may carry family rows and no shared profile at all
  # (`nested_laplace` is not a sampler and has only the ridge), so the family
  # row is applied even when the base is empty.
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
# Fill the sampler knobs a fitter was not given from the (engine, family)
# profile. A fitter declares those knobs as `NULL` formals and calls this as its
# first statement, so the table is the ONLY answer to "what is the default
# n.iter here" -- a literal formal alongside the table is a second answer, and
# for the fitters reached through `.tobs_fit_model()` it was an unreachable one
# (that entry forwards explicit values, so the formal never applied) while for
# the rest it was the live one and disagreed (gcol33/tulpaObs#188).
#
# Only knobs the fitter actually declares are touched, so one call serves
# fitters with different knob sets. `single_species = TRUE` selects the
# `.tobs_fit_model()` entry's profile.
.tobs_fill_sampler <- function(env, engine, family = NULL,
                               single_species = FALSE) {
  d <- if (isTRUE(single_species)) .tobs_single_species_defaults(engine)
       else .tobs_engine_defaults(engine, family)
  for (k in names(d)) {
    if (!exists(k, envir = env, inherits = FALSE)) next
    if (is.null(get(k, envir = env, inherits = FALSE))) {
      assign(k, d[[k]], envir = env)
    }
  }
  invisible(NULL)
}


.tobs_control_defaults <- function(control, engine, family = NULL) {
  if (is.null(control)) control <- list()
  d <- .tobs_engine_defaults(engine, family)
  for (k in names(d)) if (is.null(control[[k]])) control[[k]] <- d[[k]]
  control
}
