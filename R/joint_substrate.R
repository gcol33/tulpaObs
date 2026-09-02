# =============================================================================
# joint_substrate.R -- shared access layer for the joint nested-Laplace fits
# behind cover() (2-arm hurdle) and occu_cover() (3-arm occupancy-cover).
#
# Both families fit a `tulpa_nested_laplace_joint` object carrying the same
# substrate: per-grid weights / modes / sparse precision, an `arm_layout`, a
# shared latent field, and a working `tulpa::tulpa_posterior_draws()`. They
# differ only in (a) which slot holds the joint object, (b) the arm roster
# (occ/pos vs psi/p/pos), and (c) the field-amplitude convention -- cover()'s
# independent (sigma_occ, sigma_pos) reparam versus occu_cover()'s donor sigma
# scaled onto the cover arm by the copy coefficient alpha. This file normalizes
# those three differences ONCE into a family-agnostic draw bundle so predict()
# and the pointwise log-likelihood share a single draws -> linear-predictor
# path.
# =============================================================================


# Default outer grids for a coupled joint fit, shared by every family that
# drives `tulpa_nested_laplace_joint()` (cover, occu_cover, occu_multiscale_cover,
# the standalone occu SVC path). Single source of truth: each fitter reads them
# through the `control$*.grid` override, so the engine default and the front-door
# default cannot drift apart.
#
# The copy coefficient alpha scales the shared field onto a copied arm; 0 is on
# the grid so an uncoupled arm is reachable exactly. sigma is the field
# amplitude. bym2 additionally mixes structured and unstructured components at
# rho. Declare an outer-grid axis this package defaulted rather than one the
# user pinned. The engine's auto-recenter decides axis
# PROVENANCE, not field presence: a marked axis (or one whose nodes are exactly
# the engine's own default) may be recentred onto the hyperparameter mode,
# anything else is treated as a deliberate pin and is never moved. Because
# tulpaObs writes a grid on every joint fit -- it derives a second axis from
# the first, and hands the same vector to several blocks -- an unmarked default
# is indistinguishable from a pin and silently makes the rescue inert.
#
# The marker is an attribute, so `sort()` / `[` / `c()` / `as.numeric()` /
# `expand.grid()` all drop it: call this LAST, on the vector actually written
# into the block, after any sorting or `1 / sqrt(tau)` translation.
#
# `auto` is the "the user named nothing here" test, so a site reads
# `.tobs_mark_auto(<grid>, is.null(control$sigma.grid))`.
.tobs_mark_auto <- function(x, auto) {
  if (isTRUE(auto)) tulpa::auto_grid(x) else x
}

# `as.numeric()` on a grid whose provenance is already settled: coerce, then
# carry the source vector's own marker onto the coerced copy. Every site that
# writes `as.numeric(g)` into a block after resolving `control$*.grid %||%
# <default>()` reads this instead, so the coerce-and-lose-the-marker step has
# one spelling.
.tobs_num_auto <- function(x) {
  .tobs_mark_auto(as.numeric(x), tulpa::is_auto_grid(x))
}

#
# Each returns its vector already marked with `auto_grid()`: these functions ARE
# the layer that defaults the axis, and a caller only ever reaches them on the
# `control$*.grid %||%` fallback, so a user-supplied grid never passes through
# and never picks up the mark. A site that reshapes the result afterwards
# (`sort()`, `1 / sigma^2`, `expand.grid()`) drops the attribute and re-applies
# it with `.tobs_mark_auto()`.
#
# The copy-coefficient and field-SD axes are READ from the engine rather than
# restated here. Both fields are bound in the engine's `.NL_FAMILY_AXES` on the
# very paths these grids reach -- `alpha_grid` and `sigma_grid` on `.copy`,
# `sigma_grid` on `.joint_areal` -- so `.nl_axis_matches_default()` recognises
# a grid carrying those nodes as a default by VALUE, independently of the
# `auto_grid()` marker. Reading holds that recognition in step permanently: it
# is the belt that still classifies the axis as ours after a reshape has
# dropped the marker, which is the failure filed, and a restated copy loses the
# belt the moment the engine moves the nodes the way moved `bym2_rho`.
.tobs_default_alpha_grid <- function(n = NULL) {
  ax <- tulpa:::.nl_grid_axis("copy_alpha", n = n)
  # `n` is the caller naming the axis, so only the unresolved read carries the
  # "nothing was named here" mark.
  if (is.null(n)) tulpa::auto_grid(ax) else ax
}

# ---- the copy coefficient's outer axis ------------------------------------
#
# Two requests, one axis. NODES state the axis's coordinates, and with them the
# prior structure it carries: the atom at alpha = 0 that gives the no-coupling
# model posterior mass, and the log-spaced slab above it. A RESOLUTION states
# how finely the engine's OWN axis is read -- `n` slab nodes, atom and slab
# bounds unchanged, so the axis comes back `n + 1` nodes long and nothing about
# its structure is restated to sharpen it. One block takes one of them, and a
# block given both is refused.
#
# Both are written in the formula, on the copy whose coefficient they place:
# `share(spatial(), alpha = grid(c(...)))` states nodes, `alpha = grid(n = 9)` a
# resolution, `terms = list(<component> = ...)` either, per block, and
# `prior =` regularizes the coefficient itself.
# `control$alpha.grid[.trend]` / `alpha.n[.trend]` / `prior.alpha` are the
# lower-level spelling of the same requests, and the representation the formula
# compiles into; a fit writing both is refused in the dispatcher, where the
# knobs still carry the names the user typed (`.tobs_check_alpha_control()`).
#
# The resolution is the only way to raise this axis: it does not densify when
# the donor `sigma.grid` does, so on informative data the outer grid's
# quadrature effective sample size saturates on the copy amplitude while every
# other axis tracks the request (`NOTES_measurements.md`).
#
# A grid-integrated fit hands the ENGINE the axis rather than its nodes, so a
# fit raising the resolution never restates the structure: `alpha_n` is resolved
# against the engine's own declaration. `.tobs_alpha_axis()` returns the two
# fields a copy spec and a `field_coef` carry, exactly one of them non-NULL.
# `.tobs_default_alpha_grid(n)` is the read for the one consumer that has to
# hold the nodes itself -- the NUTS warm fit, whose sampled alpha takes the
# axis's realised span as the support of its flat prior.
#
# Stated nodes win over a resolution here: a field block with no `share()` is
# pinned (`grid = 0`, decoupled) and reaches this with the fit's `n` alongside,
# and a pinned block has no axis to resolve.
.tobs_alpha_n <- function(n) {
  n <- suppressWarnings(as.integer(n))
  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop("A copy amplitude resolution (share(alpha = grid(n = )), ",
         "control$alpha.n[.trend]) must be a single integer >= 1: it is the ",
         "number of slab nodes on the copy coefficient's axis (the atom at ",
         "alpha = 0 is carried alongside them).", call. = FALSE)
  }
  n
}

.tobs_alpha_axis <- function(grid = NULL, n = NULL) {
  if (!is.null(grid)) {
    return(list(alpha_grid = .tobs_num_auto(grid), alpha_n = NULL))
  }
  if (!is.null(n)) {
    return(list(alpha_grid = NULL, alpha_n = .tobs_alpha_n(n)))
  }
  list(alpha_grid = .tobs_num_auto(.tobs_default_alpha_grid()), alpha_n = NULL)
}

# The axis a fit's intercept field block reads off `control`, and the axis its
# weighted-trend block reads: either trend knob replaces the base axis for that
# block, and with neither the trend block rides the base axis.
#
# `[[` (exact), never `$`: `$` prefix-matches on a list, so `control$alpha.n`
# resolves to `alpha.n.trend` on a fit that sets only the trend knob.
.tobs_alpha_axis_base <- function(control) {
  .tobs_alpha_axis(control[["alpha.grid"]], control[["alpha.n"]])
}

.tobs_alpha_axis_trend <- function(control, base) {
  g <- control[["alpha.grid.trend"]]
  n <- control[["alpha.n.trend"]]
  if (is.null(g) && is.null(n)) base else .tobs_alpha_axis(g, n)
}

# Write one resolved amplitude (`.tobs_copy_amp()`) onto the pair of control
# keys its block reads. Only what the copy STATED is written: an amplitude
# stating neither nodes nor a resolution asks for the engine's default axis,
# which is exactly what composes with a `control$alpha.n[.trend]` the fit set
# alongside a bare `share(spatial())`, so writing NULL over that key (`[[<-` with
# NULL drops a list element) would discard the resolution the fit asked for.
# Nothing can leave both keys set: a share() that states either is refused
# alongside both control spellings (`.tobs_check_alpha_copy()`,
# `has_control_alpha`), and a decoupled block's stated `grid = 0` takes
# precedence over a resolution in `.tobs_alpha_axis()` -- a pinned block has no
# axis to resolve.
.tobs_alpha_control_set <- function(control, amp, grid_key, n_key) {
  if (!is.null(amp$grid)) control[[grid_key]] <- as.numeric(amp$grid)
  if (!is.null(amp$n))    control[[n_key]]    <- .tobs_alpha_n(amp$n)
  control
}

# One hyperprior on the copy coefficient reaches the engine per fit, and on a
# multi-block grid it is baked onto the FIRST block carrying a (sigma, alpha)
# axis pair (tulpa's `.joint_multi_hp_cols()`). A fit copying several blocks
# would have one block regularized and the rest not, which is neither spelling's
# meaning, so it is refused rather than applied to block 1 (gcol33/tulpa#655).
.tobs_check_alpha_prior <- function(prior, n_copied, what) {
  if (is.null(prior) || n_copied <= 1L) return(invisible(TRUE))
  stop(what, ": a prior on the copy coefficient reaches one copied block, and ",
       "this fit copies ", n_copied, " (the field's intercept and its weighted ",
       "trend blocks each carry their own amplitude). Drop the prior, or copy ",
       "one block -- terms = list(...) decouples a block by giving it ",
       "alpha = 0.", call. = FALSE)
}

# The two engine-facing shapes of one resolved axis: one copy spec's fields on
# the multi-block driver, and the pos arm's `field_coef` on the single-block
# path.
# Is a resolved amplitude axis the DECOUPLED sentinel -- the single node 0?
#
# `alpha` is the amplitude of a copy, so it is a parameter only when there IS a
# copy. `.occu_cover_apply_copy_coupling()` spells "no share() named this block"
# as an amplitude axis stating just 0, which is the same model as no coupling at
# all but still reaches the engine as a real outer axis.
#
# What the two engine paths then do with that differs, and the split is forced
# by the engine, not chosen:
#
#   * SINGLE-block backend -- the cover arm takes a plain numeric
#     `field_coef = 0` (the spelling the detection arm already uses) and no
#     amplitude axis is built at all. That backend integrates `sigma_grid` under
#     either coupling, so the field's prior is untouched.
#   * MULTI-block driver -- the copy spec STAYS, pinned at 0. That driver offers
#     the SD parameterization only to a COPIED block; a non-copied one rides its
#     precision `b<k>.tau` instead. Same implied SDs, different prior MEASURE on
#     the field's scale, and the shift is measurable (`NOTES_measurements.md`).
#     Dropping the copy to tidy the reporting would re-prior the shared field.
#
# Either way `alpha` must not be REPORTED when there is no copy: a pinned axis
# was never estimated, and emitting it put an all-zero row and column into
# vcov() (leaving it singular) and one too many in n_params.
.tobs_alpha_axis_decoupled <- function(axis) {
  if (is.null(axis)) return(FALSE)
  .tobs_copy_amp_decoupled(.tobs_copy_amp(grid = axis$alpha_grid,
                                          n = axis$alpha_n))
}

# The precision axis matching a field-SD axis. A block the engine does NOT copy
# is parameterized by its precision `b<k>.tau` rather than by `b<k>.sigma` +
# `b<k>.alpha`, so a decoupled block states its grid as tau = 1 / sigma^2. The
# auto-grid mark rides across the translation: the provenance is the SD axis the
# caller supplied, whether that was the default or stated.
.tobs_sigma_to_tau_grid <- function(sigma_grid) {
  .tobs_mark_auto(sort(1.0 / as.numeric(sigma_grid)^2),
                  tulpa::is_auto_grid(sigma_grid))
}

.tobs_alpha_copy_spec <- function(arm, block, axis) {
  list(arm = arm, block = as.integer(block),
       alpha_grid = axis$alpha_grid, alpha_n = axis$alpha_n)
}

.tobs_alpha_field_coef <- function(axis) {
  list(name = "alpha", grid = axis$alpha_grid, n = axis$alpha_n)
}

# A resolved axis's realised NODES. The grid-integrated routes hand the engine
# the axis and never need these; the NUTS warm fit does, because the sampled
# alpha takes the node set's span as the support of its flat prior and anchors
# its declared slab on the largest positive node.
.tobs_alpha_nodes <- function(axis) {
  axis$alpha_grid %||% .tobs_default_alpha_grid(axis$alpha_n)
}

# Does this `share()` STATE an amplitude -- nodes or a resolution? `alpha =
# grid(c(...))`, `alpha = grid(n = )` and a bare scalar `alpha =` all do (the
# last pins a one-node axis); `share(spatial())` with no amplitude does not, and
# asks for the engine's default axis.
.tobs_copy_states_amplitude <- function(cp) {
  if (!is.null(cp$copy_terms)) {
    return(any(vapply(cp$copy_terms,
                      function(res) !isTRUE(is.na(res$integrate)), logical(1))))
  }
  !isTRUE(is.na(cp$alpha_integrate))
}

# A `share()` that states an amplitude and `control$alpha.n[.trend]` are the same
# axis written twice; a `share()` that states none composes with the resolution
# knob.
.tobs_check_alpha_copy <- function(states_amplitude, control, what) {
  n_keys <- c("alpha.n", "alpha.n.trend")
  n_set  <- n_keys[vapply(n_keys, function(k) !is.null(control[[k]]), logical(1))]
  if (!isTRUE(states_amplitude) || length(n_set) == 0L) return(invisible(TRUE))
  stop(what, ": share(alpha = ) states the copy axis and control$", n_set[1L],
       " states how many nodes the engine's own axis is read at. Give one: ",
       "write the resolution in the formula as share(alpha = grid(n = ",
       control[[n_set[1L]]], ")) and drop control$", n_set[1L],
       ", or drop the amplitude from share().", call. = FALSE)
}

# `alpha.grid` states the axis's nodes, `alpha.n` a resolution for the engine's
# own axis; one block takes one of them. Refused here, in the dispatcher, so the
# message names the knobs as the user spelled them.
# The coupling-amplitude keys are the WIRE FORMAT the formula compiles into, not
# user surface (#295). Two front doors for one request meant a guard layer to
# keep them consistent, an asymmetry in which one composes with what, and two
# spellings a reader has to know are the same model. `share()` is the one that
# says what the model IS; these say how the engine reads an axis.
#
# Refused at the dispatcher entry, BEFORE any translation writes them, so the
# compile target still works and only the user's own spelling is rejected.
.TOBS_WIRE_ALPHA_KEYS <- c(
  alpha.grid       = "share(spatial(), alpha = grid(c(...)))",
  alpha.n          = "share(spatial(), alpha = grid(n = <k>))",
  alpha.grid.trend = "share(spatial(), terms = list(trend = grid(c(...))))",
  alpha.n.trend    = "share(spatial(), terms = list(trend = grid(n = <k>)))")

.tobs_check_alpha_control <- function(control, what) {
  hit <- intersect(names(.TOBS_WIRE_ALPHA_KEYS), names(control))
  if (!length(hit)) return(invisible(TRUE))
  stop(sprintf(paste0(
    "%s: control$%s %s not user surface -- %s the wire format the formula ",
    "compiles into. State the coupling where the model is declared:\n  %s\n",
    "A field the formula does not couple is decoupled (alpha pinned at 0); to ",
    "give one block its own amplitude, name every block with ",
    "share(spatial(), terms = list(intercept = ..., trend = ...))."),
    what, paste(hit, collapse = ", control$"),
    if (length(hit) > 1L) "are" else "is",
    if (length(hit) > 1L) "they are" else "it is",
    paste(sprintf("control$%-16s ->  %s", hit,
                  .TOBS_WIRE_ALPHA_KEYS[hit]), collapse = "\n  ")),
    call. = FALSE)
}

# ---- outer-grid adaptive knobs -------------------------------------------
#
# TWO independent mechanisms share the `adaptive.grid` prefix, and a fit can
# use either without the other:
#
#   * the post-integration REFINEMENT passes (`adaptive.grid`,
#     `.edge.thresh`, `.max.passes`), which densify a coarse axis around the
#     mode after the grid has been integrated; and
#   * the `integration = "grid_adaptive"` LATTICE BUILDER (`.cutoff`,
#     `.stride`, `.max.frac`, `.min.cells`), which evaluates a strict subset
#     of the same tensor lattice and declines back to the dense tensor when
#     the kept region would rival it.
#
# The builder knobs are forwarded UNSET so the engine keeps the single
# definition of their defaults; restating them here would fork that number.
#
# Exact `[[`, never `$`: every name in this family has `adaptive.grid` as a
# prefix, so a `$` read of the master flag returns a SUB-KNOB's value whenever
# exactly one sub-knob is set and the flag itself is not -- a fit passing only
# `adaptive.grid.edge.thresh = 0.05` would read `adaptive_grid = 0.05`.
.tobs_adaptive_grid_control <- function(control) {
  list(
    adaptive_grid             = control[["adaptive.grid"]]             %||% TRUE,
    adaptive_grid_edge_thresh = control[["adaptive.grid.edge.thresh"]] %||% 0.02,
    adaptive_grid_max_passes  = control[["adaptive.grid.max.passes"]]  %||% 1L,
    adaptive_grid_cutoff      = control[["adaptive.grid.cutoff"]],
    adaptive_grid_stride      = control[["adaptive.grid.stride"]],
    adaptive_grid_max_frac    = control[["adaptive.grid.max.frac"]],
    adaptive_grid_min_cells   = control[["adaptive.grid.min.cells"]])
}

.tobs_default_sigma_grid <- function() {
  tulpa::auto_grid(tulpa:::.nl_grid_axis("field_sd"))
}

.tobs_default_bym2_rho_grid <- function() {
  tulpa::auto_grid(c(0.25, 0.5, 0.75))
}

# Spatial-correlation axis of a proper-CAR field, shared by the cover()
# arm-specific block and the occu_cover() NUTS warm fit.
#
# These nodes are deliberately tulpaObs's own, NOT a read of the engine's
# `joint_car_rho`, which today holds the same four values. Two reasons, and they
# point the same way. The engine binds neither this axis nor the `rho_car_grid`
# field it is written to into `.NL_FAMILY_AXES`, by its own statement, precisely
# so a caller passing these nodes stays classified as a caller rather than as
# the engine's default -- so unlike the alpha and field-SD axes above there is
# no value-recognition belt to hold in step, and the equality is coincidence
# rather than contract. And the axis does not exist at this package's declared
# `Imports` floor (tulpa 0.0.136 errors on the key), so a read would have to
# fall back to these nodes anyway, making the default grid a function of which
# engine happens to be installed within one declared floor. Provenance is
# carried by the `auto_grid()` marker, which the rescue reads before it ever
# reaches the value comparison.
.tobs_default_rho_car_grid <- function() {
  tulpa::auto_grid(c(0.5, 0.8, 0.95, 0.99))
}

# Non-spatial block axes on the cover() multi-block path. The temporal
# precision axis serves ar1 and rw1 / rw2 alike, so it is one function rather
# than the same three nodes written per block type.
.tobs_default_temporal_tau_grid <- function() {
  tulpa::auto_grid(c(1, 4, 16))
}

.tobs_default_temporal_rho_grid <- function() {
  tulpa::auto_grid(c(0.3, 0.7))
}

.tobs_default_temporal_sigma_grid <- function() {
  tulpa::auto_grid(exp(seq(log(0.1), log(1), length.out = 3)))
}

.tobs_default_re_sigma_grid <- function() {
  tulpa::auto_grid(exp(seq(log(0.1), log(1.5), length.out = 3)))
}

# Field-SD axis of the standalone occu() joint path, which grids each ICAR
# block on its own precision. Four nodes rather than the coupled paths' five:
# the base outer grid is 4^n_fields cells there (16 for an intercept plus one
# trend field), and adaptive refinement densifies where the posterior piles up.
.tobs_default_occu_joint_sigma_grid <- function() {
  tulpa::auto_grid(exp(seq(log(0.15), log(3), length.out = 4)))
}


# The joint nested-Laplace object regardless of family slot: occu_cover() stores
# it at `$joint_fit`, cover() at `$joint`. NULL when neither is present (a
# non-spatial / separate-Laplace fit that carries no joint object).
.tobs_joint_fit <- function(object) {
  object$joint_fit %||% object$joint
}

# Promote the outer Pareto-k diagnostic from a joint nested-Laplace result to the
# tobs_fit top level. The joint engine attaches `pareto_k`, `pareto_k_is_ess`,
# `pareto_k_scope`, and `pareto_k_proposal_source` (the
# mode-Hessian-vs-grid-moment proposal source) to the raw result the postprocess
# wrappers nest at `$joint_fit`. A user diagnostic should read `fit$pareto_k` /
# `fit$pareto_k_proposal_source` directly rather than reach into `$joint_fit`, so
# each joint-coupled family splices the result of this into its return list.
# Returns a named list of the fields actually present (so `c(list(...),
# .tobs_promote_pareto_k(jf), list(...))` adds exactly those), or NULL when the
# joint result carries none -- diagnose.k defaults OFF, so the fields are
# NA-or-absent and this is inert unless the diagnostic was requested.
.tobs_promote_pareto_k <- function(jf) {
  if (!is.list(jf)) return(NULL)
  keys <- c("pareto_k", "pareto_k_is_ess", "pareto_k_scope",
            "pareto_k_proposal_source")
  present <- keys[keys %in% names(jf)]
  if (length(present) == 0L) return(NULL)
  # Inert unless the diagnostic actually ran: with diagnose.k OFF the engine
  # sets every field to NA and returns, so there is nothing to surface. A
  # diagnostic that ran reports a finite k-hat and its finite IS-ESS
  # (pareto_k_is_ess), so gate on a usable number rather than mirroring the
  # always-present NA placeholders.
  k    <- jf$pareto_k
  ess  <- jf$pareto_k_is_ess
  ran  <- (length(k) == 1L && is.finite(k)) ||
          (length(ess) == 1L && !is.na(ess))
  if (!ran) return(NULL)
  jf[present]
}

# Promote the outer-grid placement record from a joint nested-Laplace result to
# the tobs_fit top level. The engine reports where the outer grid ended up --
# `outer_grid_placement` ("fixed" / "auto_recentered"),
# `outer_grid_recenter_attempts`, `outer_grid_prior_added`, and
# `outer_grid_recenter_declined`, the reason a "fixed" placement stayed fixed.
# Without this promotion a caller has to reach into `$joint_fit` / `$joint` to
# learn whether the auto grid did anything, which is how an inert recenter stayed
# invisible across a whole batch. Same contract as `.tobs_promote_pareto_k()`:
# returns a named list of the fields actually present so the result splices with
# `c(list(...), .tobs_promote_outer_grid(jf), list(...))`, and NULL when the
# result carries none (a fit from an engine predating the record, or a non-joint
# path).
#
# NOT gated on the placement being "auto_recentered": a fixed placement with its
# decline reason is precisely the case worth surfacing.
.tobs_promote_outer_grid <- function(jf) {
  if (!is.list(jf)) return(NULL)
  keys <- c("outer_grid_placement", "outer_grid_recenter_attempts",
            "outer_grid_prior_added", "outer_grid_recenter_declined")
  present <- keys[keys %in% names(jf)]
  if (length(present) == 0L) return(NULL)
  jf[present]
}

# Add the outer-grid placement columns to a one-row glance data frame. Shared by
# `glance.tobs_fit()` and `glance.tobs_multiarm_fit()` -- the latter is terminal
# for `cover_fit` (class order cover_fit / tobs_multiarm_fit / tobs_fit), so it
# never reaches the former and would otherwise glance without the placement.
# Reads the promoted top-level fields first, falling back to the nested joint
# object so a fit saved before the promotion still glances. Both columns are
# written whenever the record exists, filling the absent one with NA, so a batch
# summary rbind()ing one row per species gets a rectangular frame instead of
# losing the column on whichever rows recentered.
.tobs_glance_outer_grid <- function(g, x) {
  og <- .tobs_promote_outer_grid(x) %||% .tobs_promote_outer_grid(.tobs_joint_fit(x))
  if (is.null(og)) return(g)
  g$outer_grid_placement <- og$outer_grid_placement %||% NA_character_
  g$outer_grid_recenter_declined <- og$outer_grid_recenter_declined %||% NA_character_
  g
}

# Resolve a grid amplitude axis to a per-draw vector. A multi-block fit prefixes
# the axis with its block (`b<k>.sigma`); a single-block fit uses the bare name.
# `cells` is the outer-grid cell each draw came from, so the returned length-n
# vector carries the amplitude active for that draw. A missing axis returns a
# constant `default`.
# A field block's per-draw SD, read from whichever axis it rides: `b<k>.sigma`
# when the block is copied onto another arm, `b<k>.tau` (SD = 1/sqrt(tau)) when
# it is not. Every consumer of a field amplitude goes through this, so the two
# parameterizations cannot drift apart.
.tobs_joint_field_sd <- function(theta_grid, cells, block) {
  cn <- colnames(theta_grid)
  sig_col <- sprintf("b%d.sigma", block)
  tau_col <- sprintf("b%d.tau", block)
  if (sig_col %in% cn) return(as.numeric(theta_grid[cells, sig_col]))
  if (tau_col %in% cn) return(1.0 / sqrt(as.numeric(theta_grid[cells, tau_col])))
  .tobs_joint_amp(theta_grid, cells, block, "sigma")
}

.tobs_joint_amp <- function(theta_grid, cells, block, name, default = 1) {
  cn <- colnames(theta_grid)
  j  <- match(paste0("b", block, ".", name), cn)
  if (is.na(j)) j <- match(name, cn)
  if (is.na(j)) return(rep(default, length(cells)))
  as.numeric(theta_grid[cells, j])
}

# Draw the grid-integrated joint posterior and normalize it into a family-
# agnostic bundle. Returns:
#   $n         number of draws
#   $positive  "lognormal" / "beta"
#   $cells     length-n outer-grid cell each draw came from
#   $disp      length-n per-draw positive-arm dispersion (residual SD for
#              lognormal, precision for beta), read off the `phi_pos` grid axis
#   $b         list(occ = [n x p_occ], det = [n x p_det] | NULL,
#                   pos = [n x p_pos]) of per-arm coefficient draws (in the
#              fitted, scaled design space)
#   $blocks    list of shared-field blocks, each:
#                $z       [n x n_cells] unit-variance field draws
#                $amp_occ [n] occupancy-arm field amplitude per draw
#                $amp_pos [n] positive-arm field amplitude per draw
#                $weight  NULL (intercept field) or the per-cell covariate name
#                         weighting this (trend / SVC) field
#   $n_cells   spatial-unit count
.tobs_joint_draws <- function(object, n = 1000L) {
  jf <- .tobs_joint_fit(object)
  if (is.null(jf)) {
    stop("This fit carries no joint nested-Laplace object to sample ",
         "(`$joint_fit` / `$joint`).", call. = FALSE)
  }
  layout <- jf$arm_layout
  n_arms <- layout$n_arms %||% length(layout$p)
  if (isTRUE(object$occu_only_joint)) {
    .tobs_joint_draws_occu(object, jf, layout, n)
  } else if (n_arms == 3L) {
    .tobs_joint_draws_occu_cover(object, jf, layout, n)
  } else if (isTRUE(object$armspecific)) {
    .tobs_joint_draws_cover_armspecific(object, jf, layout, n)
  } else {
    .tobs_joint_draws_cover(object, jf, layout, n)
  }
}

# occu single-arm (2-arm psi/p): one or more independent ICAR fields on the
# OCCUPANCY (psi) arm only -- no cover arm, no copy. Each block is a unit-variance
# latent z scaled on the occupancy arm by its own amplitude (b<b>.sigma); the
# detection (p) arm carries no field (amp_pos = 0). Block 1 is the unweighted
# intercept field; blocks 2.. carry the per-cell trend weight. The bundle reuses
# the cover roster slot names ("occ" for psi, "pos" empty) so the shared
# `.tobs_joint_arm_eta` accumulator and the predict path read it unchanged: the
# occupancy psi arm is "occ", and there is no positive arm.
.tobs_joint_draws_occu <- function(object, jf, layout, n) {
  tg      <- jf$theta_grid
  n_cells <- object$model$n_cells %||% object$model$n_sites
  p       <- layout$p
  bstart  <- layout$beta_start

  idx_occ <- bstart[1L] + seq_len(p[1L])
  idx_det <- bstart[2L] + seq_len(p[2L])
  starts  <- layout$field_starts %||% layout$phi_start
  n_field <- length(starts)
  field_idx <- lapply(starts, function(s0) s0 + seq_len(n_cells))

  idx   <- c(idx_occ, idx_det, unlist(field_idx))
  D     <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = n)
  cells <- attr(D, "cells")

  off  <- 0L
  take <- function(k) { v <- D[, off + seq_len(k), drop = FALSE]; off <<- off + k; v }
  b_occ <- take(p[1L]); b_det <- take(p[2L])

  cn <- colnames(tg)
  trend_cols <- object$trend_weights %||% object$trend_weight
  blocks <- lapply(seq_len(n_field), function(b) {
    z <- take(n_cells)
    # Each ICAR block grids on precision tau (axis b<b>.tau); the occupancy-arm
    # field amplitude is the SD sigma = 1 / sqrt(tau) per draw cell. The
    # detection arm carries no field (amp_pos = 0, the 0-sentinel node index at
    # fit time already excluded it).
    tau_col <- sprintf("b%d.tau", b)
    amp <- if (tau_col %in% cn) 1.0 / sqrt(as.numeric(tg[cells, tau_col]))
           else rep(1.0, length(cells))
    list(z = z, amp_occ = amp, amp_pos = rep(0, length(cells)),
         weight = if (b == 1L) NULL else trend_cols[[b - 1L]])
  })

  list(n = n, positive = NA_character_, cells = cells,
       disp = rep(1, length(cells)),
       b = list(occ = b_occ, det = b_det, pos = NULL),
       blocks = blocks, n_cells = n_cells)
}

# cover (2-arm occ/pos) with arm-specific separate latents: one or more NON-copied
# areal blocks, each placed on exactly ONE arm. Block b is stored as a
# unit-precision latent z; its amplitude on the active arm is sigma_b (=
# 1/sqrt(tau_b) for icar/car_proper, the bym2 mixed amplitude otherwise), and its
# amplitude on the OTHER arm is 0 (no cross-arm copy). The bundle's per-block
# amp_occ / amp_pos encode this directly -- one is the field amplitude, the other
# is 0 -- so the shared `.tobs_joint_arm_eta` accumulator scatters each block onto
# its own arm only. A non-intercept (slope) field also carries its per-arm weight.
.tobs_joint_draws_cover_armspecific <- function(object, jf, layout, n) {
  tg       <- jf$theta_grid
  cn       <- colnames(tg)
  positive <- object$positive %||% "lognormal"
  p        <- layout$p
  bstart   <- layout$beta_start
  meta     <- object$armspec_blocks

  idx_occ <- bstart[1L] + seq_len(p[1L])
  idx_pos <- bstart[2L] + seq_len(p[2L])

  field_starts <- layout$field_starts
  n_field <- length(field_starts %||% integer(0))
  if (n_field != length(meta)) {
    stop("internal: arm-specific cover fit has ", length(meta), " field block(s) ",
         "but the joint layout reports ", n_field, ".", call. = FALSE)
  }
  # Each block records its own n_nodes (the graph it sits on). The latent SPAN
  # is n_nodes for an intrinsic (icar / car_proper) block but 2 * n_nodes for a
  # BYM2 block, which stores the structured phi followed by the iid theta.
  n_nodes <- vapply(meta, function(m) as.integer(m$n_nodes), integer(1))
  is_bym2 <- vapply(meta, function(m) identical(m$type, "bym2"), logical(1))
  field_span <- ifelse(is_bym2, 2L * n_nodes, n_nodes)
  field_idx <- lapply(seq_len(n_field), function(b)
    field_starts[b] + seq_len(field_span[b]))

  idx <- c(idx_occ, idx_pos, unlist(field_idx))
  D   <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = n)
  cells <- attr(D, "cells")
  off <- 0L
  take <- function(k) { v <- D[, off + seq_len(k), drop = FALSE]; off <<- off + k; v }
  b_occ <- take(p[1L]); b_pos <- take(p[2L])

  blocks <- lapply(seq_len(n_field), function(b) {
    nn <- n_nodes[b]
    # Field amplitude on the outer grid: b<b>.sigma (bym2) or 1/sqrt(b<b>.tau)
    # (icar / car_proper), per draw cell. The active arm scales z by this; the
    # inactive arm by 0 (no copy).
    sig_col <- sprintf("b%d.sigma", b)
    tau_col <- sprintf("b%d.tau", b)
    if (sig_col %in% cn) {
      amp <- as.numeric(tg[cells, sig_col])
    } else if (tau_col %in% cn) {
      amp <- 1.0 / sqrt(as.numeric(tg[cells, tau_col]))
    } else {
      amp <- rep(1.0, length(cells))
    }
    if (is_bym2[b]) {
      # Reconstruct the rho-mixed unit field from the two sub-blocks, the way
      # the shared-field path does: z = sqrt(rho) * sf * phi + sqrt(1-rho) * theta.
      raw   <- take(2L * nn)
      phi   <- raw[, seq_len(nn), drop = FALSE]
      theta <- raw[, nn + seq_len(nn), drop = FALSE]
      rho   <- as.numeric(tg[cells, sprintf("b%d.rho", b)])
      sf    <- meta[[b]]$scale_factor %||% 1.0
      z <- sweep(phi, 1L, sqrt(pmax(rho, 0) + 1e-10) * sf, "*") +
           sweep(theta, 1L, sqrt(pmax(1 - rho, 0) + 1e-10), "*")
    } else {
      z <- take(nn)
    }
    slot <- meta[[b]]$slot
    amp_occ <- if (slot == 1L) amp else rep(0, length(cells))
    amp_pos <- if (slot == 2L) amp else rep(0, length(cells))
    # The block sits on ONE arm, encoded by the zero amplitude on the other arm.
    # Its node map and per-cell weight are NOT stored here: the consumer supplies
    # them via `.tobs_joint_arm_eta`'s `units` / `wfun` -- predict() the newdata
    # cell map and column, the pointwise-loglik consumer the per-observation map
    # and weight built from the fit's armspec_blocks. A non-intercept (slope)
    # field carries its covariate column name to look up.
    wt <- if (isTRUE(meta[[b]]$is_intercept)) NULL else meta[[b]]$column_name
    list(z = z, amp_occ = amp_occ, amp_pos = amp_pos, weight = wt)
  })

  list(n = n, positive = positive, cells = cells,
       disp = .tobs_joint_amp(tg, cells, 1L, "phi_pos", default = 1),
       b = list(occ = b_occ, det = NULL, pos = b_pos),
       blocks = blocks, n_cells = n_nodes[1L])
}

# occu_cover (3-arm psi/p/pos): one or more ICAR fields stored as unit-variance
# z. The occupancy arm scales each block by `sigma`, the positive arm by
# `alpha * sigma`, the detection arm not at all (the field's `field_coef = 0` on
# p). A trend fit adds weighted blocks, each with its own per-cell weight column.
.tobs_joint_draws_occu_cover <- function(object, jf, layout, n) {
  tg       <- jf$theta_grid
  positive <- object$positive %||% "lognormal"
  # Field nodes (cells), not occupancy units (sites): under group_var the shared
  # field is sized by the graph, so n_sites can exceed it. Each block holds
  # n_cells field entries; consumers map sites -> cells via model$site_cell.
  n_cells  <- object$model$n_cells %||% object$model$n_sites
  p        <- layout$p
  bstart   <- layout$beta_start

  idx_occ <- bstart[1L] + seq_len(p[1L])
  idx_det <- bstart[2L] + seq_len(p[2L])
  idx_pos <- bstart[3L] + seq_len(p[3L])
  starts  <- layout$field_starts %||% layout$phi_start
  n_field <- length(starts)
  field_idx <- lapply(starts, function(s0) s0 + seq_len(n_cells))

  # Per-arm RE latent draws: each RE block stored its latent column indices in
  # fit$re[[arm]]$latent_idx; draw them from the SAME grid-integrated posterior
  # so the BLUP offset is marginalized over the joint, not plugged in at the
  # mode. They trail the fields in the latent vector.
  re_meta <- Filter(function(r) !is.null(r$latent_idx), object$re %||% list())
  re_idx  <- lapply(re_meta, function(r) as.integer(r$latent_idx))

  idx   <- c(idx_occ, idx_det, idx_pos, unlist(field_idx), unlist(re_idx))
  D     <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = n)
  cells <- attr(D, "cells")

  off  <- 0L
  take <- function(k) {
    v <- D[, off + seq_len(k), drop = FALSE]; off <<- off + k; v
  }
  b_occ <- take(p[1L]); b_det <- take(p[2L]); b_pos <- take(p[3L])

  trend_cols  <- object$trend_weights %||% object$trend_weight
  # Per-block arm + weight labels. A shared occupancy field scales occ by
  # b<k>.sigma and cover by b<k>.alpha * b<k>.sigma; an arm-specific cover field
  # (arm == "pos") scales cover by its OWN b<k>.sigma and occ by 0 (no copy).
  # `field_specs` labels every block; older fits (no field_specs) fall back to the
  # all-shared convention.
  field_specs <- object$field_specs
  cn <- colnames(tg)
  blocks <- lapply(seq_len(n_field), function(b) {
    z     <- take(n_cells)
    spec  <- if (!is.null(field_specs) && b <= length(field_specs))
               field_specs[[b]] else NULL
    if (!is.null(spec) && identical(spec$arm, "pos")) {
      # Non-copied ICAR: amplitude is its own SD, from b<k>.sigma or 1/sqrt(b<k>.tau).
      amp <- .tobs_joint_field_sd(tg, cells, b)
      list(z = z, amp_occ = rep(0, length(cells)), amp_pos = amp,
           weight = spec$weight)
    } else {
      sigma <- .tobs_joint_field_sd(tg, cells, b)
      # No alpha axis = the block is not copied onto the pos arm, so its
      # amplitude there is 0. Defaulting to 1 would turn a decoupled fit into a
      # full-amplitude copy in the draws.
      alpha <- .tobs_joint_amp(tg, cells, b, "alpha", default = 0)
      wt <- if (!is.null(spec)) spec$weight
            else if (b == 1L) NULL else trend_cols[[b - 1L]]
      list(z = z, amp_occ = sigma, amp_pos = alpha * sigma, weight = wt)
    }
  })

  re_draws <- NULL
  if (length(re_meta) > 0L) {
    re_draws <- lapply(re_meta, function(r) {
      # [n_draws x (n_coefs * n_groups)] coefficient-major latent draws; predict
      # reshapes to per-coefficient per-group and weights by the slope covariate
      # (a random intercept is the n_coefs = 1 case, [n_draws x n_groups]).
      nc <- r$n_coefs %||% 1L; ng <- r$n_groups
      d  <- take(length(r$latent_idx))
      # Back-transform a slope coefficient's draws from the standardized covariate
      # the fit ran on to its natural units (divide coefficient block c by its
      # scale), so predict() weights by the raw covariate column in `newdata`.
      sc <- r$coef_scales %||% rep(1, nc)
      if (nc > 1L && any(sc != 1)) {
        for (cc in seq_len(nc)) {
          cols <- (cc - 1L) * ng + seq_len(ng)
          if (sc[cc] != 1) d[, cols] <- d[, cols] / sc[cc]
        }
      }
      list(arm = r$arm, var = r$var, levels = r$levels,
           n_coefs = nc, n_groups = ng,
           coef_names = r$coef_names,
           covariate_names = r$covariate_names,
           has_intercept = r$has_intercept %||% TRUE,
           draws = d)
    })
    names(re_draws) <- names(re_meta)
  }

  # Pos-arm dispersion: the `phi_pos` axis when it is integrated on the outer
  # grid (control$phi.grid.pos, or the latent path's sigma_u); otherwise the
  # dispersion the fit held FIXED in the cell-coupling spec. Falling back to a
  # bare 1 would score every spatial occu_cover fit at unit dispersion regardless
  # of the value the spec used.
  fixed_disp <- object$model$cover_pos_disp %||% 1
  list(n = n, positive = positive, cells = cells,
       disp = .tobs_joint_amp(tg, cells, 1L, "phi_pos", default = fixed_disp),
       b = list(occ = b_occ, det = b_det, pos = b_pos),
       blocks = blocks, n_cells = n_cells, re = re_draws)
}

# cover (2-arm occ/pos): a single shared field. Under the (sigma_occ, sigma_pos)
# reparam the occupancy arm scales the unit-variance field z by `sigma_occ`, the
# positive arm by `sigma_pos`. The field is stored as phi (ICAR / proper CAR) or
# phi + theta (BYM2); for BYM2 the unit-variance z is
#   z = sqrt(rho) * scale_factor * phi + sqrt(1 - rho) * theta
# reconstructed per draw from the draw's grid rho.
.tobs_joint_draws_cover <- function(object, jf, layout, n) {
  tg       <- jf$theta_grid
  cn       <- colnames(tg)
  positive <- object$positive %||% "lognormal"
  p        <- layout$p
  bstart   <- layout$beta_start

  idx_occ <- bstart[1L] + seq_len(p[1L])
  idx_pos <- bstart[2L] + seq_len(p[2L])

  # Coupled multi-block path (intercept field + one or more SVC trend fields):
  # ICAR blocks under the (sigma, alpha) per-block copy convention, axes named
  # b<k>.sigma / b<k>.alpha. The occupancy arm scales block k by b<k>.sigma, the
  # positive arm by b<k>.alpha * b<k>.sigma; block 1 is the unweighted intercept
  # field, blocks 2.. carry the per-observation trend weight.
  field_starts <- layout$field_starts
  n_field <- length(field_starts %||% integer(0))
  if (n_field > 1L) {
    n_cells <- object$n_cells %||% as.integer(field_starts[2L] - field_starts[1L])
    field_idx <- lapply(field_starts, function(s0) s0 + seq_len(n_cells))
    idx <- c(idx_occ, idx_pos, unlist(field_idx))
    D   <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = n)
    cells <- attr(D, "cells")
    off <- 0L
    take <- function(k) { v <- D[, off + seq_len(k), drop = FALSE]; off <<- off + k; v }
    b_occ <- take(p[1L]); b_pos <- take(p[2L])
    trend_cols <- object$trend_weights %||% list(object$trend_weight)
    blocks <- lapply(seq_len(n_field), function(b) {
      z     <- take(n_cells)
      sigma <- .tobs_joint_field_sd(tg, cells, b)
      # No alpha axis = the block is not copied onto the pos arm, so its
      # amplitude there is 0. Defaulting to 1 would turn a decoupled fit into a
      # full-amplitude copy in the draws.
      alpha <- .tobs_joint_amp(tg, cells, b, "alpha", default = 0)
      list(z = z, amp_occ = sigma, amp_pos = alpha * sigma,
           weight = if (b == 1L) NULL else trend_cols[[b - 1L]])
    })
    return(list(n = n, positive = positive, cells = cells,
                disp = .tobs_joint_amp(tg, cells, 1L, "phi_pos", default = 1),
                b = list(occ = b_occ, det = NULL, pos = b_pos),
                blocks = blocks, n_cells = n_cells))
  }

  phi_start   <- layout$phi_start
  theta_start <- layout$theta_start
  n_phi <- if (!is.null(theta_start)) as.integer(theta_start - phi_start)
           else as.integer(layout$n_x - phi_start)
  phi_idx   <- phi_start + seq_len(n_phi)
  has_theta <- !is.null(theta_start)
  theta_idx <- if (has_theta) theta_start + seq_len(as.integer(layout$n_x - theta_start))
               else integer(0)

  idx   <- c(idx_occ, idx_pos, phi_idx, theta_idx)
  D     <- tulpa::tulpa_posterior_draws(jf, idx = idx, n = n)
  cells <- attr(D, "cells")

  off  <- 0L
  take <- function(k) {
    v <- D[, off + seq_len(k), drop = FALSE]; off <<- off + k; v
  }
  b_occ <- take(p[1L]); b_pos <- take(p[2L])
  phi   <- take(n_phi)

  # Unit-variance field z. BYM2 mixes the structured (phi) and unstructured
  # (theta) components by the draw's grid rho; ICAR / CAR use phi directly.
  has_rho <- "rho" %in% cn
  if (has_theta && has_rho) {
    theta <- take(length(theta_idx))
    rho   <- as.numeric(tg[cells, "rho"])
    sf    <- as.numeric(attr(jf, "scale_factor") %||% 1.0)
    z <- sweep(phi, 1L, sqrt(pmax(rho, 0) + 1e-10) * sf, "*") +
         sweep(theta, 1L, sqrt(pmax(1 - rho, 0) + 1e-10), "*")
  } else {
    z <- phi
  }

  amp_occ <- .tobs_joint_amp(tg, cells, 1L, "sigma_occ")
  amp_pos <- .tobs_joint_amp(tg, cells, 1L, "sigma_pos")
  block <- list(z = z, amp_occ = amp_occ, amp_pos = amp_pos, weight = NULL)

  list(n = n, positive = positive, cells = cells,
       disp = .tobs_joint_amp(tg, cells, 1L, "phi_pos", default = 1),
       b = list(occ = b_occ, det = NULL, pos = b_pos),
       blocks = list(block), n_cells = n_phi)
}

# Pointwise-log-likelihood consumer helpers for arm-specific cover() fits: the
# per-arm fields store no node map / weight on their bundle blocks, so the loglik
# consumer -- which runs over the fit's observations -- rebuilds them from
# `object$armspec_blocks` and hands them to `.tobs_joint_arm_eta` as `units` /
# `wfun`. The predict consumer instead supplies the newdata cell map and a
# newdata-column lookup. Every block on an arm shares one per-observation node
# map (they index the same graph), so the first block on the arm carries it.
#
# Per-arm per-observation node map; `seq_len(n_rows)` when the arm has no field
# (the amplitude check in `.tobs_joint_arm_eta` then never indexes it).
.tobs_armspec_obs_units <- function(object, slot, n_rows) {
  for (m in object$armspec_blocks) {
    if (isTRUE(m$slot == slot)) return(as.integer(m$idx_active))
  }
  seq_len(n_rows)
}

# Per-arm weight lookup (column name -> per-observation weight) for the trend /
# SVC blocks on this arm; NULL when the arm carries no weighted block.
.tobs_armspec_obs_wfun <- function(object, slot) {
  lut <- list()
  for (m in object$armspec_blocks) {
    if (!isTRUE(m$slot == slot) || isTRUE(m$is_intercept)) next
    w <- if (slot == 1L) m$weight_occ else m$weight_pos
    if (!is.null(w)) lut[[m$column_name]] <- as.numeric(w)
  }
  if (length(lut) == 0L) return(NULL)
  function(nm) {
    v <- lut[[nm]]
    if (is.null(v)) {
      stop("internal: arm-specific weight column '", nm,
           "' is absent from the fit's armspec_blocks.", call. = FALSE)
    }
    v
  }
}

# Per-arm linear predictor for every draw: [nrow(X) x n]. `X` is the arm's
# design (scaled, when the family autoscales), `arm` is "occ" / "det" / "pos",
# `units` maps each design row to its spatial unit (1..n_cells), and `wfun`
# resolves a weighted field block's per-cell covariate to a length-nrow(X)
# vector (only consulted for weighted blocks). The detection arm sees no field.
# This is the single source of truth for the field accumulation across both
# families and both consumers (prediction and pointwise log-likelihood).
#
# A block's arm membership is carried entirely by its per-arm amplitude
# (`amp_occ` / `amp_pos`): a shared field scales both arms, an arm-specific field
# has zero amplitude on the arm it does not sit on (no cross-arm copy). Every
# block on an arm then reads the SAME consumer-supplied `units` / `wfun` --
# predict() passes the newdata cell map and a newdata-column lookup; the
# pointwise-log-likelihood consumer passes the per-observation node map and
# weight. Blocks store neither map nor weight themselves, so an arm-specific
# field maps correctly at predict time.
.tobs_joint_arm_eta <- function(bundle, X, arm, units, wfun = NULL) {
  B <- bundle$b[[arm]]
  if (is.null(B)) {
    stop("Arm '", arm, "' is absent from this joint fit.", call. = FALSE)
  }
  eta <- tcrossprod(X, B)                       # [nrow x n]
  if (identical(arm, "det")) return(eta)        # detection arm carries no field
  is_occ <- identical(arm, "occ")
  for (blk in bundle$blocks) {
    amp <- if (is_occ) blk$amp_occ else blk$amp_pos
    # Zero amplitude on this arm = block not on this arm; skip before any node
    # indexing, so a block sized to its own graph is never over-indexed.
    if (all(amp == 0)) next
    z_unit <- blk$z[, units, drop = FALSE]      # [n x nrow]; `units` maps each
    contr  <- t(z_unit * amp)                   # design row to its field cell
    if (!is.null(blk$weight)) {
      # Weighted (trend / SVC) field: scale each row by its per-cell covariate,
      # resolved by the consumer (predict: from newdata; loglik: per-observation).
      if (is.null(wfun)) {
        stop("Weighted field block '", blk$weight,
             "' needs a per-cell weight lookup.", call. = FALSE)
      }
      contr <- contr * wfun(blk$weight)         # row c scaled by weight[c]
    }
    eta <- eta + contr
  }
  eta
}
