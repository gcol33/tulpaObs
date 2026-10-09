# =============================================================================
# occu_cover_spatial.R - the areal-field pieces every spatial occu_cover()
# route shares: the ICAR precision Q = D - W and its Sorbye-Rue scale, unit-
# scale field draws for the simulators, and the resolver that reads the
# coupled field(s), the arm-specific fields and the occupancy-arm random
# intercept off the psi formula for the joint nested-Laplace fitter.
# =============================================================================


# ---------------------------------------------------------------------------
# Build the ICAR precision matrix Q = D - W from an adjacency.
# ---------------------------------------------------------------------------
.occu_cover_icar_Q <- function(adj) {
  n <- nrow(adj)
  D <- rowSums(adj)
  if (any(D == 0L)) {
    isolated <- which(D == 0L)
    stop(sprintf("ICAR graph has %d isolated node(s) (no neighbours): %s. ",
                 length(isolated),
                 paste(utils::head(isolated, 5L), collapse = ", ")),
         "Drop them or connect them before fitting.", call. = FALSE)
  }
  Q <- -as.matrix(adj)
  diag(Q) <- D
  Q
}


# ---------------------------------------------------------------------------
# Sorbye-Rue scaling factor: geometric mean of diag(Q^+) under the sum-to-zero
# constraint. Multiplying Q by this factor gives a precision whose
# generalized-inverse diagonal has geometric mean 1, so `sigma` in the linear
# predictor `sigma * z` is the geo-mean marginal SD of the field.
# Matches INLA's `scale.model = TRUE` (Sorbye & Rue 2014).
# ---------------------------------------------------------------------------
.occu_cover_icar_scale <- function(adj) {
  n <- nrow(adj)
  Q <- .occu_cover_icar_Q(adj)
  eig <- eigen(Q, symmetric = TRUE)
  pos <- eig$values > 1e-10
  if (!any(pos)) {
    stop("ICAR graph has no positive eigenvalues; check connectivity.",
         call. = FALSE)
  }
  V_p <- eig$vectors[, pos, drop = FALSE]
  Qinv_diag <- rowSums(V_p^2 /
                          matrix(eig$values[pos], n, sum(pos), byrow = TRUE))
  exp(mean(log(Qinv_diag)))
}


# ---------------------------------------------------------------------------
# Draw `n_draw` unit-scale ICAR fields on `adj`: f ~ N(0, Q^-) on the
# sum-to-zero constrained space, divided by sqrt(scale_q) so the geo-mean
# marginal SD is 1. `sigma * f` then carries geo-mean marginal SD `sigma`,
# which is the parameterisation the fitters and `simulate_occu_cover()` share.
# Returns an n_nodes x n_draw matrix.
# ---------------------------------------------------------------------------
.occu_cover_draw_icar_field <- function(adj, n_draw = 1L) {
  .tobs_draw_icar_unit(.occu_cover_icar_Q(adj), .occu_cover_icar_scale(adj),
                       n_draw)
}


# The grouping label and factor levels behind an occupancy-arm `re()` term. The
# resolved spec keeps the integer group codes only, so the label comes from the
# term call's own `group` argument (matched against the constructor's formals,
# not the deparsed text) and the levels from evaluating that argument in the
# term's own evaluation environment. `.tobs_index_codes()` codes a grouping
# vector through `as.factor()`, so those levels are in code order; a level count
# that disagrees with the recorded group count is left unlabelled rather than
# mismatched. Both feed the RE summary the fit reports (`fit$re$psi`, `ranef()`)
# on the same footing as the observation arms.
.occu_cover_psi_re_labels <- function(spec, formula, data) {
  out <- list(var = NULL, levels = NULL)
  cl <- tryCatch(str2lang(spec$term_call %||% ""), error = function(e) NULL)
  if (!is.call(cl)) return(out)
  mc <- tryCatch(match.call(.tobs_term_re, cl), error = function(e) NULL)
  if (is.null(mc) || is.null(mc$group)) return(out)
  out$var <- paste(deparse(mc$group), collapse = "")
  lev <- tryCatch(
    levels(as.factor(eval(mc$group,
                          .tobs_term_eval_env(data, environment(formula))))),
    error = function(e) NULL)
  if (length(lev) == as.integer(spec$n_groups)) out$levels <- lev
  out
}


# ---------------------------------------------------------------------------
# Pull the coupled spatial field(s) from the psi formula. Returns NULL when no
# spatial term is present, otherwise a list with `fe` (the fixed-effects psi
# formula) and `fields` (the ordered field specs: the one unweighted intercept
# field first, then any weighted spatially-varying-coefficient fields).
#
# A weighted areal term (`icar(graph = adj, weight = col)`) is a second coupled
# field -- a spatially-varying coefficient on `col` sharing the same areal
# graph -- the formula-DSL spelling of the trend field that `control$trend`
# also produces. The joint engine couples N such fields.
# ---------------------------------------------------------------------------
.occu_cover_spatial_fields <- function(formula, data, arm_fields = list()) {
  bind <- .tobs_bind_formulas(list(psi = formula), data)
  if (length(bind$terms) == 0L && length(arm_fields) == 0L) return(NULL)
  # `.tobs_bind_formulas` returns terms wrapped in `list(spec = ..., process = ...)`.
  spatial <- Filter(function(t) inherits(t$spec, "tobs_spatial"), bind$terms)
  # Placed arm fields (detection / positive): evaluate each lifted field call and
  # tag its spec with the arm, in the occurrence formula's environment (where the
  # occurrence field's own graph symbol resolves), then route alongside the
  # occurrence formula's own fields.
  if (length(arm_fields)) {
    extra <- lapply(arm_fields, function(af) {
      spec <- .tobs_eval_arm_field(af$call, af$arm, data, environment(formula))
      # Only the bar form carries an arm-specific field; a plain areal
      # constructor placed in an arm formula has no such spelling.
      if (!isTRUE(spec$is_bar)) {
        stop(sprintf(paste0(
          "occu_cover(): a field placed in the %s formula must use the bar form ",
          "spatial(~ 1 + w || cell, graph = adj); got %s()."),
          af$arm, spec$type %||% class(spec)[1L]), call. = FALSE)
      }
      list(spec = spec)
    })
    spatial <- c(spatial, extra)
  }
  if (length(spatial) == 0L) return(NULL)
  # A varying-coefficient bar (`spatial(~ 1 + w || node, graph = adj)`) desugars
  # in place to the intercept field + per-covariate trend fields, the same pair
  # the two-term form produces -- INDEPENDENT areal fields, each with its own
  # precision. A correlated bar (single `|`, the free-Sigma separable MCAR)
  # declares the SAME per-field design, but the intercept + coefficient fields
  # share a free cross-covariance Sigma on the occupancy arm and are copied onto
  # the cover arm together with one amplitude alpha. It is fitted as one coupled
  # MCAR block rather than independent fields, so it must be the only spatial term
  # (one field structure per fit); `correlated` flags it for the joint-coupled
  # fitter.
  specs <- list()
  correlated <- FALSE
  # Arm-specific INDEPENDENT fields, keyed by internal arm slot: "pos" (cover) and
  # "p" (detection). A single-arm `to =` bar (placement in the positive / detection
  # formula, lifted here) becomes a separate latent field on that arm ALONE, with
  # its own precision and NO cross-arm copy of the occupancy field's alpha. This is
  # the opt-in for arm-structured trends the alpha copy cannot express (the copy
  # collapses to a global slope when the shapes differ, e.g. flattening
  # delta_cover_cond). Distinct from the shared `||` desugar (which copies a
  # psi-anchored field onto the cover arm with alpha).
  armspec  <- list()
  arm_slot <- c(positive = "pos", detection = "p")
  for (t in spatial) {
    if (isTRUE(t$spec$is_bar)) {
      to <- t$spec$to %||% .tobs_cover_arms
      if (isTRUE(t$spec$correlated)) {
        if (correlated || length(specs) > 0L || length(spatial) > 1L) {
          stop(paste0(
            "occu_cover(): a correlated spatial bar (`|`) must be the only ",
            "spatial term in the psi formula (one MCAR field structure per ",
            "fit). Drop the other areal term(s), or use the INDEPENDENT ",
            "spelling `||` to combine separate per-coefficient fields."),
            call. = FALSE)
        }
        specs <- .tobs_expand_spatial_bar(t$spec, data)
        correlated <- TRUE
      } else if (length(to) == 1L) {
        if (identical(to, "presence")) {
          stop(paste0(
            "occu_cover(): there is no separate presence arm; the occupancy ",
            "field is the occurrence spatial() term. Write the field in the ",
            "`occurrence` formula, or place it in `positive` (cover) / ",
            "`detection`."), call. = FALSE)
        }
        if (!to %in% names(arm_slot)) {
          stop(sprintf(paste0(
            "occu_cover(): an arm-specific field belongs in the positive (cover) ",
            "or detection formula; got arm \"%s\"."), to), call. = FALSE)
        }
        slot <- arm_slot[[to]]
        if (!is.null(armspec[[slot]])) {
          stop(sprintf(paste0(
            "occu_cover(): at most one arm-specific field per arm (%s); ",
            "combine the coefficient fields into one bar in that arm's formula, ",
            "e.g. %s = ~ x + spatial(~ 1 + w || cell, graph = adj)."),
            to, to), call. = FALSE)
        }
        armspec[[slot]] <- .tobs_armspecific_bar_fields(t$spec, data)
      } else {
        specs <- c(specs, .cover_desugar_spatial_bar(t$spec, data))
      }
    } else {
      specs[[length(specs) + 1L]] <- t$spec
    }
  }

  # Soft guard: a bare `| / ||` RE bar whose grouping factor is also an areal
  # term's graph-node group_var is fitted as an IID random effect, not a spatial
  # field (the engine-bar-idiom papercut). Informs (message) without rejecting --
  # RE bars are legitimate -- and is silent when the bar's factor is unrelated to
  # any spatial term.
  .tobs_cover_bar_re_guard(formula, specs)

  bad <- Filter(function(s) !s$type %in% c("icar", "bym2"), specs)
  if (length(bad)) {
    stop(sprintf(
      "occu_cover() spatial path supports icar() or bym2() in the psi formula; got %s().",
      bad[[1L]]$type), call. = FALSE)
  }
  if (any(vapply(specs, function(s) identical(s$type, "bym2"), logical(1)))) {
    warning("occu_cover() reads bym2() as ICAR (rho fixed to 1).",
            call. = FALSE)
  }

  # Exactly one unweighted (intercept) field is the shared base field; the rest
  # are weighted SVC fields. The cell-coupling model has one node per cell, so
  # every coupled field shares the base graph (only the weight column differs).
  weighted <- vapply(specs, function(s) !is.null(s$weight), logical(1))
  base <- specs[!weighted]
  if (length(base) == 0L) {
    stop("occu_cover() spatial requires one unweighted intercept field ",
         "(e.g. icar(graph = adj)); only weighted SVC field(s) were given.",
         call. = FALSE)
  }
  if (length(base) > 1L) {
    stop("occu_cover() spatial supports exactly one unweighted intercept ",
         "field; got ", length(base), ". Additional coupled fields must be ",
         "weighted SVC terms, e.g. icar(graph = adj, weight = year).",
         call. = FALSE)
  }
  base_graph <- base[[1L]]$graph
  for (s in specs[weighted]) {
    if (!identical(dim(s$graph), dim(base_graph)) ||
        !all(s$graph == base_graph)) {
      stop("occu_cover() coupled fields must share the same areal graph as ",
           "the intercept field (same nodes / adjacency).", call. = FALSE)
    }
  }

  # Optional group_var maps each site (occupancy unit) to a field node, so the
  # site count can exceed the node count (e.g. site = cell-year sharing one cell
  # field). All coupled fields must name the same group_var (or none).
  gvs <- unique(Filter(Negate(is.null), lapply(specs, function(s) s$group_var)))
  if (length(gvs) > 1L) {
    stop("occu_cover() coupled fields must share a single group.var (or none).",
         call. = FALSE)
  }
  group_var <- if (length(gvs) == 1L) gvs[[1L]] else NULL

  # An optional per-group random INTERCEPT on the occupancy arm, layered on the
  # shared field (the consumer's field + per-group RE composition in the joint
  # cell-coupling engine). It joins the joint fit as a single `iid` prior block
  # whose per-cell offset rides the occupancy predictor; its variance integrates
  # on the outer grid alongside the field sigma / alpha. Scope: one
  # random-intercept term -- a scalar per group -- maps onto the one iid block. A
  # random slope or a correlated multi-coefficient block has no scalar-per-group
  # form here and errors (the non-spatial cover-hurdle EM is the route for richer
  # RE).
  re_terms <- Filter(function(t) inherits(t$spec, "tobs_re"), bind$terms)
  re_spec  <- NULL
  if (length(re_terms) > 0L) {
    if (length(re_terms) > 1L) {
      stop("occu_cover() spatial + RE supports a single random-intercept term ",
           "on the occupancy formula; got ", length(re_terms), ".", call. = FALSE)
    }
    rs <- re_terms[[1L]]$spec
    if (!identical(rs$type, "intercept")) {
      stop("occu_cover() spatial + RE supports a random INTERCEPT only ",
           "(e.g. (1 | group) / re(group)); a random slope or correlated block ",
           "is not wired on the shared-field joint engine.", call. = FALSE)
    }
    lab <- .occu_cover_psi_re_labels(rs, formula, data)
    re_spec <- list(group_idx = as.integer(rs$group_idx),
                    n_groups  = as.integer(rs$n_groups),
                    var       = lab$var, levels = lab$levels)
  }

  # Correlated (`|`) MCAR field requirements: at least one coefficient beyond
  # the intercept (a single field has no cross-covariance), intrinsic CAR only,
  # and no per-group RE (the MCAR block already spans the full coupled field
  # structure; a layered iid block is not wired with it yet).
  if (correlated) {
    if (length(specs) < 2L) {
      stop(paste0(
        "occu_cover(): a correlated spatial bar (`|`) needs at least one ",
        "coefficient beyond the intercept (e.g. spatial(~ 1 + x | cell, ",
        "graph = adj)); a single field has no cross-covariance to estimate. ",
        "Use icar()/`||` for an uncorrelated field."), call. = FALSE)
    }
    if (!all(vapply(specs, function(s) identical(s$type, "icar"), logical(1)))) {
      stop(paste0(
        "occu_cover(): a correlated spatial bar (`|`) uses the intrinsic CAR ",
        "(icar); bym2 is not supported for the MCAR field."), call. = FALSE)
    }
    if (!is.null(re_spec)) {
      stop(paste0(
        "occu_cover(): a correlated spatial bar (`|`) does not compose with a ",
        "per-group occupancy random effect; fit one or the other."),
        call. = FALSE)
    }
  }

  # Arm-specific fields (extended to the detection arm): each is a separate,
  # non-copied block on ONE arm (cover or detection); it composes with the shared
  # occupancy field (which still drives psi and, via the alpha copy,
  # delta_cover_exp) but NOT with the correlated `|` MCAR field (that already spans
  # the whole coupled structure with its own copy). Like the occu_cover shared
  # field, an arm-specific field is fitted as ICAR (rho fixed to 1); bym2 / car on
  # the bar is read as ICAR. Every arm-specific field shares the occupancy field's
  # areal graph (one node set).
  arm_label <- c(pos = "positive", p = "detection")
  for (slot in names(armspec)) {
    af <- armspec[[slot]]
    if (correlated) {
      stop(sprintf(paste0(
        "occu_cover(): an arm-specific field (to = \"%s\") does not compose with ",
        "a correlated `|` MCAR field; use one spatial structure."),
        arm_label[[slot]]), call. = FALSE)
    }
    if (!identical(af$type, "icar")) {
      warning(sprintf(paste0(
        "occu_cover() reads the arm-specific %s field as ICAR (rho fixed to 1); ",
        "bym2/car with free mixing on that arm is not yet wired."),
        arm_label[[slot]]), call. = FALSE)
      af$type <- "icar"
    }
    if (!identical(dim(af$graph), dim(base_graph)) ||
        !all(af$graph == base_graph)) {
      stop(sprintf(paste0(
        "occu_cover(): the arm-specific %s field (to = \"%s\") must share the ",
        "same areal graph as the occupancy field (same nodes / adjacency)."),
        arm_label[[slot]], arm_label[[slot]]), call. = FALSE)
    }
    armspec[[slot]] <- af
  }

  list(fe = bind$fe$psi, fields = c(base, specs[weighted]),
       armspec = armspec, pos_armspec = armspec[["pos"]],
       group_var = group_var, re = re_spec, correlated = correlated)
}
