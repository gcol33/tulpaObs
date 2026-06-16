# =============================================================================
# formula_parse.R — Split a process formula into fixed effects + structured
#                   terms.
#
# A tobs() process formula may carry structured terms (spatial / temporal /
# random-effect / latent / copy). `.tobs_parse_formula()` separates the two:
#
#   ~ elev + forest + icar(graph = adj) + re(observer)
#     |  fixed effects  |     structured terms (registry hits)     |
#
# Detection is done on the *parse tree*: `terms()` splits the top-level terms,
# and each term label is re-parsed with `str2lang()` so the call head can be
# matched against the registry by symbol — no regex on deparsed strings. The
# fixed-effects labels are reassembled into a formula for `model.matrix()`; the
# structured-term calls are evaluated with the model `data` columns in scope so
# bare symbols (`gp(lon, lat)`, `re(observer)`) resolve to columns.
#
# lme4-style bar syntax is supported as sugar over re(): `(1 | g)`, `(x | g)`,
# `(x || g)` are rewritten into the equivalent re() calls on the formula AST
# *before* terms() runs (see `.tobs_desugar_bars`), so there is one parser and
# one term type — the bar is just a second spelling of re().
# =============================================================================

.tobs_spec_classes <- c("tobs_spatial", "tobs_re", "tobs_temporal",
                        "tobs_svc", "tobs_latent", "tobs_copy")

# Is a parsed term label a registered structured term? Returns the call's head
# name if so, otherwise NA. Inspects the AST, never the string.
.tobs_match_term <- function(label, reg_names) {
  e <- tryCatch(str2lang(label), error = function(...) NULL)
  if (is.call(e) && is.symbol(e[[1L]])) {
    head <- as.character(e[[1L]])
    if (head %in% reg_names) return(head)
  }
  NA_character_
}

# ---------------------------------------------------------------------------
# lme4-style bar syntax is sugar for re(). Each bar is rewritten into the
# equivalent re() call(s) on the formula AST before terms() runs, so the
# registry parser handles them through the one re() code path. A random-slope
# block carries the group intercept (n_coefs = intercept + slopes) with a
# Cholesky-correlated covariance; `||` makes it diagonal; `0 +` drops the
# intercept; several slopes stack into one block:
#
#   (1 | g)            -> re(g, type = "intercept")
#   (x | g)            -> re(g, type = "slope", covariate = cbind(x))
#   (1 + x + z | g)    -> re(g, type = "slope", covariate = cbind(x, z))
#   (x || g)           -> re(g, type = "slope", covariate = cbind(x), correlated = FALSE)
#   (0 + x | g)        -> re(g, type = "slope", covariate = cbind(x), intercept = FALSE)
#
# Grouping may be crossed or nested (lme4 semantics): a bar expands to one
# re() per implied grouping factor, joined with `+`.
#
#   (1 | g:h)          -> re(interaction(g, h, drop = TRUE), type = "intercept")
#   (1 | g/h)          -> re(g) + re(interaction(g, h, drop = TRUE))   [both intercept]
#
# The LHS (intercept + slopes) is distributed across every grouping factor the
# RHS implies. See gcol33/tulpaObs#10.

# Recognise a parenthesised bar term `( <lhs> | <rhs> )` / `( <lhs> || <rhs> )`.
# Returns list(op, lhs, rhs) or NULL. Inspects the AST, never the string.
.tobs_bar_spec <- function(e) {
  if (!is.call(e) || !identical(e[[1L]], as.name("(")) || length(e) != 2L)
    return(NULL)
  inner <- e[[2L]]
  if (!is.call(inner) || length(inner) != 3L) return(NULL)
  op <- inner[[1L]]
  if (identical(op, as.name("|")) || identical(op, as.name("||"))) {
    return(list(op = as.character(op), lhs = inner[[2L]], rhs = inner[[3L]]))
  }
  NULL
}

# Flatten a left-associative binary operator tree `a OP b OP c` into the list
# of its atoms (a, b, c). A node that is not that operator is a single atom.
.tobs_flatten_op <- function(e, op) {
  if (is.call(e) && identical(e[[1L]], as.name(op)) && length(e) == 3L) {
    return(c(.tobs_flatten_op(e[[2L]], op), list(e[[3L]])))
  }
  list(e)
}

# Build a grouping-factor expression from grouping atoms: one atom is returned
# as-is; several fold into interaction(a, b, ..., drop = TRUE) so the effect
# groups over observed factor combinations only.
.tobs_interaction_call <- function(atoms) {
  if (length(atoms) == 1L) return(atoms[[1L]])
  as.call(c(list(as.name("interaction")), atoms, list(drop = TRUE)))
}

# Expand a bar's RHS grouping into the list of grouping-factor expressions it
# implies (lme4 semantics): crossed `a:b` is one interaction factor; nested
# `a/b/c` is the cumulative sequence a, a:b, a:b:c; a bare factor is itself.
.tobs_bar_group_terms <- function(rhs) {
  if (is.call(rhs) && identical(rhs[[1L]], as.name("/"))) {
    atoms <- .tobs_flatten_op(rhs, "/")
    return(lapply(seq_along(atoms),
                  function(i) .tobs_interaction_call(atoms[seq_len(i)])))
  }
  if (is.call(rhs) && identical(rhs[[1L]], as.name(":"))) {
    return(list(.tobs_interaction_call(.tobs_flatten_op(rhs, ":"))))
  }
  list(rhs)
}

# Convert a recognised bar into the equivalent re() call(s) (a language
# object). Multi-slope bars stack their covariates with cbind(); slope-only
# bars set intercept = FALSE; `||` sets correlated = FALSE; nested/crossed
# grouping expands to one re() per grouping factor, joined with `+`.
.tobs_bar_to_re <- function(bar) {
  groups <- .tobs_bar_group_terms(bar$rhs)

  # Read intercept / slopes off the bar LHS via the parse tree (no regex): drop
  # the LHS into a one-sided formula template and let terms() decompose it.
  lhs_f <- ~ .
  lhs_f[[2L]] <- bar$lhs
  tt      <- stats::terms(lhs_f)
  has_int <- attr(tt, "intercept") == 1L
  slopes  <- attr(tt, "term.labels")
  correlated <- identical(bar$op, "|")

  make_re <- function(group) {
    if (length(slopes) == 0L) {
      if (!has_int)
        stop("Empty random effect `(0 | g)`: nothing to estimate.", call. = FALSE)
      return(call("re", group, type = "intercept"))
    }
    # Stack all slope covariates (1 or more) into one block; cbind() preserves
    # the column names for ranef() labelling.
    cov_expr <- as.call(c(list(as.name("cbind")), lapply(slopes, str2lang)))
    re_call <- call("re", group, type = "slope", covariate = cov_expr)
    if (!has_int)    re_call[["intercept"]]  <- FALSE
    if (!correlated) re_call[["correlated"]] <- FALSE
    re_call
  }

  re_calls <- lapply(groups, make_re)
  Reduce(function(a, b) call("+", a, b), re_calls)
}

# Walk the additive structure of an expression, rewriting bars in place.
.tobs_rewrite_bars <- function(e) {
  if (is.call(e) && identical(e[[1L]], as.name("+")) && length(e) == 3L) {
    e[[2L]] <- .tobs_rewrite_bars(e[[2L]])
    e[[3L]] <- .tobs_rewrite_bars(e[[3L]])
    return(e)
  }
  bar <- .tobs_bar_spec(e)
  if (!is.null(bar)) return(.tobs_bar_to_re(bar))
  e
}

# Desugar every lme4 bar in a formula's RHS into re() calls, preserving the
# formula's class and environment.
.tobs_desugar_bars <- function(formula) {
  n <- length(formula)
  formula[[n]] <- .tobs_rewrite_bars(formula[[n]])
  formula
}

# Collect the grouping-factor expressions of every lme4 bar in a formula's RHS,
# deparsed to the strings they would index by. Walks the additive AST with the
# same bar recogniser the desugarer uses (never a regex on the deparsed formula),
# so `(1 + x | cell)` and `(x || g/h)` contribute the same grouping factors they
# expand to (`cell`; `g`, `g:h`). A bar's RHS may name several factors (crossed /
# nested), so each is returned separately. Used by the cover()/occu_cover() guard
# (gcol33/tulpaObs#62) to flag a bar grouping factor that collides with an areal
# term's graph-node `group_var`.
.tobs_collect_bar_groups <- function(formula) {
  out <- character(0)
  walk <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("+")) && length(e) == 3L) {
      walk(e[[2L]]); walk(e[[3L]]); return(invisible())
    }
    bar <- .tobs_bar_spec(e)
    if (!is.null(bar)) {
      for (g in .tobs_bar_group_terms(bar$rhs)) {
        out[[length(out) + 1L]] <<- paste(deparse(g), collapse = "")
      }
    }
  }
  walk(formula[[length(formula)]])
  out
}

# Soft guard for the cover() / occu_cover() formula papercut (gcol33/tulpaObs#62).
# A user who learns the engine's inline-MCAR bar idiom (`tulpa::spatial(graph,
# ~ 1 + x | cell)`) may carry the `| cell` spelling into a cover formula, where it
# is legitimately parsed as a random effect, not a spatial field. RE bars are
# supported and must not be rejected; but when a bar's grouping factor is ALSO the
# graph-node `group_var` of an areal term in the same formula (the strong-signal
# confusion case), emit an informative message() that the bar is being fitted as an
# IID random effect, pointing to the spatial() bar / two-term form for a spatial
# field. Suppressible (message, not warning/error) and a no-op when the bar's
# factor is unrelated to any spatial term (an unambiguous, intended RE).
#
# `formula` is the original (pre-desugar) process formula; `spatial_specs` is the
# list of parsed `tobs_spatial` specs from that formula (each carrying `group_var`
# / `label`). Both inputs are parsed objects -- the formula AST is walked, never
# the deparsed string.
.tobs_cover_bar_re_guard <- function(formula, spatial_specs) {
  if (is.null(formula) || length(spatial_specs) == 0L) return(invisible())
  bar_groups <- .tobs_collect_bar_groups(formula)
  if (length(bar_groups) == 0L) return(invisible())
  gvs <- Filter(Negate(is.null),
                lapply(spatial_specs, function(s) s$group_var))
  gvs <- unique(unlist(gvs, use.names = FALSE))
  hits <- intersect(bar_groups, gvs)
  if (length(hits) == 0L) return(invisible())
  field <- spatial_specs[[1L]]$label %||% spatial_specs[[1L]]$type %||% "icar"
  for (g in hits) {
    message(sprintf(paste0(
      "cover()/occu_cover(): the bar `| %s` is being fitted as an IID random ",
      "effect, not a spatial field, even though `%s` is also the graph-node ",
      "group_var of an areal term. For a spatial field on `%s`, use ",
      "spatial(~ ... || %s, graph = <adj>) or the two-term form ",
      "%s(graph = <adj>, group_var = \"%s\") + %s(graph = <adj>, weight = ..., ",
      "group_var = \"%s\"). Suppress with suppressMessages()."),
      g, g, g, g, field, g, field, g))
  }
  invisible()
}

#' Parse a process formula into fixed effects and structured terms
#'
#' @param formula a one- or two-sided formula for a single process (e.g. the
#'   occupancy or detection linear predictor).
#' @param data the model data frame; columns are made available when evaluating
#'   structured-term calls.
#' @param env environment in which to resolve non-data symbols (e.g. an
#'   adjacency matrix passed as `graph = adj`). Defaults to the formula's
#'   environment.
#' @return A list with `fe_formula` (the fixed-effects formula for
#'   `model.matrix`) and `terms` (a list of `tobs_*` specs).
#' @keywords internal
.tobs_parse_formula <- function(formula, data = NULL,
                                env = environment(formula)) {
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  }
  if (is.null(env)) env <- parent.frame()

  # lme4 bar syntax -> re() calls before terms() sees the formula.
  formula <- .tobs_desugar_bars(formula)

  tt        <- if (is.null(data)) stats::terms(formula, keep.order = TRUE)
               else stats::terms(formula, data = data, keep.order = TRUE)
  labels    <- attr(tt, "term.labels")
  intercept <- attr(tt, "intercept")
  reg_names <- .tobs_term_names()

  heads    <- vapply(labels, .tobs_match_term, character(1), reg_names = reg_names)
  is_spec  <- !is.na(heads)
  fe_labels <- labels[!is_spec]

  fe_formula <- stats::reformulate(
    termlabels = if (length(fe_labels)) fe_labels else "1",
    intercept  = as.logical(intercept)
  )
  # model.matrix() resolves any non-data symbols (offsets, transforms) in this
  # environment, so anchor the rebuilt formula to the original one's env.
  environment(fe_formula) <- env

  terms_list <- list()
  if (any(is_spec)) {
    reg_env  <- list2env(.tobs_terms, parent = env)
    data_env <- if (!is.null(data)) list2env(as.list(data), parent = reg_env)
                else reg_env
    for (label in labels[is_spec]) {
      call <- str2lang(label)
      spec <- tryCatch(
        eval(call, envir = data_env),
        error = function(e) stop(sprintf(
          "Could not evaluate formula term `%s`: %s",
          label, conditionMessage(e)), call. = FALSE)
      )
      if (!inherits(spec, .tobs_spec_classes)) {
        stop(sprintf("Formula term `%s` did not produce a tobs term spec.",
                     label), call. = FALSE)
      }
      spec$term_call <- label
      terms_list[[length(terms_list) + 1L]] <- spec
    }
  }

  list(fe_formula = fe_formula, terms = terms_list)
}

# Parse every process formula in a named list, returning the fixed-effects
# formulas (same names) and a flat list of structured terms, each tagged with
# the 1-based process index it appeared in and the process name. `copy` terms
# are kept in the flat list for later resolution against defining `id`s.
#
#   processes: named list like list(psi = ~ ..., p = ~ ...)
.tobs_parse_processes <- function(processes, data, env = NULL) {
  fe        <- vector("list", length(processes))
  names(fe) <- names(processes)
  terms     <- list()
  for (k in seq_along(processes)) {
    f <- processes[[k]]
    if (is.null(f)) next   # leave fe[[k]] as the pre-allocated NULL slot
    e <- if (is.null(env)) environment(f) else env
    parsed  <- .tobs_parse_formula(f, data = data, env = e)
    fe[[k]] <- parsed$fe_formula
    for (spec in parsed$terms) {
      terms[[length(terms) + 1L]] <- list(
        spec    = spec,
        process = k,
        proc_name = names(processes)[k]
      )
    }
  }
  list(fe = fe, terms = terms)
}

# Resolve the flat tagged-term list into per-class specs whose `processes`
# field lists every process index they enter. `copy(id)` merges its own
# process into the process set of the term with matching `id`.
.tobs_resolve_terms <- function(tagged) {
  if (length(tagged) == 0L) return(list())

  # Index defining (non-copy) terms by id, and gather their process sets.
  defs   <- list()           # parallel to entries that are real terms
  by_id  <- list()           # id -> index into defs
  copies <- list()           # deferred copy references

  for (t in tagged) {
    if (inherits(t$spec, "tobs_copy")) {
      copies[[length(copies) + 1L]] <- t
      next
    }
    idx <- length(defs) + 1L
    t$processes <- t$process
    defs[[idx]] <- t
    id <- t$spec$id
    if (!is.null(id)) by_id[[id]] <- idx
  }

  for (cp in copies) {
    ref <- cp$spec$ref
    target <- by_id[[ref]]
    if (is.null(target)) {
      stop(sprintf("copy(\"%s\"): no term with id = \"%s\" found in the model.",
                   ref, ref), call. = FALSE)
    }
    defs[[target]]$processes <- sort(unique(c(defs[[target]]$processes,
                                              cp$process)))
  }

  defs
}

# Bind a family's process formulas in one shot: parse each into a
# fixed-effects design formula + structured terms, then resolve copy()
# references across processes. `processes` is a named list of formulas in the
# family's natural process order (e.g. list(psi = , p = ) for single-season;
# NULL entries are allowed and keep their slot). Returns:
#   fe    — named list of fixed-effects formulas (same names/positions)
#   terms — resolved structured-term list, each tagged with `$processes`
#           (integer process indices it enters), `$proc_name`, `$spec`.
.tobs_bind_formulas <- function(processes, data, env = NULL) {
  parsed <- .tobs_parse_processes(processes, data = data, env = env)
  list(fe = parsed$fe, terms = .tobs_resolve_terms(parsed$terms))
}

# Build a visit-level (long-form) detection / cover design matrix from a
# per-visit data frame. The data frame is in unit-major order (row k belongs to
# unit (k - 1) %/% max_per_unit + 1, visit (k - 1) %% max_per_unit + 1) and must
# have exactly n_units * max_per_unit rows; cells past a unit's observed visits
# are NA-padded and zero-filled here. `arm` labels the term in the row-count
# error.
#
# The visit design is stacked onto a site-level arm that already carries an
# intercept, so `drop_intercept = TRUE` removes the visit `(Intercept)` column to
# keep the combined design full rank; every observation family uses this
# convention. Returns NULL when the formula or data is absent, or when nothing
# but the dropped intercept remains.
.tobs_build_visit_X <- function(visit_formula, visit_data,
                                n_units, max_per_unit, arm,
                                drop_intercept = TRUE) {
  if (is.null(visit_formula) || is.null(visit_data)) return(NULL)
  expected_rows <- n_units * max_per_unit
  if (nrow(visit_data) != expected_rows) {
    stop(sprintf("%s visit data must have %d rows (one row per unit-visit), got %d",
                 arm, expected_rows, nrow(visit_data)), call. = FALSE)
  }
  # Build with the intercept present so a factor term gets reference (k - 1)
  # contrast coding; the intercept column is then dropped below when stacking
  # onto the site arm. Building under a no-intercept formula would expand a
  # factor to full (k) dummy coding, collinear with the site-level intercept.
  build_formula <- if (drop_intercept) {
    stats::update(visit_formula, ~ . + 1)
  } else {
    visit_formula
  }
  mf <- stats::model.frame(build_formula, visit_data,
                           na.action = stats::na.pass)
  X <- stats::model.matrix(build_formula, mf)
  X[is.na(X)] <- 0
  if (drop_intercept) {
    int_col <- match("(Intercept)", colnames(X))
    if (!is.na(int_col)) X <- X[, -int_col, drop = FALSE]
  }
  if (ncol(X) == 0L) return(NULL)
  X
}
