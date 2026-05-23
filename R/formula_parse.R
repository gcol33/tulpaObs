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
.tobs_parse_processes <- function(processes, data, env) {
  fe        <- vector("list", length(processes))
  names(fe) <- names(processes)
  terms     <- list()
  for (k in seq_along(processes)) {
    f <- processes[[k]]
    if (is.null(f)) { fe[[k]] <- NULL; next }
    parsed  <- .tobs_parse_formula(f, data = data, env = env)
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
