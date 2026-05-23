# Standalone validation of the lme4 bar-syntax desugarer before wiring it into
# formula_parse.R. Verifies the AST rewrite produces the expected re() calls and
# that the LHS analysis (intercept? slopes?) is read off the parse tree.

# ---- candidate implementation (ported verbatim into R/formula_parse.R) ------

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

.tobs_bar_to_re <- function(bar) {
  group <- bar$rhs
  if (is.call(group) && (identical(group[[1L]], as.name("/")) ||
                         identical(group[[1L]], as.name(":")))) {
    stop(sprintf(
      "Random-effect grouping `%s` (nested/interaction) is not supported via `|`; build the factor explicitly, e.g. re(interaction(a, b)).",
      deparse(group)), call. = FALSE)
  }
  lhs_f <- ~ .
  lhs_f[[2L]] <- bar$lhs
  tt      <- stats::terms(lhs_f)
  has_int <- attr(tt, "intercept") == 1L
  slopes  <- attr(tt, "term.labels")
  correlated <- identical(bar$op, "|")

  if (length(slopes) == 0L) {
    if (!has_int) stop("Empty random effect `(0 | g)`: nothing to estimate.",
                       call. = FALSE)
    return(call("re", group, type = "intercept"))
  }
  if (length(slopes) == 1L) {
    if (!has_int) stop(sprintf(
      "Slope-only random effects `(0 + %s | g)` are not supported; the random slope always carries its group intercept. Use `(%s | g)` or call re() directly.",
      slopes, slopes), call. = FALSE)
    re_call <- call("re", group, type = "slope", covariate = str2lang(slopes))
    if (!correlated) re_call[["correlated"]] <- FALSE
    return(re_call)
  }
  stop(sprintf(
    "Multiple random slopes in one bar `(%s | g)` are not supported via `|`; split into separate terms or call re() directly.",
    paste(slopes, collapse = " + ")), call. = FALSE)
}

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

.tobs_desugar_bars <- function(formula) {
  n <- length(formula)
  formula[[n]] <- .tobs_rewrite_bars(formula[[n]])
  formula
}

# ---- checks -----------------------------------------------------------------

chk <- function(label, got, want) {
  ok <- identical(deparse(got), deparse(want))
  cat(sprintf("[%s] %s\n     got:  %s\n     want: %s\n",
              if (ok) "OK" else "FAIL", label,
              paste(deparse(got), collapse = ""),
              paste(deparse(want), collapse = "")))
  if (!ok) stop("mismatch")
}

chk("(1|g)",
    .tobs_desugar_bars(~ x + (1 | g)),
    ~ x + re(g, type = "intercept"))

chk("(x|g) correlated int+slope",
    .tobs_desugar_bars(~ x + (elev | g)),
    ~ x + re(g, type = "slope", covariate = elev))

chk("(1+x|g) == (x|g)",
    .tobs_desugar_bars(~ (1 + elev | g)),
    ~ re(g, type = "slope", covariate = elev))

chk("(x||g) uncorrelated",
    .tobs_desugar_bars(~ (elev || g)),
    ~ re(g, type = "slope", covariate = elev, correlated = FALSE))

chk("two bars + fixed effects",
    .tobs_desugar_bars(~ forest + (1 | site) + (elev || obs)),
    ~ forest + re(site, type = "intercept") +
      re(obs, type = "slope", covariate = elev, correlated = FALSE))

chk("two-sided detection formula",
    .tobs_desugar_bars(y ~ effort + (1 | observer)),
    y ~ effort + re(observer, type = "intercept"))

chk("no bars unchanged",
    .tobs_desugar_bars(~ elev + icar(graph = adj)),
    ~ elev + icar(graph = adj))

chk("transformed slope",
    .tobs_desugar_bars(~ (log(x) | g)),
    ~ re(g, type = "slope", covariate = log(x)))

# errors
err <- function(label, expr) {
  e <- tryCatch({ expr; NULL }, error = function(e) conditionMessage(e))
  cat(sprintf("[%s] %s -> %s\n", if (is.null(e)) "FAIL" else "OK", label,
              if (is.null(e)) "(no error)" else e))
  if (is.null(e)) stop("expected error")
}
err("slope-only (0+x|g)",   .tobs_desugar_bars(~ (0 + x | g)))
err("multi-slope (1+x+z|g)", .tobs_desugar_bars(~ (1 + x + z | g)))
err("nested (1|g/h)",        .tobs_desugar_bars(~ (1 | g / h)))
err("interaction (1|g:h)",   .tobs_desugar_bars(~ (1 | g : h)))

cat("\nALL PROBE CHECKS PASSED\n")
