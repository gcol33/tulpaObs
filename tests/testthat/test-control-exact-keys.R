# =============================================================================
# test-control-exact-keys.R - control knobs are read by exact key (#322).
#
# `$` on a list partial-matches. A `control$k` read where `k` is a strict prefix
# of another accepted key (`n.threads` / `n.threads.outer`, `alpha.grid` /
# `alpha.grid.trend`) returns the LONGER knob's value when only that one is set,
# and NULL when both are -- silently, and differently depending on how many of
# the pair the caller wrote. The structural test walks every function in the
# namespace so a new read of that shape fails here rather than in a fit.
# =============================================================================

test_that("cover() reads n.threads by exact key, not as a prefix of n.threads.outer", {
  f <- tulpaObs:::.cover_joint_control
  outer_only <- f(list(n.threads.outer = 8L))
  expect_identical(outer_only$n_threads, 1L)
  expect_identical(outer_only$n_threads_outer, 8L)
  expect_identical(f(list(n.threads = 4L))$n_threads, 4L)
  expect_identical(f(list())$n_threads, 1L)
})

# Every `control$key` read (on `control` itself or on `<x>$control`) and every
# `control[["key"]]` key in a function body, by walking the parse tree. A
# function taking `...` that binds a symbol to `list(...)` (`dots <- list(...)`)
# holds the same splatted knobs under that symbol, so its `$` / `[[` reads are
# tracked the same way; `dollar` entries are `"<symbol>$<key>"`.
.cek_walk <- function(e, visit) {
  empty <- list(quote(expr = ))
  walk_list <- function(lst) {
    for (i in seq_along(lst)) {
      if (identical(unname(lst[i]), empty)) next
      walk(lst[[i]])
    }
  }
  walk <- function(e) {
    if (is.call(e)) {
      visit(e)
      walk_list(as.list(e)[-1L])
    } else if (is.pairlist(e)) {
      walk_list(as.list(e))
    }
  }
  walk(e)
}

.cek_dots_symbols <- function(fun) {
  if (!"..." %in% names(formals(fun))) return(character(0))
  found <- character(0)
  .cek_walk(body(fun), function(e) {
    head <- e[[1L]]
    if ((identical(head, as.name("<-")) || identical(head, as.name("="))) &&
        length(e) == 3L && is.symbol(e[[2L]]) && is.call(e[[3L]]) &&
        identical(e[[3L]][[1L]], as.name("list")) && length(e[[3L]]) == 2L &&
        identical(e[[3L]][[2L]], as.name("...")))
      found <<- c(found, as.character(e[[2L]]))
  })
  unique(found)
}

.cek_reads <- function(fun) {
  tracked <- c("control", .cek_dots_symbols(fun))
  tracked_symbol <- function(e) {
    if (is.symbol(e) && as.character(e) %in% tracked) return(as.character(e))
    if (is.call(e) && identical(e[[1L]], as.name("$")) &&
        identical(as.character(e[[3L]]), "control")) return("control")
    NULL
  }
  acc <- list(dollar = character(0), exact = character(0))
  .cek_walk(formals(fun), function(e) NULL)
  .cek_walk(body(fun), function(e) {
    head <- e[[1L]]
    if (identical(head, as.name("$")) && length(e) == 3L) {
      sym <- tracked_symbol(e[[2L]])
      if (!is.null(sym))
        acc$dollar <<- c(acc$dollar, paste0(sym, "$", as.character(e[[3L]])))
    }
    if (identical(head, as.name("[[")) && length(e) >= 3L &&
        !is.null(tracked_symbol(e[[2L]])) && is.character(e[[3L]]))
      acc$exact <<- c(acc$exact, e[[3L]])
  })
  acc
}

.cek_key <- function(dollar) sub("^[^$]*\\$", "", dollar)

test_that("the walker tracks a `dots <- list(...)` binding like control", {
  f <- function(control, ...) {
    dots <- list(...)
    list(control$a, dots$b, dots[["c"]], other$d)
  }
  r <- .cek_reads(f)
  expect_setequal(r$dollar, c("control$a", "dots$b"))
  expect_identical(r$exact, "c")
  g <- function(control) {
    dots <- list(x = 1)
    dots$b
  }
  expect_identical(.cek_reads(g)$dollar, character(0))
})

test_that("no control$ read is a strict prefix of another accepted control key", {
  skip_on_cran()
  ns <- asNamespace("tulpaObs")
  reads <- list()
  keys  <- character(0)
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get(nm, envir = ns)
    if (!is.function(f) || is.primitive(f)) next
    r <- .cek_reads(f)
    if (length(r$dollar)) reads[[nm]] <- unique(r$dollar)
    keys <- c(keys, .cek_key(r$dollar), r$exact)
  }
  for (fam in names(ns$.tobs_family_methods)) {
    ctor <- get0(fam, envir = ns, mode = "function")
    obj  <- if (is.null(ctor)) NULL else tryCatch(ctor(), error = function(e) NULL)
    keys <- c(keys, obj$control_keys)
  }
  keys <- unique(keys[nzchar(keys)])
  # The walker must actually see reads, or an empty result proves nothing.
  expect_gt(length(reads), 10L)
  expect_true("n.threads.outer" %in% keys)

  shadowed <- character(0)
  for (fn in names(reads)) for (rd in reads[[fn]]) {
    longer <- keys[startsWith(keys, paste0(.cek_key(rd), "."))]
    if (length(longer))
      shadowed <- c(shadowed, sprintf("%s: %s prefixes %s", fn, rd,
                                      paste(longer, collapse = ", ")))
  }
  expect_identical(shadowed, character(0))
})
