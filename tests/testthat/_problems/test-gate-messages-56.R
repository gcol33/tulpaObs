# Extracted from test-gate-messages.R:56

# test -------------------------------------------------------------------------
skip_on_cran()
r_dir <- normalizePath(file.path(testthat::test_path(), "..", "..", "R"),
                         mustWork = FALSE)
skip_if_not(dir.exists(r_dir), "package sources not available")
count_directives <- function(fmt) {
    stripped <- gsub("%%", "", fmt, fixed = TRUE)
    m <- gregexpr("%[-+ #0]*[0-9*]*(\\.[0-9*]+)?[disefgGxXoObaA]", stripped)[[1]]
    if (identical(as.integer(m), -1L)) 0L else length(m)
  }
bad_directives <- function(fmt) {
    stripped <- gsub("%%", "", fmt, fixed = TRUE)
    m <- regmatches(stripped, gregexpr("%[-+ #0]*[0-9*]*(\\.[0-9*]+)?[a-zA-Z]",
                                       stripped))[[1]]
    m[!grepl("[disefgGxXoObaA]$", m)]
  }
bad <- character()
scan_call <- function(e, file) {
    if (!is.call(e)) return(invisible(NULL))
    fn <- e[[1]]
    if (is.name(fn) && identical(as.character(fn), "sprintf") && length(e) >= 2) {
      fmt <- e[[2]]
      if (is.character(fmt) && length(fmt) == 1) {
        args <- as.list(e)[-c(1, 2)]
        nms <- names(args)
        if (!is.null(nms)) args <- args[!nzchar(nms)]
        invalid <- bad_directives(fmt)
        if (count_directives(fmt) != length(args) || length(invalid)) {
          bad <<- c(bad, sprintf("%s: %s", basename(file), substr(fmt, 1, 60)))
        }
      }
    }
    for (i in seq_along(e)) {
      # A default-less formal parses to the empty symbol; forcing it errors.
      is_call_slot <- tryCatch(is.call(e[[i]]), error = function(...) FALSE)
      if (is_call_slot) scan_call(e[[i]], file)
    }
    invisible(NULL)
  }
for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
    for (e in parse(f, keep.source = FALSE)) scan_call(e, f)
  }
expect_equal(bad, character())
