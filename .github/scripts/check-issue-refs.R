#!/usr/bin/env Rscript

# Fail when an issue-tracker reference appears in source, documentation or
# tests.
#
# A code comment states what the code does or why the domain demands it. A
# tracker pointer is neither: it does not survive a renumbering, a reader of
# the source has no path to a cross-repo tracker at all, and references inside
# roxygen render into shipped help pages, so anyone installing the package
# reads internal issue numbers. Where a reference sat beside a fact worth
# keeping, the fact stays and the pointer goes; where the sentence stopped
# making sense without the number, the information was in the wrong place.
#
# NEWS.md is exempt: a changelog entry citing the issue it closes is what the
# reference is for.

PATTERN <- "gcol33/tulpa(Obs)?#[0-9]+"

roots <- c("R", "src", "man", "tests", "vignettes")
roots <- roots[dir.exists(roots)]
files <- unlist(lapply(roots, list.files, recursive = TRUE, full.names = TRUE))
files <- files[file.info(files)$isdir %in% FALSE]

hits <- list()
for (f in files) {
  lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
  idx <- grep(PATTERN, lines)
  for (i in idx) hits[[length(hits) + 1L]] <- sprintf("%s:%d: %s", f, i,
                                                      trimws(lines[i]))
}

if (length(hits)) {
  cat("Issue-tracker references found outside NEWS.md:\n\n")
  cat(paste0("  ", unlist(hits)), sep = "\n")
  cat(sprintf("\n%d reference(s). State the constraint the reference stood ",
              length(hits)),
      "for, or delete it; keep the citation in NEWS.md.\n", sep = "")
  quit(status = 1L)
}

cat(sprintf("No issue-tracker references in %s (%d files scanned).\n",
            paste(roots, collapse = ", "), length(files)))
