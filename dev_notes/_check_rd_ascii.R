# Inspect the exact non-ASCII codepoints the broad grep flagged in tobs.Rd, and
# scan all of man/ for the DANGEROUS Unicode blocks (superscripts / subscripts /
# math operators / letterlike / arrows) that actually break the LaTeX manual.
L <- readLines("man/tobs.Rd", warn = FALSE)
for (i in c(53:56, 63)) {
  cps <- utf8ToInt(L[i]); nb <- cps[cps > 127L]
  cat(sprintf("L%d: %s\n", i, L[i]))
  if (length(nb))
    cat("    >127:", paste0("U+", toupper(format(as.hexmode(nb), width = 4)),
                            collapse = " "), "\n")
}

cat("\n--- dangerous-block scan across man/ ---\n")
danger <- "[⁰-₟℀-⅏←-⇿∀-⋿]"
mans <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
any_bad <- FALSE
for (f in mans) {
  ls <- readLines(f, warn = FALSE)
  h <- grep(danger, ls, perl = TRUE)
  if (length(h)) { any_bad <- TRUE; cat(basename(f), "lines", paste(h, collapse=","), "\n") }
}
if (!any_bad) cat("No dangerous-block characters anywhere in man/.\n")
