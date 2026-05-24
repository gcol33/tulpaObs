# ASCII-scan all staged .Rd files in the nested-Laplace bundle so no forbidden
# Unicode (math ops / superscripts / arrows / letterlike) breaks the PDF manual.
files <- c(
  "man/predict.tobs_fit.Rd", "man/dot-tobs_em_nested_laplace.Rd",
  "man/dot-tobs_laplace.Rd", "man/tobs.Rd"
)
for (f in files) {
  ls <- readLines(f, warn = FALSE, encoding = "UTF-8")
  cat("===", f, "===\n")
  hit <- FALSE
  for (i in seq_along(ls)) {
    cp <- utf8ToInt(enc2utf8(ls[i]))
    bad <- which(cp > 127L)
    if (length(bad)) {
      hit <- TRUE
      for (b in bad)
        cat(sprintf("  line %d: U+%04X '%s'\n", i, cp[b], intToUtf8(cp[b])))
    }
  }
  if (!hit) cat("  clean ASCII\n")
}
