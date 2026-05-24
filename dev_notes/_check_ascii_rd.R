# Report any non-ASCII char (codepoint + glyph + line) in the regenerated Rd
# files touched this session, so forbidden Unicode (math ops / superscripts /
# arrows / letterlike) is caught before it breaks the PDF manual. Written as a
# .R file (inline -e segfaults on Windows per CLAUDE.md).
files <- c(
  "C:/Users/Gilles Colling/Documents/dev/tulpa/man/tulpa_nested_laplace.Rd",
  "C:/Users/Gilles Colling/Documents/dev/tulpaObs/man/predict.tobs_fit.Rd"
)
for (f in files) {
  ls <- readLines(f, warn = FALSE, encoding = "UTF-8")
  cat("===", basename(f), "===\n")
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
