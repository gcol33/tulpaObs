suppressMessages(devtools::document(quiet = TRUE))
cat("DOCUMENT OK\n")
# Surface any non-ASCII that crept into man/ from the roxygen edits.
mans <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
bad <- character(0)
for (f in mans) {
  lines <- readLines(f, warn = FALSE)
  hits <- grep("[^\x01-\x7f]", lines)
  if (length(hits)) bad <- c(bad, sprintf("%s:%d", basename(f), hits))
}
if (length(bad)) { cat("NON-ASCII:\n"); cat(bad, sep = "\n"); cat("\n") } else
  cat("man/ is ASCII-clean\n")
