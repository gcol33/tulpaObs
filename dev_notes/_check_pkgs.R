for (p in c("spAbundance", "unmarked", "coda")) {
  ok <- requireNamespace(p, quietly = TRUE)
  ver <- if (ok) as.character(utils::packageVersion(p)) else "NA"
  cat(sprintf("%-14s installed=%s version=%s\n", p, ok, ver))
}
