cat("=== msNMix args ===\n")
print(args(spAbundance::msNMix))
cat("\n=== msNMix help (Usage + Arguments, first lines) ===\n")
db <- tools::Rd_db("spAbundance")
rd <- db[["msNMix.Rd"]]
if (!is.null(rd)) {
  tmp <- tempfile(fileext = ".txt")
  tools::Rd2txt(rd, out = tmp)
  cat(paste(readLines(tmp), collapse = "\n"))
}
