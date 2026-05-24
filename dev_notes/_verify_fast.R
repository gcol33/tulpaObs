setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
files <- list.files("tests/testthat", pattern = "^test-.*\\.R$", full.names = TRUE)
bad <- 0L
for (f in files) {
  L <- readLines(f, warn = FALSE)
  oc <- which(grepl("^\\s*skip_on_cran\\(\\)\\s*$", L))
  for (i in oc) {
    nxt <- if (i < length(L)) L[i + 1L] else ""
    if (!grepl("^\\s*skip_if_fast\\(\\)\\s*$", nxt)) {
      cat("UNPAIRED:", basename(f), "line", i, "\n"); bad <- bad + 1L
    } else {
      ind_oc <- sub("skip.*", "", L[i]); ind_sf <- sub("skip.*", "", nxt)
      if (!identical(ind_oc, ind_sf))
        cat("INDENT-MISMATCH:", basename(f), "line", i, "\n")
    }
  }
}
cat("unpaired skip_on_cran:", bad, "\n")
cat("--- sla-cover-joint (4-space) ---\n")
L <- readLines("tests/testthat/test-sla-cover-joint.R", warn = FALSE)
oc <- which(grepl("skip_on_cran", L))[1]
cat("[", L[oc], "][", L[oc + 1L], "]\n", sep = "")
