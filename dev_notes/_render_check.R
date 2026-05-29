suppressMessages(devtools::load_all("."))
f <- commandArgs(trailingOnly = TRUE)[1]
tmp <- tempfile(fileext = ".R")
knitr::purl(f, output = tmp, quiet = TRUE)
# drop the library(tulpaObs) line (already loaded via load_all)
src <- readLines(tmp)
src <- src[!grepl("^library\\(tulpaObs\\)", src)]
writeLines(src, tmp)
cat("==== running", basename(f), "====\n")
source(tmp, echo = FALSE)
cat("\n==== OK:", basename(f), "====\n")
