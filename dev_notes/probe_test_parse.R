# Parse-check the new SLA joint test file.
#
# Confirms the file is syntactically valid R *without* running the tests
# (the implementation it targets may not be in place yet).

repo <- "C:/Users/Gilles Colling/Documents/dev/tulpaObs"
path <- file.path(repo, "tests", "testthat", "test-sla-cover-joint.R")

cat("Parsing:", path, "\n")
parsed <- parse(path)
cat("OK: parsed", length(parsed), "top-level expressions.\n")
