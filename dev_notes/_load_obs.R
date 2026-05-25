# Scratch: recompile tulpaObs C++ against the freshly-installed tulpa headers
# and load. Isolates the cross-package ABI/linkage check from test logic.
pkg <- "C:/Users/Gilles Colling/Documents/dev/tulpaObs"
Rcpp::compileAttributes(pkg)
suppressWarnings(suppressMessages(devtools::load_all(pkg, quiet = TRUE)))
cat("LOAD_OK  tulpa=", as.character(packageVersion("tulpa")),
    " tulpaObs=", as.character(packageVersion("tulpaObs")), "\n", sep = "")
