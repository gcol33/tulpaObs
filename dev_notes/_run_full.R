options(crayon.enabled = FALSE)
suppressMessages(library(devtools))
cat("== document tulpa (regenerate NAMESPACE for mcmc_diagnostics export) ==\n")
document("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE, export_all = FALSE)
cat("exported?", "mcmc_diagnostics" %in% getNamespaceExports("tulpa"),
    "/", "select_main_params" %in% getNamespaceExports("tulpa"), "\n")
load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE, export_all = FALSE)
setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
res <- as.data.frame(test(reporter = "summary", stop_on_failure = FALSE))
cat("\n=== TOTALS ===\n")
cat("failed:", sum(res$failed), " errors:", sum(res$error),
    " warnings:", sum(res$warning), " skipped:", sum(res$skipped),
    " passed:", sum(res$passed), "\n")
bad <- res[res$failed > 0 | res$error > 0 | res$warning > 0,
           c("file","context","failed","error","warning")]
if (nrow(bad)) { cat("\n=== FILES W/ FAIL|ERROR|WARN ===\n"); print(bad) } else cat("ALL GREEN\n")
