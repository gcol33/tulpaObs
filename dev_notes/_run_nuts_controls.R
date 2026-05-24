# Runner: regenerate docs, load, and exercise the renamed `visits` arg + the
# new NUTS controls / convergence diagnostics. (prefix `_` = runner.)
setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
options(crayon.enabled = FALSE)
suppressMessages(library(devtools))

cat("== document ==\n")
document(pkg = ".", quiet = TRUE)

cat("== load_all ==\n")
load_all(".", quiet = TRUE, export_all = FALSE)

cat("== targeted tests ==\n")
res <- test(filter = "nuts-controls|visit-data|issue8|methods",
            reporter = "summary", stop_on_failure = FALSE)
print(res)
