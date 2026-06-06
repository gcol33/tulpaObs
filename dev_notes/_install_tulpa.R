options(warn = 1)
cat("== install tulpa (R + NUTS/EM cpp changed) ==\n")
devtools::install("../tulpa", quick = FALSE, upgrade = FALSE, quiet = TRUE)
cat("== TULPA INSTALL DONE ==\n")
