t0 <- Sys.time()
devtools::install("../tulpa", quick = TRUE, upgrade = FALSE, quiet = TRUE)
cat("tulpa installed:", as.character(utils::packageVersion("tulpa")), "\n")
cat("elapsed:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
