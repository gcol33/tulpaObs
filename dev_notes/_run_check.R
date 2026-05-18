suppressPackageStartupMessages({
  library(devtools)
})
res <- devtools::check(args = "--no-manual", error_on = "warning",
                       quiet = FALSE)
print(res)
