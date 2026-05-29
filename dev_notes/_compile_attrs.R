# Regenerate Rcpp exports + roxygen for tulpaObs after adding the native
# single-species grouped-RE oracle (src/nmix_re_oracle.{h,cpp}).
setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
Rcpp::compileAttributes()
devtools::document()
cat("compile + document OK\n")
