# Quick R-only reinstall of tulpa (no clean_dll): picks up R-side changes
# (the K_max cap in tulpa_nmix_laplace_re) without recompiling C++.
devtools::install("C:/Users/Gilles Colling/Documents/dev/tulpa",
                  quick = TRUE, upgrade = FALSE, reload = FALSE)
cat("=== tulpa quick-installed ===\n")
