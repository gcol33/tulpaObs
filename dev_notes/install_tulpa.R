repo <- "C:/Users/Gilles Colling/Documents/dev/tulpa"
devtools::document(repo, quiet = TRUE)
devtools::install(repo, quick = TRUE, upgrade = FALSE, quiet = TRUE)
