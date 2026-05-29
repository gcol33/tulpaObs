x <- readLines("vignettes/integrated-occupancy.Rmd")
cat("em-dashes:", sum(lengths(regmatches(x, gregexpr("—", x)))), "\n")
banned <- c("delve","robust","comprehensive","seamless","leverage","crucial",
            "important to note","intricate","realm","tapestry","testament",
            "underscore","showcase","boast","navigate the","in today",
            "ever-evolving","game-chang","cutting-edge","not just","but rather")
for (b in banned) {
  n <- sum(grepl(b, x, ignore.case = TRUE))
  if (n > 0) cat("BANNED", b, ":", n, "\n")
}
# prose word count: drop code chunks, yaml, headers, fences
in_code <- FALSE; prose <- character()
for (ln in x) {
  if (grepl("^```", ln)) { in_code <- !in_code; next }
  if (in_code) next
  if (grepl("^---$|^title:|^output:|^vignette:|^\\s*%|^\\s*$", ln)) next
  prose <- c(prose, ln)
}
# strip markdown headers/table rows
prose <- prose[!grepl("^#", prose)]
prose <- prose[!grepl("^\\|", prose)]
words <- unlist(strsplit(paste(prose, collapse = " "), "\\s+"))
words <- words[nzchar(words)]
cat("prose word count (approx):", length(words), "\n")
