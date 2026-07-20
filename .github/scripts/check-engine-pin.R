#!/usr/bin/env Rscript

# Assert that the declared engine versions, the pinned Remotes tags, and the
# versions actually installed all agree, and log all three every run.
#
# tulpaObs is a thin statistical layer over the tulpa engine and links against
# its headers (LinkingTo: tulpa), so a version the package does not declare can
# change fit behaviour, not only error. DESCRIPTION states the contract twice --
# an Imports floor and an exact Remotes tag -- and those two plus the installed
# version drifted apart once already (issue #150: Imports/Remotes said 0.0.85
# while every local measurement ran on 0.0.92).
#
# Fails loudly rather than warning: an engine skew that only shows up as a NOTE
# is the failure mode this guard exists to remove.

desc <- read.dcf("DESCRIPTION")

# Engines whose version is part of the package contract. tulpa is a LinkingTo
# dependency (ABI-coupled); tulpaMesh supplies the SPDE mesh surface.
engines <- c("tulpa", "tulpaMesh")

# "tulpa (>= 0.0.93)" -> "0.0.93", for the one entry naming this engine.
floor_version <- function(field, pkg) {
  if (!field %in% colnames(desc)) return(NA_character_)
  entries <- trimws(strsplit(desc[1, field], ",")[[1]])
  hit <- entries[grepl(paste0("^", pkg, "\\b"), entries)]
  if (length(hit) != 1L) return(NA_character_)
  sub(".*\\(>=\\s*([^)]+)\\).*", "\\1", hit)
}

# "gcol33/tulpa@v0.0.93" -> "0.0.93", for the one entry naming this engine.
remote_version <- function(pkg) {
  if (!"Remotes" %in% colnames(desc)) return(NA_character_)
  entries <- trimws(strsplit(desc[1, "Remotes"], ",")[[1]])
  hit <- entries[grepl(paste0("/", pkg, "@"), entries)]
  if (length(hit) != 1L) return(NA_character_)
  sub(".*@v?", "", hit)
}

problems <- character(0)

for (pkg in engines) {
  imports <- floor_version("Imports", pkg)
  remotes <- remote_version(pkg)
  installed <- tryCatch(as.character(utils::packageVersion(pkg)),
                        error = function(e) NA_character_)

  cat(sprintf("%-10s  Imports >= %-8s  Remotes @v%-8s  installed %s\n",
              pkg, imports, remotes, installed))

  if (is.na(imports)) {
    problems <- c(problems, sprintf(
      "%s: no version floor in DESCRIPTION Imports. Declare 'Imports: %s (>= X.Y.Z)'.",
      pkg, pkg))
    next
  }
  if (is.na(remotes)) {
    problems <- c(problems, sprintf(
      "%s: no pinned tag in DESCRIPTION Remotes. Declare 'Remotes: gcol33/%s@vX.Y.Z'.",
      pkg, pkg))
    next
  }
  if (is.na(installed)) {
    problems <- c(problems, sprintf("%s: declared but not installed.", pkg))
    next
  }

  # The two declarations must agree: Remotes decides what a fresh install
  # resolves, Imports decides what the package claims to support. If the tag is
  # below the floor, a fresh install produces an unusable package; if it is
  # above, the floor understates what is actually being tested.
  if (!identical(imports, remotes)) {
    problems <- c(problems, sprintf(
      "%s: DESCRIPTION disagrees with itself -- Imports floor %s, Remotes tag v%s. Bump both together.",
      pkg, imports, remotes))
  }

  # What CI resolved must be what the package pins, or CI is not testing what
  # the package ships.
  if (package_version(installed) != package_version(remotes)) {
    problems <- c(problems, sprintf(
      "%s: installed %s but Remotes pins v%s. CI is not testing the pinned engine.",
      pkg, installed, remotes))
  }
}

if (length(problems)) {
  cat("\nEngine pin check FAILED:\n")
  cat(paste0("  - ", problems, collapse = "\n"), "\n")
  quit(status = 1L)
}

cat("\nEngine pin check passed: declared, pinned, and installed versions agree.\n")
