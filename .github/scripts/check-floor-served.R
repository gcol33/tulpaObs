#!/usr/bin/env Rscript

# Assert that every version floor DESCRIPTION declares is served, right now, by
# the repositories DESCRIPTION itself names -- so a floor is never committed
# before the version it demands exists to install.
#
# tulpa and tulpaMesh resolve from https://gcol33.r-universe.dev, declared in
# Additional_repositories. r-universe builds an engine's DEFAULT-BRANCH HEAD on
# its own poll schedule: a git tag and a GitHub release publish nothing, and R's
# resolver can install from neither. Between pushing an engine release and
# r-universe's next build there is therefore a window in which a raised floor is
# unsatisfiable, and every job that resolves dependencies fails with a solver
# error that names no version at all:
#
#   * deps::.: Can't install dependency tulpa (>= 0.0.117)
#
# This gate answers what that error does not: which floor, what is actually
# served, and from where. Run it BEFORE committing a floor bump (issue #180):
#
#   Rscript .github/scripts/check-floor-served.R
#
# Scope, against the sibling check-engine-pin.R: that one runs inside a job that
# has already resolved and installed, and asks whether what got installed is what
# the package declares. This one needs nothing installed, and asks whether the
# declaration is installable at all.
#
# Policy this encodes: DESCRIPTION currently pins tulpa by an exact Remotes tag
# (gcol33/tulpa@v<ver>) as well as an Imports floor, so a resolver that honours
# Remotes installs the tag directly and never hits this window for tulpa itself.
# The gate still matters for every OTHER hard dependency (tulpaMesh, loo, ...),
# none of which carry a Remotes tag, and for tulpa on a resolver that ignores
# Remotes and falls back to the Imports floor served from r-universe.
#
# Repositories come from DESCRIPTION rather than from options("repos"), so the
# gate measures what a user installing this package resolves, not what a
# particular workflow happened to configure.
#
# Unverifiable is reported as failure, never as a pass: a repository that cannot
# be reached leaves its packages unknown, and a floor that depends on it fails.
# A repository can only ever add supply, so one that is unreachable while every
# floor is already satisfied elsewhere cannot turn a pass into a failure -- that
# case is reported and the gate still passes.

options(timeout = 60)

fail <- function(...) {
  cat("\nFloor availability gate FAILED:\n")
  cat(paste0("  - ", c(...), collapse = "\n"), "\n", sep = "")
  quit(status = 1L)
}

if (!file.exists("DESCRIPTION")) {
  fail(sprintf(paste0("no DESCRIPTION in the working directory (%s).\n",
                      "    Run this from the package root: ",
                      "Rscript .github/scripts/check-floor-served.R"),
               normalizePath(".", winslash = "/")))
}
desc <- read.dcf("DESCRIPTION")

# Hard dependencies only. Suggests is not resolved by an install, so a floor
# there cannot make one fail.
dep_fields <- c("Depends", "Imports", "LinkingTo")

# Base packages ship with R and are never resolved from a repository; R itself
# is a floor on the interpreter, not on a package.
never_served <- c("R", rownames(installed.packages(lib.loc = .Library,
                                                   priority = "base")))

# "tulpa (>= 0.0.117)" -> package "tulpa", floor "0.0.117". An entry carrying no
# floor places no demand on any repository, so it is not gated here.
declared_floors <- function() {
  floors <- list()
  for (field in dep_fields) {
    if (!field %in% colnames(desc)) next
    entries <- trimws(strsplit(desc[1L, field], ",")[[1L]])
    for (entry in entries[nzchar(entries)]) {
      if (!grepl("\\(\\s*>=", entry)) next
      pkg <- trimws(sub("\\s*\\(.*$", "", entry))
      if (pkg %in% never_served) next
      ver <- trimws(sub(".*\\(\\s*>=\\s*([^)]+)\\).*", "\\1", entry))
      # The same package can be declared in more than one field; the highest
      # floor is the one an install has to satisfy.
      if (!is.null(floors[[pkg]]) &&
          package_version(floors[[pkg]]$floor) >= package_version(ver)) next
      floors[[pkg]] <- list(package = pkg, floor = ver, field = field)
    }
  }
  floors[order(names(floors))]
}

repositories <- function() {
  cran <- unname(getOption("repos")[["CRAN"]])
  if (is.null(cran) || !nzchar(cran) || identical(cran, "@CRAN@")) {
    cran <- "https://cloud.r-project.org"
  }
  extra <- character(0L)
  if ("Additional_repositories" %in% colnames(desc)) {
    extra <- trimws(strsplit(desc[1L, "Additional_repositories"], "[,[:space:]]+")[[1L]])
    extra <- extra[nzchar(extra)]
  }
  urls <- unique(c(cran, extra))
  data.frame(url = urls,
             kind = ifelse(urls %in% cran, "CRAN", "additional"),
             stringsAsFactors = FALSE)
}

# The source contrib index -- exactly the /src/contrib/PACKAGES the engine's
# universe is checked against by hand. Version is the same across build types,
# and source is the one index every platform serves.
query_repo <- function(url) {
  notes <- character(0L)
  index <- withCallingHandlers(
    tryCatch(
      utils::available.packages(repos = url, type = "source",
                                filters = "duplicates"),
      error = function(e) {
        notes <<- c(notes, conditionMessage(e))
        NULL
      }
    ),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (is.null(index) || nrow(index) == 0L) {
    return(list(ok = FALSE,
                versions = character(0L),
                note = if (length(notes)) notes[[1L]] else "index empty"))
  }
  list(ok = TRUE,
       versions = stats::setNames(index[, "Version"], rownames(index)),
       note = NA_character_)
}

floors <- declared_floors()
repos <- repositories()

cat("Floor availability gate -- is every declared floor installable today?\n\n")

if (!length(floors)) {
  cat("No versioned hard dependency declared; nothing to gate.\n")
  quit(status = 0L)
}

served <- lapply(repos$url, query_repo)
names(served) <- repos$url

cat("repositories\n")
for (i in seq_len(nrow(repos))) {
  hit <- served[[repos$url[i]]]
  cat(sprintf("  %-11s %-46s %s\n", repos$kind[i], repos$url[i],
              if (hit$ok) sprintf("%d packages", length(hit$versions))
              else sprintf("UNREACHABLE (%s)", hit$note)))
}

unreachable <- repos$url[!vapply(served, `[[`, logical(1L), "ok")]

# Best offer per package across every reachable repository: whichever serves the
# highest version is the one an install resolves the floor from.
best_offer <- function(pkg) {
  best <- NULL
  for (url in repos$url) {
    hit <- served[[url]]
    if (!hit$ok || !pkg %in% names(hit$versions)) next
    version <- unname(hit$versions[[pkg]])
    if (is.null(best) || package_version(version) > package_version(best$version)) {
      best <- list(url = url, version = version)
    }
  }
  best
}

name_width <- max(12L, max(nchar(names(floors))))
floor_width <- max(10L, max(nchar(vapply(floors, `[[`, character(1L), "floor"))))
row <- paste0("  %-", name_width, "s %-", floor_width, "s %-13s %s")

cat("\n")
cat(sprintf(paste0(row, "\n"), "package", "floor", "best served", "from"))

problems <- character(0L)

for (spec in floors) {
  offer <- best_offer(spec$package)

  if (is.null(offer)) {
    cat(sprintf(paste0(row, "\n"), spec$package, spec$floor, "NOT SERVED", "-"))
    problems <- c(problems, sprintf(
      paste0("%s: DESCRIPTION %s declares (>= %s), and no repository DESCRIPTION ",
             "names serves it at all.\n",
             "    Either the package is missing from every declared repository, ",
             "or the repository serving it is not declared in ",
             "Additional_repositories."),
      spec$package, spec$field, spec$floor))
    next
  }

  ok <- package_version(offer$version) >= package_version(spec$floor)
  cat(sprintf(paste0(row, "%s\n"), spec$package, spec$floor,
              offer$version, offer$url, if (ok) "" else "   <- BELOW FLOOR"))
  if (ok) next

  detail <- sprintf(
    paste0("%s: DESCRIPTION %s declares (>= %s), but the highest version served ",
           "anywhere is %s (%s).\n",
           "    Every job that resolves dependencies will fail with ",
           "\"Can't install dependency %s (>= %s)\"."),
    spec$package, spec$field, spec$floor, offer$version, offer$url,
    spec$package, spec$floor)

  if (grepl("r-universe\\.dev", offer$url)) {
    universe <- sub("^(https?://[^/]+).*$", "\\1", offer$url)
    detail <- paste0(detail, sprintf(
      paste0("\n    r-universe builds default-branch HEAD on its own poll ",
             "schedule, so no tag\n",
             "    and no GitHub release can publish %s -- only a push to the ",
             "engine's default\n",
             "    branch can, and only at the next poll. Wait for that build, ",
             "then re-run this\n",
             "    gate before committing the bump. What is built right now:\n",
             "      %s/api/packages/%s\n",
             "      %s/src/contrib/PACKAGES"),
      spec$floor, universe, spec$package, universe))
  }

  problems <- c(problems, detail)
}

# An unreachable repository can only withhold supply, so it matters exactly when
# a floor is left unsatisfied -- then the shortfall may be an outage rather than
# a missing build, and saying so is the difference between waiting and debugging.
if (length(unreachable) && length(problems)) {
  problems <- c(problems, sprintf(
    paste0("could not reach %s, so what it serves is unknown. Any shortfall ",
           "above may be that outage rather than a missing build."),
    paste(unreachable, collapse = ", ")))
}

if (length(problems)) fail(problems)

if (length(unreachable)) {
  cat(sprintf(paste0("\nNote: %s could not be reached. Every floor is served by ",
                     "a repository that could,\nso the result stands -- an ",
                     "unreachable repository can only add supply.\n"),
              paste(unreachable, collapse = ", ")))
}

cat("\nFloor availability gate passed: every declared floor is served today.\n")
