# tulpa's C++ test fixtures, reached from this suite.
#
# A few tests compare an R-side construction of this package against tulpa's
# compiled kernel through a probe entry point (an NNGP log-density twin, the
# intrinsic-ICAR log-prior). tulpa does not ship those entry points: they live
# in its tests/testthat/cpp/<name>.cpp and its own loader, load_cpp_fixture(),
# compiles one against the engine's objects from the tulpa source tree. These
# tests call that same loader on the sibling tulpa checkout this package is
# built against (`TULPA_SOURCE`, else ../tulpa next to this package's source),
# so the kernel is tulpa's own and nothing is copied here. Without a source tree
# or a C++ toolchain the block skips, as tulpa's own fixture blocks do.

.tobs_tulpa_fixture_state <- new.env(parent = emptyenv())

.tobs_tulpa_source_root <- function() {
  cands <- c(Sys.getenv("TULPA_SOURCE", ""),
             file.path(testthat::test_path("..", "..", ".."), "tulpa"))
  for (d in cands[nzchar(cands)]) {
    desc <- file.path(d, "DESCRIPTION")
    if (file.exists(desc) &&
        file.exists(file.path(d, "tests", "testthat", "helper-cpp-fixtures.R")) &&
        identical(unname(read.dcf(desc, "Package")[1L, 1L]), "tulpa")) {
      return(normalizePath(d, winslash = "/"))
    }
  }
  NULL
}

# The functions of tulpa fixture `name`, as a list. Skips in the fast tier, on
# CRAN, and where the tulpa source tree is not on disk.
tulpa_fixture <- function(name) {
  skip_if_fast()
  skip_on_cran()
  hit <- .tobs_tulpa_fixture_state[[name]]
  if (!is.null(hit)) return(hit)
  root <- .tobs_tulpa_source_root()
  skip_if(is.null(root), "tulpa source tree not found (set TULPA_SOURCE)")
  td <- file.path(root, "tests", "testthat")
  loader <- new.env(parent = globalenv())
  sys.source(file.path(td, "helper-tiers.R"), envir = loader)
  sys.source(file.path(td, "helper-cpp-fixtures.R"), envir = loader)
  out <- new.env(parent = emptyenv())
  local({
    # tulpa's loader resolves the fixture and its source root through
    # testthat::test_path(), so it runs with the test directory at tulpa's.
    testthat::local_test_directory(td)
    loader$load_cpp_fixture(name, env = out)
  })
  fns <- as.list(out)
  .tobs_tulpa_fixture_state[[name]] <- fns
  fns
}
