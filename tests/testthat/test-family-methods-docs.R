# Every family help page names the engines .tobs_family_methods accepts for it,
# as `method = "<name>"` strings under its Engines section. The table is what
# tobs() validates against, so the page is held to it in both directions: a
# method the table lists has to be named on the page, and a `method = "<x>"`
# the table lacks must not appear anywhere on it (the page would then describe
# a route the family refuses). No fit is run; the Rd is parsed from the source
# `man/` under load_all() and from the installed help database under check.

.family_rd_text <- function(name) {
  man <- testthat::test_path("..", "..", "man")
  rd <- if (dir.exists(man)) {
    tools::parse_Rd(file.path(man, paste0(name, ".Rd")))
  } else {
    tools::Rd_db("tulpaObs")[[paste0(name, ".Rd")]]
  }
  paste(unlist(rd), collapse = "")
}

.family_rd_methods <- function(txt) {
  hits <- regmatches(txt, gregexpr('method = "[a-z_]+"', txt))[[1L]]
  sort(unique(sub('^method = "([a-z_]+)"$', "\\1", hits)))
}

test_that("every family in .tobs_family_methods has a help page", {
  for (name in names(.tobs_family_methods)) {
    expect_true(nzchar(.family_rd_text(name)),
                label = sprintf("man/%s.Rd parsed to non-empty text", name))
  }
})

test_that("each family help page names exactly the engines the table lists", {
  for (name in names(.tobs_family_methods)) {
    found <- .family_rd_methods(.family_rd_text(name))
    expect_identical(
      found, sort(.tobs_family_methods[[name]]),
      label = sprintf("method strings on man/%s.Rd", name),
      expected.label = sprintf(".tobs_family_methods$%s", name))
  }
})
