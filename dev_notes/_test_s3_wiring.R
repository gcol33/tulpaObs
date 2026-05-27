suppressMessages(devtools::load_all(".", quiet = TRUE))
testthat::test_file("tests/testthat/test-ms-abun.R",
                    desc = "ms_abun S3 methods work")
