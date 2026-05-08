library(testthat)

Sys.setenv(CARF_PROJECT_ROOT = normalizePath(getwd(), mustWork = TRUE))
test_dir("tests/testthat", reporter = "summary")
