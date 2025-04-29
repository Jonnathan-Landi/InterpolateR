# Data used for create data
# Date of test creation: 2025-04-29
# Test update date: 2025-04-29

file.path <- system.file("extdata/Folds_ejs_create_data", package = "InterpolateR")
skip_on_cran()
# 1. Thest with default parameters
testthat::test_that("Create data Method with basic input", {
  out  <- create_data(file.path,
    Start_date = "2015-01-01", End_Date = "2015-03-01",
    ncores = NULL, max.na = NULL
  )
  testthat::expect_true(inherits(out, "data.table"))
})

# 2. Test with different parameters (max.na set to 10)
testthat::test_that("Create data Method with 10% NA.", {
  out  <- create_data(file.path,
                      Start_date = "2015-01-01", End_Date = "2015-03-01",
                      ncores = NULL, max.na = 10
  )
  testthat::expect_true(inherits(out$data, "data.table"))
  testthat::expect_true(inherits(out$Na_stations, "data.table"))
})

# 3. Test with paralel processing
testthat::test_that("Create data Method with parallel processing", {
  out  <- create_data(file.path,
                      Start_date = "2015-01-01", End_Date = "2015-03-01",
                      ncores = 2, max.na = NULL
  )
  testthat::expect_true(inherits(out, "data.table"))
})
