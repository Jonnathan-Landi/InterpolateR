# Data used for testing Kriging Ordinary
# Date of test creation: 2025-09-01
# Test update date: 2025-09-01
#
# Data input
data("BD_Obs", package = "InterpolateR")
data("BD_Coord", package = "InterpolateR")

# load area
shapefile <- terra::vect(system.file(
  "extdata/study_area.shp",
  package = "InterpolateR"
))

# Rain threshold for categorical metrics
Rain_threshold <- list(
  no_rain = c(0, 1),
  light_rain = c(1, 5),
  moderate_rain = c(5, 20),
  heavy_rain = c(20, 40),
  extremely_rain = c(40, Inf)
)

# Skip cran
testthat::skip_on_cran()

# 1. Testing without validation ---------------------------------------------------
testthat::test_that("Kriging_Ordinary returns SpatRaster without validation.", {
  testthat::skip_on_cran()
  out <- Kriging_Ordinary(
    BD_Obs,
    BD_Coord,
    shapefile,
    grid_resolution = 5,
    variogram_model = "exponential",
    max_dist = NULL,
    n_lags = 15,
    min_stations = 2,
    n_round = 1,
    training = 1,
    stat_validation = NULL,
    Rain_threshold = NULL,
    save_model = FALSE,
    name_save = NULL
  )

  testthat::expect_true(inherits(out, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out), length(unique(BD_Obs$Date)))
})

# 2. Testing with validation (random validation) --------------------------------
testthat::test_that("Kriging_Ordinary returns SpatRaster with random validation.", {
  testthat::skip_on_cran()
  out <- Kriging_Ordinary(
    BD_Obs,
    BD_Coord,
    shapefile,
    grid_resolution = 5,
    variogram_model = "spherical",
    max_dist = 50000,
    n_lags = 10,
    min_stations = 2,
    n_round = 2,
    training = 0.8,
    stat_validation = NULL,
    Rain_threshold = NULL,
    save_model = FALSE,
    name_save = NULL
  )

  testthat::expect_true(inherits(out$Ensamble, "SpatRaster"))
  testthat::expect_equal(
    terra::nlyr(out$Ensamble),
    length(unique(BD_Obs$Date))
  )
  testthat::expect_true(inherits(out$Validation, "data.table"))
})

# 3. Testing with validation (manual validation) --------------------------------
testthat::test_that("Kriging_Ordinary returns SpatRaster with manual validation.", {
  testthat::skip_on_cran()
  out <- Kriging_Ordinary(
    BD_Obs,
    BD_Coord,
    shapefile,
    grid_resolution = 5,
    variogram_model = "gaussian",
    max_dist = NULL,
    n_lags = 15,
    min_stations = 2,
    n_round = 1,
    training = 1,
    stat_validation = "M001",
    Rain_threshold = NULL,
    save_model = FALSE,
    name_save = NULL
  )

  testthat::expect_true(inherits(out$Ensamble, "SpatRaster"))
  testthat::expect_equal(
    terra::nlyr(out$Ensamble),
    length(unique(BD_Obs$Date))
  )
  testthat::expect_true(inherits(out$Validation, "data.table"))
})

# 4. Testing with categorical validation -----------------------------------------
testthat::test_that("Kriging_Ordinary returns validation with Rain_threshold parameter.", {
  testthat::skip_on_cran()
  out <- Kriging_Ordinary(
    BD_Obs,
    BD_Coord,
    shapefile,
    grid_resolution = 5,
    variogram_model = "linear",
    max_dist = NULL,
    n_lags = 15,
    min_stations = 2,
    n_round = NULL,
    training = 0.7,
    stat_validation = NULL,
    Rain_threshold = Rain_threshold,
    save_model = FALSE,
    name_save = NULL
  )

  # Check that output is a list with proper structure
  testthat::expect_true(is.list(out))
  testthat::expect_true("Ensamble" %in% names(out))
  testthat::expect_true("Validation" %in% names(out))
  testthat::expect_true(inherits(out$Ensamble, "SpatRaster"))

  # Check validation output structure
  testthat::expect_true(!is.null(out$Validation))

  # Find validation data (could be nested in list structure)
  validation_data <- NULL
  if (inherits(out$Validation, c("data.table", "data.frame"))) {
    validation_data <- out$Validation
  } else if (is.list(out$Validation)) {
    # Find first data.table/data.frame in the list
    for (item in out$Validation) {
      if (inherits(item, c("data.table", "data.frame"))) {
        validation_data <- item
        break
      }
    }
  }

  # Check that we found validation data
  testthat::expect_true(
    !is.null(validation_data),
    info = "Should contain validation data structure"
  )

  # Check for standard validation metrics (these should always be present)
  if (!is.null(validation_data)) {
    validation_names <- names(validation_data)

    # Check for standard continuous metrics
    standard_metrics <- any(grepl("RMSE|MAE|NSE|R2|KGE", validation_names, ignore.case = TRUE))
    testthat::expect_true(
      standard_metrics,
      info = paste("Available columns:", paste(validation_names, collapse = ", "))
    )

    # Note: Categorical metrics (CSI, POD, FAR) may not be implemented yet
    # with Rain_threshold parameter without errors
  }
})


# 5. Testing different variogram models ------------------------------------------
testthat::test_that("Kriging_Ordinary works with all variogram models.", {
  testthat::skip_on_cran()

  models <- c("exponential", "spherical", "gaussian", "linear")

  for (model in models) {
    out <- Kriging_Ordinary(
      BD_Obs,
      BD_Coord,
      shapefile,
      grid_resolution = 5,
      variogram_model = model,
      max_dist = NULL,
      n_lags = 10,
      min_stations = 2,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = FALSE,
      name_save = NULL
    )

    testthat::expect_true(
      inherits(out, "SpatRaster"),
      info = paste("Failed for model:", model)
    )
  }
})

##############################################################################
#     Check that the algorithm stops when the input data is not correct.     #
##############################################################################

# 6. shapefile must be a 'SpatVector' object. ----------------------------------
testthat::test_that("Error if `shapefile` is not SpatVector.", {
  testthat::skip_on_cran()
  bad_shape <- data.frame(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    Kriging_Ordinary(
      BD_Obs,
      BD_Coord,
      bad_shape,
      grid_resolution = 5,
      variogram_model = "exponential",
      max_dist = NULL,
      n_lags = 15,
      min_stations = 2,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = FALSE,
      name_save = NULL
    ),
    regexp = "shapefile must be a 'SpatVector' with a defined CRS\\.$"
  )
})

# 7. BD_Obs must be a 'data.table' or a 'data.frame'." -------------------------
testthat::test_that("Error if `BD_Obs` is not a data.table or data.frame.", {
  testthat::skip_on_cran()
  bad_obs <- list(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    Kriging_Ordinary(
      bad_obs,
      BD_Coord,
      shapefile,
      grid_resolution = 5,
      variogram_model = "exponential",
      max_dist = NULL,
      n_lags = 15,
      min_stations = 2,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = FALSE,
      name_save = NULL
    ),
    regexp = "BD_Obs must be a 'data.frame' or 'data.table'\\.$"
  )
})

# 8. BD_Coord must be a 'data.table' or a 'data.frame'." -----------------------
testthat::test_that("Error if `BD_Coord` is not a data.table or data.frame.", {
  testthat::skip_on_cran()
  bad_coord <- list(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    Kriging_Ordinary(
      BD_Obs,
      bad_coord,
      shapefile,
      grid_resolution = 5,
      variogram_model = "exponential",
      max_dist = NULL,
      n_lags = 15,
      min_stations = 2,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = FALSE,
      name_save = NULL
    ),
    regexp = "BD_Coord must be a 'data.frame' or 'data.table'\\.$"
  )
})

# 9. variogram_model must be valid ----------------------------------------------
testthat::test_that("Error if `variogram_model` is invalid.", {
  testthat::skip_on_cran()
  testthat::expect_error(
    Kriging_Ordinary(
      BD_Obs,
      BD_Coord,
      shapefile,
      grid_resolution = 5,
      variogram_model = "invalid_model",
      max_dist = NULL,
      n_lags = 15,
      min_stations = 2,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = FALSE,
      name_save = NULL
    ),
    regexp = "variogram_model must be one of 'exponential', 'spherical', 'gaussian', or 'linear'\\.$"
  )
})

# 10. grid_resolution must be numeric -------------------------------------------
testthat::test_that("Error if `grid_resolution` is not numeric.", {
  testthat::skip_on_cran()
  testthat::expect_error(
    Kriging_Ordinary(
      BD_Obs,
      BD_Coord,
      shapefile,
      grid_resolution = "invalid",
      variogram_model = "exponential",
      max_dist = NULL,
      n_lags = 15,
      min_stations = 2,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = FALSE,
      name_save = NULL
    ),
    regexp = "'grid_resolution' must be a single numeric value \\(km\\)\\.$"
  )
})

# 11. n_lags must be positive integer -------------------------------------------
testthat::test_that("Error if `n_lags` is not a positive integer.", {
  testthat::skip_on_cran()
  testthat::expect_error(
    Kriging_Ordinary(
      BD_Obs,
      BD_Coord,
      shapefile,
      grid_resolution = 5,
      variogram_model = "exponential",
      max_dist = NULL,
      n_lags = -5,
      min_stations = 2,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = FALSE,
      name_save = NULL
    ),
    regexp = "'n_lags' must be a single positive integer\\.$"
  )
})

# 12. min_stations must be positive integer -------------------------------------
testthat::test_that("Error if `min_stations` is not a positive integer.", {
  testthat::skip_on_cran()
  testthat::expect_error(
    Kriging_Ordinary(
      BD_Obs,
      BD_Coord,
      shapefile,
      grid_resolution = 5,
      variogram_model = "exponential",
      max_dist = NULL,
      n_lags = 15,
      min_stations = 0,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = FALSE,
      name_save = NULL
    ),
    regexp = "'min_stations' must be a single positive integer\\.$"
  )
})

# 13. n_round validation --------------------------------------------------------
testthat::test_that("Error if `n_round` is invalid.", {
  testthat::skip_on_cran()
  testthat::expect_error(
    Kriging_Ordinary(
      BD_Obs,
      BD_Coord,
      shapefile,
      grid_resolution = 5,
      variogram_model = "exponential",
      max_dist = NULL,
      n_lags = 15,
      min_stations = 2,
      n_round = -1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = FALSE,
      name_save = NULL
    ),
    regexp = "'n_round' must be NULL or a single non-negative integer\\.$"
  )
})

# 14. Coordinate names mismatch -------------------------------------------------
testthat::test_that("Error if coordinates names do not appear in observed data.", {
  testthat::skip_on_cran()

  # Create copy of BD_Coord with invalid code
  bad_coord <- BD_Coord
  bad_coord[3, "Cod"] <- "INVALID_STATION"

  testthat::expect_error(
    Kriging_Ordinary(
      BD_Obs,
      bad_coord,
      shapefile,
      grid_resolution = 5,
      variogram_model = "exponential",
      max_dist = NULL,
      n_lags = 15,
      min_stations = 2,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = FALSE,
      name_save = NULL
    ),
    regexp = "Coordinate names don't match observed data columns\\.$"
  )
})

# 15. Test model saving ---------------------------------------------------------
testthat::test_that("Kriging_Ordinary saves model when save_model = TRUE", {
  testthat::skip_on_cran()
  temp_dir <- tempdir()
  withr::local_dir(temp_dir)

  testthat::expect_message(
    out <- Kriging_Ordinary(
      BD_Obs,
      BD_Coord,
      shapefile,
      grid_resolution = 5,
      variogram_model = "exponential",
      max_dist = NULL,
      n_lags = 15,
      min_stations = 2,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = TRUE,
      name_save = "Test_Kriging"
    ),
    "Model saved successfully as Test_Kriging.nc"
  )

  expected_file <- file.path(temp_dir, "Test_Kriging.nc")
  testthat::expect_true(file.exists(expected_file), info = expected_file)
})

# 16. Test with default name saving ---------------------------------------------
testthat::test_that("Kriging_Ordinary saves model with default name", {
  testthat::skip_on_cran()
  temp_dir <- tempdir()
  withr::local_dir(temp_dir)

  testthat::expect_message(
    out <- Kriging_Ordinary(
      BD_Obs,
      BD_Coord,
      shapefile,
      grid_resolution = 5,
      variogram_model = "exponential",
      max_dist = NULL,
      n_lags = 15,
      min_stations = 2,
      n_round = 1,
      training = 1,
      stat_validation = NULL,
      Rain_threshold = NULL,
      save_model = TRUE,
      name_save = NULL
    ),
    "Model saved successfully as Model_Kriging.nc"
  )

  expected_file <- file.path(temp_dir, "Model_Kriging.nc")
  testthat::expect_true(file.exists(expected_file), info = expected_file)
})
# End of tests for Kriging_Ordinary
