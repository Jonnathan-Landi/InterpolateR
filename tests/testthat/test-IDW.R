# Data used for testing
# Date: 2025-04-29
data("BD_Obs",   package = "InterpolateR")
data("BD_Coord", package = "InterpolateR")
shapefile <- terra::vect(system.file("extdata", "study_area.shp", package = "InterpolateR"))
Rain_threshold <- list(
  no_rain       = c(0, 1),
  light_rain    = c(1, 5),
  moderate_rain = c(5, 20),
  heavy_rain    = c(20, 40),
  extremely_rain= c(40, Inf)
)

skip_on_cran()

# 1. Testing with validation ---------------------------------------------------
testthat::test_that("IDW returns SpatRaster with full validation.", {
  out <- IDW(BD_Obs, BD_Coord, shapefile,
             grid_resolution = 5, p = 2,
             n_round = 1, training = 1,
             Rain_threshold = Rain_threshold,
             stat_validation = "M004",
             save_model = FALSE)
  testthat::expect_true(inherits(out$Ensamble, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out$Ensamble), length(unique(BD_Obs$Date)))
})

# 2. Testing without validation ------------------------------------------------
testthat::test_that("IDW works without validation", {
  out <- IDW(BD_Obs, BD_Coord, shapefile,
             grid_resolution = 5, p = 2,
             n_round = 1, training = 1,
             Rain_threshold = NULL,
             stat_validation = NULL,
             save_model = FALSE)
  testthat::expect_true(inherits(out,    "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out), length(unique(BD_Obs$Date)))
})

# 3. shapefile must be a 'SpatVector' object. ----------------------------------
testthat::test_that("Error if `shapefile` is not SpatVector.", {
  bad_shape <- data.frame(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    IDW(BD_Obs, BD_Coord, bad_shape,
        grid_resolution = 5, p = 2,
        n_round = 1, training = 1,
        Rain_threshold = NULL,
        stat_validation = NULL,
        save_model = FALSE),
    regexp = "shapefile must be a 'SpatVector' object\\.$"
  )
})

# 4. BD_Obs must be a 'data.table' or a 'data.frame'." -------------------------
testthat::test_that("Error if `BD_Obs` is not a data.table or data.frame.", {
  bad_obs <- list(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    IDW(bad_obs, BD_Coord, shapefile,
        grid_resolution = 5, p = 2,
        n_round = 1, training = 1,
        Rain_threshold = NULL,
        stat_validation = NULL,
        save_model = FALSE),
    regexp = "BD_Obs must be a 'data.table' or a 'data.frame'\\.$"
  )
})

# 5. BD_Coord must be a 'data.table' or a 'data.frame'." -----------------------
testthat::test_that("Error if `BD_Coord` is not a data.table or data.frame.", {
  bad_coord <- list(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    IDW(BD_Obs, bad_coord, shapefile,
        grid_resolution = 5, p = 2,
        n_round = 1, training = 1,
        Rain_threshold = NULL,
        stat_validation = NULL,
        save_model = FALSE),
    regexp = "BD_Coord must be a 'data.table' or a 'data.frame'\\.$"
  )
})

# 6. "The names of the coordinates do not appear in the observed data." --------
testthat::test_that("Error if coordinates names do not appear in observed data.", {
  bad_coord <- BD_Coord
  bad_coord[3,1] <- "x"
  testthat::expect_error(
    IDW(BD_Obs, bad_coord, shapefile,
        grid_resolution = 5, p = 2,
        n_round = 1, training = 1,
        Rain_threshold = NULL,
        stat_validation = NULL,
        save_model = FALSE),
    regexp = "The names of the coordinates do not appear in the observed data\\.$"
  )
})

# 7. "Save the model must be a logical value." ---------------------------------
testthat::test_that("IDW saves model when save_model = TRUE", {
  temp_dir <- tempdir()
  withr::local_dir(temp_dir)
  custom_name <- "test_model_IDW"
  expect_message(
    out <- IDW(
      BD_Obs, BD_Coord, shapefile,
      grid_resolution = 5, p = 2,
      n_round = 1, training = 1,
      Rain_threshold = NULL, stat_validation = NULL,
      save_model = TRUE, name_save = custom_name
    ),
    "Model saved successfully"
  )
  expected_file <- file.path(temp_dir, paste0(custom_name, ".nc"))
  testthat::expect_true(file.exists(expected_file), info = expected_file)
  testthat::expect_true(inherits(terra::rast(expected_file), "SpatRaster"))
})
