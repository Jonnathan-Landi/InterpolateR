# Data used for testing Cressman
# Date of test creation: 2025-01-01
# Test update date: 2025-09-01

# Data input
data("BD_Obs", package = "InterpolateR")
data("BD_Coord", package = "InterpolateR")

# load area
shapefile <- terra::vect(system.file(
  "extdata/study_area.shp",
  package = "InterpolateR"
))

# Skip cran
testthat::skip_on_cran()

# 1. Testing without validation ---------------------------------------------------
testthat::test_that("Cressman returns SpatRaster without validation.", {
  testthat::skip_on_cran()
  out <- Cressman(
    BD_Obs,
    BD_Coord,
    shapefile,
    grid_resolution = 5,
    search_radius = 10,
    training = 1,
    stat_validation = NULL,
    Rain_threshold = NULL,
    save_model = FALSE
  )

  testthat::expect_true(inherits(out$`10 km`, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out$`10 km`), length(unique(BD_Obs$Date)))
})

# 2. Testing with validation (random validation) --------------------------------
testthat::test_that("Cressman returns SpatRaster with random validation.", {
  testthat::skip_on_cran()
  out <- Cressman(
    BD_Obs,
    BD_Coord,
    shapefile,
    grid_resolution = 5,
    search_radius = 10,
    training = 0.8,
    stat_validation = NULL,
    n_round = 2,
    Rain_threshold = NULL,
    save_model = FALSE
  )

  testthat::expect_true(inherits(out$Ensamble$`10 km`, "SpatRaster"))
  testthat::expect_equal(
    terra::nlyr(out$Ensamble$`10 km`),
    length(unique(BD_Obs$Date))
  )
  testthat::expect_true(inherits(out$Validation$`10 km`, "data.table"))
})

# 3. Testing with validation (manual validation) --------------------------------
testthat::test_that("Cressman returns SpatRaster with manual validation.", {
  testthat::skip_on_cran()
  out <- Cressman(
    BD_Obs,
    BD_Coord,
    shapefile,
    grid_resolution = 5,
    search_radius = 10,
    training = 1,
    stat_validation = "M001",
    Rain_threshold = NULL,
    save_model = FALSE
  )

  testthat::expect_true(inherits(out$Ensamble$`10 km`, "SpatRaster"))
  testthat::expect_equal(
    terra::nlyr(out$Ensamble$`10 km`),
    length(unique(BD_Obs$Date))
  )
  testthat::expect_true(inherits(out$Validation$`10 km`, "data.table"))
})
##############################################################################
#     Check that the algorithm stops when the input data is not correct.     #
##############################################################################
# 4. shapefile must be a 'SpatVector' object. ----------------------------------
testthat::test_that("Error if `shapefile` is not SpatVector.", {
  testthat::skip_on_cran()
  bad_shape <- data.frame(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    Cressman(
      BD_Obs,
      BD_Coord,
      bad_shape,
      grid_resolution = 5,
      search_radius = 10,
      training = 1,
      stat_validation = "M001",
      Rain_threshold = NULL,
      save_model = FALSE
    ),
    regexp = "shapefile must be a 'SpatVector' with a defined CRS\\.$"
  )
})

# 5. BD_Obs must be a 'data.table' or a 'data.frame'." -------------------------
testthat::test_that("Error if `BD_Obs` is not a data.table or data.frame.", {
  testthat::skip_on_cran()
  bad_obs <- list(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    Cressman(
      bad_obs,
      BD_Coord,
      shapefile,
      grid_resolution = 5,
      search_radius = 10,
      training = 1,
      stat_validation = "M001",
      Rain_threshold = NULL,
      save_model = FALSE
    ),
    regexp = "BD_Obs must be a 'data.table' or a 'data.frame'\\.$"
  )
})

# 6. BD_Coord must be a 'data.table' or a 'data.frame'." -----------------------
testthat::test_that("Error if `BD_Coord` is not a data.table or data.frame.", {
  testthat::skip_on_cran()
  bad_coord <- list(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    Cressman(
      BD_Obs,
      bad_coord,
      shapefile,
      grid_resolution = 5,
      search_radius = 10,
      training = 1,
      stat_validation = "M001",
      Rain_threshold = NULL,
      save_model = FALSE
    ),
    regexp = "BD_Coord must be a 'data.table' or a 'data.frame'\\.$"
  )
})

# 7. "The names of the coordinates do not appear in the observed data." --------
testthat::test_that("Error if coordinates names do not appear in observed data.", {
  testthat::skip_on_cran()

  # Crear copia de BD_Coord con un código inválido
  bad_coord <- BD_Coord
  bad_coord[3, "Cod"] <- "x"

  # Esperamos que falle con mensaje que mencione "do not appear" y el código "x"
  testthat::expect_error(
    Cressman(
      BD_Obs,
      bad_coord,
      shapefile,
      grid_resolution = 5,
      search_radius = 10,
      training = 1,
      stat_validation = "M001",
      Rain_threshold = NULL,
      save_model = FALSE
    ),
    regexp = "do not appear.*x"
  )
})

# 8. "Save the model must be a logical value." ---------------------------------
testthat::test_that("Cressman saves model when save_model = TRUE (default name)", {
  testthat::skip_on_cran()
  temp_dir <- tempdir()
  withr::local_dir(temp_dir)
  expect_message(
    out <- Cressman(
      BD_Obs,
      BD_Coord,
      shapefile,
      grid_resolution = 5,
      search_radius = 10,
      training = 1,
      stat_validation = "M001",
      Rain_threshold = NULL,
      save_model = TRUE
    ),
    "Model saved successfully"
  )
  expected_file <- file.path(temp_dir, "Radius_10 km.nc")
  testthat::expect_true(file.exists(expected_file), info = expected_file)
})
