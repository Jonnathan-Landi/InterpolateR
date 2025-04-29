# Data used for testing RFmerge
# Date of test creation: 2025-04-29
# Test update date: 2025-04-29

# Data input
data("BD_Obs", package = "InterpolateR")
data("BD_Coord", package = "InterpolateR")

# Load the covariates
cov <- list(
  MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
  CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
  DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
)

# load mask
shapefile <- terra::vect(system.file("extdata/study_area.shp", package = "InterpolateR"))

Rain_threshold = list(
  no_rain = c(0, 1),
  light_rain = c(1, 5),
  moderate_rain = c(5, 20),
  heavy_rain = c(20, 40),
  extremely_rain= c(40, Inf)
)

testthat::skip_on_cran()

# 1. Testing with validation ---------------------------------------------------
testthat::test_that("RFmerge returns SpatRaster without validation.", {
  out <- RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
    seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
    save_model = FALSE, name_save = NULL)
  
  testthat::expect_true(inherits(out, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out), length(unique(BD_Obs$Date)))
})

# 2. Testing with validation (random validation) --------------------------------
testthat::test_that("RFmerge returns SpatRaster with random validation.", {
  out <- RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
    seed = 123,  training = 0.8, stat_validation = NULL, Rain_threshold = Rain_threshold,
    save_model = FALSE, name_save = NULL)
  
  testthat::expect_true(inherits(out$Ensamble, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out$Ensamble), length(unique(BD_Obs$Date)))
  testthat::expect_true(inherits(out$Validation$gof, "data.table"))
  testthat::expect_true(inherits(out$Validation$categorical_metrics, "data.table"))
})

# 3. Testing with validation (manual validation) --------------------------------
testthat::test_that("RFmerge returns SpatRaster with manual validation.", {
  out <- RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
    seed = 123,  training = 1, stat_validation = "M004", Rain_threshold = Rain_threshold,
    save_model = FALSE, name_save = NULL)
  
  testthat::expect_true(inherits(out$Ensamble, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out$Ensamble), length(unique(BD_Obs$Date)))
  testthat::expect_true(inherits(out$Validation$gof, "data.table"))
  testthat::expect_true(inherits(out$Validation$categorical_metrics, "data.table"))
})

# 4. Testing  RFmerge returns SpatRaster with random validation without ---------
testthat::test_that("RFmerge returns SpatRaster with manual validation.", {
  out <- RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
    seed = 123,  training = 1, stat_validation = "M004", Rain_threshold = NULL,
    save_model = FALSE, name_save = NULL)
  
  testthat::expect_true(inherits(out$Ensamble, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out$Ensamble), length(unique(BD_Obs$Date)))
  testthat::expect_true(inherits(out$Validation, "data.table"))
})

##############################################################################
#     Check that the algorithm stops when the input data is not correct.     #
##############################################################################
# 5. Testing error with The cov must be a list.
testthat::test_that("The cov must be a list.", {
  bad_cov <- data.frame(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    RFmerge(BD_Obs, BD_Coord, bad_cov, mask = shapefile, n_round = 1, ntree = 2000,
      seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
      save_model = FALSE, name_save = NULL),
    regexp = "cov must be a list\\.$"
    )
  })

# 6. The cov must be of type SpatRaster.
testthat::test_that("The cov must be of type SpatRaster.", {
  bad_cov <- cov
  bad_cov$Test <- data.frame(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    RFmerge(BD_Obs, BD_Coord, bad_cov, mask = shapefile, n_round = 1, ntree = 2000,
      seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
      save_model = FALSE, name_save = NULL),
    regexp = "The cov must be of type SpatRaster\\.$"
    )
})

# 7. The extension of the cov are different (all extensions should be similar)
testthat::test_that("The extension of the cov are different.", {
  bad_cov <- cov
  terra::ext(bad_cov$MSWEP) = terra::ext(bad_cov$MSWEP) * 2
  testthat::expect_error(
    RFmerge(BD_Obs, BD_Coord, bad_cov, mask = shapefile, n_round = 1, ntree = 2000,
      seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
      save_model = FALSE, name_save = NULL),
    regexp = "^The extension of the cov are different \\(all extensions should be similar\\)\\.$"
    )
})

# 8. "The crs of the cov are different (all crs should be similar).
testthat::test_that("The crs of the cov are different.", {
  bad_cov <- cov
  terra::crs(bad_cov$MSWEP) = "EPSG:4326"
  testthat::expect_error(
    RFmerge(BD_Obs, BD_Coord, bad_cov, mask = shapefile, n_round = 1, ntree = 2000,
      seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
      save_model = FALSE, name_save = NULL),
    regexp = "^The crs of the cov are different \\(all crs should be similar\\)\\.$"
    )
})

# 9. "BD_Obs must be a 'data.table' or a 'data.frame'."
testthat::test_that("BD_Obs must be a 'data.table' or a 'data.frame'.", {
  bad_BDObs = as.matrix(BD_Obs)
  testthat::expect_error(
    RFmerge(bad_BDObs, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
      seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
      save_model = FALSE, name_save = NULL),
    regexp = "BD_Obs must be a 'data.table' or a 'data.frame'.$"
    )
})

# 10. BD_Coord must be a 'data.table' or a 'data.frame'."
testthat::test_that("BD_Coord must be a 'data.table' or a 'data.frame'.", {
  bad_BD_Coord = as.matrix(BD_Coord)
  testthat::expect_error(
    RFmerge(BD_Obs, bad_BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
      seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
      save_model = FALSE, name_save = NULL),
    regexp = "BD_Coord must be a 'data.table' or a 'data.frame'.$"
    )
})

# 11. "The names of the coordinates do not appear in the observed data."
testthat::test_that("names of the coordinates do not appear in the observed data.", {
  bad_Cords = BD_Coord
  bad_Cords[3,1] = "yy"
  testthat::expect_error(
    RFmerge(BD_Obs, bad_Cords, cov, mask = shapefile, n_round = 1, ntree = 2000,
      seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
      save_model = FALSE, name_save = NULL),
    regexp = "The names of the coordinates do not appear in the observed data.$"
    )
})

# 12. mask must be a 'SpatVector' object.
testthat::test_that("mask must be a 'SpatVector'.", {
  bad_shp = BD_Obs
  testthat::expect_error(
    RFmerge(BD_Obs, BD_Coord, cov, mask = bad_shp, n_round = 1, ntree = 2000,
      seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
      save_model = FALSE, name_save = NULL),
    regexp = "mask must be a 'SpatVector' object.$"
    )
})

# 13. Verify that all dates have at least one entry recorded
testthat::test_that("Error when BD_Obs has rows with all NA values (except Date)", {
  # Simular BD_Obs con una fila completamente NA (menos la fecha)
  BD_Obs <- data.table::data.table(
    Date = as.Date(c("2020-01-01", "2020-01-02")),
    Station1 = c(1.2, NA),
    Station2 = c(3.4, NA)
  )

  BD_Coord <- data.table::data.table(ID = c("Station1", "Station2"), x = c(1, 2), y = c(3, 4))

  # Crear una covariables falsas
  cov <- list(dummy = terra::rast(nrows = 2, ncols = 2, vals = runif(4)))
  shapefile <- terra::vect(matrix(c(1,1, 2,2, 2,1), ncol = 2, byrow = TRUE), type = "polygons")

  testthat::expect_error(
    RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile,
            n_round = 1, ntree = 10, seed = 123, training = 1,
            stat_validation = NULL, Rain_threshold = NULL,
            save_model = FALSE, name_save = NULL),
    regexp = "No data was found for the dates: 2020-01-02"
  )
})

# 14. "Save the model must be a logical value."
testthat::test_that("RFmerge saves model when save_model = TRUE", {
  temp_dir <- tempdir()
  withr::local_dir(temp_dir)
  custom_name <- "test_model_RFmerge"
  expect_message(
    out <- RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile,
      n_round = 1, ntree = 10, seed = 123, training = 1,
      stat_validation = NULL, Rain_threshold = NULL,
      save_model = T, name_save = custom_name
  ),
  "Model saved successfully"
)
  expected_file <- file.path(temp_dir, paste0(custom_name, ".nc"))
  testthat::expect_true(file.exists(expected_file), info = expected_file)
  testthat::expect_true(inherits(terra::rast(expected_file), "SpatRaster"))
})

# 15. "Save the model must be a logical value." (default name) " ------------------
testthat::test_that("RFmerge saves model when save_model = TRUE (default name)", {
  temp_dir <- tempdir()
  withr::local_dir(temp_dir)
  expect_message(
    out <- RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile,
      n_round = 1, ntree = 10, seed = 123, training = 1,
      stat_validation = NULL, Rain_threshold = NULL,
      save_model = T, name_save = NULL
  ),
    "Model saved successfully"
  )
  expected_file <- file.path(temp_dir, "Model_RFmerge.nc")
  testthat::expect_true(file.exists(expected_file), info = expected_file)
  testthat::expect_true(inherits(terra::rast(expected_file), "SpatRaster"))
})

# 15. A single layer covariate was not found. Possibly the DEM was not entered..
testthat::test_that("Single layer covariate was not found.", {
  bad_Covs = cov
  bad_Covs$DEM = c(bad_Covs$DEM, bad_Covs$DEM)
  testthat::expect_error(
    RFmerge(BD_Obs, BD_Coord, bad_Covs, mask = shapefile,
      n_round = 1, ntree = 10, seed = 123, training = 1,
      stat_validation = NULL, Rain_threshold = NULL,
      save_model = T, name_save = NULL),
    regexp = "A single layer covariate was not found. Possibly the DEM was not entered.$"
    )
})