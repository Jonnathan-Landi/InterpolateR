# Data used for testing RFplus
# Date of test creation: 2025-04-29
# Test update date: 2025-04-29

# Data input
data("BD_Obs", package = "InterpolateR")
data("BD_Coord", package = "InterpolateR")

# Load the Covariatesariates
Covariates <- list(
  MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
  CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
  DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
)

Rain_threshold = list(
  no_rain = c(0, 1),
  light_rain = c(1, 5),
  moderate_rain = c(5, 20),
  heavy_rain = c(20, 40),
  extremely_rain= c(40, Inf)
)

testthat::skip_on_cran()

# 1. Testing with validation ---------------------------------------------------
testthat::test_that("RFplus returns SpatRaster without validation.", {
  testthat::skip_on_cran()
  out <- RFplus(BD_Obs, BD_Coord, Covariates, n_round = 1, ntree = 2000,
    seed = 123, method = "none", ratio = 10, training = 1, stat_validation = NULL,
    Rain_threshold = NULL, save_model = FALSE, name_save = NULL)

  testthat::expect_true(inherits(out, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out), length(unique(BD_Obs$Date)))
})

# 2. Testing with validation (random validation) --------------------------------
testthat::test_that("RFplus returns SpatRaster with random validation.", {
  testthat::skip_on_cran()
  out <- RFplus(BD_Obs, BD_Coord, Covariates,  n_round = 1, ntree = 2000, wet.day = 0.1,
    seed = 123, training = 0.8, method = "none", ratio = 10, stat_validation = NULL,
    Rain_threshold = Rain_threshold, save_model = FALSE, name_save = NULL)

  testthat::expect_true(inherits(out$Ensamble, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out$Ensamble), length(unique(BD_Obs$Date)))
  testthat::expect_true(inherits(out$Validation$gof, "data.table"))
  testthat::expect_true(inherits(out$Validation$categorical_metrics, "data.table"))
})

# 3. Testing with validation (manual validation) --------------------------------
testthat::test_that("RFplus returns SpatRaster with manual validation.", {
  testthat::skip_on_cran()
  out <- RFplus(BD_Obs, BD_Coord, Covariates,  n_round = 1, ntree = 2000,
    seed = 123, training = 1, method = "none", ratio = 10, stat_validation = "M004",
    Rain_threshold = Rain_threshold, save_model = FALSE, name_save = NULL)

  testthat::expect_true(inherits(out$Ensamble, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out$Ensamble), length(unique(BD_Obs$Date)))
  testthat::expect_true(inherits(out$Validation$gof, "data.table"))
  testthat::expect_true(inherits(out$Validation$categorical_metrics, "data.table"))
})

# 4. Testing  RFplus returns SpatRaster with random validation without ---------
testthat::test_that("RFplus returns SpatRaster with manual validation.", {
  testthat::skip_on_cran()
  out <- RFplus(BD_Obs, BD_Coord, Covariates,  n_round = 1, ntree = 2000,
    seed = 123, training = 1, method = "none", ratio = 10, stat_validation = "M004",
    Rain_threshold = NULL, save_model = FALSE, name_save = NULL)

  testthat::expect_true(inherits(out$Ensamble, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out$Ensamble), length(unique(BD_Obs$Date)))
  testthat::expect_true(inherits(out$Validation, "data.table"))
})

# 5. Testing  RFplus with method QUANT ---------
testthat::test_that("RFplus uses the QUANT Method.", {
  testthat::skip_on_cran()
  out <- RFplus(BD_Obs, BD_Coord, Covariates,  n_round = 1, ntree = 2000,
    seed = 123, training = 1, method = "QUANT", ratio = 10, stat_validation = NULL,
    Rain_threshold = NULL, save_model = FALSE, name_save = NULL)

  testthat::expect_true(inherits(out, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out), length(unique(BD_Obs$Date)))
})

# 6. Testing  RFplus with method RQUANT ---------
testthat::test_that("RFplus uses the QUANT Method.", {
  testthat::skip_on_cran()
  out <- RFplus(BD_Obs, BD_Coord, Covariates,  n_round = 1, ntree = 2000,
    seed = 123, training = 1, method = "RQUANT", ratio = 10, stat_validation = NULL,
    Rain_threshold = NULL, save_model = FALSE, name_save = NULL)

  testthat::expect_true(inherits(out, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out), length(unique(BD_Obs$Date)))
})
##############################################################################
#     Check that the algorithm stops when the input data is not correct.     #
##############################################################################
# 5. Testing error with The Covariates must be a list.
testthat::test_that("The Covariates must be a list.", {
  testthat::skip_on_cran()
  bad_Covariates <- data.frame(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    RFplus(BD_Obs, BD_Coord, bad_Covariates,  n_round = 1, ntree = 2000,
      seed = 123, training = 1, method = "none", ratio = 10, stat_validation = NULL,
      Rain_threshold = NULL, save_model = FALSE, name_save = NULL),
    regexp = "Covariates must be a list\\.$"
    )
  })

# 6. The Covariates must be of type SpatRaster.
testthat::test_that("The Covariates must be of type SpatRaster.", {
  testthat::skip_on_cran()
  bad_Covariates <- Covariates
  bad_Covariates$Test <- data.frame(x = 1:10, y = rnorm(10))
  testthat::expect_error(
    RFplus(BD_Obs, BD_Coord, bad_Covariates, n_round = 1, ntree = 2000,
      seed = 123, training = 1, method = "none", ratio = 10, stat_validation = NULL,
      Rain_threshold = NULL, save_model = FALSE, name_save = NULL),
    regexp = "The covariates must be of type SpatRaster\\.$"
    )
})

# 7. The extension of the Covariates are different (all extensions should be similar)
testthat::test_that("The extension of the Covariates are different.", {
  testthat::skip_on_cran()
  bad_Covariates <- Covariates
  terra::ext(bad_Covariates$MSWEP) = terra::ext(bad_Covariates$MSWEP) * 2
  testthat::expect_error(
    RFplus(BD_Obs, BD_Coord, bad_Covariates,  n_round = 1, ntree = 2000,
      seed = 123, training = 1, method = "none", ratio = 10, stat_validation = NULL, Rain_threshold = NULL,
      save_model = FALSE, name_save = NULL),
    regexp = "^The extension of the covariates are different \\(all extensions should be similar\\)\\.$"
    )
})

# 8. "The crs of the Covariates are different (all crs should be similar).
testthat::test_that("The crs of the Covariates are different.", {
  testthat::skip_on_cran()
  bad_Covariates <- Covariates
  terra::crs(bad_Covariates$MSWEP) = "EPSG:4326"
  testthat::expect_error(
    RFplus(BD_Obs, BD_Coord, bad_Covariates,  n_round = 1, ntree = 2000,
      seed = 123,  training = 1, method = "none", ratio = 10, stat_validation = NULL,
      Rain_threshold = NULL, save_model = FALSE, name_save = NULL),
    regexp = "^The crs of the covariates are different \\(all crs should be similar\\)\\.$"
    )
})

# 9. "BD_Obs must be a 'data.table' or a 'data.frame'."
testthat::test_that("BD_Obs must be a 'data.table' or a 'data.frame.", {
  testthat::skip_on_cran()
  bad_BDObs = as.matrix(BD_Obs)
  testthat::expect_error(
    RFplus(bad_BDObs, BD_Coord, Covariates,  n_round = 1, ntree = 2000,
      seed = 123, training = 1, method = "none", ratio = 10, stat_validation = NULL,
      Rain_threshold = NULL, save_model = FALSE, name_save = NULL),
    regexp = "BD_Obs must be a 'data.table' or a 'data.frame'.$"
    )
})

# 10. BD_Coord must be a 'data.table' or a 'data.frame'."
testthat::test_that("BD_Coord must be a 'data.table' or a 'data.frame'.", {
  testthat::skip_on_cran()
  bad_BD_Coord = as.matrix(BD_Coord)
  testthat::expect_error(
    RFplus(BD_Obs, bad_BD_Coord, Covariates,  n_round = 1, ntree = 2000,
      seed = 123,  training = 1, method = "none", ratio = 10, stat_validation = NULL,
      Rain_threshold = NULL, save_model = FALSE, name_save = NULL),
    regexp = "BD_Coord must be a 'data.table' or a 'data.frame'.$"
    )
})

# 11. "The names of the coordinates do not appear in the observed data."
testthat::test_that("names of the coordinates do not appear in the observed data.", {
  testthat::skip_on_cran()
  bad_Cords = BD_Coord
  bad_Cords[3,1] = "yy"
  testthat::expect_error(
    RFplus(BD_Obs, bad_Cords, Covariates, n_round = 1, ntree = 2000,
      seed = 123,  training = 1, method = "none", ratio = 10, stat_validation = NULL,
      Rain_threshold = NULL, save_model = FALSE, name_save = NULL),
    regexp = "The names of the coordinates do not appear in the observed data.$"
    )
})

# 12. A single layer covariate was not found. Possibly the DEM was not entered..
testthat::test_that("Single layer covariate was not found.", {
  testthat::skip_on_cran()
  bad_Covs = Covariates
  bad_Covs$DEM = c(bad_Covs$DEM, bad_Covs$DEM)
  testthat::expect_error(
    RFplus(BD_Obs, BD_Coord, bad_Covs, n_round = 1, ntree = 2000, method = "none", ratio = 10,
      seed = 123, training = 1, stat_validation = NULL, Rain_threshold = NULL,
      save_model = FALSE, name_save = NULL),
    regexp = "A single layer covariate was not found. Possibly the DEM was not entered.$"
    )
})

# 13. Verify that all dates have at least one entry recorded
testthat::test_that("Error when BD_Obs has rows with all NA values (except Date)", {
  testthat::skip_on_cran()
  # Simular BD_Obs con una fila completamente NA (menos la fecha)
  BD_Obs <- data.table::data.table(
    Date = as.Date(c("2020-01-01", "2020-01-02")),
    Station1 = c(1.2, NA),
    Station2 = c(3.4, NA)
  )

  BD_Coord <- data.table::data.table(ID = c("Station1", "Station2"), x = c(1, 2), y = c(3, 4))

  # Crear una Covariatesariables falsas
  Covariates <- list(dummy = terra::rast(nrows = 2, ncols = 2, vals = runif(4)))
  testthat::expect_error(
    RFplus(BD_Obs, BD_Coord, Covariates,
            n_round = 1, ntree = 10, seed = 123, training = 1,
            stat_validation = NULL, Rain_threshold = NULL, method = "none", ratio = 10,
            save_model = FALSE, name_save = NULL),
    regexp = "No data was found for the dates: 2020-01-02"
  )
})

# 14. "Save the model must be a logical value."
testthat::test_that("RFplus saves model when save_model = TRUE", {
  testthat::skip_on_cran()
  temp_dir <- tempdir()
  withr::local_dir(temp_dir)
  custom_name <- "test_model_RFplus"
  expect_message(
    out <- RFplus(BD_Obs, BD_Coord, Covariates, method = "none", ratio = 10,
      n_round = 1, ntree = 10, seed = 123, training = 1,
      stat_validation = NULL, Rain_threshold = NULL,
      save_model = T, name_save = custom_name
  ),
  "Model saved successfully"
)
  expected_file <- file.path(temp_dir, paste0(custom_name, ".nc"))
  testthat::expect_true(file.exists(expected_file), info = expected_file)
})

# 15. "Save the model must be a logical value." (default name) " ---------------
testthat::test_that("RFplus saves model when save_model = TRUE (default name)", {
  testthat::skip_on_cran()
  temp_dir <- tempdir()
  withr::local_dir(temp_dir)
  expect_message(
    out <- RFplus(BD_Obs, BD_Coord, Covariates, method = "none", ratio = 10,
      n_round = 1, ntree = 10, seed = 123, training = 1,
      stat_validation = NULL, Rain_threshold = NULL,
      save_model = T, name_save = NULL
  ),
    "Model saved successfully"
  )
  expected_file <- file.path(temp_dir, "Model_RFplus.nc")
  testthat::expect_true(file.exists(expected_file), info = expected_file)
})

# 16. Test adicional para modulo de conexion
testthat::test_that("The training parameter must be between 0 and 1.", {
  testthat::skip_on_cran()
  testthat::expect_error(
    RFplus(BD_Obs, BD_Coord, Covariates, n_round = 1, ntree = 2000,
           seed = 123, training = 2, stat_validation = NULL,
           Rain_threshold = NULL, save_model = FALSE, name_save = NULL),
    regexp = "The training parameter must be between 0 and 1\\.$"
  )
})
