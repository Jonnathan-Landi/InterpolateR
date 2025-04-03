test_that("RFplus works correctly", {
  skip_on_cran()

  data("BD_Obs", package = "InterpolateR")
  data("BD_Coord", package = "InterpolateR")

  # Load the covariates
  Covariates <- list(
    MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
    CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
    DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
  )

  # validation with Rain
  Rain_threshold = list(
    no_rain = c(0, 1),
    light_rain = c(1, 5),
    moderate_rain = c(5, 20),
    heavy_rain = c(20, 40),
    extremely_rain= c(40, Inf)
  )
  ##################################################################################################
  # Testing with validation
  # Apply the RFplus with method = "none"
  model_none = RFplus(BD_Obs, BD_Coord, Covariates, n_round = 1, wet.day = 0.1,
                 ntree = 2000, seed = 123, training = 1, stat_validation = c("M004"),
                 Rain_threshold = Rain_threshold, method = "none",
                 ratio = 15, save_model = FALSE, name_save = NULL)

  # Apply the method RQUANT
  model_Rquant = RFplus(BD_Obs, BD_Coord, Covariates, n_round = 1, wet.day = 0.1,
                 ntree = 2000, seed = 123, training = 1, stat_validation = c("M004"),
                 Rain_threshold = Rain_threshold, method = "RQUANT",
                 ratio = 15, save_model = FALSE, name_save = NULL)

  # Apply the method QUANT
  model_QUANT = RFplus(BD_Obs, BD_Coord, Covariates, n_round = 1, wet.day = 0.1,
                 ntree = 2000, seed = 123, training = 1, stat_validation = c("M004"),
                 Rain_threshold = Rain_threshold, method = "QUANT",
                 ratio = 15, save_model = FALSE, name_save = NULL)

  ##################################################################################################
  # Check that the result is a raster object
  expect_true(inherits(model_none$Ensamble, "SpatRaster"))
  expect_true(inherits(model_Rquant$Ensamble, "SpatRaster"))
  expect_true(inherits(model_QUANT$Ensamble, "SpatRaster"))

  # Check that the number of layers in the raster object is equal to the number of unique dates
  expect_equal(terra::nlyr(model_none$Ensamble), length(unique(BD_Obs$Date)))
  expect_equal(terra::nlyr(model_Rquant$Ensamble), length(unique(BD_Obs$Date)))
  expect_equal(terra::nlyr(model_QUANT$Ensamble), length(unique(BD_Obs$Date)))
  ##################################################################################################
  # Testing without validation
  model_sn = RFplus(BD_Obs, BD_Coord, Covariates, n_round = 1, wet.day = 0.1,
                 ntree = 2000, seed = 123, training = 1, stat_validation = NULL,
                 Rain_threshold = NULL, method = "none",
                 ratio = 15, save_model = FALSE, name_save = NULL)

  # Check that the result is a raster object
  expect_true(inherits(model_sn$Ensamble, "SpatRaster"))
  # Check that the number of layers in the raster object is equal to the number of unique dates
  expect_equal(terra::nlyr(model_sn$Ensamble), length(unique(BD_Obs$Date)))
})

