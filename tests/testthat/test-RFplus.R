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

  # Apply the RFplus
  model = RFplus(BD_Obs, BD_Coord, Covariates, n_round = 1, wet.day = 0.1,
                 ntree = 2000, seed = 123, training = 1, stat_validation = c("M006"),
                 Rain_threshold = Rain_threshold, method = "none",
                 ratio = 15, save_model = FALSE, name_save = NULL)

  # Check that the result is a raster object
  expect_true(inherits(model$Ensamble, "SpatRaster"))
  expect_equal(terra::nlyr(model$Ensamble), length(unique(BD_Obs$Date)))
})

