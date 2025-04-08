test_that("RFmerge works correctly", {
  skip_on_cran()
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
  # Apply the RFmerge
  model_RFmerge = RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
                          seed = 123,  training = 1, stat_validation = c("M004"), Rain_threshold = Rain_threshold,
                          save_model = FALSE, name_save = NULL)

  # Check that the result is a raster object
  expect_true(inherits(model_RFmerge$Ensamble, "SpatRaster"))
  expect_equal(terra::nlyr(model_RFmerge$Ensamble), length(unique(BD_Obs$Date)))

  ##################################################################################################
  # Testing without validation
  model_Sn = RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
                          seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
                          save_model = FALSE, name_save = NULL)

  # Check that the result is a raster object
  expect_true(inherits(model_Sn, "SpatRaster"))
  expect_equal(terra::nlyr(model_Sn), length(unique(BD_Obs$Date)))

  ##################################################################################################
  # Random tesging
  model_Random = RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
                          seed = 123,  training = 0.8, stat_validation = NULL, Rain_threshold = Rain_threshold,
                          save_model = FALSE, name_save = NULL)

  # Check that the result is a raster object
  expect_true(inherits(model_Random$Ensamble, "SpatRaster"))
  expect_equal(terra::nlyr(model_Random$Ensamble), length(unique(BD_Obs$Date)))

  ##############################################################################
  #     Check that the algorithm stops when the input data is not correct.     #
  ##############################################################################
  # Cov
  cov <- data.frame(x = 1:10, y = rnorm(10))  # Aquí creamos un data.frame
  # Intentar ejecutar RFmerge con el data.frame en lugar de una lista
  resultado <- tryCatch({
    RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
            seed = 123, training = 1, stat_validation = c("M004"),
            Rain_threshold = Rain_threshold,
            save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("Parameter detected correctly: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })
})
