test_that("IDW interpolation works correctly", {
  data("BD_Obs", package = "InterpolateR")
  data("BD_Coord", package = "InterpolateR")

  # Load the study area where the interpolation will be performed.
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
  # Performing interpolation using the IDW method
  Interpolated_data <- IDW(BD_Obs, BD_Coord, shapefile, grid_resolution = 5, p = 2,
                           n_round = 1, training = 1, Rain_threshold = Rain_threshold,
                           stat_validation = "M004", save_model = FALSE, name_save = NULL)

  # Check that the result is a raster object
  expect_true(inherits(Interpolated_data$Ensamble, "SpatRaster"))
  expect_equal(terra::nlyr(Interpolated_data$Ensamble), length(unique(BD_Obs$Date)))
  ##################################################################################################
  # Testing without validation
  Interpolated_SN = IDW(BD_Obs, BD_Coord, shapefile, grid_resolution = 5, p = 2,
                           n_round = 1, training = 1, Rain_threshold = NULL,
                           stat_validation = NULL, save_model = FALSE, name_save = NULL)

  # Check that the result is a raster object
  expect_true(inherits(Interpolated_SN, "SpatRaster"))
  expect_equal(terra::nlyr(Interpolated_SN), length(unique(BD_Obs$Date)))
  ##################################################################################################
  # Testing without Rain
  Interpolated_SnRain <- IDW(BD_Obs, BD_Coord, shapefile, grid_resolution = 5, p = 2,
                           n_round = 1, training = 1, Rain_threshold = NULL,
                           stat_validation = "M004", save_model = FALSE, name_save = NULL)
  # Check that the result is a raster object
  expect_true(inherits(Interpolated_SnRain$Ensamble, "SpatRaster"))
  expect_equal(terra::nlyr(Interpolated_SnRain$Ensamble), length(unique(BD_Obs$Date)))
  ##################################################################################################
})
