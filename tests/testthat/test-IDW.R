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

  # Performing interpolation using the IDW method
  Interpolated_data <- IDW(BD_Obs, BD_Coord, shapefile, grid_resolution = 5, p = 2,
                           n_round = 1, training = 0.8, Rain_threshold = Rain_threshold,
                           stat_validation = NULL, save_model = FALSE, name_save = NULL)

  # Check that the result is a raster object
  expect_true(inherits(Interpolated_data$Ensamble, "SpatRaster"))
  expect_equal(terra::nlyr(Interpolated_data$Ensamble), length(unique(BD_Obs$Date)))
})
