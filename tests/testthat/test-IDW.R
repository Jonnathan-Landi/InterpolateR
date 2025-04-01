test_that("IDW interpolation works correctly", {
  data("BD_Obs", package = "InterpolateR")
  data("BD_Coord", package = "InterpolateR")

  # Load the study area where the interpolation will be performed.
  shapefile <- terra::vect(system.file("extdata/study_area.shp", package = "InterpolateR"))

  # Performing interpolation using the IDW method
  Interpolated_data <- IDW(BD_Obs, BD_Coord, shapefile, grid_resolution = 5, p = 2,
                           n_round = 1, training = 0.8, Rain_threshold = NULL,
                           stat_validation = NULL, save_model = FALSE, name_save = NULL)

  # Check that the result is a raster object
  expect_true(inherits(Interpolated_data$Ensamble, "SpatRaster"))
  expect_equal(terra::nlyr(Interpolated_data$Ensamble), length(unique(BD_Obs$Date)))
})
