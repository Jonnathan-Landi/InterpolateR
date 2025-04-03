test_that("Cressman Objective Analysis Method works correctly", {
  skip_on_cran()

  data("BD_Obs", package = "InterpolateR")
  data("BD_Coord", package = "InterpolateR")

  # Load the study area where the interpolation will be performed.
  shapefile = terra::vect(system.file("extdata/study_area.shp", package = "InterpolateR"))

  # validation with Rain
  Rain_threshold = list(
    no_rain = c(0, 1),
    light_rain = c(1, 5),
    moderate_rain = c(5, 20),
    heavy_rain = c(20, 40),
    extremely_rain= c(40, Inf)
  )

  # Performing interpolation using the Cressman method
  Interpolated_Cressman = Cressman(
    BD_Obs, BD_Coord, shapefile, grid_resolution = 5,
    search_radius = c(20, 10), training = 1,
    stat_validation = "M001", Rain_threshold = Rain_threshold,
    save_model = FALSE
  )

  # Results
  Radius_20 = Interpolated_Cressman$Ensamble[[1]] # Interpolated data with a 20 km radius
  Radius_10 = Interpolated_Cressman$Ensamble[[2]] # Interpolated data with a 10 km radius

  # Validation statistics
  Validation_results_20 = Interpolated_Cressman$Validation[[1]] # Validation results with a 20 km radius
  Validation_results_10 = Interpolated_Cressman$Validation[[2]] # Validation results with a 10 km radius

  # Check that the result is a raster object
  expect_true(inherits(Radius_20, "SpatRaster"))
  expect_true(inherits(Radius_10, "SpatRaster"))
  expect_equal(terra::nlyr(Radius_20), length(unique(BD_Obs$Date)))
  expect_equal(terra::nlyr(Radius_10), length(unique(BD_Obs$Date)))
})
