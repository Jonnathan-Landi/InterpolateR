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

  # Apply the RFmerge
  model_RFmerge = RFmerge(BD_Obs, BD_Coord, cov, mask = NULL, n_round = 1, ntree = 2000,
                          seed = 123,  training = 1, stat_validation = c("M006"), Rain_threshold = NULL,
                          save_model = FALSE, name_save = NULL)


  # Check that the result is a raster object
  expect_true(inherits(model_RFmerge$Ensamble, "SpatRaster"))
  expect_equal(terra::nlyr(model_RFmerge$Ensamble), length(unique(BD_Obs$Date)))
})
