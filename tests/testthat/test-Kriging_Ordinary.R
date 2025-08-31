# Data used for testing kriging ordinary
# Date of test creation: 2025-08-30
# Test update date:  2025-08-30

# Data input
data("BD_Obs", package = "InterpolateR")
data("BD_Coord", package = "InterpolateR")

# load mask
shapefile <- terra::vect(system.file(
  "extdata/study_area.shp",
  package = "InterpolateR"
))

Rain_threshold = list(
  no_rain = c(0, 1),
  light_rain = c(1, 5),
  moderate_rain = c(5, 20),
  heavy_rain = c(20, 40),
  extremely_rain = c(40, Inf)
)

# 1. Testing with validation ---------------------------------------------------
testthat::test_that("Kriging_Ordinary returns SpatRaster without validation.", {
  testthat::skip_on_cran()
  out <- Kriging_Ordinary(
    BD_Obs,
    BD_Coord,
    shapefile,
    grid_resolution = 5,
    variogram_model = "exponential",
    max_dist = 15,
    n_lags = 15,
    min_stations = 2,
    n_round = NULL,
    training = 1,
  )
  testthat::expect_true(inherits(out, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(out), length(unique(BD_Obs$Date)))
})
