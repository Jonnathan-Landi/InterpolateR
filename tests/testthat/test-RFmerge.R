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
  # Verify that cov is a list
  cov_1 <- data.frame(x = 1:10, y = rnorm(10))
  resultado <- tryCatch({
    RFmerge(BD_Obs, BD_Coord, cov_1, mask = shapefile, n_round = 1, ntree = 2000,
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

  # Verify that the cov are type SpatRaster
  cov_2 = list(
    MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
    CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
    DEM = data.frame(x = 1:10, y = rnorm(10))
  )

  resultado <- tryCatch({
    RFmerge(BD_Obs, BD_Coord, cov_2, mask = shapefile, n_round = 1, ntree = 2000,
            seed = 123, training = 1, stat_validation = c("M004"),
            Rain_threshold = Rain_threshold,
            save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("error message: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # Verify the extent of cov
  Cov_3 = list(
    MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
    CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
    DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
  )
  terra::ext(Cov_3$MSWEP) <- terra::ext(Cov_3$MSWEP) * 2
  resultado <- tryCatch({
    RFmerge(BD_Obs, BD_Coord, Cov_3, mask = shapefile, n_round = 1, ntree = 2000,
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

  # Verify the crc of cov
  Cov_4 = list(
    MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
    CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
    DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
  )
  terra::crs(Cov_4$MSWEP) <- "EPSG:4326"
  resultado <- tryCatch({
    RFmerge(BD_Obs, BD_Coord, Cov_4, mask = shapefile_e, n_round = 1, ntree = 2000,
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

  # BD_Coord can be a data.table or a data.frame
  BD_Obs_m = as.matrix(BD_Obs)
  resultado <- tryCatch({
    RFmerge(BD_Obs_m, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
            seed = 123, training = 1, stat_validation = c("M004"),
            Rain_threshold = Rain_threshold,
            save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("error message: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # BD_Coord can be a data.table or a data.frame
  BD_Coord_m <- as.matrix(BD_Coord)
  resultado <- tryCatch({
    RFmerge(BD_Obs, BD_Coord_m, cov, mask = shapefile, n_round = 1, ntree = 2000,
            seed = 123, training = 1, stat_validation = c("M004"),
            Rain_threshold = Rain_threshold,
            save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("error message: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # # Check that the coordinate names appear in the observed data
  bd_2 = BD_Coord
  bd_2[3,1] <- "aa"
  resultado <- tryCatch({
    RFmerge(BD_Obs, bd_2, cov, mask = shapefile, n_round = 1, ntree = 2000,
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


  # Check if mask is a SpatVector object
  shapefile_e <- data.frame(x = 1:10, y = rnorm(10))
  resultado <- tryCatch({
    RFmerge(BD_Obs, BD_Coord, cov, mask = shapefile_e, n_round = 1, ntree = 2000,
            seed = 123, training = 1, stat_validation = c("M004"),
            Rain_threshold = Rain_threshold,
            save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("error message: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # Verify that all dates have at least one entry recorded
  BD_na = BD_Obs
  BD_na[1, 2:11] <- NA
  resultado <- tryCatch({
    RFmerge(BD_na, BD_Coord, cov, mask = shapefile, n_round = 1, ntree = 2000,
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

  # Check if there is a DEM layer
  COV_5 = list(
    MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
    CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
    DEM = rep(terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR")),2)
  )

  resultado <- tryCatch({
    RFmerge(BD_Obs, BD_Coord, COV_5, mask = shapefile, n_round = 1, ntree = 2000,
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
  # # Verify the layers of the cov.
  # Cov_6 = list(
  #   MSWEP = rep(terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),2),
  #   CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
  #   DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
  # )
  #
  # resultado <- tryCatch({
  #   RFmerge(BD_Obs, BD_Coord, Cov_6, mask = shapefile, n_round = 1, ntree = 2000,
  #           seed = 123, training = 1, stat_validation = c("M004"),
  #           Rain_threshold = Rain_threshold,
  #           save_model = FALSE, name_save = NULL)
  # }, error = function(e) {
  #   message("Parameter detected correctly: ", e$message)
  #   return(NULL)
  # }, warning = function(w) {
  #   message("warning: ", w$message)
  #   return(NULL)
  # })

}) # end
