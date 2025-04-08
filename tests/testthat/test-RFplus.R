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
  expect_true(inherits(model_sn, "SpatRaster"))
  # Check that the number of layers in the raster object is equal to the number of unique dates
  expect_equal(terra::nlyr(model_sn), length(unique(BD_Obs$Date)))

  ##############################################################################
  #     Check that the algorithm stops when the input data is not correct.     #
  ##############################################################################
  # Verify that covariates is a list
  Cov_df <- data.frame(x = 1:10, y = rnorm(10))
  resultado <- tryCatch({
    RFplus(BD_Obs, BD_Coord, Cov_df, n_round = 1, wet.day = 0.1,
           ntree = 2000, seed = 123, training = 1, stat_validation = NULL,
           Rain_threshold = NULL, method = "none",
           ratio = 15, save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("Parameter detected correctly: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # Verify that the covariates are type SpatRaster
  cov = list(
    MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
    CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
    DEM = data.frame(x = 1:10, y = rnorm(10))
  )

  resultado <- tryCatch({
    RFplus(BD_Obs, BD_Coord, cov, n_round = 1, wet.day = 0.1,
           ntree = 2000, seed = 123, training = 1, stat_validation = NULL,
           Rain_threshold = NULL, method = "none",
           ratio = 15, save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("Parameter detected correctly: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # Verify the BD_Obs is data.table
  BD = as.matrix(BD_Obs)
  cord = as.matrix(BD_Coord)

  resultado <- tryCatch({
    RFplus(BD, BD_Coord, Covariates, n_round = 1, wet.day = 0.1,
           ntree = 2000, seed = 123, training = 1, stat_validation = NULL,
           Rain_threshold = NULL, method = "none",
           ratio = 15, save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("Parameter detected correctly: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # Verify the BD_Obs is data.table
  resultado <- tryCatch({
    RFplus(BD_Obs, cord, Covariates, n_round = 1, wet.day = 0.1,
           ntree = 2000, seed = 123, training = 1, stat_validation = NULL,
           Rain_threshold = NULL, method = "none",
           ratio = 15, save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("Parameter detected correctly: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # Verify that all dates have at least one entry recorded
  BD_na = BD_Obs
  BD_na[1, 2:11] <- NA
  resultado <- tryCatch({
    RFplus(BD_na, BD_Coord, Covariates, n_round = 1, wet.day = 0.1,
           ntree = 2000, seed = 123, training = 1, stat_validation = NULL,
           Rain_threshold = NULL, method = "none",
           ratio = 15, save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("Parameter detected correctly: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # Check that the coordinate names appear in the observed data
  bd_2 = BD_Coord
  bd_2[3,1] <- "aa"
  resultado <- tryCatch({
    RFplus(BD_Obs, bd_2, Covariates, n_round = 1, wet.day = 0.1,
           ntree = 2000, seed = 123, training = 1, stat_validation = NULL,
           Rain_threshold = NULL, method = "none",
           ratio = 15, save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("Parameter detected correctly: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # Check if there is a DEM layer
  COV_2 = list(
    MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
    CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
    DEM = rep(terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR")),2)
  )

  resultado <- tryCatch({
    RFplus(BD_Obs, BD_Coord, COV_2, n_round = 1, wet.day = 0.1,
           ntree = 2000, seed = 123, training = 1, stat_validation = NULL,
           Rain_threshold = NULL, method = "none",
           ratio = 15, save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("Parameter detected correctly: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })

  # Verify the layers of the covariates.
  Cov_3 = list(
    MSWEP = rep(terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),2),
    CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
    DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
  )

  resultado <- tryCatch({
    RFplus(BD_Obs, BD_Coord, Cov_3, n_round = 1, wet.day = 0.1,
           ntree = 2000, seed = 123, training = 1, stat_validation = NULL,
           Rain_threshold = NULL, method = "none",
           ratio = 15, save_model = FALSE, name_save = NULL)
  }, error = function(e) {
    message("Parameter detected correctly: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("warning: ", w$message)
    return(NULL)
  })


})

