#' Barnes objective analysis
#' @export
Barnes <- function(BD_Obs, BD_Coord, shapefile, grid_resolution, Kappa = 0.5, Gamma = 1, iterations = 2,
                   training = 1, stat_validation = NULL, Rain_threshold = NULL, save_model = FALSE) {
  ##############################################################################
  #                               Check input data                             #
  ##############################################################################
  # Shapefile must be a 'spatVector' object and coordinate reference system (CRS) must be defined
  if (!inherits(shapefile, "SpatVector")) stop("shapefile must be a 'SpatVector' object.")

  # BD_Obs can be a data.table or a data.frame
  if (!inherits(BD_Obs, c("data.table", "data.frame"))) stop("BD_Obs must be a 'data.table' or a 'data.frame'.")
  names(BD_Obs)[1] <- "Date"

  # BD_Coord can be a data.table or a data.frame
  if (!inherits(BD_Coord, c("data.table", "data.frame"))) stop("BD_Coord must be a 'data.table' or a 'data.frame'.")

  # Check that the coordinate names appear in the observed data
  if (!all(BD_Coord$Cod %chin% setdiff(names(BD_Obs), "Date"))) stop("The names of the coordinates do not appear in the observed data.")
  ##############################################################################
  #                          Verify if validation is to be done                #
  ##############################################################################
  names_col <- setdiff(names(BD_Obs), "Date")
  Ids <- data.table::data.table(Cod = names_col, ID = 1:length(names_col))
  if (training != 1 | !is.null(stat_validation)) {
    data_val = .select_data(BD_Obs, BD_Coord, training = training, seed = 123,
                            stat_validation = stat_validation)
    train_data <- data_val$train_data
    train_cords <- data_val$train_cords
  } else {
    train_data <- BD_Obs
    train_cords <- BD_Coord
  }
  ##############################################################################
  #                           Zone of Interpolation                            #
  ##############################################################################
  grid_resolution <- grid_resolution * 1000
  coord.ref <- terra::crs(shapefile)
  spl_layer <- rast(
    terra::ext(shapefile),
    resolution = grid_resolution,
    crs = coord.ref)
  terra::values(spl_layer) <- 0
  ##############################################################################
  #                           Data training                                    #
  ##############################################################################
  Barnes_data <- melt(
    train_data,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  )
  Barnes_data <- Ids[Barnes_data, on = "Cod"]
  Dates_extracted <- unique(Barnes_data[, Date])
  Points_Train <- merge(Barnes_data, train_cords, by = "Cod")
  setDT(Points_Train)
  Points_Train <- unique(Points_Train, by = "Cod")[, .(ID, Cod, X, Y, Z)]
  setorder(Points_Train, ID)
  Points_VectTrain <- terra::vect(Points_Train, geom = c("X", "Y"), crs = coord.ref)
  ##############################################################################
  #                          Barnes Algorithm                                  #
  ##############################################################################
  data_Barnes <- data.table::as.data.table(as.data.frame(spl_layer, xy = TRUE))
  setnames(data_Barnes, c("X", "Y", "value"))
  coords <- as.matrix(data_Barnes[, .(X, Y)])
  distancias <- as.data.table(terra::distance(terra::vect(coords, crs = crs(spl_layer)), Points_VectTrain))
  setnames(distancias, Points_VectTrain$Cod)
  data_Barnes <- cbind(data_Barnes[, .(X, Y)], distancias)

  # weight calculation
  estaciones <- Points_VectTrain$Cod
  w <- exp(-data_Barnes[, ..estaciones]^2 / Kappa)  # Gaussian weights

  barnes_logic_opt <- function(data_obs) {
    obs_values = data_obs$var[match(Points_VectTrain$Cod, data_obs$Cod)]
    obs_matrix = matrix(obs_values, nrow = nrow(data_Barnes), ncol = length(Points_VectTrain$Cod), byrow = TRUE)

    station_cols = setdiff(names(data_Barnes), c("X", "Y"))
    w = exp(-data_Barnes[, ..station_cols]^2 / gamma)  # Ponderaciones gaussianas
    numer = rowSums(w * obs_matrix, na.rm = TRUE)
    denom = rowSums(w, na.rm = TRUE)
    new_z = numer / denom

    #
    z_interp = Reduce(function(z_prev, z_new) z_prev + alpha * (z_new - z_prev),
                      x = rep(list(new_z), max_iter),
                      init = new_z,
                      accumulate = TRUE)[[max_iter]]

    # Convertir a formato raster
    rast(
      data_Barnes[, .(X, Y, value = ifelse(is.na(z_interp) | denom == 0, NA, z_interp))],
      crs = crs(spl_layer),
      type = "xyz"
    )
  }

  call_barnes_opt <- function(day) {
    data_obs <- Barnes_data[Date == as.Date(day), .(Cod, var)]
    if (nrow(data_obs) == 0) return(NULL)
    barnes_logic_opt(data_obs)
  }

  pbapply::pboptions(type = "timer", use_lb = T, style = 1, char = "=")
  message("Analysis in progress. Please wait...")
  raster_Model <- pbapply::pblapply(Dates_extracted, call_barnes_opt)
  Ensamble <- rast(raster_Model)

  ##############################################################################
  #                           Perform validation if established                #
  ##############################################################################
  if (training != 1 | !is.null(stat_validation)) {
    final_results <- .validate(
      test_cords = data_val$test_cords,
      test_data = data_val$test_data,
      crss = coord.ref,
      Ensamble = Ensamble,
      Rain_threshold = Rain_threshold
    )
  }

  ##############################################################################
  #                           Save the model if necessary                      #
  ##############################################################################
  if (save_model) {
    message("Model saved successfully")
    terra::writeCDF(
      Ensamble,
      filename = "Barnes.nc",
      overwrite = TRUE
    )
  }


}
