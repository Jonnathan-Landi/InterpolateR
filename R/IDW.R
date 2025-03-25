#' Inverse distance weighted interpolation (IDW)
#'
#' @description
#' MS-GOP is a machine learning algorithm for merging satellite-based and ground precipitation data.
#' It combines Random Forest for spatial prediction, residual modeling for bias correction, and quantile mapping for final adjustment, ensuring accurate precipitation estimates across different temporal scales
#'
IDW <- function(BD_Obs, BD_Coord, shapefile, resolution, idp = 2, n_round = 1) {
  ##############################################################################
  #                    Check input data from on-site stations                  #
  ##############################################################################
  # BD_Obs can be a data.table or a data.frame
  if (!inherits(BD_Obs, c("data.table", "data.frame"))) stop("BD_Obs must be a 'data.table' or a 'data.frame'.")

  # BD_Coord can be a data.table or a data.frame
  if (!inherits(BD_Coord, c("data.table", "data.frame"))) stop("BD_Coord must be a 'data.table' or a 'data.frame'.")

  ##############################################################################
  #                          Interpolation zone                                #
  ##############################################################################
  resolution = resolution * 1000
  bbox = ext(shapefile)
  x_min = (mean(c(bbox$xmin, bbox$xmax))) - max((bbox$xmax - bbox$xmin), (bbox$ymax - bbox$ymin)) / 2
  x_max = (mean(c(bbox$xmin, bbox$xmax))) + max((bbox$xmax - bbox$xmin), (bbox$ymax - bbox$ymin)) / 2
  y_min = (mean(c(bbox$ymin, bbox$ymax))) - max((bbox$xmax - bbox$xmin), (bbox$ymax - bbox$ymin)) / 2
  y_max = (mean(c(bbox$ymin, bbox$ymax))) + max((bbox$xmax - bbox$xmin), (bbox$ymax - bbox$ymin)) / 2

  square_polygon <- rbind(
    c(x_min, y_min),
    c(x_max, y_min),
    c(x_max, y_max),
    c(x_min, y_max),
    c(x_min, y_min)
  )

  square_vect <- vect(square_polygon, type = "polygon", crs = crs(shapefile))
  spl_layer <- rast(
    ext(square_vect),
    resolution = resolution,
    crs = crs(square_vect))
  values(spl_layer) <- 0
  ##############################################################################
  #                           Data training                                    #
  ##############################################################################
  training_data <- melt(
    BD_Obs,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  )[, ID := as.numeric(factor(Cod))]
  Dates_extracted <- unique(training_data[, Date])
  Points_Train <- merge(training_data, BD_Coord, by = "Cod")
  setDT(Points_Train)

  Points_Train <- unique(Points_Train, by = "Cod")[, .(ID, Cod, X, Y, Z)]
  setorder(Points_Train, ID)

  Points_VectTrain <- terra::vect(Points_Train, geom = c("X", "Y"), crs = crs(spl_layer))
  ##############################################################################
  #                          IDW algotithm                                     #
  ##############################################################################
  data_IDW <- as.data.table(as.data.frame(spl_layer, xy = TRUE))
  data_IDW <- data_IDW[, .(x = x, y = y)]
  coords <- as.matrix(data_IDW[, .(x, y)])
  distancias <- as.data.table(distance(vect(coords, crs = crs(spl_layer)), Points_VectTrain) / 1000)
  setnames(distancias,  Points_VectTrain$Cod)
  data_IDW = cbind(data_IDW, distancias)

  estaciones = as.character(Points_VectTrain$Cod)
  denoms <- lapply(estaciones, function(est) 1 / (data_IDW[[est]]^idp))
  denoms_dt <- setnames(as.data.table(denoms), paste0("d_", estaciones))

  idw = function(data_obs) {
    obs_values <- setNames(data_obs$var, data_obs$Cod)
    nums <- lapply(estaciones, function(est) {
      if (est %in% names(obs_values)) {
        obs_values[est] / (data_IDW[[est]]^idp)
      } else {
        rep(NA_real_, nrow(data_IDW))
      }
    })
    nums_dt <- setnames(as.data.table(nums), paste0("n_", estaciones))
    result <- data_IDW[, .(x, y)]
    result[, sum_n := rowSums(nums_dt, na.rm = TRUE)]
    result[, sum_d := rowSums(denoms_dt, na.rm = TRUE)]
    result[, value := sum_n / sum_d]
    result <- result[, .(x, y, value)]
    result <- rast(result, crs = crs(spl_layer))
    return(result)
  }


  call_idw = function(day) {
    data_obs <- training_data[Date == as.Date(day), ]
    if (sum(data_obs$var, na.rm = TRUE) == 0) {
      return(spl_layer)
    } else {
      return ((idw(data_obs)))
    }
  }

  pbapply::pboptions(type = "timer", use_lb = T, style = 1, char = "=")
  message("Analysis in progress. Please wait...")
  raster_Model <- pbapply::pblapply(Dates_extracted, function(day) {
    call_idw(day)
  })
  Ensamble <- terra::rast(raster_Model)
  if (!is.null(n_round)) Ensamble <- terra::app(Ensamble, \(x) round(x, n_round))
  names(Ensamble) <- as.character(Dates_extracted)
  return(Ensamble)
}

