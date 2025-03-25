#' Inverse distance weighted interpolation (IDW)
#'
#' @description
#' MS-GOP is a machine learning algorithm for merging satellite-based and ground precipitation data.
#' It combines Random Forest for spatial prediction, residual modeling for bias correction, and quantile mapping for final adjustment, ensuring accurate precipitation estimates across different temporal scales
#'
IDW <- function(BD_Obs, BD_Coord, shapefile, resolution, idp = 2) {
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
  data_IDW = as.data.frame(spl_layer, xy = TRUE)
  setDT(data_IDW)
  distancias = (distance(vect(data_IDW, geom = c("x", "y"), crs = "EPSG:32717"),
                         Points_VectTrain) / 1000) # Distancias en km

  weigths = function(d, idp) {
    return (1 / (d ^ idp))
  }

  for (i in seq_along(Points_VectTrain$Cod)) {
    data_IDW[, (Points_VectTrain$Cod[i]) := distancias[, i]]
  }

  estaciones = as.character(Points_VectTrain$Cod)
  for (est in estaciones) {
    data_IDW[, paste0("den_", est) := weigths(get(est), idp)]
  }

  w_cols = grep("^den_", names(data_IDW), value = TRUE) # Esto es directamente el 1 / d^idp
  stations = sub("^den_", "", w_cols)
  idw = function(data_obs) {
    for (i in seq_along(stations)) {
      station = stations[i]
      var_value = data_obs[Cod == station, var]
      data_IDW[, paste0("numer_", station) := ((var_value / get(w_cols[i])))]
      data_IDW[, paste0("denom_", station) := (1 / get(w_cols[i]))]
    }
    data_IDW[, var := rowSums(.SD, na.rm = TRUE) / rowSums(.SD, na.rm = TRUE),
             .SDcols = c(grep("^numer_", names(data_IDW), value = TRUE),
                         grep("^denom_", names(data_IDW), value = TRUE))]
    data_IDW = data_IDW[, .(x,y,var)]



  }



  idw_model = function(day) {
    data_obs <- training_data[Date == as.Date(day), ]
    if (sum(data_obs$var, na.rm = TRUE) == 0) {
      Ensamble <- spl_layer
    } else {

    } # End of if
  }



} # End of IDW function

# library(data.table)
# library(terra)
# BD_Obs = fread("C:/GitHub/InterpolateR/inst/extdata/BD_Insitu.csv")
# BD_Coord = fread("C:/GitHub/InterpolateR/inst/extdata/Cords_Insitu.csv")
# shapefile = vect("D:/ultimo_WINDOWS/SUBCUENCAS/TOMEBAMBA/CuencaRioTomebamba.shp")
# resolution = 5
# idp = 2
# plot(shapefile)
