# Declare global variables to prevent R CMD check warnings
utils::globalVariables(c("Cod", "ID", "Date", "var", "sum_n", "sum_d", "value", ".", "X", "Y", "Z", "x", "y"))
#' Inverse distance weighted interpolation (IDW)
#' @description
#' This function performs Inverse Distance Weighting (IDW), which is a spatial interpolation method that estimates unknown values using the known values of surrounding points, assigning greater influence to closer points.
#' @param BD_Obs A `data.table` or `data.frame` containing observational data with the following structure:
#'   - The first column (`Date`): A `Date` object representing the observation date.
#'   - The remaining columns: Each column corresponds to a unique ground station, where the column name is the station identifier.
#'
#'   The dataset should be structured as follows:
#'
#'   ```
#'   > BD_Obs
#'   # A data.table or data.frame with n rows (dates) and m+1 columns (stations + Date)
#'      Date        ST001  ST002  ST003  ST004  ...
#'      <date>      <dbl>  <dbl>  <dbl>  <dbl>  ...
#'   1  2015-01-01    0      0      0      0    ...
#'   2  2015-01-02    0      0      0     0.2   ...
#'   3  2015-01-03   0.1     0      0     0.1   ...
#'   ```
#'
#'   - Each station column contains numeric values representing observed measurements.
#'   - The column names (station identifiers) must be unique and match those in `BD_Coord$Cod` to ensure proper spatial referencing.
#' @param BD_Coord A `data.table` or `data.frame` containing the metadata of the ground stations. It must include the following columns:
#' #' - \code{"Cod"}:
#'    Unique identifier for each ground station.
#'
#' - \code{"X"}:
#'    Latitude of the station in UTM format.
#'
#' - \code{"Y"}:
#'    Longitude of the station in UTM format.
#' @param shapefile A shapefile defining the study area, used to constrain the interpolation to the region of interest.
#'   The shapefile must be of class `SpatVector` (from the `terra` package) and should have a UTM coordinate reference system.
#' @param resolution A numeric value indicating the resolution of the interpolation grid in kilometers (`km`).
#' @param p 'Numeric' value that controls how the influence decreases with distance. The default value is 2
#' @param n_round An integer specifying the number of decimal places to round the interpolated results.
#'   If set to `NULL`, all decimal places will be preserved. The default value is `1`.
#' @examples
#' \donttest{
#' library(InterpolateR)
#' # Load data from on-site observations
#'  data("BD_Obs", package = "InterpolateR")
#'  data("BD_Coord", package = "InterpolateR")
#'
#' # Load the study area where the interpolation is performed.
#'  shapefile <- terra::vect(system.file("extdata/study_area.shp", package = "InterpolateR"))
#'
#'  # Perform the interpolation
#'  Interpolated_data <- IDW(BD_Obs, BD_Coord, shapefile, resolution = 5, p = 2, n_round = 1)
#' }
#' @section Details:
#'  Inverse distance weighting (IDW) works as a deterministic mathematical interpolator that assumes
#'  that values closer to an unknown point are more closely related to it than those farther away. When using this method, sample points are weighted during the interpolation process so that the influence of a known point on the unknown point decreases as the distance between them increases.
#'  The IDW method calculates the value at an unknown point by a weighted average of the values of the known points, where the weights are inversely proportional to the distances between the prediction point and the known points. The basic formula defined by Shepard states that:
#'  \deqn{\hat{Z}(s_0) = \frac{\sum_{i=1}^{N} w_i Z(s_i)}{\sum_{i=1}^{N} w_i}}
#' where:
#' \describe{
#'   \item{\eqn{\hat{Z}(s_0)}}{ is the estimated value at the unknown point.}
#'   \item{\eqn{Z(s_i)}}{ are the known values.}
#'   \item{\eqn{w_i}}{ are the weights assigned to each known point.}
#'   \item{\eqn{N}}{ is the total number of known points used.}
#' }
#'
#' The weights are calculated by:
#'  \deqn{w_i = \frac{1}{d(s_0, s_i)^p}}
#' where:
#' \describe{
#'   \item{\eqn{d(s_0, s_i)}}{ is the distance between the unknown point \eqn{s_0} and the known point \eqn{s_i}.}
#'   \item{\eqn{p}}{ is the power parameter that controls how the influence decreases with distance.}
#' }
#' @return A `SpatRaster` object (from the `terra` package) where:
#'   - Each layer (`raster band`) corresponds to the interpolated values for a specific date in `BD_Obs$Date`.
#'   - The spatial resolution is defined by the `resolution` parameter (in km).
#'   - The coordinate reference system (CRS) matches that of the input `shapefile`.
#' @references Shepard, D. (1968) A two-dimensional interpolation function for irregularly-spaced data. Proceedings of the 1968 ACM National Conference, 1968, pages 517--524. DOI: 10.1145/800186.810616
#' @author Jonnathan Landi <jonnathan.landi@outlook.com>

IDW <- function(BD_Obs, BD_Coord, shapefile, resolution, p = 2, n_round = 1) {
  ##############################################################################
  #                               Check input data                             #
  ##############################################################################
  # shapefile must be a 'spatVector' object and coordinate reference system (CRS) must be defined
  if (!inherits(shapefile, "SpatVector")) stop("shapefile must be a 'SpatVector' object.")

  # BD_Obs can be a data.table or a data.frame
  if (!inherits(BD_Obs, c("data.table", "data.frame"))) stop("BD_Obs must be a 'data.table' or a 'data.frame'.")
  names(BD_Obs)[1] = "Date"

  # BD_Coord can be a data.table or a data.frame
  if (!inherits(BD_Coord, c("data.table", "data.frame"))) stop("BD_Coord must be a 'data.table' or a 'data.frame'.")

  # Check that the coordinate names appear in the observed data
  if (!all(BD_Coord$Cod %chin% setdiff(names(BD_Obs), "Date"))) stop("The names of the coordinates do not appear in the observed data.")

  # Check if Date column is present in BD_Obs
  ##############################################################################
  #                          Interpolation zone                                #
  ##############################################################################
  resolution = resolution * 1000 # Convert resolution to KM
  bbox = terra::ext(shapefile)
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

  square_vect <- terra::vect(square_polygon, type = "polygon", crs = terra::crs(shapefile))
  spl_layer <- rast(
    terra::ext(square_vect),
    resolution = resolution,
    crs = terra::crs(square_vect))
  terra::values(spl_layer) <- 0
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

  Points_VectTrain <- terra::vect(Points_Train, geom = c("X", "Y"), crs = terra::crs(spl_layer))
  ##############################################################################
  #                          IDW algotithm                                     #
  ##############################################################################
  data_IDW <- data.table::as.data.table(as.data.frame(spl_layer, xy = TRUE))
  data_IDW <- data_IDW[, .(X = x, Y = y)]
  coords <- as.matrix(data_IDW[, .(X, Y)])
  distancias <- data.table::as.data.table(terra::distance(terra::vect(coords, crs = terra::crs(spl_layer)), Points_VectTrain))
  setnames(distancias,  Points_VectTrain$Cod)
  data_IDW = cbind(data_IDW, distancias)

  estaciones = as.character(Points_VectTrain$Cod)
  denoms <- lapply(estaciones, function(est) 1 / (data_IDW[[est]]^p))
  denoms_dt <- setnames(data.table::as.data.table(denoms), paste0("d_", estaciones))

  idw = function(data_obs) {
    obs_values <- setNames(data_obs$var, data_obs$Cod)
    nums <- lapply(estaciones, function(est) {
      if (est %in% names(obs_values)) {
        obs_values[est] / (data_IDW[[est]]^p)
      } else {
        rep(NA_real_, nrow(data_IDW))
      }
    })
    nums_dt <- setnames(data.table::as.data.table(nums), paste0("n_", estaciones))
    result <- data_IDW[, .(X, Y)]
    result[, sum_n := rowSums(nums_dt, na.rm = TRUE)]
    result[, sum_d := rowSums(denoms_dt, na.rm = TRUE)]
    result[, value := sum_n / sum_d]
    result <- result[, .(X, Y, value)]
    result <- rast(result, crs = terra::crs(spl_layer))
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

#
