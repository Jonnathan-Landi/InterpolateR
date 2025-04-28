#' @title Inverse distance weighted interpolation (IDW)
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
#' - \code{"Cod"}:
#'    Unique identifier for each ground station.
#'
#' - \code{"X"}:
#'    Latitude of the station in UTM format.
#'
#' - \code{"Y"}:
#'    Longitude of the station in UTM format.
#' @param shapefile A shapefile defining the study area, used to constrain the interpolation to the region of interest.
#'   The shapefile must be of class `SpatVector` (from the `terra` package) and should have a UTM coordinate reference system.
#' @param grid_resolution A numeric value indicating the resolution of the interpolation grid in kilometers (`km`).
#' @param p 'Numeric' value that controls how the influence decreases with distance. The default value is 2
#' @param n_round An integer specifying the number of decimal places to round the interpolated results.
#'   If set to `NULL`, all decimal places will be preserved. The default value is `1`.
#' @param training Numerical value between 0 and 1 indicating the proportion of data used for model training. The remaining data are used for validation. Note that if you enter, for example, 0.8 it means that 80 % of the data will be used for training and 20 % for validation.
#' If you do not want to perform validation, set training = 1. (Default training = 1).
#' @param stat_validation A character vector specifying the names of the stations to be used for validation.
#'  This option should only be filled in when it is desired to manually enter the stations used for validation. If this parameter is NULL, and the formation is different from 1, a validation will be performed using random stations.
#'  The vector must contain the names of the stations selected by the user for validation.
#'  For example, stat_validation = c("ST001", "ST002"). (Default stat_validation = NULL).
#' @param Rain_threshold List of numerical vectors defining precipitation thresholds to classify precipitation into different categories according to its intensity.
#'  This parameter should be entered only when the validation is to include categorical metrics such as Critical Success Index (CSI), Probability of Detection (POD), False Alarm Rate (FAR), etc.
#'  Each list item should represent a category, with the category name as the list item name and a numeric vector specifying the lower and upper bounds of that category.
#'  \strong{Note:} See the "Notes" section for additional details on how to define categories, use this parameter for validation, and example configurations.
#' @param save_model Logical value indicating whether the interpolation file should be saved to disk. The default value is `FALSE`. indicating that the interpolated file should not be saved.
#'     If set to `TRUE`, be sure to set the working directory beforehand using `setwd(path)` to specify where the files should be saved.
#' @param name_save Character string indicating the name under which the interpolation raster file will be saved. By default the algorithm sets as output name: 'Model_IDW'.
#'
#' @examples
#' \donttest{
#' # Load data from on-site observations
#' data("BD_Obs", package = "InterpolateR")
#' data("BD_Coord", package = "InterpolateR")
#'
#' # Load the study area where the interpolation is performed.
#' shapefile <- terra::vect(system.file("extdata", "study_area.shp", package = "InterpolateR"))
#'
#' # Perform the interpolation
#' Interpolated_data <- IDW(BD_Obs, BD_Coord, shapefile,
#'   grid_resolution = 5, p = 2,
#'   n_round = 1, training = 0.8, Rain_threshold = NULL,
#'   stat_validation = NULL, save_model = FALSE, name_save = NULL
#' )
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
#' The `Rain_threshold` parameter is used to calculate categorical metrics such as the Critical Success Index (CSI),
#'  Probability of Detection (POD), False Alarm Rate (FAR), success ratio (SR), Hit BIAS (HB),Heidke Skill Score (HSS);
#'   Hanssen-Kuipers Discriminant (HK); Equal Threat Score (ETS) or Gilbert Skill Score.
#'   The parameter should be entered as a named list, where each item represents a category and the name of the item is the category name.
#'   The elements of each category must be a numeric vector with two values: the lower and upper limits of the category.
#'   For example:
#'  \code{Rain_threshold = list(
#'   no_rain = c(0, 1),
#'   light_rain = c(1, 5),
#'   moderate_rain = c(5, 20),
#'   heavy_rain = c(20, 40),
#'   violent_rain = c(40, Inf)
#' )}
#'
#' Precipitation values will be classified into these categories based on their intensity.
#' Users can define as many categories as necessary, or just two (e.g., "rain" vs. "no rain").
#' It is important that these categories are entered according to the study region, as each study region may have its own categories.
#' @return The return value will depend on whether validation has been performed or not.
#' If validation is not performed, the function will return a `SpatRaster` object with the interpolated values.
#' If validation is performed, the function will return a list with two elements:
#'  - `Ensamble`: A `SpatRaster` object with the interpolated values.
#'  - `Validation`: A `data.table` with the validation results, including goodness-of-fit metrics and categorical metrics (if `Rain_threshold` is provided).
#' @references Shepard, D. (1968) A two-dimensional interpolation function for irregularly-spaced data. Proceedings of the 1968 ACM National Conference, 1968, pages 517--524. DOI: 10.1145/800186.810616
#' @author Jonnathan Landi <jonnathan.landi@outlook.com>
#' @importFrom stats setNames
#' @importFrom pbapply pboptions pblapply
#' @importFrom data.table melt setDT setorder as.data.table setnames := data.table
#' @export
IDW <- function(BD_Obs, BD_Coord, shapefile, grid_resolution, p = 2,
                n_round = NULL, training = 1, stat_validation = NULL,
                Rain_threshold = NULL, save_model = FALSE, name_save = NULL) {
  ##############################################################################
  #                               Check input data                             #
  ##############################################################################
  if (!inherits(shapefile, "SpatVector")) stop("shapefile must be a 'SpatVector' object.")

  # BD_Obs can be a data.table or a data.frame
  if (!inherits(BD_Obs, c("data.table", "data.frame"))) stop("BD_Obs must be a 'data.table' or a 'data.frame'.")
  names(BD_Obs)[1] <- "Date"

  # BD_Coord can be a data.table or a data.frame
  if (!inherits(BD_Coord, c("data.table", "data.frame"))) stop("BD_Coord must be a 'data.table' or a 'data.frame'.")

  # Check that the coordinate names appear in the observed data
  if (!all(BD_Coord$Cod %in% base::setdiff(names(BD_Obs), "Date"))) stop("The names of the coordinates do not appear in the observed data.")
  ##############################################################################
  #                          Verify if validation is to be done                #
  ##############################################################################
  names_col <- base::setdiff(names(BD_Obs), "Date")
  Ids <- data.table::data.table(Cod = names_col, ID = 1:length(names_col))
  if (training != 1 | !is.null(stat_validation)) {
    data_val <- .select_data(BD_Obs, BD_Coord, training = training,
                           stat_validation = stat_validation)
    train_data <- data_val$train_data
    train_cords <- data_val$train_cords
  } else {
    message("The training parameter was not entered. The model will be trained with all the data.")
    train_data <- BD_Obs
    train_cords <- BD_Coord
  }
  ##############################################################################
  #                          Interpolation zone                                #
  ##############################################################################
  coord.ref <- terra::crs(shapefile)
  grid_resolution <- grid_resolution * 1000 # Convert resolution to KM

  spl_layer <- terra::rast(
    terra::ext(shapefile),
    resolution = grid_resolution,
    crs = coord.ref
  )

  terra::values(spl_layer) <- 0
  ##############################################################################
  #                           Data training                                    #
  ##############################################################################
  IDW_data <- data.table::melt(
    train_data,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  )

  IDW_data <- Ids[IDW_data, on = "Cod"]
  Dates_extracted <- base::unique(IDW_data[, Date])
  Points_Train <- base::merge(IDW_data, train_cords, by = "Cod")
  data.table::setDT(Points_Train)

  # Points_Train <- base::unique(Points_Train, by = "Cod")[, .(ID, Cod, X, Y, Z)]
  Points_Train <- unique(Points_Train, by = "Cod")[, .(ID, Cod, X, Y)]
  data.table::setorder(Points_Train, ID)

  Points_VectTrain <- terra::vect(Points_Train, geom = c("X", "Y"), crs = coord.ref)
  ##############################################################################
  #                          IDW algotithm                                     #
  ##############################################################################
  data_IDW <- data.table::as.data.table(terra::as.data.frame(spl_layer, xy = TRUE))
  data_IDW <- data_IDW[, .(X = x, Y = y)]
  coords <- as.matrix(data_IDW[, .(X, Y)])
  distancias <- data.table::as.data.table(terra::distance(terra::vect(coords, crs = coord.ref), Points_VectTrain))
  data.table::setnames(distancias, Points_VectTrain$Cod)
  data_IDW <- cbind(data_IDW, distancias)

  estaciones <- as.character(Points_VectTrain$Cod)
  denoms <- lapply(estaciones, function(est) 1 / (data_IDW[[est]]^p))
  denoms_dt <- data.table::setnames(data.table::as.data.table(denoms), paste0("d_", estaciones))

  idw <- function(data_obs) {
    obs_values <- stats::setNames(data_obs$var, data_obs$Cod)
    nums <- lapply(estaciones, function(est) {
      if (est %in% names(obs_values)) {
        obs_values[est] / (data_IDW[[est]]^p)
      } else {
        rep(NA_real_, nrow(data_IDW))
      }
    })

    nums_dt <- data.table::setnames(data.table::as.data.table(nums), paste0("n_", estaciones))
    result <- data_IDW[, .(X, Y)]
    result[, sum_n := rowSums(nums_dt, na.rm = TRUE)]
    result[, sum_d := rowSums(denoms_dt, na.rm = TRUE)]
    result[, value := sum_n / sum_d]
    result <- result[, .(X, Y, value)]
    result <- terra::rast(result, crs = coord.ref)
    return(result)
  }

  call_idw <- function(day) {
    data_obs <- IDW_data[Date == day, ]
    if (sum(data_obs$var, na.rm = TRUE) == 0) {
      return(spl_layer)
    } else {
      return((idw(data_obs)))
    }
  }

  pbapply::pboptions(type = "timer", use_lb = FALSE, style = 1, char = "=")
  message("Analysis in progress. Please wait...")
  raster_Model <- pbapply::pblapply(Dates_extracted, function(day) {
    call_idw(day)
  })

  Ensamble <- terra::rast(raster_Model)
  if (!is.null(n_round)) Ensamble <- terra::app(Ensamble, \(x) round(x, n_round))
  names(Ensamble) <- as.character(Dates_extracted)

  ##############################################################################
  #                           Perform validation if established                #
  ##############################################################################
  if (training != 1 | !is.null(stat_validation)) {
    test_cords <- data_val$test_cords
    test_data <- data_val$test_data
    final_results <- .validate(test_cords,  test_data, crss = coord.ref,
                              Ensamble, Rain_threshold = Rain_threshold)
  }
  ##############################################################################
  #                           Save the model if necessary                      #
  ##############################################################################
  if (save_model) {
    message("Model saved successfully")
    if (is.null(name_save)) name_save <- "Model_IDW"
    name_saving <- paste0(name_save, ".nc")
    terra::writeCDF(Ensamble, filename = name_saving, overwrite = TRUE)
  }
  ##############################################################################
  #                                      Return                                #
  ##############################################################################
  if (training != 1 | !is.null(stat_validation)) return(list(Ensamble = Ensamble, Validation = final_results))
  if (training == 1 & is.null(stat_validation)) return(Ensamble)
} # end funtion
