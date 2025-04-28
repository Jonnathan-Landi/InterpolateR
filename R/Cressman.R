#' @title Cressman Objective Analysis Method
#'
#' @description
#' The Cressman objective analysis computes values at grid points \eqn{Z_{ij}^a} (where \eqn{i} and \eqn{j} are the grid point indices for a 2D grid)
#' as the weighted average of the difference between observed values \eqn{Z_k^o} and background values interpolated to the
#' observation locations \eqn{Z_k^b} (i.e., \eqn{Z_k^o - Z_k^b}, called the observation increment) plus the background value
#' at the grid point \eqn{Z_{ij}^b}.
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
#' @param search_radius A numeric vector indicating the search radius in kilometers (`km`) for the Cressman method.
#'  \strong{Note:} See the "Notes" section for additional details on how to search radius values.
#' @param training Numerical value between 0 and 1 indicating the proportion of data used for model training. The remaining data are used for validation. Note that if you enter, for example, 0.8 it means that 80 % of the data will be used for training and 20 % for validation.
#' If you do not want to perform validation, set training = 1. (Default training = 1).
#' @param stat_validation A character vector specifying the names of the stations to be used for validation.
#'  This option should only be filled in when it is desired to manually enter the stations used for validation. If this parameter is NULL, and the formation is different from 1, a validation will be performed using random stations.
#'  The vector must contain the names of the stations selected by the user for validation.
#'  For example, stat_validation = c(“ST001”, “ST002”). (Default stat_validation = NULL).
#' @param Rain_threshold List of numerical vectors defining precipitation thresholds to classify precipitation into different categories according to its intensity.
#'  This parameter should be entered only when the validation is to include categorical metrics such as Critical Success Index (CSI), Probability of Detection (POD), False Alarm Rate (FAR), etc.
#'  Each list item should represent a category, with the category name as the list item name and a numeric vector specifying the lower and upper bounds of that category.
#'  \strong{Note:} See the "Notes" section for additional details on how to define categories, use this parameter for validation, and example configurations.
#' @param save_model Logical value indicating whether the interpolation file should be saved to disk. The default value is `FALSE`. indicating that the interpolated file should not be saved.
#'    If set to `TRUE`, be sure to set the working directory beforehand using `setwd(path)` to specify where the files should be saved.
#' @examples
#' \donttest{
#' # Load data from on-site observations
#'  data("BD_Obs", package = "InterpolateR")
#'  data("BD_Coord", package = "InterpolateR")
#'
#' # Load the study area where the interpolation is performed.
#'  shapefile <- terra::vect(system.file("extdata/study_area.shp", package = "InterpolateR"))
#'
#'  # Perform the interpolation
#'  Interpolated_Cressman <- Cressman(BD_Obs, BD_Coord, shapefile, grid_resolution = 5,
#'                                    search_radius = c(20,10), training = 1,
#'                                    stat_validation = "M001", Rain_threshold = NULL,
#'                                    save_model = FALSE)
#' # Results
#' Radius_20 = Interpolated_Cressman$Ensamble[[1]] # Interpolated data with a 20 km radius
#' Radius_10 = Interpolated_Cressman$Ensamble[[2]] # Interpolated data with a 10 km radius
#'
#' # Validation statistics
#' # Validation results with a 20 km radius
#' Validation_results_20 = Interpolated_Cressman$Validation[[1]]
#' # Validation results with a 10 km radius
#' Validation_results_10 = Interpolated_Cressman$Validation[[2]]
#'}
#' @section Details:
#' The Cressman method is defined by the following equation:
#' \deqn{Z_{ij}^a = Z_{ij}^b + \frac{\sum_{k=1}^{n} w_k (Z_k^o - Z_k^b)}{\sum_{k=1}^{n} w_k}}
#' where:
#' \describe{
#'  \item{\eqn{Z_{ij}^a}}{is the analysis value at grid point \eqn{i,j}.}
#'  \item{\eqn{Z_{ij}^b}}{is the background value at grid point \eqn{i,j}.}
#'  \item{\eqn{Z_k^o}}{is the observed value at station \eqn{k}.}
#'  \item{\eqn{Z_k^b}}{is the background value interpolated to station \eqn{k}.}
#'  \item{\eqn{w_k}}{is the weight assigned to station \eqn{k}.}
#'  \item{\eqn{n}}{is the total number of stations used.}
#'  }
#' The weight \eqn{w_k} is a function of the distance \eqn{r = \sqrt{(x_{ij} - x_k)^2 + (y_{ij} - y_k)^2}}
#' between the individual observation \eqn{k} and grid point \eqn{(i, j)}. \eqn{R} is the influence radius.
#' Beyond the influence radius, the weight is set to zero. \eqn{R} is therefore often referred to as
#' the cut-off radius.
#' @note
#' The `search_radius` parameter defines the influence range for the Cressman interpolation method.
#' It determines the maximum distance (in kilometers) within which observational data points contribute
#' to the interpolated value at a given location. A larger radius results in smoother interpolated fields
#' but may oversmooth local variations, while a smaller radius preserves finer details but may introduce noise.
#'
#' The Cressman method typically applies an iterative approach, where the search radius is progressively reduced
#' to refine the interpolation. Each iteration recalculates interpolated values with a smaller radius,
#' allowing a better representation of small-scale features in the dataset.
#'
#' The `search_radius` should be defined as a numeric vector representing the influence range in kilometers (`km`)
#' for each interpolation iteration. For example, setting `search_radius = c(50, 20, 10)` means the first iteration
#' considers a 50 km influence radius, the second iteration uses 20 km, and the final iteration refines
#' the interpolation with a 10 km radius. The number of elements in `search_radius` determines the total number of iterations.
#'
#'  The `Rain_threshold` parameter is used to calculate categorical metrics such as the Critical Success Index (CSI),
#'  Probability of Detection (POD), False Alarm Rate (FAR), success ratio (SR), Hit BIAS (HB),Heidke Skill Score (HSS);
#'   Hanssen-Kuipers Discriminant (HK); Equal Threat Score (ETS) or Gilbert Skill Score.
#'   The parameter should be entered as a named list, where each item represents a category and the name of the item is the category name.
#'   The elements of each category must be a numeric vector with two values: the lower and upper limits of the category.
#'   For example:
#' \code{Rain_threshold = list(
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
#' @return
#' The return value depends on whether validation has been performed.
#'
#' - **Without validation:** The function returns a `list`, where each element is a `SpatRaster` object containing the interpolated values for a specific search radius defined in `search_radius`. The number of elements in this list matches the length of `search_radius`.
#'
#' - **With validation:** The function returns a named `list` with two elements:
#'   - **`Ensamble`**: A `list` where each element corresponds to a `SpatRaster` object containing the interpolated values for a specific search radius in `search_radius`.
#'   - **`Validation`**: A `list` where each element is a `data.table` containing the validation results for the corresponding interpolated `SpatRaster`. Each `data.table` incluye métricas de bondad de ajuste como RMSE, MAE y Kling-Gupta Efficiency (KGE), junto con métricas categóricas si se proporciona `Rain_threshold`.
#'
#' The number of elements in both the `Ensamble` and `Validation` lists matches the length of `search_radius`, ensuring that each interpolation result has an associated validation dataset.
#' @references
#' Cressman, G. P., 1959: An operational objective analysis system. Mon. Wea. Rev., 87, 367-374, doi:10.1175/1520-0493(1959)087%3C0367:AOOAS%3E2.0.CO;2.
#' @author Jonnathan Landi <jonnathan.landi@outlook.com>
#' @export
Cressman <- function(BD_Obs, BD_Coord, shapefile, grid_resolution, search_radius,
                     training = 1, stat_validation = NULL, Rain_threshold = NULL,
                     save_model = FALSE) {
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
  if (!all(BD_Coord$Cod %in% base::setdiff(names(BD_Obs), "Date"))) stop("The names of the coordinates do not appear in the observed data.")
    ##############################################################################
  #                          Verify if validation is to be done                #
  ##############################################################################
  names_col <- base::setdiff(names(BD_Obs), "Date")
  Ids <- data.table::data.table(Cod = names_col, ID = 1:length(names_col))
  if (training != 1 | !is.null(stat_validation)) {
    data_val = .select_data(BD_Obs, BD_Coord, training = training,
                           stat_validation = stat_validation)
    train_data <- data_val$train_data
    train_cords <- data_val$train_cords
  } else {
    message("The training parameter was not entered. The model will be trained with all the data.")
    train_data <- BD_Obs
    train_cords <- BD_Coord
  }
    ##############################################################################
  #                           Zone of Interpolation                            #
  ##############################################################################
  grid_resolution <- grid_resolution * 1000
  coord.ref <- terra::crs(shapefile)
  spl_layer <- terra::rast(
    terra::ext(shapefile),
    resolution = grid_resolution,
    crs = coord.ref)
  terra::values(spl_layer) <- 0
    ##############################################################################
  #                           Data training                                    #
  ##############################################################################
  Cressman_data <- data.table::melt(
    train_data,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  )
  Cressman_data <- Ids[Cressman_data, on = "Cod"]
  Dates_extracted <- unique(Cressman_data[, Date])
  Points_Train <- merge(Cressman_data, train_cords, by = "Cod")
  data.table::setDT(Points_Train)

  Points_Train <- unique(Points_Train, by = "Cod")[, .(ID, Cod, X, Y, Z)]
  data.table::setorder(Points_Train, ID)

  Points_VectTrain <- terra::vect(Points_Train, geom = c("X", "Y"), crs = coord.ref)

  ##############################################################################
  #                          Cressman algotithm                                #
  ##############################################################################
  data_Crsmn <- data.table::as.data.table(as.data.frame(spl_layer, xy = TRUE))
  data.table::setnames(data_Crsmn, c("X", "Y", "value"))
  coords <- as.matrix(data_Crsmn[, .(X, Y)])

  distancias <- data.table::as.data.table(terra::distance(terra::vect(coords, crs = terra::crs(spl_layer)), Points_VectTrain))
  data.table::setnames(distancias, Points_VectTrain$Cod)
  data_Crsmn <- cbind(data_Crsmn[, .(X, Y)], distancias)

  # Pre-calcular matrices de pesos para cada radio
  search_radius <- search_radius * 1000
  estaciones <- Points_VectTrain$Cod

  weights <- lapply(search_radius, function(R) {
    R_sq <- R^2
    w <- (R_sq - data_Crsmn[, estaciones, with = FALSE]^2) / (R_sq + data_Crsmn[, estaciones, with = FALSE]^2)
    w[data_Crsmn[, estaciones, with = FALSE] > R] <- 0
    w
  })
  names(weights) <- as.character(search_radius)

  # Función optimizada para cálculo de Cressman
  crsmn_logic_opt <- function(data_obs) {
    obs_values <- data_obs$var[match(estaciones, data_obs$Cod)]
    obs_matrix <- matrix(obs_values, nrow = nrow(data_Crsmn), ncol = length(estaciones), byrow = TRUE)

    lapply(names(weights), function(r) {
      w <- as.matrix(weights[[r]])
      numer <- rowSums(w * obs_matrix, na.rm = TRUE)
      denom <- rowSums(w, na.rm = TRUE)
      result <- numer / denom

      # Crear raster directamente desde data.table
      terra::rast(
        data_Crsmn[, .(X, Y, value = ifelse(is.na(result) | denom == 0, NA, result))],
        crs = terra::crs(spl_layer),
        type = "xyz"
      )
    }) |> stats::setNames(names(weights))
  }

  # Función de llamada optimizada
  call_crsm_opt <- function(day) {
    data_obs <- Cressman_data[Date == day, .(Cod, var)]
    if (nrow(data_obs) == 0) return(NULL)
    crsmn_logic_opt(data_obs)
  }

  # Procesamiento paralelizable
  pbapply::pboptions(type = "timer", use_lb = F, style = 1, char = "=")
  message("Analysis in progress. Please wait...")
  raster_Model <- pbapply::pblapply(Dates_extracted, call_crsm_opt)

  # Crear stacks por radio
  Ensamble <- lapply(names(weights), function(r) {
    terra::rast(lapply(raster_Model, \(x) x[[r]]))
  }) |> stats::setNames(names(weights))

  ##############################################################################
  #                           Perform validation if established                #
  ##############################################################################
  if (training != 1 | !is.null(stat_validation)) {
    final_results <- lapply(names(Ensamble), function(name) {
      .validate(
        test_cords = data_val$test_cords,
        test_data = data_val$test_data,
        crss = coord.ref,
        Ensamble = Ensamble[[name]],
        Rain_threshold = Rain_threshold
      )
    }) |> stats::setNames(names(Ensamble))
  }
  ##############################################################################
  #                           Save the model if necessary                      #
  ##############################################################################
  if (save_model) {
    message("Model saved successfully")
    names_saving <- paste0("Radius_", names(Ensamble), ".nc")
    invisible(
      lapply(seq_along(Ensamble), function(i) {
        terra::writeCDF(
          Ensamble[[i]],
          filename = names_saving[i],
          overwrite = TRUE
        )
      })
    )
  }
  ##############################################################################
  #                                      Return                                #
  ##############################################################################
  if (training != 1 | !is.null(stat_validation)) return(list(Ensamble = Ensamble, Validation = final_results))
  if (training == 1 & is.null(stat_validation)) return(Ensamble)

} # End function cressman