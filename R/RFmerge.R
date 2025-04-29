#' @title Merging of satellite datasets with ground observations using Random Forest
#'
#' @description
#' RFmerge is a methodology developed by Baez-Villanueva et al. (2020) for the fusion of satellite precipitation datasets with ground-based observations, with the objective of improving the accuracy and spatial representativeness of the data.
#' This package implements RFmerge using Random Forest as a machine learning technique to correct biases and adjust the distribution of satellite products to in situ measurements.
#' In addition, it allows the integration of multiple sources of information, including geographic and environmental variables, optimizing the interpolation and spatial extrapolation of precipitation in data-limited regions.
#' Unlike previous implementations, this package has been optimized to improve computational efficiency and reduce processing times by incorporating advanced data manipulation techniques with `data.table`.
#'
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
#'
#' - \code{"Z"}:
#'   Altitude of the station in meters.
#' @param cov A list of cov used as independent variables in the RFmerge. Each covariate should be a
#'   `SpatRaster` object (from the `terra` package) and can represent satellite-derived weather variables or a Digital
#'    Elevation Model (DEM). All cov should have the same number of layers (bands), except for the DEM, which must have only one layer.
#'
#' @param mask A shapefile defining the study area.
#'   If provided, the shapefile must be of class `SpatVector` (from the `terra` package) with a UTM coordinate reference system.
#'   When specified, a spatial mask is applied to ensure that the final precipitation estimates are restricted to the defined study area.
#'   Defaults to `NULL`, meaning no spatial mask is applied.
#' @param training Numerical value between 0 and 1 indicating the proportion of data used for model training. The remaining data are used for validation. Note that if you enter, for example, 0.8 it means that 80 % of the data will be used for training and 20 % for validation.
#' If you do not want to perform validation, set training = 1. (Default training = 1).
#' @param stat_validation A character vector specifying the names of the stations to be used for validation.
#'  This option should only be filled in when it is desired to manually enter the stations used for validation. If this parameter is NULL, and the formation is different from 1, a validation will be performed using random stations.
#'  The vector must contain the names of the stations selected by the user for validation.
#'  For example, stat_validation = c(“ST001”, “ST002”). (Default stat_validation = NULL).
#' @param seed Integer for setting the random seed to ensure reproducibility of results (default: 123).
#' @param ntree Numeric indicating the maximum number trees to grow in the Random Forest algorithm. The default value is set to 2000.
#' This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. If this value is too low, the prediction may be biased.
#' @param n_round An integer specifying the number of decimal places to round the interpolated results.
#'   If set to `NULL`, all decimal places will be preserved. The default value is `1`.
#' @param Rain_threshold List of numerical vectors defining precipitation thresholds to classify precipitation into different categories according to its intensity.
#'  This parameter should be entered only when the validation is to include categorical metrics such as Critical Success Index (CSI), Probability of Detection (POD), False Alarm Rate (FAR), etc.
#'  Each list item should represent a category, with the category name as the list item name and a numeric vector specifying the lower and upper bounds of that category.
#'  \strong{Note:} See the "Notes" section for additional details on how to define categories, use this parameter for validation, and example configurations.
#' @param save_model Logical value indicating whether the interpolation file should be saved to disk. The default value is `FALSE`. indicating that the interpolated file should not be saved.
#'     If set to `TRUE`, be sure to set the working directory beforehand using `setwd(path)` to specify where the files should be saved.
#' @param name_save Character string indicating the name under which the interpolation raster file will be saved. By default the algorithm sets as output name: 'Model_RFmerge'.
#' as the code will internally assign it.
#' @references Baez-Villanueva, O. M.; Zambrano-Bigiarini, M.; Beck, H.; McNamara, I.; Ribbe, L.; Nauditt, A.; Birkel, C.; Verbist, K.; Giraldo-Osorio, J.D.; Thinh, N.X. (2020). RF-MEP: a novel Random Forest method for merging gridded precipitation products and ground-based measurements, Remote Sensing of Environment, 239, 111610. doi:10.1016/j.rse.2019.111606.
#' @examples
#' \donttest{
#' # Load data from on-site observations
#'  data("BD_Obs", package = "InterpolateR")
#'  data("BD_Coord", package = "InterpolateR")
#'
#' # Load the cov
#' cov <- list(
#'  MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
#'  CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
#'  DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
#'  )
#'
#'  # Apply the RFmerge
#'  model_RFmerge = RFmerge(BD_Obs, BD_Coord, cov, mask = NULL, n_round = 1, ntree = 2000,
#'                          seed = 123,  training = 0.8, stat_validation = NULL,
#'                          Rain_threshold = NULL, save_model = FALSE, name_save = NULL)
#'
#' # Visualize the results
#' # Precipitation results within the study area
#' modelo_rainfall = model_RFmerge$Ensamble
#'
#' # Validation statistic results
#' # goodness-of-fit metrics
#' metrics_gof = model_RFmerge$Validation$gof
#'
#' # categorical metrics
#' metrics_cat = model_RFmerge$Validation$categorical_metrics
#' }
#' @return If a value other than 1 is set (point to pixel validation is performed), a list containing two elemeentis returned:
#'
#' \strong{Ensamble:}
#' A `SpatRaster` object containing the bias-corrected layers for each time step. The number of layers
#' corresponds to the number of dates for which the correction is applied. This represents the corrected satellite data adjusted for bias.
#'
#' \strong{Validation:}
#' A list containing the statistical results obtained from the validation process. This list includes:
#'
#' - \code{gof}:
#'   A data table with goodness-of-fit metrics such as Kling-Gupta Efficiency (KGE), Nash-Sutcliffe Efficiency (NSE), Percent Bias (PBIAS), Root Mean Square Error (RMSE), and Pearson Correlation Coefficient (CC). These metrics assess the overall performance of the bias correction process.
#'
#' - \code{categorical_metrics}:
#'   A data frame containing categorical evaluation metrics such as Probability of Detection (POD), Success Ratio (SR), False Alarm Rate (FAR), Critical Success Index (CSI), and Hit Bias (HB). These metrics evaluate the classification performance of rainfall event predictions based on user-defined precipitation thresholds.
#'
#' If training is set to 1 (No validation is performed) only the Assembly mentioned above is returned.
#' @section Details:
#' The `Rain_threshold` parameter is used to calculate categorical metrics such as the Critical Success Index (CSI),
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
#' @export
RFmerge = function(BD_Obs, BD_Coord, cov, mask = NULL, n_round = NULL, ntree = 2000,
                   seed = 123,  training = 1, stat_validation = NULL, Rain_threshold = NULL,
                   save_model = FALSE, name_save = NULL) {
  ##############################################################################
  #                               Check input data                             #
  ##############################################################################
  # Verify that cov is a list
  if (!inherits(cov, "list")) stop("cov must be a list.")

  # Verify that the cov are type SpatRaster
  if (!all(sapply(cov, function(x) inherits(x, "SpatRaster")))) stop("The cov must be of type SpatRaster.")

  # Verify the extent of cov
  ext_list = lapply(cov, terra::ext)
  if (!all(sapply(ext_list, function(x) x == ext_list[[1]]))) stop("The extension of the cov are different (all extensions should be similar).")

  # Verify the crc of cov
  if (length(unique(vapply(cov, terra::crs, character(1)))) > 1) stop("The crs of the cov are different (all crs should be similar).")
  if (!inherits(BD_Obs, c("data.table", "data.frame"))) stop("BD_Obs must be a 'data.table' or a 'data.frame'.")
  names(BD_Obs)[1] = "Date"

  # BD_Coord can be a data.table or a data.frame
  if (!inherits(BD_Coord, c("data.table", "data.frame"))) stop("BD_Coord must be a 'data.table' or a 'data.frame'.")

  # Check that the coordinate names appear in the observed data
  if (!all(BD_Coord$Cod %in% base::setdiff(names(BD_Obs), "Date"))) stop("The names of the coordinates do not appear in the observed data.")

  # Check if mask is a SpatVector object
  if (!is.null(mask) && !inherits(mask, "SpatVector")) stop("mask must be a 'SpatVector' object.")

  # Verify that all dates have at least one entry recorded
  Dates_NA <- BD_Obs[apply(BD_Obs[, .SD, .SDcols = -1], 1, function(x) all(is.na(x))), Date]
  if (length(Dates_NA) > 0) stop(paste0("No data was found for the dates: ", paste(Dates_NA, collapse = ", ")))
  ##############################################################################
  #               Verify that there is a DEM and manage DEM layers.            #
  ##############################################################################
  # Check if there is a DEM layer
  nlyr_covs <- sapply(cov, function(x) terra::nlyr(x))
  index_dem <- which(nlyr_covs == 1)
  if (length(index_dem) == 0) stop("A single layer covariate was not found. Possibly the DEM was not entered.")

  # Replicating the DEM at covariate scale
  nlyrs_tots <- which(nlyr_covs != 1)
  nlyr_rep <- nlyr_covs[nlyrs_tots[1]]
  DEM <- cov[[index_dem]]
  ##############################################################################
  #                          Verify if validation is to be done                #
  ##############################################################################
  if (training != 1 | !is.null(stat_validation)) {
    data_val = .select_data(BD_Obs, BD_Coord, training = training,
                            stat_validation = stat_validation)
    train_data = data_val$train_data
    train_cords = data_val$train_cords
  } else {
    message("The training parameter was not entered. The model will be trained with all the data.")
    train_data <- BD_Obs
    train_cords <- BD_Coord
  }
  ##############################################################################
  #                         Prepare data for training                          #
  ##############################################################################
  # Layer to sample
  Sample_lyrs <- DEM[[1]] * 0

  # Data for training
  training_data <- data.table::melt(
    train_data,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  )

  # Date of the data
  Dates_extracted <- unique(training_data[, Date])
  Points_Train <- merge(training_data, train_cords, by = "Cod")
  data.table::setDT(Points_Train)

  Points_Train <- unique(Points_Train, by = "Cod")[, .(Cod, X, Y, Z)]
  Points_VectTrain <- terra::vect(Points_Train, geom = c("X", "Y"), crs = terra::crs(Sample_lyrs))

  # Calculate the Distance Euclidean
  distance_ED <- stats::setNames(lapply(1:nrow(Points_VectTrain), function(i) {
    terra::distance(Sample_lyrs, Points_VectTrain[i, ], rasterize = FALSE)
  }), Points_VectTrain$Cod)
  ##############################################################################
  #                                 RF merge method                            #
  ##############################################################################
  day_COV <- list(
    DEM = DEM,
    distance_ED = terra::rast(distance_ED)
  )

  day_COV = terra::rast(day_COV)
  data_cov = lapply(day_COV, terra::extract, y = Points_VectTrain) |>
    Reduce(\(x, y) merge(x, y, by = "ID", all = TRUE), x = _) |>
    (\(d) {
      data.table::setDT(d)
    })()

  data_cov$DEM <- Points_VectTrain$Z
  data_cov$Cod <- Points_VectTrain$Cod
  data_cov$ID = NULL

  # Generate data for prediction
  name_covs <- names(cov)[nlyr_covs != 1]
  data_simSat = lapply(name_covs, function (name) {
    raster = cov[[name]]
    dt = data.table::data.table(terra::extract(raster, y = Points_VectTrain))
    dt[, Cod := Points_VectTrain$Cod]
    dt[,ID := NULL]
    dt = data.table::transpose(dt, make.names = "Cod")
    dt = cbind(
      Date = Dates_extracted,
      dt
    )

    dt <- data.table::melt(
      dt,
      id.vars = "Date",
      variable.name = "Cod",
      value.name = name
    )
    return(dt)
  })

  dt_sim =  Reduce(function(x, y) merge(x, y, by = c("Date", "Cod"), all = TRUE), data_simSat)
  dt_merged = merge(dt_sim, data_cov, by = "Cod", all.x = TRUE)
  dt_final = merge(dt_merged, training_data[, .(Date, Cod, var)], by = c("Date", "Cod"), all.x = TRUE)

  # Module for Ffmerge refactoring

  RF_Modelmerge = function(date_P, Cov_pred_combined) {
    train_RF = dt_final[Date == date_P, ]
    if (train_RF[, sum(var, na.rm = TRUE)] == 0) return(Sample_lyrs)
    features_ff <- base::setdiff(names(train_RF), c("Cod", "Date"))
    set.seed(seed)

    Model_P1 <- randomForest::randomForest(
      var ~ .,
      data = train_RF[, ..features_ff],
      ntree = ntree,
      na.action = stats::na.omit # REVISAR ESTO
     ) |>
      suppressWarnings()

    return(terra::predict(Cov_pred_combined, Model_P1, na.rm = TRUE))
  }

  # Run the model
  Cov_pred <- cov[!grepl("DEM", names(cov))]
  pbapply::pboptions(type = "timer", use_lb = F, style = 1, char = "=")
  message("Analysis in progress. Please wait...")
  raster_Model <- pbapply::pblapply(Dates_extracted, function(date_P) {
    Cov_pred <- lapply(Cov_pred, function(x) x[[match(date_P, Dates_extracted)]])
    Cov_pred_combined <- c(Cov_pred, list(day_COV))
    Cov_pred_combined <- terra::rast(Cov_pred_combined)
    ff = base::setdiff(names(dt_final), c("Cod", "Date", "var"))
    names(Cov_pred_combined) = ff
    RF_Modelmerge(date_P, Cov_pred_combined)
  })

  Ensamble <- terra::rast(raster_Model)

  if (!is.null(n_round)) Ensamble <- terra::app(Ensamble, \(x) round(x, n_round))

  if (!is.null(mask)) {
    Ensamble <- terra::mask(Ensamble, mask)
  }
  ##############################################################################
  #                           Perform validation if established                #
  ##############################################################################
  if (training != 1 | !is.null(stat_validation)) {
    test_cords = data_val$test_cords
    test_data = data_val$test_data
    final_results <- .validate(test_cords, test_data, crss = terra::crs(Sample_lyrs),
                              Ensamble, Rain_threshold = Rain_threshold)
  }
  ##############################################################################
  #                           Save the model if necessary                      #
  ##############################################################################
  if (save_model) {
    message("Model saved successfully")
    if (is.null(name_save)) name_save = "Model_RFmerge"
    name_saving <- paste0(name_save, ".nc")
    terra::writeCDF(Ensamble, filename = name_saving, overwrite=TRUE)
  }
  ##############################################################################
  #                                      Return                                #
  ##############################################################################
  if (training != 1 | !is.null(stat_validation)) return(list(Ensamble = Ensamble, Validation = final_results))
  if (training == 1 & is.null(stat_validation)) return(Ensamble)
} # End Rfmerge function
