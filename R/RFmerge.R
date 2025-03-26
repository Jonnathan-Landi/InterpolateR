#' Merging of satellite datasets with ground observations using RFmerge
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
#' #' - \code{"Cod"}:
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
#' @param cov A list of covariates used as independent variables in the RFmerge. Each covariate should be a
#'   `SpatRaster` object (from the `terra` package) and can represent satellite-derived weather variables or a Digital
#'    Elevation Model (DEM). All covariates should have the same number of layers (bands), except for the DEM, which must have only one layer.
#'
#' @param mask A shapefile defining the study area.
#'   If provided, the shapefile must be of class `SpatVector` (from the `terra` package) with a UTM coordinate reference system.
#'   When specified, a spatial mask is applied to ensure that the final precipitation estimates are restricted to the defined study area.
#'   Defaults to `NULL`, meaning no spatial mask is applied.
#' @param training Numerical value between 0 and 1 indicating the proportion of data used for model training. The remaining data are used for validation. Note that if you enter, for example, 0.8 it means that 80 % of the data will be used for training and 20 % for validation.
#' If you do not want to perform validation, set training = 1. (Default training = 1).
#' @param seed Integer for setting the random seed to ensure reproducibility of results (default: 123).
#' @param ntree Numeric indicating the maximum number trees to grow in the Random Forest algorithm. The default value is set to 2000.
#' This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. If this value is too low, the prediction may be biased.
#' @param n_round An integer specifying the number of decimal places to round the interpolated results.
#'   If set to `NULL`, all decimal places will be preserved. The default value is `1`.
#' @param Rain_threshold
#' A list of numeric vectors that define the precipitation thresholds for classifying rainfall events into different categories based on intensity.
#' Each element of the list should represent a category, with the category name as the list element's name and a numeric vector specifying the lower and upper bounds for that category.
#' @param save_model Logical value indicating whether the corrected raster layers should be saved to disk. The default is `FALSE`.
#'    If set to `TRUE`, make sure to set the working directory beforehand using `setwd(path)` to specify where the files should be saved.
#' @param name_save Character string. Base name for output file (default: NULL). The output file will be saved as "Model_RFplus.nc".
#' If you set a different name, make sure you do not set the ".nc" format,
#' as the code will internally assign it.
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
RFmerge = function(BD_Obs, BD_Coord, cov, mask = NULL, training, seed = 123,
                   ntree = 2000, n_round = NULL, Rain_threshold = list(no_rain = c(0, 1)),
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
  if (!all(BD_Coord$Cod %chin% setdiff(names(BD_Obs), "Date"))) stop("The names of the coordinates do not appear in the observed data.")

  # Check if mask is a SpatVector object
  if (!is.null(mask) && !inherits(mask, "SpatVector")) stop("mask must be a 'SpatVector' object.")

  # Verify that all dates have at least one entry recorded
  Dates_NA <- BD_Obs[apply(BD_Obs[, .SD, .SDcols = -1], 1, function(x) all(is.na(x))), Date]
  if (length(Dates_NA) > 0) stop(paste0("No data was found for the dates: ", paste(Dates_NA, collapse = ", ")))

  # Check that the coordinate names appear in the observed data
  if (!all(BD_Coord$Cod %chin% setdiff(names(BD_Obs), "Date"))) stop("The names of the coordinates do not appear in the observed data.")
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
  DEM <- terra::rast(replicate(nlyr_rep, DEM))
  cov[[index_dem]] <- DEM

  # Verify the layers of the cov.
  if (length(unique(sapply(cov, function(x) terra::nlyr(x)))) > 1) stop("The number of covariate layers does not match. Check the input data.")

  ##############################################################################
  #                          Verify if validation is to be done                #
  ##############################################################################
  if (training != 1) {
    message(paste("The training parameter has been introduced. The model will be trained with:", (training * 100), "%", "data and validated with:", (100 - (training * 100)), "%"))

    # Verify if the number entered in training is valid
    if (!(training %between% c(0, 1))) stop("The training parameter must be between 0 and 1.")

    # Exclude Date column and split remaining columns
    set.seed(seed)
    columns <- setdiff(names(BD_Obs), "Date")

    #  Randomly select the % of the training columns
    train_columns <- sample(columns, size = floor(training * length(columns)))

    # Data train
    train_data <- BD_Obs[, .SD, .SDcols = c("Date", train_columns)]

    # Data test
    test_data <- BD_Obs[, .SD, .SDcols = setdiff(names(BD_Obs), train_columns)]

    # Function used to validate data
    evaluation_metrics <- function(data, rain_thresholds) {
      ##############################################################################
      #                       metrics of goodness of fit                           #
      ##############################################################################
      gof = data.table(
        MAE = round(hydroGOF::mae(data$Sim, data$Obs, na.rm = T), 3),
        CC = round(hydroGOF::rSpearman(data$Sim, data$Obs, na.rm = T), 3),
        RMSE = round(hydroGOF::rmse(data$Sim, data$Obs, na.rm = T), 3),
        KGE = round(hydroGOF::KGE(data$Sim, data$Obs, na.rm = T), 3),
        NSE = round(hydroGOF::NSE(data$Sim, data$Obs, na.rm = T), 3),
        PBIAS = round(hydroGOF::pbias(data$Sim, data$Obs, na.rm = T), 3)
      )
      ##############################################################################
      #                       metrics of categorical                               #
      ##############################################################################
      create_threshold_categories <- function(rain_thresholds) {
        cat_names <- names(rain_thresholds)
        cat_min_values <- sapply(rain_thresholds, function(x) if(length(x) == 2) x[1] else x)

        # Sort categories by threshold value
        sorted_indices <- order(cat_min_values)
        cat_names <- cat_names[sorted_indices]

        # Create vectors for thresholds and categories
        thresholds <- c(sapply(rain_thresholds[cat_names], function(x) if(length(x) == 2) x[1] else x), Inf)

        return(list(thresholds = thresholds, categories = cat_names))
      }

      # Calculate performance metrics for each precipitation category
      calculate_category_metrics <- function(dt, category) {
        filtered_data <- dt[observed == category | estimated == category]

        # Calculate metrics
        hits <- filtered_data[observed == category & estimated == category, .N]
        misses <- filtered_data[observed == category & estimated != category, .N]
        false_alarms <- filtered_data[estimated == category & observed != category, .N]
        correct_negatives <- dt[observed != category & estimated != category, .N]

        # Calculate indices, handling zero denominators
        POD <- ifelse((hits + misses) > 0, hits / (hits + misses), NA) # Probability of Detection
        SR <- ifelse((hits + false_alarms) > 0, 1 - (false_alarms / (hits + false_alarms)), NA) # Success Ratio
        CSI <- ifelse((hits + misses + false_alarms) > 0, hits / (hits + misses + false_alarms), NA) #Critical Success Index
        HB <- ifelse((hits + misses) > 0, (hits + false_alarms) / (hits + misses), NA) # Hit BIAS
        FAR <- ifelse((hits + false_alarms) > 0, false_alarms / (hits + false_alarms), NA) # False Alarm Rate
        HK <- POD - (false_alarms / (false_alarms + correct_negatives)) # Hanssen-Kuipers Discriminant
        HSS <- ifelse((hits + misses)*(misses + correct_negatives) +
                        (hits + false_alarms)*(false_alarms + correct_negatives) != 0, (2 * (hits * correct_negatives - misses * false_alarms)) / (hits + misses)*(misses + correct_negatives) +
                        (hits + false_alarms)*(false_alarms + correct_negatives), NA) # Heidke Skill Score
        a_random <- ( (hits + false_alarms) * (hits + misses) ) /
          (hits + misses + false_alarms + correct_negatives)
        ETS <- ifelse(hits + misses + false_alarms - a_random != 0, hits - a_random /
                        hits + misses + false_alarms - a_random, NA) # Equitable Threat Score

        # Return results as data.table
        return(data.table(
          Category = category,
          POD = POD,
          SR = SR,
          CSI = CSI,
          HB = HB,
          FAR = FAR,
          HK = HK,
          HSS = HSS,
          ETS = ETS
        ))
      }

      dt <- copy(data)
      setkey(dt, Date)
      threshold_info <- create_threshold_categories(rain_thresholds)

      # Classify observed and estimated precipitation
      dt[, observed := ifelse(is.na(Obs), NA_character_,
                              as.character(cut(Obs, breaks = threshold_info$thresholds,
                                               labels = threshold_info$categories, right = FALSE)))]
      dt[, estimated := ifelse(is.na(Sim), NA_character_,
                               as.character(cut(Sim, breaks = threshold_info$thresholds,
                                                labels = threshold_info$categories, right = FALSE)))]

      # Create simplified data.table with just the categories
      category_data <- dt[, .(observed, estimated)]

      # Calculate metrics for each category and combine results
      categorical_metrics <- rbindlist(lapply(threshold_info$categories, function(cat) {
        calculate_category_metrics(category_data, cat)
      }))

      return(list(gof = gof, categorical_metrics = categorical_metrics))
    }

    # Coordinates of the training data
    train_cords <- BD_Coord[Cod %chin% train_columns, ]

    # Coordinates of the test data
    test_cords <- BD_Coord[Cod %chin% setdiff(names(BD_Obs), train_columns), ]

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
  training_data <- melt(
    train_data,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  )[, ID := as.numeric(factor(Cod))]

  # Date of the data
  Dates_extracted <- unique(training_data[, Date])
  Points_Train <- merge(training_data, train_cords, by = "Cod")
  setDT(Points_Train)

  Points_Train <- unique(Points_Train, by = "Cod")[, .(ID, Cod, X, Y, Z)]
  setorder(Points_Train, ID)

  Points_VectTrain <- terra::vect(Points_Train, geom = c("X", "Y"), crs = crs(Sample_lyrs))

  # Calculate the Distance Euclidean
  distance_ED <- setNames(lapply(1:nrow(Points_VectTrain), function(i) {
    terra::distance(Sample_lyrs, Points_VectTrain[i, ], rasterize = FALSE)
  }), Points_VectTrain$Cod)

  difference_altitude <- setNames(lapply(1:nrow(Points_VectTrain), function(i) {
    z_station = Points_VectTrain$Z[i]
    diff_alt = DEM[[1]] - z_station
    return(diff_alt)
  }), Points_VectTrain$Cod)
  ##############################################################################
  #                                 RF merge method                            #
  ##############################################################################
  RF_Modelmerge = function(day_COV, fecha) {
    names(day_COV) = sapply(day_COV, names)

    data_obs <- training_data[Date == as.Date(fecha), ]

    if (data_obs[, sum(var, na.rm = TRUE)] == 0) return(Sample_lyrs) # if the sum of var is 0, I assume that there is no precipitation in the whole basin.

    # If the sum of var not is 0, I assume that there is precipitation in the whole basin.
      points_EstTrain <- merge(
        data_obs[, .(ID, Cod)],
        Points_Train[, .(Cod, X, Y, Z)],
        by = "Cod"
      )[order(ID)] |>
        terra::vect(geom = c("X", "Y"), crs = crs(Sample_lyrs))

        add_rasters <- function(lyr, pattern) {
          r = terra::rast(get(lyr)[points_EstTrain$Cod])  # Obtiene el raster basado en la capa `lyr`
          names(r) = paste0(pattern, "_", seq_along(points_EstTrain$Cod))  # Usa `pattern` en los nombres
          r
        }

        day_COV$dist_ED <- add_rasters("distance_ED", "dist_ED")
        day_COV$diff_alt <- add_rasters("difference_altitude", "diff_alt")

          data_cov = lapply(day_COV, terra::extract, y = points_EstTrain) |>
            Reduce(\(x, y) merge(x, y, by = "ID", all = TRUE), x = _) |>
            (\(d) {
              setDT(d)
              d[, DEM := points_EstTrain$Z[match(ID, points_EstTrain$ID)]]
              d
            })()

            dt.train = merge(
              data_obs[, .(ID, var)],
              data_cov,
              by = "ID"
            )
            features <- setdiff(names(dt.train), "ID")
            cov_Sat <- terra::rast(day_COV)
                Model_P1 <- randomForest::randomForest(
                  var ~ .,
                  data = dt.train[, ..features],
                  ntree = ntree,
                  na.action = na.omit
                ) |>
                  suppressWarnings()

                return(terra::predict(cov_Sat, Model_P1, na.rm = TRUE))
  }

  # Run the model
  pbapply::pboptions(type = "timer", use_lb = T, style = 1, char = "=")
  message("Analysis in progress: Stage 1 of 2. Please wait...")
  raster_Model <- pbapply::pblapply(Dates_extracted, function(fecha) {
    day_COV <- lapply(cov, function(x) x[[match(fecha, Dates_extracted)]])
    RF_Modelmerge(day_COV, fecha)
  })

  Ensamble <- terra::rast(raster_Model)
  if (!is.null(n_round)) Ensamble <- terra::app(Ensamble, \(x) round(x, n_round))

  if (!is.null(mask)) {
    Ensamble <- terra::mask(Ensamble, mask)
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
  if (training == 1) return(Ensamble)
  if (training != 1) return(list(Ensamble = Ensamble, Validation = final_results))

}
