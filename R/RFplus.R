# Declare global variables to prevent R CMD check warnings
utils::globalVariables(c("Cod", "ID", "X", "Y", "Z","Date", "Obs", "sim", "residuals",
                         "var", ".", ".SD", "observed","estimated", "copy", ".N", "Sim"))
#' Machine learning algorithm for fusing ground and satellite precipitation data.
#'
#' @description
#' MS-GOP is a machine learning algorithm for merging satellite-based and ground precipitation data.
#' It combines Random Forest for spatial prediction, residual modeling for bias correction, and quantile mapping for final adjustment, ensuring accurate precipitation estimates across different temporal scales
#'
#' @details
#' The `RFplus` method implements a three-step approach:
#'
#' - \strong{Base Prediction}:
#'   Random Forest model is trained using satellite data and covariates.
#'
#' - \strong{Residual Correction}:
#'   A second Random Forest model is used to correct the residuals from the base prediction.
#'
#' - \strong{Distribution Adjustment}:
#'   Quantile mapping (QUANT or RQUANT) is applied to adjust the distribution of satellite data to match the observed data distribution.
#'
#' The final result combines all three steps, correcting the biases while preserving the outliers, and improving the accuracy of satellite-derived data such as precipitation and temperature.
#'
#' @param BD_Insitu `data.table` containing the ground truth measurements (dependent variable) used to train the RFplus model.
#'   Each column represents a ground station, and station identifiers are stored as column names. The class of `BD_Insitu`
#'   must be `data.table`. Each row represents a time step with measurements of the corresponding station.
#' @param Cords_Insitu `data.table` containing metadata for the ground stations. Must include the following columns:
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
#'
#' @param Covariates A list of covariates used as independent variables in the RFplus model. Each covariate should be a
#'   `SpatRaster` object (from the `terra` package) and can represent satellite-derived weather variables or a Digital
#'    Elevation Model (DEM). All covariates should have the same number of layers (bands), except for the DEM, which must have only one layer.
#' @param n_round Numeric indicating the number of decimal places to round the corrected values. If `n_round` is set to `NULL`, no rounding is applied.
#' @param wet.day Numeric value indicating the threshold for wet day correction. Values below this threshold will be set to zero.
#'   - `wet.day = FALSE`: No correction is applied (default).
#'   - For wet day correction, provide a numeric threshold (e.g., `wet.day = 0.1`).
#' @param ntree Numeric indicating the maximum number trees to grow in the Random Forest algorithm. The default value is set to 2000. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. If this value is too low, the prediction may be biased.
#' @param seed Integer for setting the random seed to ensure reproducibility of results (default: 123).
#' @param training Numerical value between 0 and 1 indicating the proportion of data used for model training. The remaining data are used for validation. Note that if you enter, for example, 0.8 it means that 80 % of the data will be used for training and 20 % for validation.
#' If you do not want to perform validation, set training = 1. (Default training = 1).
#' @param Rain_threshold
#' A list of numeric vectors that define the precipitation thresholds for classifying rainfall events into different categories based on intensity.
#' Each element of the list should represent a category, with the category name as the list element's name and a numeric vector specifying the lower and upper bounds for that category.
#'
#' \strong{Note:} See the "Notes" section for additional details on how to define categories, use this parameter for validation, and example configurations.
#' @param method
#' A character string specifying the quantile mapping method used for distribution adjustment. Options are:
#'
#' - \code{"RQUANT"}:
#'   Robust quantile mapping to adjust satellite data distribution to observed data.
#'
#' - \code{"QUANT"}:
#'   Standard quantile mapping.
#'
#' - \code{"none"}:
#'   No distribution adjustment is applied. Only Random Forest-based bias correction and residual correction are performed.
#' @param ratio integer Maximum search radius (in kilometers) for the quantile mapping setting using the nearest station. (default = 15 km)
#' @param save_model Logical value indicating whether the corrected raster layers should be saved to disk. The default is `FALSE`.
#'    If set to `TRUE`, make sure to set the working directory beforehand using `setwd(path)` to specify where the files should be saved.
#' @param name_save Character string. Base name for output file (default: NULL). The output file will be saved as "Model_RFplus.nc".
#' If you set a different name, make sure you do not set the ".nc" format,
#' as the code will internally assign it.
#' @param ... Additional arguments to pass to the underlying methods (e.g., for model tuning or future extensions).
#' @examples
#' \donttest{
#' # Load the libraries
#' library(terra)
#' library(data.table)
#'
#' # Load the data
#'  data("BD_Insitu", package = "RFplus")
#'  data("Cords_Insitu", package = "RFplus")
#'
#' # Load the covariates
#' Covariates <- list(
#'  MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "RFplus")),
#'  CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "RFplus")),
#'  DEM = terra::rast(system.file("extdata/DEM.nc", package = "RFplus"))
#'  )
#'
#'  # Apply the RFplus bias correction model
#' model = RFplus(BD_Insitu, Cords_Insitu, Covariates, n_round = 1, wet.day = 0.1,
#'         ntree = 2000, seed = 123, training = 1,
#'         Rain_threshold = list(no_rain = c(0, 1), light_rain = c(1, 5)),
#'         method = "RQUANT", ratio = 10, save_model = FALSE, name_save = NULL)
#' # Visualize the results
#'
#' # Precipitation results within the study area
#' modelo_rainfall = model$Ensamble
#'
#' # Validation statistic results
#' # goodness-of-fit metrics
#' metrics_gof = model$Validation$gof
#'
#' # categorical metrics
#' metrics_cat = model$Validation$categorical_metrics
#'
#' # Note: In the above example we used 80% of the data for training and 20% for # model validation.
#' }
#'
#' @section Notes:
#' The `Rain_threshold` parameter is used to classify precipitation events based on intensity into different categories. For example:
#'
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
#'
#' This parameter is required only when `training` is not equal to 1, as it is needed to calculate performance metrics such as the Probability of Detection (POD), False Alarm Rate (FAR), and Critical Success Index (CSI).
#'
#' @return A list containing two elements:
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
#' @author
#'  Jonnathan Augusto landi Bermeo, jonnathan.landi@outlook.com
#' @rdname RFplus
#' @export
RFplus <- function(BD_Insitu, Cords_Insitu, Covariates,...) {
  UseMethod("RFplus")
}

#' @rdname RFplus
#' @export
RFplus.default <- function(BD_Insitu, Cords_Insitu, Covariates, n_round = NULL, wet.day = FALSE,
                           ntree = 2000, seed = 123, training = 1, Rain_threshold = list(no_rain = c(0, 1)),
                           method = c("RQUANT","QUANT","none"), ratio = 15, save_model = FALSE, name_save = NULL, ...) {

  if (!inherits(BD_Insitu, "data.table")) stop("BD_Insitu must be a data.table.")

  RFplus.data.table(BD_Insitu = BD_Insitu, Cords_Insitu = Cords_Insitu, Covariates = Covariates,
                    n_round = n_round, wet.day = wet.day, ntree = ntree, seed = seed, training = training,
                    Rain_threshold = Rain_threshold,method = method, ratio = ratio, save_model = save_model,
                    name_save = name_save, ...)
}

#' @rdname RFplus
#' @export
RFplus.data.table <- function(BD_Insitu, Cords_Insitu, Covariates, n_round = NULL, wet.day = FALSE,
                              ntree = 2000, seed = 123, training = 1, Rain_threshold = list(no_rain = c(0, 1)), method = c("RQUANT","QUANT","none"),
                              ratio = 15, save_model = FALSE, name_save = NULL, ...) {

  ##############################################################################
  #                      Check the input data of the covariates                #
  ##############################################################################
  # Verify that covariates is a list
  if (!inherits(Covariates, "list")) stop("Covariates must be a list.")

  # Verify that the covariates are type SpatRaster
  if (!all(sapply(Covariates, function(x) inherits(x, "SpatRaster")))) stop("The covariates must be of type SpatRaster.")

  # Verify the extent of covariates
  ext_list = lapply(Covariates, terra::ext)
  if (!all(sapply(ext_list, function(x) x == ext_list[[1]]))) stop("The extension of the covariates are different (all extensions should be similar).")

  # Verify the crc of covariates
  if (length(unique(vapply(Covariates, terra::crs, character(1)))) > 1) stop("The crs of the covariates are different (all crs should be similar).")

  ##############################################################################
  #                    Check input data from on-site stations                  #
  ##############################################################################
  # Verify the BD_Insitu is data.table
  if (!inherits(BD_Insitu, "data.table")) stop("The data of the on-site stations should be a data.table.")

  # Verify the columns of Cords_Insitu
  if (!inherits(Cords_Insitu, "data.table")) stop("The coordinate data of the on-site stations must be a data.table.")

  # Check if there is a Date column
  date_column <- names(BD_Insitu)[which(sapply(BD_Insitu, function(x) inherits(x, c("Date", "IDate", "POSIXct"))))]
  if (length(date_column) == 0) stop("The Date column was not found in the observed data.")

  # Change the column name to match the full code
  if (date_column != "Date") setnames(BD_Insitu, date_column, "Date")

  # Verify that all dates have at least one entry recorded
  Dates_NA <- BD_Insitu[apply(BD_Insitu[, .SD, .SDcols = -1], 1, function(x) all(is.na(x))), Date]
  if (length(Dates_NA) > 0) stop(paste0("No data was found for the dates: ", paste(Dates_NA, collapse = ", ")))

  # Check that the coordinate names appear in the observed data
  if (!all(Cords_Insitu$Cod %chin% setdiff(names(BD_Insitu), "Date"))) stop("The names of the coordinates do not appear in the observed data.")

  ##############################################################################
  #                Checking the input parameters for quantile mapping          #
  ##############################################################################
  # Verify the method
  method <- match.arg(method, choices = c("RQUANT", "QUANT", "none"))

  ##############################################################################
  #               Verify that there is a DEM and manage DEM layers.            #
  ##############################################################################
  # Check if there is a DEM layer
  nlyr_covs <- sapply(Covariates, function(x) terra::nlyr(x))
  index_dem <- which(nlyr_covs == 1)
  if (length(index_dem) == 0) stop("A single layer covariate was not found. Possibly the DEM was not entered.")

  # Replicating the DEM at covariate scale
  nlyrs_tots <- which(nlyr_covs != 1)
  nlyr_rep <- nlyr_covs[nlyrs_tots[1]]
  DEM <- Covariates[[index_dem]]
  DEM <- terra::rast(replicate(nlyr_rep, DEM))
  Covariates[[index_dem]] <- DEM

  # Verify the layers of the covariates.
  if (length(unique(sapply(Covariates, function(x) terra::nlyr(x)))) > 1) stop("The number of covariate layers does not match. Check the input data.")

  ##############################################################################
  #                          Verify if validation is to be done                #
  ##############################################################################
  if (training != 1) {
    message(paste("The training parameter has been introduced. The model will be trained with:", (training * 100), "%", "data and validated with:", (100 - (training * 100)), "%"))

    # Verify if the number entered in training is valid
    if (!(training %between% c(0, 1))) stop("The training parameter must be between 0 and 1.")

    # Exclude Date column and split remaining columns
    set.seed(seed)
    columns <- setdiff(names(BD_Insitu), "Date")

    #  Randomly select the % of the training columns
    train_columns <- sample(columns, size = floor(training * length(columns)))

    # Data train
    train_data <- BD_Insitu[, .SD, .SDcols = c("Date", train_columns)]

    # Data test
    test_data <- BD_Insitu[, .SD, .SDcols = setdiff(names(BD_Insitu), train_columns)]

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
    train_cords <- Cords_Insitu[Cod %chin% train_columns, ]

    # Coordinates of the test data
    test_cords <- Cords_Insitu[Cod %chin% setdiff(names(BD_Insitu), train_columns), ]

  } else {
    message("The training parameter was not entered. The model will be trained with all the data.")
    train_data <- BD_Insitu
    train_cords <- Cords_Insitu
  }
  ##############################################################################
  #                         Prepare data for training                          #
  ##############################################################################
  # Layer to sample
  Sample_lyrs <- DEM[[1]]
  layer_0 = Covariates$DEM[[1]] * 0 # Layer when there is no precipitation.

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

  # Calculate altitude difference
  difference_altitude <- setNames(lapply(1:nrow(Points_VectTrain), function(i) {
    Covariates$DEM[[1]] - terra::extract(Covariates$DEM[[1]], Points_VectTrain[i, ])[, 2]
  }), Points_VectTrain$Cod)

  # Calculate the difference in altitude
  difference_altitude <- setNames(lapply(1:nrow(Points_VectTrain), function(i) {
    z_station = Points_VectTrain$Z[i]
    diff_alt = Sample_lyrs - z_station
    return(diff_alt)
  }), Points_VectTrain$Cod)
  ##############################################################################
  ##############################################################################
  #                    Progressive correction methodology                      #
  ##############################################################################
  ##############################################################################

  # Model of the Random Forest for the progressive correction 1 y 2 ------------
  RF_Modelplus = function(day_COV, fecha) {

    for (i in seq_along(day_COV)) {
      names(day_COV)[i] <- names(day_COV[[i]])
    }

    data_obs <- training_data[Date == as.Date(fecha), ]

    points_EstTrain = merge(
      data_obs[, .(ID, Cod)],
      Points_Train[, .(Cod, X, Y)],
      by = "Cod"
    )[order(ID)]

    points_EstTrain <- terra::vect(points_EstTrain, geom = c("X", "Y"), crs = crs(Sample_lyrs))

    # Covariates extras
    day_COV$dist_ED <- terra::rast(distance_ED[points_EstTrain$Cod])
    names(day_COV$dist_ED) = paste("dist_ED_", seq_along(points_EstTrain$Cod), sep = "")

    day_COV$diff_alt <- terra::rast(difference_altitude[points_EstTrain$Cod])
    names(day_COV$diff_alt) = paste("diff_alt_", seq_along(points_EstTrain$Cod), sep = "")

    #                   Training the model Random Forest (1)                   #
    data_cov <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
                       lapply(day_COV, terra::extract, y = points_EstTrain))

    dt.train <- merge(
      data_obs[, .(ID, Date, var)],
      data_cov,
      by = "ID"
    )

    dt.train <- dt.train[, setdiff(names(dt.train), "Date"), with = FALSE]

    # If the sum of var is 0, I assume that there is no precipitation in the whole basin.
    if (sum(dt.train$var, na.rm = TRUE) == 0) {
      Ensamble <- layer_0 # No precipitation
    } else {
      set.seed(seed)
      Model_P1 <- suppressWarnings(
        randomForest::randomForest(
          formula = var ~ .,
          data = dt.train[, setdiff(names(dt.train), "ID"), with = FALSE],
          ntree = ntree,
          na.action = stats::na.omit
        )
      )

      #                 Training the model Random Forest (2)                   #
      # Prediction of the Model 1
      val_P1 <- dt.train[, .(ID, Obs = var, sim = predict(Model_P1, .SD, na.rm = TRUE)), .SDcols = !"ID"]
      # val_p3 <- dt.train[, .(ID, Obs = var, sim = predict(Model_P1, .SD, na.rm = TRUE)), .SDcols = !"ID"]
      val_P1[, residuals := Obs - sim]

      dt.train_resi <- data.table(
        residuals = val_P1$residuals,
        dt.train[, setdiff(names(dt.train), c("ID", "var")), with = FALSE]
      )

      # Train Model 2
      Model_P2 <- suppressWarnings(randomForest::randomForest(
        formula = residuals ~ .,
        data = dt.train_resi,
        ntree = ntree,
        na.action = stats::na.omit
      )
      )
      # Create the corrected model 1
      cov_Sat <- terra::rast(day_COV)
      Ensamble <- predict(cov_Sat, Model_P1, na.rm = TRUE, fun = predict) +
        predict(cov_Sat, Model_P2, na.rm = TRUE, fun = predict)
      # Extra operations in case they have been established
      if (!is.null(n_round)) Ensamble <- terra::app(Ensamble, function(x) round(x, n_round))
      if (wet.day != FALSE) Ensamble <- terra::app(Ensamble, function(x) ifelse(x < wet.day, 0, x))
    } # End of the if statement
    return(Ensamble)
  }

  # Run the model
  pbapply::pboptions(type = "timer", use_lb = T, style = 1, char = "=")
  message("Analysis in progress: Stage 1 of 2. Please wait...")
  raster_Model <- pbapply::pblapply(Dates_extracted, function(fecha) {
    day_COV <- lapply(Covariates, function(x) x[[match(fecha, Dates_extracted)]])
    RF_Modelplus(day_COV, fecha)
  })

  Ensamble <- terra::rast(raster_Model)
  if (!is.null(n_round)) Ensamble <- terra::app(Ensamble, \(x) round(x, n_round))

  # Model of the QM or QDM correction ------------------------------------------
  if (method == "none") {
    message("Analysis completed, QUANT or RQUANT correction phase not applied.")
  } else if (method %in% c("RQUANT","QUANT")) {
    message(paste0("Analysis in progress: Stage 2 of 2. Correction by: ", method, ". Please wait..."))

    data_CM <- data.table(terra::extract(Ensamble, Points_VectTrain))
    data_CM[, ID := as.numeric(as.character(ID))]

    names_train = unique(training_data[, .(ID, Cod)])
    setkey(names_train, ID)

    data_CM[, ID := names_train[data_CM, on = "ID", Cod]]
    data_CM <- na.omit(data_CM)

    names <- as.character(data_CM$ID)

    data_CM <- data.table(t(data_CM[, -1]))
    colnames(data_CM) <- names
    data_CM <- data.table(data.table(Date = Dates_extracted), data_CM)

    common_columns <- setdiff(names(data_CM), c("ID", "Date"))

    res_interpolation <- lapply(common_columns, function(col) {
      dt <- merge(
        train_data[, .(Date, Obs = get(col))],
        data_CM[, .(Date, Sim = get(col))],
        by = "Date",
        all = FALSE
      )

      dt[, c("Obs", "Sim") := lapply(.SD, as.numeric), .SDcols = c("Obs", "Sim")]

      return(na.omit(dt))
    })

    names(res_interpolation) <- common_columns
    data_complete <- data.table(terra::as.data.frame(Ensamble, xy = TRUE))
    setnames(data_complete, new = c("x", "y", as.character(Dates_extracted)))

    points <- train_cords[Cod %chin% names, ]
    points <- terra::vect(points, geom = c("X", "Y"), crs = crs(Sample_lyrs))
    dat_final <- data.table()

    process_data <- function(method_fun, doQmap_fun) {
      cuantiles <- lapply(res_interpolation, function(x) method_fun(x$Obs, x$Sim, method = method, wet.day = wet.day))
      message("Applying correction method. This may take a while...")

      dat_final <- pblapply(seq_len(nrow(data_complete)), function(i) {
        x <- data_complete[i, x]
        y <- data_complete[i, y]

        Points_VectTrain

        distances <- terra::distance(
          terra::vect(data.table(x, y), geom = c("x", "y"), crs = crs(Sample_lyrs)),
          points, unit = "km"
        )

        distances <- data.table(dist = as.vector(distances), Cod = points$Cod)

        if (any(distances$dist <= ratio)) {
          name <- distances[which.min(distances$dist), Cod]
          data <- data.table(Sim = t(data_complete[i, -c("x", "y")]))
          data_corregido <- doQmap_fun(data$Sim.V1, cuantiles[[name]])
          data_sat <- cbind(data_complete[i, c("x", "y")], t(data_corregido))
          colnames(data_sat) <- c("x", "y", as.character(Dates_extracted))

          return(data_sat)
        } else {
          return(data_complete[i])
        }
      })

      rbindlist(dat_final)
    }

    # Apply the correction method
    if (method == "QUANT") {
      dat_final = process_data(method_fun = fitQmapQUANT, doQmap_fun = doQmapQUANT)
    } else {
      dat_final <- process_data(method_fun = fitQmapRQUANT, doQmap_fun = doQmapRQUANT)
    }

    Ensamble <- terra::rast(dat_final, crs = crs(Sample_lyrs))
    if (!is.null(n_round)) Ensamble <- terra::app(Ensamble, \(x) round(x, n_round))
    message("Analysis completed.")
  }

  ##############################################################################
  #                           Perform validation if established                #
  ##############################################################################
  if (training != 1) {
    message("Validation process in progress. Please wait.")
    test_cords$ID <- seq_len(nrow(test_cords))
    Points_VectTest <- terra::vect(test_cords, geom = c("X", "Y"), crs = crs(Sample_lyrs))

    data_validation = data.table(terra::extract(Ensamble, Points_VectTest))
    setkey(data_validation, ID)

    testing_data <- melt(
      test_data,
      id.vars = "Date",
      variable.name = "Cod",
      value.name = "var"
    )[, ID := as.numeric(Cod)]

    names_test = unique(testing_data[, .(ID, Cod)])
    setkey(names_test, ID)

    data_validation$ID = names_test[data_validation, on = "ID", Cod]
    data_validation <- data.table(t(data_validation[, -1, with = FALSE]))
    setnames(data_validation, new = as.character(names_test$Cod))

    data_validation <- data.table(testing_data[, .(Date)], data_validation)

    # Store in a list
    common_columns = setdiff(names(data_validation), "Date")

    res_validate <- lapply(common_columns, function(col) {
      merge(test_data[, .(Date, Obs = get(col))], data_validation[, .(Date, Sim = get(col))], by = "Date", all = FALSE)
    })
    names(res_validate) <- common_columns
    validation_results <- lapply(res_validate, function(x) evaluation_metrics(x, rain_thresholds = Rain_threshold))

    combine_results <- function(results, metric_type) {
      lapply(seq_along(results), function(i) {
        results[[i]][[metric_type]]$ID = names(results)[i]
        return(results[[i]][[metric_type]])
      })
    }

    gof_final_results <- data.table::rbindlist(combine_results(validation_results, "gof"))
    categorical_metrics_final_results <- data.table::rbindlist(combine_results(validation_results, "categorical_metrics"))
    final_results <- list(gof = gof_final_results, categorical_metrics = categorical_metrics_final_results)
  } else {
    final_results <- NULL
  }
  ##############################################################################
  #                           Save the model if necessary                      #
  ##############################################################################
  if (save_model) {
    message("Model saved successfully")
    if (is.null(name_save)) name_save = "Model_RFplus"
    name_saving <- paste0(name_save, ".nc")
    terra::writeCDF(Ensamble, filename = name_saving, overwrite=TRUE)
  }
  return(list(Ensamble = Ensamble, Validation = final_results))
}
