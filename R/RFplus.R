#' @title Machine learning algorithm for fusing ground and satellite precipitation data.
#'
#' @description
#' MS-GOP (RFplus) is a machine learning algorithm for merging satellite-based and ground precipitation data.
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
#' @param BD_Coord `data.table` containing metadata for the ground stations. Must include the following columns:
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
#' @param stat_validation A character vector specifying the names of the stations to be used for validation.
#'  This option should only be filled in when it is desired to manually enter the stations used for validation. If this parameter is NULL, and the formation is different from 1, a validation will be performed using random stations.
#'  The vector must contain the names of the stations selected by the user for validation.
#'  For example, stat_validation = c(“ST001”, “ST002”). (Default stat_validation = NULL).
#' @param Rain_threshold List of numerical vectors defining precipitation thresholds to classify precipitation into different categories according to its intensity.
#'  This parameter should be entered only when the validation is to include categorical metrics such as Critical Success Index (CSI), Probability of Detection (POD), False Alarm Rate (FAR), etc.
#'  Each list item should represent a category, with the category name as the list item name and a numeric vector specifying the lower and upper bounds of that category.
#'  \strong{Note:} See the "Notes" section for additional details on how to define categories, use this parameter for validation, and example configurations.
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
#' @param save_model Logical value indicating whether the interpolation file should be saved to disk. The default value is `FALSE`. indicating that the interpolated file should not be saved.
#'     If set to `TRUE`, be sure to set the working directory beforehand using `setwd(path)` to specify where the files should be saved.
#' @param name_save Character string indicating the name under which the interpolation raster file will be saved. By default the algorithm sets as output name: 'Model_RFplus'.
#' @examples
#' \donttest{
#' # Load the data
#'  data("BD_Obs", package = "InterpolateR")
#'  data("BD_Coord", package = "InterpolateR")
#'
#' # Load the covariates
#' Covariates <- list(
#'  MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
#'  CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
#'  DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
#'  )
#'
#'  # Apply the RFplus bias correction model
#' model = RFplus(BD_Obs, BD_Coord, Covariates, n_round = 1, wet.day = 0.1,
#'         ntree = 2000, seed = 123, training = 0.8,
#'         Rain_threshold = list(no_rain = c(0, 1), light_rain = c(1, 5)),
#'         method = "RQUANT", ratio = 10, save_model = FALSE, name_save = NULL)
#'
#' # Visualize the results
#' # Precipitation results within the study area
#' modelo_rainfall = model$Ensamble
#'
#' # Validation statistic results
#' # goodness-of-fit metrics
#' metrics_gof = model$Validation$gof
#'
#' # categorical metrics
#' metrics_cat = model$Validation$categorical_metrics
#' }
#'
#' @section Notes:
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
#' @export
RFplus <- function(
  BD_Obs,
  BD_Coord,
  Covariates,
  n_round = NULL,
  wet.day = FALSE,
  ntree = 2000,
  seed = 123,
  training = 1,
  stat_validation = NULL,
  Rain_threshold = NULL,
  method = c("RQUANT", "QUANT", "none"),
  ratio = 15,
  save_model = FALSE,
  name_save = NULL
) {
  ##############################################################################
  #                      Check the input data of the covariates                #
  ##############################################################################
  # Verify that covariates is a list
  if (!inherits(Covariates, "list")) stop("Covariates must be a list.")

  # Verify that the covariates are type SpatRaster
  if (!all(sapply(Covariates, function(x) inherits(x, "SpatRaster"))))
    stop("The covariates must be of type SpatRaster.")

  # Verify the extent of covariates
  ext_list = lapply(Covariates, terra::ext)
  if (!all(sapply(ext_list, function(x) x == ext_list[[1]])))
    stop(
      "The extension of the covariates are different (all extensions should be similar)."
    )

  # Verify the crc of covariates
  if (length(unique(vapply(Covariates, terra::crs, character(1)))) > 1)
    stop("The crs of the covariates are different (all crs should be similar).")

  ##############################################################################
  #                    Check input data from on-site stations                  #
  ##############################################################################
  # Verify the BD_Obs is data.table
  if (!inherits(BD_Obs, c("data.table", "data.frame")))
    stop("BD_Obs must be a 'data.table' or a 'data.frame'.")
  names(BD_Obs)[1] <- "Date"

  # Verify the columns of BD_Coord
  if (!inherits(BD_Coord, c("data.table", "data.frame")))
    stop("BD_Coord must be a 'data.table' or a 'data.frame'.")

  # Verify that all dates have at least one entry recorded
  Dates_NA <- BD_Obs[
    apply(BD_Obs[, .SD, .SDcols = -1], 1, function(x) all(is.na(x))),
    Date
  ]
  if (length(Dates_NA) > 0)
    stop(paste0(
      "No data was found for the dates: ",
      paste(Dates_NA, collapse = ", ")
    ))

  # Check that the coordinate names appear in the observed data
  if (!all(BD_Coord$Cod %in% setdiff(names(BD_Obs), "Date")))
    stop("The names of the coordinates do not appear in the observed data.")

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
  if (length(index_dem) == 0)
    stop(
      "A single layer covariate was not found. Possibly the DEM was not entered."
    )

  # Identify the DEM layer
  nlyrs_tots <- which(nlyr_covs != 1)
  nlyr_rep <- nlyr_covs[nlyrs_tots[1]]
  DEM <- Covariates[[index_dem]]
  ##############################################################################
  #                          Verify if validation is to be done                #
  ##############################################################################
  if (training != 1 | !is.null(stat_validation)) {
    data_val = .select_data(
      BD_Obs,
      BD_Coord,
      training = training,
      stat_validation = stat_validation
    )
    train_data = data_val$train_data
    train_cords = data_val$train_cords
  } else {
    message(
      "The training parameter was not entered. The model will be trained with all the data."
    )
    train_data <- BD_Obs
    train_cords <- BD_Coord
  }
  ##############################################################################
  #                         Prepare data for training                          #
  ##############################################################################
  # Layer to sample
  Sample_lyrs <- DEM * 0

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
  Points_VectTrain <- terra::vect(
    Points_Train,
    geom = c("X", "Y"),
    crs = terra::crs(Sample_lyrs)
  )

  # Calculate the Distance Euclidean
  distance_ED <- stats::setNames(
    lapply(1:nrow(Points_VectTrain), function(i) {
      terra::distance(DEM, Points_VectTrain[i, ], rasterize = FALSE)
    }),
    Points_VectTrain$Cod
  )

  difference_altitude <- stats::setNames(
    lapply(1:nrow(Points_VectTrain), function(i) {
      z_station = Points_VectTrain$Z[i]
      diff_alt = DEM[[1]] - z_station
      return(diff_alt)
    }),
    Points_VectTrain$Cod
  )

  ##############################################################################
  #                    Progressive correction methodology                      #
  ##############################################################################
  day_COV <- list(
    DEM = DEM,
    distance_ED = terra::rast(distance_ED),
    difference_altitude = terra::rast(difference_altitude)
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
  name_covs <- names(Covariates)[nlyr_covs != 1]
  data_simSat = lapply(name_covs, function(name) {
    raster = Covariates[[name]]
    dt = data.table::data.table(terra::extract(raster, y = Points_VectTrain))
    dt[, Cod := Points_VectTrain$Cod]
    dt[, ID := NULL]
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

  dt_sim = Reduce(
    function(x, y) merge(x, y, by = c("Date", "Cod"), all = TRUE),
    data_simSat
  )
  dt_merged = merge(dt_sim, data_cov, by = "Cod", all.x = TRUE)
  dt_final = merge(
    dt_merged,
    training_data[, .(Date, Cod, var)],
    by = c("Date", "Cod"),
    all.x = TRUE
  )
  ##############################################################################
  RF_Modelplus = function(date_P, Cov_pred_combined) {
    train_RF = dt_final[Date == date_P, ]
    if (train_RF[, sum(var, na.rm = TRUE)] == 0) return(Sample_lyrs)

    features_ff <- setdiff(names(train_RF), c("Cod", "Date"))
    set.seed(seed)
    Model_P1 <- randomForest::randomForest(
      var ~ .,
      data = train_RF[, ..features_ff],
      ntree = ntree,
      na.action = stats::na.omit
    ) |>
      suppressWarnings()

    val_RF = train_RF[,
      .(Cod, Obs = var, sim = stats::predict(Model_P1, .SD)),
      .SDcols = features_ff
    ]
    val_RF[, residuals := Obs - sim]

    # Model post-correction
    dt.train_resi <- data.table(
      residuals = val_RF$residuals,
      train_RF[,
        setdiff(names(train_RF), c("Cod", "Date", "var")),
        with = FALSE
      ]
    )

    set.seed(seed)
    Model_P2 <- randomForest::randomForest(
      residuals ~ .,
      data = dt.train_resi,
      ntree = ntree,
      na.action = stats::na.omit
    ) |>
      suppressWarnings()

    Pred1 = terra::predict(Cov_pred_combined, model = Model_P1, na.rm = TRUE)
    Pred2 = terra::predict(Cov_pred_combined, model = Model_P2, na.rm = TRUE)
    Ensamble = Pred1 + Pred2
    return(Ensamble)
  }

  Cov_pred <- Covariates[!grepl("DEM", names(Covariates))]
  pbapply::pboptions(type = "timer", use_lb = F, style = 1, char = "=")
  message("Analysis in progress. Please wait...")
  raster_Model <- pbapply::pblapply(Dates_extracted, function(date_P) {
    Cov_pred <- lapply(
      Cov_pred,
      function(x) x[[match(date_P, Dates_extracted)]]
    )
    Cov_pred_combined <- c(Cov_pred, list(day_COV))
    Cov_pred_combined <- terra::rast(Cov_pred_combined)
    ff = setdiff(names(dt_final), c("Cod", "Date", "var"))
    names(Cov_pred_combined) = ff
    RF_Modelplus(date_P, Cov_pred_combined)
  })

  Ensamble <- terra::rast(raster_Model)
  # Model of the QM or QDM correction ------------------------------------------
  if (method == "none") {
    message("Analysis completed, QUANT or RQUANT correction phase not applied.")
  } else if (method %in% c("RQUANT", "QUANT")) {
    message(paste0(
      "Analysis in progress: Stage 2 of 2. Correction by: ",
      method,
      ". Please wait..."
    ))

    data_CM <- data.table::data.table(terra::extract(
      Ensamble,
      Points_VectTrain
    ))
    data_CM[, ID := as.numeric(as.character(ID))]

    names_train = unique(training_data[, .(Cod)])
    names_train[, ID := seq_len(nrow(names_train))]
    data.table::setkey(names_train, ID)

    data_CM[, ID := names_train[data_CM, on = "ID", Cod]]
    data_CM <- stats::na.omit(data_CM)

    names <- as.character(data_CM$ID)

    data_CM <- data.table::data.table(t(data_CM[, -1]))
    colnames(data_CM) <- names
    data_CM <- data.table::data.table(
      data.table::data.table(Date = Dates_extracted),
      data_CM
    )

    common_columns <- setdiff(names(data_CM), c("ID", "Date"))

    res_interpolation <- lapply(common_columns, function(col) {
      dt <- merge(
        train_data[, .(Date, Obs = get(col))],
        data_CM[, .(Date, Sim = get(col))],
        by = "Date",
        all = FALSE
      )

      dt[,
        c("Obs", "Sim") := lapply(.SD, as.numeric),
        .SDcols = c("Obs", "Sim")
      ]

      return(stats::na.omit(dt))
    })

    names(res_interpolation) <- common_columns
    data_complete <- data.table::data.table(terra::as.data.frame(
      Ensamble,
      xy = TRUE
    ))
    data.table::setnames(
      data_complete,
      new = c("x", "y", as.character(Dates_extracted))
    )

    points <- train_cords[Cod %in% names, ]
    points <- terra::vect(
      points,
      geom = c("X", "Y"),
      crs = terra::crs(Sample_lyrs)
    )
    dat_final <- data.table::data.table()

    process_data <- function(method_fun, doQmap_fun) {
      cuantiles <- lapply(
        res_interpolation,
        function(x) method_fun(x$Obs, x$Sim, method = method, wet.day = wet.day)
      )
      message("Applying correction method. This may take a while...")

      dat_final <- pblapply(seq_len(nrow(data_complete)), function(i) {
        x <- data_complete[i, x]
        y <- data_complete[i, y]

        Points_VectTrain

        distances <- terra::distance(
          terra::vect(
            data.table::data.table(x, y),
            geom = c("x", "y"),
            crs = terra::crs(Sample_lyrs)
          ),
          points,
          unit = "km"
        )

        distances <- data.table::data.table(
          dist = as.vector(distances),
          Cod = points$Cod
        )

        if (any(distances$dist <= ratio)) {
          name <- distances[which.min(distances$dist), Cod]
          data <- data.table::data.table(
            Sim = t(data_complete[i, -c("x", "y")])
          )
          data_corregido <- doQmap_fun(data$Sim.V1, cuantiles[[name]])
          data_sat <- cbind(data_complete[i, c("x", "y")], t(data_corregido))
          colnames(data_sat) <- c("x", "y", as.character(Dates_extracted))

          return(data_sat)
        } else {
          return(data_complete[i])
        }
      })

      data.table::rbindlist(dat_final)
    }

    # Apply the correction method
    if (method == "QUANT") {
      dat_final = process_data(
        method_fun = qmap::fitQmapQUANT,
        doQmap_fun = qmap::doQmapQUANT
      )
    } else {
      dat_final <- process_data(
        method_fun = qmap::fitQmapRQUANT,
        doQmap_fun = qmap::doQmapRQUANT
      )
    }

    Ensamble <- terra::rast(dat_final, crs = terra::crs(Sample_lyrs))
    message("Analysis completed.")
  }

  ##############################################################################
  #                           Perform validation if established                #
  ##############################################################################
  if (training != 1 | !is.null(stat_validation)) {
    test_cords = data_val$test_cords
    test_data = data_val$test_data
    final_results <- .validate(
      test_cords,
      test_data,
      crss = terra::crs(Sample_lyrs),
      Ensamble,
      Rain_threshold = Rain_threshold
    )
  }

  if (!is.null(n_round))
    Ensamble <- terra::app(Ensamble, \(x) round(x, n_round))
  if (wet.day != FALSE)
    Ensamble <- terra::app(
      Ensamble,
      \(x) data.table::fifelse(x < wet.day, 0, x)
    )
  names(Ensamble) <- as.character(Dates_extracted)
  ##############################################################################
  #                           Save the model if necessary                      #
  ##############################################################################
  if (save_model) {
    message("Model saved successfully")
    if (is.null(name_save)) name_save = "Model_RFplus"
    name_saving <- paste0(name_save, ".nc")
    terra::writeCDF(Ensamble, filename = name_saving, overwrite = TRUE)
  }

  if (training != 1 | !is.null(stat_validation))
    return(list(Ensamble = Ensamble, Validation = final_results))
  if (training == 1 & is.null(stat_validation)) return(Ensamble)
} # End Rfplus functio
