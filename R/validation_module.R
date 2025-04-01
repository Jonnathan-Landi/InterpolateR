# This function prepares the files for validation
.select_data = function(BD_Obs, BD_Coord, training, seed, stat_validation = NULL) {
  if (is.null(stat_validation)) {
    message(paste("Random validation mode.", (training * 100), "% of the stations will be used for training and", (100 - (training * 100)), "% for validation."))
    if (!(training %between% c(0, 1))) stop("The training parameter must be between 0 and 1.")
    set.seed(seed)
    columns <- setdiff(names(BD_Obs), "Date")
    train_columns <- sample(columns, size = floor(training * length(columns)))
  } else {
    message(paste("Manual validation mode. Stations:", paste(stat_validation, collapse = ", "), "will be used for validation."))
    train_columns <-  setdiff(names(BD_Obs),  c("Date", stat_validation))
  }

  # Data train
  train_data <- BD_Obs[, .SD, .SDcols = c("Date", train_columns)]

  # Data test
  test_data <- BD_Obs[, .SD, .SDcols = setdiff(names(BD_Obs), train_columns)]

  # Coordinates of the training data
  train_cords <- BD_Coord[Cod %chin% train_columns, ]

  # Coordinates of the test data
  test_cords <- BD_Coord[Cod %chin% setdiff(names(BD_Obs), train_columns), ]
  return(list(train_data = train_data, test_data = test_data, train_cords = train_cords, test_cords = test_cords))
}

# This function contains the evaluation metrics
metrics <- function(data, Rain_threshold) {
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
  if (!is.null(Rain_threshold )) {
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

    setkey(data, Date)
    threshold_info <- create_threshold_categories(Rain_threshold)

    # Classify observed and estimated precipitation
    data[, observed := ifelse(is.na(Obs), NA_character_,
                              as.character(cut(Obs, breaks = threshold_info$thresholds,
                                               labels = threshold_info$categories, right = FALSE)))]
    data[, estimated := ifelse(is.na(Sim), NA_character_,
                               as.character(cut(Sim, breaks = threshold_info$thresholds,
                                                labels = threshold_info$categories, right = FALSE)))]

    # Create simplified data.table with just the categories
    category_data <- data[, .(observed, estimated)]

    # Calculate metrics for each category and combine results
    categorical_metrics <- rbindlist(lapply(threshold_info$categories, function(cat) {
      calculate_category_metrics(category_data, cat)
    }))

    return(list(gof = gof, categorical_metrics = categorical_metrics))
  } else {
    return(gof)
  }
}

# This function validates the model
.validate = function(test_cords, test_data, crss, Ensamble, Rain_threshold) {
  message("Validation process in progress. Please wait.")
  test_cords$ID <- seq_len(nrow(test_cords))
  Points_VectTest <- terra::vect(test_cords, geom = c("X", "Y"), crs = crss)

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

  common_columns = setdiff(names(data_validation), "Date")

  res_validate <- lapply(common_columns, function(col) {
    merge(test_data[, .(Date, Obs = get(col))], data_validation[, .(Date, Sim = get(col))], by = "Date", all = FALSE)
  })

  names(res_validate) <- common_columns
  validation_results <- lapply(res_validate, function(x) metrics(x, Rain_threshold  = Rain_threshold))
  if (!is.null(Rain_threshold)) {

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
    final_results <- data.table::rbindlist(validation_results, idcol = "ID")
  }
  return(final_results)
}

