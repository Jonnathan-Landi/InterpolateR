#' @title Ordinary Kriging interpolation
#' @description
#' This function performs Ordinary Kriging, which is a spatial interpolation method that provides the best linear unbiased estimator
#' for unknown values based on the spatial correlation structure of the observed data. Unlike IDW, Kriging incorporates the spatial
#' autocorrelation through variogram modeling to determine optimal weights.
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
#' @param variogram_model Character string specifying the variogram model to fit. Options are:
#'   - `"exponential"` (default): Exponential variogram model
#'   - `"spherical"`: Spherical variogram model
#'   - `"gaussian"`: Gaussian variogram model
#'   - `"linear"`: Linear variogram model
#' @param max_dist Numeric value specifying the maximum distance (in meters) for variogram calculation and kriging prediction.
#'   If `NULL` (default), it will be set to half the maximum distance between stations.
#' @param min_stations Integer specifying the minimum number of stations required for kriging interpolation.
#'   If fewer stations are available (due to NAs), the function will return a constant field or use fallback methods.
#'   Default is 2, which is more permissive for sparse precipitation data.
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
#' @param name_save Character string indicating the name under which the interpolation raster file will be saved. By default the algorithm sets as output name: 'Model_Kriging'.
#'
#' @param n_lags Integer specifying the number of lag bins to use in the empirical variogram calculation. Default is 15.
#'
#' @examples
#' \donttest{
#' # Load data from on-site observations
#' data("BD_Obs", package = "InterpolateR")
#' data("BD_Coord", package = "InterpolateR")
#'
#' # Load the study area where the interpolation is performed.
#' shp_path = system.file("extdata", "study_area.shp", package = "InterpolateR")
#' shapefile = terra::vect(shp_path)
#'
#' # Perform the interpolation
#' Interpolated_data <- Kriging_Ordinary(BD_Obs, BD_Coord, shapefile,
#'   grid_resolution = 5, variogram_model = "linear",
#'   max_dist = NULL, n_lags = 15, n_round = 1, training = 0.8,
#'   Rain_threshold = NULL, stat_validation = NULL,
#'   save_model = FALSE, name_save = NULL)
#' }
#'
#' @section Details:
#'  Ordinary Kriging is a geostatistical interpolation technique that provides the best linear unbiased estimator (BLUE)
#'  for spatial data. Unlike deterministic methods like IDW, Kriging incorporates the spatial correlation structure of the
#'  data through variogram modeling.
#'
#'  The Ordinary Kriging estimator is defined as:
#'  \deqn{\hat{Z}(s_0) = \sum_{i=1}^{n} \lambda_i Z(s_i)}
#' where:
#' \describe{
#'   \item{\eqn{\hat{Z}(s_0)}}{ is the estimated value at the unknown point.}
#'   \item{\eqn{Z(s_i)}}{ are the known values.}
#'   \item{\eqn{\lambda_i}}{ are the Kriging weights that sum to 1.}
#'   \item{\eqn{n}}{ is the total number of known points used.}
#' }
#'
#' The Kriging weights are obtained by solving the Kriging system:
#'  \deqn{\begin{bmatrix} \Gamma & \mathbf{1} \\ \mathbf{1}^T & 0 \end{bmatrix} \begin{bmatrix} \boldsymbol{\lambda} \\ \mu \end{bmatrix} = \begin{bmatrix} \boldsymbol{\gamma} \\ 1 \end{bmatrix}}
#' where:
#' \describe{
#'   \item{\eqn{\Gamma}}{ is the variogram matrix between observation points.}
#'   \item{\eqn{\boldsymbol{\gamma}}}{ is the variogram vector between prediction and observation points.}
#'   \item{\eqn{\mu}}{ is the Lagrange multiplier ensuring unbiasedness.}
#' }
#'
#' The variogram models available are:
#' \describe{
#'   \item{Exponential}{ \eqn{\gamma(h) = \sigma^2(1 - \exp(-3h/a))} for \eqn{h > 0}}
#'   \item{Spherical}{ \eqn{\gamma(h) = \sigma^2(1.5h/a - 0.5(h/a)^3)} for \eqn{h \leq a}, \eqn{\sigma^2} for \eqn{h > a}}
#'   \item{Gaussian}{ \eqn{\gamma(h) = \sigma^2(1 - \exp(-3(h/a)^2))} for \eqn{h > 0}}
#'   \item{Linear}{ \eqn{\gamma(h) = \sigma^2 \cdot h/a} for \eqn{h \leq a}, \eqn{\sigma^2} for \eqn{h > a}}
#' }
#'
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
#' @return The return value will depend on whether validation has been performed or not.
#' If validation is not performed, the function will return a `SpatRaster` object with the interpolated values.
#' If validation is performed, the function will return a list with two elements:
#'  - `Ensamble`: A `SpatRaster` object with the interpolated values.
#'  - `Validation`: A `data.table` with the validation results, including goodness-of-fit metrics and categorical metrics (if `Rain_threshold` is provided).
#' @references
#' Matheron, G. (1963). Principles of geostatistics. Economic Geology, 58(8), 1246-1266.
#' Cressie, N. (1993). Statistics for Spatial Data. John Wiley & Sons.
#' @author Marco Mogro <marcov.mogro@ucuenca.edu.ec>
#' @importFrom stats setNames optimize
#' @importFrom pbapply pboptions pblapply
#' @importFrom data.table melt setDT setorder as.data.table setnames := data.table
#' @export
#'

Kriging_Ordinary <- function(
    BD_Obs,
    BD_Coord,
    shapefile,
    grid_resolution,
    variogram_model = c("exponential", "spherical", "gaussian", "linear"),
    max_dist = NULL,
    n_lags = 15,
    min_stations = 2,
    n_round = NULL,
    training = 1,
    stat_validation = NULL,
    Rain_threshold = NULL,
    save_model = FALSE,
    name_save = NULL
) {
  ##############################################################################
  #                               Check input data                             #
  ##############################################################################
  # Shapefile must be a 'spatVector' object and coordinate reference system (CRS) must be defined
  if (!inherits(shapefile, "SpatVector")) {
    stop("shapefile must be a 'SpatVector' with a defined CRS.")
  }
  # BD_Obs must be a data.frame or data.table
  if (!inherits(BD_Obs, c("data.table", "data.frame"))) {
    stop("BD_Obs must be a 'data.frame' or 'data.table'.")
  }
  # BD_Coord must be a data.frame or data.table
  if (!inherits(BD_Coord, c("data.table", "data.frame"))) {
    stop("BD_Coord must be a 'data.frame' or 'data.table'.")
  }
  # Variogram model must be one of the specified options
  if (!variogram_model %in% c("exponential", "spherical", "gaussian", "linear")) {
    stop("variogram_model must be one of 'exponential', 'spherical', 'gaussian', or 'linear'.")
  }
  # grid_resolution must be numeric and positive
  if (!is.numeric(grid_resolution) || length(grid_resolution) != 1L) {
    stop("'grid_resolution' must be a single numeric value (km).")
  }
  # n_lags must be a positive integer
  if (!is.numeric(n_lags) || length(n_lags) != 1L || n_lags <= 0 || n_lags != floor(n_lags)) {
    stop("'n_lags' must be a single positive integer.")
  }
  # min_stations must be a positive integer
  if (!is.numeric(min_stations) || length(min_stations) != 1L || min_stations <= 0 || min_stations != floor(min_stations)) {
    stop("'min_stations' must be a single positive integer.")
  }
  # n_round must be NULL or a positive integer
  if (!is.null(n_round) && (!is.numeric(n_round) || length(n_round) != 1L || n_round < 0 || n_round != floor(n_round))) {
    stop("'n_round' must be NULL or a single non-negative integer.")
  }

  # Convert BD_Obs and BD_Coord to data.table if they are not already
  data.table::setDT(BD_Obs)
  data.table::setDT(BD_Coord)

  # Input validation (consolidated)
  stopifnot(
    inherits(shapefile, "SpatVector"),
    inherits(BD_Obs, c("data.table", "data.frame")),
    inherits(BD_Coord, c("data.table", "data.frame")),
    variogram_model %in% c("exponential", "spherical", "gaussian", "linear")
  )

  names(BD_Obs)[1] <- "Date"
  if (!all(BD_Coord$Cod %in% setdiff(names(BD_Obs), "Date"))) {
    stop("Coordinate names don't match observed data columns.")
  }

  # Setup training/validation data
  names_col <- setdiff(names(BD_Obs), "Date")
  Ids <- data.table::data.table(Cod = names_col, ID = seq_along(names_col))

  if (training != 1 || !is.null(stat_validation)) {
    data_val <- .select_data(BD_Obs, BD_Coord, training, stat_validation)
    train_data <- data_val$train_data
    train_cords <- data_val$train_cords
    message("Using training subset of data.")
  } else {
    train_data <- BD_Obs
    train_cords <- BD_Coord
    message("Training with all available data.")
  }

  # Setup interpolation grid
  coord.ref <- terra::crs(shapefile)
  spl_layer <- terra::rast(
    terra::ext(shapefile),
    resolution = grid_resolution * 1000,
    crs = coord.ref
  )
  terra::values(spl_layer) <- 0

  # Prepare kriging data structure
  Kriging_data <- data.table::melt(
    train_data,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  )
  Kriging_data <- Ids[Kriging_data, on = "Cod"]
  Dates_extracted <- unique(Kriging_data$Date)

  Points_Train <- unique(
    merge(Kriging_data, train_cords, by = "Cod"),
    by = "Cod"
  )[, .(ID, Cod, X, Y)]
  data.table::setorder(Points_Train, ID)
  Points_VectTrain <- terra::vect(
    Points_Train,
    geom = c("X", "Y"),
    crs = coord.ref
  )

  # Set maximum distance for variogram
  max_dist <- max_dist %||%
    (max(terra::distance(Points_VectTrain, Points_VectTrain)) / 2)

  # Variogram model functions (consolidated)
  variogram_models <- function(h, nugget, sill, range, model) {
    switch(
      model,
      "exponential" = nugget + sill * (1 - exp(-3 * h / range)),
      "spherical" = ifelse(
        h <= range,
        nugget + sill * (1.5 * h / range - 0.5 * (h / range)^3),
        nugget + sill
      ),
      "gaussian" = nugget + sill * (1 - exp(-3 * (h / range)^2)),
      "linear" = ifelse(h <= range, nugget + sill * h / range, nugget + sill)
    )
  }

  # Empirical variogram calculation (optimized)
  calc_empirical_variogram <- function(coords, values, n_lags, max_dist) {
    # Handle edge cases efficiently
    valid_idx <- !is.na(values)
    if (sum(valid_idx) < 2 || length(unique(values[valid_idx])) == 1) {
      lag_distances <- seq(
        max_dist / (2 * n_lags),
        max_dist - max_dist / (2 * n_lags),
        length.out = n_lags
      )
      return(data.frame(distance = lag_distances, gamma = rep(0.01, n_lags)))
    }

    coords <- coords[valid_idx, , drop = FALSE]
    values <- values[valid_idx]
    n <- length(values)
    lag_size <- max_dist / n_lags

    # Vectorized distance and gamma calculation
    idx_pairs <- expand.grid(i = 1:(n - 1), j = 2:n)
    idx_pairs <- idx_pairs[idx_pairs$i < idx_pairs$j, ]

    distances <- sqrt(rowSums(
      (coords[idx_pairs$i, ] - coords[idx_pairs$j, ])^2
    ))
    gamma_values <- 0.5 * (values[idx_pairs$i] - values[idx_pairs$j])^2

    # Filter by max distance and bin
    valid_pairs <- distances <= max_dist
    if (sum(valid_pairs) == 0) {
      lag_distances <- seq(lag_size / 2, max_dist - lag_size / 2, by = lag_size)
      return(data.frame(
        distance = lag_distances,
        gamma = rep(0.01, length(lag_distances))
      ))
    }

    lag_bins <- cut(
      distances[valid_pairs],
      breaks = seq(0, max_dist, by = lag_size),
      include.lowest = TRUE
    )
    empirical_gamma <- tapply(
      gamma_values[valid_pairs],
      lag_bins,
      mean,
      na.rm = TRUE
    )
    lag_distances <- seq(lag_size / 2, max_dist - lag_size / 2, by = lag_size)
    empirical_gamma[is.na(empirical_gamma)] <- 0.01

    data.frame(distance = lag_distances, gamma = as.numeric(empirical_gamma))
  }

  # Variogram fitting (streamlined)
  fit_variogram <- function(emp_variogram, model) {
    max_gamma <- max(emp_variogram$gamma, na.rm = TRUE)
    max_dist_vario <- max(emp_variogram$distance, na.rm = TRUE)

    # Handle flat variograms
    if (max_gamma <= 0.02) {
      return(list(
        nugget = 0.01,
        sill = 0.01,
        range = max_dist_vario / 3,
        model = model,
        is_flat = TRUE
      ))
    }

    # Optimization setup
    objective <- function(params) {
      if (any(params <= 0) || params[1] + params[2] > max_gamma * 2) {
        return(Inf)
      }
      pred_gamma <- variogram_models(
        emp_variogram$distance,
        params[1],
        params[2],
        params[3],
        model
      )
      sum((emp_variogram$gamma - pred_gamma)^2, na.rm = TRUE)
    }

    # Initial parameters and optimization
    init_params <- c(
      max(min(emp_variogram$gamma, na.rm = TRUE), 0.001),
      max(max_gamma - min(emp_variogram$gamma, na.rm = TRUE), 0.001),
      max_dist_vario / 3
    )

    result <- tryCatch(
      {
        stats::optim(
          init_params,
          objective,
          method = "L-BFGS-B",
          lower = c(0.001, 0.001, max_dist_vario / 100),
          upper = c(max_gamma, max_gamma, max_dist_vario)
        )
      },
      error = function(e) list(par = init_params)
    )

    list(
      nugget = result$par[1],
      sill = result$par[2],
      range = result$par[3],
      model = model,
      is_flat = FALSE
    )
  }

  # Setup prediction grid and distance matrices
  data_Kriging <- data.table::as.data.table(terra::as.data.frame(
    spl_layer,
    xy = TRUE
  ))[, .(X = x, Y = y)]
  coords_pred <- as.matrix(data_Kriging)
  coords_obs <- as.matrix(Points_Train[, .(X, Y)])

  dist_obs_obs <- as.matrix(terra::distance(Points_VectTrain, Points_VectTrain))
  dist_pred_obs <- as.matrix(terra::distance(
    terra::vect(coords_pred, crs = coord.ref),
    Points_VectTrain
  ))

  # Main kriging function (consolidated and optimized)
  perform_kriging <- function(data_obs, variogram_params) {
    obs_values <- stats::setNames(data_obs$var, data_obs$Cod)
    available_values <- obs_values[!is.na(obs_values)]
    available_stations <- names(available_values)

    # Handle insufficient data cases
    if (length(available_stations) < 2) {
      result <- data_Kriging[, .(
        X,
        Y,
        value = if (length(available_values) > 0) mean(available_values) else 0
      )]
      return(terra::rast(result, crs = coord.ref))
    }

    # Handle constant values
    if (length(unique(available_values)) == 1) {
      result <- data_Kriging[, .(X, Y, value = unique(available_values)[1])]
      return(terra::rast(result, crs = coord.ref))
    }

    station_indices <- which(Points_Train$Cod %in% available_stations)
    n_stations <- length(station_indices)

    # Use IDW for flat variograms or ill-conditioned matrices
    use_idw <- (!is.null(variogram_params$is_flat) && variogram_params$is_flat)

    if (!use_idw) {
      # Build gamma matrix
      gamma_matrix <- matrix(
        variogram_params$nugget,
        n_stations + 1,
        n_stations + 1
      )

      for (i in 1:n_stations) {
        for (j in 1:n_stations) {
          if (i != j) {
            h <- dist_obs_obs[station_indices[i], station_indices[j]]
            gamma_matrix[i, j] <- variogram_models(
              h,
              variogram_params$nugget,
              variogram_params$sill,
              variogram_params$range,
              variogram_params$model
            )
          }
        }
      }

      gamma_matrix[n_stations + 1, 1:n_stations] <- 1
      gamma_matrix[1:n_stations, n_stations + 1] <- 1
      gamma_matrix[n_stations + 1, n_stations + 1] <- 0

      # Check matrix condition
      use_idw <- (det(gamma_matrix[1:n_stations, 1:n_stations]) == 0 ||
                    kappa(gamma_matrix[1:n_stations, 1:n_stations]) > 1e12)
    }

    # Prediction loop (optimized)
    if (use_idw) {
      # Inverse Distance Weighting fallback
      predictions <- apply(
        dist_pred_obs[, station_indices, drop = FALSE],
        1,
        function(distances) {
          distances[distances == 0] <- 1e-10
          weights <- 1 / (distances^2)
          sum(weights * available_values) / sum(weights)
        }
      )
    } else {
      # Kriging predictions
      predictions <- apply(
        dist_pred_obs[, station_indices, drop = FALSE],
        1,
        function(distances) {
          gamma_vector <- c(
            variogram_models(
              distances,
              variogram_params$nugget,
              variogram_params$sill,
              variogram_params$range,
              variogram_params$model
            ),
            1
          )

          tryCatch(
            {
              weights <- solve(gamma_matrix, gamma_vector)
              max(0, sum(weights[1:n_stations] * available_values))
            },
            error = function(e) mean(available_values, na.rm = TRUE)
          )
        }
      )
    }

    result <- data_Kriging[, .(X, Y)]
    result[, value := predictions]
    terra::rast(result, crs = coord.ref)
  }

  # Daily kriging function (streamlined)
  process_day <- function(day) {
    data_obs <- Kriging_data[Date == day & !is.na(var)]

    if (nrow(data_obs) < 2 || all(data_obs$var == 0)) {
      return(spl_layer)
    }

    available_stations <- intersect(data_obs$Cod, Points_Train$Cod)
    if (length(available_stations) < 2) {
      return(spl_layer)
    }

    if (length(unique(data_obs$var)) == 1) {
      constant_raster <- spl_layer
      terra::values(constant_raster) <- unique(data_obs$var)[1]
      return(constant_raster)
    }

    # Calculate variogram and perform kriging
    station_indices <- which(Points_Train$Cod %in% available_stations)
    day_coords <- coords_obs[station_indices, , drop = FALSE]
    day_values <- data_obs$var[match(available_stations, data_obs$Cod)]

    tryCatch(
      {
        emp_variogram <- calc_empirical_variogram(
          day_coords,
          day_values,
          n_lags,
          max_dist
        )
        variogram_params <- fit_variogram(emp_variogram, variogram_model)
        perform_kriging(data_obs, variogram_params)
      },
      error = function(e) {
        fallback_raster <- spl_layer
        terra::values(fallback_raster) <- max(0, mean(day_values, na.rm = TRUE))
        fallback_raster
      }
    )
  }

  # Execute kriging for all dates
  pbapply::pboptions(type = "timer", use_lb = FALSE, style = 6, char = "=")
  message("Kriging analysis in progress. Please wait...")
  raster_Model <- pbapply::pblapply(Dates_extracted, process_day)

  Ensamble <- terra::rast(raster_Model)
  if (!is.null(n_round)) {
    Ensamble <- terra::app(Ensamble, \(x) round(x, n_round))
  }
  names(Ensamble) <- as.character(Dates_extracted)

  # Validation if required
  if (training != 1 || !is.null(stat_validation)) {
    final_results <- .validate(
      data_val$test_cords,
      data_val$test_data,
      coord.ref,
      Ensamble,
      Rain_threshold
    )
  }

  # Save model if requested
  if (save_model) {
    name_saving <- paste0(name_save %||% "Model_Kriging", ".nc")
    terra::writeCDF(Ensamble, filename = name_saving, overwrite = TRUE)
    message("Model saved successfully as ", name_saving)
  }

  # Return results
  if (training != 1 || !is.null(stat_validation)) {
    return(list(Ensamble = Ensamble, Validation = final_results))
  } else {
    return(Ensamble)
  }
}

