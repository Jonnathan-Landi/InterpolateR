#' Create a Unified Observation Dataset in the BD_Obs Format from Multiple CSV Files
#' This function constructs a unified dataset (`BD_Obs` structure) by merging multiple CSV files,
#' each containing in-situ observations from different stations. The function standardizes the format
#' required by downstream interpolation or bias correction algorithms by aligning all station data
#' into a single `data.table`, with dates as rows and station identifiers as columns.
#'
#' Each input CSV file must contain exactly two columns: the first with dates (`Date`) and the second
#' with the in-situ measurements of the variable to be interpolated.
#'
#' @param file.path `character`. Path to the folder containing the CSV files. Each file should represent a single station and be named using the station ID (e.g., `M001.csv`). Each file must have exactly two columns: a date column and a column of in-situ observations for the variable to be interpolated.
#' @param Start_date `Date`. Start date of the period to be included in the merged dataset. If the CSV files cover different date ranges, this defines the initial bound of the common time window for merging.
#' @param End_Date `Date`. End date of the period to be included in the merged dataset. This sets the upper bound of the time window to consider across all files.
#' @param ncores `integer`. Number of processing cores to be used when reading and merging CSV files in parallel. If If you want to perform the procedure without parallelization, set `ncores = NULL`. The default is `NULL`.
#' @param max.na `numeric`, optional. Maximum acceptable percentage of missing values per station (from 0 to 100). Stations exceeding this threshold will be excluded. If `NULL`, no filtering is performed. Default is `NULL`.
#'
#' @examples
#' \donttest{
#' # Example usage
#' file.path <- system.file("extdata/Folds_ejs_create_data", package = "InterpolateR")
#'
#' # Create a data with all stations
#' data <- create_data(file.path, Start_date = "2015-01-01", End_Date = "2015-03-01", ncores = NULL)
#' }
#' @return
#' If `max.na` is `NULL`, the function returns a `data.table` structured in the `BD_Obs` format,
#' where the first column contains the dates and the remaining columns correspond to individual stations.
#' This format preserves the full dataset without filtering for missing values.
#'
#' If `max.na` is not `NULL`, the function returns a named list containing:
#' \describe{
#'   \item{`data`}{A `data.table` in the `BD_Obs` format that includes only stations with a percentage
#'                of missing values less than or equal to `max.na`.}
#'   \item{`Na_stations`}{A `data.table` summarizing the percentage of missing values for each station,
#'                        useful for assessing data quality and supporting decisions about station selection.}
#' }
#'
#' @author Jonnathan Augusto landi Bermeo, jonnathan.landi@outlook.com
#' @importFrom data.table fread setnames data.table
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#'
#' @export
create_data <- function(
  file.path,
  Start_date,
  End_Date,
  ncores = NULL,
  max.na = NULL
) {
  file_list <- list.files(
    path = file.path,
    pattern = "*.csv",
    full.names = TRUE
  )
  station_names <- paste0(gsub("\\.csv$", "", basename(file_list)))

  read_data <- function(file) {
    name <- gsub("\\.csv$", "", basename(file))
    est <- data.table::fread(file)
    data.table::setnames(est, c("Date", name))
    return(est)
  }

  if (is.null(ncores)) {
    data_base <- lapply(file_list, read_data)
  } else {
    future::plan(multisession, workers = ncores)
    data_base <- future.apply::future_lapply(file_list, read_data)
    future::plan(sequential)
  }

  data_base <- lapply(data_base, function(dt) {
    dt[Date >= Start_date & Date <= End_Date]
  })

  data_base <- Reduce(
    function(x, y) merge(x, y, by = "Date", all = TRUE),
    data_base
  )
  if (is.null(max.na)) {
    return(data_base)
  }

  # Calculation of void percentage
  na_percentages <- data_base[,
    lapply(.SD, function(x) 100 * sum(is.na(x)) / .N),
    .SDcols = !c("Date")
  ]
  estat_valid <- names(na_percentages)[na_percentages < max.na]

  # Select stations that meet the criteria
  data_base <- data_base[, c("Date", estat_valid), with = FALSE]
  return(list(data = data_base, Na_stations = na_percentages))
}
