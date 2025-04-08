create_data = function(file.path, Start_date, End_Date, ncores = 5, variable, max.na = NULL) {
  file_list <- list.files(path = file.path, pattern = "*.csv", full.names = TRUE)
  station_names <- paste0(gsub("\\.csv$", "", basename(file_list)))

  read_data = function(file){
    name = gsub("\\.csv$", "", basename(file))
    est = fread(file)
    est = est[, .(TIMESTAMP, get(variable))]
    setnames(est, c("Date", name))
    return(est)
  }

  plan(multisession, workers = ncores)
  data_base = future_lapply(file_list, read_data)
  plan(sequential)

  data_base = lapply(data_base, function(dt) {
    dt[Date >= Start_date & Date <= End_Date]
  })

  data_base = Reduce(function(x, y) merge(x, y, by = "Date", all = TRUE), data_base)
  if (is.null(max.na)) return(data_base)

  # Calculation of void percentage
  na_percentages <- data_base[, lapply(.SD, function(x) 100 * sum(is.na(x)) / .N), .SDcols = !c("Date")]
  estat_valid <- names(na_percentages)[na_percentages < max.na]

  # Select stations that meet the criteria
  data_base = data_base[, c("Date", estat_valid), with = FALSE]
  return(list(data = data_base, Na_stations = na_percentages))
}
