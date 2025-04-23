test_that("Create data Method works correctly", {
  skip_on_cran()
  file.path = system.file("extdata/Folds_ejs_create_data", package = "InterpolateR")

  # Test with default parameters
  data = create_data(file.path, Start_date = "2015-01-01", End_Date = "2015-03-01",
                     ncores = NULL)
  print(head(data))

  # Test with different parameters
  data_2 = create_data(file.path, Start_date = "2015-01-01", End_Date = "2015-03-01",
                       ncores = NULL, max.na = 10)

  # Visualize the data
  print(head(data_2$data))

  # Visualize percentage of missing data
  print(head(data_2$Na_stations))

})
