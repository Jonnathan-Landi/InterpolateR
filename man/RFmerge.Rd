% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RFmerge.R
\name{RFmerge}
\alias{RFmerge}
\title{Merging of satellite datasets with ground observations using Random Forest}
\usage{
RFmerge(
  BD_Obs,
  BD_Coord,
  cov,
  mask = NULL,
  n_round = NULL,
  ntree = 2000,
  seed = 123,
  training = 1,
  stat_validation = NULL,
  Rain_threshold = NULL,
  save_model = FALSE,
  name_save = NULL
)
}
\arguments{
\item{BD_Obs}{A \code{data.table} or \code{data.frame} containing observational data with the following structure:
\itemize{
\item The first column (\code{Date}): A \code{Date} object representing the observation date.
\item The remaining columns: Each column corresponds to a unique ground station, where the column name is the station identifier.
}

The dataset should be structured as follows:

\if{html}{\out{<div class="sourceCode">}}\preformatted{> BD_Obs
# A data.table or data.frame with n rows (dates) and m+1 columns (stations + Date)
   Date        ST001  ST002  ST003  ST004  ...
   <date>      <dbl>  <dbl>  <dbl>  <dbl>  ...
1  2015-01-01    0      0      0      0    ...
2  2015-01-02    0      0      0     0.2   ...
3  2015-01-03   0.1     0      0     0.1   ...
}\if{html}{\out{</div>}}
\itemize{
\item Each station column contains numeric values representing observed measurements.
\item The column names (station identifiers) must be unique and match those in \code{BD_Coord$Cod} to ensure proper spatial referencing.
}}

\item{BD_Coord}{A \code{data.table} or \code{data.frame} containing the metadata of the ground stations. It must include the following columns:
\itemize{
\item \code{"Cod"}:
Unique identifier for each ground station.
\item \code{"X"}:
Latitude of the station in UTM format.
\item \code{"Y"}:
Longitude of the station in UTM format.
\item \code{"Z"}:
Altitude of the station in meters.
}}

\item{cov}{A list of cov used as independent variables in the RFmerge. Each covariate should be a
\code{SpatRaster} object (from the \code{terra} package) and can represent satellite-derived weather variables or a Digital
Elevation Model (DEM). All cov should have the same number of layers (bands), except for the DEM, which must have only one layer.}

\item{mask}{A shapefile defining the study area.
If provided, the shapefile must be of class \code{SpatVector} (from the \code{terra} package) with a UTM coordinate reference system.
When specified, a spatial mask is applied to ensure that the final precipitation estimates are restricted to the defined study area.
Defaults to \code{NULL}, meaning no spatial mask is applied.}

\item{n_round}{An integer specifying the number of decimal places to round the interpolated results.
If set to \code{NULL}, all decimal places will be preserved. The default value is \code{1}.}

\item{ntree}{Numeric indicating the maximum number trees to grow in the Random Forest algorithm. The default value is set to 2000.
This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. If this value is too low, the prediction may be biased.}

\item{seed}{Integer for setting the random seed to ensure reproducibility of results (default: 123).}

\item{training}{Numerical value between 0 and 1 indicating the proportion of data used for model training. The remaining data are used for validation. Note that if you enter, for example, 0.8 it means that 80 \% of the data will be used for training and 20 \% for validation.
If you do not want to perform validation, set training = 1. (Default training = 1).}

\item{stat_validation}{A character vector specifying the names of the stations to be used for validation.
This option should only be filled in when it is desired to manually enter the stations used for validation. If this parameter is NULL, and the formation is different from 1, a validation will be performed using random stations.
The vector must contain the names of the stations selected by the user for validation.
For example, stat_validation = c(“ST001”, “ST002”). (Default stat_validation = NULL).}

\item{Rain_threshold}{List of numerical vectors defining precipitation thresholds to classify precipitation into different categories according to its intensity.
This parameter should be entered only when the validation is to include categorical metrics such as Critical Success Index (CSI), Probability of Detection (POD), False Alarm Rate (FAR), etc.
Each list item should represent a category, with the category name as the list item name and a numeric vector specifying the lower and upper bounds of that category.
\strong{Note:} See the "Notes" section for additional details on how to define categories, use this parameter for validation, and example configurations.}

\item{save_model}{Logical value indicating whether the interpolation file should be saved to disk. The default value is \code{FALSE}. indicating that the interpolated file should not be saved.
If set to \code{TRUE}, be sure to set the working directory beforehand using \code{setwd(path)} to specify where the files should be saved.}

\item{name_save}{Character string indicating the name under which the interpolation raster file will be saved. By default the algorithm sets as output name: 'Model_RFmerge'.
as the code will internally assign it.}
}
\value{
If a value other than 1 is set (point to pixel validation is performed), a list containing two elemeentis returned:

\strong{Ensamble:}
A \code{SpatRaster} object containing the bias-corrected layers for each time step. The number of layers
corresponds to the number of dates for which the correction is applied. This represents the corrected satellite data adjusted for bias.

\strong{Validation:}
A list containing the statistical results obtained from the validation process. This list includes:
\itemize{
\item \code{gof}:
A data table with goodness-of-fit metrics such as Kling-Gupta Efficiency (KGE), Nash-Sutcliffe Efficiency (NSE), Percent Bias (PBIAS), Root Mean Square Error (RMSE), and Pearson Correlation Coefficient (CC). These metrics assess the overall performance of the bias correction process.
\item \code{categorical_metrics}:
A data frame containing categorical evaluation metrics such as Probability of Detection (POD), Success Ratio (SR), False Alarm Rate (FAR), Critical Success Index (CSI), and Hit Bias (HB). These metrics evaluate the classification performance of rainfall event predictions based on user-defined precipitation thresholds.
}

If training is set to 1 (No validation is performed) only the Assembly mentioned above is returned.
}
\description{
RFmerge is a methodology developed by Baez-Villanueva et al. (2020) for the fusion of satellite precipitation datasets with ground-based observations, with the objective of improving the accuracy and spatial representativeness of the data.
This package implements RFmerge using Random Forest as a machine learning technique to correct biases and adjust the distribution of satellite products to in situ measurements.
In addition, it allows the integration of multiple sources of information, including geographic and environmental variables, optimizing the interpolation and spatial extrapolation of precipitation in data-limited regions.
Unlike previous implementations, this package has been optimized to improve computational efficiency and reduce processing times by incorporating advanced data manipulation techniques with \code{data.table}.
}
\section{Details}{

The \code{Rain_threshold} parameter is used to calculate categorical metrics such as the Critical Success Index (CSI),
Probability of Detection (POD), False Alarm Rate (FAR), success ratio (SR), Hit BIAS (HB),Heidke Skill Score (HSS);
Hanssen-Kuipers Discriminant (HK); Equal Threat Score (ETS) or Gilbert Skill Score.
The parameter should be entered as a named list, where each item represents a category and the name of the item is the category name.
The elements of each category must be a numeric vector with two values: the lower and upper limits of the category.
For example:
\code{Rain_threshold = list(
  no_rain = c(0, 1),
  light_rain = c(1, 5),
  moderate_rain = c(5, 20),
  heavy_rain = c(20, 40),
  violent_rain = c(40, Inf)
)}

Precipitation values will be classified into these categories based on their intensity.
Users can define as many categories as necessary, or just two (e.g., "rain" vs. "no rain").
It is important that these categories are entered according to the study region, as each study region may have its own categories.
}

\examples{
\donttest{
# Load data from on-site observations
 data("BD_Obs", package = "InterpolateR")
 data("BD_Coord", package = "InterpolateR")

# Load the cov
cov <- list(
 MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
 CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
 DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
 )

 # Apply the RFmerge
 model_RFmerge = RFmerge(BD_Obs, BD_Coord, cov, mask = NULL, n_round = 1, ntree = 2000,
                         seed = 123,  training = 0.8, stat_validation = NULL,
                         Rain_threshold = NULL, save_model = FALSE, name_save = NULL)

# Visualize the results
# Precipitation results within the study area
modelo_rainfall = model_RFmerge$Ensamble

# Validation statistic results
# goodness-of-fit metrics
metrics_gof = model_RFmerge$Validation$gof

# categorical metrics
metrics_cat = model_RFmerge$Validation$categorical_metrics
}
}
\references{
Baez-Villanueva, O. M.; Zambrano-Bigiarini, M.; Beck, H.; McNamara, I.; Ribbe, L.; Nauditt, A.; Birkel, C.; Verbist, K.; Giraldo-Osorio, J.D.; Thinh, N.X. (2020). RF-MEP: a novel Random Forest method for merging gridded precipitation products and ground-based measurements, Remote Sensing of Environment, 239, 111610. doi:10.1016/j.rse.2019.111606.
}
