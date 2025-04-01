---
title: "A Comprehensive Toolkit for Fast and Efficient Spatial Interpolation"
author: "Jonnathan Landi"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{InterpolateR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduction

`InterpolateR` is one of the most efficient R packages for spatial interpolation. It leverages vectorized operations to maximize computational efficiency while maintaining high accuracy in interpolation tasks. The package integrates both traditional and advanced interpolation methods, making it a versatile tool for a wide range of spatial analysis applications.

Among the traditional interpolation techniques included are:

-   **Inverse Distance Weighting (IDW)**: A widely used method that assigns weights to observations based on their inverse distance.

-   **Cressman Objective Analysis Method (Cressman)**: A powerful approach that iteratively refines interpolated values based on surrounding observations.

Additionally, `InterpolateR` incorporates advanced machine learning-based interpolation techniques:

-   **RFmerge**: A Random Forest-based method designed to merge multiple precipitation products and ground observations to enhance interpolation accuracy.

-   **RFplus**: A novel approach that combines Random Forest with Quantile Mapping to refine interpolated estimates and minimize bias.

The package is continuously evolving, with new interpolation methods being added in future releases. Its highly optimized implementation ensures fast computations, making it an excellent choice for researchers and practitioners working with large spatial datasets.

# Installation

Install the development version from GitHub

```{r}
devtools::install_github("Jonnathan-Landi/InterpolateR")
```

# Getting Started

`InterpolateR` is designed to maintain a consistent data structure across all functions. The input datasets BD_Obs and BD_Coord must follow specific formats to ensure smooth functionality.

## Structure of BD_Obs

`BD_Obs` can be a `data.frame` or `data.table`, and its structure should follow this format:

| Date       | Est_1 | Est_2 | Est_3 | ... | Est_n |
|------------|-------|-------|-------|-----|-------|
| 2024-01-01 | 12.4  | 15.3  | 10.2  | ... | 18.1  |
| 2024-01-02 | 14.1  | 16.8  | 11.0  | ... | 19.5  |
| 2024-01-03 | 13.5  | 14.9  | 10.8  | ... | 16.4  |

Each column (except the first one) represents a different season (Est_1, Est_2, etc.), while the rows correspond to different time stamps (dates).

## Structure of BD_Coord

`BD_Coord` must have the following structure:

|       |        |         |      |
|-------|--------|---------|------|
| Cod   | X      | Y       | Z    |
| Est_1 | 720227 | 9680901 | 2668 |
| Est_2 | 732260 | 9682763 | 2527 |
| Est_3 | 714574 | 967587  | 2956 |
| ...   | ...    | ...     | ...  |
| Est_n | 734574 | 9683686 | 2622 |

-   X and Y correspond to the longitude and latitude of each station in UTM.

-   Z is the altitude at which the station is located.

-   Cod is the identifier that matches the column names in BD_Obs.

For correct operation, make sure that the first column of Data_Obs contains the date of each observation. For BD_Coords it is important to structure the data in this way, i.e. Cod, X,Y,Z. Also, make sure that the Cod you assign matches the corresponding column in BD_Obs.

By ensuring that BD_Obs and BD_Coord follow these structures, you can seamlessly use any interpolation method provided by InterpolateR.

# Available Interpolation Methods

## Inverse distance weighting (IDW)

Inverse Distance Weighting (IDW) is a deterministic interpolation method that assumes that values closer to an unknown point have a greater influence than those farther away. This method assigns weights to sample points so that their influence decreases as the distance to the unknown point increases.

IDW estimates the value at an unknown point using a weighted average of the known values, where the weights are inversely proportional to the distances between the prediction point and the known points. The basic formula, as defined by Shepard, is:

$$
\hat{Z}(s_0) = \frac{\sum_{i=1}^{N} w_i Z(s_i)}{\sum_{i=1}^{N} w_i}
$$

where:

-   $\hat{Z}(s_0)$ is the estimated value at the unknown point.
-   $Z(s_i)$ are the known values.
-   $w_i$ are the weights assigned to each known point.
-   $N$ is the total number of known points used.

The weights are calculated using:

$$
w_i = \frac{1}{d(s_0, s_i)^p}
$$

where:

-   $d(s_0, s_i)$ is the distance between the unknown point $s_0$ and the known point $s_i$.
-   $p$ is the power parameter that controls how the influence decreases with distance.

### Application of the method in R

```{r}
library(InterpolateR)
# Load data from on-site observations
data("BD_Obs", package = "InterpolateR")
data("BD_Coord", package = "InterpolateR")

# Load the study area where the interpolation is performed.
shapefile <- terra::vect(system.file("extdata/study_area.shp", package = "InterpolateR"))

# Perform the interpolation
Interpolated_data <- IDW(BD_Obs, BD_Coord, shapefile, grid_resolution = 5, p = 2, n_round = 1, training = 1, Rain_threshold = NULL, stat_validation = "M001", save_model = FALSE, name_save = NULL)

# Summary of model performance
print(Interpolated_data$Validation)
```

Since the algorithms are designed to work in the UTM system, make sure that grid_resolution is entered in Km.

## Cressman Objective Analysis Method (Cressman)

The Cressman objective analysis computes values at grid points $Z_{ij}^a$ (where $i$ and $j$ are the grid point indices for a 2D grid) as the weighted average of the difference between observed values $Z_k^o$ and background values interpolated to the observation locations $Z_k^b$ (i.e., $Z_k^o - Z_k^b$, called the observation increment) plus the background value at the grid point $Z_{ij}^b$.

The Cressman method is defined by the following equation:

$$
Z_{ij}^a = Z_{ij}^b + \frac{\sum_{k=1}^{n} w_k (Z_k^o - Z_k^b)}{\sum_{k=1}^{n} w_k}
$$

where:

-   $Z_{ij}^a$ is the analysis value at grid point $i,j$.

-   $Z_{ij}^b$ is the background value at grid point $i,j$.

-   $Z_k^o$ is the observed value at station $k$.

-   $Z_k^b$ is the background value interpolated to station $k$.

-   $w_k$ is the weight assigned to station $k$.

-   $n$ is the total number of stations used.

The weight $w_k$ is a relation of influence radius $R$ and the distance $r$ between the observation point and the grid point. The weight is defined as: $$
w_k = \frac{R^2 - r^2}{R^2 + r^2}
$$

where:

-   $r = \sqrt{(x_{ij} - x_k)^2 + (y_{ij} - y_k)^2}$ is the distance between the individual observation $k$ and grid point $(i, j)$.
-   $R$ is the influence radius.
-   Beyond the influence radius, the weight is set to zero. $R$ is therefore often referred to as the cut-off radius.

### Search Radius Considerations

The `search_radius` parameter defines the influence range for the Cressman interpolation method. It determines the maximum distance (in kilometers) within which observational data points contribute to the interpolated value at a given location. A larger radius results in smoother interpolated fields but may oversmooth local variations, while a smaller radius preserves finer details but may introduce noise.

The Cressman method typically applies an iterative approach, where the search radius is progressively reduced to refine the interpolation. Each iteration recalculates interpolated values with a smaller radius, allowing a better representation of small-scale features in the dataset.

### Application of the method in R

```{r}
# Load data from on-site observations
data("BD_Obs", package = "InterpolateR")
data("BD_Coord", package = "InterpolateR")

# Load the study area where the interpolation is performed.
shapefile <- terra::vect(system.file("extdata/study_area.shp", package = "InterpolateR"))

# Perform the interpolation
Interpolated_Cressman <- Cressman(BD_Obs, BD_Coord, shapefile, grid_resolution = 5, search_radius = c(20,10), training = 1,stat_validation = "M001", Rain_threshold = NULL, save_model = FALSE)

# Summary of model performance
print(Interpolated_Cressman$Validation)
```

### Important considerations

-   The `search_radius` should be defined as a numeric vector representing the influence range in kilometers (`km`) for each interpolation iteration. For example, setting: `search_radius = c(50, 20, 10)`.
-   Since the algorithm works in the UTM coordinate system, make sure that the `grid_resolution` is entered in km.

## RFplus

The RFplus package implements a novel spatial extrapolation and bias correction framework, integrating Random Forest (RF) and Quantile Mapping (QM) in a multi-stage process to improve the accuracy of satellite precipitation estimates. The methodology consists of three key stages:

1.  **Spatial Extrapolation of Precipitation:** The first stage employs a Random Forest model to extrapolate the spatial distribution of precipitation. The model is trained using in-situ measurements as the response variable and a diverse set of satellite precipitation products and environmental covariates as predictors. This approach enables the generation of an initial precipitation field that extends observed precipitation patterns across unmonitored regions with high spatial flexibility, allowing applications at different temporal scales (e.g., daily, monthly, or annual).

2.  **Residual Correction through a Secondary RF Model:** To enhance predictive accuracy, a second Random Forest model is trained to estimate residual errors from the initial predictions. The residuals are defined as the difference between observed and modeled precipitation values at station locations. By modeling these residuals as a function of the same covariates used in the first stage, systematic biases are identified and corrected iteratively. The corrected precipitation estimate is obtained by summing the residual predictions to the initial RF-based precipitation estimates, leading to a refined precipitation product with reduced bias and improved spatial coherence.

3.  **Bias Adjustment via Non-Parametric Quantile Mapping (QM):** In the third stage, a nonparametric quantile mapping (QM) is applied to adapt the distribution of each time series to the in situ observations of the nearest station. The QM correction will be applied to those pixels that meet the proximity criterion, which states that only pixels within a predefined radius of influence (e.g., ≤15 km) are QM corrected.

The RFplus package is designed to be highly adaptable and can be utilized across various satellite precipitation products and geographic regions. Although initially developed for precipitation bias correction, its methodology is applicable to other environmental variables such as temperature, wind speed, and soil moisture. This versatility makes RFplus a powerful tool for enhancing the accuracy of remote sensing-based estimations across diverse environmental conditions.

### Application of the method in R

```{r}
# Load the in-situ data and the coordinates
data("BD_Obs", package = "InterpolateR")
data("BD_Coord", package = "InterpolateR")

# Load the covariates
Covariates <- list(
  MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
  CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
  DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
   )

# 2. Apply de the model
Model_rfplus = RFplus(BD_Obs, BD_Coord, Covariates, n_round = 1, wet.day = 0.1,ntree = 2000, seed = 123, training = 1, stat_validation = c("M001"), Rain_threshold = NULL, method = "none", ratio = 15, save_model = FALSE, name_save = NULL)

# 3. Summary of model performance
print(Model_rfplus$Validation)

```

## RFmerge

RFmerge is a methodology developed by Baez-Villanueva et al. (2020) for the fusion of satellite precipitation datasets with ground-based observations, with the objective of improving the accuracy and spatial representativeness of the data. This package implements RFmerge using Random Forest as a machine learning technique to correct biases and adjust the distribution of satellite products to in situ measurements.

### Application of the method in R

```{r}
# Load the in-situ data and the coordinates
data("BD_Obs", package = "InterpolateR")
data("BD_Coord", package = "InterpolateR")

# Load the covariates
cov <- list(
  MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "InterpolateR")),
   CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "InterpolateR")),
   DEM = terra::rast(system.file("extdata/DEM.nc", package = "InterpolateR"))
   )

# Apply the RFmerge
model_RFmerge = RFmerge(BD_Obs, BD_Coord, cov, mask = NULL, n_round = 1, ntree = 2000, seed = 123,  training = 1, stat_validation = c("M001"), Rain_threshold = NULL, save_model = FALSE, name_save = NULL)

# 3. Summary of model performance
print(model_RFmerge$Validation)

```

# Validation methods

In all interpolation methods it is possible to perform a validation of the results at point to pixel level. To perform this validation process, several parameters such as training, stat_validation and Rain_threshold must be correctly configured. These parameters allow to define how the data will be divided for training and validation.

## Key parameters for validation

1.  **training:** This parameter defines the proportion of stations to be used to train the model. The value must be between 0 and 1. For example, if you want to use 80% of the stations for training and the remaining 20% for validation, the training parameter should be set to 0.8. This means that 80% of the stations will be used to train the interpolation model, while the remaining 20% will be used to validate it. If the training value is 1, all stations will be used for training and no validation will be performed.

2.  **stat_validation:** This parameter allows to perform a manual validation, i.e. to select the specific stations to be used for validation. Instead of the selection of the stations being random, you can explicitly set the stations you want to use for this purpose. For example, if you want to use the stations “Est_1”, “Est_5” and “Est_7” for the validation, you must set stat_validation as follows: stat_validation = c(“Est_1”, “Est_5”, “Est_7”)

3.  **Rain_threshold:** The Rain_threshold parameter is used exclusively when performing precipitation validation. This parameter is critical for calculating categorical metrics that help evaluate the performance of the interpolation model, such as: Critical Success Index (CSI); Probability of Detection (POD); False Alarm Rate (FAR); Success Ratio (SR), etc.

    These metrics are essential for assessing the quality of precipitation prediction by comparing estimates with actual precipitation observations.

    **Parameter format**

    When using the Rain_Threshold parameter, it should be entered in named list format, where each element represents a precipitation category. The element name should be the category name, and the values associated with each category should be a numeric vector with two values: the lower limit and the upper limit of the category. For example:

    ```{r}
    Rain_threshold = list(
      no_rain = c(0, 1),
      light_rain = c(1, 5),
      moderate_rain = c(5, 20),
      heavy_rain = c(20, 40),
      violent_rain = c(40, Inf)
    )
    ```

    By entering these values, InterpolateR will classify the precipitation values into these categories and can calculate the corresponding validation metrics to evaluate model performance.

## Goodness-of-fit and categorical metrics calculated by InterpolateR

InterpolateR calculates two main types of metrics to assess the accuracy of interpolation predictions: Goodness-of-fit metrics and categorical metrics.

### Goodness-of-fit metrics

These metrics evaluate the overall agreement between predicted and observed values. Metrics calculated by InterpolateR include:

-   Mean absolute error (MAE)
-   Spearman's coefficient (CC)
-   Root mean square error (RMSE)
-   Kling-Gupta Efficiency (KGE)
-   Nash-Sutcliffe Efficiency (NSE)
-   Percentage bias (PBIAS)

### Categorical metrics

These metrics are used to evaluate the model's ability to predict extreme precipitation events by ranking values according to their intensity. Categorical metrics calculated by InterpolateR include:

-   Critical Success Index (CSI).
-   Probability of detection (POD)
-   False Alarm Rate (FAR)
-   Success Rate (SR)
-   Hit Bias (HB)
-   Heidke Skill Score (HSS)
-   Hanssen-Kuipers Discriminant (HK)
-   Equal Threat Score (ETS)
-   Gilbert Skill Score

These metrics allow for a more detailed and specific assessment of the quality of precipitation predictions based on intensity categories.
