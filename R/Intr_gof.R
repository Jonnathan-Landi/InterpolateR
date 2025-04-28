# This script contains the codes for the model evaluation metrics
# Note that these are internal functions that will only be used by the exported functions.
##############################################################################
#                             Goodness-of-fit metrics                        #
##############################################################################
# Mean Absolute Error (MAE)
.mae <- function(sim, obs) {
  valid <- stats::complete.cases(cbind(obs, sim))
  obs <- obs[valid]
  sim <- sim[valid]
  if (length(obs) == 0) return(NA_real_)
  mean(abs(sim - obs))
}

# Mean Squared Error (MSE)
.mse <- function(sim, obs) {
  valid <- stats::complete.cases(cbind(obs, sim))
  obs <- obs[valid]
  sim <- sim[valid]
  if (length(obs) == 0) return(NA_real_)
  mean((sim - obs)^2)
}

# Spearman's Rank Correlation
.rspearman <- function(sim, obs) {
  valid <- stats::complete.cases(cbind(obs, sim))
  obs <- obs[valid]
  sim <- sim[valid]
  if (length(obs) < 2) return(NA_real_)
  stats::cor(obs, sim, method = "spearman")
}

# Root Mean Squared Error (RMSE)
.rmse <- function(sim, obs) {
  valid <- stats::complete.cases(cbind(obs, sim))
  obs <- obs[valid]
  sim <- sim[valid]
  if (length(obs) == 0) return(NA_real_)
  sqrt(mean((sim - obs)^2))
}

# Kling-Gupta Efficiency (KGE)
.kge <- function(sim, obs) {
  valid <- stats::complete.cases(cbind(obs, sim))
  obs <- obs[valid]
  sim <- sim[valid]
  n <- length(obs)
  if (n < 2) return(NA_real_)
  
  mean_obs <- mean(obs)
  if (mean_obs == 0) return(NA_real_)
  
  sd_obs <- stats::sd(obs)
  if (sd_obs == 0) return(NA_real_)
  
  r <- stats::cor(obs, sim, method = "pearson")
  if (is.na(r)) return(NA_real_)
  
  alpha <- stats::sd(sim)/sd_obs
  beta <- mean(sim)/mean_obs
  
  1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
}

# Nash-Sutcliffe Efficiency (NSE)
.nse <- function(sim, obs) {
  valid <- stats::complete.cases(cbind(obs, sim))
  obs <- obs[valid]
  sim <- sim[valid]
  if (length(obs) == 0) return(NA_real_)
  
  mean_obs <- mean(obs)
  denominator <- sum((obs - mean_obs)^2)
  if (denominator == 0) return(NA_real_)
  
  numerator <- sum((sim - obs)^2)
  1 - (numerator/denominator)
}

# Percent Bias (PBIAS)
.pbias <- function(sim, obs) {
  valid <- stats::complete.cases(cbind(obs, sim))
  obs <- obs[valid]
  sim <- sim[valid]
  if (length(obs) == 0) return(NA_real_)
  
  sum_obs <- sum(obs)
  if (sum_obs == 0) return(NA_real_)
  
  100 * sum(sim - obs)/sum_obs
}