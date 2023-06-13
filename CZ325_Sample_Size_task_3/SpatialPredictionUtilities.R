
# CZ325 Deliverable 3a: Evaluation of sample size to assess managed-area-level trends
# Utilities for spatial prediction
# this prediction function is based on the spatial_mean_RE.stan model
predict_mean <- function(global_mean, alpha, rho, obs_gp_function,
                       obs_distances, pred_distances,
                       obs_to_pred_distances){
  
  # Estimated covariance kernel for the fitted points
  K_obs <- alpha ^ 2 * exp(-0.5 * ((obs_distances / rho) ^ 2)) +
    diag(1e-4, dim(obs_distances)[1])
  
  # Cross-covariance between prediction points and fitted points
  K_new <- alpha ^ 2 * exp(-0.5 * ((obs_to_pred_distances / rho) ^ 2))
  
  # Estimated covariance kernel for prediction points
  K_star <- alpha ^ 2 * exp(-0.5 * ((pred_distances / rho) ^ 2)) +
    diag(1e-4, dim(pred_distances)[1])
  
  # Estimated mean for prediction points
  t(K_new) %*% solve(K_obs, obs_gp_function + global_mean) +
    # One realisation of uncertainty for prediction points
    MASS::mvrnorm(1, mu = rep(0, dim(pred_distances)[1]),
                  Sigma = K_star - t(K_new) %*% solve(K_obs, K_new))
}

predict_covariate_mean <- function(site_mean, alpha, rho, obs_gp_function,
                         obs_distances, pred_distances,
                         obs_to_pred_distances){
  
  # Estimated covariance kernel for the fitted points
  K_obs <- alpha ^ 2 * exp(-0.5 * ((obs_distances / rho) ^ 2)) +
    diag(1e-4, dim(obs_distances)[1])
  
  # Cross-covariance between prediction points and fitted points
  K_new <- alpha ^ 2 * exp(-0.5 * ((obs_to_pred_distances / rho) ^ 2))
  
  # Estimated covariance kernel for prediction points
  K_star <- alpha ^ 2 * exp(-0.5 * ((pred_distances / rho) ^ 2)) +
    diag(1e-4, dim(pred_distances)[1])
  
  # Estimated mean for prediction points
  t(K_new) %*% solve(K_obs, site_mean) +
    # One realisation of uncertainty for prediction points
    MASS::mvrnorm(1, mu = rep(0, dim(pred_distances)[1]),
                  Sigma = K_star - t(K_new) %*% solve(K_obs, K_new))
}


my_calc_point_distances = function(coordinates, coordinates2){
  if(missing(coordinates2)){
    distances <- pointDistance(coordinates, lonlat = FALSE)
    distances[upper.tri(distances)] = t(distances)[upper.tri(distances)]
  } else {
    distances <- pointDistance(coordinates, coordinates2, lonlat = FALSE)
  }
  
  # Convert meters distance to kilometers
  #distances <- distances / 1000
  
  return(distances)
}

