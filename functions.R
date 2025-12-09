# Dog Fouling & Toxocara ABM (Ben McAnoy)
# Script 1/3. Requires: tidyverse, ggplot2, spatstat

# ---- Functions ----

library(MASS) # for mvrnorm (multivariate normal)

# Generate correlated attributes for N agents
# rho: Correlation coefficient (-1 to 1).
# uses negative (e.g., -0.5) if high cleanup = low infection risk.
draw_correlated_attributes <- function(n, cleanup_dist, infection_prob, rho = -0.5) {

  # generate correlated standard normal variables
  # Sigma is the covariance matrix
  Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
  mu <- c(0, 0)

  # Z is a matrix with n rows and 2 columns
  Z <- mvrnorm(n, mu = mu, Sigma = Sigma)

  # convert to uniform Distributions [0, 1] using CDF
  U <- pnorm(Z)

  # map U[,1] to Cleanup Propensity
  # assumes 'cleanup_dist' is the vector of probabilities
  # sort then pick values based on the rank of U[,1]
  sorted_cleanup <- sort(cleanup_dist)
  #map 0-1 to index 1-Length
  idx <- ceiling(U[, 1] * length(sorted_cleanup))
  #ensures idx is within bounds (1 to length)
  idx <- pmax(1, pmin(length(sorted_cleanup), idx))
  assigned_cleanup <- sorted_cleanup[idx]

  # U[,2] to Infection Status
  # If rho is NEGATIVE (e.g. -0.5), High Cleanup (High U1) pairs with Low U2.
  # High Cleanup = healthy (Not Infected).
  # So "Low U2" = "healthy".
  # Therefore, "High U2" = "infected".#

  # Threshold is the top X% of the distribution
  threshold <- 1 - infection_prob
  assigned_infected <- U[, 2] > threshold

  return(tibble(cleanup_propensity = assigned_cleanup, infected = assigned_infected))
}

# (Keep your existing fouling_weight_surface and sample_location functions here)
# draw cleanup propensity from distribution
draw_cleanup_propensity <- function(n) {
  sample(cleanup_probs, n, replace = TRUE)
}

# fouling probability surface for a given entry point
# returns matrix of weights (not normalised)

fouling_weight_surface <- function(entry_x, entry_y, on_lead = TRUE) {
  # distance from entry and path already precomputed
  # entrance-based weight (exponential decay)
  w_ent <- exp(- entrance_lambda * dist_to_entrance)
  # path-based weight: Weibull for distance from path
  if (on_lead) {
    w_path <- pweibull(dist_to_path, shape = weibull_shape, scale = weibull_scale)
  } else {
    # flatten: increase scale to broaden distribution
    w_path <- pweibull(dist_to_path, shape = weibull_shape,
      scale = weibull_scale * flatten_factor_offlead
    )
  }
  # multiplicative combination
  w <- w_ent * w_path
  if (all(w == 0)) w <- matrix(1, GRID_SIZE, GRID_SIZE)
  return(w)
}

# sample fouling location matrix indices given weight surface
sample_location_from_weights <- function(weight_mat) {
  w <- as.numeric(weight_mat)
  w[w < 0] <- 0
  if (all(w == 0)) w <- rep(1, length(w))
  p <- w / sum(w)
  idx <- sample(seq_along(p), size = 1, prob = p)
  # convert index to x,y
  x <- ((idx - 1) %% GRID_SIZE) + 1
  y <- floor((idx - 1) / GRID_SIZE) + 1
  list(x = x, y = y)
}

# proximity to a bin (min Euclidean distance)
# dist_to_nearest_bin <- function(x,y) {
# min(sqrt((bins$x - x)^2 + (bins$y - y)^2))
# }

# check if within bin radius
# is_within_bin_radius <- function(x,y, radius = bin_radius) {
# dist_to_nearest_bin(x,y) <= radius
# }




