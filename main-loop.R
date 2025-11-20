# Dog Fouling & Toxocara ABM (Ben McAnoy)
# Single-file script. Requires: tidyverse, ggplot2, spatstat 

library(tidyverse)
library(ggplot2)
library(spatstat)

set.seed(123) # pseudo-random start

# ---- Input Parameters ---- 

GRID_SIZE <- 100 # 100 x 100 grid
SIM_DAYS <- 90  # simulation length
N_daily <- 200 # number of dog-owner pairs entering per day
prop_infected <- 0.15 # proportion of dogs shedding eggs
degradation_time <- 7 # faeces visual degradation time (days)
burn_in <- 10 # days to run before collecting (0 = none)

# Fouling probability parameters
entrance_lambda <- 0.03 # exponential decay rate with distance from entrance
weibull_shape <- 2  # weibull shape for distance from path (on-lead)
weibull_scale <- 3  # weibull scale
flatten_factor_offlead <- 3 # multiply scale to flatten distribution off-lead

# Exponential calculation
lambda <- 2
exp_sim <- rexp(n = 1000, rate = lambda)

u <- runif(1000)
exp_inv <- -log(u)/lambda
par(mar =c(4,4,1,1))
x <- seq(0, max(exp_sim), len = 500)
hist(exp_inv, probability = T, ylim= c(0, lambda), main = "", xlab = "t", ylab = "density")
lines(x, dexp(x, rate = lambda))

# Weibull calculation
w <- rweibull(n = 1000, shape = 1.5, scale = 2)

hist(w,probability = TRUE, ylim = c(0, 1.5/2), main = "", xlab = "t", ylab = "density")
x <- seq(0, max(w), length.out = 500)
lines(x, dweibull(x, shape = 1.5, scale = 2))



# Cleanup propensity distribution (bimodal)
cleanup_probs <- c(rep(0.9, 0.8 * 1000), rep(0.3, 0.2 * 1000)) # draw from this

# Bins effect radius
# bin_radius <- 5 # cells within which cleanup propensity increases

# Movement parameters
walk_duration_mean <- 30 # number of steps (cells) per visit (not used strictly; we'll use path length)
on_lead_prob <- 0.6 # fraction of dogs on lead (drawn per dog)
# offlead_boldness_mean <- 0.3 # when off-lead, probability to deviate at each step

# Visualization options
plot_final_heatmaps <- TRUE

# ---- Park Layout ----

# grid coordinates
coords <- expand.grid(x = 1:GRID_SIZE, y = 1:GRID_SIZE)

# create a simple set of entrances (on edge): 4 entrances roughly in middle of each side
entrances <- tibble(
  id = 1:4,
  x = c(1, GRID_SIZE, round(GRID_SIZE/2), round(GRID_SIZE/2)),
  y = c(round(GRID_SIZE/2), round(GRID_SIZE/2), 1, GRID_SIZE)
)

# create a few linear paths (represented as lists of cells)
# defines two crossing paths: horizontal and vertical through centre, plus a winding path
center <- round(GRID_SIZE / 2)
path_cells <- tibble()
path_cells <- bind_rows(path_cells,
                        tibble(x = 1:GRID_SIZE, y = center),  # horizontal main path
                        tibble(x = center, y = 1:GRID_SIZE))  # vertical main path

# winding path 
t <- seq(0, 1, length.out = GRID_SIZE)
wx <- round((GRID_SIZE/4) + (GRID_SIZE/2) * (t + 0.2 * sin(4 * pi * t)))
wy <- round((GRID_SIZE/4) + (GRID_SIZE/2) * (rev(t) + 0.2 * cos(3 * pi * t)))
path_cells <- bind_rows(path_cells, tibble(x = pmin(pmax(wx,1),GRID_SIZE), y = pmin(pmax(wy,1),GRID_SIZE)))

path_cells <- distinct(path_cells)

# Plot out the park 
plot(path_cells$x, path_cells$y,type = "l",asp = 1, xlab = "x", ylab = "y", main = "Park path",col = "blue", lwd = 2)

points(path_cells$x, path_cells$y, pch = 16, cex = 0.3) # individual points 


# bins/disposal points (a few near entrances and central)
bins <- tibble(id = 1:5,
               x = c(2, GRID_SIZE-1, center-3, center+4, center),
               y = c(center, center, 2, GRID_SIZE-1, center+6))

# mark path distance raster (distance from nearest path cell)
# precompute distance from path for every cell
compute_euclid_dist <- function(x1,y1,x2,y2) sqrt((x1-x2)^2 + (y1-y2)^2)
path_matrix <- matrix(FALSE, nrow = GRID_SIZE, ncol = GRID_SIZE)
for(i in 1:nrow(path_cells)) path_matrix[path_cells$x[i], path_cells$y[i]] <- TRUE

dist_to_path <- matrix(NA, GRID_SIZE, GRID_SIZE)
for(x in 1:GRID_SIZE){
  for(y in 1:GRID_SIZE){
    pts <- path_cells
    # compute nearest distance to any path cell (euclidean)
    dist_to_path[x,y] <- min((pts$x - x)^2 + (pts$y - y)^2) ^ 0.5
  }
}

# precompute distance from nearest entrance
dist_to_entrance <- matrix(NA, GRID_SIZE, GRID_SIZE)
for(x in 1:GRID_SIZE){
  for(y in 1:GRID_SIZE){
    dist_to_entrance[x,y] <- min(((entrances$x - x)^2 + (entrances$y - y)^2) ^ 0.5)
  }
}

# contamination layer i.e. whether a cell has ever received infected eggs (biary)
contamination_layer <- matrix(FALSE, GRID_SIZE, GRID_SIZE)

# ---- Object states (Faeces) ----

# data.frame of faeces deposits (each row = one deposit). Will grow/shrink.
# columns: id, x, y, age (days), contains_eggs (logical), degradation_time, persistent_contamination (logical)
faeces <- tibble(
  id = integer(),
  x = integer(),
  y = integer(),
  age = integer(),
  contains_eggs = logical(),
  degradation_time = integer(),
  persistent_contamination = logical()
)

# daily summary storage
daily_summary <- tibble(day = integer(),
                        total_deposited = integer(),
                        total_removed = integer(),
                        total_visible = integer(),
                        contaminated_cells = integer(),
                        mean_dist_from_entrance = numeric(),
                        var_dist_from_entrance = numeric(),
                        mean_dist_from_path = numeric(),
                        var_dist_from_path = numeric(),
                        cleanup_success_rate = numeric())

# unique id generator for faeces
next_faeces_id <- 1L

# ---- Functions ----

# draw cleanup propensity from distribution
draw_cleanup_propensity <- function(n) {
  sample(cleanup_probs, n, replace = TRUE)
}

# fouling probability surface for a given entry point
# returns matrix of weights (not normalised)

fouling_weight_surface <- function(entry_x, entry_y, on_lead = TRUE) {
  # distance from entry and path already precomputed
  # entrance-based weight (exponential decay)
  w_ent <- exp(-entrance_lambda * dist_to_entrance)
  # path-based weight: Weibull for distance from path
  if (on_lead) {
    w_path <- pweibull(dist_to_path, shape = weibull_shape, scale = weibull_scale)
  } else {
    # flatten: increase scale to broaden distribution
    w_path <- pweibull(dist_to_path, shape = weibull_shape, scale = weibull_scale * flatten_factor_offlead)
  }
  # multiplicative combination
  w <- w_ent * w_path
  if(all(w == 0)) w <- matrix(1, GRID_SIZE, GRID_SIZE)
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

# ---- Main sim loop ----

for(day in 1:SIM_DAYS) {
  # 1. Agent entry: generate N_daily dog-owner pairs
  dogs <- tibble(
    dog_id = 1:N_daily,
    entry_id = sample(entrances$id, N_daily, replace = TRUE),
    entry_x = entrances$x[match(entry_id, entrances$id)],
    entry_y = entrances$y[match(entry_id, entrances$id)],
    on_lead = runif(N_daily) < on_lead_prob,
    cleanup_propensity = draw_cleanup_propensity(N_daily),
    fouling_urge = runif(N_daily) < 0.5,  # could be parameterized; here 50% have urge
    infected = runif(N_daily) < prop_infected
  )
  
  # movement & fouling phase: for each dog, simulate a simple path-following walk
  # For simplicity each dog will "simulate" a sequence of steps along a path towards the opposite edge:
  deposited_today <- tibble()  # collects today's new faeces
  removed_today <- 0L
  visible_today <- 0L
  dist_from_entrance_list <- numeric()
  dist_from_path_list <- numeric()
  cleanup_successes <- 0L
  fouling_events <- 0L
  
  for(i in seq_len(nrow(dogs))) {
    d <- dogs[i, ]
    # choose a route: simple straight to central or follow one of the paths
    # make dogs follow the nearest big path if close enough, otherwise go straight toward park center then exit
    # once fouling location chosen per visit if fouling_urge==TRUE
    if(!d$fouling_urge) next
    
    # compute fouling weight surface for this dog's lead state and entry
    weight_surf <- fouling_weight_surface(d$entry_x, d$entry_y, on_lead = d$on_lead)
    
    # random draw of location
    loc <- sample_location_from_weights(weight_surf)
    fx <- loc$x; fy <- loc$y
    fouling_events <- fouling_events + 1
    
    # record distances for stats
    dist_from_entrance <- sqrt((fx - d$entry_x)^2 + (fy - d$entry_y)^2)
    dist_from_path <- dist_to_path[fx, fy]
    dist_from_entrance_list <- c(dist_from_entrance_list, dist_from_entrance)
    dist_from_path_list <- c(dist_from_path_list, dist_from_path)
    
    # Create faeces deposit (owner may clean up)
    # modify cleanup propensity if near bin
    cleanup_p <- d$cleanup_propensity
    if (is_within_bin_radius(fx, fy, bin_radius)) cleanup_p <- min(1, cleanup_p + 0.15) # bins increase cleanup prob
    
    cleaned <- runif(1) < cleanup_p
    if(cleaned) {
      removed_today <- removed_today + 1L
      cleanup_successes <- cleanup_successes + 1L
    } else {
      # create deposit
      new_row <- tibble(
        id = next_faeces_id,
        x = fx,
        y = fy,
        age = 0L,
        contains_eggs = d$infected,
        degradation_time = degradation_time,
        persistent_contamination = FALSE
      )
      faeces <- bind_rows(faeces, new_row)
      next_faeces_id <- next_faeces_id + 1L
      visible_today <- visible_today + 1L
      
      # if infected, mark contamination layer persistent
      if(d$infected) {
        contamination_layer[fx, fy] <- TRUE
        # mark that this deposit left persistent contamination (depending on model)
        faeces$persistent_contamination[nrow(faeces)] <- TRUE
      }
    }
  } # end dogs loop
  
  # 3. Cleanup Phase: (we already modelled immediate cleanup at time of fouling above)
  # 4. Degradation/Contamination Update:
  if(nrow(faeces) > 0) {
    faeces <- faeces %>% 
      mutate(age = age + 1L)
    # remove those exceeding degradation_time (visual disappears)
    # But persistent_contamination remains in contamination_layer if contains_eggs
    to_remove <- which(faeces$age > faeces$degradation_time)
    if(length(to_remove) > 0) {
      # For removed visual deposits we keep contamination layer (already set on deposit creation)
      faeces <- faeces[-to_remove, ]
    }
  }
  
  # 5. Data collection at end of day
  total_deposited_today <- fouling_events
  total_removed_today <- removed_today
  total_visible_now <- nrow(faeces)
  contaminated_cells_now <- sum(contamination_layer)
  mean_e <- if(length(dist_from_entrance_list)>0) mean(dist_from_entrance_list) else NA
  var_e <- if(length(dist_from_entrance_list)>1) var(dist_from_entrance_list) else NA
  mean_p <- if(length(dist_from_path_list)>0) mean(dist_from_path_list) else NA
  var_p <- if(length(dist_from_path_list)>1) var(dist_from_path_list) else NA
  cleanup_rate <- if(total_deposited_today>0) cleanup_successes/total_deposited_today else NA
  
  daily_summary <- bind_rows(daily_summary,
                             tibble(day = day,
                                    total_deposited = total_deposited_today,
                                    total_removed = total_removed_today,
                                    total_visible = total_visible_now,
                                    contaminated_cells = contaminated_cells_now,
                                    mean_dist_from_entrance = mean_e,
                                    var_dist_from_entrance = var_e,
                                    mean_dist_from_path = mean_p,
                                    var_dist_from_path = var_p,
                                    cleanup_success_rate = cleanup_rate))
  
  # optional: simple progress
  if(day %% 10 == 0) message("Completed day ", day, " / ", SIM_DAYS)
}

# -----------------------
# OUTPUTS & PLOTS
# -----------------------
# Print summary
print(glimpse(daily_summary))
print(tail(daily_summary, 10))

# Final heatmap of visible faeces frequency & contamination
# Create dataframe for plotting frequency of visible faeces by cell
visible_counts <- faeces %>%
  group_by(x,y) %>%
  summarise(visible = n(), .groups = "drop")

# contamination layer to dataframe
cont_df <- as_tibble(which(contamination_layer, arr.ind = TRUE)) %>%
  rename(x = row, y = col) %>%
  mutate(contaminated = TRUE)

# for plotting a heatmap of historical contamination counts, compute how many times contaminated (binary so it's 1)
cont_plot_df <- expand.grid(x = 1:GRID_SIZE, y = 1:GRID_SIZE) %>%
  as_tibble() %>%
  left_join(cont_df, by = c("x","y")) %>%
  mutate(contaminated = if_else(is.na(contaminated), FALSE, contaminated))

# visible heatmap
if(plot_final_heatmaps) {
  p1 <- cont_plot_df %>%
    ggplot(aes(x = x, y = y, fill = contaminated)) +
    geom_raster() +
    scale_fill_manual(values = c("white", "red")) +
    coord_fixed() +
    labs(title = "Persistent contamination (cells that have ever received infected faeces)") +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank())
  print(p1)
  
  if(nrow(visible_counts) > 0) {
    p2 <- visible_counts %>%
      ggplot(aes(x = x, y = y, fill = visible)) +
      geom_raster() +
      coord_fixed() +
      labs(title = "Visible faeces counts (present at simulation end)") +
      theme_minimal() +
      theme(axis.text = element_blank(), axis.title = element_blank())
    print(p2)
  } else {
    message("No visible faeces at end of simulation to plot.")
  }
}

# Example summary plots: daily counts
p3 <- daily_summary %>%
  pivot_longer(cols = c(total_deposited, total_removed, total_visible), names_to = "metric") %>%
  ggplot(aes(x = day, y = value, color = metric)) +
  geom_line() +
  labs(title = "Daily fouling/removal/visible counts", y = "count") +
  theme_minimal()
print(p3)

p4 <- daily_summary %>%
  ggplot(aes(x = day, y = contaminated_cells)) +
  geom_line() +
  labs(title = "Number of contaminated cells over time", y = "contaminated cells") +
  theme_minimal()
print(p4)

# Save results to list for later analysis
results <- list(
  params = list(GRID_SIZE = GRID_SIZE, SIM_DAYS = SIM_DAYS, N_daily = N_daily,
                prop_infected = prop_infected, degradation_time = degradation_time),
  faeces = faeces,
  contamination_layer = contamination_layer,
  daily_summary = daily_summary,
  park = list(entrances = entrances, paths = path_cells, bins = bins)
)

# Optionally save to disk
saveRDS(results, file = "dog_fouling_abm_results.rds")
message("Simulation complete. Results stored in `results` object.")

