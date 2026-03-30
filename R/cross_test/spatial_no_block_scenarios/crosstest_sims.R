# Purpose: To simulate data for cross-testing models (conditioned simulation w/o selectivity blocks)
# Notes: Runs an MSE loop for one year, but anything prior to n_years can essentially
# be used as a self / cross test without any feedback. A bit slower, but easy enough for our purposes here.
# Creator: Matthew LH. Cheng
# Date: 12/9/25

# Setup -------------------------------------------------------------------
library(here)
library(SPoRC)
library(tidyverse)
library(future.apply)
library(progressr)
library(furrr)

# Load in functions
source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "single_region_em.R"))

# Read in model output
om_values <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_NoBlock_model_results.RDS"))

# Condition OM ------------------------------------------------------------

# Extract model components from sablefish model
data <- om_values$data
rep <- om_values$rep
parameters <- om_values$parameters
mapping <- om_values$mapping
sd_rep <- om_values$sd_rep

# Define operating model parameters
closed_loop_yrs <- 1      # Years to project forward
n_years <- length(data$years)  # number of years
burnin_years <- 1:n_years  # Historical conditioning period
n_sims <- 100              # Number of replicate simulations
n_regions <- 5             # number of regions
n_fish_fleets <- 2         # number of fishery fleets
n_srv_fleets <- 3          # number of survey fleets

# get fleet allocation by region
fleet_allocation <- array(NA, dim = c(n_regions, n_fish_fleets))
fleet_allocation[,1] <- c(0.5,0.75,0.8,0.8,0.95) # from fmp
fleet_allocation[,2] <- 1 - fleet_allocation[,1] # trawl gear allocation

# Condition closed-loop simulations, mean recruitment
sim_list_lowsamp <- condition_closed_loop_simulations(
  closed_loop_yrs = closed_loop_yrs,
  n_sims = n_sims,
  data = data,
  parameters = parameters,
  mapping = mapping,
  sd_rep = sd_rep,
  rep = rep,
  random = NULL,
  recruitment_opt = 'resample_from_input',
  # setup variances
  ISS_FishAgeComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_FishLenComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_SrvAgeComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ISS_SrvLenComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
  ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets)),
  n_tags_rel_input = rep(2e3, nrow(data$tag_release_indicator))
)

# Condition closed-loop simulations, mean recruitment with time-varying movement
# Modify report file to incorporate time-varying movement
tv_moverep <- rep
tv_moverep$Movement[,,20:39,,] <- tv_moverep$Movement[5:1,5:1,1:20,,,drop = FALSE] # reverse movement (turn 1 to 5 rates into 5 to 1 rates)
tv_moverep$Rec[,]

sim_list_lowsamp_tvmove <- condition_closed_loop_simulations(
  closed_loop_yrs = closed_loop_yrs,
  n_sims = n_sims,
  data = data,
  parameters = parameters,
  mapping = mapping,
  sd_rep = sd_rep,
  rep = tv_moverep,
  random = NULL,
  recruitment_opt = 'resample_from_input',
  # setup variances
  ISS_FishAgeComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_FishLenComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_SrvAgeComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ISS_SrvLenComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
  ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets)),
  n_tags_rel_input = rep(2e3, nrow(data$tag_release_indicator))
)

# Condition closed-loop simulations, GOA recruitment crash, with no movement
dt_rep <- rep
dt_rep$Movement[,,,,] <- diag(1, data$n_regions) # no movement
dt_rep$Rec[3:5,33:65] <- dt_rep$Rec[3:5,33:65] * 0.1

sim_list_dt <- condition_closed_loop_simulations(
  closed_loop_yrs = closed_loop_yrs,
  n_sims = n_sims,
  data = data,
  parameters = parameters,
  mapping = mapping,
  sd_rep = sd_rep,
  rep = dt_rep,
  random = NULL,
  recruitment_opt = 'resample_from_input',
  # setup variances
  ISS_FishAgeComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_FishLenComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_SrvAgeComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ISS_SrvLenComps = array(30, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
  ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets)),
  n_tags_rel_input = rep(2e3, nrow(data$tag_release_indicator))
)


# Condition closed-loop simulations, mean recruitment (high sample sizes)
# To simulate more tags or resimulate tags with new sample sizes
data$tag_release_indicator <- expand.grid(
  regions = 1:5, tag_yrs = 1:(n_years + closed_loop_yrs)
)

sim_list_highsamp <- condition_closed_loop_simulations(
  closed_loop_yrs = closed_loop_yrs,
  n_sims = n_sims,
  data = data,
  parameters = parameters,
  mapping = mapping,
  sd_rep = sd_rep,
  rep = rep,
  random = NULL,
  recruitment_opt = 'resample_from_input',
  # setup variances
  ISS_FishAgeComps = array(200, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_FishLenComps = array(200, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_SrvAgeComps = array(200, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ISS_SrvLenComps = array(200, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
  ObsSrvIdx_SE = array(0.1, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets)),
  n_tags_rel_input = rep(2e3, nrow(data$tag_release_indicator))
)


# Run MSEs to get "data" ----------------------------------------------------------------
# Single-region low sample size
sim_env_lowsamp <- Setup_sim_env(sim_list = sim_list_lowsamp)
sim_env_lowsamp <- add_aggregated_obj_to_simenv(sim_env = sim_env_lowsamp)
sim_env_lowsamp <- run_single_rg_closedloop_parallel(sim_env = sim_env_lowsamp, n_sims = n_sims,
                                                     fleet_allocation = fleet_allocation,
                                                     lls_design_type = "historical", n_cores = 8)
saveRDS(sim_env_lowsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "spt_rand_OM_lowsamp.RDS"))

# Single-region low sample size (tv movement)
sim_env_lowsamp_tvmove <- Setup_sim_env(sim_list = sim_list_lowsamp_tvmove)
sim_env_lowsamp_tvmove <- add_aggregated_obj_to_simenv(sim_env = sim_env_lowsamp_tvmove)
sim_env_lowsamp_tvmove <- run_single_rg_closedloop_parallel(sim_env = sim_env_lowsamp_tvmove, n_sims = n_sims,
                                                     fleet_allocation = fleet_allocation,
                                                     lls_design_type = "historical", n_cores = 8)
saveRDS(sim_env_lowsamp_tvmove, here("outputs", "cross_test", "spatial_no_block_scenarios", "spt_rand_OM_lowsamp_tvmove_sens.RDS"))


# Divergent trends (use single-region model to initialize)
sim_env_dt <- Setup_sim_env(sim_list = sim_list_dt)
sim_env_dt <- add_aggregated_obj_to_simenv(sim_env = sim_env_dt)
sim_env_dt <- run_single_rg_closedloop_parallel(sim_env = sim_env_dt, n_sims = n_sims,
                                                fleet_allocation = fleet_allocation,
                                                lls_design_type = "historical", n_cores = 15)
saveRDS(sim_env_dt, here("outputs", "cross_test", "spatial_no_block_scenarios", "spt_rand_OM_lowsamp_dt.RDS"))


# Single-region high sample size
sim_env_highsamp <- Setup_sim_env(sim_list = sim_list_highsamp)
sim_env_highsamp <- add_aggregated_obj_to_simenv(sim_env = sim_env_highsamp)
sim_env_highsamp <- run_single_rg_closedloop_parallel(sim_env = sim_env_highsamp, n_sims = n_sims,
                                                      fleet_allocation = fleet_allocation,
                                                      lls_design_type = "historical", n_cores = 8)
saveRDS(sim_env_highsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "spt_rand_OM_highsamp.RDS"))
