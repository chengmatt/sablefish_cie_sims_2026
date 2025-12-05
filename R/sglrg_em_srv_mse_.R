# Purpose: To run a five region OM, with a single region EM
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 11/26/25

# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)
library(furrr)
library(progressr)

# Read in model output
om_values <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_model_results.RDS"))

source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "single_region_em.R"))

# Condition OM ------------------------------------------------------------

# Extract model components from sablefish model
data <- om_values$data
rep <- om_values$rep
parameters <- om_values$parameters
mapping <- om_values$mapping
sd_rep <- om_values$sd_rep

# Define operating model parameters
closed_loop_yrs <- 30      # Years to project forward
n_years <- length(data$years)  # number of years
burnin_years <- 1:n_years  # Historical conditioning period
n_sims <- 100              # Number of replicate simulations
assess_freq <- 1           # Assessment frequency
data_yr_freq <- 1          # Data collection frequency
n_regions <- 5             # number of regions
n_fish_fleets <- 2         # number of fishery fleets
n_srv_fleets <- 3          # number of survey fleets

# get fleet allocation by region
fleet_allocation <- array(NA, dim = c(n_regions, n_fish_fleets))
fleet_allocation[,1] <- c(0.5,0.75,0.8,0.8,0.95) # from fmp
fleet_allocation[,2] <- 1 - fleet_allocation[,1] # trawl gear allocation

# Condition closed-loop simulations, random recruitment
sim_list_rand <- condition_closed_loop_simulations(
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
  ISS_FishAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_FishLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
  ISS_SrvAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ISS_SrvLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
  ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
  ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets))
)

# Get constant B40 for comparison
global_spr <- Get_Reference_Points(data = om_values$data,
                                   rep = om_values$rep,
                                   SPR_x = 0.4,
                                   type = 'multi_region',
                                   what = 'global_SPR',
                                   calc_rec_st_yr = 20,
                                   rec_age = 2
)


# Run MSEs ----------------------------------------------------------------

# Single-region, current design
sim_env_current <- Setup_sim_env(sim_list = sim_list_rand)
sim_env_current <- add_aggregated_obj_to_simenv(sim_env = sim_env_current)
sim_env_current <- run_single_rg_closedloop_parallel(sim_env = sim_env_current, n_sims = n_sims,
                                                     fleet_allocation = fleet_allocation, lls_design_type = "current", n_cores = 7)
saveRDS(sim_env_current, here("outputs", "mse_results", "single_region_rand_current.RDS"))

# Single-region, historical design
sim_env_hist <- Setup_sim_env(sim_list = sim_list_rand)
sim_env_hist <- add_aggregated_obj_to_simenv(sim_env = sim_env_hist)
sim_env_hist <- run_single_rg_closedloop_parallel(sim_env = sim_env_hist, n_sims = n_sims, fleet_allocation = fleet_allocation,
                                                  lls_design_type = "historical", n_cores = 7)
saveRDS(sim_env_hist, here("outputs", "mse_results", "single_region_rand_hist.RDS"))

# Single-region, sample all regions
sim_env_all <- Setup_sim_env(sim_list = sim_list_rand)
sim_env_all <- add_aggregated_obj_to_simenv(sim_env = sim_env_all)
sim_env_all <- run_single_rg_closedloop_parallel(sim_env = sim_env_hist, n_sims = n_sims,
                                                 fleet_allocation = fleet_allocation,lls_design_type =  "all", n_cores = 7)
saveRDS(sim_env_all, here("outputs", "mse_results", "single_region_rand_all.RDS"))


