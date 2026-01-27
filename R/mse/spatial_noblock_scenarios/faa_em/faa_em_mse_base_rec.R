# Purpose: To run a five region OM, with a FAA EM (base recruitment)
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 11/26/25

# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)
library(furrr)
library(progressr)

# Read in model output
om_values <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_NoBlock_model_results.RDS"))

source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "faa_em.R"))

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
n_sims <- 5              # Number of replicate simulations
n_regions <- 5             # number of regions
n_fish_fleets <- 2         # number of fishery fleets
n_srv_fleets <- 3          # number of survey fleets

# get fleet allocation by region
fleet_allocation <- array(NA, dim = c(n_regions, n_fish_fleets))
fleet_allocation[,1] <- c(0.5,0.75,0.8,0.8,0.95) # from fmp
fleet_allocation[,2] <- 1 - fleet_allocation[,1] # trawl gear allocation

# Condition closed-loop simulations, random recruitment
sim_list <- condition_closed_loop_simulations(
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
sim_env_current <- Setup_sim_env(sim_list = sim_list)
sim_env_current <- add_aggregated_obj_to_simenv(sim_env = sim_env_current, type = 'faa', faa_n_fish_fleets = 4, faa_n_srv_fleets = 4)
sim_env_current <- run_faa_closedloop_parallel(sim_env = sim_env_current,
                                               n_sims = n_sims,
                                               fleet_allocation = fleet_allocation,
                                               srv_idx_se =  0.2,
                                               age_lag = 1,
                                               lls_design_type = 'current',
                                               srv_wgt = 'numbers',
                                               fish_wgt = 'numbers',
                                               faa_n_fish_fleets = 4,
                                               faa_n_srv_fleets = 4,
                                               fish_sel_model = c("gamma_Fleet_1",   "logist1_Fleet_2",
                                                                  "logist1_Fleet_3", "gamma_Fleet_4"),
                                               srv_sel_model = c('gamma_Fleet_1', "logist1_Fleet_2",
                                                                 "logist1_Fleet_3", "logist1_Fleet_4"),
                                               fish_selex_prior = expand.grid(region = 1, fleet = 1:4, block = 1, sex = 1:2, par = 1:2, mu =  3.5, sd = 2),
                                               srv_selex_prior = expand.grid(region = 1, fleet = 1:4, block = 1, sex = 1:2, par = 1:2, mu = 3.5, sd = 2),
                                               n_cores = 7)

saveRDS(sim_env_current, here("outputs", "mse_results", "spatial_noblock_scenarios", "faa_base.RDS"))
