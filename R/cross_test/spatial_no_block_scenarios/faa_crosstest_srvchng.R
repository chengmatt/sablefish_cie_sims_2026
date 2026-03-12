# Purpose: To self test the FAA EM prior to the MSE year (No Selex Blocks)
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
source(here("R", "functions", "faa_em.R"))

# Read in OMs
srvchng <- readRDS(here("outputs", "mse_results", "spatial_noblock_scenarios", "five_region_base.RDS"))

# Simulation dimensions
y <- 95 # terminal year prior to MSE (no feedback)
n_sims <- 100 # number of sims

# Setup EMs ---------------------------------------------------------------
ems_grid <- expand.grid(
  # Fishery structure
  fishery_structure = c(
    "1F_1T", # Single fleet for trawl and fishery
    "3F_1T", # Three fixed-gear, 1 trawl gear
    "3F_2T", # Three fixed gear, 2 trawl gear (BS and AI + GOA Trawl)
    "3F_3T" # three fixed-gear, three trawl gear
  ),
  # Fishery selex options
  fishery_selex = c(
    "All_Logist", # All logistic
    "BS_Gamma", # BS Gamma
    "BS_AI_Gamma", # BS, AI Gamma
    "All_Gamma" # All Gamma
  ),
  # Survey selex options
  survey_selex = c(
    "All_Logist", # All logistic
    "BS_Gamma", # BS Gamma
    "BS_AI_Gamma", # BS, AI Gamma
    "All_Gamma" # All Gamma
  ),
  stringsAsFactors = FALSE
)

# Remove invalid combinations:
# 1F_1T can only use All_Logist or All_Gamma
ems_grid <- ems_grid[!(ems_grid$fishery_structure == "1F_1T" &
                         ems_grid$fishery_selex %in% c("BS_Gamma", "BS_AI_Gamma")), ]
ems_grid$model_id <- 1:nrow(ems_grid)

# Get fishery structure specs
get_fish_structure <- function(struct_name) {
  specs <- list(
    "1F_1T" = list(
      n_fish_fleets = 2,
      n_srv_fleets = 4,
      fish_blocks = c('none_Fleet_1', 'none_Fleet_2'),
      fish_prior_blocks = list("1" = 1, "2" = 1)
    ),
    "3F_1T" = list(
      n_fish_fleets = 4,
      n_srv_fleets = 4,
      fish_blocks = paste("none_Fleet_", 1:4, sep = ''),
      fish_prior_blocks = list("1" = 1, "2" = 1, "3" = 1, "4" = 1)
    ),
    "3F_2T" = list(
      n_fish_fleets = 5,
      n_srv_fleets = 4,
      fish_blocks = paste("none_Fleet_", 1:5, sep = ''),
      fish_prior_blocks = list("1" = 1, "2" = 1, "3" = 1, "4" = 1, "5" = 1)
    ),
    "3F_3T" = list(
      n_fish_fleets = 6,
      n_srv_fleets = 4,
      fish_blocks = paste("none_Fleet_", 1:6, sep = ''),
      fish_prior_blocks = list("1" = 1, "2" = 1, "3" = 1, "4" = 1, "5" = 1, "6" = 1)
    )
  )
  return(specs[[struct_name]])
}

# Get fishery selectivity models
get_fish_selex <- function(selex_name, n_fleets) {
  specs <- list(
    "All_Logist" = list(
      "2" = c("logist1_Fleet_1", "gamma_Fleet_2"),
      "4" = c(paste("logist1_Fleet_", 1:3, sep = ''), "gamma_Fleet_4"),
      "5" = c(paste("logist1_Fleet_", 1:3, sep = ''), paste("gamma_Fleet_", 4:5, sep = '')),
      "6" = c(paste("logist1_Fleet_", 1:3, sep = ''), paste("gamma_Fleet_", 4:6, sep = ''))
    ),
    "All_Gamma" = list(
      "2" = c("gamma_Fleet_1", "gamma_Fleet_2"),
      "4" = c(paste("gamma_Fleet_", 1:3, sep = ''), "gamma_Fleet_4"),
      "5" = c(paste("gamma_Fleet_", 1:3, sep = ''), paste("gamma_Fleet_", 4:5, sep = '')),
      "6" = c(paste("gamma_Fleet_", 1:3, sep = ''), paste("gamma_Fleet_", 4:6, sep = ''))
    ),
    "BS_Gamma" = list(
      "4" = c(paste("gamma_Fleet_", 1, sep = ''), paste("logist1_Fleet_", 2:3, sep = ''), "gamma_Fleet_4"),
      "5" = c(paste("gamma_Fleet_", 1, sep = ''), paste("logist1_Fleet_", 2:3, sep = ''), paste("gamma_Fleet_", 4:5, sep = '')),
      "6" = c(paste("gamma_Fleet_", 1, sep = ''), paste("logist1_Fleet_", 2:3, sep = ''), paste("gamma_Fleet_", 4:6, sep = ''))
    ),
    "BS_AI_Gamma" = list(
      "4" = c(paste("gamma_Fleet_", 1:2, sep = ''), paste("logist1_Fleet_", 3, sep = ''), "gamma_Fleet_4"),
      "5" = c(paste("gamma_Fleet_", 1:2, sep = ''), paste("logist1_Fleet_", 3, sep = ''), paste("gamma_Fleet_", 4:5, sep = '')),
      "6" = c(paste("gamma_Fleet_", 1:2, sep = ''), paste("logist1_Fleet_", 3, sep = ''), paste("gamma_Fleet_", 4:6, sep = ''))
    )
  )
  return(specs[[selex_name]][[as.character(n_fleets)]])
}

# Get survey selectivity models
get_survey_selex <- function(selex_name) {
  specs <- list(
    "All_Logist" = c("logist1_Fleet_1", "logist1_Fleet_2", "logist1_Fleet_3", "logist1_Fleet_4"),
    "BS_Gamma" = c("gamma_Fleet_1", "logist1_Fleet_2", "logist1_Fleet_3", "logist1_Fleet_4"),
    "BS_AI_Gamma" = c("gamma_Fleet_1", "gamma_Fleet_2", "logist1_Fleet_3", "logist1_Fleet_4"),
    "All_Gamma" = c("gamma_Fleet_1", "gamma_Fleet_2", "gamma_Fleet_3", "logist1_Fleet_4")
  )
  return(specs[[selex_name]])
}

ems <- list()

for (i in 1:nrow(ems_grid)) {

  # Get specs for this row
  fish_struct <- get_fish_structure(ems_grid$fishery_structure[i])
  fish_selex <- get_fish_selex(ems_grid$fishery_selex[i], fish_struct$n_fish_fleets)
  surv_selex <- get_survey_selex(ems_grid$survey_selex[i])

  # Build the model
  ems[[i]] <- list(
    faa_n_fish_fleets = fish_struct$n_fish_fleets,
    faa_n_srv_fleets = fish_struct$n_srv_fleets,

    # Fishery selectivity
    fish_sel_model = fish_selex,
    fish_fixed_sel_pars_spec = rep("est_all", fish_struct$n_fish_fleets),
    fish_selex_prior = build_selex_prior(
      fleets_blocks = fish_struct$fish_prior_blocks,
      sex_par = expand.grid(sex = 1:2, par = 1:2),
      region = 1,
      mu = 3.5,
      sd = 2
    ),

    # Survey selectivity
    srv_sel_model = surv_selex,
    srv_fixed_sel_pars_spec = rep("est_all", 4),
    srv_selex_prior = build_selex_prior(
      fleets_blocks = list(
        "1" = 1,
        "2" = 1,
        "3" = 1,
        "4" = 1
      ),
      sex_par = expand.grid(sex = 1:2, par = 1:2),
      region = 1,
      mu = 3.5,
      sd = 2
    )
  )
}


# Current design  -----------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 6)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_srvchng <- future_lapply(1:n_sims, function(i) {

    j <- 19 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

    # get faa data (need to iterate to build up)
    for(y1 in 65:y) {
      asmt_list <- faa_em(sim_env = srvchng,
                          y = y1,
                          sim = i,
                          srv_idx_se = 0.2,
                          age_lag = 1,
                          lls_design_type = 'current',
                          faa_n_fish_fleets = ems[[j]]$faa_n_fish_fleets,
                          faa_n_srv_fleets = ems[[j]]$faa_n_srv_fleets,
                          srv_wgt = 'numbers',
                          fish_wgt = 'numbers',
                          fish_sel_blocks = paste('none_Fleet_',1:ems[[j]]$faa_n_fish_fleets, sep = ''),
                          fish_sel_model = ems[[j]]$fish_sel_model,
                          fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                          fish_selex_prior = ems[[j]]$fish_selex_prior,
                          srv_sel_blocks =  paste('none_Fleet_',1:ems[[j]]$faa_n_srv_fleets, sep = ''),
                          srv_sel_model = ems[[j]]$srv_sel_model,
                          srv_fixed_sel_pars_spec = ems[[j]]$srv_fixed_sel_pars_spec,
                          srv_selex_prior = ems[[j]]$srv_selex_prior
      )
    }

    asmt_list$par$ln_fish_fixed_sel_pars[] <- log(5)
    asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

    model <- tryCatch(
      {
        # fit model
        model <- fit_model(
          asmt_list$data,
          asmt_list$par,
          asmt_list$map,
          NULL,
          2,
          silent = F
        )

        model$data <- asmt_list$data    # save data
        model$sd_rep <- sdreport(model)
        model
      },
      error = function(e) {
        NA
      }
    )

    p()
    model
  })
})

saveRDS(model_list_srvchng, here("outputs", "cross_test", "spatial_no_block_scenarios", "final_faa_crosstest_srvchng_current.RDS"))

# All design  -----------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_srvchng <- future_lapply(1:n_sims, function(i) {

    j <- 19 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

    # get faa data (need to iterate to build up)
    for(y1 in 65:y) {
      asmt_list <- faa_em(sim_env = srvchng,
                          y = y1,
                          sim = i,
                          srv_idx_se = 0.2,
                          age_lag = 1,
                          lls_design_type = 'all',
                          faa_n_fish_fleets = ems[[j]]$faa_n_fish_fleets,
                          faa_n_srv_fleets = ems[[j]]$faa_n_srv_fleets,
                          srv_wgt = 'numbers',
                          fish_wgt = 'numbers',
                          fish_sel_blocks = paste('none_Fleet_',1:ems[[j]]$faa_n_fish_fleets, sep = ''),
                          fish_sel_model = ems[[j]]$fish_sel_model,
                          fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                          fish_selex_prior = ems[[j]]$fish_selex_prior,
                          srv_sel_blocks =  paste('none_Fleet_',1:ems[[j]]$faa_n_srv_fleets, sep = ''),
                          srv_sel_model = ems[[j]]$srv_sel_model,
                          srv_fixed_sel_pars_spec = ems[[j]]$srv_fixed_sel_pars_spec,
                          srv_selex_prior = ems[[j]]$srv_selex_prior
      )
    }

    asmt_list$par$ln_fish_fixed_sel_pars[] <- log(5)
    asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

    model <- tryCatch(
      {
        # fit model
        model <- fit_model(
          asmt_list$data,
          asmt_list$par,
          asmt_list$map,
          NULL,
          2,
          silent = F
        )

        model$data <- asmt_list$data    # save data
        model$sd_rep <- sdreport(model)
        model
      },
      error = function(e) {
        NA
      }
    )

    p()
    model
  })
})

saveRDS(model_list_srvchng, here("outputs", "cross_test", "spatial_no_block_scenarios", "final_faa_crosstest_srvchng_all.RDS"))
