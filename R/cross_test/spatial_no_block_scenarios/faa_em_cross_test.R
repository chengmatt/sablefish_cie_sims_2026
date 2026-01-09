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
lowsamp <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_lowsamp.RDS"))
lowsamp_movesense <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_lowsamp_tvmove_sens.RDS"))

# Simulation dimensions
y <- 65 # terminal year prior to MSE (no feedback)
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


# Run Models (First Half) --------------------------------------------------------------

handlers(global = TRUE)  # progress bar
plan(multisession, workers = 5)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    # container
    out_i <- vector("list", length(ems))

    for(j in 1:(length(ems) / 2)) {

      # get faa data
      asmt_list <- faa_em(sim_env = lowsamp,
                          y = y,
                          sim = i,
                          srv_idx_se = 0.2,
                          age_lag = 1,
                          lls_design_type = 'historical',
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

      asmt_list$par$ln_fish_fixed_sel_pars[] <- log(8)
      asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

      # fit model
      model <- tryCatch({
        # fit model
        model <- fit_model(
          asmt_list$data,
          asmt_list$par,
          asmt_list$map,
          NULL,
          3,
          silent = F
        )

        plot(model$rep$fish_sel[1,1,,1,4], ylim = c(0,1))

        model$data <- asmt_list$data # save data
        model$sd_rep <- sdreport(model)
        model
      }, error = function(e) {
        NA
      })

      out_i[[j]] <- model # output results

    }
    p()
    out_i
  })
})

# Sim, EM list dimension
saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_firsthalf.RDS"))

# Run Models (Second Half) --------------------------------------------------------------

handlers(global = TRUE)  # progress bar
plan(multisession, workers = 5)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    # container
    out_i <- vector("list", length(ems))

    for(j in ((length(ems) / 2) + 1):length(ems)) {

      # get faa data
      asmt_list <- faa_em(sim_env = lowsamp,
                          y = y,
                          sim = i,
                          srv_idx_se = 0.2,
                          age_lag = 1,
                          lls_design_type = 'historical',
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

      asmt_list$par$ln_fish_fixed_sel_pars[] <- log(10)
      asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

      # fit model
      model <- tryCatch({
        # fit model
        model <- fit_model(
          asmt_list$data,
          asmt_list$par,
          asmt_list$map,
          NULL,
          3,
          silent = F
        )
        model$data <- asmt_list$data # save data
        model$sd_rep <- sdreport(model)
        model
      }, error = function(e) {
        NA
      })
      out_i[[j]] <- model # output results
    }
    p()
    out_i
  })
})

# Sim, EM list dimension
saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_secondhalf.RDS"))

# Run Models (Finalized Param FAA W Francis)  -----------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    j <- 19 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

    # get faa data
    asmt_list <- faa_em(sim_env = lowsamp,
                        y = y,
                        sim = i,
                        srv_idx_se = 0.2,
                        age_lag = 1,
                        lls_design_type = 'historical',
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

    asmt_list$par$ln_fish_fixed_sel_pars[] <- log(10)
    asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

    model <- tryCatch(
      {
        # fit model
        model_francis <- run_francis(
          asmt_list$data,
          asmt_list$par,
          asmt_list$map,
          NULL,
          n_francis_iter = 10,
          newton_loops = 2
        )

        model <- model_francis$obj      # return model object
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

saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_finalFAA_francis.RDS"))

# Run Models (Finalized Param FAA w/o Francis - Sensitivity)  -----------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    j <- 19 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

    # get faa data
    asmt_list <- faa_em(sim_env = lowsamp_movesense,
                        y = y,
                        sim = i,
                        srv_idx_se = 0.2,
                        age_lag = 1,
                        lls_design_type = 'historical',
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

    asmt_list$par$ln_fish_fixed_sel_pars[] <- log(10)
    asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

    # fit model
    model <- tryCatch({
      # fit model
      model <- fit_model(
        asmt_list$data,
        asmt_list$par,
        asmt_list$map,
        NULL,
        3,
        silent = F
      )
      model$data <- asmt_list$data # save data
      model$sd_rep <- sdreport(model)
      model
    }, error = function(e) {
      NA
    })

    p()
    model
  })
})

saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_movesense_finalFAA.RDS"))

# Run Models (Finalized Param FAA w Francis - Sensitivity)  -----------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    j <- 19 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

    # get faa data
    asmt_list <- faa_em(sim_env = lowsamp_movesense,
                        y = y,
                        sim = i,
                        srv_idx_se = 0.2,
                        age_lag = 1,
                        lls_design_type = 'historical',
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

    asmt_list$par$ln_fish_fixed_sel_pars[] <- log(10)
    asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

    model <- tryCatch(
      {
        # fit model
        model_francis <- run_francis(
          asmt_list$data,
          asmt_list$par,
          asmt_list$map,
          NULL,
          n_francis_iter = 10,
          newton_loops = 2
        )

        model <- model_francis$obj      # return model object
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

saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_movesense_finalFAA_francis.RDS"))

# Run FAA Model (w/o Francis Retro) ---------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  retros_wo_francis <- future_lapply(1:n_sims, function(i) {

    j <- 19 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

    # get faa data
    asmt_list <- faa_em(sim_env = lowsamp,
                        y = y,
                        sim = i,
                        srv_idx_se = 0.2,
                        age_lag = 1,
                        lls_design_type = 'historical',
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

    asmt_list$par$ln_fish_fixed_sel_pars[] <- log(10)
    asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

    retros <- tryCatch(
      {
        # fit model via retros
        retros <- do_retrospective(
          8,
          data = asmt_list$data,
          parameters = asmt_list$par,
          mapping = asmt_list$map,
          random = NULL,
          do_par = FALSE,
          newton_loops = 3,
          do_francis = FALSE
        )

        retros$Sim <- i  # name retros
        retros
      },
      error = function(e) {
        NA
      }
    )

    p()
    retros

  })
})

retros_wo_francis <- retros_wo_francis[!is.na(retros_wo_francis)]
retros_wo_francis_df <- data.table::rbindlist(retros_wo_francis)
write.csv(retros_wo_francis_df, here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_finalized_cross_test_retro.csv"))

# Run FAA Model (w/ Francis Retro) ---------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 7)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  retros_w_francis <- future_lapply(1:n_sims, function(i) {

    j <- 19 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

    # get faa data
    asmt_list <- faa_em(sim_env = lowsamp,
                        y = y,
                        sim = i,
                        srv_idx_se = 0.2,
                        age_lag = 1,
                        lls_design_type = 'historical',
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

    asmt_list$par$ln_fish_fixed_sel_pars[] <- log(10)
    asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

    retros <- tryCatch(
      {
        # fit model via retros
        retros <- do_retrospective(
          n_retro = 8,
          data = asmt_list$data,
          parameters = asmt_list$par,
          mapping = asmt_list$map,
          random = NULL,
          do_par = FALSE,
          newton_loops = 3,
          do_francis = TRUE,
          n_francis_iter = 10
        )

        retros$Sim <- i  # name retros
        retros
      },
      error = function(e) {
        NA
      }
    )

    p()
    retros

  })
})

retros_w_francis <- retros_w_francis[!is.na(retros_w_francis)]
retros_w_francis_df <- data.table::rbindlist(retros_w_francis)
write.csv(retros_w_francis_df, here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_finalized_cross_test_retro_francis.csv"))

# Process Results ---------------------------------------------------------
model_list_all <- vector('list', 100)
model_list_first <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_firsthalf.RDS"))
model_list_second <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_secondhalf.RDS"))
finalized_model_francis <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_finalFAA_francis.RDS"))

# Combine models
for(i in 1:100) {
  for(j in 1:nrow(ems_grid)) {
    if(j <= (nrow(ems_grid) / 2)) model_list_all[[i]][[j]] <- model_list_first[[i]][[j]]
    else model_list_all[[i]][[j]] <- model_list_second[[i]][[j]]
  }
}

# Input finalized francis model in as well
for(i in 1:100) model_list_all[[i]][[nrow(ems_grid) + 1]] <- finalized_model_francis[[i]]

# storage containers
ssb_store <- array(NA, dim = c(65, 100, nrow(ems_grid) + 1))
rec_store <- array(NA, dim = c(65, 100, nrow(ems_grid) + 1))
dep_store <- array(NA, dim = c(65, 100, nrow(ems_grid) + 1))

for(mod in 1:(nrow(ems_grid) + 1)) {
  for(i in 1:100) {
    tryCatch({
      if(length(model_list_all[[i]][[mod]]) > 1) {
        ssb_store[,i,mod] <- (t(model_list_all[[i]][[mod]]$rep$SSB) - colSums(lowsamp$SSB[,-66,1])) / colSums(lowsamp$SSB[,-66,1])
        rec_store[,i,mod] <- (t(model_list_all[[i]][[mod]]$rep$Rec) - colSums(lowsamp$Rec[,-66,1])) / colSums(lowsamp$Rec[,-66,1])
        depletion_om <- colSums(lowsamp$SSB[,-66,1]) / sum(lowsamp$SSB[,1,1])
        depletion_em <- t(model_list_all[[i]][[mod]]$rep$SSB) / t(model_list_all[[i]][[mod]]$rep$SSB)[1]
        dep_store[,i,mod] <- (t(depletion_em) - depletion_om) / depletion_om
      }
    }, error = function(e) NULL)
  }
}

# Add one more row
ems_grid <- rbind(ems_grid, data.frame(fishery_structure = '3F_2T', fishery_selex = 'BS_Gamma', survey_selex = 'BS_Gamma', model_id = nrow(ems_grid) + 1))

# Summarized results (SSB)
ssb_sum_results <- reshape2::melt(ssb_store) %>%
  group_by(Var1, Var3) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  left_join(ems_grid, by = c("Var3" = "model_id"))

# Summarized results (Rec)
rec_sum_results <- reshape2::melt(rec_store) %>%
  group_by(Var1, Var3) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  left_join(ems_grid, by = c("Var3" = "model_id"))

# Summarized results (Depletion)
dep_sum_results <- reshape2::melt(dep_store) %>%
  group_by(Var1, Var3) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  left_join(ems_grid, by = c("Var3" = "model_id"))

ems_grid %>%
  filter(!str_detect(fishery_selex, 'All_Gamma|All_Logist'),
         !str_detect(survey_selex, "All_Gamma|All_Logist"),
         fishery_structure == '3F_2T') %>% view()

# Only look at 3 fixed gear fleets, 2 trawl fleet variants for now
ssb_sum_results %>%
  filter(fishery_structure == '3F_2T') %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

# Removing all fishery selex = gamma
ssb_sum_results %>%
  filter(fishery_structure == '3F_2T',
         !str_detect(fishery_selex, 'All_Gamma')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

# Removing all survey selex = gamma
ssb_sum_results %>%
  filter(fishery_structure == '3F_2T',
         !str_detect(fishery_selex, 'All_Gamma'),
         !str_detect(survey_selex, "All_Gamma")) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

# Removing all fishery selex = logist
ssb_sum_results %>%
  filter(fishery_structure == '3F_2T',
         !str_detect(fishery_selex, 'All_Gamma|All_Logist'),
         !str_detect(survey_selex, "All_Gamma")) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

# Removing all survey selex = logist
ssb_sum_results %>%
  filter(fishery_structure == '3F_2T',
         !str_detect(fishery_selex, 'All_Gamma|All_Logist'),
         !str_detect(survey_selex, "All_Gamma|All_Logist")) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

# Only look at 3 fixed gear fleets, 2 trawl fleet variants for now
rec_sum_results %>%
  filter(fishery_structure == '3F_2T') %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

# Removing all fishery selex = gamma
rec_sum_results %>%
  filter(fishery_structure == '3F_2T',
         !str_detect(fishery_selex, 'All_Gamma')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

# Removing all survey selex = gamma
rec_sum_results %>%
  filter(fishery_structure == '3F_2T',
         !str_detect(fishery_selex, 'All_Gamma'),
         !str_detect(survey_selex, "All_Gamma")) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

# Removing all fishery selex = logist
rec_sum_results %>%
  filter(fishery_structure == '3F_2T',
         !str_detect(fishery_selex, 'All_Gamma|All_Logist'),
         !str_detect(survey_selex, "All_Gamma")) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

# Removing all survey selex = logist
rec_sum_results %>%
  filter(fishery_structure == '3F_2T',
         !str_detect(fishery_selex, 'All_Gamma|All_Logist'),
         !str_detect(survey_selex, "All_Gamma|All_Logist")) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

# Optimal seems to be BS Gamma for both fishery and survey ...
# Look at how that changes with fleet structure
ssb_sum_results %>%
  filter(fishery_selex == 'BS_Gamma',
         survey_selex == 'BS_Gamma'
  ) %>%
  mutate(fishery_structure = ifelse(Var3 == 57, '3F_2T_Francis', fishery_structure)) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr,
             group = fishery_structure, color = fishery_structure, fill = fishery_structure)) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

rec_sum_results %>%
  filter(fishery_selex == 'BS_Gamma',
         survey_selex == 'BS_Gamma'
         ) %>%
  mutate(fishery_structure = ifelse(Var3 == 57, '3F_2T_Francis', fishery_structure)) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr,
             group = fishery_structure, color = fishery_structure, fill = fishery_structure)) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

dep_sum_results %>%
  filter(fishery_selex == 'BS_Gamma',
         survey_selex == 'BS_Gamma'
  ) %>%
  mutate(fishery_structure = ifelse(Var3 == 57, '3F_2T_Francis', fishery_structure)) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr,
             group = fishery_structure, color = fishery_structure, fill = fishery_structure, group = Var3)) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Depletion RE')


# Sensitivity Run Results ---------------------------------------------------------
model_list_lowsamp <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_movesense_finalFAA.RDS"))
model_list_lowsamp_francis <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_movesense_finalFAA_francis.RDS"))

ssb_em_results <- array(NA, dim = c(2, length(1:y), n_sims))
rec_em_results <- array(NA, dim = c(2, length(1:y), n_sims))
ssb_om_results <- array(NA, dim = c(length(1:y), n_sims))
rec_om_results <- array(NA, dim = c(length(1:y), n_sims))

for(i in 1:n_sims) {
  # get aggregated OM results
  ssb_om_results[,i] <- colSums(lowsamp_movesense$SSB[,1:y,i])
  rec_om_results[,i] <- colSums(lowsamp_movesense$Rec[,1:y,i])

  # get EM results
  ssb_em_results[1,,i] <- t(model_list_lowsamp[[i]]$rep$SSB)
  rec_em_results[1,,i] <- t(model_list_lowsamp[[i]]$rep$Rec)
  if(length(model_list_lowsamp_francis[[i]]) > 1) ssb_em_results[2,,i] <- t(model_list_lowsamp_francis[[i]]$rep$SSB)
  if(length(model_list_lowsamp_francis[[i]]) > 1) rec_em_results[2,,i] <- t(model_list_lowsamp_francis[[i]]$rep$Rec)
} # end i


par(mfrow = c(2,2))
# SSB
plot(apply((ssb_em_results[1,,] - ssb_om_results) / ssb_om_results, 1, median, na.rm = T), type = 'l', ylim = c(-0.25, 0.25),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in SSB')
lines(apply((ssb_em_results[2,,] - ssb_om_results) / ssb_om_results, 1, median, na.rm = T), col = 'red', lwd = 3, lty = 2)
abline(h = 0, lty = 2, lwd = 3)
# SSB abs
plot(apply(abs((ssb_em_results[1,,] - ssb_om_results) / ssb_om_results), 1, median, na.rm = T), type = 'l', ylim = c(0, 0.25),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in SSB')
lines(apply(abs((ssb_em_results[2,,] - ssb_om_results) / ssb_om_results), 1, median, na.rm = T), col = 'red', lwd = 3, lty = 2)

# Rec
plot(apply((rec_em_results[1,,] - rec_om_results) / rec_om_results, 1, median, na.rm = T), type = 'l', ylim = c(-1, 1),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in Rec')
lines(apply((rec_em_results[2,,] - rec_om_results) / rec_om_results, 1, median, na.rm = T), col = 'red', lwd = 3, lty = 2)
abline(h = 0, lty = 2, lwd = 3)
# Rec abs
plot(apply(abs((rec_em_results[1,,] - rec_om_results) / rec_om_results), 1, median, na.rm = T), type = 'l', ylim = c(0, 1),
     col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in Rec')
lines(apply(abs((rec_em_results[2,,] - rec_om_results) / rec_om_results), 1, median, na.rm = T), col = 'red', lwd = 3, lty = 2)

# Retrospective Results ---------------------------------------------------
retro_nofrancis <- read_csv(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_finalized_cross_test_retro.csv"))
retro_wfrancis <- read_csv(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_finalized_cross_test_retro_francis.csv"))

# Summarize info
retro_nofrancis_sum <- retro_nofrancis %>%
  group_by(Region, Year, peel, Type) %>%
  summarize(lwr_95 = quantile(value, 0.025),
            upr_95 = quantile(value, 0.975),
            value = median(value)) %>%
  mutate(type = 'No Francis') %>%
  ungroup()


# Summarize info
retro_wfrancis_sum <- retro_wfrancis %>%
  group_by(Region, Year, peel, Type) %>%
  summarize(lwr_95 = quantile(value, 0.025),
            upr_95 = quantile(value, 0.975),
            value = median(value)) %>%
  mutate(type = 'W Francis') %>%
  ungroup()

# get relative difference
retro_nofrancis_rd <- get_retrospective_relative_difference(retro_nofrancis_sum) %>% mutate(type = 'No Francis')
retro_francis_rd <- get_retrospective_relative_difference(retro_wfrancis_sum)%>% mutate(type = 'W Francis')

# get mohn's rho
mohns_rho_nofrancis <- retro_nofrancis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
  dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'No Francis')
mohns_rho_wfrancis <- retro_francis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
  dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'W Francis')

# bind
retro_sgl_region <- rbind(retro_wfrancis_sum, retro_nofrancis_sum)
mohns_rho <- rbind(mohns_rho_nofrancis, mohns_rho_wfrancis)
ret_df <- rbind(retro_nofrancis_rd, retro_francis_rd)

# Generally look the same
ggplot() +
  geom_hline(yintercept = 0,
             lty = 2,
             lwd = 1.3) +
  geom_line(
    ret_df,
    mapping = aes(
      x = Year,
      y = rd,
      group = as.numeric(peel),
      color = as.numeric(peel)
    ),
    lwd = 1.3
  ) +
  geom_point(
    ret_df %>% dplyr::filter(peel ==
                               max(Year) - Year),
    mapping = aes(
      x = Year,
      y = rd,
      group = as.numeric(peel),
      fill = as.numeric(peel)
    ),
    pch = 21,
    size = 6
  ) +
  geom_text(
    mohns_rho,
    mapping = aes(
      x = -Inf,
      y = Inf,
      label = paste("Mohns Rho:", round(rho, 4))
    ),
    hjust = -0.3,
    vjust = 3,
    size = 5
  )+
  guides(color = guide_colourbar(barwidth = 15, barheight = 1.3)) +
  labs(x = "Year",
       y = "Relative Difference from Terminal Year",
       color = "Retrospective Year",
       fill = "Retrospective Year") +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  facet_grid(type ~ Type, scales = "free") +
  theme_sablefish() +
  theme(legend.position = "top")

ggplot() +
  geom_line(
    retro_sgl_region %>%
      filter(peel != 0),
    mapping = aes(
      x = Year,
      y = value,
      group = peel,
      color = peel
    ),
    lwd = 1
  ) +
  geom_line(
    retro_sgl_region %>%
      filter(peel == 0),
    mapping = aes(x = Year, y = value),
    lty = 2,
    lwd = 1
  ) +
  geom_ribbon(
    retro_sgl_region %>%
      filter(peel != 0),
    mapping = aes(
      x = Year,
      y = value,
      group = peel,
      ymin = lwr_95,
      ymax = upr_95,
      fill = peel
    ),
    color = NA, alpha = 0.2
  ) +
  geom_text(mohns_rho,
            mapping = aes(x = -Inf, y = Inf, label = paste("Mohns Rho:", round(rho, 4))),
            hjust = -0.3, vjust = 3, size = 5) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +

  coord_cartesian(ylim = c(0, NA)) +
  facet_grid(Type ~ type, scales = "free_y") +
  theme_sablefish() +
  labs(x = "Year", y = "Value", color = "Peel")



# Compare Retro to True Values --------------------------------------------
# Input true values into peel == 0 for evaluation
true_ssb <- colSums(lowsamp$SSB[,-66,1])
true_rec <- colSums(lowsamp$Rec[,-66,1])
true_values <- c(true_ssb, true_rec)

retro_nofrancis[retro_nofrancis$peel == 0, ]$value <- rep(true_values, length(unique(retro_nofrancis$Sim)))
retro_wfrancis[retro_wfrancis$peel == 0, ]$value <- rep(true_values, length(unique(retro_wfrancis$Sim)))

# Summarize info
retro_nofrancis_sum <- retro_nofrancis %>%
  group_by(Region, Year, peel, Type) %>%
  summarize(lwr_95 = quantile(value, 0.025),
            upr_95 = quantile(value, 0.975),
            value = median(value)) %>%
  mutate(type = 'No Francis') %>%
  ungroup()


# Summarize info
retro_wfrancis_sum <- retro_wfrancis %>%
  group_by(Region, Year, peel, Type) %>%
  summarize(lwr_95 = quantile(value, 0.025),
            upr_95 = quantile(value, 0.975),
            value = median(value)) %>%
  mutate(type = 'W Francis') %>%
  ungroup()

# get relative difference
retro_nofrancis_rd <- get_retrospective_relative_difference(retro_nofrancis_sum) %>% mutate(type = 'No Francis')
retro_francis_rd <- get_retrospective_relative_difference(retro_wfrancis_sum)%>% mutate(type = 'W Francis')

# get mohn's rho
mohns_rho_nofrancis <- retro_nofrancis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
  dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'No Francis')
mohns_rho_wfrancis <- retro_francis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
  dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'W Francis')

# bind
retro_sgl_region <- rbind(retro_wfrancis_sum, retro_nofrancis_sum)
mohns_rho <- rbind(mohns_rho_nofrancis, mohns_rho_wfrancis)
ret_df <- rbind(retro_nofrancis_rd, retro_francis_rd)
# Generally look the same
ggplot() +
  geom_hline(yintercept = 0,
             lty = 2,
             lwd = 1.3) +
  geom_line(
    ret_df,
    mapping = aes(
      x = Year,
      y = rd,
      group = as.numeric(peel),
      color = as.numeric(peel)
    ),
    lwd = 1.3
  ) +
  geom_point(
    ret_df %>% dplyr::filter(peel ==
                               max(Year) - Year),
    mapping = aes(
      x = Year,
      y = rd,
      group = as.numeric(peel),
      fill = as.numeric(peel)
    ),
    pch = 21,
    size = 6
  ) +
  geom_text(
    mohns_rho,
    mapping = aes(
      x = -Inf,
      y = Inf,
      label = paste("Mohns Rho:", round(rho, 4))
    ),
    hjust = -0.3,
    vjust = 3,
    size = 5
  )+
  guides(color = guide_colourbar(barwidth = 15, barheight = 1.3)) +
  labs(x = "Year",
       y = "Relative Difference from Terminal Year",
       color = "Retrospective Year",
       fill = "Retrospective Year") +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  facet_grid(type ~ Type, scales = "free") +
  theme_sablefish() +
  theme(legend.position = "top")

ggplot() +
  geom_line(
    retro_sgl_region %>%
      filter(peel != 0),
    mapping = aes(
      x = Year,
      y = value,
      group = peel,
      color = peel
    ),
    lwd = 1
  ) +
  geom_line(
    retro_sgl_region %>%
      filter(peel == 0),
    mapping = aes(x = Year, y = value),
    lty = 2,
    lwd = 1
  ) +
  geom_ribbon(
    retro_sgl_region %>%
      filter(peel != 0),
    mapping = aes(
      x = Year,
      y = value,
      group = peel,
      ymin = lwr_95,
      ymax = upr_95,
      fill = peel
    ),
    color = NA, alpha = 0.2
  ) +
  geom_text(mohns_rho,
            mapping = aes(x = -Inf, y = Inf, label = paste("Mohns Rho:", round(rho, 4))),
            hjust = -0.3, vjust = 3, size = 5) +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +
  coord_cartesian(ylim = c(0, NA)) +
  facet_grid(Type ~ type, scales = "free_y") +
  theme_sablefish() +
  labs(x = "Year", y = "Value", color = "Peel")



