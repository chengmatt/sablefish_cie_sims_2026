# Purpose: To crash the GOA and see how single reigon and FAA models respond
# Creator: Matthew LH. Cheng
# 3/30/26


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
source(here("R", "functions", "single_region_em.R"))

# Load in OM
om <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "spt_rand_OM_lowsamp_dt.RDS"))

# Simulation dimensions
y <- 65 # terminal year prior to MSE (no feedback)
n_sims <- 100 # number of sims

# Run single-region EM  -----------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list <- future_lapply(1:n_sims, function(i) {

    # get single region data
    asmt_list <- single_region_em(
      sim_env = om,
      y = y,
      sim = i,
      srv_idx_se = 0.2,
      age_lag = 1,
      lls_design_type = "historical",
      srv_wgt = 'numbers',
      fish_wgt = 'numbers'
    )

    # fit model
    model <- fit_model(
      asmt_list$data,
      asmt_list$par,
      asmt_list$map,
      NULL,
      2,
      silent = TRUE
    )

    model$data <- asmt_list$data # save data
    model$sd_rep <- sdreport(model)

    p()
    model
  })
})

# Sim, EM list dimension
saveRDS(model_list, here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_sens_dt.RDS"))

# Quick look at bias
sgl_mod_list <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_sens_dt.RDS"))
ssb_sgl_rg_re <- array(NA, dim = c(y, n_sims))
rec_sgl_rg_re <- array(NA, dim = c(y, n_sims))

for(i in 1:100) {
  if(length(sgl_mod_list[[i]]) > 1) {
    ssb_sgl_rg_re[,i] <- (t(sgl_mod_list[[i]]$rep$SSB) - colSums(om$SSB[,1:65,i]) )/ colSums(om$SSB[,1:65,i])
    rec_sgl_rg_re[,i] <- (t(sgl_mod_list[[i]]$rep$Rec) - colSums(om$Rec[,1:65,i]) )/ colSums(om$Rec[,1:65,i])
  }
}

plot(apply(ssb_sgl_rg_re, 1, median, na.rm = T), type = 'l', ylim = c(-1, 1), main = 'Single Region SSB Error')
abline(h = 0, lty = 2, lwd = 2)

plot(apply(rec_sgl_rg_re, 1, median, na.rm = T), type = 'l', ylim = c(-1, 1))
abline(h = 0, lty = 2, lwd = 2)

for(i in 1:100) {
  if(i == 1) plot(sgl_mod_list[[i]]$rep$fish_sel[1,1,,1,1], ylim = c(0,1))
  else lines(sgl_mod_list[[i]]$rep$fish_sel[1,1,,1,1], ylim = c(0,1))
}

# Run FAA EM --------------------------------------------------------------
# setup EM grid and then filter to EM 19
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



# Run FAA EM --------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list <- future_lapply(1:n_sims, function(i) {

    j <- 19 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

    # get faa data
    asmt_list <- faa_em(sim_env = om,
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

saveRDS(model_list, here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_sens_dt.RDS"))

# Quick look at bias
faa_mod_list <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_sens_dt.RDS"))
ssb_faa_re <- array(NA, dim = c(y, n_sims))
rec_faa_re <- array(NA, dim = c(y, n_sims))

for(i in 1:100) {
  if(length(faa_mod_list[[i]]) > 1) {
    ssb_faa_re[,i] <- (t(faa_mod_list[[i]]$rep$SSB) - colSums(om$SSB[,1:65,i]) )/ colSums(om$SSB[,1:65,i])
    rec_faa_re[,i] <- (t(faa_mod_list[[i]]$rep$Rec) - colSums(om$Rec[,1:65,i]) )/ colSums(om$Rec[,1:65,i])
  }
}

plot(apply(ssb_faa_re, 1, median, na.rm = T), type = 'l', ylim = c(-1, 1), main = 'FAA SSB Error')
abline(h = 0, lty = 2, lwd = 2)

plot(apply(rec_faa_re, 1, median, na.rm = T), type = 'l', ylim = c(-1, 1))
abline(h = 0, lty = 2, lwd = 2)

# Run FAA Comp Weighting EM --------------------------------------------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list <- future_lapply(1:n_sims, function(i) {

    j <- 19 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

    # get faa data
    asmt_list <- faa_em(sim_env = om,
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

    # reset weights at 20 and then weight by catch
    asmt_list$data$ISS_FishAgeComps[,,,1:3] <- 100 * asmt_list$data$ObsCatch[,,1:3] / rowSums(asmt_list$data$ObsCatch[,,1:3])
    asmt_list$data$ISS_SrvAgeComps[,,,1:3] <- 100 * asmt_list$data$ObsSrvIdx[,,1:3] / rowSums(asmt_list$data$ObsSrvIdx[,,1:3])

    # fit model
    model <- tryCatch({
      # fit model
      model <- run_francis(
        asmt_list$data,
        asmt_list$par,
        asmt_list$map,
        NULL,
        n_francis_iter = 8,
        newton_loops = 3
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

saveRDS(model_list, here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_sens_dt_comp_wt.RDS"))

# Quick look at bias
faa_comp_wt_mod_list <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_sens_dt_comp_wt.RDS"))
ssb_faa_comp_wt_re <- array(NA, dim = c(y, n_sims))
rec_faa_comp_wt_re <- array(NA, dim = c(y, n_sims))

for(i in 1:100) {
  if(length(faa_comp_wt_mod_list[[i]]) > 1) {
    ssb_faa_comp_wt_re[,i] <- (t(faa_comp_wt_mod_list[[i]]$rep$SSB) - colSums(om$SSB[,1:65,i]) )/ colSums(om$SSB[,1:65,i])
    rec_faa_comp_wt_re[,i] <- (t(faa_comp_wt_mod_list[[i]]$rep$Rec) - colSums(om$Rec[,1:65,i]) )/ colSums(om$Rec[,1:65,i])
  }
}

plot(apply(ssb_faa_comp_wt_re, 1, median, na.rm = T), type = 'l', ylim = c(-1, 1))
abline(h = 0, lty = 2, lwd = 2)

plot(apply(rec_faa_comp_wt_re, 1, median, na.rm = T), type = 'l', ylim = c(-1, 1))
abline(h = 0, lty = 2, lwd = 2)
