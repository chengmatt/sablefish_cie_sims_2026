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
  fishery_structure = c("1F_1T", "3F_1T", "3F_2T", "3F_3T"),
  fishery_selex     = c("All_Logist", "BS_Gamma", "BS_AI_Gamma", "All_Gamma"),
  survey_selex      = c("All_Logist", "BS_Gamma", "BS_AI_Gamma", "All_Gamma"),
  fish_blocks_type  = c("none", "blocks"),
  srv_blocks_type   = c("none", "blocks"),
  stringsAsFactors  = FALSE
)

# Remove invalid combinations:
# 1F_1T can only use All_Logist or All_Gamma for fishery selex
ems_grid <- ems_grid[!(ems_grid$fishery_structure == "1F_1T" &
                         ems_grid$fishery_selex %in% c("BS_Gamma", "BS_AI_Gamma")), ]
ems_grid$model_id <- seq_len(nrow(ems_grid))


# Helper: base fleet-structure specs
get_fish_structure <- function(struct_name) {
  specs <- list(
    "1F_1T" = list(n_fish_fleets = 2, n_srv_fleets = 4),
    "3F_1T" = list(n_fish_fleets = 4, n_srv_fleets = 4),
    "3F_2T" = list(n_fish_fleets = 5, n_srv_fleets = 4),
    "3F_3T" = list(n_fish_fleets = 6, n_srv_fleets = 4)
  )
  specs[[struct_name]]
}

# Helper: fishery selectivity models
get_fish_selex <- function(selex_name, n_fleets) {
  specs <- list(
    "All_Logist" = list(
      "2" = c("logist1_Fleet_1", "gamma_Fleet_2"),
      "4" = c(paste0("logist1_Fleet_", 1:3), "gamma_Fleet_4"),
      "5" = c(paste0("logist1_Fleet_", 1:3), paste0("gamma_Fleet_", 4:5)),
      "6" = c(paste0("logist1_Fleet_", 1:3), paste0("gamma_Fleet_", 4:6))
    ),
    "All_Gamma" = list(
      "2" = c("gamma_Fleet_1", "gamma_Fleet_2"),
      "4" = c(paste0("gamma_Fleet_", 1:3), "gamma_Fleet_4"),
      "5" = c(paste0("gamma_Fleet_", 1:3), paste0("gamma_Fleet_", 4:5)),
      "6" = c(paste0("gamma_Fleet_", 1:3), paste0("gamma_Fleet_", 4:6))
    ),
    "BS_Gamma" = list(
      "4" = c("gamma_Fleet_1",  paste0("logist1_Fleet_", 2:3), "gamma_Fleet_4"),
      "5" = c("gamma_Fleet_1",  paste0("logist1_Fleet_", 2:3), paste0("gamma_Fleet_", 4:5)),
      "6" = c("gamma_Fleet_1",  paste0("logist1_Fleet_", 2:3), paste0("gamma_Fleet_", 4:6))
    ),
    "BS_AI_Gamma" = list(
      "4" = c(paste0("gamma_Fleet_",   1:2), "logist1_Fleet_3", "gamma_Fleet_4"),
      "5" = c(paste0("gamma_Fleet_",   1:2), "logist1_Fleet_3", paste0("gamma_Fleet_", 4:5)),
      "6" = c(paste0("gamma_Fleet_",   1:2), "logist1_Fleet_3", paste0("gamma_Fleet_", 4:6))
    )
  )
  specs[[selex_name]][[as.character(n_fleets)]]
}

# Helper: survey selectivity models
get_survey_selex <- function(selex_name) {
  specs <- list(
    "All_Logist"  = paste0("logist1_Fleet_", 1:4),
    "BS_Gamma"    = c("gamma_Fleet_1",  paste0("logist1_Fleet_", 2:4)),
    "BS_AI_Gamma" = c(paste0("gamma_Fleet_", 1:2), paste0("logist1_Fleet_", 3:4)),
    "All_Gamma"   = c(paste0("gamma_Fleet_", 1:3), "logist1_Fleet_4")
  )
  specs[[selex_name]]
}

# Helper: fishery selectivity blocks
# Fixed-gear fleets (1:n_fixed) get two time blocks; trawl fleets get "none".
# Returns fish_sel_blocks vector and matching prior_blocks list.
get_fish_blocks <- function(struct_name, blocks_type) {
  n_total <- get_fish_structure(struct_name)$n_fish_fleets
  n_fixed <- switch(struct_name, "1F_1T" = 1, 3)
  n_trawl <- n_total - n_fixed

  if (blocks_type == "none") {
    blk_vec   <- paste0("none_Fleet_", seq_len(n_total))
    prior_blk <- setNames(as.list(rep(1, n_total)), as.character(seq_len(n_total)))
  } else {
    blk_vec <- c(
      as.vector(sapply(seq_len(n_fixed), function(f) c(
        paste0("Block_1_Year_1-55_Fleet_",       f),
        paste0("Block_2_Year_56-terminal_Fleet_", f)
      ))),
      paste0("none_Fleet_", (n_fixed + 1):n_total)
    )
    # fixed-gear fleets -> 2 blocks each, trawl fleets -> 1 block each
    prior_blk <- setNames(
      as.list(c(rep(2, n_fixed), rep(1, n_trawl))),
      as.character(seq_len(n_total))
    )
  }
  list(fish_sel_blocks = blk_vec, fish_prior_blocks = prior_blk)
}

# Helper: survey selectivity blocks ─────────────────────────────────────────
# Surveys fleets 1-3 get two time blocks; fleet 4 gets none.
get_srv_blocks <- function(blocks_type) {
  if (blocks_type == "none") {
    blk_vec   <- paste0("none_Fleet_", 1:4)
    prior_blk <- setNames(as.list(rep(1, 4)), as.character(1:4))
  } else {
    blk_vec <- c(
      "Block_1_Year_1-55_Fleet_1",       "Block_2_Year_56-terminal_Fleet_1",
      "Block_1_Year_1-55_Fleet_2",       "Block_2_Year_56-terminal_Fleet_2",
      "Block_1_Year_1-55_Fleet_3",       "Block_2_Year_56-terminal_Fleet_3",
      "none_Fleet_4"
    )
    prior_blk <- list("1" = 2, "2" = 2, "3" = 2, "4" = 1)
  }
  list(srv_sel_blocks = blk_vec, srv_prior_blocks = prior_blk)
}

# Build EM list
ems <- vector("list", nrow(ems_grid))

for (i in seq_len(nrow(ems_grid))) {

  fish_struct  <- get_fish_structure(ems_grid$fishery_structure[i])
  fish_selex   <- get_fish_selex(ems_grid$fishery_selex[i], fish_struct$n_fish_fleets)
  surv_selex   <- get_survey_selex(ems_grid$survey_selex[i])
  fish_blk     <- get_fish_blocks(ems_grid$fishery_structure[i], ems_grid$fish_blocks_type[i])
  srv_blk      <- get_srv_blocks(ems_grid$srv_blocks_type[i])

  ems[[i]] <- list(
    faa_n_fish_fleets = fish_struct$n_fish_fleets,
    faa_n_srv_fleets  = fish_struct$n_srv_fleets,

    # Fishery
    fish_sel_blocks          = fish_blk$fish_sel_blocks,
    fish_sel_model           = fish_selex,
    fish_fixed_sel_pars_spec = rep("est_all", fish_struct$n_fish_fleets),
    fish_selex_prior         = build_selex_prior(
      fleets_blocks = fish_blk$fish_prior_blocks,
      sex_par       = expand.grid(sex = 1:2, par = 1:2),
      region        = 1,
      mu            = 3.5,
      sd            = 2
    ),

    # Survey
    srv_sel_blocks           = srv_blk$srv_sel_blocks,
    srv_sel_model            = surv_selex,
    srv_fixed_sel_pars_spec  = rep("est_all", 4),
    srv_selex_prior          = build_selex_prior(
      fleets_blocks = srv_blk$srv_prior_blocks,
      sex_par       = expand.grid(sex = 1:2, par = 1:2),
      region        = 1,
      mu            = 3.5,
      sd            = 2
    )
  )
}


# Current design  -----------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 13)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_srvchng <- future_lapply(1:n_sims, function(i) {

    j <- 187 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

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
                          fish_sel_blocks = ems[[j]]$fish_sel_blocks,
                          fish_sel_model = ems[[j]]$fish_sel_model,
                          fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                          fish_selex_prior = ems[[j]]$fish_selex_prior,
                          srv_sel_blocks     = ems[[j]]$srv_sel_blocks,
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
plan(multisession, workers = 13)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_srvchng <- future_lapply(1:n_sims, function(i) {

    j <- 187 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

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
                          fish_sel_blocks = ems[[j]]$fish_sel_blocks,
                          fish_sel_model = ems[[j]]$fish_sel_model,
                          fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                          fish_selex_prior = ems[[j]]$fish_selex_prior,
                          srv_sel_blocks     = ems[[j]]$srv_sel_blocks,
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
