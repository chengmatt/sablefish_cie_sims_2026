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


# Run FAA EMs -----------------------------------------------------------------
split_names <- c("firstpart", "secondpart", "thirdpart", "fourthpart", "fifthpart", "sixthpart")
n <- nrow(ems_grid)
n_parts <- length(split_names)
splits <- lapply(1:n_parts, function(i) seq(ceiling((i-1)*n/n_parts) + 1, ceiling(i*n/n_parts)))

for (k in seq_along(split_names)) {

  handlers(global = TRUE)
  plan(multisession, workers = 13)
  options(future.globals.maxSize = 15e9)

  with_progress({
    p <- progressor(steps = n_sims)

    model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

      out_i <- vector("list", length(ems))

      for (j in splits[[k]]) {

        asmt_list <- faa_em(
          sim_env = lowsamp, y = y, sim = i,
          srv_idx_se = 0.2, age_lag = 1, lls_design_type = 'historical',
          faa_n_fish_fleets       = ems[[j]]$faa_n_fish_fleets,
          faa_n_srv_fleets        = ems[[j]]$faa_n_srv_fleets,
          srv_wgt = 'numbers', fish_wgt = 'numbers',
          fish_sel_blocks         = ems[[j]]$fish_sel_blocks,
          fish_sel_model          = ems[[j]]$fish_sel_model,
          fish_fixed_sel_pars_spec= ems[[j]]$fish_fixed_sel_pars_spec,
          fish_selex_prior        = ems[[j]]$fish_selex_prior,
          srv_sel_blocks          = ems[[j]]$srv_sel_blocks,
          srv_sel_model           = ems[[j]]$srv_sel_model,
          srv_fixed_sel_pars_spec = ems[[j]]$srv_fixed_sel_pars_spec,
          srv_selex_prior         = ems[[j]]$srv_selex_prior
        )

        if (ems[[j]]$faa_n_fish_fleets != 2) {
          asmt_list$par$ln_fish_fixed_sel_pars[] <- log(8)
          asmt_list$par$ln_srv_fixed_sel_pars[]  <- log(2)
        } else {
          asmt_list$par$ln_fish_fixed_sel_pars[] <- log(3)
          asmt_list$par$ln_srv_fixed_sel_pars[]  <- log(1)
        }

        out_i[[j]] <- tryCatch({
          model <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 3, silent = FALSE)
          model$data   <- asmt_list$data
          model$sd_rep <- sdreport(model)
          model
        }, error = function(e) NA)
      }

      p()
      out_i
    })
  })

  saveRDS(
    model_list_lowsamp,
    here("outputs", "cross_test", "spatial_no_block_scenarios",
         paste0("faa_crosstest_lowsamp_", split_names[k], ".RDS"))
  )
}


# Run Models (Finalized Param FAA W Francis)  -----------------------------
handlers(global = TRUE)  # progress bar
plan(multisession, workers = 13)
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
                        fish_sel_blocks = ems[[j]]$fish_sel_blocks,
                        fish_sel_model = ems[[j]]$fish_sel_model,
                        fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                        fish_selex_prior = ems[[j]]$fish_selex_prior,
                        srv_sel_blocks     = ems[[j]]$srv_sel_blocks,
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
plan(multisession, workers = 13)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    j <- 187 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

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
                        fish_sel_blocks = ems[[j]]$fish_sel_blocks,
                        fish_sel_model = ems[[j]]$fish_sel_model,
                        fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                        fish_selex_prior = ems[[j]]$fish_selex_prior,
                        srv_sel_blocks     = ems[[j]]$srv_sel_blocks,
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
plan(multisession, workers = 13)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  model_list_lowsamp <- future_lapply(1:n_sims, function(i) {

    j <- 187 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

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
                        fish_sel_blocks = ems[[j]]$fish_sel_blocks,
                        fish_sel_model = ems[[j]]$fish_sel_model,
                        fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                        fish_selex_prior = ems[[j]]$fish_selex_prior,
                        srv_sel_blocks     = ems[[j]]$srv_sel_blocks,
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
plan(multisession, workers = 13)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  retros_wo_francis <- future_lapply(1:n_sims, function(i) {

    j <- 187 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

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
                        fish_sel_blocks = ems[[j]]$fish_sel_blocks,
                        fish_sel_model = ems[[j]]$fish_sel_model,
                        fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                        fish_selex_prior = ems[[j]]$fish_selex_prior,
                        srv_sel_blocks     = ems[[j]]$srv_sel_blocks,
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
plan(multisession, workers = 13)
options(future.globals.maxSize = 15e9)

# loop through
with_progress({
  p <- progressor(steps = n_sims)

  retros_w_francis <- future_lapply(1:n_sims, function(i) {

    j <- 187 # 3F_1T, BS_Gamm (Srv), BS_Gamma (Srv)

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
                        fish_sel_blocks = ems[[j]]$fish_sel_blocks,
                        fish_sel_model = ems[[j]]$fish_sel_model,
                        fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                        fish_selex_prior = ems[[j]]$fish_selex_prior,
                        srv_sel_blocks     = ems[[j]]$srv_sel_blocks,
                        srv_sel_model = ems[[j]]$srv_sel_model,
                        srv_fixed_sel_pars_spec = ems[[j]]$srv_fixed_sel_pars_spec,
                        srv_selex_prior = ems[[j]]$srv_selex_prior
    )

    asmt_list$par$ln_fish_fixed_sel_pars[] <- log(3)
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
n <- nrow(ems_grid)
split_names <- c("firstpart", "secondpart", "thirdpart", "fourthpart", "fifthpart", "sixthpart")
n_parts <- length(split_names)
splits <- lapply(1:n_parts, function(i) seq(ceiling((i-1)*n/n_parts) + 1, ceiling(i*n/n_parts)))

ssb_store <- array(NA, dim = c(65, 100, n + 1))
rec_store <- array(NA, dim = c(65, 100, n + 1))
dep_store <- array(NA, dim = c(65, 100, n + 1))
conv_df <- data.frame()

depletion_om <- colSums(lowsamp$SSB[,-66,1]) / sum(lowsamp$SSB[,1,1])
ssb_om <- colSums(lowsamp$SSB[,-66,1])
rec_om <- colSums(lowsamp$Rec[,-66,1])


for(part in 1:n_parts) {
  part_data <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios",
                            paste0("faa_crosstest_lowsamp_", split_names[part], ".RDS")))
  idx <- splits[[part]]

  for(mod in idx) {
    for(i in 1:100) {
      tryCatch({
        if(length(part_data[[i]][[mod]]) > 1) {

          ssb_em <- as.vector(t(part_data[[i]][[mod]]$rep$SSB))
          rec_em <- as.vector(t(part_data[[i]][[mod]]$rep$Rec))

          ssb_store[, i, mod] <- (ssb_em - ssb_om) / ssb_om
          rec_store[, i, mod] <- (rec_em - rec_om) / rec_om

          depletion_em <- ssb_em / ssb_em[1]
          dep_store[, i, mod] <- (depletion_em - depletion_om) / depletion_om

          # get convergence
          tmp_conv <- data.frame(model = mod, sim = i, pd = part_data[[i]][[mod]]$sd_rep$pdHess,
                                 grad = max(abs(part_data[[i]][[mod]]$sd_rep$gradient.fixed)))
          conv_df <- rbind(tmp_conv, conv_df)
        }
      }, error = function(e) NULL)
    }
  }

  rm(part_data)
  gc()
}

conv_df %>%
  mutate(conv = ifelse(pd == TRUE & grad <= 1e-5, T, F)) %>%
  group_by(model) %>%
  summarize(sum = sum(conv)) %>%
  filter(model %in% 1:56
         # model %in% c(22, 1, 8, 5, 19)
         ) %>%
  left_join(ems_grid, by = c("model" = "model_id")) %>% view()

# which model has the minimum mean absolute bias?
reshape2::melt(ssb_store) %>%
  drop_na() %>%
  left_join(ems_grid, by = c("Var3" = "model_id")) %>%
  filter(fish_blocks_type == 'none', srv_blocks_type == 'none') %>%
  group_by(Var3) %>%
  summarize(median_abs = (median(abs(value)))) %>% view()

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

ssb_sum_results %>%
  drop_na() %>%
  ggplot(aes(x = fishery_structure, y = median)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  facet_grid(survey_selex ~ fishery_selex, scales = 'free',
             labeller = labeller(
               survey_selex = \(x) paste("Survey:", x),
               fishery_selex = \(x) paste("Fishery:", x)
             )) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Fishery Fleet Structure', y = 'Spawning Stock Biomass Relative Error')

ssb_sum_results %>%
  filter(fishery_selex != 'All_Gamma' & survey_selex != 'All_Gamma') %>%
  ggplot(aes(x = fishery_structure, y = median)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  facet_grid(survey_selex~fishery_selex, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Fishery Fleet Structure', y = 'Spawning Stock Biomass Relative Error')

ssb_sum_results %>%
  filter((fishery_selex != 'All_Gamma' &
           survey_selex != 'All_Gamma'),
         fishery_structure != '1F_1T') %>%
  ggplot(aes(x = fishery_structure, y = median)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  facet_grid(survey_selex~fishery_selex, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Fishery Fleet Structure', y = 'Spawning Stock Biomass Relative Error')

ssb_sum_results %>%
  filter((fishery_selex != 'All_Gamma' &
            survey_selex != 'All_Gamma'),
         fishery_structure != '1F_1T',
         (fishery_selex != 'All_Logist' &
            survey_selex != 'All_Logist')) %>%
  ggplot(aes(x = fishery_structure, y = median)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  facet_grid(survey_selex~fishery_selex, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Fishery Fleet Structure', y = 'Spawning Stock Biomass Relative Error')

# 3F1T models seem to do better in estimaitng inittal pop scale
ssb_sum_results %>%
  filter(Var3 %in% c(22, 1, 8, 5, 19)) %>%
  left_join(ems_grid %>%
              filter(model_id %in% c(22, 1, 8, 5, 19)) %>%
              select(model_id, fishery_structure, fishery_selex, survey_selex),
            by = c("Var3" = 'model_id')) %>%
  mutate(model_name = paste(
    fishery_structure.x, paste("FS:", fishery_selex.x), paste("SS:", survey_selex.x),
    sep = '_'
  )) %>%
  ggplot(aes(x = Var1 + 1959, y = median, ymin = lwr, ymax = upr,
             fill = factor(model_name), color = factor(model_name))) +
  geom_ribbon(alpha = 0.1, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.25, 0.25)) +
  theme_bw(base_size = 15) +
  labs(x = 'Year',
       y = 'Spawning Stock Biomass Relative Error',
       color = 'Model', fill = 'Model')

# # Sensitivity Run Results ---------------------------------------------------------
# model_list_lowsamp <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_movesense_finalFAA.RDS"))
# model_list_lowsamp_francis <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_movesense_finalFAA_francis.RDS"))
#
# ssb_em_results <- array(NA, dim = c(2, length(1:y), n_sims))
# rec_em_results <- array(NA, dim = c(2, length(1:y), n_sims))
# ssb_om_results <- array(NA, dim = c(length(1:y), n_sims))
# rec_om_results <- array(NA, dim = c(length(1:y), n_sims))
#
# for(i in 1:n_sims) {
#   # get aggregated OM results
#   ssb_om_results[,i] <- colSums(lowsamp_movesense$SSB[,1:y,i])
#   rec_om_results[,i] <- colSums(lowsamp_movesense$Rec[,1:y,i])
#
#   # get EM results
#   ssb_em_results[1,,i] <- t(model_list_lowsamp[[i]]$rep$SSB)
#   rec_em_results[1,,i] <- t(model_list_lowsamp[[i]]$rep$Rec)
#   if(length(model_list_lowsamp_francis[[i]]) > 1) ssb_em_results[2,,i] <- t(model_list_lowsamp_francis[[i]]$rep$SSB)
#   if(length(model_list_lowsamp_francis[[i]]) > 1) rec_em_results[2,,i] <- t(model_list_lowsamp_francis[[i]]$rep$Rec)
# } # end i
#
#
# par(mfrow = c(2,2))
# # SSB
# plot(apply((ssb_em_results[1,,] - ssb_om_results) / ssb_om_results, 1, median, na.rm = T), type = 'l', ylim = c(-0.25, 0.25),
#      col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in SSB')
# lines(apply((ssb_em_results[2,,] - ssb_om_results) / ssb_om_results, 1, median, na.rm = T), col = 'red', lwd = 3, lty = 2)
# abline(h = 0, lty = 2, lwd = 3)
# # SSB abs
# plot(apply(abs((ssb_em_results[1,,] - ssb_om_results) / ssb_om_results), 1, median, na.rm = T), type = 'l', ylim = c(0, 0.25),
#      col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in SSB')
# lines(apply(abs((ssb_em_results[2,,] - ssb_om_results) / ssb_om_results), 1, median, na.rm = T), col = 'red', lwd = 3, lty = 2)
#
# # Rec
# plot(apply((rec_em_results[1,,] - rec_om_results) / rec_om_results, 1, median, na.rm = T), type = 'l', ylim = c(-1, 1),
#      col = 'black', lwd = 3, xlab = 'Year', ylab = 'Relative Error in Rec')
# lines(apply((rec_em_results[2,,] - rec_om_results) / rec_om_results, 1, median, na.rm = T), col = 'red', lwd = 3, lty = 2)
# abline(h = 0, lty = 2, lwd = 3)
# # Rec abs
# plot(apply(abs((rec_em_results[1,,] - rec_om_results) / rec_om_results), 1, median, na.rm = T), type = 'l', ylim = c(0, 1),
#      col = 'black', lwd = 3, xlab = 'Year', ylab = 'Median Absolute Relative Error in Rec')
# lines(apply(abs((rec_em_results[2,,] - rec_om_results) / rec_om_results), 1, median, na.rm = T), col = 'red', lwd = 3, lty = 2)
#
# # Retrospective Results ---------------------------------------------------
# retro_nofrancis <- read_csv(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_finalized_cross_test_retro.csv"))
# retro_wfrancis <- read_csv(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_finalized_cross_test_retro_francis.csv"))
#
# # Summarize info
# retro_nofrancis_sum <- retro_nofrancis %>%
#   group_by(Region, Year, peel, Type) %>%
#   summarize(lwr_95 = quantile(value, 0.025),
#             upr_95 = quantile(value, 0.975),
#             value = median(value)) %>%
#   mutate(type = 'No Francis') %>%
#   ungroup()
#
#
# # Summarize info
# retro_wfrancis_sum <- retro_wfrancis %>%
#   group_by(Region, Year, peel, Type) %>%
#   summarize(lwr_95 = quantile(value, 0.025),
#             upr_95 = quantile(value, 0.975),
#             value = median(value)) %>%
#   mutate(type = 'W Francis') %>%
#   ungroup()
#
# # get relative difference
# retro_nofrancis_rd <- get_retrospective_relative_difference(retro_nofrancis_sum) %>% mutate(type = 'No Francis')
# retro_francis_rd <- get_retrospective_relative_difference(retro_wfrancis_sum)%>% mutate(type = 'W Francis')
#
# # get mohn's rho
# mohns_rho_nofrancis <- retro_nofrancis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
#   dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'No Francis')
# mohns_rho_wfrancis <- retro_francis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
#   dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'W Francis')
#
# # bind
# retro_sgl_region <- rbind(retro_wfrancis_sum, retro_nofrancis_sum)
# mohns_rho <- rbind(mohns_rho_nofrancis, mohns_rho_wfrancis)
# ret_df <- rbind(retro_nofrancis_rd, retro_francis_rd)
#
# # Generally look the same
# ggplot() +
#   geom_hline(yintercept = 0,
#              lty = 2,
#              lwd = 1.3) +
#   geom_line(
#     ret_df,
#     mapping = aes(
#       x = Year,
#       y = rd,
#       group = as.numeric(peel),
#       color = as.numeric(peel)
#     ),
#     lwd = 1.3
#   ) +
#   geom_point(
#     ret_df %>% dplyr::filter(peel ==
#                                max(Year) - Year),
#     mapping = aes(
#       x = Year,
#       y = rd,
#       group = as.numeric(peel),
#       fill = as.numeric(peel)
#     ),
#     pch = 21,
#     size = 6
#   ) +
#   geom_text(
#     mohns_rho,
#     mapping = aes(
#       x = -Inf,
#       y = Inf,
#       label = paste("Mohns Rho:", round(rho, 4))
#     ),
#     hjust = -0.3,
#     vjust = 3,
#     size = 5
#   )+
#   guides(color = guide_colourbar(barwidth = 15, barheight = 1.3)) +
#   labs(x = "Year",
#        y = "Relative Difference from Terminal Year",
#        color = "Retrospective Year",
#        fill = "Retrospective Year") +
#   scale_color_viridis_c() +
#   scale_fill_viridis_c() +
#   facet_grid(type ~ Type, scales = "free") +
#   theme_sablefish() +
#   theme(legend.position = "top")
#
# ggplot() +
#   geom_line(
#     retro_sgl_region %>%
#       filter(peel != 0),
#     mapping = aes(
#       x = Year,
#       y = value,
#       group = peel,
#       color = peel
#     ),
#     lwd = 1
#   ) +
#   geom_line(
#     retro_sgl_region %>%
#       filter(peel == 0),
#     mapping = aes(x = Year, y = value),
#     lty = 2,
#     lwd = 1
#   ) +
#   geom_ribbon(
#     retro_sgl_region %>%
#       filter(peel != 0),
#     mapping = aes(
#       x = Year,
#       y = value,
#       group = peel,
#       ymin = lwr_95,
#       ymax = upr_95,
#       fill = peel
#     ),
#     color = NA, alpha = 0.2
#   ) +
#   geom_text(mohns_rho,
#             mapping = aes(x = -Inf, y = Inf, label = paste("Mohns Rho:", round(rho, 4))),
#             hjust = -0.3, vjust = 3, size = 5) +
#   scale_color_viridis_c() +
#   scale_fill_viridis_c() +
#
#   coord_cartesian(ylim = c(0, NA)) +
#   facet_grid(Type ~ type, scales = "free_y") +
#   theme_sablefish() +
#   labs(x = "Year", y = "Value", color = "Peel")
#
#
#
# # Compare Retro to True Values --------------------------------------------
# # Input true values into peel == 0 for evaluation
# true_ssb <- colSums(lowsamp$SSB[,-66,1])
# true_rec <- colSums(lowsamp$Rec[,-66,1])
# true_values <- c(true_ssb, true_rec)
#
# retro_nofrancis[retro_nofrancis$peel == 0, ]$value <- rep(true_values, length(unique(retro_nofrancis$Sim)))
# retro_wfrancis[retro_wfrancis$peel == 0, ]$value <- rep(true_values, length(unique(retro_wfrancis$Sim)))
#
# # Summarize info
# retro_nofrancis_sum <- retro_nofrancis %>%
#   group_by(Region, Year, peel, Type) %>%
#   summarize(lwr_95 = quantile(value, 0.025),
#             upr_95 = quantile(value, 0.975),
#             value = median(value)) %>%
#   mutate(type = 'No Francis') %>%
#   ungroup()
#
#
# # Summarize info
# retro_wfrancis_sum <- retro_wfrancis %>%
#   group_by(Region, Year, peel, Type) %>%
#   summarize(lwr_95 = quantile(value, 0.025),
#             upr_95 = quantile(value, 0.975),
#             value = median(value)) %>%
#   mutate(type = 'W Francis') %>%
#   ungroup()
#
# # get relative difference
# retro_nofrancis_rd <- get_retrospective_relative_difference(retro_nofrancis_sum) %>% mutate(type = 'No Francis')
# retro_francis_rd <- get_retrospective_relative_difference(retro_wfrancis_sum)%>% mutate(type = 'W Francis')
#
# # get mohn's rho
# mohns_rho_nofrancis <- retro_nofrancis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
#   dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'No Francis')
# mohns_rho_wfrancis <- retro_francis_rd %>% dplyr::filter(ifelse(Type == 'Recruitment', peel == (max(Year) - 1 -  Year), peel == (max(Year) -  Year))) %>%
#   dplyr::group_by(Type, Region) %>% dplyr::summarize(rho = mean(rd)) %>% mutate(type = 'W Francis')
#
# # bind
# retro_sgl_region <- rbind(retro_wfrancis_sum, retro_nofrancis_sum)
# mohns_rho <- rbind(mohns_rho_nofrancis, mohns_rho_wfrancis)
# ret_df <- rbind(retro_nofrancis_rd, retro_francis_rd)
# # Generally look the same
# ggplot() +
#   geom_hline(yintercept = 0,
#              lty = 2,
#              lwd = 1.3) +
#   geom_line(
#     ret_df,
#     mapping = aes(
#       x = Year,
#       y = rd,
#       group = as.numeric(peel),
#       color = as.numeric(peel)
#     ),
#     lwd = 1.3
#   ) +
#   geom_point(
#     ret_df %>% dplyr::filter(peel ==
#                                max(Year) - Year),
#     mapping = aes(
#       x = Year,
#       y = rd,
#       group = as.numeric(peel),
#       fill = as.numeric(peel)
#     ),
#     pch = 21,
#     size = 6
#   ) +
#   geom_text(
#     mohns_rho,
#     mapping = aes(
#       x = -Inf,
#       y = Inf,
#       label = paste("Mohns Rho:", round(rho, 4))
#     ),
#     hjust = -0.3,
#     vjust = 3,
#     size = 5
#   )+
#   guides(color = guide_colourbar(barwidth = 15, barheight = 1.3)) +
#   labs(x = "Year",
#        y = "Relative Difference from Terminal Year",
#        color = "Retrospective Year",
#        fill = "Retrospective Year") +
#   scale_color_viridis_c() +
#   scale_fill_viridis_c() +
#   facet_grid(type ~ Type, scales = "free") +
#   theme_sablefish() +
#   theme(legend.position = "top")
#
# ggplot() +
#   geom_line(
#     retro_sgl_region %>%
#       filter(peel != 0),
#     mapping = aes(
#       x = Year,
#       y = value,
#       group = peel,
#       color = peel
#     ),
#     lwd = 1
#   ) +
#   geom_line(
#     retro_sgl_region %>%
#       filter(peel == 0),
#     mapping = aes(x = Year, y = value),
#     lty = 2,
#     lwd = 1
#   ) +
#   geom_ribbon(
#     retro_sgl_region %>%
#       filter(peel != 0),
#     mapping = aes(
#       x = Year,
#       y = value,
#       group = peel,
#       ymin = lwr_95,
#       ymax = upr_95,
#       fill = peel
#     ),
#     color = NA, alpha = 0.2
#   ) +
#   geom_text(mohns_rho,
#             mapping = aes(x = -Inf, y = Inf, label = paste("Mohns Rho:", round(rho, 4))),
#             hjust = -0.3, vjust = 3, size = 5) +
#   scale_color_viridis_c() +
#   scale_fill_viridis_c() +
#   coord_cartesian(ylim = c(0, NA)) +
#   facet_grid(Type ~ type, scales = "free_y") +
#   theme_sablefish() +
#   labs(x = "Year", y = "Value", color = "Peel")
#
#
#
