# Purpose: To self test the FAA EM prior to the MSE year
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
lowsamp <- readRDS(here("outputs", 'cross_test', "spatial_blocked_scenarios", "spt_rand_OM_lowsamp.RDS"))

# Simulation dimensions
y <- 65 # terminal year prior to MSE (no feedback)
n_sims <- 100 # number of sims

# Setup EMs ---------------------------------------------------------------
ems_grid <- expand.grid(

  # Fishery structure
  fishery_structure = c(
    "Single_FixedGear_SingleTrawl", # single fleet aggregated across all areas
    "Three_FixedGear_SingleTrawl", # three fixed-gear fleets, trawl aggregated across all areas
    "Three_FixedGear_Three_Trawl" # three fixed-gear, three trawl gear
  ),
  # Fishery selex options
  fishery_selex = c(
    "FixedGear_Logistic_Trawl_Gamma", # fixed-gear logistic, trawl gamma (only applicable for single fleet)
    "All_FixedGear_Logistic_Trawl_Gamma", # multi FAA, all logistic, trawl gamma
    "BS_FixedGear_Gamma_Others_Logistic_Trawl_Gamma", # BS Gamma, all other fixed gears are logistic, gamma trawl
    "BS_AI_FixedGear_Gamma_GOA_Logistic_Trawl_Gamma", # BS, AI Gamma, GOA Logistic, Gamma trawl
    "All_Gamma" # all gamma
  ),
  # Survey block options
  survey_blocks = c(
    "NoBlocks_Domestic", # no blocks for domestic fleet
    "TwoBlocks_Domestic" # two blocks for domestic fleet
  ),
  # survey selex options
  survey_selex = c(
    "All_Logistic", # all logistic selex
    "BS_Gamma_Others_Logistic", # bs domestic gamma, all others logistic
    "BS_AI_Gamma_Others_Logistic", # bs ai gamma, others logistic
    "BS_AI_GOA_Gamma_JP_Logistic" # everything in domestic is gamma, jp is logistic
  ),
  stringsAsFactors = FALSE
)

# Remove invalid combinations:
# Structure 1 (Single FixedGear) can only use: FixedGear_Logistic_Trawl_Gamma
ems_grid <- ems_grid[!(ems_grid$fishery_structure == "Single_FixedGear_SingleTrawl" &
                         ems_grid$fishery_selex %in% c("All_FixedGear_Logistic_Trawl_Gamma",
                                                       "BS_FixedGear_Gamma_Others_Logistic_Trawl_Gamma",
                                                       "BS_AI_FixedGear_Gamma_GOA_Logistic_Trawl_Gamma",
                                                       "All_Gamma")), ]

# Structure 2 and 3 don't use FixedGear_Logistic_Trawl_Gamma
ems_grid <- ems_grid[!(ems_grid$fishery_structure != "Single_FixedGear_SingleTrawl" &
                         ems_grid$fishery_selex == ("FixedGear_Logistic_Trawl_Gamma")), ]

ems_grid$model_id <- 1:nrow(ems_grid)
ems_grid <- ems_grid[, c("model_id", "fishery_structure", "fishery_selex", "survey_blocks", "survey_selex")]

# Get fishery structure specs
get_fish_structure <- function(struct_name) {
  specs <- list(
    "Single_FixedGear_SingleTrawl" = list(
      n_fish_fleets = 2,
      n_srv_fleets = 4,
      fish_blocks = c("Block_1_Year_1-35_Fleet_1", "Block_2_Year_36-56_Fleet_1",
                      "Block_3_Year_57-terminal_Fleet_1", "none_Fleet_2"),
      fish_prior_blocks = list("1" = 1:3, "2" = 1)
    ),
    "Three_FixedGear_SingleTrawl" = list(
      n_fish_fleets = 4,
      n_srv_fleets = 4,
      fish_blocks = c("Block_1_Year_1-35_Fleet_1", "Block_2_Year_36-terminal_Fleet_1",
                      "Block_1_Year_1-35_Fleet_2", "Block_2_Year_36-terminal_Fleet_2",
                      "Block_1_Year_1-35_Fleet_3", "Block_2_Year_36-56_Fleet_3",
                      "Block_3_Year_57-terminal_Fleet_3", "none_Fleet_4"),
      fish_prior_blocks = list("1" = 1:2, "2" = 1:2, "3" = 1:3, "4" = 1)
    ),
    "Three_FixedGear_Three_Trawl" = list(
      n_fish_fleets = 6,
      n_srv_fleets = 4,
      fish_blocks = c("Block_1_Year_1-35_Fleet_1", "Block_2_Year_36-terminal_Fleet_1",
                      "Block_1_Year_1-35_Fleet_2", "Block_2_Year_36-terminal_Fleet_2",
                      "Block_1_Year_1-35_Fleet_3", "Block_2_Year_36-56_Fleet_3",
                      "Block_3_Year_57-terminal_Fleet_3", "none_Fleet_4", "none_Fleet_5", "none_Fleet_6"),
      fish_prior_blocks = list("1" = 1:2, "2" = 1:2, "3" = 1:3, "4" = 1, "5" = 1, "6" = 1)
    )
  )
  return(specs[[struct_name]])
}

# Get fishery selectivity models
get_fish_selex <- function(selex_name, n_fleets) {
  specs <- list(
    # For Single FixedGear structure (2 fleets: 1 fixed, 1 trawl)
    "FixedGear_Logistic_Trawl_Gamma" = list(
      "2" = c("logist1_Fleet_1", "gamma_Fleet_2")
    ),
    "All_Gamma" = list(
      "4" = c("gamma_Fleet_1", "gamma_Fleet_2", "gamma_Fleet_3", "gamma_Fleet_4"),
      "6" = c("gamma_Fleet_1", "gamma_Fleet_2", "gamma_Fleet_3",
              "gamma_Fleet_4", "gamma_Fleet_5", "gamma_Fleet_6")
    ),
    # For Three FixedGear structures (4 or 6 fleets: 3 fixed regions + 1 or 3 trawl)
    "All_FixedGear_Logistic_Trawl_Gamma" = list(
      "4" = c("logist1_Fleet_1", "logist1_Fleet_2", "logist1_Fleet_3", "gamma_Fleet_4"),
      "6" = c("logist1_Fleet_1", "logist1_Fleet_2", "logist1_Fleet_3",
              "gamma_Fleet_4", "gamma_Fleet_5", "gamma_Fleet_6")
    ),
    "BS_FixedGear_Gamma_Others_Logistic_Trawl_Gamma" = list(
      "4" = c("gamma_Fleet_1", "logist1_Fleet_2", "logist1_Fleet_3", "gamma_Fleet_4"),
      "6" = c("gamma_Fleet_1", "logist1_Fleet_2", "logist1_Fleet_3",
              "gamma_Fleet_4", "gamma_Fleet_5", "gamma_Fleet_6")
    ),
    "BS_AI_FixedGear_Gamma_GOA_Logistic_Trawl_Gamma" = list(
      "4" = c("gamma_Fleet_1", "gamma_Fleet_2", "logist1_Fleet_3", "gamma_Fleet_4"),
      "6" = c("gamma_Fleet_1", "gamma_Fleet_2", "logist1_Fleet_3",
              "gamma_Fleet_4", "gamma_Fleet_5", "gamma_Fleet_6")
    )
  )
  return(specs[[selex_name]][[as.character(n_fleets)]])
}

# Get survey blocks
get_survey_blocks <- function(blocks_name) {
  specs <- list(
    "NoBlocks_Domestic" = list(
      blocks = c("none_Fleet_1", "none_Fleet_2", "none_Fleet_3", "none_Fleet_4"),
      n_blocks = 1
    ),
    "TwoBlocks_Domestic" = list(
      blocks = c("Block_1_Year_1-56_Fleet_1", "Block_2_Year_57-terminal_Fleet_1",
                 "Block_1_Year_1-56_Fleet_2", "Block_2_Year_57-terminal_Fleet_2",
                 "Block_1_Year_1-56_Fleet_3", "Block_2_Year_57-terminal_Fleet_3",
                 "none_Fleet_4"),
      n_blocks = 2
    )
  )
  return(specs[[blocks_name]])
}

# Get survey selectivity models
get_survey_selex <- function(selex_name) {
  specs <- list(
    "All_Logistic" = c("logist1_Fleet_1", "logist1_Fleet_2", "logist1_Fleet_3", "logist1_Fleet_4"),
    "BS_Gamma_Others_Logistic" = c("gamma_Fleet_1", "logist1_Fleet_2", "logist1_Fleet_3", "logist1_Fleet_4"),
    "BS_AI_Gamma_Others_Logistic" = c("gamma_Fleet_1", "gamma_Fleet_2", "logist1_Fleet_3", "logist1_Fleet_4"),
    "BS_AI_GOA_Gamma_JP_Logistic" = c("gamma_Fleet_1", "gamma_Fleet_2", "gamma_Fleet_3", "logist1_Fleet_4")
  )
  return(specs[[selex_name]])
}

ems <- list()

for (i in 1:nrow(ems_grid)) {

  # Get specs for this row
  fish_struct <- get_fish_structure(ems_grid$fishery_structure[i])
  fish_selex <- get_fish_selex(ems_grid$fishery_selex[i], fish_struct$n_fish_fleets)
  surv_blocks <- get_survey_blocks(ems_grid$survey_blocks[i])
  surv_selex <- get_survey_selex(ems_grid$survey_selex[i])

  # Build the model
  ems[[i]] <- list(
    faa_n_fish_fleets = fish_struct$n_fish_fleets,
    faa_n_srv_fleets = fish_struct$n_srv_fleets,

    # Fishery selectivity
    fish_sel_blocks = fish_struct$fish_blocks,
    fish_sel_model = fish_selex,
    fish_fixed_sel_pars_spec = rep("est_all", fish_struct$n_fish_fleets),
    fish_selex_prior = build_selex_prior(
      fleets_blocks = fish_struct$fish_prior_blocks,
      sex_par = expand.grid(sex = 1:2, par = 1:2),
      region = 1,
      mu = 1,
      sd = 3
    ),

    # Survey selectivity
    srv_sel_blocks = surv_blocks$blocks,
    srv_sel_model = surv_selex,
    srv_fixed_sel_pars_spec = rep("est_all", 4),
    srv_selex_prior = build_selex_prior(
      fleets_blocks = list(
        "1" = 1:surv_blocks$n_blocks,
        "2" = 1:surv_blocks$n_blocks,
        "3" = 1:surv_blocks$n_blocks,
        "4" = 1
      ),
      sex_par = expand.grid(sex = 1:2, par = 1:2),
      region = 1,
      mu = 1,
      sd = 3
    )
  )
}


# Run Models (First Half) --------------------------------------------------------------

handlers(global = TRUE)  # progress bar
plan(multisession, workers = 6)
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
                          fish_sel_blocks = ems[[j]]$fish_sel_blocks,
                          fish_sel_model = ems[[j]]$fish_sel_model,
                          fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                          fish_selex_prior = ems[[j]]$fish_selex_prior,
                          srv_sel_blocks = ems[[j]]$srv_sel_blocks,
                          srv_sel_model = ems[[j]]$srv_sel_model,
                          srv_fixed_sel_pars_spec = ems[[j]]$srv_fixed_sel_pars_spec,
                          srv_selex_prior = ems[[j]]$srv_selex_prior
      )

      asmt_list$par$ln_fish_fixed_sel_pars[,,,,1] <- log(4)
      asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

      # fit model
      model <- tryCatch({
        # fit model
        model <- fit_model(
          asmt_list$data,
          asmt_list$par,
          asmt_list$map,
          NULL,
          2,
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
saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_blocked_scenarios", "faa_crosstest_lowsamp_1.RDS"))

# Run Models (Second Half) --------------------------------------------------------------

handlers(global = TRUE)  # progress bar
plan(multisession, workers = 5)
options(future.globals.maxSize = 20e9)

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
                          fish_sel_blocks = ems[[j]]$fish_sel_blocks,
                          fish_sel_model = ems[[j]]$fish_sel_model,
                          fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
                          fish_selex_prior = ems[[j]]$fish_selex_prior,
                          srv_sel_blocks = ems[[j]]$srv_sel_blocks,
                          srv_sel_model = ems[[j]]$srv_sel_model,
                          srv_fixed_sel_pars_spec = ems[[j]]$srv_fixed_sel_pars_spec,
                          srv_selex_prior = ems[[j]]$srv_selex_prior
      )

      asmt_list$par$ln_fish_fixed_sel_pars[,,,,1] <- log(4)
      asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

      # fit model
      model <- tryCatch({
        # fit model
        model <- fit_model(
          asmt_list$data,
          asmt_list$par,
          asmt_list$map,
          NULL,
          2,
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
saveRDS(model_list_lowsamp, here("outputs", "cross_test", "spatial_blocked_scenarios", "faa_crosstest_lowsamp_2.RDS"))


# Visualize Results -------------------------------------------------------

# Read in models
model_list_lowsamp1 <- readRDS(here("outputs", "cross_test", "spatial_blocked_scenarios", "faa_crosstest_lowsamp_1.RDS"))
model_list_lowsamp2 <- readRDS(here("outputs", "cross_test", "spatial_blocked_scenarios", "faa_crosstest_lowsamp_2.RDS"))
model_list_lowsamp_all <- vector('list', 100)

# Combine models
for(i in 1:100) {
  for(j in 1:nrow(ems_grid)) {
    if(j <= (nrow(ems_grid) / 2)) model_list_lowsamp_all[[i]][[j]] <- model_list_lowsamp1[[i]][[j]]
    else model_list_lowsamp_all[[i]][[j]] <- model_list_lowsamp2[[i]][[j]]
  }
}

# Store models
ssb_store <- array(NA, c(65, 100, nrow(ems_grid)))
rec_store <- array(NA, c(65, 100, nrow(ems_grid)))
for(mod in 1:72) {
  for(i in 1:100) {
    if(length(model_list_lowsamp_all[[i]][[mod]]) > 1) {
      ssb_store[,i,mod] <- (t(model_list_lowsamp_all[[i]][[mod]]$rep$SSB) - colSums(lowsamp$SSB[,-66,1])) / colSums(lowsamp$SSB[,-66,1])
      rec_store[,i,mod] <- (t(model_list_lowsamp_all[[i]][[mod]]$rep$Rec) - colSums(lowsamp$Rec[,-66,1])) / colSums(lowsamp$Rec[,-66,1])
    }
  }
}

ems_grid %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl|Three_FixedGear_Three_Trawl'),
         !str_detect(fishery_selex, 'All_Gamma|All_FixedGear_Logistic_Trawl_Gamma'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(survey_blocks, 'NoBlocks_Domestic')) %>% view()

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

ssb_sum_results %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

ssb_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

ssb_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

ssb_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

ssb_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

ssb_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma|All_FixedGear_Logistic_Trawl_Gamma')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

ssb_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma|All_FixedGear_Logistic_Trawl_Gamma'),
         !str_detect(survey_blocks, 'NoBlocks')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

ssb_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma|All_FixedGear_Logistic_Trawl_Gamma'),
         !str_detect(survey_blocks, 'NoBlocks')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr,
             group = Var3, color = factor(Var3), fill = factor(Var3))) +
  geom_line(lwd = 1.3) +
  # geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  # facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')

ssb_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl|Three_FixedGear_Three_Trawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma|All_FixedGear_Logistic_Trawl_Gamma'),
         !str_detect(survey_blocks, 'NoBlocks')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr,
             group = Var3, color = factor(Var3), fill = factor(Var3))) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_grid(survey_selex~fishery_selex, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'SSB RE')


rec_sum_results %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

rec_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

rec_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

rec_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

rec_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

rec_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma|All_FixedGear_Logistic_Trawl_Gamma')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

rec_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma|All_FixedGear_Logistic_Trawl_Gamma'),
         !str_detect(survey_blocks, 'NoBlocks')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr, group = Var3)) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

rec_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma|All_FixedGear_Logistic_Trawl_Gamma'),
         !str_detect(survey_blocks, 'NoBlocks')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr,
             group = Var3, color = factor(Var3), fill = factor(Var3))) +
  geom_line(lwd = 1.3) +
  # geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  # facet_wrap(~Var3, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

rec_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl|Three_FixedGear_Three_Trawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma|All_FixedGear_Logistic_Trawl_Gamma'),
         !str_detect(survey_blocks, 'NoBlocks')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr,
             group = Var3, color = factor(Var3), fill = factor(Var3))) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_grid(survey_selex~fishery_selex, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')

rec_sum_results %>%
  filter(!str_detect(fishery_structure, 'Single_FixedGear_SingleTrawl|Three_FixedGear_Three_Trawl'),
         !str_detect(survey_selex, 'All_Logistic|BS_AI_GOA_Gamma_JP_Logistic'),
         !str_detect(fishery_selex, 'All_Gamma|All_FixedGear_Logistic_Trawl_Gamma'),
         !str_detect(survey_blocks, 'NoBlocks')) %>%
  ggplot(aes(x = Var1, y = median, ymin = lwr, ymax = upr,
             group = Var3, color = factor(Var3), fill = factor(Var3))) +
  geom_line(lwd = 1.3) +
  geom_ribbon(alpha = 0.35, color = NA) +
  geom_hline(yintercept = 0, lty = 2, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  facet_grid(survey_selex~fishery_selex, scales = 'free') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'top') +
  labs(x = 'Yr', y = 'Rec RE')



# lines(colSums(lowsamp$Rec[,-66,sim]), lty = 2, col = 'red', lwd = 3)
# lines(colSums(lowsamp$SSB[,-66,sim]), lty = 2, col = 'red', lwd = 3)
#
# plot(store[,,3])
#
# #
#
# mod1 <- 5
# mod2 <- 51
# par(mfrow = c(2,2))
# for(i in 1:100) {
#   if(i == 1) plot(model_list_lowsamp_all[[i]][[mod1]]$rep$srv_sel[1,60,,1,1], type = 'l', ylim = c(0,1), main = 'BS Fleet (Model 49)', ylab = 'sel')
#   else lines(model_list_lowsamp_all[[i]][[mod1]]$rep$srv_sel[1,60,,1,1], type = 'l')
# }
#
# for(i in 1:100) {
#   if(i == 1) plot(model_list_lowsamp_all[[i]][[mod1]]$rep$srv_sel[1,60,,1,2], type = 'l', ylim = c(0,1), main = 'AI Fleet (Model 49)', ylab = 'sel')
#   else lines(model_list_lowsamp_all[[i]][[mod1]]$rep$srv_sel[1,60,,1,2], type = 'l')
# }
#
# for(i in 1:100) {
#   if(i == 1) plot(model_list_lowsamp_all[[i]][[mod2]]$rep$srv_sel[1,60,,1,1], type = 'l', ylim = c(0,1), main = 'BS Fleet (Model 51)', ylab = 'sel')
#   else lines(model_list_lowsamp_all[[i]][[mod2]]$rep$srv_sel[1,60,,1,1], type = 'l')
# }
#
# for(i in 1:100) {
#   if(i == 1) plot(model_list_lowsamp_all[[i]][[mod2]]$rep$srv_sel[1,60,,1,2], type = 'l', ylim = c(0,1), main = 'AI Fleet (Model 51)', ylab = 'sel')
#   else lines(model_list_lowsamp_all[[i]][[mod2]]$rep$srv_sel[1,60,,1,2], type = 'l')
# }
#
#
# mod1 <- 3
# mod2 <- 51
# plot(apply(model_list_lowsamp_all[[1]][[mod1]]$rep$SrvIAA[,,,,1], 2, sum))
# lines(apply(model_list_lowsamp_all[[1]][[mod2]]$rep$SrvIAA[,,,,1], 2, sum))
#
# plot(apply(model_list_lowsamp_all[[1]][[mod1]]$rep$SrvIAA[,,,,2], 2, sum))
# lines(apply(model_list_lowsamp_all[[1]][[mod2]]$rep$SrvIAA[,,,,2], 2, sum))
#
#
# plot(apply(model_list_lowsamp_all[[1]][[mod1]]$rep$SrvIAA[,,,,3], 2, sum))
# lines(apply(model_list_lowsamp_all[[1]][[mod2]]$rep$SrvIAA[,,,,3], 2, sum))

#
# j <- 49
# # get faa data
# asmt_list <- faa_em(sim_env = lowsamp,
#                     y = y,
#                     sim = i,
#                     srv_idx_se = 0.2,
#                     age_lag = 1,
#                     lls_design_type = 'historical',
#                     faa_n_fish_fleets = ems[[j]]$faa_n_fish_fleets,
#                     faa_n_srv_fleets = ems[[j]]$faa_n_srv_fleets,
#                     srv_wgt = 'numbers',
#                     fish_wgt = 'numbers',
#                     fish_sel_blocks = ems[[j]]$fish_sel_blocks,
#                     fish_sel_model = ems[[j]]$fish_sel_model,
#                     fish_fixed_sel_pars_spec = ems[[j]]$fish_fixed_sel_pars_spec,
#                     fish_selex_prior = ems[[j]]$fish_selex_prior,
#                     srv_sel_blocks = ems[[j]]$srv_sel_blocks,
#                     srv_sel_model = ems[[j]]$srv_sel_model,
#                     srv_fixed_sel_pars_spec = ems[[j]]$srv_fixed_sel_pars_spec,
#                     srv_selex_prior = ems[[j]]$srv_selex_prior
# )
#
# retros <- do_retrospective(
#   8, data = model_list_lowsamp_all[[i]][[51]]$data,
#   parameters = asmt_list$par,
#   mapping = asmt_list$map,
#   random = NULL, do_par = F, newton_loops = 3, do_francis = T, n_francis_iter = 5
# )
#
# get_retrospective_plot(retros, 2)
