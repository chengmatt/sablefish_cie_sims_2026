  # Purpose: To condition a 5-region OM for sablefish
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
  closed_loop_yrs <- 1      # Years to project forward
  n_years <- length(data$years)  # number of years
  burnin_years <- 1:n_years  # Historical conditioning period
  n_sims <- 1               # Number of replicate simulations
  assess_freq <- 1           # Assessment frequency
  data_yr_freq <- 1          # Data collection frequency
  n_regions <- 5             # number of regions
  n_fish_fleets <- 2         # number of fishery fleets
  n_srv_fleets <- 3          # number of survey fleets

  # get fleet allocation by region
  fleet_allocation <- array(NA, dim = c(n_regions, n_fish_fleets))
  fleet_allocation[,1] <- c(0.5,0.75,0.8,0.8,0.95) # from fmp
  fleet_allocation[,2] <- 1 - fleet_allocation[,1] # trawl gear allocation

  # Global SPR
  global_spr <- Get_Reference_Points(data = om_values$data,
                                     rep = om_values$rep,
                                     SPR_x = 0.4,
                                     type = 'multi_region',
                                     what = 'global_SPR',
                                     calc_rec_st_yr = 20,
                                     rec_age = 2
                                     )

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
    ISS_FishAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
    ISS_FishLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_fish_fleets, n_sims)),
    ISS_SrvAgeComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
    ISS_SrvLenComps = array(20, dim = c(n_regions, n_years + closed_loop_yrs, om_values$data$n_sexes, n_srv_fleets, n_sims)),
    ObsFishIdx_SE = array(NA, dim = c(n_regions, n_years + closed_loop_yrs, n_fish_fleets)),
    ObsSrvIdx_SE = array(0.2, dim = c(n_regions, n_years + closed_loop_yrs, n_srv_fleets))
  )

  # Setup Closed Loop -------------------------------------------------------

  # Single-region, current design
  sim_env_current <- Setup_sim_env(sim_list_rand)
  sim_env_current <- add_aggregated_obj_to_simenv(sim_env_current)
  sim_env_current <- run_single_rg_closedloop_parallel(sim_env_current, n_sims, fleet_allocation, "current", 8)
  saveRDS(sim_env, here("outputs", "mse_results", "single_region_rand_current.RDS"))

  # Single-region, historical design
  sim_env_hist <- Setup_sim_env(sim_list_rand)
  sim_env_hist <- add_aggregated_obj_to_simenv(sim_env_hist)
  sim_env_hist <- run_single_rg_closedloop_parallel(sim_env_hist, n_sims, fleet_allocation, "historical", 8)
  saveRDS(sim_env_hist, here("outputs", "mse_results", "single_region_rand_hist.RDS"))

  # Single-region, sample all regions
  sim_env_all <- Setup_sim_env(sim_list_rand)
  sim_env_all <- add_aggregated_obj_to_simenv(sim_env_all)
  sim_env_all <- run_single_rg_closedloop_parallel(sim_env = sim_env_all, n_sims, fleet_allocation, lls_design_type = "all", 8)
  saveRDS(sim_env_all, here("outputs", "mse_results", "single_region_rand_all.RDS"))

  sim_env_current <- readRDS( here("outputs", "mse_results", "single_region_rand_current.RDS"))
  sim_env_hist <- readRDS( here("outputs", "mse_results", "single_region_rand_all.RDS"))
  sim_env_all <- readRDS( here("outputs", "mse_results", "single_region_highregime_current.RDS"))

  scenarios <- list(sim_env_current, sim_env_hist, sim_env_all)
  res <- array(NA, c(3, 95,100))
  for(j in 1:length(scenarios)) {
    for(i in 1:100) {
      res[j,,i] <- (t(scenarios[[j]]$models[[i]]$rep$SSB) -
                    colSums(scenarios[[j]]$SSB[,,i]) ) / colSums(scenarios[[j]]$SSB[,,i])

    }
  }

  reshape2::melt(res) %>%
    group_by(Var2, Var1) %>%
    summarize(value = median(value)) %>%
    ggplot(aes(x = Var2, y = value, color = factor(Var1))) +
    geom_line()

  lines(apply(res, 1, median), col = 'red', lty = 2, lwd = 3)
  abline(h = 0, lty = 2, lwd = 3, col = 'blue')
  plot(  apply(res, 1, median)[-c(1:65)], col = 'red', lty = 2, lwd = 3)

  reshape2::melt(sim_env_hist$SSB) %>%
    mutate(period = ifelse(Var2 <= 65, 'H', 'P')) %>%
    filter(period != 'H' | Var3 == 1, value != 0) %>%
    group_by(Var1, Var2, period) %>%
    summarize(median = median(value),
              lwr = quantile(value, 0.025),
              upr = quantile(value, 0.975)) %>%
    ggplot(aes(x = Var2, y = median,  ymin = lwr, ymax = upr)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, color = NA) +
    facet_wrap(~paste("Region", Var1), nrow = 1, scales = 'free') +
    geom_vline(xintercept = 65, lty = 2) +
    labs(x = 'Year', y = 'SSB') +
    ylim(0,NA) +
    theme_bw(base_size = 20) +
    theme(legend.position = 'none')

  reshape2::melt(sim_env_current$SSB) %>%
    mutate(period = ifelse(Var2 <= 65, 'H', 'P')) %>%
    filter(period != 'H' | Var3 == 1, value != 0) %>%
    group_by(Var1, Var2, period) %>%
    summarize(median = median(value),
              lwr = quantile(value, 0.025),
              upr = quantile(value, 0.975)) %>%
    ggplot(aes(x = Var2, y = median,  ymin = lwr, ymax = upr)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, color = NA) +
    facet_wrap(~paste("Region", Var1), nrow = 1, scales = 'free') +
    geom_vline(xintercept = 65, lty = 2) +
    labs(x = 'Year', y = 'SSB') +
    ylim(0,NA) +
    theme_bw(base_size = 20) +
    theme(legend.position = 'none')

  reshape2::melt(sim_env_hist$SSB) %>%
    group_by(Var2, Var3) %>%
    mutate(value = sum(value)) %>%
    mutate(period = ifelse(Var2 <= 65, 'H', 'P')) %>%
    filter(period != 'H' | Var3 == 1, value != 0) %>%
    group_by(Var2, period) %>%
    summarize(median = median(value),
              lwr = quantile(value, 0.025),
              upr = quantile(value, 0.975)) %>%
    ggplot(aes(x = Var2, y = median, ymin = lwr, ymax = upr)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, color = NA) +
    geom_vline(xintercept = 65, lty = 2) +
    # geom_hline(yintercept = sum(global_spr$b_ref_pt)) +
    labs(x = 'Year', y = 'SSB') +
    theme_bw(base_size = 20) +
    ylim(0, NA) +
    theme(legend.position = 'none')

  reshape2::melt(sim_env_current$SSB) %>%
    group_by(Var2, Var3) %>%
    mutate(value = sum(value)) %>%
    mutate(period = ifelse(Var2 <= 65, 'H', 'P')) %>%
    filter(period != 'H' | Var3 == 1, value != 0) %>%
    group_by(Var2, period) %>%
    summarize(median = median(value),
              lwr = quantile(value, 0.025),
              upr = quantile(value, 0.975)) %>%
    ggplot(aes(x = Var2, y = median, ymin = lwr, ymax = upr)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, color = NA) +
    geom_vline(xintercept = 65, lty = 2) +
    # geom_hline(yintercept = sum(global_spr$b_ref_pt)) +
    labs(x = 'Year', y = 'SSB') +
    theme_bw(base_size = 20) +
    ylim(0, NA) +
    theme(legend.position = 'none')

  reshape2::melt(sim_env_all$SSB) %>%
    group_by(Var2, Var3) %>%
    mutate(value = sum(value)) %>%
    mutate(period = ifelse(Var2 <= 65, 'H', 'P')) %>%
    filter(period != 'H' | Var3 == 1, value != 0) %>%
    group_by(Var2, period) %>%
    summarize(median = median(value),
              lwr = quantile(value, 0.025),
              upr = quantile(value, 0.975)) %>%
    ggplot(aes(x = Var2, y = median, ymin = lwr, ymax = upr)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, color = NA) +
    geom_vline(xintercept = 65, lty = 2) +
    # geom_hline(yintercept = sum(global_spr$b_ref_pt)) +
    labs(x = 'Year', y = 'SSB') +
    theme_bw(base_size = 20) +
    ylim(0, NA) +
    theme(legend.position = 'none')

  reshape2::melt(sim_env_all$Rec) %>%
    mutate(period = ifelse(Var2 <= 65, 'H', 'P')) %>%
    filter(period != 'H' | Var3 == 1, value != 0) %>%
    group_by(Var1, Var2, period) %>%
    summarize(median = median(value),
              lwr = quantile(value, 0.025),
              upr = quantile(value, 0.975)) %>%
    ggplot(aes(x = Var2, y = median, ymin = lwr, ymax = upr)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, color = NA) +
    facet_wrap(~paste("Region", Var1), nrow = 1) +
    geom_vline(xintercept = 65, lty = 2) +
    labs(x = 'Year', y = 'Recruitment') +
    theme_bw(base_size = 20) +
    theme(legend.position = 'none')

  reshape2::melt(sim_env_hist$Rec) %>%
    group_by(Var2, Var3) %>%
    mutate(value = sum(value)) %>%
    mutate(period = ifelse(Var2 <= 65, 'H', 'P')) %>%
    filter(period != 'H' | Var3 == 1, value != 0) %>%
    group_by(Var2, period) %>%
    summarize(median = median(value),
              lwr = quantile(value, 0.025),
              upr = quantile(value, 0.975)) %>%
    ggplot(aes(x = Var2, y = median, ymin = lwr, ymax = upr)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, color = NA) +
    geom_vline(xintercept = 65, lty = 2) +
    labs(x = 'Year', y = 'Recruitment') +
    theme_bw(base_size = 20) +
    ylim(0, NA) +
    theme(legend.position = 'none')

  reshape2::melt(sim_env_hist$TrueCatch) %>%
    mutate(period = ifelse(Var2 <= 65, 'H', 'P')) %>%
    filter(period != 'H' | Var4 == 1, value != 0) %>%
    group_by(Var1, Var2, Var3, period) %>%
    summarize(median = median(value),
              lwr = quantile(value, 0.025),
              upr = quantile(value, 0.975)) %>%
    ggplot(aes(x = Var2, y = median, ymin = lwr, ymax = upr)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, color = NA) +
    facet_grid(paste("Fleet", Var3) ~paste("Region", Var1)) +
    geom_vline(xintercept = 65, lty = 2) +
    labs(x = 'Year', y = 'Catch') +
    theme_bw(base_size = 20) +
    scale_color_manual(values = c('black', 'black')) +
    theme(legend.position = 'none')

  reshape2::melt(sim_env_hist$TrueCatch) %>%
    group_by(Var2, Var4) %>%
    mutate(value = sum(value)) %>%
    mutate(period = ifelse(Var2 <= 65, 'H', 'P')) %>%
    filter(period != 'H' | Var4 == 1, value != 0) %>%
    group_by(Var2, period) %>%
    summarize(median = median(value),
              lwr = quantile(value, 0.025),
              upr = quantile(value, 0.975)) %>%
    ggplot(aes(x = Var2, y = median,  ymin = lwr, ymax = upr)) +
    geom_line() +
    geom_ribbon(alpha = 0.5, color = NA) +
    geom_vline(xintercept = 65, lty = 2) +
    labs(x = 'Year', y = 'Catch') +
    theme_bw(base_size = 20) +
    theme(legend.position = 'none')

  reshape2::melt(sim_env_hist$Fmort) %>%
    mutate(period = ifelse(Var2 <= 65, 'H', 'P')) %>%
    filter(period != 'H' | Var4 == 1) %>%
    ggplot(aes(x = Var2, y = value, group = Var4)) +
    geom_line() +
    facet_grid(paste("Fleet", Var3) ~paste("Region", Var1)) +
    geom_vline(xintercept = 65, lty = 2) +
    labs(x = 'Year', y = 'F') +
    theme_bw(base_size = 20) +
    theme(legend.position = 'none')


  truesrv <- array(apply(sim_env_current$TrueSrvIdx, 2:4, sum), dim = c(1, dim(apply(sim_env_current$TrueSrvIdx, 2:4, sum))))
  reshape2::melt((sim_env_current$Agg_TrueSrvIdx - truesrv) / truesrv) %>%
    filter(Var3 == 1) %>%
    group_by(Var2) %>%
    summarize(median = median(value),
              lwr = quantile(value, 0.025),
              upr = quantile(value, 0.975)) %>%
    ggplot(aes(x = Var2, y = median, ymin = lwr, ymax = upr)) +
    geom_line() +
    geom_ribbon(alpha = 0.3) +
    labs(y = 'Relative Error in Aggregated Index')
