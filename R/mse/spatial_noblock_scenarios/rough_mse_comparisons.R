# Purpose: To roughly compare MSE performance across models
# Creator: Matthew LH. Cheng
# 1/22/25


# Setup -------------------------------------------------------------------

library(here)
library(tidyverse)
library(SPoRC)

source(here("R", "functions", "mse_functions.R"))

# Read in MSEs (base recruitment)
file_dir <- here("outputs", "mse_results", "spatial_noblock_scenarios")
sgl_rg <- readRDS(here(file_dir, "single_region_base.RDS"))
faa <- readRDS(here(file_dir, "faa_base.RDS"))
five_rg <- readRDS(here(file_dir, 'five_region_base.RDS'))

# Read in and munge three region model results
files <- list.files(file_dir, full.names = TRUE)
three_rg_files <- files[stringr::str_detect(files, "three_region_base")]
all_three_rg_results <- lapply(three_rg_files, readRDS)
three_rg <- setNames(lapply(names(all_three_rg_results[[1]]), function(nm) {
  slices <- lapply(all_three_rg_results, `[[`, nm)
  if(is.array(slices[[1]])) {
    abind::abind(slices, along = length(dim(slices[[1]])))
  } else {
    do.call(c, slices)
  }
}), names(all_three_rg_results[[1]]))


# Dimensions
n_sims <- 150
st_yr <- 65
end_yr <- 95

# combine to a list
mse_list <- list(sgl_rg, faa, three_rg, five_rg)

# Storage containers
ssb_results <- array(NA, dim = c(length(mse_list), dim(sgl_rg$SSB)))
rec_results <- array(NA, dim = c(length(mse_list), dim(sgl_rg$SSB)))
catch_results <- array(NA, dim = c(length(mse_list), dim(sgl_rg$TrueCatch)))

# Loop through to get results
for(i in 1:length(mse_list)) {

  # # Remove non-converged runs / incomplete runs for three region spatial model
  if(i == 3) {
    # get valid simulations
    zero_ssb = apply((three_rg$SSB), 3, sum)
    zero_ssb = which(zero_ssb == 0)
    non_zero_ssb = (1:250)[-zero_ssb]
    null_models = sapply(three_rg$models, is.null)
    non_null_models <- which(null_models == F)
    valid_sims = non_zero_ssb[non_zero_ssb %in% non_null_models]
    mse_list[[i]]$SSB <- mse_list[[i]]$SSB[,,valid_sims]
    mse_list[[i]]$Rec <- mse_list[[i]]$Rec[,,valid_sims]
    mse_list[[i]]$TrueCatch <- mse_list[[i]]$TrueCatch[,,,valid_sims]

    # input into array
    ssb_results[i,,,1:length(valid_sims)] <- mse_list[[i]]$SSB[,,1:length(valid_sims)]
    rec_results[i,,,1:length(valid_sims)] <- mse_list[[i]]$Rec[,,1:length(valid_sims)]
    catch_results[i,,,,1:length(valid_sims)] <- mse_list[[i]]$TrueCatch[,,,1:length(valid_sims)]

  } else {
    ssb_results[i,,,] <- mse_list[[i]]$SSB[,,1:n_sims]
    rec_results[i,,,] <- mse_list[[i]]$Rec[,,1:n_sims]
    catch_results[i,,,,] <- mse_list[[i]]$TrueCatch[,,,1:n_sims]
  }
}

# Munge results
# regional ssb
ssb_rg_df <- reshape2::melt(ssb_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  mutate(region = case_when(
    region == 1 ~ 'BS',
    region == 2 ~ 'AI',
    region == 3 ~ 'WGOA',
    region == 4 ~ 'CGOA',
    region == 5 ~ 'EGOA'
  ),
  region = factor(region, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA')),
  model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# global ssb
ssb_global_df <- reshape2::melt(ssb_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  mutate(model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# regional recruitment
rec_rg_df <- reshape2::melt(rec_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  mutate(region = case_when(
    region == 1 ~ 'BS',
    region == 2 ~ 'AI',
    region == 3 ~ 'WGOA',
    region == 4 ~ 'CGOA',
    region == 5 ~ 'EGOA'
  ),
  region = factor(region, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA')),
  model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# global recruitment
rec_global_df <- reshape2::melt(rec_results) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  mutate(model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# regional catch
catch_rg_df <- reshape2::melt(catch_results) %>%
  rename(model = Var1, region = Var2, year = Var3, fleet = Var4, sim = Var5) %>%
  group_by(model, region, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  mutate(region = case_when(
    region == 1 ~ 'BS',
    region == 2 ~ 'AI',
    region == 3 ~ 'WGOA',
    region == 4 ~ 'CGOA',
    region == 5 ~ 'EGOA'
  ),
  region = factor(region, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA')),
  model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# global catch
catch_global_df <- reshape2::melt(catch_results) %>%
  rename(model = Var1, region = Var2, year = Var3, fleet = Var4, sim = Var5) %>%
  group_by(model, year, sim) %>%
  summarize(value = sum(value)) %>%
  group_by(model, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T)) %>%
  mutate(
  model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  ),
  period = ifelse(year <= st_yr, "Historical", 'Projection'))

# Terminal Year Bias
ssb_bias_results <- array(NA, dim = c(3, 95, 150))
rec_bias_results <- array(NA, dim = c(3, 95, 150))

# Loop through to get results
for(i in 1:3) { # sgl, faa, three_rg
  converged_runs <- dim(mse_list[[i]]$SSB)[3]
  for(j in 1:converged_runs) { # 150 sims
    if((sum(mse_list[[i]]$models[[j]]$rep$SSB)) != 0) ssb_bias_results[i,,j] <- (colSums(mse_list[[i]]$models[[j]]$rep$SSB) - colSums(mse_list[[i]]$SSB[,,j])) / colSums(mse_list[[i]]$SSB[,,j])
    if((sum(mse_list[[i]]$models[[j]]$rep$Rec)) != 0) rec_bias_results[i,,j] <- (colSums(mse_list[[i]]$models[[j]]$rep$Rec) - colSums(mse_list[[i]]$Rec[,,j])) / colSums(mse_list[[i]]$Rec[,,j])
   } # end j loop
} # end i loop

reshape2::melt(ssb_bias_results) %>%
  group_by(Var2, Var1) %>%
  mutate(median = median(value, na.rm = T)) %>%
  ggplot() +
  geom_line(aes(x = Var2 + 1959, y = value, group = Var3)) +
  geom_line(aes(x = Var2 + 1959, y = median, group = Var1), col = 'green', lwd = 3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 2, col = 'red') +
  facet_wrap(~Var1, labeller = labeller(Var1 = c("1" = "sgl", "2" = "faa", "3" = 'three_rg'))) +
  labs(x = 'Year', y = 'SSB Bias')

reshape2::melt(rec_bias_results) %>%
  group_by(Var2, Var1) %>%
  mutate(median = median(value, na.rm = T)) %>%
  ggplot() +
  geom_line(aes(x = Var2 + 1959, y = value, group = Var3)) +
  geom_line(aes(x = Var2 + 1959, y = median, group = Var1), col = 'green', lwd = 3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 2, col = 'red') +
  facet_wrap(~Var1, labeller = labeller(Var1 = c("1" = "sgl", "2" = "faa", "3" = 'three_rg'))) +
  labs(x = 'Year', y = 'Rec Bias')


# Compare Time Series -------------------------------------------------------------------

# Regional SSB
ggplot() +
  geom_line(ssb_rg_df %>% filter(year <= st_yr), mapping = aes(x = year, y = median), color = 'black', lwd = 1.3) +
  geom_line(ssb_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1.3) +
  # geom_ribbon(ssb_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.25, color = NA) +
  facet_wrap(~region, nrow = 1, scales = 'free') +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 65, lty = 2) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'SSB', fill = 'Model', color = 'Model') +
  theme_sablefish()

# Global SSB
ggplot() +
  geom_line(ssb_global_df %>% filter(year <= st_yr), mapping = aes(x = year, y = median), color = 'black', lwd = 1.3) +
  geom_line(ssb_global_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  # geom_ribbon(ssb_global_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)),
              # alpha = 0.25, color = NA) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 65, lty = 2) +
  theme_bw(base_size = 18) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'SSB', fill = 'Model', color = 'Model')

# Regional Catch
ggplot() +
  geom_line(catch_rg_df %>% filter(year <= st_yr), mapping = aes(x = year, y = median), color = 'black', lwd = 1.3) +
  geom_line(catch_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1.3) +
  # geom_ribbon(catch_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.25, color = NA) +
  facet_wrap(~region, nrow = 1, scales = 'free') +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 65, lty = 2) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Catch', fill = 'Model', color = 'Model') +
  theme_sablefish()

# Global Catch
ggplot() +
  geom_line(catch_global_df %>% filter(year <= st_yr), mapping = aes(x = year, y = median), color = 'black', lwd = 1.3) +
  geom_line(catch_global_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1) +
  # geom_ribbon(catch_global_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)),
              # alpha = 0.25, color = NA) +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 65, lty = 2) +
  theme_bw(base_size = 18) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Catch', fill = 'Model', color = 'Model')

# Regional Recruitment
ggplot() +
  geom_line(rec_rg_df %>% filter(year <= st_yr), mapping = aes(x = year, y = median), color = 'black', lwd = 1.3) +
  geom_line(rec_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, color = factor(model)), lwd = 1.3) +
  geom_ribbon(rec_rg_df %>% filter(year >= st_yr), mapping = aes(x = year, y = median, ymin = lwr, ymax = upr, color = factor(model), fill = factor(model)), alpha = 0.25, color = NA) +
  facet_wrap(~region, nrow = 1, scales = 'free') +
  coord_cartesian(ylim = c(0, NA)) +
  geom_vline(xintercept = 65, lty = 2) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Recruitment', fill = 'Model', color = 'Model') +
  theme_sablefish()


# Compare Apportionment ---------------------------------------------------

# Figure out apportionment for each model
get_apportionment <- function(sim_env, lls_design_type, y, sim, rolling_avg_yrs = 5) {
  # get survey apportionment (survey biomass)
  yr_lag <- rolling_avg_yrs - 1 # how many years to rollying avg across
  tmp_Srv_BiomIA <- sim_env$SrvIAA[,((y - yr_lag): y),,,1,sim] * sim_env$WAA_srv[,((y - yr_lag): y),,,1,sim]
  tmp_srvidx_biom <- apply(tmp_Srv_BiomIA, 1:2, sum) # get biomass
  apportionment <- get_single_region_survey_apportionment(feedback_start_yr = sim_env$feedback_start_yr,
                                                          n_yrs = sim_env$n_yrs,
                                                          y = y,
                                                          srv_idx = tmp_srvidx_biom, # using true survey biomass values to do apportionment
                                                          rolling_avg_yrs = rolling_avg_yrs,
                                                          lls_design_type = lls_design_type)
  return(apportionment)
}

# loop through to get apportionment
apportionment <- array(NA, dim = c(4, 5, length(st_yr:end_yr), n_sims))
for(y in st_yr:end_yr) {
  for(sim in 1:n_sims) {
    apportionment[1,,y - st_yr + 1,sim] <- rowSums(sgl_rg$TrueCatch[,y,,sim]) / sum(sgl_rg$TrueCatch[,y,,sim])
    apportionment[2,,y - st_yr + 1,sim] <- rowSums(faa$TrueCatch[,y,,sim]) / sum(faa$TrueCatch[,y,,sim])
    apportionment[3,,y - st_yr + 1,sim] <- rowSums(three_rg$TrueCatch[,y,,sim]) / sum(three_rg$TrueCatch[,y,,sim])
    apportionment[4,,y - st_yr + 1,sim] <- rowSums(five_rg$TrueCatch[,y,,sim]) / sum(five_rg$TrueCatch[,y,,sim])
  }
}

# Plot apportionment
apportionment_df <- reshape2::melt(apportionment) %>%
  rename(model = Var1, region = Var2, year = Var3, sim = Var4) %>%
  mutate(region = case_when(
    region == 1 ~ 'BS',
    region == 2 ~ 'AI',
    region == 3 ~ 'WGOA',
    region == 4 ~ 'CGOA',
    region == 5 ~ 'EGOA'
  ),
  region = factor(region, levels = c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA')),
  model = case_when(
    model == 1 ~ 'sgl',
    model == 2 ~ 'faa',
    model == 3 ~ 'three_rg',
    model == 4 ~ 'five_rg'
  )) %>%
  filter(year != 1)

# summarize apportionment
approtionment_sum <- apportionment_df %>%
  group_by(region, year, model) %>%
  group_by(model, region, year) %>%
  summarize(median = median(value, na.rm = T),
            lwr = quantile(value, 0.025, na.rm = T),
            upr = quantile(value, 0.975, na.rm = T))

ggplot(approtionment_sum, aes(x = year + 65, y = median,
                              ymin = lwr, ymax = upr,
                              color = factor(model), fill = factor(model))) +
  geom_line(lwd = 1) +
  facet_wrap(~region, nrow = 1) +
  theme_bw(base_size = 15) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(x = 'Year', y = 'Apportionment', fill = 'Model', color = 'Model')


