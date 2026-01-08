# Purpose: To compare the simulation results of a single-region, FAA, and three-region model
# Creator: Matthew LH. Cheng
# Date 12/24/25

conv_check <- function(sd_rep) {
  if(sd_rep$pdHess == T &  max(abs(sd_rep$gradient.fixed)) <= 0.01) TRUE
  else FALSE
}

# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)

# Read in models
sgl_rg <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp.RDS"))
sgl_rg_francis <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp_francis.RDS"))
faa <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_firsthalf.RDS")) # Model 20
faa_francis <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_finalFAA_francis.RDS"))
three_rg_low <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_lowsamp.RDS"))
three_rg_high <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_highsamp.RDS"))

# Read in OM
oms <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_lowsamp.RDS"))

# Reverse order of FAA models
faa_model <- vector("list", 100)
for(i in 1:100) faa_model[[i]] <- faa[[i]][[19]]

# Dimensions
n_yrs <- 65
n_sims <- 100
n_ages <- 30
n_sexes <- 2

# Process Results ---------------------------------------------------------
model_names <- data.frame(model = 1:6, model_name =  c("sgl", "sgl_f", "faa", "faa_f", "three_rg_low", "three_rg_high"))

# Storage containers
agg_ssb_store <- array(NA, dim = c(nrow(model_names), n_yrs, n_sims))
agg_rec_store <- array(NA, dim = c(nrow(model_names), n_yrs, n_sims))
agg_dep_store <- array(NA, dim = c(nrow(model_names), n_yrs, n_sims))
agg_ssb_om <- array(NA, dim = c(n_yrs, n_sims))
agg_rec_om <- array(NA, dim = c(n_yrs, n_sims))
agg_dep_om <- array(NA, dim = c(n_yrs, n_sims))
conv_store <- array(NA, dim = c(nrow(model_names), n_sims))

spt_ssb_store <- array(NA, dim = c(3, n_yrs, n_sims, 2))
spt_rec_store <- array(NA, dim = c(3, n_yrs, n_sims, 2))
spt_dep_store <- array(NA, dim = c(3, n_yrs, n_sims, 2))
spt_ssb_om <- array(NA, dim = c(3, n_yrs, n_sims))
spt_rec_om <- array(NA, dim = c(3, n_yrs, n_sims))
spt_dep_om <- array(NA, dim = c(3, n_yrs, n_sims))

faa_fish_sel <- array(NA, dim = c(n_ages, n_sexes, 4, 2, n_sims))
faa_srv_sel <- array(NA, dim = c(n_ages, n_sexes, 4, 2, n_sims))

for(i in 1:n_sims) {

  # Input OM values
  agg_ssb_om[,i] <- colSums(oms$SSB[,1:n_yrs,i])
  agg_rec_om[,i] <- colSums(oms$Rec[,1:n_yrs,i])
  agg_dep_om[,i] <- colSums(oms$SSB[,1:n_yrs,i]) / sum(oms$SSB[,1,i])

  # SSB
  agg_ssb_store[1,,i] <- sgl_rg[[i]]$rep$SSB
  agg_ssb_store[2,,i] <- sgl_rg_francis[[i]]$rep$SSB
  if(length(faa_model[[i]]) > 1) agg_ssb_store[3,,i] <- faa_model[[i]]$rep$SSB
  if(length(faa_francis[[i]]) > 1)  agg_ssb_store[4,,i] <- faa_francis[[i]]$rep$SSB
  if(length(three_rg_low[[i]]) > 1) agg_ssb_store[5,,i] <- colSums(three_rg_low[[i]]$rep$SSB)
  if(length(three_rg_high[[i]]) > 1) agg_ssb_store[6,,i] <- colSums(three_rg_high[[i]]$rep$SSB)

  # Recruitment
  agg_rec_store[1,,i] <- sgl_rg[[i]]$rep$Rec
  agg_rec_store[2,,i] <- sgl_rg_francis[[i]]$rep$Rec
  if(length(faa_model[[i]]) > 1) agg_rec_store[3,,i] <- faa_model[[i]]$rep$Rec
  if(length(faa_francis[[i]]) > 1)  agg_rec_store[4,,i] <- faa_francis[[i]]$rep$Rec
  if(length(three_rg_low[[i]]) > 1) agg_rec_store[5,,i] <- colSums(three_rg_low[[i]]$rep$Rec)
  if(length(three_rg_high[[i]]) > 1) agg_rec_store[6,,i] <- colSums(three_rg_high[[i]]$rep$Rec)

  # Depletion
  agg_dep_store[1,,i] <- sgl_rg[[i]]$rep$SSB / sgl_rg[[i]]$rep$SSB[1]
  agg_dep_store[2,,i] <- sgl_rg_francis[[i]]$rep$SSB / sgl_rg_francis[[i]]$rep$SSB[1]
  if(length(faa_model[[i]]) > 1) agg_dep_store[3,,i] <- faa_model[[i]]$rep$SSB / faa_model[[i]]$rep$SSB[1]
  if(length(faa_francis[[i]]) > 1) agg_dep_store[4,,i] <- faa_francis[[i]]$rep$SSB / faa_francis[[i]]$rep$SSB[1]
  if(length(three_rg_low[[i]]) > 1) agg_dep_store[5,,i] <- colSums(three_rg_low[[i]]$rep$SSB) / colSums(three_rg_low[[i]]$rep$SSB)[1]
  if(length(three_rg_high[[i]]) > 1) agg_dep_store[6,,i] <- colSums(three_rg_high[[i]]$rep$SSB) / colSums(three_rg_high[[i]]$rep$SSB)[1]

  # OM values from three region model
  spt_ssb_om[1,,i] <- oms$SSB[1,1:n_yrs,i]
  spt_ssb_om[2,,i] <- oms$SSB[2,1:n_yrs,i]
  spt_ssb_om[3,,i] <- colSums(oms$SSB[3:5,1:n_yrs,i])
  spt_rec_om[1,,i] <- oms$Rec[1,1:n_yrs,i]
  spt_rec_om[2,,i] <- oms$Rec[2,1:n_yrs,i]
  spt_rec_om[3,,i] <- colSums(oms$Rec[3:5,1:n_yrs,i])
  spt_dep_om[1,,i] <- oms$SSB[1,1:n_yrs,i] / oms$SSB[1,1,i]
  spt_dep_om[2,,i] <- oms$SSB[2,1:n_yrs,i] / oms$SSB[2,1,i]
  spt_dep_om[3,,i] <- colSums(oms$SSB[3:5,1:n_yrs,i]) / sum(oms$SSB[3:5,1,i])

  # Estimates from three region model (low)
  if(length(three_rg_low[[i]]) > 1) {
    spt_ssb_store[,,i,1] <- three_rg_low[[i]]$rep$SSB[,1:n_yrs]
    spt_rec_store[,,i,1] <- three_rg_low[[i]]$rep$Rec[,1:n_yrs]
    spt_dep_store[,,i,1] <- three_rg_low[[i]]$rep$SSB[,1:n_yrs] / three_rg_low[[i]]$rep$SSB[,1]
  }

  # Estimates from three region model (high)
  if(length(three_rg_high[[i]]) > 1) {
    spt_ssb_store[,,i,2] <- three_rg_high[[i]]$rep$SSB[,1:n_yrs]
    spt_rec_store[,,i,2] <- three_rg_high[[i]]$rep$Rec[,1:n_yrs]
    spt_dep_store[,,i,2] <- three_rg_high[[i]]$rep$SSB[,1:n_yrs] / three_rg_high[[i]]$rep$SSB[,1]
  }

  # Get FAA model selex
  if(length(faa_model[[i]]) > 1) faa_fish_sel[,,,1,i] <- faa_model[[i]]$rep$fish_sel[1,1,,,]
  if(length(faa_model[[i]]) > 1) faa_srv_sel[,,,1,i] <- faa_model[[i]]$rep$srv_sel[1,1,,,]
  # if(length(faa_francis[[i]]) > 1) faa_fish_sel[,,,2,i] <- faa_francis[[i]]$rep$fish_sel[1,1,,,]
  # if(length(faa_francis[[i]]) > 1) faa_srv_sel[,,,2,i] <- faa_francis[[i]]$rep$srv_sel[1,1,,,]

} # end i loop


# Check convergence
for(i in 1:n_sims) {
  conv_store[1,i] <- conv_check(sd_rep = sgl_rg[[i]]$sd_rep)
  conv_store[2,i] <- conv_check(sd_rep = sgl_rg_francis[[i]]$sd_rep)
  conv_store[3,i] <- conv_check(sd_rep = faa_model[[i]]$sd_rep)
  if(length(faa_francis[[i]]) > 1) conv_store[4,i] <- conv_check(sd_rep = faa_francis[[i]]$sd_rep)
  if(length(three_rg_low[[i]]) > 1)  conv_store[5,i] <- conv_check(sd_rep = three_rg_low[[i]]$sd_rep)
  if(length(three_rg_high[[i]]) > 1) conv_store[6,i] <- conv_check(sd_rep = three_rg_high[[i]]$sd_rep)
}

sqrt(diag(three_rg_low[[40]]$sd_rep$cov.fixed))
which(conv_store[5,] == F)

# figure out convergence rates
apply(conv_store, 1, sum, na.rm = T)

# OM munging
agg_ssb_om_df <- reshape2::melt(agg_ssb_om) %>%
  rename(Year = Var1, Sim = Var2, OM = value)
agg_rec_om_df <- reshape2::melt(agg_rec_om) %>%
  rename(Year = Var1, Sim = Var2, OM = value)
agg_dep_om_df <- reshape2::melt(agg_dep_om) %>%
  rename(Year = Var1, Sim = Var2, OM = value)
spt_ssb_om_df <- reshape2::melt(spt_ssb_om) %>%
  rename(Region = Var1, Year = Var2, Sim = Var3, OM = value)
spt_rec_om_df <- reshape2::melt(spt_rec_om) %>%
  rename(Region = Var1, Year = Var2, Sim = Var3, OM = value)
spt_dep_om_df <- reshape2::melt(spt_dep_om) %>%
  rename(Region = Var1, Year = Var2, Sim = Var3, OM = value)

# Munge dataframes into the right format
agg_ssb_df <- reshape2::melt(agg_ssb_store) %>%
  left_join(model_names, by = c("Var1" = 'model')) %>%
  rename(Model = Var1, Year = Var2, Sim = Var3) %>%
  left_join(agg_ssb_om_df, by = c("Year", "Sim")) %>%
  mutate(RE = (value - OM) / OM,
         Abs_RE = abs(RE))

agg_rec_df <- reshape2::melt(agg_rec_store) %>%
  left_join(model_names, by = c("Var1" = 'model')) %>%
  rename(Model = Var1, Year = Var2, Sim = Var3) %>%
  left_join(agg_rec_om_df, by = c("Year", "Sim")) %>%
  mutate(RE = (value - OM) / OM,
         Abs_RE = abs(RE))

agg_dep_df <- reshape2::melt(agg_dep_store) %>%
  left_join(model_names, by = c("Var1" = 'model')) %>%
  rename(Model = Var1, Year = Var2, Sim = Var3) %>%
  left_join(agg_dep_om_df, by = c("Year", "Sim")) %>%
  mutate(RE = (value - OM) / OM,
         Abs_RE = abs(RE))

spt_ssb_df <- reshape2::melt(spt_ssb_store) %>%
  rename(Region = Var1, Year = Var2, Sim = Var3, Data = Var4) %>%
  left_join(spt_ssb_om_df, by = c("Region", "Year", "Sim")) %>%
  mutate(RE = (value - OM) / OM,
         Abs_RE = abs(RE),
         Data = ifelse(Data == 1, 'low', 'high'))

spt_rec_df <- reshape2::melt(spt_rec_store) %>%
  rename(Region = Var1, Year = Var2, Sim = Var3, Data = Var4) %>%
  left_join(spt_rec_om_df, by = c("Region", "Year", "Sim")) %>%
  mutate(RE = (value - OM) / OM,
         Abs_RE = abs(RE),
         Data = ifelse(Data == 1, 'low', 'high'))

spt_dep_df <- reshape2::melt(spt_dep_store) %>%
  rename(Region = Var1, Year = Var2, Sim = Var3, Data = Var4) %>%
  left_join(spt_dep_om_df, by = c("Region", "Year", "Sim")) %>%
  mutate(RE = (value - OM) / OM,
         Abs_RE = abs(RE),
         Data = ifelse(Data == 1, 'low', 'high'))

# Summarize selectivities
faa_fish_sel_df <- reshape2::melt(faa_fish_sel) %>%
  rename(Age = Var1, Sex = Var2, Fleet = Var3, Francis = Var4, Sim = Var5) %>%
  mutate(
    Francis = ifelse(Francis == 1, 'No Francis', 'W Francis'),
    Sex = ifelse(Sex == 1, 'F', 'M'),
    Fleet = case_when(
      Fleet == 1 ~ 'BS_Fix',
      Fleet == 2 ~ 'AI_Fix',
      Fleet == 3 ~ 'GOA_Fix',
      Fleet == 4 ~ 'BS_Trwl',
      Fleet == 5 ~ 'AI+GOA_Trwl'
    ),
    Fleet = factor(Fleet, levels = c('BS_Fix', 'AI_Fix', 'GOA_Fix', 'BS_Trwl', 'AI+GOA_Trwl'))
  )

faa_srv_sel_df <- reshape2::melt(faa_srv_sel) %>%
  rename(Age = Var1, Sex = Var2, Fleet = Var3, Francis = Var4, Sim = Var5) %>%
  mutate(
    Francis = ifelse(Francis == 1, 'No Francis', 'W Francis'),
    Sex = ifelse(Sex == 1, 'F', 'M'),
    Fleet = case_when(
      Fleet == 1 ~ 'BS_Dom',
      Fleet == 2 ~ 'AI_Dom',
      Fleet == 3 ~ 'GOA_Dom',
      Fleet == 4 ~ 'BS+AI+GOA_JP'
    ),
    Fleet = factor(Fleet, levels = c('BS_Dom', 'AI_Dom', 'GOA_Dom', 'BS+AI+GOA_JP'))
  )

# Plots -------------------------------------------------------------------
### FAA Selectivities -------------------------------------------------------
faa_fish_sel_df %>%
  filter(Sex == 'F', Fleet == 'BS_Trwl') %>%
  drop_na() %>%
  ggplot(aes(x = Age, y = value)) +
  geom_line() +
  facet_wrap(~Sim)

faa_fish_sel_df %>%
  filter(Sex == 'M') %>%
  ggplot(aes(x = Age, y = value, group = Sim)) +
  geom_line() +
  facet_grid(Francis~Fleet)

faa_srv_sel_df %>%
  filter(Sex == 'F') %>%
  ggplot(aes(x = Age, y = value, group = Sim)) +
  geom_line() +
  facet_grid(Francis~Fleet)

faa_srv_sel_df %>%
  filter(Sex == 'M') %>%
  ggplot(aes(x = Age, y = value, group = Sim)) +
  geom_line() +
  facet_grid(Francis~Fleet)

### SSB ---------------------------------------------------------------------

# Absolute
ggplot() +
  geom_line(agg_ssb_df, mapping =  aes(x = Year, y = value, group = Sim)) +
  geom_line(agg_ssb_df %>% filter(Sim == 1), mapping =  aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_wrap(~model_name) +
  ylim(0, NA) +
  theme_bw(base_size = 15) +
  labs(y = 'SSB')

# Relative Error
ggplot() +
  geom_line(agg_ssb_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = model_name), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  # coord_cartesian(ylim = c(-0.4, 0.4)) +
  theme_bw(base_size = 15)

ggplot() +
  geom_ribbon(agg_ssb_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, ymin = lwr, ymax = upr, fill = model_name), lwd = 1.3, alpha = 0.25) +
  geom_line(agg_ssb_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = model_name), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  theme_bw(base_size = 15)


# Absolute Relative Error
ggplot() +
  geom_line(agg_ssb_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(Abs_RE, 0.025, na.rm = T),
                        upr = quantile(Abs_RE, 0.975, na.rm = T),
                        median = median(Abs_RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = model_name),  lwd = 1.3) +
  theme_bw(base_size = 15)

# Absolute
ggplot() +
  geom_line(spt_ssb_df, mapping =  aes(x = Year, y = value, group = Sim)) +
  geom_line(spt_ssb_df %>% filter(Sim == 1), mapping =  aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_grid(Region~Data, scales = 'free') +
  ylim(0, NA) +
  theme_bw(base_size = 15) +
  labs(y = 'SSB')

# Relative Error
ggplot() +
  geom_line(spt_ssb_df %>%
              group_by(Year, Region, Data) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = factor(Region)), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  theme_bw(base_size = 15) +
  facet_grid(Region~Data, scales = 'free')

# Absolute Relative Error
ggplot() +
  geom_line(spt_ssb_df %>%
              group_by(Year, Region) %>%
              summarize(lwr = quantile(Abs_RE, 0.025, na.rm = T),
                        upr = quantile(Abs_RE, 0.975, na.rm = T),
                        median = median(Abs_RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = factor(Region)),  lwd = 1.3) +
  theme_bw(base_size = 15)


### Recruitment -------------------------------------------------------------
# Absolute
ggplot() +
  geom_line(agg_rec_df, mapping =  aes(x = Year, y = value, group = Sim)) +
  geom_line(agg_rec_df %>% filter(Sim == 1), mapping =  aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_wrap(~model_name) +
  ylim(0, NA) +
  theme_bw(base_size = 15)

# Relative Error
ggplot() +
  geom_line(agg_rec_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = model_name), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  # coord_cartesian(ylim = c(-1, 1)) +
  theme_bw(base_size = 15)

# Absolute Relative Error
ggplot() +
  geom_line(agg_rec_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(Abs_RE, 0.025, na.rm = T),
                        upr = quantile(Abs_RE, 0.975, na.rm = T),
                        median = median(Abs_RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = model_name),  lwd = 1.3) +
  theme_bw(base_size = 15)

# Absolute
ggplot() +
  geom_line(spt_rec_df, mapping =  aes(x = Year, y = value, group = Sim)) +
  geom_line(spt_rec_df %>% filter(Sim == 1), mapping =  aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_grid(Region~Data, scales = 'free') +
  ylim(0, NA) +
  theme_bw(base_size = 15)

# Relative Error
ggplot() +
  geom_line(spt_rec_df %>%
              group_by(Year, Region, Data) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = factor(Region)), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  # coord_cartesian(ylim = c(-4, 4)) +
  theme_bw(base_size = 15) +
  facet_grid(Region~Data, scales = 'free')

# Absolute Relative Error
ggplot() +
  geom_line(spt_rec_df %>%
              group_by(Year, Region) %>%
              summarize(lwr = quantile(Abs_RE, 0.025, na.rm = T),
                        upr = quantile(Abs_RE, 0.975, na.rm = T),
                        median = median(Abs_RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = factor(Region)),  lwd = 1.3) +
  theme_bw(base_size = 15)


### Depletion ---------------------------------------------------------------
# Absolute
ggplot() +
  geom_line(agg_dep_df, mapping =  aes(x = Year, y = value, group = Sim)) +
  geom_line(agg_dep_df %>% filter(Sim == 1), mapping =  aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_wrap(~model_name) +
  ylim(0, NA) +
  theme_bw(base_size = 15)

# Relative Error
ggplot() +
  geom_line(agg_dep_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = model_name), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  # coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme_bw(base_size = 15)

# Absolute Relative Error
ggplot() +
  geom_line(agg_dep_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(Abs_RE, 0.025, na.rm = T),
                        upr = quantile(Abs_RE, 0.975, na.rm = T),
                        median = median(Abs_RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = model_name),  lwd = 1.3) +
  theme_bw(base_size = 15)

# Absolute
ggplot() +
  geom_line(spt_dep_df, mapping =  aes(x = Year, y = value, group = Sim)) +
  geom_line(spt_dep_df %>% filter(Sim == 1), mapping =  aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_grid(Region~Data, scales = 'free') +
  ylim(0, NA) +
  theme_bw(base_size = 15)

# Relative Error
ggplot() +
  geom_line(spt_dep_df %>%
              group_by(Year, Region, Data) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = factor(Region)), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  # coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme_bw(base_size = 15) +
  facet_grid(Region~Data, scales = 'free')

# Absolute Relative Error
ggplot() +
  geom_line(spt_dep_df %>%
              group_by(Year, Region) %>%
              summarize(lwr = quantile(Abs_RE, 0.025, na.rm = T),
                        upr = quantile(Abs_RE, 0.975, na.rm = T),
                        median = median(Abs_RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = factor(Region)),  lwd = 1.3) +
  theme_bw(base_size = 15)


# Movement ----------------------------------------------------------------
# Read in five region model
om_values <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_NoBlock_model_results.RDS"))

move_df <- data.frame()
for(i in 1:100) {
  high_move <- reshape2::melt(three_rg_high[[i]]$rep$Movement) %>% mutate(Sim = i, Type = 'high')
  if(length(three_rg_low[[i]]) > 1) low_move <- reshape2::melt(three_rg_low[[i]]$rep$Movement) %>% mutate(Sim = i, Type = 'low')
  move_df <- rbind(move_df, high_move, low_move)
}

# Melt the OM movement array
mov_long <- reshape2::melt(om_values$rep$Movement)

# 1. Collapse destinations 3:5 → 3 for from = 3
mov_3 <- mov_long %>%
  filter(from == 3) %>%
  mutate(to_collapsed = ifelse(to %in% 3:5, 3, to)) %>%
  group_by(from, to_collapsed, years, ages, sexes) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  rename(to = to_collapsed)

# 2. Keep origins 1 and 2
mov_12 <- mov_long %>%
  filter(from %in% 1:2, to %in% 1:3)   # keep to = 3 if you want to back it out

# 3. Back out destination 3 for from = 1 and 2
mov_12_backed_out <- mov_12 %>%
  filter(to != 3) %>%                  # keep only known destinations
  group_by(from, years, ages, sexes) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  mutate(to = 3,
         value = pmax(0, 1 - value))   # remainder goes to to = 3

# 4. Combine everything
mov_final <- bind_rows(
  mov_12 %>% filter(to != 3),  # from 1/2 → 1/2
  mov_12_backed_out,            # from 1/2 → 3
  mov_3                         # from 3 → 1/2/3 (collapsed)
)
# combine
move_joined_df <- move_df %>%
  left_join(mov_final %>% rename(true = value),
            by = c("from", 'to', 'years', 'ages', 'sexes'))

ggplot() +
  geom_line(move_joined_df %>% filter(years == 1, sexes == 1, Type == 'high'),
            mapping = aes(x = ages, y = value, group = interaction(Sim, Type), color = Type),
            alpha = 0.3) +
  geom_line(move_joined_df %>% filter(years == 1, sexes == 1),
            mapping = aes(x = ages, y = true), lty = 2, color = 'black', lwd = 1.3) +
  facet_grid(paste("to", to) ~paste("from", from))

# Reference Points --------------------------------------------------------

# Read in five region model to figure out global F40
om_values <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_NoBlock_model_results.RDS"))

# get true five region key quantities (since conditioned on these values)
true_five_rg_f40 <- get_key_quants(data = list(om_values$data),
                                   rep = list(om_values$rep),
                                   reference_points_opt = list(SPR_x = 0.4,
                                                               t_spwn = 0,
                                                               sex_ratio_f = 0.5,
                                                               calc_rec_st_yr = 20,
                                                               rec_age = 2,
                                                               type = "multi_region",
                                                               what = 'global_SPR'),
                                   proj_model_opt = list(n_proj_yrs = 2,
                                                         HCR_function = HCR_threshold,
                                                         recruitment_opt = 'mean_rec',
                                                         fmort_opt = 'HCR_global',
                                                         n_avg_yrs = 1),
                                   1)[[1]]

# Compute differences in F40
f40_values <- array(NA, dim = c(nrow(model_names), n_sims))
b40_values <- array(NA, dim = c(nrow(model_names), n_sims))
abc_values <- array(NA, dim = c(nrow(model_names), n_sims))

for(sim in 1:n_sims) {

  # Single region model
  sgl_key_quants <- get_key_quants(data = list(sgl_rg[[sim]]$data),
                                   rep = list(sgl_rg[[sim]]$rep),
                                   reference_points_opt = list(SPR_x = 0.4,
                                                               t_spwn = 0,
                                                               sex_ratio_f = 0.5,
                                                               calc_rec_st_yr = 20,
                                                               rec_age = 2,
                                                               type = "single_region",
                                                               what = 'SPR'),
                                   proj_model_opt = list(n_proj_yrs = 2,
                                                         HCR_function = HCR_threshold,
                                                         recruitment_opt = 'mean_rec',
                                                         fmort_opt = 'HCR',
                                                         n_avg_yrs = 1),
                                   1)[[1]]

  f40_values[1,sim] <- sgl_key_quants$F_Ref_Pt
  b40_values[1,sim] <- sgl_key_quants$B_Ref_Pt
  abc_values[1,sim] <- sgl_key_quants$Catch_Advice

  # Single region model w/ francis
  sgl_francis_key_quants <- get_key_quants(data = list(sgl_rg_francis[[sim]]$data),
                                   rep = list(sgl_rg_francis[[sim]]$rep),
                                   reference_points_opt = list(SPR_x = 0.4,
                                                               t_spwn = 0,
                                                               sex_ratio_f = 0.5,
                                                               calc_rec_st_yr = 20,
                                                               rec_age = 2,
                                                               type = "single_region",
                                                               what = 'SPR'),
                                   proj_model_opt = list(n_proj_yrs = 2,
                                                         HCR_function = HCR_threshold,
                                                         recruitment_opt = 'mean_rec',
                                                         fmort_opt = 'HCR',
                                                         n_avg_yrs = 1),
                                   1)[[1]]

  f40_values[2,sim] <- sgl_francis_key_quants$F_Ref_Pt
  b40_values[2,sim] <- sgl_francis_key_quants$B_Ref_Pt
  abc_values[2,sim] <- sgl_francis_key_quants$Catch_Advice

  # FAA model
  faa_key_quants <- get_key_quants(data = list(faa_model[[sim]]$data),
                                           rep = list(faa_model[[sim]]$rep),
                                           reference_points_opt = list(SPR_x = 0.4,
                                                                       t_spwn = 0,
                                                                       sex_ratio_f = 0.5,
                                                                       calc_rec_st_yr = 20,
                                                                       rec_age = 2,
                                                                       type = "single_region",
                                                                       what = 'SPR'),
                                           proj_model_opt = list(n_proj_yrs = 2,
                                                                 HCR_function = HCR_threshold,
                                                                 recruitment_opt = 'mean_rec',
                                                                 fmort_opt = 'HCR',
                                                                 n_avg_yrs = 1),
                                           1)[[1]]

  f40_values[3,sim] <- faa_key_quants$F_Ref_Pt
  b40_values[3,sim] <- faa_key_quants$B_Ref_Pt
  abc_values[3,sim] <- faa_key_quants$Catch_Advice


  # FAA model w/ francis
  if(length(faa_francis[[sim]]) > 1) {
    faa_francis_key_quants <- get_key_quants(data = list(faa_francis[[sim]]$data),
                                     rep = list(faa_francis[[sim]]$rep),
                                     reference_points_opt = list(SPR_x = 0.4,
                                                                 t_spwn = 0,
                                                                 sex_ratio_f = 0.5,
                                                                 calc_rec_st_yr = 20,
                                                                 rec_age = 2,
                                                                 type = "single_region",
                                                                 what = 'SPR'),
                                     proj_model_opt = list(n_proj_yrs = 2,
                                                           HCR_function = HCR_threshold,
                                                           recruitment_opt = 'mean_rec',
                                                           fmort_opt = 'HCR',
                                                           n_avg_yrs = 1),
                                     1)[[1]]

    f40_values[4,sim] <- faa_francis_key_quants$F_Ref_Pt
    b40_values[4,sim] <- faa_francis_key_quants$B_Ref_Pt
    abc_values[4,sim] <- faa_francis_key_quants$Catch_Advice
  }

  # Low three region model
  if(length(three_rg_low[[sim]]) > 1) {
    three_rg_low_key_quants <- get_key_quants(data = list(three_rg_low[[sim]]$data),
                                             rep = list(three_rg_low[[sim]]$rep),
                                             reference_points_opt = list(SPR_x = 0.4,
                                                                         t_spwn = 0,
                                                                         sex_ratio_f = 0.5,
                                                                         calc_rec_st_yr = 20,
                                                                         rec_age = 2,
                                                                         type = "multi_region",
                                                                         what = 'global_SPR'),
                                             proj_model_opt = list(n_proj_yrs = 2,
                                                                   HCR_function = HCR_threshold,
                                                                   recruitment_opt = 'mean_rec',
                                                                   fmort_opt = 'HCR_global',
                                                                   n_avg_yrs = 1),
                                             1)[[1]]

    f40_values[5,sim] <- unique(three_rg_low_key_quants$F_Ref_Pt)
    b40_values[5,sim] <- sum(three_rg_low_key_quants$B_Ref_Pt)
    abc_values[5,sim] <- sum(three_rg_low_key_quants$Catch_Advice)
  }

  # High three region model
  three_rg_high_key_quants <- get_key_quants(data = list(three_rg_high[[sim]]$data),
                                            rep = list(three_rg_high[[sim]]$rep),
                                            reference_points_opt = list(SPR_x = 0.4,
                                                                        t_spwn = 0,
                                                                        sex_ratio_f = 0.5,
                                                                        calc_rec_st_yr = 20,
                                                                        rec_age = 2,
                                                                        type = "multi_region",
                                                                        what = 'global_SPR'),
                                            proj_model_opt = list(n_proj_yrs = 2,
                                                                  HCR_function = HCR_threshold,
                                                                  recruitment_opt = 'mean_rec',
                                                                  fmort_opt = 'HCR_global',
                                                                  n_avg_yrs = 1),
                                            1)[[1]]

  f40_values[6,sim] <- unique(three_rg_high_key_quants$F_Ref_Pt)
  b40_values[6,sim] <- sum(three_rg_high_key_quants$B_Ref_Pt)
  abc_values[6,sim] <- sum(three_rg_high_key_quants$Catch_Advice)

}

# Munge dataframes
f40_df <- reshape2::melt(f40_values) %>%
  left_join(model_names, by = c("Var1" = 'model')) %>%
  rename(Model = Var1, Sim = Var2) %>%
  mutate(model_name = factor(model_name, levels = c("sgl", "sgl_f",
                                                    "faa", "faa_f",
                                                    "three_rg_low", "three_rg_high")))
b40_df <- reshape2::melt(b40_values) %>%
  left_join(model_names, by = c("Var1" = 'model')) %>%
  rename(Model = Var1, Sim = Var2) %>%
  mutate(model_name = factor(model_name, levels = c("sgl", "sgl_f",
                                                    "faa", "faa_f",
                                                    "three_rg_low", "three_rg_high")))

abc_df <- reshape2::melt(abc_values) %>%
  left_join(model_names, by = c("Var1" = 'model')) %>%
  rename(Model = Var1, Sim = Var2) %>%
  mutate(model_name = factor(model_name, levels = c("sgl", "sgl_f",
                                                    "faa", "faa_f",
                                                    "three_rg_low", "three_rg_high")))

# Absolute
ggplot(f40_df %>%
         group_by(model_name) %>%
         summarize(median = median(value, na.rm = T),
                   lwr = quantile(value, 0.025, na.rm = T),
                   upr = quantile(value, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = unique(true_five_rg_f40$F_Ref_Pt), lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'F40') +
  theme_bw()

# Relative Error
ggplot(f40_df %>%
         mutate(RE = (value - unique(true_five_rg_f40$F_Ref_Pt)) / unique(true_five_rg_f40$F_Ref_Pt)) %>%
         group_by(model_name) %>%
         summarize(median = median(RE, na.rm = T),
                   lwr = quantile(RE, 0.025, na.rm = T),
                   upr = quantile(RE, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'F40 Relative Error') +
  theme_bw()

# Absolute
ggplot(b40_df %>%
         group_by(model_name) %>%
         summarize(median = median(value, na.rm = T),
                   lwr = quantile(value, 0.025, na.rm = T),
                   upr = quantile(value, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = sum(true_five_rg_f40$B_Ref_Pt), lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'B40') +
  theme_bw()

# Relative Error
ggplot(b40_df %>%
         mutate(RE = (value - sum(true_five_rg_f40$B_Ref_Pt)) / sum(true_five_rg_f40$B_Ref_Pt)) %>%
         group_by(model_name) %>%
         summarize(median = median(RE, na.rm = T),
                   lwr = quantile(RE, 0.025, na.rm = T),
                   upr = quantile(RE, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'b40 Relative Error') +
  theme_bw()

# Absolute
ggplot(abc_df %>%
         group_by(model_name) %>%
         summarize(median = median(value, na.rm = T),
                   lwr = quantile(value, 0.025, na.rm = T),
                   upr = quantile(value, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = sum(true_five_rg_f40$Catch_Advice), lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'abc') +
  theme_bw()

# Relative Error
ggplot(abc_df %>%
         mutate(RE = (value - sum(true_five_rg_f40$Catch_Advice)) / sum(true_five_rg_f40$Catch_Advice)) %>%
         group_by(model_name) %>%
         summarize(median = median(RE, na.rm = T),
                   lwr = quantile(RE, 0.025, na.rm = T),
                   upr = quantile(RE, 0.975, na.rm = T)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = 'Model', y = 'abc Relative Error') +
  theme_bw()

# Other Diagnostic Plots -------------------------------------------------------------
# Survey Index Fits
sim <- 5
sgl_rg[[sim]]$data$ObsSrvIdx[which(sgl_rg[[sim]]$data$UseSrvIdx == 0)] <- NA
get_idx_fits_plot(list(sgl_rg[[sim]]$data), list(sgl_rg[[sim]]$rep), 1)

# Composition Fits
comp_prop <- get_comp_prop(data = sgl_rg[[sim]]$data,
                           rep = sgl_rg[[sim]]$rep,
                           age_labels = 1:30,
                           len_labels = seq(41, 99, 2),
                           year_labels = 1960:2024)

# get one step ahead fishery ages
fishages <- get_osa(obs_mat = comp_prop$Obs_FishAge_mat, # observed fishery age compositions
                    exp_mat = comp_prop$Pred_FishAge_mat, # predicted fishery age compositions
                    N = array(sgl_rg[[sim]]$data$ISS_FishAgeComps[,,1,1], c(1,65)),
                    years = list(
                      which(sgl_rg[[sim]]$data$UseFishAgeComps[1,,1] == 1)
                      # which(three_rg_high[[sim]]$data$UseFishAgeComps[2,,1] == 1),
                      # which(three_rg_high[[sim]]$data$UseFishAgeComps[3,,1] == 1)
                    ), # years with fishery ages
                    fleet = 1, # fleet
                    bins = 1:30, # age bins
                    comp_type = 2, # composition type (age-specific)
                    bin_label = "Ages" # bin labels
)

resid_plot <- SPoRC::plot_resids(osa_results = fishages)
resid_plot[[1]]
resid_plot[[2]]

# get one step ahead survey ages
srvages <- get_osa(obs_mat = comp_prop$Obs_SrvAge_mat, # observed fishery age compositions
                    exp_mat = comp_prop$Pred_SrvAge_mat, # predicted fishery age compositions
                    N = three_rg_high[[sim]]$data$ISS_SrvAgeComps[,,1,1],
                    years = list(
                      which(three_rg_high[[sim]]$data$UseSrvAgeComps[1,,1] == 1),
                      which(three_rg_high[[sim]]$data$UseSrvAgeComps[2,,1] == 1),
                      which(three_rg_high[[sim]]$data$UseSrvAgeComps[3,,1] == 1)
                    ), # years with fishery ages
                    fleet = 1, # fleet
                    bins = 1:30, # age bins
                    comp_type = 2, # composition type (age-specific)
                    bin_label = "Ages" # bin labels
)

resid_plot <- SPoRC::plot_resids(osa_results = srvages)
resid_plot[[1]]
resid_plot[[2]]
