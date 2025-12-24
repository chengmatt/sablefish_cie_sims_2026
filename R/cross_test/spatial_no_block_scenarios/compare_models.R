# Purpose: To compare the simulation results of a single-region, FAA, and three-region model
# Creator: Matthew LH. Cheng
# Date 12/24/25


# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)

# Read in models
sgl_rg <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp.RDS"))
sgl_rg_francis <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "single_region_crosstest_lowsamp_francis.RDS"))
faa <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_firsthalf.RDS")) # Model 20
faa_francis <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "faa_crosstest_lowsamp_finalFAA_francis.RDS"))
three_rg <- readRDS(here("outputs", "cross_test", "spatial_no_block_scenarios", "three_region_crosstest_lowsamp.RDS"))

# Read in OM
oms <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_lowsamp.RDS"))

# Reverse order of FAA models
faa_model <- vector("list", 100)
for(i in 1:100) faa_model[[i]] <- faa[[i]][[20]]

# Dimensions
n_yrs <- 65
n_sims <- 100

# Process Results ---------------------------------------------------------
model_names <- data.frame(model = 1:5, model_name =  c("sgl", "sgl_f", "faa", "faa_f", "three_rg"))

# Storage containers
agg_ssb_store <- array(NA, dim = c(nrow(model_names), n_yrs, n_sims))
agg_rec_store <- array(NA, dim = c(nrow(model_names), n_yrs, n_sims))
agg_dep_store <- array(NA, dim = c(nrow(model_names), n_yrs, n_sims))
agg_ssb_om <- array(NA, dim = c(n_yrs, n_sims))
agg_rec_om <- array(NA, dim = c(n_yrs, n_sims))
agg_dep_om <- array(NA, dim = c(n_yrs, n_sims))

spt_ssb_store <- array(NA, dim = c(3, n_yrs, n_sims))
spt_rec_store <- array(NA, dim = c(3, n_yrs, n_sims))
spt_dep_store <- array(NA, dim = c(3, n_yrs, n_sims))
spt_ssb_om <- array(NA, dim = c(3, n_yrs, n_sims))
spt_rec_om <- array(NA, dim = c(3, n_yrs, n_sims))
spt_dep_om <- array(NA, dim = c(3, n_yrs, n_sims))

for(i in 1:n_sims) {

  # Input OM values
  agg_ssb_om[,i] <- colSums(oms$SSB[,1:n_yrs,i])
  agg_rec_om[,i] <- colSums(oms$Rec[,1:n_yrs,i])
  agg_dep_om[,i] <- colSums(oms$SSB[,1:n_yrs,i]) / sum(oms$SSB[,1,i])

  # SSB
  agg_ssb_store[1,,i] <- sgl_rg[[i]]$rep$SSB
  agg_ssb_store[2,,i] <- sgl_rg_francis[[i]]$rep$SSB
  if(length(faa_model[[i]]) > 1) agg_ssb_store[3,,i] <- faa_model[[i]]$rep$SSB
  if(!is.null(faa_francis[[i]])) agg_ssb_store[4,,i] <- faa_francis[[i]]$rep$SSB
  if(length(three_rg[[i]]) > 1) agg_ssb_store[5,,i] <- colSums(three_rg[[i]]$rep$SSB)

  # Recruitment
  agg_rec_store[1,,i] <- sgl_rg[[i]]$rep$Rec
  agg_rec_store[2,,i] <- sgl_rg_francis[[i]]$rep$Rec
  if(length(faa_model[[i]]) > 1) agg_rec_store[3,,i] <- faa_model[[i]]$rep$Rec
  if(!is.null(faa_francis[[i]])) agg_rec_store[4,,i] <- faa_francis[[i]]$rep$Rec
  if(length(three_rg[[i]]) > 1) agg_rec_store[5,,i] <- colSums(three_rg[[i]]$rep$Rec)

  # Depletion
  agg_dep_store[1,,i] <- sgl_rg[[i]]$rep$SSB / sgl_rg[[i]]$rep$SSB[1]
  agg_dep_store[2,,i] <- sgl_rg_francis[[i]]$rep$SSB / sgl_rg_francis[[i]]$rep$SSB[1]
  if(length(faa_model[[i]]) > 1) agg_dep_store[3,,i] <- faa_model[[i]]$rep$SSB / faa_model[[i]]$rep$SSB[1]
  if(!is.null(faa_francis[[i]])) agg_dep_store[4,,i] <- faa_francis[[i]]$rep$SSB / faa_francis[[i]]$rep$SSB[1]
  if(length(three_rg[[i]]) > 1) agg_dep_store[5,,i] <- colSums(three_rg[[i]]$rep$SSB) / colSums(three_rg[[i]]$rep$SSB)[1]

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

  # Estimates from three region model
  if(length(three_rg[[i]]) > 1) {
    spt_ssb_store[,,i] <- three_rg[[i]]$rep$SSB[,1:n_yrs]
    spt_rec_store[,,i] <- three_rg[[i]]$rep$Rec[,1:n_yrs]
    spt_dep_store[,,i] <- three_rg[[i]]$rep$SSB[,1:n_yrs] / three_rg[[i]]$rep$SSB[,1]
  }

} # end i loop

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
  rename(Region = Var1, Year = Var2, Sim = Var3) %>%
  left_join(spt_ssb_om_df, by = c("Region", "Year", "Sim")) %>%
  mutate(RE = (value - OM) / OM,
         Abs_RE = abs(RE))

spt_rec_df <- reshape2::melt(spt_rec_store) %>%
  rename(Region = Var1, Year = Var2, Sim = Var3) %>%
  left_join(spt_rec_om_df, by = c("Region", "Year", "Sim")) %>%
  mutate(RE = (value - OM) / OM,
         Abs_RE = abs(RE))

spt_dep_df <- reshape2::melt(spt_dep_store) %>%
  rename(Region = Var1, Year = Var2, Sim = Var3) %>%
  left_join(spt_dep_om_df, by = c("Region", "Year", "Sim")) %>%
  mutate(RE = (value - OM) / OM,
         Abs_RE = abs(RE))

# Plots -------------------------------------------------------------------
### SSB ---------------------------------------------------------------------

# Absolute
ggplot() +
  geom_line(agg_ssb_df, mapping =  aes(x = Year, y = value, group = Sim)) +
  geom_line(agg_ssb_df %>% filter(Sim == 1), mapping =  aes(x = Year, y = OM), lty = 2, color = 'red', lwd = 1.3) +
  facet_wrap(~model_name) +
  ylim(0, NA) +
  theme_bw(base_size = 15)

# Relative Error
ggplot() +
  geom_line(agg_ssb_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = model_name), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  theme_bw(base_size = 15)

ggplot() +
  geom_ribbon(agg_ssb_df %>%
              group_by(Year, model_name) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, ymin = lwr, ymax = upr, fill = model_name), lwd = 1.3, alpha = 0.25) +
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
  facet_wrap(~Region, scales = 'free') +
  ylim(0, NA) +
  theme_bw(base_size = 15)

# Relative Error
ggplot() +
  geom_line(spt_ssb_df %>%
              group_by(Year, Region) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = factor(Region)), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  theme_bw(base_size = 15)

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
  coord_cartesian(ylim = c(-3, 3)) +
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
  facet_wrap(~Region, scales = 'free') +
  ylim(0, NA) +
  theme_bw(base_size = 15)

# Relative Error
ggplot() +
  geom_line(spt_rec_df %>%
              group_by(Year, Region) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = factor(Region)), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  coord_cartesian(ylim = c(-4, 4)) +
  theme_bw(base_size = 15)

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
  coord_cartesian(ylim = c(-0.5, 0.5)) +
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
  facet_wrap(~Region, scales = 'free') +
  ylim(0, NA) +
  theme_bw(base_size = 15)

# Relative Error
ggplot() +
  geom_line(spt_dep_df %>%
              group_by(Year, Region) %>%
              summarize(lwr = quantile(RE, 0.025, na.rm = T),
                        upr = quantile(RE, 0.975, na.rm = T),
                        median = median(RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = factor(Region)), lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3, col = 'black') +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme_bw(base_size = 15)

# Absolute Relative Error
ggplot() +
  geom_line(spt_dep_df %>%
              group_by(Year, Region) %>%
              summarize(lwr = quantile(Abs_RE, 0.025, na.rm = T),
                        upr = quantile(Abs_RE, 0.975, na.rm = T),
                        median = median(Abs_RE, na.rm = T)),
            mapping =  aes(x = Year, y = median, color = factor(Region)),  lwd = 1.3) +
  theme_bw(base_size = 15)


