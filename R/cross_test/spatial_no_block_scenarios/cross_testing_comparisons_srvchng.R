# Purpose: To compare the simulation results of a single-region, FAA, and three-region model under future survey designs.
# Creator: Matthew LH. Cheng
# Date 3/24/26

# Helper functions --------------------------------------------------------
conv_check <- function(sd_rep) {
  if(sd_rep$pdHess == T &  max(abs(sd_rep$gradient.fixed)) <= 0.01) TRUE
  else FALSE
}

melt_om <- function(arr, is_spatial = FALSE) {
  df <- reshape2::melt(arr)
  if (is_spatial)
    rename(df, Region = Var1, Year = Var2, Sim = Var3, OM = value)
  else
    rename(df, Year = Var1, Sim = Var2, OM = value)
}

melt_em_srv <- function(store, om_df, is_spt = FALSE) {
  df <- reshape2::melt(store)
  if (is_spt) {
    df <- rename(df, Region = Var1, Year = Var2, Sim = Var3, Survey = Var4)
    join_keys <- c("Region", "Year", "Sim")
  } else {
    df <- df %>%
      left_join(srv_model_names, by = c("Var1" = "model")) %>%
      rename(Model = Var1, Year = Var2, Sim = Var3, Survey = Var4)
    join_keys <- c("Year", "Sim")
  }
  df %>%
    mutate(Survey = survey_designs[Survey]) %>%
    left_join(om_df, by = join_keys) %>%
    mutate(RE = (value - OM) / OM, Abs_RE = abs(RE))
}

converged <- function(mod) length(mod) > 1
get_agg  <- function(mod, qty, multi) {
  x <- mod$rep[[qty]]
  if (multi) colSums(x) else x
}

HCR_threshold <- function(x, frp, brp, alpha = 0.05) {
  stock_status <- x / brp

  if (stock_status >= 1) {
    f <- frp
  } else if (stock_status > alpha) {
    f <- frp * (stock_status - alpha) / (1 - alpha)
  } else {
    f <- 0
  }

  return(f)
}

melt_ref_pt_srv <- function(arr, true_col) {
  reshape2::melt(arr) %>%
    rename(Model = Var1, Sim = Var2, Survey = Var3) %>%
    mutate(model_name = srv_model_names$model_name[Model],
           Survey     = survey_designs[Survey]) %>%
    left_join(true_rp, by = "Sim") %>%
    mutate(RE = (value - .data[[true_col]]) / .data[[true_col]])
}

summarise_re <- function(df, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      lo  = quantile(RE, 0.025, na.rm = TRUE),
      med = median(RE,          na.rm = TRUE),
      hi  = quantile(RE, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}
# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)

# Survey Design Comparison (1960 - 202x) ------------------------------------------------
base_path <- here("outputs", "cross_test", "spatial_no_block_scenarios")

survey_designs   <- c("all", "current", "historical")
n_surveys        <- length(survey_designs)

# Dimensions
n_yrs <- 95
n_sims <- 100
n_ages <- 30
n_sexes <- 2
n_rg <- 3

srv_model_names  <- data.frame(model = 1:3, model_name = c("sgl", "faa", "three_rg"))
n_models_srv     <- nrow(srv_model_names)
multi_rg_srv     <- c(FALSE, FALSE, TRUE)

file_patterns <- c(
  sgl      = "single_region_crosstest_srvchng_%s.RDS",
  faa      = "final_faa_crosstest_srvchng_%s.RDS",
  three_rg = "three_region_crosstest_srvchng_%s.RDS"
)


# Process Results ---------------------------------------------------------
### Time Series -------------------------------------------------------------
# mods_srv[[model]][[survey]][[sim]]
mods_srv <- lapply(file_patterns, function(pat) {
  lapply(survey_designs, function(srv)
    readRDS(file.path(base_path, sprintf(pat, srv)))
  ) |> setNames(survey_designs)
}) |> setNames(names(file_patterns))

srvchng_om <- readRDS(here("outputs", "mse_results", "spatial_noblock_scenarios", "five_region_base.RDS"))

# OM storage (invariant across survey designs)
agg_ssb_om_srv <- array(NA, dim = c(n_yrs, n_sims))
agg_rec_om_srv <- array(NA, dim = c(n_yrs, n_sims))
agg_dep_om_srv <- array(NA, dim = c(n_yrs, n_sims))
spt_ssb_om_srv <- array(NA, dim = c(n_rg, n_yrs, n_sims))
spt_rec_om_srv <- array(NA, dim = c(n_rg, n_yrs, n_sims))
spt_dep_om_srv <- array(NA, dim = c(n_rg, n_yrs, n_sims))

for (i in 1:n_sims) {
  agg_ssb_om_srv[, i]   <- colSums(srvchng_om$SSB[, 1:n_yrs, i])
  agg_rec_om_srv[, i]   <- colSums(srvchng_om$Rec[, 1:n_yrs, i])
  agg_dep_om_srv[, i]   <- agg_ssb_om_srv[, i] / agg_ssb_om_srv[1, i]
  spt_ssb_om_srv[1, , i] <- srvchng_om$SSB[1, 1:n_yrs, i]
  spt_ssb_om_srv[2, , i] <- srvchng_om$SSB[2, 1:n_yrs, i]
  spt_ssb_om_srv[3, , i] <- colSums(srvchng_om$SSB[3:5, 1:n_yrs, i])
  spt_rec_om_srv[1, , i] <- srvchng_om$Rec[1, 1:n_yrs, i]
  spt_rec_om_srv[2, , i] <- srvchng_om$Rec[2, 1:n_yrs, i]
  spt_rec_om_srv[3, , i] <- colSums(srvchng_om$Rec[3:5, 1:n_yrs, i])
  for (r in 1:n_rg)
    spt_dep_om_srv[r, , i] <- spt_ssb_om_srv[r, , i] / spt_ssb_om_srv[r, 1, i]
}

# EM storage: Survey is the 4th dimension
agg_ssb_srv <- array(NA, dim = c(n_models_srv, n_yrs, n_sims, n_surveys))
agg_rec_srv <- array(NA, dim = c(n_models_srv, n_yrs, n_sims, n_surveys))
agg_dep_srv <- array(NA, dim = c(n_models_srv, n_yrs, n_sims, n_surveys))
spt_ssb_srv <- array(NA, dim = c(n_rg, n_yrs, n_sims, n_surveys))
spt_rec_srv <- array(NA, dim = c(n_rg, n_yrs, n_sims, n_surveys))
spt_dep_srv <- array(NA, dim = c(n_rg, n_yrs, n_sims, n_surveys))
conv_srv    <- array(NA, dim = c(n_models_srv, n_sims, n_surveys))

for (s in seq_along(survey_designs)) {
  srv <- survey_designs[s]
  for (i in 1:n_sims) {
    mods <- lapply(names(mods_srv), function(nm) mods_srv[[nm]][[srv]][[i]])
    # Aggregate quantities across all models
    for (m in seq_along(mods)) {
      if (!converged(mods[[m]])) next
      conv_srv[m, i, s] <- conv_check(mods[[m]]$sd_rep)
      ssb <- get_agg(mods[[m]], "SSB", multi_rg_srv[m])
      rec <- get_agg(mods[[m]], "Rec", multi_rg_srv[m])
      agg_ssb_srv[m, , i, s] <- ssb
      agg_rec_srv[m, , i, s] <- rec
      agg_dep_srv[m, , i, s] <- ssb / ssb[1]
    }
    # Spatial quantities for three_rg only
    three_mod <- mods_srv[["three_rg"]][[srv]][[i]]
    if (converged(three_mod)) {
      ssb <- three_mod$rep$SSB[, 1:n_yrs]
      spt_ssb_srv[, , i, s] <- ssb
      spt_rec_srv[, , i, s] <- three_mod$rep$Rec[, 1:n_yrs]
      spt_dep_srv[, , i, s] <- ssb / three_mod$rep$SSB[, 1]
    }
  }
}

# OM dfs
agg_ssb_om_srv_df <- melt_om(agg_ssb_om_srv)
agg_rec_om_srv_df <- melt_om(agg_rec_om_srv)
agg_dep_om_srv_df <- melt_om(agg_dep_om_srv)
spt_ssb_om_srv_df <- melt_om(spt_ssb_om_srv, is_spatial = TRUE)
spt_rec_om_srv_df <- melt_om(spt_rec_om_srv, is_spatial = TRUE)
spt_dep_om_srv_df <- melt_om(spt_dep_om_srv, is_spatial = TRUE)

# EM dfs
agg_ssb_srv_df <- melt_em_srv(agg_ssb_srv, agg_ssb_om_srv_df)
agg_rec_srv_df <- melt_em_srv(agg_rec_srv, agg_rec_om_srv_df)
agg_dep_srv_df <- melt_em_srv(agg_dep_srv, agg_dep_om_srv_df)
spt_ssb_srv_df <- melt_em_srv(spt_ssb_srv, spt_ssb_om_srv_df, is_spt = TRUE)
spt_rec_srv_df <- melt_em_srv(spt_rec_srv, spt_rec_om_srv_df, is_spt = TRUE)
spt_dep_srv_df <- melt_em_srv(spt_dep_srv, spt_dep_om_srv_df, is_spt = TRUE)

# Aggregate summaries (Year x model_name x Survey)
agg_ssb_srv_sum <- summarise_re(agg_ssb_srv_df, Year, model_name, Survey)
agg_rec_srv_sum <- summarise_re(agg_rec_srv_df, Year, model_name, Survey)
agg_dep_srv_sum <- summarise_re(agg_dep_srv_df, Year, model_name, Survey)

# Spatial summaries (Year x Region x Survey)
spt_ssb_srv_sum <- summarise_re(spt_ssb_srv_df, Year, Region, Survey) %>% mutate(
             Region = case_when(
             Region == 1 ~ "BS+AI+WGOA",
             Region == 2 ~ "CGOA",
             Region == 3 ~ "EGOA"
           )
)
spt_rec_srv_sum <- summarise_re(spt_rec_srv_df, Year, Region, Survey) %>% mutate(
  Region = case_when(
    Region == 1 ~ "BS+AI+WGOA",
    Region == 2 ~ "CGOA",
    Region == 3 ~ "EGOA"
  )
)
spt_dep_srv_sum <- summarise_re(spt_dep_srv_df, Year, Region, Survey) %>% mutate(
  Region = case_when(
    Region == 1 ~ "BS+AI+WGOA",
    Region == 2 ~ "CGOA",
    Region == 3 ~ "EGOA"
  )
)


### Reference Points and Catch Advice ---------------------------------------
ref_pt_opt_sgl   <- list(SPR_x = 0.4, t_spwn = 0, sex_ratio_f = 0.5,
                         calc_rec_st_yr = 20, rec_age = 2,
                         type = "single_region", what = "SPR")
ref_pt_opt_multi <- list(SPR_x = 0.4, t_spwn = 0, sex_ratio_f = 0.5,
                         calc_rec_st_yr = 20, rec_age = 2,
                         type = "multi_region", what = "global_SPR")
proj_opt_sgl   <- list(n_proj_yrs = 2, HCR_function = HCR_threshold,
                       recruitment_opt = "mean_rec", fmort_opt = "HCR",        n_avg_yrs = 1)
proj_opt_multi <- list(n_proj_yrs = 2, HCR_function = HCR_threshold,
                       recruitment_opt = "mean_rec", fmort_opt = "HCR_global", n_avg_yrs = 1)

# True OM reference points per sim
true_rp <- lapply(1:n_sims, function(i) {
  data_obj <- list(
    ages             = 1:srvchng_om$n_ages,
    years            = 1:n_yrs,
    n_fish_fleets    = srvchng_om$n_fish_fleets,
    n_regions        = srvchng_om$n_regions,
    n_sexes          = srvchng_om$n_sexes,
    rec_dd           = srvchng_om$rec_dd,
    rec_lag          = srvchng_om$rec_lag,
    do_recruits_move = srvchng_om$do_recruits_move,
    WAA              = array(srvchng_om$WAA[, 1:n_yrs, , , i],
                             dim = c(srvchng_om$n_regions, n_yrs, srvchng_om$n_ages, srvchng_om$n_sexes)),
    WAA_fish         = array(srvchng_om$WAA_fish[, 1:n_yrs, , , , i],
                             dim = c(srvchng_om$n_regions, n_yrs, srvchng_om$n_ages, srvchng_om$n_sexes, srvchng_om$n_fish_fleets)),
    MatAA            = array(srvchng_om$MatAA[, 1:n_yrs, , , i],
                             dim = c(srvchng_om$n_regions, n_yrs, srvchng_om$n_ages, srvchng_om$n_sexes))
  )
  rep_obj <- list(
    NAA            = array(srvchng_om$NAA[, 1:n_yrs, , , i],
                           dim = c(srvchng_om$n_regions, n_yrs, srvchng_om$n_ages, srvchng_om$n_sexes)),
    NAA0           = array(srvchng_om$NAA0[, 1:n_yrs, , , i],
                           dim = c(srvchng_om$n_regions, n_yrs, srvchng_om$n_ages, srvchng_om$n_sexes)),
    Fmort          = array(srvchng_om$Fmort[, 1:n_yrs, , i],
                           dim = c(srvchng_om$n_regions, n_yrs, srvchng_om$n_fish_fleets)),
    fish_sel       = array(srvchng_om$fish_sel[, 1:n_yrs, , , , i, drop = FALSE],
                           dim = c(srvchng_om$n_regions, n_yrs, srvchng_om$n_ages, srvchng_om$n_sexes, srvchng_om$n_fish_fleets)),
    natmort        = array(srvchng_om$natmort[, 1:n_yrs, , , i],
                           dim = c(srvchng_om$n_regions, n_yrs, srvchng_om$n_ages, srvchng_om$n_sexes)),
    sexratio       = array(srvchng_om$sexratio[, 1:n_yrs, , i],
                           dim = c(srvchng_om$n_regions, n_yrs, srvchng_om$n_sexes)),
    h_trans        = srvchng_om$h[, n_yrs, i],
    R0             = sum(srvchng_om$R0[, n_yrs, i]),
    Rec_trans_prop = srvchng_om$R0[, n_yrs, i] / sum(srvchng_om$R0[, n_yrs, i]),
    Rec            = array(srvchng_om$Rec[, 1:n_yrs, i],
                           dim = c(srvchng_om$n_regions, n_yrs)),
    SSB            = array(srvchng_om$SSB[, 1:n_yrs, i],
                           dim = c(srvchng_om$n_regions, n_yrs)),
    Movement       = array(srvchng_om$Movement[, , 1:n_yrs, , , i],
                           dim = c(srvchng_om$n_regions, srvchng_om$n_regions, n_yrs, srvchng_om$n_ages, srvchng_om$n_sexes))
  )
  kq <- get_key_quants(
    data = list(data_obj), rep = list(rep_obj),
    reference_points_opt = ref_pt_opt_multi,
    proj_model_opt       = proj_opt_multi,
    model_names = "OM")[[1]]
  data.frame(Sim = i, true_f40 = unique(kq$F_Ref_Pt),
             true_b40 = sum(kq$B_Ref_Pt), true_abc = sum(kq$Catch_Advice))
}) |> bind_rows()

# EM reference points across survey designs
f40_srv <- array(NA, dim = c(n_models_srv, n_sims, n_surveys))
b40_srv <- array(NA, dim = c(n_models_srv, n_sims, n_surveys))
abc_srv <- array(NA, dim = c(n_models_srv, n_sims, n_surveys))

for (s in seq_along(survey_designs)) {
  srv <- survey_designs[s]
  for (i in 1:n_sims) {
    mods <- lapply(names(mods_srv), function(nm) mods_srv[[nm]][[srv]][[i]])
    for (m in seq_along(mods)) {
      if (!converged(mods[[m]])) next
      multi <- multi_rg_srv[m]
      kq <- get_key_quants(
        data = list(mods[[m]]$data), rep = list(mods[[m]]$rep),
        reference_points_opt = if (multi) ref_pt_opt_multi else ref_pt_opt_sgl,
        proj_model_opt       = if (multi) proj_opt_multi   else proj_opt_sgl,
        model_names = srv_model_names$model_name[m])[[1]]
      f40_srv[m, i, s] <- if (multi) unique(kq$F_Ref_Pt)  else kq$F_Ref_Pt
      b40_srv[m, i, s] <- if (multi) sum(kq$B_Ref_Pt)     else kq$B_Ref_Pt
      abc_srv[m, i, s] <- if (multi) sum(kq$Catch_Advice) else kq$Catch_Advice
    }
  }
}

f40_srv_df <- melt_ref_pt_srv(f40_srv, "true_f40")
b40_srv_df <- melt_ref_pt_srv(b40_srv, "true_b40")
abc_srv_df <- melt_ref_pt_srv(abc_srv, "true_abc")



### Movement ----------------------------------------------------------------
true_move_vals <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_NoBlock_model_results.RDS")) # read in true movement
mov_long <- reshape2::melt(true_move_vals$rep$Movement)

# Collapse destinations 3:5 → 3 for from = 3
mov_3 <- mov_long %>%
  filter(from == 3) %>%
  mutate(to_collapsed = ifelse(to %in% 3:5, 3, to)) %>%
  group_by(from, to_collapsed, years, ages, sexes) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  rename(to = to_collapsed)

# Keep origins 1 and 2
mov_12 <- mov_long %>% filter(from %in% 1:2, to %in% 1:3)

# Back out destination 3 for from = 1 and 2
mov_12_backed_out <- mov_12 %>%
  filter(to != 3) %>%
  group_by(from, years, ages, sexes) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  mutate(to = 3, value = pmax(0, 1 - value))

# Combine everything
mov_final <- bind_rows(
  mov_12 %>% filter(to != 3),
  mov_12_backed_out,
  mov_3
)

# Extract movement estimates across survey designs
move_df <- data.frame()
for (s in seq_along(survey_designs)) {
  srv <- survey_designs[s]
  for (i in 1:n_sims) {
    mod <- mods_srv[["three_rg"]][[srv]][[i]]
    if (!converged(mod)) next
    move_df <- rbind(move_df,
                     reshape2::melt(mod$rep$Movement) %>%
                       mutate(Sim = i, Survey = srv)
    )
  }
}

# join
move_joined_df <- move_df %>%
  left_join(mov_final %>% rename(true = value),
            by = c("from", "to", "years", "ages", "sexes"))

region_labels <- c("1" = "BS+AI+WGOA", "2" = "CGOA", "3" = "EGOA")

# summarize
move_ribbon <- move_joined_df %>%
  filter(years == 1, sexes == 1) %>%
  group_by(from, to, ages, Survey) %>%
  summarise(
    lo   = quantile(value, 0.025, na.rm = TRUE),
    med  = quantile(value, 0.500, na.rm = TRUE),
    hi   = quantile(value, 0.975, na.rm = TRUE),
    true = first(true),
    .groups = "drop"
  ) %>%
  mutate(from = region_labels[as.character(from)],
         to   = region_labels[as.character(to)])


### Model Convergence -------------------------------------------------------

conv_srv_df <- reshape2::melt(conv_srv) %>%
  rename(Model = Var1, Sim = Var2, Survey = Var3) %>%
  mutate(
    model_name = srv_model_names$model_name[Model],
    Survey     = survey_designs[Survey]
  ) %>%
  group_by(model_name, Survey) %>%
  summarise(
    n_converged = sum(value, na.rm = TRUE),
    n_total     = sum(!is.na(value)),
    conv_rate   = n_converged / n_total,
    .groups     = "drop"
  )

fig_dir_ct <- here(base_path, "figs")
dir.create(fig_dir_ct, showWarnings = FALSE, recursive = TRUE)

#### Aggregate Bias ----------------------------------------------------------
ggplot(agg_ssb_srv_sum, aes(x = Year, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~model_name) +
  theme_bw(base_size = 15) +
  labs(y = "SSB Relative Error")
ggsave(file.path(fig_dir_ct, 'agg_ssb_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

ggplot(agg_rec_srv_sum, aes(x = Year, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~model_name) +
  theme_bw(base_size = 15) +
  labs(y = "Recruitment Relative Error")
ggsave(file.path(fig_dir_ct, 'agg_rec_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

ggplot(agg_dep_srv_sum, aes(x = Year, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~model_name) +
  theme_bw(base_size = 15) +
  labs(y = "Depletion Relative Error")
ggsave(file.path(fig_dir_ct, 'agg_dep_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

#### Spatial Model Bias ------------------------------------------------------
ggplot(spt_ssb_srv_sum, aes(x = Year, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~Region) +
  theme_bw(base_size = 15) +
  labs(y = "SSB Relative Error")
ggsave(file.path(fig_dir_ct, 'spt_ssb_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

ggplot(spt_rec_srv_sum, aes(x = Year, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~Region) +
  theme_bw(base_size = 15) +
  labs(y = "Recruitment Relative Error")
ggsave(file.path(fig_dir_ct, 'spt_rec_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

ggplot(spt_dep_srv_sum, aes(x = Year, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~Region) +
  theme_bw(base_size = 15) +
  labs(y = "Depletion Relative Error")
ggsave(file.path(fig_dir_ct, 'spt_dep_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

### Reference Points and Catch Advice ----------------------------------------
ggplot(f40_srv_df %>%
         group_by(model_name) %>%
         summarise(median = median(RE, na.rm = TRUE),
                   lwr    = quantile(RE, 0.025, na.rm = TRUE),
                   upr    = quantile(RE, 0.975, na.rm = TRUE)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = "Model", y = "F40 RE") +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'f40_bias_srvchng.png'), width = 6, height = 4, dpi = 300)

ggplot(b40_srv_df %>%
         group_by(model_name) %>%
         summarise(median = median(RE, na.rm = TRUE),
                   lwr    = quantile(RE, 0.025, na.rm = TRUE),
                   upr    = quantile(RE, 0.975, na.rm = TRUE)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = "Model", y = "B40 RE") +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'b40_bias_srvchng.png'), width = 6, height = 4, dpi = 300)

ggplot(abc_srv_df %>%
         group_by(model_name) %>%
         summarise(median = median(RE, na.rm = TRUE),
                   lwr    = quantile(RE, 0.025, na.rm = TRUE),
                   upr    = quantile(RE, 0.975, na.rm = TRUE)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = "Model", y = "ABC RE") +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'abc_bias_srvchng.png'), width = 6, height = 4, dpi = 300)

### Spatial Model Movement Bias ----------------------------------------------
ggplot(move_ribbon, aes(x = ages, fill = Survey, color = Survey)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3, color = NA) +
  geom_line(aes(y = med)) +
  geom_line(aes(y = true), lty = 2, color = "black", lwd = 1.3) +
  ggh4x::facet_grid2(to ~ from, scales = "free", independent = "all") +
  theme_bw(base_size = 15)
ggsave(file.path(fig_dir_ct, 'movement_bias_srvchng.png'), width = 10, height = 8, dpi = 300)

