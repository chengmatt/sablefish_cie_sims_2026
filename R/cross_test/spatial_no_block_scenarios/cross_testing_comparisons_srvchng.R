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

melt_om_naa <- function(arr, is_spatial = FALSE) {
  df <- reshape2::melt(arr)
  if (is_spatial)
    rename(df, Region = Var1, Year = Var2, Age = Var3, Sex = Var4, Sim = Var5, OM = value)
  else
    rename(df, Year = Var1, Age = Var2, Sex = Var3, Sim = Var4, OM = value)
}

melt_em_srv_naa <- function(store, om_df, is_spt = FALSE) {
  df <- reshape2::melt(store)
  if (is_spt) {
    df <- rename(df, Region = Var1, Year = Var2, Age = Var3, Sex = Var4, Sim = Var5, Survey = Var6)
    join_keys <- c("Region", "Year", "Age", "Sex", "Sim")
  } else {
    df <- df %>%
      left_join(srv_model_names, by = c("Var1" = "model")) %>%
      rename(Model = Var1, Year = Var2, Age = Var3, Sex = Var4, Sim = Var5, Survey = Var6)
    join_keys <- c("Year", "Age", "Sex", "Sim")
  }
  df %>%
    mutate(Survey = survey_designs[Survey]) %>%
    left_join(om_df, by = join_keys) %>%
    mutate(RE = (value - OM) / OM, Abs_RE = abs(RE))
}

get_agg_naa <- function(mod, multi) {
  naa <- mod$rep$NAA
  apply(naa[, 1:n_yrs, , ,drop = F], c(2, 3, 4), sum)
}


get_om_caa_sgl <- function(i) {
  caa <- srvchng_om$CAA[, 1:n_yrs, , , , i]
  out <- array(NA, dim = c(n_yrs, n_ages, n_sexes, 2))
  out[, , , 1] <- apply(caa[, , , , 1], c(2, 3, 4), sum)  # all regions, fixed
  out[, , , 2] <- apply(caa[, , , , 2], c(2, 3, 4), sum)  # all regions, trawl
  out
}

get_om_caa_faa <- function(i) {
  caa <- srvchng_om$CAA[, 1:n_yrs, , , , i]
  out <- array(NA, dim = c(n_yrs, n_ages, n_sexes, 4))
  out[, , , 1] <- caa[1,  , , , 1]  # BS fixed
  out[, , , 2] <- caa[1,  , , , 2]  # AI fixed
  out[, , , 3] <- apply(caa[3:5, , , , 1, drop = FALSE],               c(2, 3, 4), sum)  # GOA fixed
  out[, , , 4] <- apply(caa[,   , , , 2, drop = FALSE],                c(2, 3, 4), sum)  # all trawl
  out
}

get_om_caa_threerg <- function(i) {
  caa <- srvchng_om$CAA[, 1:n_yrs, , , , i]
  out <- array(NA, dim = c(n_rg, n_yrs, n_ages, n_sexes, 2))
  out[1, , , , 1] <- caa[1,  , , , 1]
  out[2, , , , 1] <- caa[2,  , , , 1]
  out[3, , , , 1] <- apply(caa[3:5, , , , 1], c(2, 3, 4), sum)
  out[1, , , , 2] <- caa[1,  , , , 2]
  out[2, , , , 2] <- caa[2,  , , , 2]
  out[3, , , , 2] <- apply(caa[3:5, , , , 2], c(2, 3, 4), sum)
  out
}

get_em_caa <- function(mod, m) {
  caa <- mod$rep$CAA  # [region, yr, age, sex, fleet]
  n_f <- n_fleets_by_model[m]
  if (m == 3) {
    caa[, 1:n_yrs, , , ]
  } else {
    # sgl/faa: sum over region dim (should be 1 anyway), return [yr, age, sex, fleet]
    apply(caa[, 1:n_yrs, , , , drop = FALSE], c(2, 3, 4, 5), sum)
  }
}

make_caa_store <- function(n_f, spatial = FALSE) {
  if (spatial)
    array(NA, dim = c(n_rg, n_yrs, n_ages, n_sexes, n_f, n_sims, n_surveys))
  else
    array(NA, dim = c(n_yrs, n_ages, n_sexes, n_f, n_sims, n_surveys))
}

make_caa_re <- function(em_arr, om_arr, model_label, spatial = FALSE) {
  em_df <- reshape2::melt(em_arr)
  om_df <- reshape2::melt(om_arr)
  if (spatial) {
    em_df <- rename(em_df, Region = Var1, Year = Var2, Age = Var3, Sex = Var4,
                    Fleet = Var5, Sim = Var6, Survey = Var7)
    om_df <- rename(om_df, Region = Var1, Year = Var2, Age = Var3, Sex = Var4,
                    Fleet = Var5, Sim = Var6, OM = value)
    join_keys <- c("Region", "Year", "Age", "Sex", "Fleet", "Sim")
  } else {
    em_df <- rename(em_df, Year = Var1, Age = Var2, Sex = Var3,
                    Fleet = Var4, Sim = Var5, Survey = Var6)
    om_df <- rename(om_df, Year = Var1, Age = Var2, Sex = Var3,
                    Fleet = Var4, Sim = Var5, OM = value)
    join_keys <- c("Year", "Age", "Sex", "Fleet", "Sim")
  }
  em_df %>%
    mutate(Survey = survey_designs[Survey]) %>%
    left_join(om_df, by = join_keys) %>%
    mutate(RE = (value - OM) / OM, model_name = model_label)
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
             Region == 1 ~ "BS",
             Region == 2 ~ "AI",
             Region == 3 ~ "GOA"
           )
)
spt_rec_srv_sum <- summarise_re(spt_rec_srv_df, Year, Region, Survey) %>% mutate(
  Region = case_when(
    Region == 1 ~ "BS",
    Region == 2 ~ "AI",
    Region == 3 ~ "GOA"
  )
)
spt_dep_srv_sum <- summarise_re(spt_dep_srv_df, Year, Region, Survey) %>%
  mutate(
  Region = case_when(
    Region == 1 ~ "BS",
    Region == 2 ~ "AI",
    Region == 3 ~ "GOA"
  )
)


### Numbers at Age ----------------------------------------------------------
# OM NAA storage
agg_naa_om_srv <- array(NA, dim = c(n_yrs, n_ages, n_sexes, n_sims))
spt_naa_om_srv <- array(NA, dim = c(n_rg, n_yrs, n_ages, n_sexes, n_sims))

for (i in 1:n_sims) {
  agg_naa_om_srv[, , , i]    <- apply(srvchng_om$NAA[, 1:n_yrs, , , i], c(2, 3, 4), sum)
  spt_naa_om_srv[1, , , , i] <- srvchng_om$NAA[1, 1:n_yrs, , , i]
  spt_naa_om_srv[2, , , , i] <- srvchng_om$NAA[2, 1:n_yrs, , , i]
  spt_naa_om_srv[3, , , , i] <- apply(srvchng_om$NAA[3:5, 1:n_yrs, , , i], c(2, 3, 4), sum)
}

# EM NAA storage
agg_naa_srv <- array(NA, dim = c(n_models_srv, n_yrs, n_ages, n_sexes, n_sims, n_surveys))
spt_naa_srv <- array(NA, dim = c(n_rg,         n_yrs, n_ages, n_sexes, n_sims, n_surveys))

for (s in seq_along(survey_designs)) {
  srv <- survey_designs[s]
  for (i in 1:n_sims) {
    mods <- lapply(names(mods_srv), function(nm) mods_srv[[nm]][[srv]][[i]])
    for (m in seq_along(mods)) {
      if (!converged(mods[[m]]))           next
      if (!isTRUE(conv_srv[m, i, s]))      next
      agg_naa_srv[m, , , , i, s] <- get_agg_naa(mods[[m]], multi_rg_srv[m])
    }
    three_mod <- mods_srv[["three_rg"]][[srv]][[i]]
    if (converged(three_mod) && isTRUE(conv_srv[3, i, s]))
      spt_naa_srv[, , , , i, s] <- three_mod$rep$NAA[, 1:n_yrs, , ]
  }
}

# Melt to long format
agg_naa_om_srv_df <- melt_om_naa(agg_naa_om_srv)
spt_naa_om_srv_df <- melt_om_naa(spt_naa_om_srv, is_spatial = TRUE)

agg_naa_srv_df <- melt_em_srv_naa(agg_naa_srv, agg_naa_om_srv_df)
spt_naa_srv_df <- melt_em_srv_naa(spt_naa_srv, spt_naa_om_srv_df, is_spt = TRUE)

# Summaries
agg_naa_srv_sum <- agg_naa_srv_df %>%
  group_by(model_name, Survey, Year, Age) %>%
  summarise(med = median(RE, na.rm = TRUE), .groups = "drop") %>%
  mutate(Year = Year + 1959, Age = Age + 1,Cohort = Year - Age)

agg_naa_srv_sum_abs <- agg_naa_srv_df %>%
  group_by(model_name, Survey, Age) %>%
  summarise(OM = median(OM, na.rm = TRUE),
            value = median(value, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(Age = Age + 1)

spt_naa_srv_sum <- spt_naa_srv_df %>%
  mutate(Region = case_when(
    Region == 1 ~ "BS",
    Region == 2 ~ "AI",
    Region == 3 ~ "GOA"
  )) %>%
  group_by(Region, Survey, Year, Age) %>%
  summarise(med = median(RE, na.rm = TRUE), .groups = "drop") %>%
  mutate(Year = Year + 1959, Age = Age + 1,Cohort = Year - Age)


### Catch at Age ------------------------------------------------------------
n_fleets_by_model <- c(2, 4, 2)  # sgl, faa, three_rg

# storage
caa_om_sgl   <- array(NA, dim = c(n_yrs, n_ages, n_sexes, 2, n_sims))
caa_om_faa   <- array(NA, dim = c(n_yrs, n_ages, n_sexes, 4, n_sims))
caa_om_3rg   <- array(NA, dim = c(n_rg, n_yrs, n_ages, n_sexes, 2, n_sims))
caa_em_sgl   <- make_caa_store(2)
caa_em_faa   <- make_caa_store(4)
caa_em_3rg   <- make_caa_store(2, spatial = TRUE)

# get om values
for (i in 1:n_sims) {
  caa_om_sgl[, , , , i] <- get_om_caa_sgl(i)
  caa_om_faa[, , , , i] <- get_om_caa_faa(i)
  caa_om_3rg[, , , , , i] <- get_om_caa_threerg(i)
}

# get em values
for (s in seq_along(survey_designs)) {
  srv <- survey_designs[s]
  for (i in 1:n_sims) {
    # SGL (m=1)
    mod <- mods_srv[["sgl"]][[srv]][[i]]
    if (converged(mod) && isTRUE(conv_srv[1, i, s]))
      caa_em_sgl[, , , , i, s] <- get_em_caa(mod, 1)
    # FAA (m=2)
    mod <- mods_srv[["faa"]][[srv]][[i]]
    if (converged(mod) && isTRUE(conv_srv[2, i, s]))
      caa_em_faa[, , , , i, s] <- get_em_caa(mod, 2)
    # three_rg (m=3)
    mod <- mods_srv[["three_rg"]][[srv]][[i]]
    if (converged(mod) && isTRUE(conv_srv[3, i, s]))
      caa_em_3rg[, , , , , i, s] <- get_em_caa(mod, 3)
  }
}

# CAA for single region model
caa_re_sgl <- make_caa_re(caa_em_sgl, caa_om_sgl, "sgl") %>%
  mutate(Year = Year + 1959)

# aggregated CAA for single region model
agg_caa_sgl_sum <- caa_re_sgl %>%
  group_by(model_name, Survey, Year, Age) %>%
  summarize(OM = mean(OM),
            value = mean(value))

# CAA for faa model
caa_re_faa <- make_caa_re(caa_em_faa, caa_om_faa, "faa") %>%
  mutate(Year = Year + 1959)

# aggregated CAA for single region model
agg_caa_faa_sum <- caa_re_faa %>%
  group_by(model_name, Survey, Year, Age) %>%
  summarize(OM = mean(OM),
            value = mean(value))

# CAA for faa model
caa_re_faa <- make_caa_re(caa_em_faa, caa_om_faa, "faa") %>%
  mutate(Year = Year + 1959)

# aggregated CAA for faa model
agg_caa_faa_sum <- caa_re_faa %>%
  group_by(model_name, Survey, Year, Age) %>%
  summarize(OM = mean(OM, na.rm = T),
            value = mean(value, na.rm = T))

# CAA for faa model by fleet
agg_caa_faa_fleet_sum <- caa_re_faa %>%
  group_by(model_name, Survey, Age, Fleet) %>%
  summarize(OM = mean(OM, na.rm = T),
            value = mean(value, na.rm = T))

# CAA for three region model
caa_re_3rg <- make_caa_re(caa_em_3rg, caa_om_3rg, "three_rg", spatial = TRUE) %>%
  mutate(Year = Year + 1959)

# aggregated CAA for 3_rg model
agg_caa_3rg_sum <- caa_re_3rg %>%
  group_by(model_name, Survey, Year, Age) %>%
  summarize(OM = mean(OM, na.rm = T),
            value = mean(value, na.rm = T))

# CAA for 3rg model by fleet
agg_caa_3rg_rg_sum <- caa_re_3rg %>%
  group_by(model_name, Survey, Age, Region) %>%
  summarize(OM = mean(OM, na.rm = T),
            value = mean(value, na.rm = T))


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

region_labels <- c("1" = "BS", "2" = "AI", "3" = "GOA")

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

fig_dir_ct <- here('outputs', "cross_test", "figs")
dir.create(fig_dir_ct, showWarnings = FALSE, recursive = TRUE)

### Aggregate Bias ----------------------------------------------------------
ggplot(agg_ssb_srv_sum, aes(x = Year + 1959, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~model_name) +
  theme_bw(base_size = 15) +
  labs(y = "SSB Relative Error", x = 'Year')
ggsave(file.path(fig_dir_ct, 'agg_ssb_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

ggplot(agg_rec_srv_sum, aes(x = Year + 1959, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~model_name) +
  theme_bw(base_size = 15) +
  labs(y = "Recruitment Relative Error", x = 'Year')
ggsave(file.path(fig_dir_ct, 'agg_rec_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

ggplot(agg_rec_srv_sum %>% filter(Year >= 50, Survey == 'current'),
       aes(x = Year + 1959, y = med)) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  geom_vline(xintercept = c(2018, 2025), lty = 2, color = 'black', lwd = 1 ) +
  facet_wrap(~model_name) +
  theme_bw(base_size = 15) +
  labs(y = "Recruitment Relative Error", x = 'Year')
ggsave(file.path(fig_dir_ct, 'agg_rec_bias_current_zoom.png'), width = 10, height = 4, dpi = 300)

ggplot(agg_dep_srv_sum, aes(x = Year + 1959, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~model_name) +
  theme_bw(base_size = 15) +
  labs(y = "Depletion Relative Error", x = 'Year')
ggsave(file.path(fig_dir_ct, 'agg_dep_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

### Spatial Model Bias ------------------------------------------------------
ggplot(spt_ssb_srv_sum, aes(x = Year + 1959, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~Region) +
  geom_vline(xintercept = c(2018, 2025), lty = 2, color = 'black', lwd = 1 ) +
  theme_bw(base_size = 15) +
  labs(y = "SSB Relative Error", x = 'Year')
ggsave(file.path(fig_dir_ct, 'spt_ssb_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

ggplot(spt_rec_srv_sum, aes(x = Year + 1959, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~Region) +
  theme_bw(base_size = 15) +
  labs(y = "Recruitment Relative Error", x = 'Year')
ggsave(file.path(fig_dir_ct, 'spt_rec_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

ggplot(spt_dep_srv_sum, aes(x = Year + 1959, y = med, ymin = lo, ymax = hi,
                            color = Survey, fill = Survey)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line(lwd = 1.3) +
  geom_vline(xintercept = c(2018, 2025), lty = 2, color = 'black', lwd = 1 ) +
  geom_hline(yintercept = 0, lty = 2, color = "black", lwd = 1) +
  facet_wrap(~Region) +
  theme_bw(base_size = 15) +
  labs(y = "Depletion Relative Error", x = 'Year')
ggsave(file.path(fig_dir_ct, 'spt_dep_bias_srvchng.png'), width = 10, height = 4, dpi = 300)

### Reference Points and Catch Advice ----------------------------------------
ggplot(f40_srv_df %>%
         group_by(model_name, Survey) %>%
         summarise(median = median(RE, na.rm = TRUE),
                   lwr    = quantile(RE, 0.025, na.rm = TRUE),
                   upr    = quantile(RE, 0.975, na.rm = TRUE)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr, color = Survey)) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = "Model", y = "F40 RE") +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'f40_bias_srvchng.png'), width = 6, height = 4, dpi = 300)

ggplot(b40_srv_df %>%
         group_by(model_name, Survey) %>%
         summarise(median = median(RE, na.rm = TRUE),
                   lwr    = quantile(RE, 0.025, na.rm = TRUE),
                   upr    = quantile(RE, 0.975, na.rm = TRUE)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr, color = Survey)) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  labs(x = "Model", y = "B40 RE") +
  theme_bw()
ggsave(file.path(fig_dir_ct, 'b40_bias_srvchng.png'), width = 6, height = 4, dpi = 300)

ggplot(abc_srv_df %>%
         group_by(model_name, Survey) %>%
         summarise(median = median(RE, na.rm = TRUE),
                   lwr    = quantile(RE, 0.025, na.rm = TRUE),
                   upr    = quantile(RE, 0.975, na.rm = TRUE)),
       aes(x = model_name, y = median, ymin = lwr, ymax = upr, color = Survey)) +
  geom_pointrange(position = position_dodge(width = 0.3)) +
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


### Numbers at Age ----------------------------------------------------------
# Average NAA
ggplot() +
  geom_line(agg_naa_srv_sum_abs,
            mapping = aes(x = Age, y = OM), color = 'black') +
  geom_line(agg_naa_srv_sum_abs,
            mapping = aes(x = Age, y = value, color = model_name)) +
  theme_bw(base_size = 14) +
  facet_wrap(~Survey) +
  labs(y = 'NAA')

# Bias in NAA across time
ggplot(agg_naa_srv_sum %>% filter(Survey == 'all'), aes(x = Year, y = med, color = model_name)) +
  geom_line() +
  facet_wrap(~Age, scales = 'free') +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  theme_bw(base_size = 14)

# Spatial bias
ggplot(spt_naa_srv_sum %>% filter(Survey == 'all'), aes(x = Year, y = med, color = Region)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2, lwd = 1.3) +
  facet_wrap(~Age, scales = 'free') +
  theme_bw(base_size = 14)


### Catch at Age ------------------------------------------------------------

# SGL CAA
ggplot() +
  geom_line(agg_caa_sgl_sum, mapping = aes(x = Age, y = OM)) +
  geom_line(agg_caa_sgl_sum, mapping = aes(x = Age, y = value, color = Survey)) +
  facet_wrap(~Year, scales = 'free') +
  labs(y = 'CAA')

# SGL CAA aggregated
ggplot() +
  geom_line(agg_caa_sgl_sum %>% group_by(Age, Survey) %>% summarize(OM = mean(OM)), mapping = aes(x = Age, y = OM)) +
  geom_line(agg_caa_sgl_sum %>% group_by(Age, Survey) %>% summarize(value = mean(value)),
            mapping = aes(x = Age, y = value, color = Survey)) +
  labs(y = 'CAA')

# faa CAA
ggplot() +
  geom_line(agg_caa_faa_sum, mapping = aes(x = Age, y = OM)) +
  geom_line(agg_caa_faa_sum, mapping = aes(x = Age, y = value, color = Survey)) +
  facet_wrap(~Year, scales = 'free') +
  labs(y = 'CAA')

# faa CAA aggregated
ggplot() +
  geom_line(agg_caa_faa_sum %>% group_by(Age, Survey) %>% summarize(OM = mean(OM)), mapping = aes(x = Age, y = OM)) +
  geom_line(agg_caa_faa_sum %>% group_by(Age, Survey) %>% summarize(value = mean(value)),
            mapping = aes(x = Age, y = value, color = Survey)) +
  labs(y = 'CAA')

# faa CAA by fleet
ggplot() +
  geom_line(agg_caa_faa_fleet_sum, mapping = aes(x = Age, y = OM)) +
  geom_line(agg_caa_faa_fleet_sum, mapping = aes(x = Age, y = value, color = Survey)) +
  labs(y = 'CAA') +
  facet_wrap(~Fleet, scales = 'free')

# three_region CAA
ggplot() +
  geom_line(agg_caa_3rg_sum, mapping = aes(x = Age, y = OM)) +
  geom_line(agg_caa_3rg_sum, mapping = aes(x = Age, y = value, color = Survey)) +
  facet_wrap(~Year, scales = 'free') +
  labs(y = 'CAA')

# three_region aggregated
ggplot() +
  geom_line(agg_caa_3rg_sum %>% group_by(Age, Survey) %>% summarize(OM = mean(OM)), mapping = aes(x = Age, y = OM)) +
  geom_line(agg_caa_3rg_sum %>% group_by(Age, Survey) %>% summarize(value = mean(value)),
            mapping = aes(x = Age, y = value, color = Survey)) +
  labs(y = 'CAA')

# three_region by region
ggplot() +
  geom_line(agg_caa_3rg_rg_sum, mapping = aes(x = Age, y = OM)) +
  geom_line(agg_caa_3rg_rg_sum, mapping = aes(x = Age, y = value, color = Survey)) +
  labs(y = 'CAA') +
  facet_wrap(~Region, scales = 'free')



### Survey Stuff ------------------------------------------------------------
# get mean observed survey indices
mean_ObsSrvIdx_sgl <- lapply(survey_designs, function(srv) {
  vals <- lapply(seq_len(n_sims), function(i) {
    mod <- mods_srv[["sgl"]][[srv]][[i]]
    if (converged(mod)) mod$data$ObsSrvIdx else NULL
  })
  vals <- Filter(Negate(is.null), vals)
  Reduce("+", vals) / length(vals)
})
names(mean_ObsSrvIdx_sgl) <- survey_designs

# get true index
true_agg_srv <- apply(apply(srvchng_om$TrueSrvIdx, 2:4, sum), 1:2, mean)
true_df <- reshape2::melt(true_agg_srv, varnames = c("Year", "Fleet"), value.name = "value") %>%
  mutate(type = "True")

obs_df <- reshape2::melt(mean_ObsSrvIdx_sgl) %>%
  rename(Region = Var1, Year = Var2, Fleet = Var3, Survey = L1) %>%
  group_by(Survey, Year, Fleet) %>%
  summarize(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(type = "Observed", Year = Year + 1959)

true_df <- reshape2::melt(true_agg_srv, varnames = c("Year", "Fleet"), value.name = "value") %>%
  mutate(type = "True", Survey = "True", Year = Year + 1959)

# plot absolute differences
plot_df <- bind_rows(obs_df, true_df)
ggplot(plot_df %>% filter(Fleet == 1), aes(x = Year, y = value, color = Survey, linetype = type)) +
  geom_line(lwd = 1) +
  theme_bw()

# Plot relative error
re_df <- obs_df %>%
  left_join(true_df %>% select(Year, Fleet, value) %>% rename(true_val = value),
            by = c("Year", "Fleet")) %>%
  mutate(value = (value - true_val) / true_val)

ggplot(re_df %>% filter(Fleet == 1), aes(x = Year, y = value, color = Survey)) +
  geom_line(lwd = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(2018, 2025), lty = 2, color = 'black', lwd = 1 ) +
  labs(y = "Mean Relative Error in Survey Index") +
  theme_bw()

# Single region model starts off engativeyl biased because of carry over values
# from previous years prior to "observing" large recruitment events in
# other regions. It then switches to positive bias
# later on following new carry over values from large recruitment events

# This trend kind of manifests in all the models across cross-testing and MSE runs.
# Start of negatively biased in the first closed loop year and it turns to positive bias
# thereafter.   less pronounced trend in the faa likely because u dont fit those data, but the same trend still manifests because u dont have data in years in the bs, ai where there are high rec so it reverts closer to mean (hence goes to neg bias), and then as recruitment goes back to normal levels, slight
# positive upward bias because it takes a couple of years for the
# mean recruitment equilibriate back?
# kinda similar but different issue for the three_rg model but compounded w/ model instability etc

### General Conclusions -----------------------------------------------------

# FAA model generally predicts NAA well
# sgl model overpredicts NAA the second most
# three_rg model overpredicts NAA the most (most likely due to data availiability)

# Under different survey change scenarios ... current is better than historical and all b/c less regions sampled.
# W/ more regions sampled, more conflicting signals that you need to reconcile via recruitment and fishing mortality in
# a single population age-structure.

# For the sgl region model, similar story thought flipping back and forth BSAI and GOA is better than historical, because
# historical has GOA dominating, but the pop scale has switched from not necessarily being a GOA dominated system, so better
# to average two large areas rather tahn bias towards the GOA.

# For three_rg, genearlly as expected, more regions sampled = less bias because able to acocmodate that spatial structure there.

# For the FAA model, generally able to get the population scale correct because
# there really isn't any aggregation going on, and the catchabilities from the surveys
# are therefore time-invariant and less mis-specified. However, some biases still manifest, primarily
# because of FAA can't really account for time-varying availiability w/o time-varying selectivity (
# i.e., if a new recruitment event pops up in one given region). It can only account for spatial averages.

# For the sgl model, increasing pop scale bias, because q is time-varying from carry over values.

# For the three_rg model, bias primarily from overparameterization and trying to reconcile recruitment and movement.
# New recruits coming in changes movement, vice versa, and end up overpredicting biomass due to recruitment being so
# variable and uncertain. Also likely some auto-correlated spatially-varying recruitment process that is not being accounted for.
# Some bias likely from aggregating the GOA too. However, we did show with good data, a lot of these biases go away.


