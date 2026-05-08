# Purpose: To self test the single region EM and FAA model w/ different aggregaiton schemes
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
source(here("R", "functions", "single_region_em.R"))
source(here("R", "functions", "faa_em.R"))

# Read in OMs
highsamp <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_highsamp.RDS"))
lowsamp <- readRDS(here("outputs", 'cross_test', "spatial_no_block_scenarios", "spt_rand_OM_lowsamp.RDS"))

# Simulation dimensions
y <- 65 # terminal year prior to MSE (no feedback)
n_sims <- 100 # number of sims


# FAA model setup ---------------------------------------------------------

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


# Perfect data single region -----------------------------------------------------------------
n_sims <- 50
i <- 1
# get single region data
asmt_list <- single_region_em(
  sim_env = highsamp,
  y = y,
  sim = 1,
  srv_idx_se = 0.01,
  age_lag = 1,
  lls_design_type = "all",
  srv_wgt = 'numbers',
  fish_wgt = 'numbers'
)

# True Comps
ae <- asmt_list$data$AgeingError
ae_kron <- function(yr) kronecker(diag(2), ae[yr,,])
caa_sum <- apply(highsamp$CAA, 2:6, sum)[1:y,,,,i]
saa_sum <- apply(highsamp$SrvIAA, 2:6, sum)[1:y,,,,i]
for(yr in 1:y) for(f in 1:2) asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- as.vector(caa_sum[yr,,,f]) %*% ae_kron(yr)
for(yr in 1:y) for(f in 1:3) asmt_list$data$ObsSrvAgeComps[1,yr,,,f]  <- as.vector(saa_sum[yr,,,f]) %*% ae_kron(yr)
asmt_list$data$ObsFishLenComps[] <- apply(highsamp$CAL, 2:6, sum)[1:y,,,,i]
asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]
asmt_list$data$ObsCatch[]  <- apply(highsamp$TrueCatch,  2:4, sum)[1:y,,i]
asmt_list$data$ISS_FishAgeComps[] <- 1e3
asmt_list$data$ISS_SrvAgeComps[]  <- 1e3
asmt_list$data$ISS_FishLenComps[] <- 1e3
asmt_list$data$ObsSrvIdx_SE[] <- 0.01
asmt_list$par$ln_sigmaC[] <- log(0.01)
asmt_list$data$UseSrvIdx

sgl_noerror <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
sgl_noerror$data <- asmt_list$data
sgl_noerror$sd_rep <- sdreport(sgl_noerror)

get_idx_fits_plot(list(asmt_list$data), list(sgl_noerror$rep),  1)

plot((t(sgl_noerror$rep$SSB) - (colSums(highsamp$SSB[,1:65,1]))) /
       (colSums(highsamp$SSB[,1:65,1])) * 100, type = 'l',
     ylim = c(-50, 50), ylab = "SSB Relative Error"
     )
abline(h = 0)

plot((t(sgl_noerror$rep$Total_Biom) - (colSums(highsamp$Total_Biom[,1:65,1]))) /
       (colSums(highsamp$Total_Biom[,1:65,1])) * 100, type = 'l'
     # ylim = c(-50, 50)
     )
abline(h = 0)

plot((t(sgl_noerror$rep$Rec) - (colSums(highsamp$Rec[,1:65,1]))) / (colSums(highsamp$Rec[,1:65,1])) * 100,
     type = 'l'
     # ylim = c(-50, 50)
     )
abline(h = 0, v = 60)

plot(
  (apply(sgl_noerror$rep$NAA[,1:65,,], 1, sum) - apply(highsamp$NAA[,1:65,,,1], 2, sum)) / apply(highsamp$NAA[,1:65,,,1], 2, sum),
  type = 'l'
)


plot(t(sgl_noerror$rep$Total_Biom))
lines((colSums(highsamp$Total_Biom[,-66,1])))

plot(t(sgl_noerror$rep$Rec))
lines((colSums(highsamp$Rec[,-66,1])))


# Perfect data single region - no age lag -----------------------------------------------------------------

asmt_list <- single_region_em(
  sim_env = highsamp,
  y = y,
  sim = 1,
  srv_idx_se = 0.1,
  age_lag = 0,
  lls_design_type = "all",
  srv_wgt = 'numbers',
  fish_wgt = 'numbers'
)

ae <- asmt_list$data$AgeingError
ae_kron <- function(yr) kronecker(diag(2), ae[yr,,])
caa_sum <- apply(highsamp$CAA, 2:6, sum)[1:y,,,,i]
saa_sum <- apply(highsamp$SrvIAA, 2:6, sum)[1:y,,,,i]
for(yr in 1:y) for(f in 1:2) asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- as.vector(caa_sum[yr,,,f]) %*% ae_kron(yr)
for(yr in 1:y) for(f in 1:3) asmt_list$data$ObsSrvAgeComps[1,yr,,,f]  <- as.vector(saa_sum[yr,,,f]) %*% ae_kron(yr)
asmt_list$data$ObsFishLenComps[] <- apply(highsamp$CAL, 2:6, sum)[1:y,,,,i]
asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]
asmt_list$data$ObsCatch[]  <- apply(highsamp$TrueCatch,  2:4, sum)[1:y,,i]
asmt_list$data$ISS_FishAgeComps[] <- 1e3
asmt_list$data$ISS_SrvAgeComps[]  <- 1e3
asmt_list$data$ISS_FishLenComps[] <- 1e3
asmt_list$data$ObsSrvIdx_SE[] <- 0.01
asmt_list$par$ln_sigmaC[] <- log(0.01)

sgl_noerror_noagelag <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
sgl_noerror_noagelag$data <- asmt_list$data
sgl_noerror_noagelag$sd_rep <- sdreport(sgl_noerror_noagelag)

plot((t(sgl_noerror_noagelag$rep$SSB) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, type = 'l', ylim = c(-50, 50))
abline(h = 0)

# Error in Index - single region -----------------------------------------------------------------
handlers(global = TRUE)
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

with_progress({
  p <- progressor(steps = n_sims)
  sgl_idx_error <- future_lapply(1:n_sims, function(i) {
    asmt_list <- single_region_em(
      sim_env = highsamp, y = y, sim = i, srv_idx_se = 0.2,
      age_lag = 1, lls_design_type = "historical",
      srv_wgt = 'numbers', fish_wgt = 'numbers'
    )
    ae <- asmt_list$data$AgeingError
    ae_kron <- function(yr) kronecker(diag(2), ae[yr,,])
    caa_sum <- apply(highsamp$CAA, 2:6, sum)[1:y,,,,i]
    saa_sum <- apply(highsamp$SrvIAA, 2:6, sum)[1:y,,,,i]
    for(yr in 1:y) for(f in 1:2) asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- as.vector(caa_sum[yr,,,f]) %*% ae_kron(yr)
    for(yr in 1:y) for(f in 1:3) asmt_list$data$ObsSrvAgeComps[1,yr,,,f]  <- as.vector(saa_sum[yr,,,f]) %*% ae_kron(yr)
    asmt_list$data$ObsFishLenComps[] <- apply(highsamp$CAL, 2:6, sum)[1:y,,,,i]
    asmt_list$data$ObsCatch[] <- apply(highsamp$TrueCatch, 2:4, sum)[1:y,,i]
    asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i] * exp(rnorm(length(apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]), 0, 0.2))
    asmt_list$data$ISS_FishAgeComps[] <- 1e3
    asmt_list$data$ISS_SrvAgeComps[]  <- 1e3
    asmt_list$data$ISS_FishLenComps[] <- 1e3
    asmt_list$data$ObsSrvIdx_SE[] <- 0.2
    asmt_list$par$ln_sigmaC[] <- log(0.01)
    out <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
    out$data <- asmt_list$data; out$sd_rep <- sdreport(out); p(); out
  })
})
saveRDS(sgl_idx_error, here("outputs", "cross_test", "spatial_no_block_scenarios", "sgl_idx_error_aggregation_test.RDS"))

plot((t(apply(simplify2array(lapply(sgl_idx_error, \(x) x$rep$SSB)), 1:2, median)) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-20,20))
abline(h = 0)

# Error in Index - single region - no age lag -----------------------------------------------------------------
handlers(global = TRUE)
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

with_progress({
  p <- progressor(steps = n_sims)
  sgl_idx_error_noagelag <- future_lapply(1:n_sims, function(i) {
    asmt_list <- single_region_em(
      sim_env = highsamp, y = y, sim = i, srv_idx_se = 0.2,
      age_lag = 0, lls_design_type = "historical",
      srv_wgt = 'numbers', fish_wgt = 'numbers'
    )
    ae <- asmt_list$data$AgeingError
    ae_kron <- function(yr) kronecker(diag(2), ae[yr,,])
    caa_sum <- apply(highsamp$CAA, 2:6, sum)[1:y,,,,i]
    saa_sum <- apply(highsamp$SrvIAA, 2:6, sum)[1:y,,,,i]
    for(yr in 1:y) for(f in 1:2) asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- as.vector(caa_sum[yr,,,f]) %*% ae_kron(yr)
    for(yr in 1:y) for(f in 1:3) asmt_list$data$ObsSrvAgeComps[1,yr,,,f]  <- as.vector(saa_sum[yr,,,f]) %*% ae_kron(yr)
    asmt_list$data$ObsFishLenComps[] <- apply(highsamp$CAL, 2:6, sum)[1:y,,,,i]
    asmt_list$data$ObsCatch[] <- apply(highsamp$TrueCatch, 2:4, sum)[1:y,,i]
    asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i] * exp(rnorm(length(apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]), 0, 0.2))
    asmt_list$data$ISS_FishAgeComps[] <- 1e3
    asmt_list$data$ISS_SrvAgeComps[]  <- 1e3
    asmt_list$data$ISS_FishLenComps[] <- 1e3
    asmt_list$data$ObsSrvIdx_SE[] <- 0.2
    asmt_list$par$ln_sigmaC[] <- log(0.01)
    out <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
    out$data <- asmt_list$data; out$sd_rep <- sdreport(out); p(); out
  })
})
saveRDS(sgl_idx_error_noagelag, here("outputs", "cross_test", "spatial_no_block_scenarios", "sgl_idx_error_noagelag_aggregation_test.RDS"))

plot((t(apply(simplify2array(lapply(sgl_idx_error_noagelag, \(x) x$rep$SSB)), 1:2, median)) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-20,20))
abline(h = 0)

# Error in comps only - single region -----------------------------------------------------------------
handlers(global = TRUE)
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

with_progress({
  p <- progressor(steps = n_sims)
  sgl_error_comps <- future_lapply(1:n_sims, function(i) {
    asmt_list <- single_region_em(
      sim_env = highsamp, y = y, sim = i, srv_idx_se = 0.2,
      age_lag = 1, lls_design_type = "historical",
      srv_wgt = 'numbers', fish_wgt = 'numbers'
    )
    caa_sum <- apply(highsamp$CAA,    2:6, sum)[1:y,,,,i]
    cal_sum <- apply(highsamp$CAL,    2:6, sum)[1:y,,,,i]
    saa_sum <- apply(highsamp$SrvIAA, 2:6, sum)[1:y,,,,i]
    ae <- asmt_list$data$AgeingError
    ae_kron <- function(yr) kronecker(diag(2), ae[yr,,])
    for(yr in 1:65) for(f in 1:2) asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- as.vector(rmultinom(1, 500, as.vector(caa_sum[yr,,,f]))) %*% ae_kron(yr)
    for(yr in 1:65) for(f in 1:2) asmt_list$data$ObsFishLenComps[1,yr,,,f] <- rmultinom(1, 500, as.vector(cal_sum[yr,,,f]))
    for(yr in 1:65) for(f in 1:3) asmt_list$data$ObsSrvAgeComps[1,yr,,,f]  <- as.vector(rmultinom(1, 500, as.vector(saa_sum[yr,,,f]))) %*% ae_kron(yr)
    asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]
    asmt_list$data$ObsCatch[]  <- apply(highsamp$TrueCatch,  2:4, sum)[1:y,,i]
    asmt_list$data$ISS_FishAgeComps[] <- 500
    asmt_list$data$ISS_SrvAgeComps[]  <- 500
    asmt_list$data$ISS_FishLenComps[] <- 500
    asmt_list$data$ObsSrvIdx_SE[] <- 0.01
    asmt_list$par$ln_sigmaC[] <- log(0.01)
    out <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
    out$data <- asmt_list$data; out$sd_rep <- sdreport(out); p(); out
  })
})
saveRDS(sgl_error_comps, here("outputs", "cross_test", "spatial_no_block_scenarios", "sgl_error_comps_aggregation_test.RDS"))

plot((t(apply(simplify2array(lapply(sgl_error_comps, \(x) x$rep$SSB)), 1:2, median)) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-20,20))
abline(h = 0)

# Error in comps only - single region - no age lag -----------------------------------------------------------------
handlers(global = TRUE)
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

with_progress({
  p <- progressor(steps = n_sims)
  sgl_error_comps_noagelag <- future_lapply(1:n_sims, function(i) {
    asmt_list <- single_region_em(
      sim_env = highsamp, y = y, sim = i, srv_idx_se = 0.2,
      age_lag = 0, lls_design_type = "historical",
      srv_wgt = 'numbers', fish_wgt = 'numbers'
    )
    caa_sum <- apply(highsamp$CAA,    2:6, sum)[1:y,,,,i]
    cal_sum <- apply(highsamp$CAL,    2:6, sum)[1:y,,,,i]
    saa_sum <- apply(highsamp$SrvIAA, 2:6, sum)[1:y,,,,i]
    ae <- asmt_list$data$AgeingError
    ae_kron <- function(yr) kronecker(diag(2), ae[yr,,])
    for(yr in 1:65) for(f in 1:2) asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- as.vector(rmultinom(1, 500, as.vector(caa_sum[yr,,,f]))) %*% ae_kron(yr)
    for(yr in 1:65) for(f in 1:2) asmt_list$data$ObsFishLenComps[1,yr,,,f] <- rmultinom(1, 500, as.vector(cal_sum[yr,,,f]))
    for(yr in 1:65) for(f in 1:3) asmt_list$data$ObsSrvAgeComps[1,yr,,,f]  <- as.vector(rmultinom(1, 500, as.vector(saa_sum[yr,,,f]))) %*% ae_kron(yr)
    asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]
    asmt_list$data$ObsCatch[]  <- apply(highsamp$TrueCatch,  2:4, sum)[1:y,,i]
    asmt_list$data$ISS_FishAgeComps[] <- 500
    asmt_list$data$ISS_SrvAgeComps[]  <- 500
    asmt_list$data$ISS_FishLenComps[] <- 500
    asmt_list$data$ObsSrvIdx_SE[] <- 0.01
    asmt_list$par$ln_sigmaC[] <- log(0.01)
    out <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
    out$data <- asmt_list$data; out$sd_rep <- sdreport(out); p(); out
  })
})
saveRDS(sgl_error_comps_noagelag, here("outputs", "cross_test", "spatial_no_block_scenarios", "sgl_error_comps_noagelag_aggregation_test.RDS"))

plot((t(apply(simplify2array(lapply(sgl_error_comps_noagelag, \(x) x$rep$SSB)), 1:2, median)) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-20,20))
abline(h = 0)

# Error in comps only (sim then agg) - single region -----------------------------------------------------------------
handlers(global = TRUE)
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

with_progress({
  p <- progressor(steps = n_sims)
  sgl_error_comps_simthenagg <- future_lapply(1:n_sims, function(i) {
    asmt_list <- single_region_em(
      sim_env = highsamp, y = y, sim = i, srv_idx_se = 0.2,
      age_lag = 1, lls_design_type = "historical",
      srv_wgt = 'numbers', fish_wgt = 'numbers'
    )
    n_per_reg <- 100
    get_strat_iss <- function(reg_draws, w_norm) {
      regional_n <- colSums(reg_draws)
      regional_p <- sweep(reg_draws, 2, regional_n, "/")
      pooled_p   <- rowSums(sweep(regional_p, 2, w_norm, "*"))
      sum(pooled_p * (1 - pooled_p)) / sum(w_norm^2 * (1 / regional_n) * colSums(regional_p * (1 - regional_p)))
    }
    for(yr in 1:65) for(f in 1:2) {
      w <- apply(highsamp$CAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) as.vector(rmultinom(1, n_per_reg, as.vector(highsamp$CAA[r,yr,,,f,i]))) %*% kronecker(diag(2), asmt_list$data$AgeingError[yr,,])), ncol = length(active))
      agg_props <- (reg_draws / n_per_reg) %*% w_norm; iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_FishAgeComps[1,yr,1,f] <- iss
    }
    for(yr in 1:65) for(f in 1:2) {
      w <- apply(highsamp$CAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) rmultinom(1, n_per_reg, as.vector(highsamp$CAL[r,yr,,,f,i]))), ncol = length(active))
      agg_props <- rowSums(reg_draws / n_per_reg * w_norm); iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsFishLenComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_FishLenComps[1,yr,1,f] <- iss
    }
    for(yr in 1:65) for(f in 1:3) {
      w <- apply(highsamp$SrvIAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) as.vector(rmultinom(1, n_per_reg, as.vector(highsamp$SrvIAA[r,yr,,,f,i]))) %*% kronecker(diag(2), asmt_list$data$AgeingError[yr,,])), ncol = length(active))
      agg_props <- (reg_draws / n_per_reg) %*% w_norm; iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsSrvAgeComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_SrvAgeComps[1,yr,1,f] <- iss
    }
    asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]
    asmt_list$data$ObsCatch[]  <- apply(highsamp$TrueCatch,  2:4, sum)[1:y,,i]
    asmt_list$data$ObsSrvIdx_SE[] <- 0.01
    asmt_list$par$ln_sigmaC[] <- log(0.01)
    out <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
    out$data <- asmt_list$data; out$sd_rep <- sdreport(out); p(); out
  })
})
saveRDS(sgl_error_comps_simthenagg, here("outputs", "cross_test", "spatial_no_block_scenarios", "sgl_error_comps_simthenagg_aggregation_test.RDS"))

plot((t(apply(simplify2array(lapply(sgl_error_comps_simthenagg, \(x) x$rep$SSB)), 1:2, median)) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-20,20))
abline(h = 0)

# Error in comps only (sim then agg) - single region - no age lag -----------------------------------------------------------------
handlers(global = TRUE)
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

with_progress({
  p <- progressor(steps = n_sims)
  sgl_error_comps_simthenagg_noagelag <- future_lapply(1:n_sims, function(i) {
    asmt_list <- single_region_em(
      sim_env = highsamp, y = y, sim = i, srv_idx_se = 0.2,
      age_lag = 0, lls_design_type = "historical",
      srv_wgt = 'numbers', fish_wgt = 'numbers'
    )
    n_per_reg <- 100
    get_strat_iss <- function(reg_draws, w_norm) {
      regional_n <- colSums(reg_draws)
      regional_p <- sweep(reg_draws, 2, regional_n, "/")
      pooled_p   <- rowSums(sweep(regional_p, 2, w_norm, "*"))
      sum(pooled_p * (1 - pooled_p)) / sum(w_norm^2 * (1 / regional_n) * colSums(regional_p * (1 - regional_p)))
    }
    for(yr in 1:65) for(f in 1:2) {
      w <- apply(highsamp$CAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) as.vector(rmultinom(1, n_per_reg, as.vector(highsamp$CAA[r,yr,,,f,i]))) %*% kronecker(diag(2), asmt_list$data$AgeingError[yr,,])), ncol = length(active))
      agg_props <- (reg_draws / n_per_reg) %*% w_norm; iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_FishAgeComps[1,yr,1,f] <- iss
    }
    for(yr in 1:65) for(f in 1:2) {
      w <- apply(highsamp$CAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) rmultinom(1, n_per_reg, as.vector(highsamp$CAL[r,yr,,,f,i]))), ncol = length(active))
      agg_props <- rowSums(reg_draws / n_per_reg * w_norm); iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsFishLenComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_FishLenComps[1,yr,1,f] <- iss
    }
    for(yr in 1:65) for(f in 1:3) {
      w <- apply(highsamp$SrvIAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) as.vector(rmultinom(1, n_per_reg, as.vector(highsamp$SrvIAA[r,yr,,,f,i]))) %*% kronecker(diag(2), asmt_list$data$AgeingError[yr,,])), ncol = length(active))
      agg_props <- (reg_draws / n_per_reg) %*% w_norm; iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsSrvAgeComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_SrvAgeComps[1,yr,1,f] <- iss
    }
    asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]
    asmt_list$data$ObsCatch[]  <- apply(highsamp$TrueCatch,  2:4, sum)[1:y,,i]
    asmt_list$data$ObsSrvIdx_SE[] <- 0.01
    asmt_list$par$ln_sigmaC[] <- log(0.01)
    out <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
    out$data <- asmt_list$data; out$sd_rep <- sdreport(out); p(); out
  })
})
saveRDS(sgl_error_comps_simthenagg_noagelag, here("outputs", "cross_test", "spatial_no_block_scenarios", "sgl_error_comps_simthenagg_noagelag_aggregation_test.RDS"))

plot((t(apply(simplify2array(lapply(sgl_error_comps_simthenagg_noagelag, \(x) x$rep$SSB)), 1:2, median)) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-20,20))
abline(h = 0)

# Error in all - single region -----------------------------------------------------------------
handlers(global = TRUE)
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

with_progress({
  p <- progressor(steps = n_sims)
  sgl_errorall <- future_lapply(1:n_sims, function(i) {
    asmt_list <- single_region_em(
      sim_env = highsamp, y = y, sim = i, srv_idx_se = 0.2,
      age_lag = 1, lls_design_type = "historical",
      srv_wgt = 'numbers', fish_wgt = 'numbers'
    )
    caa_sum <- apply(highsamp$CAA,    2:6, sum)[1:y,,,,i]
    cal_sum <- apply(highsamp$CAL,    2:6, sum)[1:y,,,,i]
    saa_sum <- apply(highsamp$SrvIAA, 2:6, sum)[1:y,,,,i]
    ae <- asmt_list$data$AgeingError
    ae_kron <- function(yr) kronecker(diag(2), ae[yr,,])
    for(yr in 1:65) for(f in 1:2) asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- as.vector(rmultinom(1, 500, as.vector(caa_sum[yr,,,f]))) %*% ae_kron(yr)
    for(yr in 1:65) for(f in 1:2) asmt_list$data$ObsFishLenComps[1,yr,,,f] <- rmultinom(1, 500, as.vector(cal_sum[yr,,,f]))
    for(yr in 1:65) for(f in 1:3) asmt_list$data$ObsSrvAgeComps[1,yr,,,f]  <- as.vector(rmultinom(1, 500, as.vector(saa_sum[yr,,,f]))) %*% ae_kron(yr)
    asmt_list$data$ObsCatch[] <- apply(highsamp$TrueCatch, 2:4, sum)[1:y,,i]
    asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i] * exp(rnorm(length(apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]), 0, 0.2))
    asmt_list$data$ISS_FishAgeComps[] <- 500
    asmt_list$data$ISS_SrvAgeComps[]  <- 500
    asmt_list$data$ISS_FishLenComps[] <- 500
    asmt_list$data$ObsSrvIdx_SE[] <- 0.2
    asmt_list$par$ln_sigmaC[] <- log(0.01)
    out <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
    out$data <- asmt_list$data; out$sd_rep <- sdreport(out); p(); out
  })
})
saveRDS(sgl_errorall, here("outputs", "cross_test", "spatial_no_block_scenarios", "sgl_errorall_aggregation_test.RDS"))

plot((t(apply(simplify2array(lapply(sgl_errorall, \(x) x$rep$SSB)), 1:2, median)) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-20,20))
abline(h = 0)

# Error in all - single region - no age lag -----------------------------------------------------------------
handlers(global = TRUE)
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

with_progress({
  p <- progressor(steps = n_sims)
  sgl_errorall_noagelag <- future_lapply(1:n_sims, function(i) {
    asmt_list <- single_region_em(
      sim_env = highsamp, y = y, sim = i, srv_idx_se = 0.2,
      age_lag = 0, lls_design_type = "historical",
      srv_wgt = 'numbers', fish_wgt = 'numbers'
    )
    caa_sum <- apply(highsamp$CAA,    2:6, sum)[1:y,,,,i]
    cal_sum <- apply(highsamp$CAL,    2:6, sum)[1:y,,,,i]
    saa_sum <- apply(highsamp$SrvIAA, 2:6, sum)[1:y,,,,i]
    ae <- asmt_list$data$AgeingError
    ae_kron <- function(yr) kronecker(diag(2), ae[yr,,])
    for(yr in 1:65) for(f in 1:2) asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- as.vector(rmultinom(1, 500, as.vector(caa_sum[yr,,,f]))) %*% ae_kron(yr)
    for(yr in 1:65) for(f in 1:2) asmt_list$data$ObsFishLenComps[1,yr,,,f] <- rmultinom(1, 500, as.vector(cal_sum[yr,,,f]))
    for(yr in 1:65) for(f in 1:3) asmt_list$data$ObsSrvAgeComps[1,yr,,,f]  <- as.vector(rmultinom(1, 500, as.vector(saa_sum[yr,,,f]))) %*% ae_kron(yr)
    asmt_list$data$ObsCatch[] <- apply(highsamp$TrueCatch, 2:4, sum)[1:y,,i]
    asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i] * exp(rnorm(length(apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]), 0, 0.2))
    asmt_list$data$ISS_FishAgeComps[] <- 500
    asmt_list$data$ISS_SrvAgeComps[]  <- 500
    asmt_list$data$ISS_FishLenComps[] <- 500
    asmt_list$data$ObsSrvIdx_SE[] <- 0.2
    asmt_list$par$ln_sigmaC[] <- log(0.01)
    out <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
    out$data <- asmt_list$data; out$sd_rep <- sdreport(out); p(); out
  })
})
saveRDS(sgl_errorall_noagelag, here("outputs", "cross_test", "spatial_no_block_scenarios", "sgl_errorall_noagelag_aggregation_test.RDS"))

plot((t(apply(simplify2array(lapply(sgl_errorall_noagelag, \(x) x$rep$SSB)), 1:2, median)) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-20,20))
abline(h = 0)

# Error in all (sim then agg) - single region -----------------------------------------------------------------
handlers(global = TRUE)
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

with_progress({
  p <- progressor(steps = n_sims)
  sgl_errorall_simthenagg <- future_lapply(1:n_sims, function(i) {
    asmt_list <- single_region_em(
      sim_env = highsamp, y = y, sim = i, srv_idx_se = 0.2,
      age_lag = 1, lls_design_type = "historical",
      srv_wgt = 'numbers', fish_wgt = 'numbers'
    )
    n_per_reg <- 100
    get_strat_iss <- function(reg_draws, w_norm) {
      regional_n <- colSums(reg_draws)
      regional_p <- sweep(reg_draws, 2, regional_n, "/")
      pooled_p   <- rowSums(sweep(regional_p, 2, w_norm, "*"))
      sum(pooled_p * (1 - pooled_p)) / sum(w_norm^2 * (1 / regional_n) * colSums(regional_p * (1 - regional_p)))
    }
    for(yr in 1:65) for(f in 1:2) {
      w <- apply(highsamp$CAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) as.vector(rmultinom(1, n_per_reg, as.vector(highsamp$CAA[r,yr,,,f,i]))) %*% kronecker(diag(2), asmt_list$data$AgeingError[yr,,])), ncol = length(active))
      agg_props <- (reg_draws / n_per_reg) %*% w_norm; iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_FishAgeComps[1,yr,1,f] <- iss
    }
    for(yr in 1:65) for(f in 1:2) {
      w <- apply(highsamp$CAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) rmultinom(1, n_per_reg, as.vector(highsamp$CAL[r,yr,,,f,i]))), ncol = length(active))
      agg_props <- rowSums(reg_draws / n_per_reg * w_norm); iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsFishLenComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_FishLenComps[1,yr,1,f] <- iss
    }
    for(yr in 1:65) for(f in 1:3) {
      w <- apply(highsamp$SrvIAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) as.vector(rmultinom(1, n_per_reg, as.vector(highsamp$SrvIAA[r,yr,,,f,i]))) %*% kronecker(diag(2), asmt_list$data$AgeingError[yr,,])), ncol = length(active))
      agg_props <- (reg_draws / n_per_reg) %*% w_norm; iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsSrvAgeComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_SrvAgeComps[1,yr,1,f] <- iss
    }
    asmt_list$data$ObsCatch[] <- apply(highsamp$TrueCatch, 2:4, sum)[1:y,,i]
    asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i] * exp(rnorm(length(apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]), 0, 0.2))
    asmt_list$data$ObsSrvIdx_SE[] <- 0.2
    asmt_list$par$ln_sigmaC[] <- log(0.01)
    out <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
    out$data <- asmt_list$data; out$sd_rep <- sdreport(out); p(); out
  })
})
saveRDS(sgl_errorall_simthenagg, here("outputs", "cross_test", "spatial_no_block_scenarios", "sgl_errorall_simthenagg_aggregation_test.RDS"))

plot((t(apply(simplify2array(lapply(sgl_errorall_simthenagg, \(x) x$rep$SSB)), 1:2, median)) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-20,20))
abline(h = 0)

# Error in all (sim then agg) - single region - no age lag -----------------------------------------------------------------
handlers(global = TRUE)
plan(multisession, workers = 15)
options(future.globals.maxSize = 15e9)

with_progress({
  p <- progressor(steps = n_sims)
  sgl_errorall_simthenagg_noagelag <- future_lapply(1:n_sims, function(i) {
    asmt_list <- single_region_em(
      sim_env = highsamp, y = y, sim = i, srv_idx_se = 0.2,
      age_lag = 0, lls_design_type = "historical",
      srv_wgt = 'numbers', fish_wgt = 'numbers'
    )
    n_per_reg <- 100
    get_strat_iss <- function(reg_draws, w_norm) {
      regional_n <- colSums(reg_draws)
      regional_p <- sweep(reg_draws, 2, regional_n, "/")
      pooled_p   <- rowSums(sweep(regional_p, 2, w_norm, "*"))
      sum(pooled_p * (1 - pooled_p)) / sum(w_norm^2 * (1 / regional_n) * colSums(regional_p * (1 - regional_p)))
    }
    for(yr in 1:65) for(f in 1:2) {
      w <- apply(highsamp$CAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) as.vector(rmultinom(1, n_per_reg, as.vector(highsamp$CAA[r,yr,,,f,i]))) %*% kronecker(diag(2), asmt_list$data$AgeingError[yr,,])), ncol = length(active))
      agg_props <- (reg_draws / n_per_reg) %*% w_norm; iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsFishAgeComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_FishAgeComps[1,yr,1,f] <- iss
    }
    for(yr in 1:65) for(f in 1:2) {
      w <- apply(highsamp$CAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) rmultinom(1, n_per_reg, as.vector(highsamp$CAL[r,yr,,,f,i]))), ncol = length(active))
      agg_props <- rowSums(reg_draws / n_per_reg * w_norm); iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsFishLenComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_FishLenComps[1,yr,1,f] <- iss
    }
    for(yr in 1:65) for(f in 1:3) {
      w <- apply(highsamp$SrvIAA[,yr,,,f,i], 1, sum); active <- which(w > 0); w_norm <- w[active] / sum(w[active])
      reg_draws <- matrix(sapply(active, function(r) as.vector(rmultinom(1, n_per_reg, as.vector(highsamp$SrvIAA[r,yr,,,f,i]))) %*% kronecker(diag(2), asmt_list$data$AgeingError[yr,,])), ncol = length(active))
      agg_props <- (reg_draws / n_per_reg) %*% w_norm; iss <- get_strat_iss(reg_draws, w_norm)
      asmt_list$data$ObsSrvAgeComps[1,yr,,,f] <- agg_props * iss
      asmt_list$data$ISS_SrvAgeComps[1,yr,1,f] <- iss
    }
    asmt_list$data$ObsCatch[] <- apply(highsamp$TrueCatch, 2:4, sum)[1:y,,i]
    asmt_list$data$ObsSrvIdx[] <- apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i] * exp(rnorm(length(apply(highsamp$TrueSrvIdx, 2:4, sum)[1:y,,i]), 0, 0.2))
    asmt_list$data$ObsSrvIdx_SE[] <- 0.2
    asmt_list$par$ln_sigmaC[] <- log(0.01)
    out <- fit_model(asmt_list$data, asmt_list$par, asmt_list$map, NULL, 2, silent = FALSE)
    out$data <- asmt_list$data; out$sd_rep <- sdreport(out); p(); out
  })
})
saveRDS(sgl_errorall_simthenagg_noagelag, here("outputs", "cross_test", "spatial_no_block_scenarios", "sgl_errorall_simthenagg_noagelag_aggregation_test.RDS"))

plot((t(apply(simplify2array(lapply(sgl_errorall_simthenagg_noagelag, \(x) x$rep$SSB)), 1:2, median)) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-20,20))
abline(h = 0)

# Single Region Model Summaries -------------------------------------------
med_ssb <- function(lst) as.vector(t(apply(simplify2array(lapply(lst, \(x) x$rep$SSB)), 1:2, median)))

true_ssb <- colSums(highsamp$SSB[,-66,1])
years    <- 1:65

re <- list(
  "No Error"                          = as.vector(t(sgl_noerror$rep$SSB)),
  "No Error (no age lag)"             = as.vector(t(sgl_noerror_noagelag$rep$SSB)),
  "Index Error"                       = med_ssb(sgl_idx_error),
  "Index Error (no age lag)"          = med_ssb(sgl_idx_error_noagelag),
  "Comp Error"                        = med_ssb(sgl_error_comps),
  "Comp Error (no age lag)"           = med_ssb(sgl_error_comps_noagelag),
  "Comp Error (Sim\u2192Agg)"         = med_ssb(sgl_error_comps_simthenagg),
  "Comp Error (Sim\u2192Agg, no lag)" = med_ssb(sgl_error_comps_simthenagg_noagelag),
  "All Error"                         = med_ssb(sgl_errorall),
  "All Error (no age lag)"            = med_ssb(sgl_errorall_noagelag),
  "All Error (Sim\u2192Agg)"          = med_ssb(sgl_errorall_simthenagg),
  "All Error (Sim\u2192Agg, no lag)"  = med_ssb(sgl_errorall_simthenagg_noagelag)
)
re <- lapply(re, function(x) (x - true_ssb) / true_ssb * 100)

cols <- rep(c("black", "steelblue", "firebrick", "darkorange", "purple", "darkgreen"), each = 2)
ltys <- rep(c(1, 2), 6)

par(mar = c(4, 4.5, 1.5, 1))
plot(NULL, xlim = range(years), ylim = c(-50, 50),
     xlab = "Year", ylab = "Relative Error (%)", las = 1, bty = 'l')
abline(h = 0,         lty = 2, col = "grey60")
abline(h = c(-20, 20), lty = 3, col = "grey80")
for(j in seq_along(re)) lines(years, re[[j]], col = cols[j], lwd = 2, lty = ltys[j])
legend("topright", legend = names(re), col = cols, lwd = 2, lty = ltys,
       bty = 'n', cex = 0.65, ncol = 2)

# Perfect data FAA - Equal Weighting -----------------------------------------------------------------
# needs fully time-varying survey selectivity to get to 0 bias
j <- 19
i <- 1
asmt_list <- faa_em(
  sim_env = highsamp,
  y = y,
  sim = 1,
  srv_idx_se = 0.01,
  age_lag = 0,
  lls_design_type = 'all',
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


asmt_list$par$ln_fish_fixed_sel_pars[,,,,1] <- log(10)
# asmt_list$par$ln_fish_fixed_sel_pars[,,,,4] <- log(5)
asmt_list$par$ln_srv_fixed_sel_pars[] <- log(2)

# True composition data
# asmt_list$data$AgeingError[] <- aperm(replicate(65, diag(1, 30)), c(3,1,2))
ae <- asmt_list$data$AgeingError
ae_kron <- function(yr) kronecker(diag(2), ae[yr,,])
n_ages <- dim(highsamp$CAA)[3]
n_sexes <- dim(highsamp$CAA)[4]

ae <- asmt_list$data$AgeingError  # [yr, n_age, n_age]

apply_ae_yas <- function(x_yas) {
  # x_yas: [yr, age, sex]
  out <- x_yas
  for(yr in 1:dim(x_yas)[1])
    for(s in 1:dim(x_yas)[3])
      out[yr,,s] <- as.vector(x_yas[yr,,s]) %*% ae[yr,,]
  out
}

# Fishery age comps
asmt_list$data$ObsFishAgeComps[,,,,1] <- apply_ae_yas(highsamp$CAA[1,1:y,,,1,i])
asmt_list$data$ObsFishAgeComps[,,,,2] <- apply_ae_yas(highsamp$CAA[2,1:y,,,1,i])
asmt_list$data$ObsFishAgeComps[,,,,3] <- apply_ae_yas(apply(highsamp$CAA[3:5,1:y,,,,i,drop=F], 2:6, sum)[1:y,,,1,i])

# Length comps - no AE
asmt_list$data$ObsFishLenComps[] <- apply(highsamp$CAL, 2:6, sum)[1:y,,,,i]

# Survey age comps
asmt_list$data$ObsSrvAgeComps[,,,,1] <- apply_ae_yas(highsamp$SrvIAA[1,1:y,,,1,i])
asmt_list$data$ObsSrvAgeComps[,,,,2] <- apply_ae_yas(highsamp$SrvIAA[2,1:y,,,1,i])
asmt_list$data$ObsSrvAgeComps[,,,,3] <- apply_ae_yas(apply(highsamp$SrvIAA[3:5,1:y,,,,i,drop=F], 2:6, sum)[1:y,,,1,i])
asmt_list$data$ObsSrvAgeComps[,,,,4] <- apply_ae_yas(apply(highsamp$SrvIAA, 2:6, sum)[1:y,,,3,i])

# # True survey index data
asmt_list$data$ObsSrvIdx[,,1] <- highsamp$TrueSrvIdx[1,1:y,1,i]
asmt_list$data$ObsSrvIdx[,,2] <- highsamp$TrueSrvIdx[2,1:y,1,i]
asmt_list$data$ObsSrvIdx[,,3] <- colSums(highsamp$TrueSrvIdx[3:5,1:y,1,i])
asmt_list$data$ObsSrvIdx[,,4] <- colSums(highsamp$TrueSrvIdx[,1:y,3,i])
#
# # # True catch
asmt_list$data$ObsCatch[,,1] <- highsamp$TrueCatch[1,1:y,1,i]
asmt_list$data$ObsCatch[,,2] <- highsamp$TrueCatch[2,1:y,1,i]
asmt_list$data$ObsCatch[,,3] <- colSums(highsamp$TrueCatch[3:5,1:y,1,i])
asmt_list$data$ObsCatch[,,4] <- colSums(highsamp$TrueCatch[,1:y,2,i])

# # Change uncertainty
asmt_list$data$ISS_FishAgeComps[] <- 1e3 / 3
asmt_list$data$ISS_SrvAgeComps[]  <- 1e3 / 3
asmt_list$data$ISS_FishLenComps[] <- 1e3
asmt_list$data$ObsSrvIdx_SE[] <- 0.01
asmt_list$par$ln_sigmaC[] <- log(0.01)

# fit model
faa_eq <- SPoRC::fit_model(
  asmt_list$data,
  asmt_list$par,
  asmt_list$map,
  NULL,
  newton_loops = 3
  # 3,
  # silent = FALSE, do_optim = T
)

faa_eq$data <- asmt_list$data
faa_eq$sd_rep <- sdreport(faa_eq)

# Look at SSB comparison
plot((t(faa_eq$rep$SSB) - (colSums(highsamp$SSB[,-66,1]))) / (colSums(highsamp$SSB[,-66,1])) * 100, ylim = c(-50,50), type = 'l', ylab = 'SSB Relative Error')
abline(h = 0)

