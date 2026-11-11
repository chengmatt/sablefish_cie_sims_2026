# Purpose: What if the historical LLS had followed the current alternating design? Refits the single region EM with,
#          from 1997 (model year 38), BS + AI surveyed in odd years and the GOA in even years (1990-1996 as historical).
#          Index: the unsurveyed half is carried over from the previous year, with the same observation noise as the
#          as-run index. Comps: only the regions surveyed that year (survey-numbers weights, stratified ISS).
#          Also fits the alternating index alone (comps as run) to isolate the carryover effect. Compare with `hist`
#          and `full_idx_partial_comps` in single_region_comps_test.RDS (same build and seed, so catch, fishery data
#          and noise are identical).
# Creator: Matthew LH. Cheng
# Date: 11/11/26
# Notes: Requires SPoRC at commit 72b8426a (pre-seasonality); see README.

# Setup -------------------------------------------------------------------
library(here)
library(SPoRC)
library(tidyverse)
library(parallel)

source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "single_region_em.R"))
source(here("R", "cross_test", "spatial_no_block_scenarios", "comps_vs_index_test", "index_helpers.R"))   # carry_hist(), carry_alt()

sims <- 1:100     # simulations to run
ncores <- 3       # forked workers (mclapply)
Y <- 65           # terminal year (2024)

out_dir <- here("outputs", "cross_test", "spatial_no_block_scenarios", "comps_vs_index_test")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

om <- readRDS(here("outputs", "mse_results", "spatial_no_block_scenarios", "five_region_base.RDS"))

# Helpers -----------------------------------------------------------------
# index: same observation noise as the as-run index, applied to the alternating carryover truth
alt_index <- function(a, i) {
  tr <- om$TrueSrvIdx[, 1:Y, , i]
  mult <- a$data$ObsSrvIdx[1, 1:Y, ] / apply(carry_hist(tr), 2:3, sum)
  a$data$ObsSrvIdx[1, 1:Y, ] <- apply(carry_alt(tr), 2:3, sum) * mult
  a
}
# comps: from 1997, only the regions surveyed that year (same weighting / ISS as agg_data_to_single_rg)
alt_comps <- function(a, i) {
  sd <- simulation_data_to_SPoRC(om, Y, i)
  for (yr in 38:Y) {
    regs <- if (bsai_year(yr)) 1:2 else 3:5
    w <- om$TrueSrvIdx[regs, yr, 1, i]; w <- w / sum(w)
    rp <- sd$ObsSrvAgeComps[regs, yr, , , 1, drop = FALSE]
    rp <- array(rp, dim = c(length(regs), dim(rp)[3], dim(rp)[4]))      # region x age x sex
    n_r <- apply(rp, 1, sum); if (any(n_r <= 0)) next
    p_r <- sweep(rp, 1, n_r, "/")
    agg <- apply(sweep(p_r, 1, w, "*"), 2:3, sum)
    a$data$ObsSrvAgeComps[1, yr, , , 1] <- agg / sum(agg)
    a$data$ISS_SrvAgeComps[1, yr, 1, 1] <- sum(agg * (1 - agg)) / sum(w^2 * (1 / n_r) * apply(p_r * (1 - p_r), 1, sum))
  }
  a
}

fitit <- function(a) {
  f <- tryCatch(fit_model(a$data, a$par, a$map, NULL, 2, silent = TRUE), error = function(e) NULL)
  if (is.null(f)) NULL else list(SSB = as.vector(f$rep$SSB[1, 1:Y]), NAA = apply(f$rep$NAA[1, 1:Y, , ], 1:2, sum),
                                 maxgrad = max(abs(f$gr())))
}

# Run ---------------------------------------------------------------------
run_one <- function(i) {
  set.seed(1000 + i)   # same build as `hist` in run_single_region_comps_test.R
  H <- single_region_em(sim_env = om, y = Y, sim = i, srv_idx_se = 0.2, age_lag = 1, lls_design_type = "historical",
                        srv_wgt = "numbers", fish_wgt = "numbers")
  out <- list(alt_both  = fitit(alt_comps(alt_index(H, i), i)),   # alternating index and comps (current design)
              alt_index = fitit(alt_index(H, i)))                 # alternating index only, comps as run (carryover effect)
  message("sim ", i, " done")
  out
}

res <- mclapply(sims, run_one, mc.cores = ncores)
names(res) <- sims
saveRDS(res, file.path(out_dir, "single_region_alt_design.RDS"))
