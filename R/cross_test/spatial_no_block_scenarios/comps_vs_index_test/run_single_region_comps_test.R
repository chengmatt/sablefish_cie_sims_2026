# Purpose: Separate the effect of the LLS age comps from the carryover LLS index on single region EM bias.
#          Fits 8 single region variants per simulation (same observation noise within a sim), crossing the
#          LLS index (carryover vs. complete) with the LLS age comps (missing BS or AI, GOA only, BS only, complete).
# Creator: Matthew LH. Cheng
# Date: 11/10/26
# Notes: Requires SPoRC at commit 72b8426a (pre-seasonality) to match the five_region_base OM

# Setup -------------------------------------------------------------------
library(here)
library(SPoRC)
library(tidyverse)
library(parallel)

source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "single_region_em.R"))

sims <- 1:100     # simulations to run
ncores <- 3       # forked workers (mclapply)
Y <- 65           # terminal year (2024), i.e., an assessment run with data through 2024

out_dir <- here("outputs", "cross_test", "spatial_no_block_scenarios", "comps_vs_index_test")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

om <- readRDS(here("outputs", "mse_results", "spatial_no_block_scenarios", "five_region_base.RDS"))

# Helpers -----------------------------------------------------------------
# build single region data; the seed is reset per build so both designs share the same catch / index noise
build <- function(i, design) {
  set.seed(1000 + i)
  single_region_em(sim_env = om, y = Y, sim = i, srv_idx_se = 0.2, age_lag = 1,
                   lls_design_type = design, srv_wgt = "numbers", fish_wgt = "numbers")
}

# replace the LLS age comps with comps aggregated from a subset of regions (same survey-numbers weighting
# and stratified ISS as agg_data_to_single_rg); regs = 3:5 is GOA only, regs = 1 is BS only
region_comps <- function(a, i, regs) {
  sd <- simulation_data_to_SPoRC(om, Y, i)
  for (yr in 1:Y) {
    w <- om$TrueSrvIdx[regs, yr, 1, i]; if (sum(w) <= 0) next; w <- w / sum(w)
    rp <- sd$ObsSrvAgeComps[regs, yr, , , 1, drop = FALSE]
    rp <- array(rp, dim = c(length(regs), dim(rp)[3], dim(rp)[4]))      # region x age x sex
    n_r <- apply(rp, 1, sum); if (any(n_r <= 0)) next
    p_r <- sweep(rp, 1, n_r, "/")
    agg <- apply(sweep(p_r, 1, w, "*"), 2:3, sum)
    iss <- sum(agg * (1 - agg)) / sum(w^2 * (1 / n_r) * apply(p_r * (1 - p_r), 1, sum))
    a$data$ObsSrvAgeComps[1, yr, , , 1] <- agg / sum(agg)
    a$data$ISS_SrvAgeComps[1, yr, 1, 1] <- iss
  }
  a
}

summ <- function(fit) list(SSB = as.vector(fit$rep$SSB[1, 1:Y]),
                           NAA = apply(fit$rep$NAA[1, 1:Y, , ], 1:2, sum),   # year x age (sexes summed)
                           maxgrad = max(abs(fit$gr())))
fitit <- function(a) {
  f <- tryCatch(fit_model(a$data, a$par, a$map, NULL, 2, silent = TRUE), error = function(e) NULL)
  if (is.null(f)) NULL else summ(f)
}

# Run ---------------------------------------------------------------------
run_one <- function(i) {
  H  <- build(i, "historical")   # carryover index, comps missing BS or AI each year
  Al <- build(i, "all")          # complete index, complete comps
  swap_comps <- function(a, from) { for (k in c("ObsSrvAgeComps", "ISS_SrvAgeComps", "UseSrvAgeComps")) a$data[[k]] <- from$data[[k]]; a }
  swap_index <- function(a, from) { for (k in c("ObsSrvIdx", "UseSrvIdx", "ObsSrvIdx_SE")) a$data[[k]] <- from$data[[k]]; a }
  out <- list(
    hist                   = fitit(H),                          # carryover index + comps missing BS or AI (as run)
    full_idx_partial_comps = fitit(swap_index(H, Al)),          # complete index + comps missing BS or AI
    goa_comps_carry_idx    = fitit(region_comps(H, i, 3:5)),    # carryover index + GOA-only comps
    goa_comps_full_idx     = fitit(region_comps(Al, i, 3:5)),   # complete index + GOA-only comps
    bs_comps_carry_idx     = fitit(region_comps(H, i, 1)),      # carryover index + BS-only comps
    bs_comps_full_idx      = fitit(region_comps(Al, i, 1)),     # complete index + BS-only comps
    carry_idx_full_comps   = fitit(swap_comps(H, Al)),          # carryover index + complete comps
    all                    = fitit(Al)                          # complete index + complete comps
  )
  message("sim ", i, " done")
  out
}

res <- mclapply(sims, run_one, mc.cores = ncores)
names(res) <- sims
saveRDS(res, file.path(out_dir, "single_region_comps_test.RDS"))
