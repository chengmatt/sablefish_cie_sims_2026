# Purpose: To valiidate and understand potential sources of bias from aggregating
# spatial outputs to a single-region
# Creator: Matthew LH. Cheng
# Date: 12/8/25

# Setup -------------------------------------------------------------------
library(here)
library(SPoRC)
library(tidyverse)

# Read in runs
model_list_lowsamp <- readRDS(here("outputs", "cross_test", "single_region_crosstest_lowsamp.RDS"))
lowsamp <- readRDS(here("outputs", 'cross_test', "spt_rand_OM_lowsamp.RDS"))

# Reverse model list to match format
n_sims <- length(model_list_lowsamp)
n_models <- length(model_list_lowsamp[[1]])
model_list_lowsamp_reversed <- lapply(1:n_models, function(model) {
  lapply(1:n_sims, function(sim) {
    model_list_lowsamp[[sim]][[model]]
  })
})

em_names <- c("srv_biomass_fish_biomass", "srv_numbers_fish_biomass", 'srv_biomass_fish_numbers',
              'srv_numbers_fish_numbers', 'eqwt')


# Load in functions
source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "single_region_em.R"))

# Compare Assessments ------------------------------------------------------
sgl <- readRDS("/Users/matthewcheng/Desktop/PostDoc/2025-Sablefish-SAFE-Dev/Sept PT Model Runs/25_12_Drop_TS_Upd_M/25_12_Drop_TS_Upd_M_model_results_v2.RDS")
spt <- readRDS("/Users/matthewcheng/Desktop/PostDoc/2025-Sablefish-SAFE-Dev/Sept PT Model Runs/25_15_Spatial/Spatial_MltRel/Spatial_MltRel_model_results.RDS")

# SSB
plot(t(sgl$rep$SSB), ylim=c(0,300), type = 'l')
lines(colSums(spt$rep$SSB), lty = 2)
plot((t(sgl$rep$SSB) - colSums(spt$rep$SSB)) / colSums(spt$rep$SSB))

# Recruitment
plot(t(sgl$rep$Rec), ylim=c(0,100), type = 'l')
lines(colSums(spt$rep$Rec), lty = 2)
plot((t(sgl$rep$Rec) - colSums(spt$rep$Rec)) / colSums(spt$rep$Rec), type = 'l')

# Compare Selectivities ---------------------------------------------------

### Fishery -----------------------------------------------------------------

# define iterations
gears <- list(
  list(idx = 1, name = "Fixed-Gear"),
  list(idx = 2, name = "Trawl-Gear")
)

# blocks
blocks <- list(
  list(idx = 1, year = 1),
  list(idx = 2, year = 40),
  list(idx = 3, year = 65)
)

# sexes
sexes <- list(
  list(idx = 1, label = "F"),
  list(idx = 2, label = "M")
)

# get selex
models_list <- list(model_list_lowsamp_reversed[[1]], model_list_lowsamp_reversed[[2]],
                    model_list_lowsamp_reversed[[3]], model_list_lowsamp_reversed[[4]],
                    model_list_lowsamp_reversed[[5]])

# loop through all
par(mfrow = c(2,4))
for (gear in gears) {
  for (block in blocks) {
    for (sex in sexes) {

      # Extract selectivity from all models in each set, median across sims
      sel_list <- lapply(models_list, function(model_set) {
        sel_matrix <- do.call(rbind, lapply(model_set, function(m) {
          m$rep$fish_sel[1, block$year, , sex$idx, gear$idx]
        }))
        apply(sel_matrix, 2, median, na.rm = TRUE)
      })


      # Get true values
      true_sel <- lowsamp$fish_sel[1, block$year, , sex$idx, gear$idx, 1]

      # Plot
      if(!(gear$idx == 2 && block$idx %in% c(2,3))) {
        plot(true_sel, ylim = c(0, 1), col = 'black', lty = 1, type = 'l',
             main = sprintf('%s Block %d - %s', gear$name, block$idx, sex$label),
             ylab = 'Median Sel', xlab = 'Age', lwd = 4)
        lines(sel_list[[1]], col = 'green', lwd = 4, lty = 2)  # srv_biomass_fish_biomass
        lines(sel_list[[2]], col = 'blue', lwd = 4, lty = 3)   # srv_numbers_fish_biomass
        lines(sel_list[[3]], col = 'purple', lty = 4, lwd = 4)             # srv_biomass_fish_numbers
        lines(sel_list[[4]], col = 'orange', lty = 4, lwd = 4)             # srv_numbers_fish_numbers
        lines(sel_list[[5]], col = 'brown', lty = 4, lwd = 4)             # eqwt

      }
    }
  }
}

# Add legend
legend('bottomleft',
       legend = c('total at age', em_names),
       col = c('black', c('green', 'blue', 'purple', 'orange', 'brown')),
       lty = c(1, c(2:4)),
       lwd = 4,
       cex = 0.8)


### Survey ------------------------------------------------------

# define iterations
gears <- list(
  list(idx = 1, name = "Domestic"),
  list(idx = 3, name = "JP")
)

# blocks
blocks <- list(
  list(idx = 1, year = 1),
  list(idx = 2, year = 40),
  list(idx = 3, year = 65)
)

# sexes
sexes <- list(
  list(idx = 1, label = "F"),
  list(idx = 2, label = "M")
)

# get selex
models_list <- list(model_list_lowsamp_reversed[[1]], model_list_lowsamp_reversed[[2]],
                    model_list_lowsamp_reversed[[3]], model_list_lowsamp_reversed[[4]],
                    model_list_lowsamp_reversed[[5]])

# loop through all
par(mfrow = c(2,4))
for (gear in gears) {
  for (block in blocks) {
    for (sex in sexes) {

      # Extract selectivity from all models in each set, median across sims
      sel_list <- lapply(models_list, function(model_set) {
        sel_matrix <- do.call(rbind, lapply(model_set, function(m) {
          m$rep$srv_sel[1, block$year, , sex$idx, gear$idx]
        }))
        apply(sel_matrix, 2, median, na.rm = TRUE)
      })


      # Get true values
      true_sel <- lowsamp$srv_sel[1, block$year, , sex$idx, gear$idx, 1]

      # Plot
      if(!(gear$idx == 3 && block$idx %in% c(2,3))) {
        plot(true_sel, ylim = c(0, 1), col = 'black', lty = 1, type = 'l',
             main = sprintf('%s Block %d - %s', gear$name, block$idx, sex$label),
             ylab = 'Median Sel', xlab = 'Age', lwd = 4)
        lines(sel_list[[1]], col = 'green', lwd = 4, lty = 2)  # srv_biomass_fish_biomass
        lines(sel_list[[2]], col = 'blue', lwd = 4, lty = 3)   # srv_numbers_fish_biomass
        lines(sel_list[[3]], col = 'purple', lty = 4, lwd = 4)             # srv_biomass_fish_numbers
        lines(sel_list[[4]], col = 'orange', lty = 4, lwd = 4)             # srv_numbers_fish_numbers
        lines(sel_list[[5]], col = 'brown', lty = 4, lwd = 4)             # eqwt
      }
    }
  }
}

# Add legend
legend('bottomright',
       legend = c('total at age', em_names),
       col = c('black', c('green', 'blue', 'purple', 'orange', 'brown')),
       lty = c(1, c(2:4)),
       lwd = 4,
       cex = 0.8)


# Index Fits --------------------------------------------------------------
sim <- 1
model_list_lowsamp_reversed[[1]][[sim]]$data$ObsSrvIdx[model_list_lowsamp_reversed[[1]][[sim]]$data$UseSrvIdx == 0] <- NA
get_idx_fits_plot(list(model_list_lowsamp_reversed[[1]][[sim]]$data, model_list_lowsamp_reversed[[1]][[sim]]$data,
                       model_list_lowsamp_reversed[[1]][[sim]]$data, model_list_lowsamp_reversed[[1]][[sim]]$data,
                       model_list_lowsamp_reversed[[1]][[sim]]$data),
                  list(model_list_lowsamp_reversed[[1]][[sim]]$rep, model_list_lowsamp_reversed[[2]][[sim]]$rep,
                       model_list_lowsamp_reversed[[3]][[sim]]$rep, model_list_lowsamp_reversed[[4]][[sim]]$rep,
                       model_list_lowsamp_reversed[[4]][[sim]]$rep), em_names) +
  theme(legend.position = 'right')

get_catch_fits_plot(list(model_list_lowsamp_reversed[[1]][[sim]]$data, model_list_lowsamp_reversed[[1]][[sim]]$data,
                       model_list_lowsamp_reversed[[1]][[sim]]$data, model_list_lowsamp_reversed[[1]][[sim]]$data,
                       model_list_lowsamp_reversed[[1]][[sim]]$data),
                  list(model_list_lowsamp_reversed[[1]][[sim]]$rep, model_list_lowsamp_reversed[[2]][[sim]]$rep,
                       model_list_lowsamp_reversed[[3]][[sim]]$rep, model_list_lowsamp_reversed[[4]][[sim]]$rep,
                       model_list_lowsamp_reversed[[4]][[sim]]$rep), em_names) +
  theme(legend.position = 'right')



# Compare Numbers at Age ----------------------------------------------------
# Function to extract and normalize naa
get_norm_naa <- function(model_obj) {
  naa <- apply(model_obj$rep$NAA[,,,], c(2,3), mean)
  naa_vec <- as.vector(naa)
  # naa_vec / sum(naa_vec)
}

# Calculate average across all sims for each model
n_models <- length(model_list_lowsamp_reversed)
n_sims <- length(model_list_lowsamp_reversed[[1]])

# get true aggregated numbers
true_naa <- apply(lowsamp$NAA[,1:65,,,], 2:5, sum)
true_naa_mean <- apply(true_naa, c(2,3), mean)
colors <- c('green', 'blue', 'purple', 'orange', 'brown')

par(mfrow = c(1,2))
plot(as.vector(true_naa_mean), type = 'l', lty = 2, lwd = 4, ylab = 'Aggregated Numbers at Age')
for (model in 1:n_models) {
  # Average across all sims
  naa_all_sims <- lapply(1:n_sims, function(sim) {
    get_norm_naa(model_list_lowsamp_reversed[[model]][[sim]])
  })

  # Take mean across sims
  naa_mean <- Reduce("+", naa_all_sims) / n_sims

  lines(naa_mean, lwd = 4, type = 'l', col = colors[model])
}

plot(rep(0, length(as.vector(true_naa_mean))), type = 'l', lty = 2, lwd = 4, ylab = 'Aggregated Numbers at Age', ylim = c(-0.5, 0.5))
for (model in 1:n_models) {
  # Average across all sims
  naa_all_sims <- lapply(1:n_sims, function(sim) {
    get_norm_naa(model_list_lowsamp_reversed[[model]][[sim]])
  })

  # Take mean across sims
  naa_mean <- Reduce("+", naa_all_sims) / n_sims

  lines((naa_mean - as.vector(true_naa_mean)) / as.vector(true_naa_mean), lwd = 4, type = 'l', col = colors[model])
}


# Add legend
legend('topright',
       legend = c('total at age', em_names),
       col = c('black', c('green', 'blue', 'purple', 'orange', 'brown')),
       lty = c(1, c(2:4)),
       lwd = 4,
       cex = 0.85)

# Function to extract and normalize caa
get_norm_caa <- function(model_obj) {
  caa <- apply(model_obj$rep$CAA[,,,,1], c(2,3), mean)
  caa_vec <- as.vector(caa)
  # caa_vec / sum(caa_vec)
}

# get true aggregated numbers
true_caa <- apply(lowsamp$CAA[,1:65,,,1,], 2:5, sum)
true_caa_mean <- apply(true_caa, c(2,3), mean)
colors <- c('green', 'blue', 'purple', 'orange', 'brown')

par(mfrow = c(1,2))
plot(as.vector(true_caa_mean), type = 'l', lty = 2, lwd = 4, ylab = 'Aggregated Catch at Age')
for (model in 1:n_models) {
  # Average across all sims
  caa_all_sims <- lapply(1:n_sims, function(sim) {
    get_norm_caa(model_list_lowsamp_reversed[[model]][[sim]])
  })

  # Take mean across sims
  caa_mean <- Reduce("+", caa_all_sims) / n_sims

  lines(caa_mean, lwd = 4, type = 'l', col = colors[model])
}

plot(rep(0, length(as.vector(true_caa_mean))), type = 'l', lty = 2, lwd = 4, ylab = 'Aggregated Catch at Age', ylim = c(-0.5, 0.5))
for (model in 1:n_models) {
  # Average across all sims
  caa_all_sims <- lapply(1:n_sims, function(sim) {
    get_norm_caa(model_list_lowsamp_reversed[[model]][[sim]])
  })

  # Take mean across sims
  caa_mean <- Reduce("+", caa_all_sims) / n_sims

  lines((caa_mean - as.vector(true_caa_mean)) / as.vector(true_caa_mean), lwd = 4, type = 'l', col = colors[model])
}


# Add legend
legend('bottomright',
       legend = c('total at age', em_names),
       col = c('black', c('green', 'blue', 'purple', 'orange', 'brown')),
       lty = c(1, c(2:4)),
       lwd = 4,
       cex = 0.85)

# Function to extract and normalize srvaa
get_norm_srvaa <- function(model_obj) {
  srvaa <- apply(model_obj$rep$SrvIAA[,,,,1], c(2,3), mean)
  srvaa_vec <- as.vector(srvaa)
  # srvaa_vec / sum(srvaa_vec)
}

# get true aggregated numbers
true_srvaa <- apply(lowsamp$SrvIAA[,1:65,,,1,], 2:5, sum)
true_srvaa_mean <- apply(true_srvaa, c(2,3), mean)
colors <- c('green', 'blue', 'purple', 'orange', 'brown')

par(mfrow = c(1,2))
plot(as.vector(true_srvaa_mean), type = 'l', lty = 2, lwd = 4, ylab = 'Aggregated Survey at Age')
for (model in 1:n_models) {
  # Average across all sims
  srvaa_all_sims <- lapply(1:n_sims, function(sim) {
    get_norm_srvaa(model_list_lowsamp_reversed[[model]][[sim]])
  })

  # Take mean across sims
  srvaa_mean <- Reduce("+", srvaa_all_sims) / n_sims

  lines(srvaa_mean, lwd = 4, type = 'l', col = colors[model])
}

plot(rep(0, length(as.vector(true_srvaa_mean))), type = 'l', lty = 2, lwd = 4, ylab = 'Aggregated Survey at Age', ylim = c(-0.5, 0.5))
for (model in 1:n_models) {
  # Average across all sims
  srvaa_all_sims <- lapply(1:n_sims, function(sim) {
    get_norm_srvaa(model_list_lowsamp_reversed[[model]][[sim]])
  })

  # Take mean across sims
  srvaa_mean <- Reduce("+", srvaa_all_sims) / n_sims

  lines((srvaa_mean - as.vector(true_srvaa_mean)) / as.vector(true_srvaa_mean), lwd = 4, type = 'l', col = colors[model])
}


# Add legend
legend('bottomright',
       legend = c('total at age', em_names),
       col = c('black', c('green', 'blue', 'purple', 'orange', 'brown')),
       lty = c(1, c(2:4)),
       lwd = 4,
       cex = 0.85)

# Function to extract and normalize caa
get_norm_caa <- function(model_obj) {
  caa <- apply(model_obj$data$ObsFishAgeComps[,1:65,,,1], c(2,3), mean)
  caa_vec <- as.vector(caa)
  caa_vec / sum(caa_vec)
}

par(mfrow = c(1,1))

# get true aggregated numbers
true_caa <- apply(lowsamp$CAA[,1:65,,,1,], 2:5, sum)
true_caa_mean <- apply(true_caa, c(2,3), mean)
true_caa_mean <- true_caa_mean / sum(true_caa_mean)
colors <- c('green', 'blue', 'purple', 'orange', 'brown')

plot(as.vector(true_caa_mean), type = 'l', lty = 2, lwd = 4, ylab = 'Aggregated Obs Fishery at Age', col = 'black')
for (model in 1:n_models) {
  # Average across all sims
  caa_all_sims <- lapply(1:n_sims, function(sim) {
    get_norm_caa(model_list_lowsamp_reversed[[model]][[sim]])
  })

  # Take mean across sims
  caa_mean <- Reduce("+", caa_all_sims) / n_sims
  lines(caa_mean, lwd = 4, type = 'l', col = colors[model])
}

# Add legend
legend('topright',
       legend = c('total at age', em_names),
       col = c('black', c('green', 'blue', 'purple', 'orange', 'brown')),
       lty = c(1, c(2:4)),
       lwd = 4,
       cex = 0.85)

# Function to extract and normalize srvaa
get_norm_srvaa <- function(model_obj) {
  srvaa <- apply(model_obj$data$ObsSrvAgeComps[,1:64,,,1], c(2,3), mean)
  srvaa_vec <- as.vector(srvaa)
  srvaa_vec / sum(srvaa_vec)
}

# get true aggregated numbers
true_srvaa <- apply(lowsamp$SrvIAA[,1:65,,,1,], 2:5, sum)
true_srvaa_mean <- apply(true_srvaa, c(2,3), mean)
true_srvaa_mean <- true_srvaa_mean / sum(true_srvaa_mean)
colors <- c('green', 'blue', 'purple', 'orange', 'brown')

plot(as.vector(true_srvaa_mean), type = 'l', lty = 2, lwd = 4, ylab = 'Aggregated Obs Survey at Age', col = 'black')
for (model in 1:n_models) {
  # Average across all sims
  srvaa_all_sims <- lapply(1:n_sims, function(sim) {
    get_norm_srvaa(model_list_lowsamp_reversed[[model]][[sim]])
  })

  # Take mean across sims
  srvaa_mean <- Reduce("+", srvaa_all_sims) / n_sims
  lines(srvaa_mean, lwd = 4, type = 'l', col = colors[model])
}

# Add legend
legend('topright',
       legend = c('total at age', em_names),
       col = c('black', c('green', 'blue', 'purple', 'orange', 'brown')),
       lty = c(1, c(2:4)),
       lwd = 4,
       cex = 0.85)

# Function to extract and normalize caa
get_norm_caa <- function(model_obj) {
  caa <- apply(model_obj$data$ObsFishAgeComps[,1:65,,,1] * model_obj$data$WAA_fish[,1:65,,,1], c(2,3), mean)
  caa_vec <- as.vector(caa)
  caa_vec / sum(caa_vec)
}

par(mfrow = c(1,1))

# get true aggregated numbers
true_caa <- apply(lowsamp$CAA[,1:65,,,1,] * lowsamp$WAA_fish[,1:65,,,1,], 2:5, sum)
true_caa_mean <- apply(true_caa, c(2,3), mean)
true_caa_mean <- true_caa_mean / sum(true_caa_mean)
colors <- c('green', 'blue', 'purple', 'orange', 'brown')

plot(as.vector(true_caa_mean), type = 'l', lty = 2, lwd = 4, ylab = 'Aggregated Obs Fishery at Age Biom', col = 'black')
for (model in 1:n_models) {
  # Average across all sims
  caa_all_sims <- lapply(1:n_sims, function(sim) {
    get_norm_caa(model_list_lowsamp_reversed[[model]][[sim]])
  })

  # Take mean across sims
  caa_mean <- Reduce("+", caa_all_sims) / n_sims
  lines(caa_mean, lwd = 4, type = 'l', col = colors[model])
}

# Add legend
legend('topright',
       legend = c('total at age', em_names),
       col = c('black', c('green', 'blue', 'purple', 'orange', 'brown')),
       lty = c(1, c(2:4)),
       lwd = 4,
       cex = 0.85)

# Function to extract and normalize srvaa
get_norm_srvaa <- function(model_obj) {
  srvaa <- apply(model_obj$data$ObsSrvAgeComps[,1:64,,,1] * model_obj$data$WAA_srv[,1:64,,,1], c(2,3), mean)
  srvaa_vec <- as.vector(srvaa)
  srvaa_vec / sum(srvaa_vec)
}

# get true aggregated numbers
true_srvaa <- apply(lowsamp$SrvIAA[,1:65,,,1,] * lowsamp$WAA_srv[,1:65,,,1,], 2:5, sum)
true_srvaa_mean <- apply(true_srvaa, c(2,3), mean)
true_srvaa_mean <- true_srvaa_mean / sum(true_srvaa_mean)
colors <- c('green', 'blue', 'purple', 'orange', 'brown')

plot(as.vector(true_srvaa_mean), type = 'l', lty = 2, lwd = 4, ylab = 'Aggregated Obs Survey at Age Biom', col = 'black')
for (model in 1:n_models) {
  # Average across all sims
  srvaa_all_sims <- lapply(1:n_sims, function(sim) {
    get_norm_srvaa(model_list_lowsamp_reversed[[model]][[sim]])
  })

  # Take mean across sims
  srvaa_mean <- Reduce("+", srvaa_all_sims) / n_sims
  lines(srvaa_mean, lwd = 4, type = 'l', col = colors[model])
}

# Add legend
legend('topright',
       legend = c('total at age', em_names),
       col = c('black', c('green', 'blue', 'purple', 'orange', 'brown')),
       lty = c(1, c(2:4)),
       lwd = 4,
       cex = 0.85)

# Unweighted tends to underpredict recruits more, and weighted tends to
# overpredict other ages a bit.

# From CAA, it seems like weighted approaches do not remove as many intermediate aged
# fish, which results in more numbers for SrvIAA. For unweighted approach, seems like
# selectivity removes more fish for the fishery.

# From comparing what the average observed fishery age comps tell us about removals-at-age,
# it seems like it indicates that the unweighted approach removes relatively more intermediate
# aged fish. The weighted appraoches seem to remove less intermedaite aged fish, likely
# explaining the positive bias

# Likwise, in the observed survey age compositions, it seems the intermediate
# age groups are always underepresented, while for the weighted, the
# intermediate age groups are a bit more overepresented (leading to positive bias)
