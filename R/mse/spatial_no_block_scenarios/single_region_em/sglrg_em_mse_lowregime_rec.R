# Purpose: To run a five region OM, with a single region EM (low egime recruitment)
# Creator: Matthew LH. Cheng (UAF - CFOS)
# Date: 11/26/25

# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(tidyverse)
library(furrr)
library(progressr)

# Read in model output
om_values <- readRDS(here("data", "spatial_outputs", "Spatial_MltRel_NoBlock_model_results.RDS"))

source(here("R", "functions", "mse_functions.R"))
source(here("R", "functions", "single_region_em.R"))
