# data

Inputs to the simulations. Nothing in this project writes here.

`spatial_outputs/` holds the conditioned operating model, i.e. fitted results of the multi-region spatial sablefish assessment (produced outside this repo). Each RDS is a list with `data`, `rep`, `parameters`, `mapping`, and `sd_rep`, which the OM conditioning scripts unpack directly:

- `Spatial_MltRel_NoBlock_model_results.RDS`: spatial model WITHOUT survey selectivity blocks. This conditions everything currently active (`R/cross_test` and `R/mse` under `spatial_no_block_scenarios`).
- `Spatial_MltRel_model_results.RDS`: spatial model WITH selectivity blocks. Conditioned the blocked-era work now archived in `R/old/`; still read by a couple of `R/scratch` scripts.

`.Rapp.history` is a stray R console history file, safe to delete.
