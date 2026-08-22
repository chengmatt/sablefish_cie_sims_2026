# R/old

Superseded scripts from the blocked era, when the OM was conditioned on the spatial sablefish model WITH selectivity blocks. The current workflows in `R/cross_test` and `R/mse` replaced these after the OM switched to the no-selectivity-block conditioning (`Spatial_MltRel_NoBlock_model_results.RDS`).

- `spatial_blocked_scenarios/single_region_em/`: single-region EM closed loops (high regime, low regime, random recruitment) with survey feedback.
- `spatial_blocked_selex_scenarios/`: the blocked-era OM simulation generator and single-region/FAA cross tests.

Their results live in `outputs/cross_test/old/` and `outputs/mse_results/old/`.

Frozen: paths inside these scripts still point at pre-archive locations and are not kept up to date. To rerun one, update its input/output paths first.
