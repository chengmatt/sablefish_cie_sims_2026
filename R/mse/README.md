# R/mse

Closed-loop simulations: the five-region OM runs forward with annual feedback (assess with the EM, project, set ABC via the threshold HCR, apportion, harvest), so EMs are compared as part of a full management procedure rather than as one-off fits.

Layout within `spatial_no_block_scenarios/`: one folder per EM, and inside each, one script per recruitment scenario:

- `*_base_rec.R`: base (mean) recruitment
- `*_highregime_rec.R`: high recruitment regime
- `*_lowregime_rec.R`: low recruitment regime

The EM folders are `single_region_em/`, `three_rg_em/`, `five_rg_em/`, and `faa_em/`. Note the Purpose headers inside several scripts say "single region EM" from copy-paste; the folder name and the `run_*_closedloop_parallel()` call are what the script actually runs. All scripts use the same five-region OM; only the EM differs.

The FAA scripts each run twice: once with the standard survey selectivity and once with a time-blocked survey selectivity EM variant, saved with a `_tb` tag.

Results go to `outputs/mse_results/spatial_no_block_scenarios/` (gitignored) as `<em>_<scenario>_<i>.RDS`, where `i` is a batch index over chunks of simulations (batch counts differ by EM; the five-region EM saves a single file per scenario).

`mse_comparisons_base.R`, `mse_comparisons_highregimerec.R`, and `mse_comparisons_lowregimerec.R` read those batches and write the comparison figures to `outputs/mse_results/figs/`, one script per recruitment scenario, with the scenario as the figure suffix.
