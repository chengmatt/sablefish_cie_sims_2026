# R/scratch

Active exploratory scripts. Nothing here is part of the main cross test or MSE workflows, and nothing sources these.

- `animation_testing.R`: gganimate experiments animating EM vs OM trajectories over assessment years. Reads a `three_region_test.RDS` that no longer exists in outputs, so it needs repointing before rerunning.
- `faa_selex_curves.R`: quick look at what the FAA fleet selectivity curves look like.
- `srvdesign_testing_plots.R`: plotting for the survey-design-change results (reads the `srvchng` cross test fits and five-region MSE output).
- `test_survey_mse.R`: sandbox for conditioning a five-region OM and testing survey MSE machinery.
- `three_rg_em_mse_base_testing.R`: interactive testing version of the three-region MSE base run.
- `validate_sglrg_inputs.R`: digs into sources of bias when aggregating spatial OM output to a single region.
