# Copy to config_real.R (ignored by git) and edit. Paths are relative to analysis/.
params_prep <- list(
  data_dir          = "../data/raw",
  mali_patients     = "patients_mali.csv",
  mali_baseline     = "children_vaccination_mali.csv",
  mali_followup     = "children_vaccination_followup_mali.csv",
  mali_verification = "hpv_verification_mali.csv",
  civ_patients      = "patients_rci.csv",
  civ_baseline      = "children_vaccination_rci.csv",
  civ_followup      = "children_vaccination_followup_r.csv",
  civ_verification  = "hpv_verification_rci.csv",
  out_dir           = "../outputs/real")
params_analysis <- list(data = "../outputs/real/analysis_dataset.rds", out_dir = "../outputs/real")
