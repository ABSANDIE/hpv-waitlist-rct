# =============================================================================
# run_all.R - reproduce every output from scratch.
#   Rscript run_all.R          full run on synthetic data (50 imputations, ~20-40 min)
#   Rscript run_all.R quick    quick check (5 imputations, 50 bootstrap replicates)
#   Rscript run_all.R real     real CHT extracts, using config_real.R (never committed)
# Run from the repository root.
# =============================================================================
args  <- commandArgs(trailingOnly = TRUE)
mode  <- if (length(args)) args[1] else "full"
stopifnot(mode %in% c("full", "quick", "real"))
params_prep <- list(); params_analysis <- list()
if (mode == "real") {
  if (!file.exists("config_real.R")) stop("Copy config_real_template.R to config_real.R and edit the paths.")
  source("config_real.R")            # defines params_prep (and optionally params_analysis)
} else {
  if (!file.exists("data/synthetic/mali_patients.csv")) {
    message("Generating synthetic data ...")
    system2("Rscript", c("R/00_simulate_test_data.R", "data/synthetic", "0.15"))
  }
  if (mode == "quick") params_analysis <- list(mi_m = 5, boot_B = 50)
}
# Real-data reports go to outputs/real/, next to the real-data tables, so they never
# overwrite (or get mistaken for) the synthetic reports in outputs/
out <- if (mode == "real") "outputs/real" else "outputs"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
rmarkdown::render("analysis/01_prepare_analysis_dataset.Rmd", params = params_prep, output_dir = out,
                  envir = new.env(), quiet = TRUE)
rmarkdown::render("analysis/02_prespecified_analysis.Rmd",   params = params_analysis, output_dir = out,
                  envir = new.env(), quiet = TRUE)
message("Done. Reports and tables are in ", out, "/")
