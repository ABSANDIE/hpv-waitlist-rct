# =============================================================================
# export_public_summaries.R
# Builds app_public/data/summaries.rds, the only data file of the public demo
# app, from the SYNTHETIC outputs (outputs/analysis_dataset.rds). It contains
# aggregated results only (counts, percentages, model estimates, survival-curve
# points). The public app is a demonstration of the analysis on invented data:
# this script refuses to run on real trial data.
#
# Run from the repository root, after Rscript run_all.R (or run_all.R quick):
#   Rscript app_public/export_public_summaries.R
# Takes a few minutes (about 40 model fits).
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(survival); library(ggplot2)
})
src <- "outputs"                       # synthetic outputs only; never outputs/real
min_cell <- 10
source("R/hpv_functions.R"); source("R/hpv_analysis_functions.R"); source("app/app_helpers.R")
stopifnot(file.exists(file.path(src, "analysis_dataset.rds")))
message("Reading ", file.path(src, "analysis_dataset.rds"))
D <- load_analysis_data(file.path(src, "analysis_dataset.rds")); an <- D$an
# Guard: the synthetic generator names girls "Mali-girl-00001", "CIV-girl-00001", ...
if (!all(grepl("^(Mali|CIV)-girl-[0-9]+$", an$girl_id)))
  stop("outputs/analysis_dataset.rds is not the synthetic dataset. The public app is for synthetic data only: ",
       "run `Rscript run_all.R quick` to rebuild the synthetic outputs, then try again.")

pops <- names(pop_label); strata <- names(country_label)
grid <- function(f, a, b) setNames(map(a, function(x) setNames(map(b, function(y) f(x, y)), b)), a)
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

message("Flow, baseline and crude coverage")
S <- list()
S$meta <- list(synthetic = TRUE, created = Sys.Date(), knots = D$knots,
               n_mitt = sum(an$pop_itt), min_cell = min_cell)
S$flow <- flow_counts(an) %>% mutate(Step = as.character(Step))
S$causes <- grid(function(p, s) {
  d <- select_pop(an, p, s)
  list(`1` = unknown_causes(d, 1), `2` = unknown_causes(d, 2))
}, pops, strata)
S$table1 <- grid(function(p, s) table1(select_pop(an, p, s)), pops, strata)
S$crude <- grid(function(p, s) {
  d <- select_pop(an, p, s)
  map(setNames(names(outcomes), names(outcomes)), ~ crude_by_country(outcome_data(d, .x), outcomes[[.x]]$var))
}, pops, strata)

message("Adjusted models (24 fits)")
S$adj <- grid(function(p, s) {
  d <- select_pop(an, p, s)
  map(setNames(names(outcomes), names(outcomes)), function(o) {
    message("  ", p, " / ", s, " / ", o)
    r <- quiet(adjusted_row(outcome_data(d, o), outcomes[[o]]$var))
    list(summary = r$summary, diag = r$diag)              # the fitted model object is not kept
  })
}, pops, strata)

message("Time to vaccination")
S$tte <- grid(function(p, s) {
  d <- select_pop(an, p, s)
  map(c(`1` = 1, `2` = 2), function(ph) {
    r <- quiet(tte_analysis(d, ph, s))
    list(rmst = r$surv$rmst, cum = r$surv$cum, dRMST = r$surv$dRMST, logrank_p = r$surv$logrank_p,
         cox = r$cox$summary, period = if (!is.null(r$period)) list(summary = r$period$summary,
                                                                   interaction_p = r$period$interaction_p),
         n = r$n, events = r$events, tau = r$tau)
  })
}, pops, strata)
# Kaplan-Meier curves: one point per day with an event (cumulative proportion only)
S$km <- grid(function(p, s) {
  d <- select_pop(an, p, s)
  map(c(`1` = 1, `2` = 2), ~ km_curves(d, .x) %>% select(country, arm, time, cum, lo, hi))
}, pops, strata)

message("Subgroups")
S$sg <- setNames(map(pops, function(p) {
  map(c(cov_p1 = "cov_p1", upt_p1 = "upt_p1"), function(o) {
    d <- select_pop(an, p, "Overall") %>% outcome_data(o)
    list(sg = quiet(run_subgroups(d) %>% flag_subgroups()), overall = quiet(fit_glmer(d, "HPV_vaccinated")$summary))
  })
}), pops)

message("Coverage by CHW zone (anonymized, small cells suppressed)")
zone_key <- an %>% distinct(country, chw_zone) %>% arrange(country, as.character(chw_zone)) %>%
  group_by(country) %>% mutate(zone = sprintf("Zone %02d", row_number())) %>% ungroup()
S$zone <- grid(function(p, s) {
  d <- select_pop(an, p, s)
  map(setNames(names(outcomes), names(outcomes)), function(o)
    zone_coverage(outcome_data(d, o), outcomes[[o]]$var) %>%
      left_join(zone_key, by = c("country", "chw_zone")) %>%
      select(country, zone, arm, n, vacc, cov, lo, hi) %>%
      filter(n >= min_cell))
}, pops, strata)

message("Report tables")
S$report <- compact(map(report_tables, function(f) {
  p <- file.path(src, f); if (file.exists(p)) read_csv(p, show_col_types = FALSE) else NULL
}))

# ---- Safety check: nothing girl-level ---------------------------------------
walk_tbls <- function(x) if (is.data.frame(x)) list(x) else if (is.list(x)) flatten(map(x, walk_tbls)) else list()
tbls <- walk_tbls(S)
bad <- c("girl_id", "uuid", "patient_id", "c50_patient_uuid", "chw_zone", "date_of_birth", "dob")
hit <- keep(tbls, ~ any(tolower(names(.x)) %in% bad))
if (length(hit)) stop("A summary table contains an identifier column: ", paste(unique(unlist(map(hit, names))), collapse = ", "))
if (max(map_int(tbls, nrow)) > 0.2 * nrow(an)) stop("A summary table is too large to be an aggregate: check before publishing.")
message(sprintf("%d summary tables; largest has %d rows (the dataset has %d girls).",
                length(tbls), max(map_int(tbls, nrow)), nrow(an)))

dir.create("app_public/data", showWarnings = FALSE, recursive = TRUE)
saveRDS(S, "app_public/data/summaries.rds")
message("Written: app_public/data/summaries.rds (synthetic data)")
