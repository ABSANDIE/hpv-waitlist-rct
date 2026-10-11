# =============================================================================
# hpv_functions.R
# Shared helpers for the HPV community follow-up trial (Mali, Cote d'Ivoire).
# Implements the derivations of the registered SAP v2.1, mapped to the CHT
# fields as documented in SAP v2.5 ("Derivation from CHT Fields").
# Sourced by 01_prepare_analysis_dataset.Rmd and 02_prespecified_analysis.Rmd.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(purrr)
})

# ---- Study constants (SAP) ---------------------------------------------
study_dates <- list(
  launch        = c(Mali = as.Date("2026-04-08"), CIV = as.Date("2026-04-09")),
  phase1_end    = as.Date("2026-07-01"),   # primary endline: status on this date
  crossover     = as.Date("2026-07-02"),   # control arm starts the intervention
  verif_close_1 = as.Date("2026-08-30"),   # original close of verification window (sensitivity)
  window_end    = as.Date("2026-10-01")    # verification of phase 1 status closes; phase 2 endline
)
# National HPV eligibility (age in completed years): lower inclusive, upper exclusive.
# A girl "ages out" on the day she reaches the upper limit.
eligibility <- list(
  Mali = c(lower = 10, upper = 11),
  CIV  = c(lower = 9,  upper = 15)
)
# Known number of girls allocated to the intervention arm by the age-out rule
expected_ageout <- c(Mali = 1945, CIV = 165)

# ---- Parsing helpers --------------------------------------------------------
# Dates arrive as "2026-04-08", "2026-07-01T00:00:00.000-00:00" or
# "2026-04-08 08:17:14.34 UTC": keep the first 10 characters.
to_date <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA
  suppressWarnings(as.Date(substr(x, 1, 10)))
}
# CHT booleans and yes/no selects: yes/no, true/false, TRUE/FALSE, VRAI/FAUX, oui/non, 1/0
to_bool <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(x))
  out[x %in% c("yes", "true", "vrai", "oui", "1", "t")] <- TRUE
  out[x %in% c("no", "false", "faux", "non", "0", "f")] <- FALSE
  out
}
norm_sex <- function(x) {
  x <- tolower(trimws(as.character(x)))
  case_when(x %in% c("female", "femme", "f", "fille") ~ "female",
            x %in% c("male", "homme", "m", "garcon", "garçon") ~ "male",
            TRUE ~ NA_character_)
}
norm_arm <- function(x) {
  x <- tolower(trimws(as.character(x)))
  case_when(str_detect(x, "interv") ~ "Intervention",
            str_detect(x, "contr") ~ "Control",
            TRUE ~ NA_character_)
}
# Column accessor that returns NA when a column is absent (forms differ slightly
# between countries and versions).
col_or_na <- function(df, name) if (name %in% names(df)) df[[name]] else rep(NA, nrow(df))
add_years <- function(d, y) {           # calendar-exact anniversaries (29 Feb -> 1 Mar)
  as.Date(ifelse(is.na(d), NA,
          as.character(lubridate::add_with_rollback(d, lubridate::years(y), roll_to_first = TRUE))))
}

read_table <- function(path) {
  if (is.null(path) || is.na(path) || !file.exists(path)) return(NULL)
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    as.data.frame(readxl::read_excel(path, guess_max = 1e5))
  } else {
    as.data.frame(readr::read_csv(path, col_types = cols(.default = col_character()),
                                  progress = FALSE, show_col_types = FALSE))
  }
}

# ---- Assessment records: one row per form submission ------------------------
# hpv: "yes" (vaccinated), "no" (not vaccinated), NA (status not obtained)
# vacc_date: documented HPV vaccination date (NA if not documented)
# reason flags explain why status was not obtained.

records_baseline <- function(df, country) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  assess_date <- coalesce(to_date(col_or_na(df, "visited_date")),
                          if_else(to_bool(col_or_na(df, "s_assess_today")) %in% FALSE,
                                  to_date(col_or_na(df, "s_assess_date")), as.Date(NA)),
                          to_date(col_or_na(df, "reported")))
  any_vacc <- to_bool(col_or_na(df, "s_child_received_vaccination"))
  card     <- to_bool(col_or_na(df, "s_child_has_vaccination_card"))
  hpv_card <- to_bool(col_or_na(df, "s_vaccine_HPV_1"))
  # A girl who never received any vaccination is not HPV-vaccinated: the form
  # does not ask about HPV in that case, so the HPV field is blank but status is known.
  hpv <- case_when(any_vacc %in% FALSE ~ "no",
                   hpv_card %in% TRUE  ~ "yes",
                   hpv_card %in% FALSE ~ "no",
                   TRUE ~ NA_character_)
  tibble(country = country, girl_id = df$c50_patient_uuid, form = "baseline",
         record_uuid = col_or_na(df, "uuid"), reported = to_date(col_or_na(df, "reported")),
         assess_date = assess_date, hpv = hpv,
         vacc_date = if_else(hpv %in% "yes", to_date(col_or_na(df, "s_vaccination_date_HPV_1")), as.Date(NA)),
         no_card = card %in% FALSE, info_refused = FALSE, moved = FALSE, info_unavailable = FALSE,
         chw_zone_form = col_or_na(df, "c30_supervisor_area_uuid"),
         age_days_form = suppressWarnings(as.numeric(col_or_na(df, "patient_age_in_days"))))
}

records_followup <- function(df, country, no_info_values = c("no", "non", "info_not_available")) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  assess_date <- coalesce(if_else(to_bool(col_or_na(df, "s_assess_today")) %in% FALSE,
                                  to_date(col_or_na(df, "s_followup_date")), as.Date(NA)),
                          to_date(col_or_na(df, "visited_date")),
                          to_date(col_or_na(df, "reported")))
  hpv_card <- to_bool(col_or_na(df, "s_vaccine_HPV_1"))
  hpv_done <- to_bool(col_or_na(df, "vaccine_HPV_1_done"))
  session  <- to_bool(col_or_na(df, "s_child_had_vaccination_session"))
  infos    <- tolower(trimws(as.character(col_or_na(df, "s_can_we_have_vaccination_infos"))))
  no_info  <- infos %in% no_info_values
  hpv <- case_when(hpv_card %in% TRUE | hpv_done %in% TRUE ~ "yes",
                   hpv_card %in% FALSE ~ "no",
                   session %in% FALSE & !no_info ~ "no",
                   TRUE ~ NA_character_)
  tibble(country = country, girl_id = df$c50_patient_uuid, form = "followup",
         record_uuid = col_or_na(df, "uuid"), reported = to_date(col_or_na(df, "reported")),
         assess_date = assess_date, hpv = hpv,
         vacc_date = if_else(hpv %in% "yes", to_date(col_or_na(df, "s_vaccination_date_HPV_1")), as.Date(NA)),
         no_card = FALSE, info_refused = FALSE, moved = FALSE, info_unavailable = no_info,
         chw_zone_form = col_or_na(df, "c30_supervisor_area_uuid"),
         age_days_form = suppressWarnings(as.numeric(col_or_na(df, "patient_age_in_days"))))
}

records_verification <- function(df, country) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  assess_date <- coalesce(to_date(col_or_na(df, "assessment_date")),
                          if_else(to_bool(col_or_na(df, "s_reported_today")) %in% FALSE,
                                  to_date(col_or_na(df, "s_reported_date")), as.Date(NA)),
                          to_date(col_or_na(df, "reported")))
  avail    <- to_bool(coalesce(col_or_na(df, "vaccination_info_available"),
                               col_or_na(df, "s_vaccination_info_available")))
  received <- to_bool(coalesce(col_or_na(df, "hpv_vaccine_received"),
                               col_or_na(df, "s_hpv_vaccine_received")))
  hpv <- case_when(avail %in% FALSE ~ NA_character_,
                   received %in% TRUE ~ "yes",
                   received %in% FALSE ~ "no",
                   TRUE ~ NA_character_)
  tibble(country = country, girl_id = df$c50_patient_uuid, form = "verification",
         record_uuid = col_or_na(df, "uuid"), reported = to_date(col_or_na(df, "reported")),
         assess_date = assess_date, hpv = hpv,
         vacc_date = if_else(hpv %in% "yes",
                             to_date(coalesce(col_or_na(df, "hpv_vaccine_date"),
                                              col_or_na(df, "s_hpv_vaccine_date"))), as.Date(NA)),
         no_card = FALSE,
         info_refused     = to_bool(col_or_na(df, "info_refused_by_parents")) %in% TRUE,
         moved            = to_bool(col_or_na(df, "patient_must_be_transferred")) %in% TRUE,
         info_unavailable = to_bool(col_or_na(df, "info_unavailable_other")) %in% TRUE,
         chw_zone_form = col_or_na(df, "c30_supervisor_area_uuid"),
         age_days_form = NA_real_)
}

# Data-quality rules applied to the record table (counts are reported by the QC section)
clean_records <- function(rec, launch_by_country, first_vaccine_date = as.Date("2019-01-01")) {
  rec %>%
    mutate(launch = launch_by_country[country],
           qc_assess_missing   = is.na(assess_date),
           qc_assess_outside   = !is.na(assess_date) & (assess_date < launch - 30 | assess_date > study_dates$window_end),
           qc_vacc_after_visit = !is.na(vacc_date) & !is.na(assess_date) & vacc_date > assess_date,
           qc_vacc_implausible = !is.na(vacc_date) & (vacc_date < first_vaccine_date | vacc_date > study_dates$window_end),
           # A vaccination date after the visit, or outside the programme period, is not
           # usable; the dose is kept but treated as undated (vaccinated by the visit date).
           vacc_date = if_else(qc_vacc_after_visit | qc_vacc_implausible, as.Date(NA), vacc_date)) %>%
    filter(!is.na(girl_id), !qc_assess_missing, !qc_assess_outside)
}

# ---- Girl-level status at a cut-off date -------------------------------------
# Returns one row per girl:
#   vacc        : 1 vaccinated by `cutoff`, 0 known not vaccinated on `cutoff`, NA unknown
#   first_dose  : earliest documented dose date (any time)
#   undated_by  : earliest visit reporting a dose without a usable date
#   last_unvacc : last date on/before `cutoff` the girl was documented as unvaccinated
#   status_date : date of the visit that established `vacc` (assess_date_phase1 in the SAP)
#   undated_after: a dose reported without a usable date that cannot be counted:
#                  at a visit after `cutoff`, or at any visit when undated_counts = FALSE
#                  (-> cause 2 if status unknown)
# `undated_counts`: FALSE (registered SAP v2.1): a dose without a documented date
#   gives unknown status. TRUE (sensitivity): a dose reported without a date at a
#   visit on/before `cutoff` counts as vaccinated by that visit.
# `known_from`: earliest assessment date that can establish "not vaccinated on cutoff".
#   Phase 1 (SAP): assessments dated 2 July - 1 October 2026 -> known_from = 2 July.
#   Phase 2 (SAP): status at the last visit on/before 1 October -> known_from = launch.
status_at <- function(rec, cutoff, known_from, window_end = study_dates$window_end,
                      undated_counts = FALSE) {
  rec %>%
    filter(assess_date <= window_end) %>%
    mutate(kf = if (length(known_from) == 1) known_from else known_from[country],
           yes_dated_by  = hpv %in% "yes" & !is.na(vacc_date) & vacc_date <= cutoff,
           yes_undat_by  = undated_counts & hpv %in% "yes" & is.na(vacc_date) & assess_date <= cutoff,
           est_no        = (hpv %in% "no" & assess_date >= kf) |
                           (hpv %in% "yes" & !is.na(vacc_date) & vacc_date > cutoff & assess_date >= kf)) %>%
    group_by(country, girl_id) %>%
    summarise(
      first_dose  = suppressWarnings(min(vacc_date[hpv %in% "yes"], na.rm = TRUE)),
      undated_by  = suppressWarnings(min(assess_date[hpv %in% "yes" & is.na(vacc_date)], na.rm = TRUE)),
      known_no    = any(est_no),
      last_unvacc = suppressWarnings(max(assess_date[hpv %in% "no" & assess_date <= cutoff], na.rm = TRUE)),
      date_yes    = suppressWarnings(min(assess_date[yes_dated_by | yes_undat_by], na.rm = TRUE)),
      date_no     = suppressWarnings(min(assess_date[est_no], na.rm = TRUE)),
      undated_after = any(hpv %in% "yes" & is.na(vacc_date) & (assess_date > cutoff | !undated_counts)),
      conflict    = any(yes_dated_by) & any(hpv %in% "no" & assess_date > cutoff),
      n_contacts  = n(),
      last_visit  = max(assess_date),
      .groups = "drop") %>%
    mutate(across(c(first_dose, undated_by, last_unvacc, date_yes, date_no),
                  ~ if_else(is.finite(.x), .x, as.Date(NA))),
           vacc = case_when(!is.na(first_dose) & first_dose <= cutoff ~ 1L,
                            undated_counts & !is.na(undated_by) & undated_by <= cutoff ~ 1L,
                            known_no ~ 0L,
                            TRUE ~ NA_integer_),
           status_date = case_when(vacc == 1L ~ date_yes, vacc == 0L ~ date_no))
}

# ---- Restricted cubic spline basis (Harrell), default 3 knots at 10/50/90th pct
rcs_knots <- function(x, nk = 3) {
  q <- switch(as.character(nk), "3" = c(.10, .50, .90), "4" = c(.05, .35, .65, .95),
              "5" = c(.05, .275, .50, .725, .95))
  unname(quantile(x, q, na.rm = TRUE))
}
rcs_basis <- function(x, knots) {
  k <- length(knots); tk <- knots[k]; tk1 <- knots[k - 1]
  out <- sapply(seq_len(k - 2), function(j) {
    pmax(x - knots[j], 0)^3 - pmax(x - tk1, 0)^3 * (tk - knots[j]) / (tk - tk1) +
      pmax(x - tk, 0)^3 * (tk1 - knots[j]) / (tk - tk1)
  })
  out <- as.matrix(out) / (tk - knots[1])^2       # Harrell's scaling
  cbind(x, out)
}

# ---- Small statistics helpers ------------------------------------------------
wilson_ci <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2); p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(est = p, lo = centre - half, hi = centre + half)
}
# Newcombe hybrid score interval (method 10) for p1 - p0
newcombe_diff <- function(x1, n1, x0, n0, conf = 0.95) {
  a <- wilson_ci(x1, n1, conf); b <- wilson_ci(x0, n0, conf)
  d <- a[["est"]] - b[["est"]]
  c(diff = d,
    lo = d - sqrt((a[["est"]] - a[["lo"]])^2 + (b[["hi"]] - b[["est"]])^2),
    hi = d + sqrt((a[["hi"]] - a[["est"]])^2 + (b[["est"]] - b[["lo"]])^2))
}
smd <- function(x, g) {          # standardized mean difference (continuous or binary)
  x1 <- x[g == 1]; x0 <- x[g == 0]
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) /
    sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
}
fmt_ci <- function(est, lo, hi, digits = 2, pct = FALSE) {
  k <- if (pct) 100 else 1
  sprintf(paste0("%.", digits, "f (%.", digits, "f to %.", digits, "f)"), k * est, k * lo, k * hi)
}
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
