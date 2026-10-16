# =============================================================================
# 00_simulate_test_data.R
# Generates SYNTHETIC CHT extracts (patients, baseline, follow-up, verification)
# for Mali and Cote d'Ivoire, with the same field names as the real tables.
# Purpose: test 01_prepare_analysis_dataset.Rmd and 02_prespecified_analysis.Rmd
# without touching real outcome data. Nothing here is real.
#
# Usage:  Rscript 00_simulate_test_data.R <out_dir> [scale]
#   scale = fraction of the real sample size (default 0.15)
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(purrr); library(readr) })
args  <- commandArgs(trailingOnly = TRUE)
out   <- if (length(args) >= 1) args[1] else "synthetic_data"
scale <- if (length(args) >= 2) as.numeric(args[2]) else 0.15
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(20260408)

d <- function(x) format(as.Date(x), "%Y-%m-%d")
uuid <- function(n, prefix) sprintf("%s-%06d", prefix, sample.int(1e6, n))

sim_country <- function(country, launch, lower, upper, n_rand_i, n_rand_c, n_ageout,
                        p_prior, rate_c, hr, n_zone, n_site) {
  phase1_end <- as.Date("2026-07-01"); cross <- as.Date("2026-07-02"); wend <- as.Date("2026-10-01")
  zones <- sprintf("%s-zone-%02d", country, seq_len(n_zone))
  sites <- sprintf("%s-site-%d", country, seq_len(n_site))
  zone_site <- setNames(sample(sites, n_zone, TRUE), zones)
  zone_eff  <- setNames(rnorm(n_zone, 0, 0.35), zones)

  n <- n_rand_i + n_rand_c + n_ageout
  arm <- c(rep("intervention", n_rand_i), rep("control", n_rand_c), rep("intervention", n_ageout))
  ageout_alloc <- c(rep(FALSE, n_rand_i + n_rand_c), rep(TRUE, n_ageout))
  # days from launch to the upper age limit (age-out date)
  days_left <- ifelse(ageout_alloc, runif(n, 1, 84),
                      runif(n, 86, 365 * (upper - lower) + 170))
  ageout_date <- launch + round(days_left)
  dob <- as.Date(sapply(ageout_date, function(x) as.character(seq(x, by = sprintf("-%d years", upper), length.out = 2)[2])))
  sex <- rep("female", n)
  # a few records whose current DOB or sex is ineligible (eligibility review)
  bad <- sample(which(!ageout_alloc), max(3, round(0.004 * n)))
  sex[bad[seq_along(bad) %% 2 == 0]] <- "male"
  dob[bad[seq_along(bad) %% 2 == 1]] <- dob[bad[seq_along(bad) %% 2 == 1]] - 365 * 3
  zone <- sample(zones, n, TRUE)
  id   <- sprintf("%s-girl-%05d", country, seq_len(n))

  # truth: prior vaccination, then post-launch vaccination time (days since launch)
  prior <- runif(n) < p_prior
  prior_date <- launch - round(runif(n, 30, 900))
  lp <- log(rate_c) + log(hr) * (arm == "intervention") + zone_eff[zone] + 0.3 * (days_left < 200)
  t1 <- rexp(n, exp(lp))                                   # before crossover
  t2 <- 84 + rexp(n, exp(log(rate_c) + log(hr) + zone_eff[zone]))  # control after crossover
  tvac <- ifelse(arm == "intervention", t1, ifelse(t1 <= 84, t1, t2))
  vac_date <- launch + ceiling(tvac)
  vac_date[vac_date >= ageout_date] <- NA                  # cannot be vaccinated once aged out
  vac_date[prior] <- prior_date[prior]
  vac_date[!is.na(vac_date) & vac_date > wend + 200] <- NA
  dead  <- runif(n) < 0.003
  muted <- runif(n) < 0.02
  mute_date <- launch + round(runif(n, 10, 170))

  patients <- tibble(
    uuid = id, name = "SYNTHETIC", sexe = ifelse(sex == "female", "Féminin", "Masculin"),
    etl_sex = ifelse(sex == "female", "female", "male"),
    date_of_birth = ifelse(runif(n) < 0.97, d(dob), NA),
    etl_dob_combine = d(dob), etl_dob_approx = ifelse(runif(n) < 0.05, d(dob), NA),
    date_of_death = ifelse(dead, paste0(d(mute_date), " 00:00:00 UTC"), NA),
    muted = ifelse(muted | dead, paste0(d(mute_date), "T10:00:00"), NA),
    c50_family_muted = NA, supervisor_area_uuid = zone, site_uuid = zone_site[zone],
    site_name = zone_site[zone], in_vaccination_trial = "true",
    vaccination_trial_branch = arm)

  status_by <- function(i, date) !is.na(vac_date[i]) && vac_date[i] <= date
  active    <- function(i, date) !((muted[i] | dead[i]) && date >= mute_date[i]) && date < ageout_date[i] + 20

  # baseline: ~80% of girls, around launch
  bl <- map_dfr(sample(seq_len(n), round(0.8 * n)), function(i) {
    vd <- launch + sample(-5:25, 1)
    any_v <- if (status_by(i, vd)) TRUE else runif(1) < 0.7
    card  <- runif(1) < 0.65
    hpv   <- if (!any_v) NA else if (!card) NA else if (status_by(i, vd)) "yes" else "no"
    tibble(uuid = uuid(1, "bl"), reported = paste0(d(vd), " 09:00:00 UTC"), c50_patient_uuid = id[i],
           c30_supervisor_area_uuid = zone[i], patient_age_in_days = as.numeric(vd - dob[i]),
           patient_sex_en = ifelse(sex[i] == "female", "Female", "Male"),
           hpv_study_branch_label_en = ifelse(arm[i] == "intervention", "Intervention", "Control"),
           s_assess_today = "yes", s_assess_date = NA, visited_date = d(vd),
           s_child_received_vaccination = ifelse(any_v, "yes", "no"),
           s_child_has_vaccination_card = ifelse(any_v, ifelse(card, "yes", "no"), NA),
           s_vaccine_HPV_1 = hpv,
           s_vaccination_date_HPV_1 = ifelse(identical(hpv, "yes") && runif(1) < 0.95, d(vac_date[i]), NA))
  })
  bl <- bind_rows(bl, bl[sample(nrow(bl), round(0.015 * nrow(bl))), ] %>% mutate(uuid = uuid(n(), "bldup")))

  # follow-up visits: more frequent in the intervention arm before crossover
  fu <- map_dfr(seq_len(n), function(i) {
    k1 <- if (arm[i] == "intervention") rpois(1, 2.2) else rpois(1, 0.3)
    k2 <- rpois(1, 1.1)
    dates <- sort(c(launch + sample(3:83, k1, TRUE), cross + sample(0:90, k2, TRUE)))
    dates <- dates[vapply(dates, function(x) active(i, x), logical(1))]
    if (!length(dates)) return(NULL)
    map_dfr(dates, function(vd) {
      info <- runif(1) < 0.92
      v <- status_by(i, vd)
      tibble(uuid = uuid(1, "fu"), reported = paste0(d(vd + sample(0:1, 1)), " 10:00:00 UTC"),
             c50_patient_uuid = id[i], c30_supervisor_area_uuid = zone[i],
             patient_age_in_days = as.numeric(vd - dob[i]),
             s_assess_today = "no", s_followup_date = d(vd), visited_date = d(vd),
             s_can_we_have_vaccination_infos = ifelse(info, "yes", "no"),
             s_child_had_vaccination_session = ifelse(info, ifelse(v, "yes", "no"), NA),
             s_vaccine_HPV_1 = ifelse(info, ifelse(v, "yes", "no"), NA),
             vaccine_HPV_1_done = ifelse(info & v, "true", "false"),
             s_vaccination_date_HPV_1 = ifelse(info & v & runif(1) < 0.93, d(vac_date[i]), NA))
    })
  })

  # verification round: girls without an in-window assessment
  seen <- unique(fu$c50_patient_uuid[as.Date(fu$s_followup_date) >= cross])
  cand <- which(!(id %in% seen))
  ver <- map_dfr(cand[runif(length(cand)) < 0.75], function(i) {
    vd <- cross + sample(5:88, 1)
    if (!active(i, vd) && !(muted[i] || dead[i])) return(NULL)
    u <- runif(1)
    moved <- (muted[i] || dead[i]) || u < 0.03
    refused <- !moved && u < 0.08
    other <- !moved && !refused && u < 0.15
    avail <- !(moved || refused || other)
    v <- status_by(i, vd)
    tibble(uuid = uuid(1, "ver"), reported = paste0(d(vd), " 11:00:00 UTC"), c50_patient_uuid = id[i],
           c30_supervisor_area_uuid = zone[i], s_reported_today = "yes", s_reported_date = NA,
           assessment_date = d(vd), vaccination_info_available = ifelse(avail, "true", "false"),
           hpv_vaccine_received = ifelse(avail, ifelse(v, "true", "false"), NA),
           hpv_vaccine_date = ifelse(avail & v & runif(1) < 0.9, d(vac_date[i]), NA),
           patient_must_be_transferred = ifelse(moved, "true", "false"),
           info_refused_by_parents = ifelse(refused, "true", "false"),
           info_unavailable_other = ifelse(other, "true", "false"))
  })
  ageout <- tibble(c50_patient_uuid = id[ageout_alloc], ao_group = "age_out")
  list(patients = patients, baseline = bl, followup = fu, verification = ver, ageout = ageout)
}

r <- function(x) max(20, round(x * scale))
sets <- list(
  Mali = sim_country("Mali", as.Date("2026-04-08"), 10, 11, r(3420), r(3458), r(1945),
                     p_prior = 0.05, rate_c = 0.0006, hr = 1.8, n_zone = 30, n_site = 4),
  CIV  = sim_country("CIV",  as.Date("2026-04-09"),  9, 15, r(4020), r(4180), r(165),
                     p_prior = 0.90, rate_c = 0.0015, hr = 1.6, n_zone = 30, n_site = 4))
for (cn in names(sets)) for (tb in names(sets[[cn]]))
  write_csv(sets[[cn]][[tb]], file.path(out, sprintf("%s_%s.csv", tolower(cn), tb)), na = "")
cat("Synthetic tables written to", normalizePath(out), "\n")
