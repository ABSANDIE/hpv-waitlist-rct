# =============================================================================
# app_helpers.R
# Helpers for the results explorer (app/app.R). They reuse the model functions
# of R/hpv_functions.R and R/hpv_analysis_functions.R, and mirror the code of
# analysis/02_prespecified_analysis.Rmd so that the app gives the same numbers
# as the report. The report remains the reference: if the two ever disagree,
# the report is right and this file needs fixing.
# =============================================================================

country_label <- c(Overall = "Both countries", Mali = "Mali", CIV = "C\u00f4te d'Ivoire")
pop_label     <- c(itt = "All girls (mITT)", randomized = "Randomized girls")
pal <- c(Control = "#6b6a66", Intervention = "#2a78d6")

# ---- Data -------------------------------------------------------------------
# Same preparation as the "data" chunk of 02: factors, then the age spline with
# knots fixed once on the whole mITT population (whatever the filters later).
load_analysis_data <- function(path) {
  an <- readRDS(path) %>%
    mutate(country = factor(country, levels = c("Mali", "CIV")),
           arm = factor(arm, levels = c("Control", "Intervention")),
           chw_zone = factor(chw_zone), site = factor(site))
  knots <- rcs_knots(an$months_remaining[an$pop_itt])
  list(an = add_age_spline(an, knots), knots = knots)
}

select_pop <- function(an, population, stratum) {
  d <- if (population == "randomized") filter(an, pop_randomized) else filter(an, pop_itt)
  if (stratum != "Overall") d <- d %>% filter(country == stratum) %>% droplevels()
  d
}

# ---- Outcomes ---------------------------------------------------------------
outcomes <- list(
  cov_p1 = list(label = "Coverage at 12 weeks (1 July 2026): primary outcome",
                var = "HPV_vaccinated", uptake = FALSE, phase = 1),
  upt_p1 = list(label = "Uptake at 12 weeks among girls unvaccinated at launch",
                var = "HPV_vaccinated", uptake = TRUE, phase = 1),
  cov_p2 = list(label = "Coverage on 1 October 2026 (phase 2, exploratory)",
                var = "HPV_vaccinated_phase2", uptake = FALSE, phase = 2),
  upt_p2 = list(label = "Uptake by 1 October among girls unvaccinated at launch (phase 2)",
                var = "HPV_vaccinated_phase2", uptake = TRUE, phase = 2))
outcome_choices <- setNames(names(outcomes), vapply(outcomes, `[[`, "", "label"))
outcome_data <- function(d, key) if (outcomes[[key]]$uptake) d %>% filter(prior_vacc %in% 0L) else d

# ---- Participant flow (chunk "consort") --------------------------------------
flow_counts <- function(an) {
  an %>% group_by(country, arm) %>%
    summarise(Allocated = n(),
              `Allocated by the age-out rule` = sum(ageout),
              `Excluded at the eligibility review` = sum(eligible_final == 0),
              `Analysed (mITT)` = sum(pop_itt),
              `Randomized girls analysed` = sum(pop_randomized),
              `Phase 1 status known` = sum(pop_itt & status_known_phase1 == 1),
              `Unknown: 1 deactivated, deceased or moved` = sum(pop_itt & unknown_reason %in% "1"),
              `Unknown: 2 no information source` = sum(pop_itt & unknown_reason %in% "2"),
              `Unknown: 3a aged out` = sum(pop_itt & unknown_reason %in% "3a"),
              `Unknown: 3c not reached in time` = sum(pop_itt & unknown_reason %in% "3c"),
              `Unknown: 4 parents refused (v2.5 scheme only)` = sum(pop_itt & unknown_reason %in% "4"),
              `Phase 2 status known` = sum(pop_itt & status_known_phase2 == 1), .groups = "drop") %>%
    pivot_longer(-c(country, arm), names_to = "Step") %>%
    mutate(col = paste(country_label[as.character(country)], arm, sep = ", ")) %>%
    select(Step, col, value) %>%
    pivot_wider(names_from = col, values_from = value) %>%
    mutate(Step = factor(Step, levels = unique(Step)))
}
unknown_causes <- function(d, phase = 1) {
  v <- if (phase == 1) "unknown_reason" else "unknown_reason_phase2"
  lab <- c("1" = "1 Deactivated, deceased or moved", "2" = "2 No information source",
           "3a" = "3a Aged out", "3c" = "3c Not reached in time", "4" = "4 Parents refused")
  d %>% group_by(country, arm) %>%
    mutate(N = n()) %>% filter(!is.na(.data[[v]])) %>%
    count(country, arm, N, cause = lab[.data[[v]]]) %>% ungroup() %>%
    mutate(pct = n / N)
}

# ---- Table 1 (chunk "table1") ------------------------------------------------
table1 <- function(d) {
  g <- d$arm01
  bin <- function(x, lab) tibble(Characteristic = lab,
    Intervention = sprintf("%s (%.1f%%)", format(sum(x[g == 1], na.rm = TRUE), big.mark = ","), 100 * mean(x[g == 1], na.rm = TRUE)),
    Control = sprintf("%s (%.1f%%)", format(sum(x[g == 0], na.rm = TRUE), big.mark = ","), 100 * mean(x[g == 0], na.rm = TRUE)),
    SMD = smd(as.numeric(x), g))
  cont <- function(x, lab) {
    f <- function(v) sprintf("%.1f (%.1f); %.1f [%.1f\u2013%.1f]", mean(v, na.rm = TRUE), sd(v, na.rm = TRUE),
                             median(v, na.rm = TRUE), quantile(v, .25, na.rm = TRUE), quantile(v, .75, na.rm = TRUE))
    tibble(Characteristic = lab, Intervention = f(x[g == 1]), Control = f(x[g == 0]), SMD = smd(x, g))
  }
  bind_rows(
    tibble(Characteristic = "N", Intervention = format(sum(g == 1), big.mark = ","),
           Control = format(sum(g == 0), big.mark = ","), SMD = NA),
    cont(d$age_years, "Age at launch, years: mean (SD); median [IQR]"),
    cont(d$months_remaining, "Months remaining in eligibility window: mean (SD); median [IQR]"),
    bin(d$ageout, "Allocated by the age-out rule"),
    bin(d$prior_vacc == 1, "HPV-vaccinated before launch"),
    bin(is.na(d$prior_vacc), "Vaccination status at launch unknown"),
    bin(d$dob_approx, "Date of birth approximated"),
    map_dfr(levels(droplevels(d$site)), ~ bin(d$site == .x, paste("Site:", .x)))) %>%
    mutate(SMD = round(SMD, 2))
}

# ---- Crude coverage by country (Table 2 rows) --------------------------------
crude_by_country <- function(d, var) {
  cs <- intersect(c("Mali", "CIV"), unique(as.character(d$country)))
  if (length(cs) > 1) cs <- c("Overall", cs)
  map_dfr(cs, function(s) {
    dd <- if (s == "Overall") d else filter(d, country == s)
    u <- unadjusted(dd, var)
    if (!nrow(u)) return(tibble())
    u %>% mutate(stratum = s, .before = 1)
  })
}

# ---- Adjusted model (Table 3 row) ---------------------------------------------
adjusted_row <- function(d, var) {
  r <- fit_glmer(d, var)
  list(summary = r$summary, diag = tryCatch(diagnostics(r$fit), error = function(e) tibble()), fit = r$fit)
}

# ---- Time to vaccination (chunks "tte" and "phase2_tte") ---------------------
tte_analysis <- function(d, phase, stratum) {
  d <- d %>% filter(prior_vacc %in% 0L, !is.na(age_s1))
  if (phase == 1) {
    tau <- if (stratum == "Overall") 84 else unique(d$tau_p1)[1]
    s <- surv_by_arm(d, tau); cx <- fit_cox(d)
    per <- NULL
  } else {
    tau <- if (stratum == "Overall") 176 else unique(d$tau_p2)[1]
    s <- surv_by_arm(d, tau, "time_p2", "event_p2", weeks = c(12, 18, 25))
    cx <- fit_cox(d, "time_p2", "event_p2")
    sp <- split_periods(d, as.list(d$tau_p1), "time_p2", "event_p2")
    per <- fit_cox_period(sp, c("Before 2 July", "After 2 July"))
  }
  list(surv = s, cox = cx, period = per, n = nrow(d),
       events = sum(d[[if (phase == 1) "event_p1" else "event_p2"]]), tau = tau, data = d)
}
km_curves <- function(d, phase) {
  d <- d %>% filter(prior_vacc %in% 0L)
  tm <- if (phase == 1) "time_p1" else "time_p2"; ev <- if (phase == 1) "event_p1" else "event_p2"
  map_dfr(intersect(c("Mali", "CIV"), unique(as.character(d$country))), function(cn) {
    dd <- filter(d, country == cn)
    km_df(survfit(as.formula(sprintf("Surv(%s, %s) ~ arm", tm, ev)), data = dd)) %>%
      mutate(country = factor(country_label[cn], levels = country_label))
  })
}
plot_km <- function(df, phase) {
  tmax <- if (phase == 1) 84 else 176
  p <- ggplot(df, aes(time, cum, colour = arm)) +
    geom_step(linewidth = .8) + facet_wrap(~ country) +
    scale_colour_manual(values = pal) +
    scale_y_continuous(labels = scales::percent) + coord_cartesian(xlim = c(0, tmax)) +
    labs(x = "Days since launch", y = "Cumulative proportion vaccinated", colour = NULL)
  if (phase == 2) p <- p + geom_vline(xintercept = 84.5, linetype = 2, colour = "grey40") +
    annotate("text", x = 84.5, y = Inf, label = " control arm starts (2 July)", hjust = 0, vjust = 1.5, size = 3.5)
  p
}

# ---- Subgroups (chunks "subgroup_fun", "forest_cov", "forest_upt") -----------
run_subgroups <- function(dat, outcome = "HPV_vaccinated", progress = function(...) NULL) {
  sub_country <- function(d, s) d %>% filter(country == s) %>% droplevels()
  progress("Country subgroups")
  a <- map_dfr(intersect(c("Mali", "CIV"), unique(as.character(dat$country))),
               ~ fit_glmer(sub_country(dat, .x), outcome, label = .x)$summary %>% mutate(group = "Country"))
  progress("Age groups in C\u00f4te d'Ivoire")
  b <- if ("CIV" %in% dat$country) map_dfr(c("9-11", "12-14"), function(a) {
    rng <- if (a == "9-11") 9:11 else 12:14
    fit_glmer(dat %>% filter(country == "CIV", age_years %in% rng) %>% droplevels(), outcome,
              label = paste0("CIV, age ", a, " years"))$summary %>% mutate(group = "Age (CIV)")
  }) else tibble()
  progress("Sites")
  c <- map_dfr(levels(droplevels(dat$site)), ~ fit_glmer(dat %>% filter(site == .x) %>% droplevels(), outcome,
                                                         label = paste("Site", .x))$summary %>% mutate(group = "Site"))
  bind_rows(a, b, c) %>% mutate(group = factor(group, levels = c("Country", "Age (CIV)", "Site")))
}
flag_subgroups <- function(sg, min_events = 5) {
  sg %>% mutate(plotted = !is.na(OR) & events_int >= min_events & events_ctl >= min_events &
                  !is.na(se_logOR) & se_logOR > 0.01 & se_logOR < 3,
                why_not = case_when(plotted ~ NA_character_,
                                    is.na(OR) ~ "not estimable",
                                    events_int < min_events | events_ctl < min_events ~
                                      sprintf("fewer than %d vaccinated girls in an arm", min_events),
                                    TRUE ~ "unstable model fit"))
}
plot_subgroups <- function(sg, overall) {
  d <- sg %>% filter(plotted) %>%
    mutate(lab = sprintf("%s  (n = %s; %d/%d vaccinated)", analysis, format(n, big.mark = ","), events_int, events_ctl),
           ci_txt = fmt_ci(OR, OR_lo, OR_hi))
  if (!nrow(d)) return(NULL)
  lim <- c(min(d$OR_lo, overall$OR_lo, 1), max(d$OR_hi, overall$OR_hi, 1)) * c(0.8, 1.25)
  brks <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)
  ggplot(d, aes(OR, reorder(lab, OR))) +
    geom_vline(xintercept = 1, colour = "grey60") +
    geom_vline(xintercept = overall$OR, linetype = 2, colour = "grey40") +
    geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0, colour = pal[2], linewidth = .7) +
    geom_point(aes(size = n), colour = pal[2], show.legend = FALSE) + scale_size_area(max_size = 3.5) +
    geom_text(aes(x = lim[2], label = ci_txt), hjust = -0.15, size = 3.5, colour = "grey25") +
    scale_x_log10(breaks = brks, labels = format(brks, drop0trailing = TRUE)) +
    coord_cartesian(xlim = lim, clip = "off") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = sprintf("Adjusted OR (log scale) with 95%% CI; dashed line = overall estimate, OR %s",
                     fmt_ci(overall$OR, overall$OR_lo, overall$OR_hi)), y = NULL) +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold"),
          strip.placement = "outside", panel.spacing.y = unit(8, "pt"),
          plot.margin = margin(5.5, 130, 5.5, 5.5))
}

# ---- Coverage by CHW zone (descriptive) --------------------------------------
zone_coverage <- function(d, var) {
  d %>% filter(!is.na(.data[[var]])) %>%
    group_by(country, site, chw_zone, arm) %>%
    summarise(n = n(), vacc = sum(.data[[var]]), .groups = "drop") %>%
    mutate(cov = vacc / n,
           lo = map2_dbl(vacc, n, ~ wilson_ci(.x, .y)[["lo"]]),
           hi = map2_dbl(vacc, n, ~ wilson_ci(.x, .y)[["hi"]]))
}

# ---- Report tables written by 02 --------------------------------------------
report_tables <- c(
  "Table 1. Baseline characteristics" = "table1_baseline.csv",
  "Table 2. Coverage by arm" = "table2_coverage.csv",
  "Table 3. Primary and supplementary analyses" = "table3_primary.csv",
  "Table 4. Time to vaccination" = "table4_secondary_tte.csv",
  "Table 6. Unknown status by cause" = "table6_unknown_by_cause.csv",
  "Table 7. Sensitivity analyses (summary)" = "table7_sensitivity_summary.csv",
  "Multiple imputation" = "sens_multiple_imputation.csv",
  "Manski bounds" = "sens_manski_bounds.csv",
  "Kaplan-Meier coverage at 12 weeks" = "sens_km_coverage_12wk.csv",
  "Figure 5a data. Subgroups, coverage" = "fig5a_subgroups_coverage.csv",
  "Figure 5b data. Subgroups, uptake" = "fig5b_subgroups_uptake.csv",
  "Table 8a. Phase 2 coverage" = "table8a_phase2_coverage.csv",
  "Table 8b. Phase 2 Kaplan-Meier coverage" = "table8b_phase2_km_coverage.csv",
  "Table 8c. Phase 2 GEE" = "table8c_phase2_gee.csv",
  "Table 8d. Phase 2 time to vaccination" = "table8d_phase2_tte.csv",
  "Table 8e. Historical comparison" = "table8e_historical.csv",
  "Figure 3 data. CONSORT counts" = "fig3_consort_counts.csv")

plot_sensitivity <- function(t7) {
  d <- t7 %>% filter(!is.na(OR)) %>% mutate(scenario = factor(scenario, levels = rev(unique(scenario))),
                                            primary = row_number() == 1)
  ggplot(d, aes(OR, scenario)) +
    geom_vline(xintercept = 1, colour = "grey60") +
    geom_vline(xintercept = d$OR[1], linetype = 2, colour = "grey40") +
    geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi), height = 0, colour = pal[2], linewidth = .7) +
    geom_point(aes(shape = primary), colour = pal[2], size = 2.4, show.legend = FALSE) +
    scale_shape_manual(values = c(`TRUE` = 15, `FALSE` = 16)) +
    scale_x_log10(breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)) +
    labs(x = "Adjusted OR (log scale) with 95% CI; dashed line = primary estimate; square = primary", y = NULL)
}

# Format helpers for display
pct1 <- function(x) sprintf("%.1f%%", 100 * x)
