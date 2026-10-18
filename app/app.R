# =============================================================================
# HPV community follow-up trial: results explorer (Shiny)
#
# Reads the girl-level dataset made by analysis/01_prepare_analysis_dataset.Rmd
# (analysis_dataset.rds) and the tables written by 02_prespecified_analysis.Rmd,
# and lets you look at the pre-specified analyses by population, country and
# outcome. Models are refitted with the same functions as the report.
#
# Run from the repository root, after run_all.R:
#   shiny::runApp("app")
# For local use only: with real data the app shows girl-level derived results.
# Do not deploy it to a public server.
# =============================================================================
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(DT); library(ggplot2)
  library(dplyr); library(tidyr); library(purrr); library(readr); library(survival)
})
# shiny::runApp("app") sets the working directory to app/; also allow running from the root
root <- if (file.exists("../R/hpv_functions.R")) ".." else "."
source(file.path(root, "R/hpv_functions.R"))
source(file.path(root, "R/hpv_analysis_functions.R"))
source(file.path(root, "app/app_helpers.R"))

theme_set(theme_minimal(base_size = 13) +
            theme(panel.grid.minor = element_blank(), legend.position = "top",
                  strip.text = element_text(face = "bold", size = 12)))

# ---- Data sources: every outputs folder that contains analysis_dataset.rds ---
find_sources <- function() {
  extra <- Sys.getenv("HPV_OUTPUTS")                    # optional: another folder
  cand <- c(file.path(root, "outputs/real"), file.path(root, "outputs"), if (nzchar(extra)) extra)
  cand <- cand[file.exists(file.path(cand, "analysis_dataset.rds"))]
  lab <- ifelse(grepl("outputs/real$", cand), "Real trial data (outputs/real)",
         ifelse(grepl("outputs$", cand), "Synthetic data (outputs)", cand))
  setNames(normalizePath(cand), lab)
}
sources <- find_sources()
if (!length(sources)) stop("No analysis_dataset.rds found. Run run_all.R (or run_all.R real) first.")

# Simple cache shared by all sessions: model fits are slow, inputs are few
cache <- new.env()
cached <- function(key, expr) {
  if (!exists(key, envir = cache, inherits = FALSE)) assign(key, force(expr), envir = cache)
  get(key, envir = cache)
}

dt <- function(x, digits = NULL, page = 25, ...) {
  out <- datatable(x, rownames = FALSE, options = list(pageLength = page, scrollX = TRUE, dom = "tip"), ...)
  num <- names(x)[vapply(x, is.double, logical(1))]
  if (length(num) && !is.null(digits)) out <- formatRound(out, num, digits)
  out
}
note <- function(...) div(class = "text-muted small mt-2", ...)
comma <- function(x) formatC(x, format = "d", big.mark = ",")
med <- function(x) ifelse(is.na(x), "not reached", format(x))
country_f <- function(x) factor(country_label[as.character(x)], levels = country_label)
stat_tile <- function(title, est, lo, hi, sub = NULL, primary = FALSE) {
  div(class = paste("p-3 rounded flex-fill", if (primary) "text-white bg-primary" else "border"),
      style = "min-width: 170px;",
      div(class = "small", title),
      div(style = "font-size: 1.9rem; font-weight: 600; line-height: 1.2;", sprintf("%.2f", est)),
      div(class = "small", sprintf("95%% CI %.2f to %.2f", lo, hi)),
      if (!is.null(sub)) div(class = "small opacity-75", sub))
}

# =============================================================================
ui <- page_navbar(
  title = "HPV trial results explorer",
  fillable = FALSE,
  theme = bs_theme(version = 5, primary = "#2a78d6", base_font = font_collection("system-ui", "sans-serif")),
  sidebar = sidebar(
    width = 290,
    selectInput("source", "Data", choices = sources),
    uiOutput("source_badge"),
    radioButtons("pop", "Population", choices = setNames(names(pop_label), pop_label)),
    radioButtons("stratum", "Country", choices = setNames(names(country_label), country_label)),
    note("Under the SAP, the estimate in randomized girls takes precedence for causal",
         "interpretation if the two populations diverge. The age-out girls are all in the",
         "intervention arm, so crude comparisons pooling both countries are confounded."),
    note(textOutput("knots_txt"))
  ),

  nav_panel("Participant flow",
    layout_columns(col_widths = c(12, 12), fill = FALSE,
      card(card_header("CONSORT counts by country and arm"), DTOutput("flow_tbl", fill = FALSE)),
      card(card_header("Causes of unknown status (share of the analysed population)"),
           radioButtons("cause_phase", NULL, inline = TRUE,
                        choices = c("Phase 1 (1 July)" = 1, "Phase 2 (1 October)" = 2)),
           plotOutput("cause_plot", height = 420),
           note("Follows the population and country chosen on the left.")))),

  nav_panel("Baseline",
    card(card_header(textOutput("t1_title")), DTOutput("t1_tbl", fill = FALSE),
         note("SMD = standardized mean difference (intervention vs control)."))),

  nav_panel("Coverage and uptake",
    selectInput("outcome", "Outcome", choices = outcome_choices, width = "520px"),
    layout_columns(col_widths = c(5, 7), fill = FALSE,
      card(card_header("Crude coverage by arm (Wilson 95% CI)"), plotOutput("crude_plot", height = 330)),
      card(card_header("Adjusted analysis (pre-specified model)"),
           uiOutput("adj_box"),
           note("Mixed-effects logistic regression: arm + country + age spline (3 knots) +",
                "random intercept for CHW supervision zone. ARD and RR by marginal standardization."))),
    card(card_header("Crude comparison (Table 2 layout)"), DTOutput("crude_tbl", fill = FALSE))),

  nav_panel("Time to vaccination",
    radioButtons("tte_phase", NULL, inline = TRUE,
                 choices = c("12-week time frame (secondary outcome)" = 1,
                             "Full follow-up to 1 October (phase 2, exploratory)" = 2)),
    card(card_header("Cumulative proportion vaccinated, girls unvaccinated at launch (Kaplan-Meier)"),
         plotOutput("km_plot", height = 380)),
    layout_columns(col_widths = c(6, 6), fill = FALSE,
      card(card_header("RMST and hazard ratio"), uiOutput("tte_box")),
      card(card_header("Cumulative proportion vaccinated by arm"), DTOutput("tte_cum", fill = FALSE)))),

  nav_panel("Subgroups",
    radioButtons("sg_outcome", NULL, inline = TRUE,
                 choices = c("Coverage at 12 weeks (Figure 5a)" = "cov_p1",
                             "Uptake at 12 weeks, unvaccinated at launch (Figure 5b, post hoc)" = "upt_p1")),
    card(card_header("Subgroup-specific adjusted ORs"), plotOutput("sg_plot", height = 520),
         note("Exploratory; 95% CIs without multiplicity adjustment. Uses both countries of the",
              "chosen population. Subgroups with fewer than 5 vaccinated girls in an arm or an",
              "unstable fit are listed below and not plotted.")),
    card(card_header("Subgroups not plotted"), DTOutput("sg_tbl", fill = FALSE))),

  nav_panel("CHW zones",
    selectInput("zone_outcome", "Outcome", choices = outcome_choices, width = "520px"),
    card(card_header("Crude coverage by CHW supervision zone and arm"),
         plotOutput("zone_plot", height = 700),
         note("Descriptive. Zones are sorted by control-arm coverage; point size shows the number of girls",
              "with known status. The primary model accounts for zones with a random intercept."))),

  nav_panel("Sensitivity analyses",
    card(card_header("Table 7. Sensitivity analyses of the primary outcome (from the report)"),
         plotOutput("sens_plot", height = 620), DTOutput("sens_tbl", fill = FALSE),
         note("Read from table7_sensitivity_summary.csv written by 02_prespecified_analysis.Rmd.",
              "A material change is an OR differing by more than 10% or a change in significance."))),

  nav_panel("Report tables",
    selectInput("rep_tbl", "Table", choices = report_tables, width = "520px"),
    card(DTOutput("rep_view", fill = FALSE), note(textOutput("rep_path")))),

  nav_spacer(),
  nav_item(tags$a("NCT07807202", href = "https://clinicaltrials.gov/study/NCT07807202", target = "_blank"))
)

# =============================================================================
server <- function(input, output, session) {
  src <- reactive(input$source)
  key0 <- reactive(paste(src(), file.mtime(file.path(src(), "analysis_dataset.rds"))))

  dat <- reactive({
    cached(paste("data", key0()), {
      withProgress(message = "Loading the analysis dataset", load_analysis_data(file.path(src(), "analysis_dataset.rds")))
    })
  })
  pop_data <- reactive(select_pop(dat()$an, input$pop, input$stratum))
  sel_label <- reactive(paste(pop_label[input$pop], country_label[input$stratum], sep = ", "))

  output$source_badge <- renderUI({
    if (grepl("real", names(sources)[sources == src()]))
      div(class = "alert alert-warning py-2 small", strong("Real trial data."),
          " For the study team only: do not share screenshots or exports outside the team.")
    else div(class = "alert alert-info py-2 small", strong("Synthetic data."),
             " Results say nothing about the trial.")
  })
  output$knots_txt <- renderText(sprintf("Age spline knots (months remaining): %s",
                                         paste(round(dat()$knots, 2), collapse = ", ")))

  # ---- Flow -------------------------------------------------------------------
  output$flow_tbl <- renderDT({
    f <- flow_counts(dat()$an) %>% mutate(Step = as.character(Step), across(-Step, comma))
    dt(f, page = 20)
  })
  output$cause_plot <- renderPlot(res = 96, {
    u <- unknown_causes(pop_data(), as.integer(input$cause_phase))
    validate(need(nrow(u) > 0, "No girls with unknown status."))
    ggplot(u, aes(pct, arm, fill = cause)) + geom_col(width = .6) +
      facet_wrap(~ country_f(country), ncol = 1) +
      scale_x_continuous(labels = scales::percent) +
      scale_fill_manual(values = c("#2a78d6", "#e3a33b", "#8a5cc2", "#5aa469", "#c9473e")) +
      labs(x = "Share of girls with unknown status", y = NULL, fill = NULL) +
      guides(fill = guide_legend(nrow = 2))
  })

  # ---- Baseline ----------------------------------------------------------------
  output$t1_title <- renderText(paste("Table 1. Baseline characteristics:", sel_label()))
  output$t1_tbl <- renderDT(dt(table1(pop_data()), page = 30))

  # ---- Coverage ---------------------------------------------------------------
  cov_data <- reactive(outcome_data(pop_data(), input$outcome))
  crude <- reactive(crude_by_country(cov_data(), outcomes[[input$outcome]]$var))
  output$crude_plot <- renderPlot(res = 96, {
    cr <- crude()
    validate(need(nrow(cr) > 0, "Not estimable."))
    d <- bind_rows(cr %>% transmute(stratum, arm = "Intervention", cov = cov_int, lo = cov_int_lo, hi = cov_int_hi),
                   cr %>% transmute(stratum, arm = "Control", cov = cov_ctl, lo = cov_ctl_lo, hi = cov_ctl_hi)) %>%
      mutate(stratum = factor(country_label[stratum], levels = country_label),
             arm = factor(arm, levels = c("Control", "Intervention")))
    ggplot(d, aes(stratum, cov, fill = arm)) +
      geom_col(position = position_dodge(.75), width = .7) +
      geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(.75), width = .15, colour = "grey20") +
      geom_text(aes(label = pct1(cov), y = pmax(cov, .08)), position = position_dodge(.75), vjust = 1.6, size = 3.6, colour = "white") +
      scale_fill_manual(values = pal) + scale_y_continuous(labels = scales::percent, expand = expansion(c(0, .12))) +
      labs(x = NULL, y = "Vaccinated (girls with known status)", fill = NULL)
  })
  output$crude_tbl <- renderDT({
    cr <- crude()
    validate(need(nrow(cr) > 0, "Not estimable."))
    dt(cr %>% transmute(Country = country_label[stratum],
                        Intervention = sprintf("%s/%s, %s", comma(vacc_int), comma(n_int), fmt_pp(cov_int, cov_int_lo, cov_int_hi)),
                        Control = sprintf("%s/%s, %s", comma(vacc_ctl), comma(n_ctl), fmt_pp(cov_ctl, cov_ctl_lo, cov_ctl_hi)),
                        `Difference, pp (95% CI)` = fmt_pp(diff, diff_lo, diff_hi),
                        `RR (95% CI)` = fmt_ci(RR, RR_lo, RR_hi), p = fmt_p(p), Test = test))
  })
  adj <- reactive({
    k <- paste("adj", key0(), input$pop, input$stratum, input$outcome)
    cached(k, withProgress(message = "Fitting the mixed-effects model", detail = "a few seconds",
                           adjusted_row(cov_data(), outcomes[[input$outcome]]$var)))
  })
  output$adj_box <- renderUI({
    a <- adj(); s <- a$summary
    if (is.null(s$OR) || is.na(s$OR)) return(p("Not estimable in this selection."))
    dg <- a$diag
    tagList(
      div(class = "d-flex flex-wrap gap-2",
        stat_tile("Adjusted OR", s$OR, s$OR_lo, s$OR_hi, paste("p =", fmt_p(s$p)), primary = TRUE),
        stat_tile("Risk difference, points", 100 * s$ARD, 100 * s$ARD_lo, 100 * s$ARD_hi, "marginal standardization"),
        stat_tile("Risk ratio", s$RR, s$RR_lo, s$RR_hi, paste("n =", comma(s$n)))),
      tags$table(class = "table table-sm small mt-2",
        tags$tr(tags$td("Vaccinated, intervention / control"),
                tags$td(sprintf("%s/%s and %s/%s", comma(s$events_int), comma(s$n_int), comma(s$events_ctl), comma(s$n_ctl)))),
        tags$tr(tags$td("ICC (CHW zone)"), tags$td(sprintf("%.3f", s$ICC))),
        if (nrow(dg)) tags$tr(tags$td("AUC; Hosmer-Lemeshow"),
                              tags$td(sprintf("%.3f; \u03c7\u00b2 %.1f (df %d), p %s", dg$AUC, dg$HL_chisq, dg$HL_df, fmt_p(dg$HL_p)))),
        if (nrow(dg)) tags$tr(tags$td("Max adjusted GVIF"), tags$td(sprintf("%.2f", dg$max_GVIF_adj))),
        if (nzchar(s$note)) tags$tr(tags$td("Note"), tags$td(s$note))))
  })

  # ---- Time to vaccination ------------------------------------------------------
  tte <- reactive({
    k <- paste("tte", key0(), input$pop, input$stratum, input$tte_phase)
    cached(k, withProgress(message = "Fitting survival models",
                           tte_analysis(pop_data(), as.integer(input$tte_phase), input$stratum)))
  })
  output$km_plot <- renderPlot(res = 96, plot_km(km_curves(pop_data(), as.integer(input$tte_phase)), as.integer(input$tte_phase)))
  output$tte_box <- renderUI({
    r <- tte(); s <- r$surv; cx <- r$cox$summary
    rm <- s$rmst
    rows <- list(
      c("Girls (vaccinated during follow-up)", sprintf("%s (%s)", comma(r$n), comma(r$events))),
      c("Horizon \u03c4, days", r$tau),
      c("RMST, intervention / control (days unvaccinated)",
        sprintf("%.2f / %.2f", rm$rmean[rm$arm == "Intervention"], rm$rmean[rm$arm == "Control"])),
      c("Median days to vaccination, intervention / control",
        paste(med(rm$median[rm$arm == "Intervention"]), "/", med(rm$median[rm$arm == "Control"]))),
      c("Log-rank p (stratified by country)", fmt_p(s$logrank_p)),
      c("Cox adjusted HR (robust SE by zone)", sprintf("%s, p %s", fmt_ci(cx$HR, cx$HR_lo, cx$HR_hi), fmt_p(cx$p))),
      c("Proportional hazards test for arm, p", fmt_p(cx$PH_p_arm)))
    if (!is.null(r$period)) {
      ps <- r$period$summary
      rows <- c(rows, list(
        c("HR before 2 July", fmt_ci(ps$HR[1], ps$HR_lo[1], ps$HR_hi[1])),
        c("HR after 2 July", fmt_ci(ps$HR[2], ps$HR_lo[2], ps$HR_hi[2])),
        c("Arm \u00d7 period p", fmt_p(r$period$interaction_p))))
    }
    tagList(
      stat_tile("\u0394RMST, days (intervention \u2212 control)", s$dRMST$dRMST, s$dRMST$lo, s$dRMST$hi,
                "negative = vaccinated earlier in the intervention arm", primary = TRUE),
      tags$table(class = "table table-sm small mt-2",
                 map(rows, ~ tags$tr(tags$td(.x[1]), tags$td(.x[2])))),
      if (input$stratum == "Overall")
        note("The pooled estimate is not standardized to the country mix here; the report's Table 4 gives the",
             "country-standardized \u0394RMST."))
  })
  output$tte_cum <- renderDT({
    cm <- tte()$surv$cum
    dt(cm %>% transmute(Week = round(day / 7), Day = day, Arm = arm, `Vaccinated, % (95% CI)` = fmt_pp(cum_vacc, lo, hi)) %>%
         arrange(Week, Arm), page = 10)
  })

  # ---- Subgroups --------------------------------------------------------------
  sg <- reactive({
    k <- paste("sg", key0(), input$pop, input$sg_outcome)
    d <- select_pop(dat()$an, input$pop, "Overall") %>% outcome_data(input$sg_outcome)
    cached(k, withProgress(message = "Fitting subgroup models", value = 0, {
      res <- run_subgroups(d, "HPV_vaccinated", progress = function(x) incProgress(.3, detail = x)) %>% flag_subgroups()
      list(sg = res, overall = fit_glmer(d, "HPV_vaccinated")$summary)
    }))
  })
  output$sg_plot <- renderPlot(res = 96, {
    p <- plot_subgroups(sg()$sg, sg()$overall)
    validate(need(!is.null(p), "No subgroup has a stable estimate."))
    p
  })
  output$sg_tbl <- renderDT(dt(sg()$sg %>% filter(!plotted) %>%
                                 select(Subgroup = analysis, n, n_int, n_ctl, events_int, events_ctl, Reason = why_not)))

  # ---- Zones --------------------------------------------------------------------
  output$zone_plot <- renderPlot(res = 96, {
    d <- outcome_data(pop_data(), input$zone_outcome)
    z <- zone_coverage(d, outcomes[[input$zone_outcome]]$var)
    validate(need(nrow(z) > 0, "No girls with known status."))
    ord <- z %>% filter(arm == "Control") %>% arrange(cov) %>% pull(chw_zone) %>% as.character()
    z <- z %>% mutate(chw_zone = factor(chw_zone, levels = unique(c(ord, as.character(chw_zone)))))
    ggplot(z, aes(cov, chw_zone, colour = arm)) +
      geom_line(aes(group = chw_zone), colour = "grey75") +
      geom_point(aes(size = n), alpha = .9) + scale_size_area(max_size = 4) +
      scale_colour_manual(values = pal) + scale_x_continuous(labels = scales::percent) +
      facet_grid(country_f(country) ~ ., scales = "free_y", space = "free_y") +
      labs(x = "Vaccinated (girls with known status)", y = "CHW supervision zone", colour = NULL, size = "Girls") +
      theme(axis.text.y = element_text(size = 8))
  })

  # ---- Report tables -----------------------------------------------------------
  read_rep <- function(f) {
    p <- file.path(src(), f)
    validate(need(file.exists(p), sprintf("%s not found: render 02_prespecified_analysis.Rmd for this data source first.", f)))
    read_csv(p, show_col_types = FALSE)
  }
  output$sens_plot <- renderPlot(res = 96, plot_sensitivity(read_rep("table7_sensitivity_summary.csv")))
  output$sens_tbl <- renderDT({
    t7 <- read_rep("table7_sensitivity_summary.csv")
    dt(t7 %>% transmute(Scenario = scenario, n, `OR (95% CI)` = ifelse(is.na(OR), "", fmt_ci(OR, OR_lo, OR_hi)),
                        p = fmt_p(p), `Difference, pp (95% CI)` = ifelse(is.na(ARD), "", fmt_pp(ARD, ARD_lo, ARD_hi)),
                        `OR change vs primary` = ifelse(is.na(OR), "", sprintf("%+.1f%%", 100 * (OR / t7$OR[1] - 1)))),
       page = 30)
  })
  output$rep_view <- renderDT(dt(read_rep(input$rep_tbl), digits = 3, page = 25))
  output$rep_path <- renderText(file.path(src(), input$rep_tbl))
}

shinyApp(ui, server)
