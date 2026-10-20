# =============================================================================
# HPV community follow-up trial: public demo of the results explorer (Shiny)
#
# Reads ONLY data/summaries.rds, made by export_public_summaries.R from the
# SYNTHETIC data. It shows how the trial is analysed; every number is invented.
# Real trial results are never published through this app.
# Locally:  shiny::runApp("app_public")
# =============================================================================
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(DT); library(ggplot2); library(dplyr)
})
source("public_helpers.R")
S <- readRDS("data/summaries.rds")
if (!isTRUE(S$meta$synthetic)) stop("This public app only runs on summaries of the synthetic data.")

theme_set(theme_minimal(base_size = 13) +
            theme(panel.grid.minor = element_blank(), legend.position = "top",
                  strip.text = element_text(face = "bold", size = 12)))
dt <- function(x, digits = NULL, page = 25) {
  out <- datatable(x, rownames = FALSE, options = list(pageLength = page, scrollX = TRUE, dom = "tip"))
  num <- names(x)[vapply(x, is.double, logical(1))]
  if (length(num) && !is.null(digits)) out <- formatRound(out, num, digits)
  out
}
note <- function(...) div(class = "text-muted small mt-2", ...)

ui <- page_navbar(
  title = "HPV trial analysis: demo on synthetic data",
  fillable = FALSE,
  theme = bs_theme(version = 5, primary = "#2a78d6", base_font = font_collection("system-ui", "sans-serif")),
  sidebar = sidebar(
    width = 290,
    div(class = "alert alert-warning py-2 small", strong("Synthetic data."),
        " Every number in this app is invented: it shows how the trial is analysed, not its results."),
    radioButtons("pop", "Population", choices = setNames(names(pop_label), pop_label)),
    radioButtons("stratum", "Country", choices = setNames(names(country_label), country_label)),
    note("Girls expected to age out of eligibility within three months were allocated to the",
         "intervention arm without randomization. Crude comparisons pooling both countries",
         "are therefore confounded; the estimate in randomized girls takes precedence for",
         "causal interpretation."),
    note(sprintf("Synthetic summaries generated %s.", format(S$meta$created, "%d %B %Y")))
  ),

  nav_panel("About",
    card(card_header("About this demonstration"),
      div(class = "alert alert-warning", strong("All results shown here are invented."),
          " They come from synthetic data with the same structure as the trial data. They say nothing",
          " about the effect of the intervention. The trial results will be reported in the peer-reviewed publication."),
      p("The app demonstrates the pre-specified analysis of a pragmatic, individually randomized trial",
        "with a waitlist crossover. The trial tests community health worker (CHW) follow-up visits,",
        "supported by a digital package, to improve HPV vaccination of adolescent girls in Mali and",
        "C\u00f4te d'Ivoire (Muso Health)."),
      tags$ul(
        tags$li("Primary outcome: HPV vaccination status on 1 July 2026 (12 weeks after launch)."),
        tags$li("Secondary outcome: time to HPV vaccination within 12 weeks."),
        tags$li("Phase 2 (exploratory): status on 1 October 2026, after the control arm received the intervention from 2 July."),
        tags$li("Primary model: mixed-effects logistic regression with arm, country, age (spline) and a random intercept for CHW supervision zone.")),
      p("Registration: ", tags$a("ClinicalTrials.gov NCT07807202", href = "https://clinicaltrials.gov/study/NCT07807202", target = "_blank"),
        ". Analysis code, synthetic data generator and statistical analysis plan: ",
        tags$a("github.com/ABSANDIE/hpv-waitlist-rct", href = "https://github.com/ABSANDIE/hpv-waitlist-rct", target = "_blank"), "."))),

  nav_panel("Participant flow",
    card(card_header("CONSORT counts by country and arm"), DTOutput("flow_tbl", fill = FALSE)),
    card(card_header("Causes of unknown status (share of the analysed population)"),
         radioButtons("cause_phase", NULL, inline = TRUE, choices = c("Phase 1 (1 July)" = "1", "Phase 2 (1 October)" = "2")),
         plotOutput("cause_plot", height = 420))),

  nav_panel("Baseline",
    card(card_header(textOutput("t1_title")), DTOutput("t1_tbl", fill = FALSE),
         note("SMD = standardized mean difference (intervention vs control)."))),

  nav_panel("Coverage and uptake",
    selectInput("outcome", "Outcome", choices = outcome_choices, width = "520px"),
    layout_columns(col_widths = breakpoints(sm = c(12, 12), lg = c(5, 7)), fill = FALSE,
      card(card_header("Crude coverage by arm (Wilson 95% CI)"), plotOutput("crude_plot", height = 330)),
      card(card_header("Adjusted analysis (pre-specified model)"), uiOutput("adj_box"),
           note("ARD and RR by marginal standardization."))),
    card(card_header("Crude comparison"), DTOutput("crude_tbl", fill = FALSE))),

  nav_panel("Time to vaccination",
    radioButtons("tte_phase", NULL, inline = TRUE,
                 choices = c("12-week time frame (secondary outcome)" = "1",
                             "Full follow-up to 1 October (phase 2, exploratory)" = "2")),
    card(card_header("Cumulative proportion vaccinated, girls unvaccinated at launch (Kaplan-Meier)"),
         plotOutput("km_plot", height = 380)),
    layout_columns(col_widths = c(6, 6), fill = FALSE,
      card(card_header("RMST and hazard ratio"), uiOutput("tte_box")),
      card(card_header("Cumulative proportion vaccinated by arm"), DTOutput("tte_cum", fill = FALSE)))),

  nav_panel("Subgroups",
    radioButtons("sg_outcome", NULL, inline = TRUE,
                 choices = c("Coverage at 12 weeks" = "cov_p1",
                             "Uptake at 12 weeks, unvaccinated at launch (post hoc)" = "upt_p1")),
    card(card_header("Subgroup-specific adjusted ORs (both countries)"), plotOutput("sg_plot", height = 520),
         note("Exploratory; 95% CIs without multiplicity adjustment. Subgroups with fewer than 5",
              "vaccinated girls in an arm or an unstable fit are not plotted."))),

  nav_panel("CHW zones",
    selectInput("zone_outcome", "Outcome", choices = outcome_choices, width = "520px"),
    card(card_header("Crude coverage by CHW supervision zone and arm"), plotOutput("zone_plot", height = 700),
         note("Descriptive. Zones are anonymized and sorted by control-arm coverage; cells with fewer than",
              S$meta$min_cell, "girls are not shown."))),

  nav_panel("Sensitivity analyses",
    card(card_header("Sensitivity analyses of the primary outcome"),
         plotOutput("sens_plot", height = 620), DTOutput("sens_tbl", fill = FALSE),
         note("A material change is an OR differing by more than 10% or a change in significance."))),

  nav_panel("Tables",
    selectInput("rep_tbl", "Table", choices = names(S$report), width = "520px"),
    card(DTOutput("rep_view", fill = FALSE)))
)

server <- function(input, output, session) {
  pick <- function(x) x[[input$pop]][[input$stratum]]
  sel_label <- reactive(paste(pop_label[input$pop], country_label[input$stratum], sep = ", "))

  output$flow_tbl <- renderDT(dt(S$flow %>% mutate(across(-Step, comma)), page = 20))
  output$cause_plot <- renderPlot(res = 96, {
    u <- pick(S$causes)[[input$cause_phase]]
    validate(need(nrow(u) > 0, "No girls with unknown status."))
    ggplot(u, aes(pct, arm, fill = cause)) + geom_col(width = .6) +
      facet_wrap(~ country_f(country), ncol = 1) + scale_x_continuous(labels = scales::percent) +
      scale_fill_manual(values = c("#2a78d6", "#e3a33b", "#8a5cc2", "#5aa469", "#c9473e")) +
      labs(x = "Share of girls with unknown status", y = NULL, fill = NULL) + guides(fill = guide_legend(nrow = 2))
  })

  output$t1_title <- renderText(paste("Baseline characteristics:", sel_label()))
  output$t1_tbl <- renderDT(dt(pick(S$table1), page = 30))

  crude <- reactive(pick(S$crude)[[input$outcome]])
  output$crude_plot <- renderPlot(res = 96, {
    cr <- crude(); validate(need(nrow(cr) > 0, "Not estimable."))
    d <- bind_rows(cr %>% transmute(stratum, arm = "Intervention", cov = cov_int, lo = cov_int_lo, hi = cov_int_hi),
                   cr %>% transmute(stratum, arm = "Control", cov = cov_ctl, lo = cov_ctl_lo, hi = cov_ctl_hi)) %>%
      mutate(stratum = factor(country_label[stratum], levels = country_label),
             arm = factor(arm, levels = c("Control", "Intervention")))
    pd <- position_dodge(.8)
    ggplot(d, aes(cov, stratum, fill = arm)) +
      geom_col(position = pd, width = .75) +
      geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", position = pd, width = .2, colour = "grey20") +
      geom_text(aes(x = hi, label = pct1(cov), group = arm), position = pd, hjust = -0.2, size = 3.4, colour = "grey20") +
      scale_fill_manual(values = pal, breaks = c("Intervention", "Control")) +
      scale_x_continuous(labels = scales::percent, limits = c(0, 1.18), breaks = seq(0, 1, .25)) +
      scale_y_discrete(limits = rev) +
      labs(y = NULL, x = "Vaccinated (known status)", fill = NULL) +
      theme(panel.grid.major.y = element_blank())
  })
  output$crude_tbl <- renderDT({
    cr <- crude(); validate(need(nrow(cr) > 0, "Not estimable."))
    dt(cr %>% transmute(Country = country_label[stratum],
                        Intervention = sprintf("%s/%s, %s", comma(vacc_int), comma(n_int), fmt_pp(cov_int, cov_int_lo, cov_int_hi)),
                        Control = sprintf("%s/%s, %s", comma(vacc_ctl), comma(n_ctl), fmt_pp(cov_ctl, cov_ctl_lo, cov_ctl_hi)),
                        `Difference, pp (95% CI)` = fmt_pp(diff, diff_lo, diff_hi),
                        `RR (95% CI)` = fmt_ci(RR, RR_lo, RR_hi), p = fmt_p(p)))
  })
  output$adj_box <- renderUI({
    a <- pick(S$adj)[[input$outcome]]; s <- a$summary; dg <- a$diag
    if (is.null(s$OR) || is.na(s$OR)) return(p("Not estimable in this selection."))
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
                              tags$td(sprintf("%.3f; χ² %.1f (df %d), p %s", dg$AUC, dg$HL_chisq, dg$HL_df, fmt_p(dg$HL_p))))))
  })

  tte <- reactive(pick(S$tte)[[input$tte_phase]])
  output$km_plot <- renderPlot(res = 96, plot_km(pick(S$km)[[input$tte_phase]], as.integer(input$tte_phase)))
  output$tte_box <- renderUI({
    r <- tte(); rm <- r$rmst; cx <- r$cox
    rows <- list(
      c("Girls (vaccinated during follow-up)", sprintf("%s (%s)", comma(r$n), comma(r$events))),
      c("Horizon τ, days", r$tau),
      c("RMST, intervention / control (days unvaccinated)",
        sprintf("%.2f / %.2f", rm$rmean[rm$arm == "Intervention"], rm$rmean[rm$arm == "Control"])),
      c("Median days to vaccination, intervention / control",
        paste(med(rm$median[rm$arm == "Intervention"]), "/", med(rm$median[rm$arm == "Control"]))),
      c("Log-rank p (stratified by country)", fmt_p(r$logrank_p)),
      c("Cox adjusted HR (robust SE by zone)", sprintf("%s, p %s", fmt_ci(cx$HR, cx$HR_lo, cx$HR_hi), fmt_p(cx$p))),
      c("Proportional hazards test for arm, p", fmt_p(cx$PH_p_arm)))
    if (!is.null(r$period)) {
      ps <- r$period$summary
      rows <- c(rows, list(c("HR before 2 July", fmt_ci(ps$HR[1], ps$HR_lo[1], ps$HR_hi[1])),
                           c("HR after 2 July", fmt_ci(ps$HR[2], ps$HR_lo[2], ps$HR_hi[2])),
                           c("Arm × period p", fmt_p(r$period$interaction_p))))
    }
    tagList(
      stat_tile("ΔRMST, days (intervention − control)", r$dRMST$dRMST, r$dRMST$lo, r$dRMST$hi,
                "negative = vaccinated earlier in the intervention arm", primary = TRUE),
      tags$table(class = "table table-sm small mt-2", lapply(rows, function(x) tags$tr(tags$td(x[1]), tags$td(x[2])))))
  })
  output$tte_cum <- renderDT(dt(tte()$cum %>% transmute(Week = round(day / 7), Day = day, Arm = arm,
                                                        `Vaccinated, % (95% CI)` = fmt_pp(cum_vacc, lo, hi)) %>%
                                  arrange(Week, Arm), page = 10))

  output$sg_plot <- renderPlot(res = 96, {
    x <- S$sg[[input$pop]][[input$sg_outcome]]
    p <- plot_subgroups(x$sg, x$overall)
    validate(need(!is.null(p), "No subgroup has a stable estimate."))
    p
  })

  output$zone_plot <- renderPlot(res = 96, {
    z <- pick(S$zone)[[input$zone_outcome]]
    validate(need(nrow(z) > 0, "No zone has enough girls with known status."))
    z <- z %>% mutate(key = paste(country, zone))
    ord <- z %>% filter(arm == "Control") %>% arrange(cov) %>% pull(key)
    z <- z %>% mutate(key = factor(key, levels = unique(c(ord, key))))
    ggplot(z, aes(cov, key, colour = arm)) +
      geom_line(aes(group = key), colour = "grey75") +
      geom_point(aes(size = n), alpha = .9) + scale_size_area(max_size = 4) +
      scale_colour_manual(values = pal) + scale_x_continuous(labels = scales::percent) +
      scale_y_discrete(labels = function(k) sub("^\\S+ ", "", k)) +
      facet_grid(country_f(country) ~ ., scales = "free_y", space = "free_y") +
      labs(x = "Vaccinated (girls with known status)", y = "CHW supervision zone", colour = NULL, size = "Girls") +
      theme(axis.text.y = element_text(size = 8))
  })

  t7 <- S$report[["Table 7. Sensitivity analyses (summary)"]]
  output$sens_plot <- renderPlot(res = 96, { validate(need(!is.null(t7), "Not available.")); plot_sensitivity(t7) })
  output$sens_tbl <- renderDT({
    validate(need(!is.null(t7), "Not available."))
    dt(t7 %>% transmute(Scenario = scenario, n, `OR (95% CI)` = ifelse(is.na(OR), "", fmt_ci(OR, OR_lo, OR_hi)),
                        p = fmt_p(p), `Difference, pp (95% CI)` = ifelse(is.na(ARD), "", fmt_pp(ARD, ARD_lo, ARD_hi)),
                        `OR change vs primary` = ifelse(is.na(OR), "", sprintf("%+.1f%%", 100 * (OR / t7$OR[1] - 1)))),
       page = 30)
  })
  output$rep_view <- renderDT(dt(S$report[[input$rep_tbl]], digits = 3))
}

shinyApp(ui, server)
