# =============================================================================
# public_helpers.R - labels, formatting and plots for the public app.
# Self-contained copies of functions in R/hpv_functions.R and app/app_helpers.R,
# so that the deployed folder needs nothing else from the repository.
# =============================================================================
country_label <- c(Overall = "Both countries", Mali = "Mali", CIV = "Côte d'Ivoire")
pop_label     <- c(itt = "All girls (mITT)", randomized = "Randomized girls")
pal <- c(Control = "#6b6a66", Intervention = "#2a78d6")
outcome_choices <- c(
  "Coverage at 12 weeks (1 July 2026): primary outcome" = "cov_p1",
  "Uptake at 12 weeks among girls unvaccinated at launch" = "upt_p1",
  "Coverage on 1 October 2026 (phase 2, exploratory)" = "cov_p2",
  "Uptake by 1 October among girls unvaccinated at launch (phase 2)" = "upt_p2")

fmt_ci <- function(est, lo, hi, digits = 2, pct = FALSE) {
  k <- if (pct) 100 else 1
  sprintf(paste0("%.", digits, "f (%.", digits, "f to %.", digits, "f)"), k * est, k * lo, k * hi)
}
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_pp <- function(e, l, h) ifelse(is.na(e), "n.e.", fmt_ci(e, l, h, 2, pct = TRUE))
pct1   <- function(x) sprintf("%.1f%%", 100 * x)
comma  <- function(x) formatC(x, format = "d", big.mark = ",")
med    <- function(x) ifelse(is.na(x), "not reached", format(x))
country_f <- function(x) factor(country_label[as.character(x)], levels = country_label)

plot_km <- function(df, phase) {
  tmax <- if (phase == 1) 84 else 176
  df$country <- factor(as.character(df$country), levels = country_label)
  p <- ggplot(df, aes(time, cum, colour = arm)) +
    geom_step(linewidth = .8) + facet_wrap(~ country) +
    scale_colour_manual(values = pal) +
    scale_y_continuous(labels = scales::percent) + coord_cartesian(xlim = c(0, tmax)) +
    labs(x = "Days since launch", y = "Cumulative proportion vaccinated", colour = NULL)
  if (phase == 2) p <- p + geom_vline(xintercept = 84.5, linetype = 2, colour = "grey40") +
    annotate("text", x = 84.5, y = Inf, label = " control arm starts (2 July)", hjust = 0, vjust = 1.5, size = 3.5)
  p
}

plot_subgroups <- function(sg, overall) {
  d <- sg %>% filter(plotted) %>%
    mutate(lab = sprintf("%s  (n = %s; %d/%d vaccinated)", analysis, comma(n), events_int, events_ctl),
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

stat_tile <- function(title, est, lo, hi, sub = NULL, primary = FALSE) {
  div(class = paste("p-3 rounded flex-fill", if (primary) "text-white bg-primary" else "border"),
      style = "min-width: 170px;",
      div(class = "small", title),
      div(style = "font-size: 1.9rem; font-weight: 600; line-height: 1.2;", sprintf("%.2f", est)),
      div(class = "small", sprintf("95%% CI %.2f to %.2f", lo, hi)),
      if (!is.null(sub)) div(class = "small opacity-75", sub),
      div(class = "small opacity-75 fst-italic", "synthetic data"))
}


