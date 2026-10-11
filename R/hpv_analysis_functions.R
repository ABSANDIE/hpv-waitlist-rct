# =============================================================================
# hpv_analysis_functions.R
# Model-fitting helpers for 02_prespecified_analysis.Rmd (registered SAP v2.1).
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(lme4); library(survival)
  library(marginaleffects)
})

ctrl_glmer <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# Add the age spline (months remaining in the eligibility window at launch).
# Knots are fixed once on the full ITT population and reused in every model.
add_age_spline <- function(dat, knots) {
  b <- rcs_basis(dat$months_remaining, knots)
  dat$age_s1 <- b[, 1]; dat$age_s2 <- b[, 2]
  dat
}

# Right-hand side of the pre-specified model; drops country when only one level
rhs_primary <- function(dat, extra = NULL, interaction = FALSE) {
  terms <- c(if (interaction) "arm01 * country" else "arm01",
             if (!interaction && n_distinct(dat$country) > 1) "country",
             "age_s1", "age_s2", extra)
  paste(terms, collapse = " + ")
}

# ---- Mixed-effects logistic model: OR, ARD and RR (marginal standardization) --
fit_glmer <- function(dat, outcome, extra = NULL, interaction = FALSE, label = NULL) {
  dat <- dat[!is.na(dat[[outcome]]), , drop = FALSE]
  dat$y <- dat[[outcome]]
  out <- tibble(analysis = label %||% outcome, n = nrow(dat),
                n_int = sum(dat$arm01 == 1), n_ctl = sum(dat$arm01 == 0),
                events_int = sum(dat$y[dat$arm01 == 1]), events_ctl = sum(dat$y[dat$arm01 == 0]))
  if (nrow(dat) < 20 || length(unique(dat$y)) < 2 || length(unique(dat$arm01)) < 2)
    return(list(fit = NULL, summary = out %>% mutate(note = "not estimable")))
  f <- as.formula(paste("y ~", rhs_primary(dat, extra, interaction), "+ (1 | chw_zone)"))
  fit <- tryCatch(glmer(f, data = dat, family = binomial, control = ctrl_glmer), error = function(e) e)
  note <- ""
  if (inherits(fit, "error")) {          # fall back to a GLM with cluster-robust SE
    note <- paste("glmer failed:", conditionMessage(fit), "- GLM used")
    fit <- glm(update(f, . ~ . - (1 | chw_zone)), data = dat, family = binomial)
  } else if (isSingular(fit)) note <- "singular fit (zone variance ~ 0)"
  co <- if (inherits(fit, "merMod")) summary(fit)$coefficients else lmtest_like(fit, dat$chw_zone)
  b <- co["arm01", 1]; se <- co["arm01", 2]
  me_d <- tryCatch(avg_comparisons(fit, variables = "arm01"), error = function(e) NULL)
  me_r <- tryCatch(avg_comparisons(fit, variables = "arm01", comparison = "lnratioavg",
                                   transform = exp), error = function(e) NULL)
  icc <- if (inherits(fit, "merMod")) { v <- as.numeric(VarCorr(fit)$chw_zone); v / (v + pi^2 / 3) } else NA
  out <- out %>% mutate(
    OR = exp(b), OR_lo = exp(b - 1.96 * se), OR_hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)),
    ARD = me_d$estimate %||% NA, ARD_lo = me_d$conf.low %||% NA, ARD_hi = me_d$conf.high %||% NA,
    RR = me_r$estimate %||% NA, RR_lo = me_r$conf.low %||% NA, RR_hi = me_r$conf.high %||% NA,
    ICC = icc, logOR = b, se_logOR = se, note = note)
  list(fit = fit, summary = out)
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# GLM with cluster-robust (sandwich, CR0) SE - fallback only
lmtest_like <- function(fit, cluster) {
  X <- model.matrix(fit); r <- residuals(fit, type = "working") * fit$weights
  bread <- solve(crossprod(X * sqrt(fit$weights)))
  meat <- crossprod(rowsum(X * r, cluster))
  V <- bread %*% meat %*% bread
  cbind(Estimate = coef(fit), SE = sqrt(diag(V)))
}

# ---- Unadjusted comparison from the 2x2 table ---------------------------------
unadjusted <- function(dat, outcome) {
  d <- dat[!is.na(dat[[outcome]]), ]
  x1 <- sum(d[[outcome]][d$arm01 == 1]); n1 <- sum(d$arm01 == 1)
  x0 <- sum(d[[outcome]][d$arm01 == 0]); n0 <- sum(d$arm01 == 0)
  if (n1 == 0 || n0 == 0) return(tibble())
  tab <- matrix(c(x1, n1 - x1, x0, n0 - x0), 2, byrow = TRUE)
  expct <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
  tst <- if (any(expct < 5)) fisher.test(tab) else chisq.test(tab, correct = FALSE)
  a <- x1 + .5; b <- n1 - x1 + .5; c <- x0 + .5; dd <- n0 - x0 + .5     # Haldane for OR
  lor <- log(a * dd / (b * c)); se <- sqrt(1 / a + 1 / b + 1 / c + 1 / dd)
  rr <- (x1 / n1) / (x0 / n0); se_rr <- sqrt(1 / max(x1, .5) - 1 / n1 + 1 / max(x0, .5) - 1 / n0)
  nd <- newcombe_diff(x1, n1, x0, n0)
  w1 <- wilson_ci(x1, n1); w0 <- wilson_ci(x0, n0)
  tibble(n_int = n1, vacc_int = x1, cov_int = w1[1], cov_int_lo = w1[2], cov_int_hi = w1[3],
         n_ctl = n0, vacc_ctl = x0, cov_ctl = w0[1], cov_ctl_lo = w0[2], cov_ctl_hi = w0[3],
         diff = nd[1], diff_lo = nd[2], diff_hi = nd[3],
         RR = rr, RR_lo = exp(log(rr) - 1.96 * se_rr), RR_hi = exp(log(rr) + 1.96 * se_rr),
         OR = exp(lor), OR_lo = exp(lor - 1.96 * se), OR_hi = exp(lor + 1.96 * se),
         test = if (grepl("Fisher", tst$method)) "Fisher" else "Chi-squared",
         p = tst$p.value)
}

# ---- Model diagnostics: AUC, Hosmer-Lemeshow (10 groups), VIF -----------------
diagnostics <- function(fit) {
  if (is.null(fit)) return(tibble())
  y <- if (inherits(fit, "merMod")) getME(fit, "y") else fit$y
  p <- fitted(fit)
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  auc <- (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  g <- cut(p, unique(quantile(p, seq(0, 1, .1))), include.lowest = TRUE)
  obs <- tapply(y, g, sum); ex <- tapply(p, g, sum); nn <- tapply(y, g, length)
  hl <- sum((obs - ex)^2 / (ex * (1 - ex / nn))); dfhl <- length(nn) - 2
  X <- if (inherits(fit, "merMod")) getME(fit, "X") else model.matrix(fit)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  X <- X[, apply(X, 2, var) > 0, drop = FALSE]
  # Generalized VIF (Fox & Monette 1992); the two spline columns form one term
  grp <- ifelse(grepl("^age_s", colnames(X)), "age (spline)", colnames(X))
  R <- cor(X); dR <- det(R)
  vif <- sapply(unique(grp), function(g) {
    i <- grp == g
    if (all(i)) return(1)
    gv <- det(R[i, i, drop = FALSE]) * det(R[!i, !i, drop = FALSE]) / dR
    gv^(1 / (2 * sum(i)))
  })
  tibble(AUC = auc, HL_chisq = hl, HL_df = dfhl, HL_p = pchisq(hl, dfhl, lower.tail = FALSE),
         max_GVIF_adj = max(vif), `GVIF^(1/(2df))` = paste(sprintf("%s=%.2f", names(vif), vif), collapse = "; "))
}

# ---- Kaplan-Meier coverage: pi + (1 - pi) * (1 - S(tau)) --------------------
km_coverage <- function(dat, tau, time = "time_p1", event = "event_p1", B = 0) {
  one <- function(d) {
    known <- d[!is.na(d$prior_vacc), ]
    pi_hat <- mean(known$prior_vacc == 1); n_pi <- nrow(known)
    s <- d[d$prior_vacc %in% 0L, ]
    if (nrow(s) == 0) return(c(cov = pi_hat, var = pi_hat * (1 - pi_hat) / n_pi))
    fit <- survfit(Surv(s[[time]], s[[event]]) ~ 1)
    sm <- summary(fit, times = tau, extend = TRUE)
    S <- sm$surv; vS <- sm$std.err^2
    c(cov = pi_hat + (1 - pi_hat) * (1 - S),
      var = (1 - pi_hat)^2 * vS + S^2 * pi_hat * (1 - pi_hat) / n_pi)
  }
  a1 <- one(dat[dat$arm01 == 1, ]); a0 <- one(dat[dat$arm01 == 0, ])
  res <- tibble(cov_int = a1[["cov"]], var_int = a1[["var"]], cov_ctl = a0[["cov"]], var_ctl = a0[["var"]])
  bt <- if (B > 0) replicate(B, {   # bootstrap: resample girls within arm
    b1 <- dat[dat$arm01 == 1, ][sample(sum(dat$arm01 == 1), replace = TRUE), ]
    b0 <- dat[dat$arm01 == 0, ][sample(sum(dat$arm01 == 0), replace = TRUE), ]
    one(b1)[["cov"]] - one(b0)[["cov"]]
  }) else NULL
  attr(res, "boot") <- bt
  res
}
# By country (country-specific tau), plus an overall estimate standardized to
# the country distribution of the population (same weights in both arms)
km_coverage_std <- function(dat, tau_var, time, event, B = 0) {
  cs <- intersect(c("Mali", "CIV"), unique(as.character(dat$country)))
  per <- map(cs, function(cn) { d <- dat[dat$country == cn, ]
                                km_coverage(d, unique(d[[tau_var]])[1], time, event, B) })
  w <- map_dbl(cs, ~ mean(dat$country == .x))
  fin <- function(r, bt) {
    d <- r$cov_int - r$cov_ctl; se <- sqrt(r$var_int + r$var_ctl)
    out <- r %>% transmute(cov_int, cov_ctl, diff = d, diff_lo = d - 1.96 * se, diff_hi = d + 1.96 * se)
    if (!is.null(bt)) out <- out %>% mutate(boot_lo = quantile(bt, .025), boot_hi = quantile(bt, .975))
    out
  }
  rows <- map2_dfr(per, cs, ~ fin(.x, attr(.x, "boot")) %>%
                     mutate(stratum = .y, tau = unique(dat[[tau_var]][dat$country == .y])[1], .before = 1))
  if (length(cs) > 1) {
    ov <- tibble(cov_int = sum(w * map_dbl(per, "cov_int")), var_int = sum(w^2 * map_dbl(per, "var_int")),
                 cov_ctl = sum(w * map_dbl(per, "cov_ctl")), var_ctl = sum(w^2 * map_dbl(per, "var_ctl")))
    bt <- if (B > 0) Reduce(`+`, map2(per, w, ~ .y * attr(.x, "boot"))) else NULL
    rows <- bind_rows(fin(ov, bt) %>% mutate(stratum = "Overall (country-standardized)", tau = NA, .before = 1), rows)
  }
  rows
}

# ---- Survival summaries: cumulative proportion, median, RMST -----------------
surv_by_arm <- function(dat, tau, time = "time_p1", event = "event_p1", weeks = c(4, 8, 12)) {
  f <- as.formula(sprintf("Surv(%s, %s) ~ arm", time, event))
  fit <- survfit(f, data = dat)
  tab <- summary(fit, rmean = tau)$table
  rm <- tibble(arm = sub("arm=", "", rownames(tab)), rmean = tab[, "rmean"], se_rmean = tab[, "se(rmean)"],
               median = tab[, "median"], median_lo = tab[, "0.95LCL"], median_hi = tab[, "0.95UCL"])
  times <- pmin(weeks * 7, tau)
  cp <- summary(fit, times = times, extend = TRUE)
  cum <- tibble(arm = sub("arm=", "", as.character(cp$strata)), day = cp$time,
                cum_vacc = 1 - cp$surv, lo = 1 - cp$upper, hi = 1 - cp$lower)
  d <- rm$rmean[rm$arm == "Intervention"] - rm$rmean[rm$arm == "Control"]
  se <- sqrt(sum(rm$se_rmean^2))
  drm <- tibble(tau = tau, dRMST = d, lo = d - 1.96 * se, hi = d + 1.96 * se, p = 2 * pnorm(-abs(d / se)))
  if (requireNamespace("survRM2", quietly = TRUE)) {      # Uno et al. (2014), as in the SAP
    r2 <- tryCatch(survRM2::rmst2(dat[[time]], dat[[event]], dat$arm01, tau = tau), error = function(e) NULL)
    if (!is.null(r2)) drm <- drm %>% mutate(dRMST_survRM2 = r2$unadjusted.result[1, 1],
                                            lo_survRM2 = r2$unadjusted.result[1, 2],
                                            hi_survRM2 = r2$unadjusted.result[1, 3])
  }
  strata_term <- if (n_distinct(dat$country) > 1) " + strata(country)" else ""
  lr <- survdiff(as.formula(sprintf("Surv(%s, %s) ~ arm%s", time, event, strata_term)), data = dat)
  list(fit = fit, rmst = rm, cum = cum, dRMST = drm,
       logrank_p = pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE))
}

# Cox model with cluster-robust SE at the CHW-zone level, and PH checks
fit_cox <- function(dat, time = "time_p1", event = "event_p1") {
  rhs <- paste(c("arm01", if (n_distinct(dat$country) > 1) "country", "age_s1", "age_s2"), collapse = " + ")
  f <- as.formula(sprintf("Surv(%s, %s) ~ %s", time, event, rhs))
  fit <- coxph(f, data = dat, cluster = chw_zone)
  s <- summary(fit)$coefficients
  zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
  list(fit = fit, zph = zph,
       summary = tibble(HR = exp(s["arm01", "coef"]),
                        HR_lo = exp(s["arm01", "coef"] - 1.96 * s["arm01", "robust se"]),
                        HR_hi = exp(s["arm01", "coef"] + 1.96 * s["arm01", "robust se"]),
                        p = s["arm01", "Pr(>|z|)"],
                        PH_p_arm = if (!is.null(zph)) zph$table["arm01", "p"] else NA,
                        PH_p_global = if (!is.null(zph)) zph$table["GLOBAL", "p"] else NA))
}

# Split follow-up at girl-specific cut points (days) and fit arm x period
split_periods <- function(dat, cuts, time, event) {
  # cuts: list (one element per girl) of cut points in days; girls with time 0 contribute nothing
  n <- nrow(dat); tm <- dat[[time]]; ev <- dat[[event]]
  cl <- tibble(row = rep(seq_len(n), lengths(cuts)), b = unlist(cuts)) %>%
    filter(b > 0, b < tm[row])
  sp <- bind_rows(tibble(row = seq_len(n), b = 0), cl, tibble(row = seq_len(n), b = tm)) %>%
    distinct() %>% arrange(row, b) %>% group_by(row) %>%
    mutate(tstart = lag(b), tstop = b) %>% filter(!is.na(tstart)) %>%
    mutate(period = row_number(), ev = if_else(row_number() == n(), ev[row], 0)) %>% ungroup()
  cbind(dat[sp$row, setdiff(names(dat), c("tstart", "tstop", "ev", "period")), drop = FALSE],
        as.data.frame(sp[, c("tstart", "tstop", "ev", "period")]))
}
fit_cox_period <- function(sp, labels) {
  sp$period <- factor(labels[sp$period], levels = labels)
  rhs <- paste(c("arm01:period", if (n_distinct(sp$country) > 1) "country", "age_s1", "age_s2"), collapse = " + ")
  fit <- coxph(as.formula(paste("Surv(tstart, tstop, ev) ~", rhs)), data = sp, cluster = chw_zone)
  s <- summary(fit)$coefficients
  rr <- grep("^arm01:period", rownames(s))
  out <- tibble(period = sub("arm01:period", "", rownames(s)[rr]), HR = exp(s[rr, "coef"]),
                HR_lo = exp(s[rr, "coef"] - 1.96 * s[rr, "robust se"]),
                HR_hi = exp(s[rr, "coef"] + 1.96 * s[rr, "robust se"]), p = s[rr, "Pr(>|z|)"])
  # Wald test that the HRs are equal across periods (interaction)
  b <- coef(fit)[rr]; V <- vcov(fit)[rr, rr]
  L <- cbind(-1, diag(length(rr) - 1))
  w <- if (length(rr) > 1) as.numeric(t(L %*% b) %*% solve(L %*% V %*% t(L)) %*% (L %*% b)) else NA
  list(fit = fit, summary = out,
       interaction_p = if (length(rr) > 1) pchisq(w, length(rr) - 1, lower.tail = FALSE) else NA)
}

# ---- Rubin's rules on the log-OR ---------------------------------------------
pool_rubin <- function(est, se) {
  m <- length(est); qbar <- mean(est); ubar <- mean(se^2); b <- var(est)
  tt <- ubar + (1 + 1 / m) * b
  df <- if (b > 0) (m - 1) * (1 + ubar / ((1 + 1 / m) * b))^2 else Inf
  q <- qt(.975, df)
  tibble(m = m, OR = exp(qbar), OR_lo = exp(qbar - q * sqrt(tt)), OR_hi = exp(qbar + q * sqrt(tt)),
         p = 2 * pt(-abs(qbar / sqrt(tt)), df), fmi = (1 + 1 / m) * b / tt)
}

# Approximate-Bayesian logistic imputation with a delta shift (tipping point)
# Draws coefficients from N(beta_hat, V), adds `delta` to the linear predictor
# of girls flagged in `shift`, and draws y ~ Bernoulli.
impute_delta <- function(dat, formula_rhs, shift, delta, m) {
  obs <- !is.na(dat$y)
  fit <- glm(as.formula(paste("y ~", formula_rhs)), data = dat[obs, ], family = binomial)
  X <- model.matrix(as.formula(paste("~", formula_rhs)), data = dat)
  X <- X[, names(coef(fit))[!is.na(coef(fit))], drop = FALSE]
  bh <- coef(fit)[colnames(X)]; V <- vcov(fit)[colnames(X), colnames(X)]
  lapply(seq_len(m), function(k) {
    bk <- MASS::mvrnorm(1, bh, V)
    lp <- drop(X[!obs, , drop = FALSE] %*% bk) + delta * shift[!obs]
    d <- dat; d$y[!obs] <- rbinom(sum(!obs), 1, plogis(lp)); d
  })
}

# ---- Manski bounds ----------------------------------------------------------
# For girls with unknown status whose cause is in `causes`: worst case for the
# intervention (intervention = 0, control = 1) and best case (reverse).
manski_data <- function(dat, outcome, reason, causes, case = c("worst", "best")) {
  case <- match.arg(case)
  d <- dat; y <- d[[outcome]]
  fill <- is.na(y) & d[[reason]] %in% causes
  y[fill] <- if (case == "worst") as.integer(d$arm01[fill] == 0) else as.integer(d$arm01[fill] == 1)
  d$y_bound <- y
  d
}

# ---- Formatting --------------------------------------------------------------
fmt_or <- function(s) if (is.null(s$OR) || is.na(s$OR)) "n.e." else fmt_ci(s$OR, s$OR_lo, s$OR_hi, 2)
fmt_pp <- function(e, l, h) ifelse(is.na(e), "n.e.", fmt_ci(e, l, h, 2, pct = TRUE))
km_df <- function(fit) {
  s <- summary(fit, censored = FALSE)
  st <- if (is.null(s$strata)) "all" else as.character(s$strata)
  d <- tibble(strata = st, time = s$time, cum = 1 - s$surv, lo = 1 - s$upper, hi = 1 - s$lower)
  start <- d %>% distinct(strata) %>% mutate(time = 0, cum = 0, lo = 0, hi = 0)
  bind_rows(start, d) %>% arrange(strata, time) %>% mutate(arm = sub(".*arm=", "", strata))
}
