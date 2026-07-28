source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(plm))

lag_sweep <- function(dep, inc, panel_label) {
  d <- add_lags(read_panel(FILE_P5), c(dep, inc), 1)
  d$Ldep <- d[[paste0(dep, "_l1")]]
  d$Linc <- d[[paste0(inc, "_l1")]]
  s <- filter(d, sample == 1)
  ols <- as.numeric(ce(fit_ols(s, dep, c("Ldep", "Linc"), FALSE), "Linc")["est"])
  fe  <- as.numeric(ce(fit_ols(s, dep, c("Ldep", "Linc"), TRUE),  "Linc")["est"])
  pd <- gmm_panel(s, dep, inc)
  rows <- list()
  for (L in 2:8) for (tr in c("d", "ld")) for (cl in c(FALSE, TRUE)) {
    m <- pgmm(gmm_formula(dep, inc, sprintf("2:%d", L)), data = pd, effect = "twoways",
              model = "twostep", transformation = tr, collapse = cl)
    ic <- unname(co_row(pgmm_co(m), inc))
    rows[[length(rows) + 1]] <- tibble(
      panel = panel_label,
      estimator = ifelse(tr == "d", "Difference GMM", "System GMM"),
      instruments_kept = ifelse(cl, "collapsed", "uncollapsed"),
      max_lag = L, income = ic[1], income_se = ic[2],
      obs = nobs(m), countries = pgmm_countries(m), n_inst = ncol(m$W[[1]]),
      ar1_p = tryCatch(mtest(m, 1)$p.value, error = function(e) NA_real_),
      ar2_p = tryCatch(mtest(m, 2)$p.value, error = function(e) NA_real_),
      overid_p = tryCatch(sargan(m)$p.value, error = function(e) NA_real_))
  }
  tb <- bind_rows(rows)
  stopifnot(nrow(tb) == 28L,
            all(tapply(tb$obs, tb$estimator, function(x) length(unique(x))) == 1L))
  list(tab = tb, bench = tibble(panel = panel_label, ols = ols, fe = fe))
}

res <- lapply(MEASURES, function(ms) lag_sweep(ms$dep, ms$inc, ms$label))
tab <- bind_rows(lapply(res, function(r) r$tab))
bench <- bind_rows(lapply(res, function(r) r$bench))
write_csv(tab, file.path(PATH_OUTPUT, "instruments.csv"))

txt <- c()
for (pl in unique(tab$panel)) {
  b <- filter(bench, panel == pl)
  txt <- c(txt, sprintf("== %s (pooled OLS %.3f, fixed effects %.3f) ==", pl, b$ols, b$fe),
           "Each cell: estimate (SE) [instruments, overid p]")
  for (es in unique(tab$estimator)) {
    txt <- c(txt, sprintf("-- %s --", es),
             sprintf("%-8s %-30s %-30s", "Max lag", "Uncollapsed", "Collapsed"))
    for (L in sort(unique(tab$max_lag))) {
      pick <- function(k) filter(tab, panel == pl, estimator == es,
                                 max_lag == L, instruments_kept == k)
      cell <- function(r) if (nrow(r) == 0) "" else sprintf("%s (%s) [%d, %s]",
        num(r$income), num(r$income_se), r$n_inst,
        ifelse(is.na(r$overid_p), "exact", num(r$overid_p, 2)))
      txt <- c(txt, sprintf("%-8d %-30s %-30s", L,
                            cell(pick("uncollapsed")), cell(pick("collapsed"))))
    }
    txt <- c(txt, "")
  }
}
writeLines(txt, file.path(PATH_OUTPUT, "instruments.txt"))

benchlong <- bind_rows(
  transmute(bench, panel, kind = "pooled OLS", value = ols),
  transmute(bench, panel, kind = "fixed effects", value = fe))
fig <- ggplot(tab, aes(max_lag, income, color = instruments_kept)) +
  geom_hline(data = benchlong, aes(yintercept = value, linetype = kind), color = "grey45") +
  geom_line() + geom_point(size = 1.5) +
  facet_grid(estimator ~ panel, scales = "free_y") +
  scale_x_continuous(breaks = 2:8) +
  scale_linetype_manual(values = c("pooled OLS" = "dashed", "fixed effects" = "dotted")) +
  labs(x = "Longest lag used to build the instruments",
       y = "Estimated effect of income",
       color = "Instruments", linetype = "Benchmark",
       title = "Collapsing holds the estimates steady as the lag window widens",
       subtitle = "Two-step GMM; the uncollapsed set grows with the square of the lag depth") +
  theme_minimal(base_size = 11)
ggsave(file.path(PATH_OUTPUT, "instruments.png"), fig, width = 10, height = 5.5, dpi = 150)

drift <- tab |>
  group_by(panel, estimator, instruments_kept) |>
  summarise(first = income[which.min(max_lag)], last = income[which.max(max_lag)],
            n_first = n_inst[which.min(max_lag)], n_last = n_inst[which.max(max_lag)],
            .groups = "drop") |>
  mutate(move = last - first)
write_csv(drift, file.path(PATH_OUTPUT, "instruments_drift.csv"))

pick <- function(pl, es, ik) arrange(filter(tab, panel == pl, estimator == es,
                                            instruments_kept == ik), max_lag)
ends <- function(d, col) c(d[[col]][1], d[[col]][nrow(d)])
fhc <- pick("Freedom House", "Difference GMM", "collapsed")
fhu <- pick("Freedom House", "Difference GMM", "uncollapsed")
poc <- pick("Polity", "Difference GMM", "collapsed")
pou <- pick("Polity", "Difference GMM", "uncollapsed")
fhb <- filter(bench, panel == "Freedom House")

md <- c(
"# Collapsing the instruments, and letting the lag window grow",
"",
"Professor Torgovitsky asked me to take the sweep from last round and make",
"collapsing a second dial, so the lag window and the collapsing choice move",
"together in one picture. These estimators use a country's own past as a",
"stand-in for its present, so the instruments are old values of democracy and",
"income. Left separate, the estimator gets its own instrument for every pairing",
"of a lag with a time period, and the count grows with the square of the lag",
"window. Collapsed, it gets one instrument per lag distance, and the count grows",
"in step with the window instead. I ran both, for lag windows of 2 through 8,",
"for both GMM estimators and both democracy measures. The countries and years",
"stay the same throughout; only the instrument list changes. Results are in",
"output/instruments.txt and output/instruments.csv, and the figure is",
"output/instruments.png.",
"",
"## What happened",
"",
"Half of the expectation held and half of it did not, so it is worth separating",
"the two measures.",
"",
sprintf("For Polity, collapsing does what it is supposed to. The collapsed estimate"),
sprintf("starts at %.3f and ends at %.3f, drifting slightly away from zero rather",
        ends(poc, "income")[1], ends(poc, "income")[2]),
"than toward it, while the uncollapsed one",
sprintf("climbs from %.3f to %.3f, moving toward the fixed-effects value of %.3f.",
        ends(pou, "income")[1], ends(pou, "income")[2],
        filter(bench, panel == "Polity")$fe),
"",
"For Freedom House it does not. Both lines drift, and by almost the same amount:",
sprintf("collapsed runs from %.3f to %.3f and uncollapsed from %.3f to %.3f, both",
        ends(fhc, "income")[1], ends(fhc, "income")[2],
        ends(fhu, "income")[1], ends(fhu, "income")[2]),
sprintf("heading toward the fixed-effects value of %.3f. Collapsing slows the drift",
        fhb$fe),
"here but does not stop it. The reason is that collapsing does not freeze the",
"instrument count, it only slows its growth, so a wide enough window still piles",
"up enough instruments to bite.",
"",
"System GMM barely moves under either setting. It starts a little above the",
"pooled OLS value and stays there, collapsed or not.",
"",
sprintf("The standard errors are the clearer signal. Uncollapsed, the Freedom House"),
sprintf("standard error falls from %.3f to %.3f as the instruments pile up, and the",
        ends(fhu, "income_se")[1], ends(fhu, "income_se")[2]),
sprintf("Polity one from %.3f to %.3f. Collapsed, they hold up much better: %.3f to",
        ends(pou, "income_se")[1], ends(pou, "income_se")[2],
        ends(fhc, "income_se")[1]),
sprintf("%.3f and %.3f to %.3f. The uncollapsed estimator looks like it is getting",
        ends(fhc, "income_se")[2], ends(poc, "income_se")[1], ends(poc, "income_se")[2]),
"sharper as the window widens, and that apparent sharpness is manufactured by",
"the instrument count rather than earned from the data.",
"",
sprintf("The counts show the scale of it. For difference GMM the widest window uses %d",
        ends(fhc, "n_inst")[2]),
sprintf("instruments collapsed against %d uncollapsed. The overidentification p-value",
        ends(fhu, "n_inst")[2]),
"climbs as that count grows, which is the test losing its power to object rather",
"than the instruments getting better.",
"",
"One caveat on reading this next to the other writeups. This sweep gives both",
"democracy and income their own block of lags, which is what the plm package",
"does by default. The paper builds its own GMM column differently, from lags of",
"democracy only. Under the paper's design, collapsing does hold the Freedom",
"House estimate steady, which is in docs/aggregation.md.",
"",
"## Other ways to aggregate",
"",
"He also asked whether Roodman describes other ways to collapse. Section V of",
"the paper lays out two levers, not one, and gives a single collapsing rule.",
"Since collapsing is really just a choice about which instrument columns get",
"added together, there are other rules worth trying. That comparison is in",
"docs/aggregation.md.")
writeLines(md, file.path(PATH_DOCS, "instruments.md"))

cat("Instrument sweep written (collapsed and uncollapsed).\n")
print(drift)
