source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(plm))

lag_sweep <- function(dep, inc, panel_label) {
  d <- prep_dynamic(FILE_P5, dep, inc)
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
      overid_df = tryCatch(unname(sargan(m)$parameter), error = function(e) NA_real_),
      overid_p = tryCatch({
        sg <- sargan(m)
        if (is.na(sg$parameter) || sg$parameter < 1) NA_real_ else sg$p.value
      }, error = function(e) NA_real_))
  }
  tb <- bind_rows(rows)
  stopifnot(nrow(tb) == 28L,
            all(tapply(tb$obs, tb$estimator, function(x) length(unique(x))) == 1L),
            all(tapply(tb$countries, tb$estimator, function(x) length(unique(x))) == 1L))
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
           "Each cell: estimate (SE) [instruments, overid p]",
           "'exact id' means there is nothing left over to test")
  for (es in unique(tab$estimator)) {
    txt <- c(txt, sprintf("-- %s --", es),
             sprintf("%-8s %-30s %-30s", "Max lag", "Uncollapsed", "Collapsed"))
    for (L in sort(unique(tab$max_lag))) {
      pick <- function(k) filter(tab, panel == pl, estimator == es,
                                 max_lag == L, instruments_kept == k)
      cell <- function(r) if (nrow(r) == 0) "" else sprintf("%s (%s) [%d, %s]",
        num(r$income), num(r$income_se), r$n_inst,
        ifelse(is.na(r$overid_p), "exact id", num(r$overid_p, 2)))
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
       title = "Collapsing steadies Polity, but not Freedom House",
       subtitle = "Two-step GMM; collapsing cuts the instrument count by roughly two thirds") +
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
"of a lag with a time period. Collapsed, it gets one instrument per lag",
"distance, so the count stays much smaller. I ran both, for lag windows of 2",
"through 8, for both GMM estimators and both democracy measures. Results are in",
"output/instruments.txt and output/instruments.csv, and the figure is",
"output/instruments.png.",
"",
"Two notes on reading it. Within each estimator the sample is fixed, so only the",
"instrument list changes down a column, but difference and system GMM do not use",
sprintf("the same sample as each other: for Freedom House they use %d and %d countries",
        filter(tab, panel == "Freedom House", estimator == "Difference GMM")$countries[1],
        filter(tab, panel == "Freedom House", estimator == "System GMM")$countries[1]),
"respectively, because system GMM can also use the level equations. Their levels",
"are therefore not strictly comparable, only their trends. Also, the instrument",
"count grows quickly with the window but not without limit. For difference GMM",
sprintf("it runs %s across windows of 2 to 8,",
        paste(arrange(filter(tab, panel == "Freedom House", estimator == "Difference GMM",
                             instruments_kept == "uncollapsed"), max_lag)$n_inst,
              collapse = ", ")),
"flattening out once the window reaches back as far as the panel goes.",
"",
"## What happened",
"",
"Half of the expectation held and half of it did not, so it is worth separating",
"the two measures.",
"",
"For Polity, collapsing does what it is supposed to. The collapsed estimate",
sprintf("starts at %.3f and ends at %.3f, drifting slightly away from zero rather",
        ends(poc, "income")[1], ends(poc, "income")[2]),
sprintf("than toward it, while the uncollapsed one climbs from %.3f to %.3f, moving",
        ends(pou, "income")[1], ends(pou, "income")[2]),
sprintf("toward the fixed-effects value of %.3f.", filter(bench, panel == "Polity")$fe),
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
"The standard errors are the clearer signal. Uncollapsed, the Freedom House",
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
sprintf("instruments collapsed against %d uncollapsed.", ends(fhu, "n_inst")[2]),
"",
"The overidentification test broadly loses its bite as the count grows, which is",
"the test running out of power rather than the instruments improving, but it is",
"not a clean trend and I do not want to oversell it. Uncollapsed for Freedom",
sprintf("House it goes %s.",
        paste(num(arrange(filter(tab, panel == "Freedom House", estimator == "Difference GMM",
                                 instruments_kept == "uncollapsed"), max_lag)$overid_p, 2),
              collapse = ", ")),
"Collapsed, the narrowest window is exactly identified and leaves nothing to",
sprintf("test; from there the p-value rises most of the way, %s,",
        paste(num(arrange(filter(tab, panel == "Freedom House", estimator == "Difference GMM",
                                 instruments_kept == "collapsed",
                                 max_lag >= 3, max_lag <= 7), max_lag)$overid_p, 2),
              collapse = ", ")),
sprintf("before dropping back sharply to %s at the widest window, which is also where",
        num(filter(tab, panel == "Freedom House", estimator == "Difference GMM",
                   instruments_kept == "collapsed", max_lag == 8)$overid_p, 2)),
"that estimate makes its largest jump. I do not have a clean account of that",
"last point, and it is the weakest link in the Freedom House story above.",
"",
"The test that matters most for whether these instruments are allowed at all is",
"the second-order serial correlation check, and it passes everywhere: the",
sprintf("smallest p-value across all %d fits is %.2f, so nothing here suggests the lag-2",
        nrow(tab), min(tab$ar2_p, na.rm = TRUE)),
"instruments are invalid.",
"",
"One caveat on reading this next to docs/aggregation.md, which runs the same",
"idea through the hand-built estimator and reaches a tidier conclusion. The two",
"differ in three ways at once: this sweep uses two-step GMM through the plm",
"package and that one uses the one-step estimator written for the replication,",
"the samples differ (this one keeps fewer countries), and the levels differ",
"noticeably as a result. I checked whether the instrument design was the",
"explanation and it is not, since the symmetric design stays flat there too. So",
"the honest summary is that the Freedom House drift under collapsing shows up",
"with the two-step plm estimator on its sample and not with the one-step",
"estimator on the paper's sample, and I have not isolated which of those two",
"differences is doing the work.",
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
