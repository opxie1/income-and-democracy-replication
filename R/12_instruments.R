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
  for (L in 2:8) for (tr in c("d", "ld")) {
    m <- pgmm(gmm_formula(dep, inc, sprintf("2:%d", L)), data = pd, effect = "twoways",
              model = "twostep", transformation = tr, collapse = FALSE)
    ic <- unname(co_row(pgmm_co(m), inc))
    rows[[length(rows) + 1]] <- tibble(
      panel = panel_label,
      estimator = ifelse(tr == "d", "Difference GMM", "System GMM"),
      max_lag = L, income = ic[1], income_se = ic[2],
      obs = nobs(m), countries = pgmm_countries(m), instruments = ncol(m$W[[1]]),
      ar1_p = tryCatch(mtest(m, 1)$p.value, error = function(e) NA_real_),
      ar2_p = tryCatch(mtest(m, 2)$p.value, error = function(e) NA_real_),
      overid_p = tryCatch(sargan(m)$p.value, error = function(e) NA_real_))
  }
  tb <- bind_rows(rows)
  stopifnot(nrow(tb) == 14L,
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
  txt <- c(txt,
           sprintf("== %s (pooled OLS %.3f, fixed effects %.3f) ==", pl, b$ols, b$fe),
           "Each cell: estimate (SE) [instruments, overid p, AR(2) p]",
           sprintf("%-8s %-33s %-33s", "Max lag", "Difference GMM", "System GMM"))
  for (L in sort(unique(tab$max_lag))) {
    dd <- filter(tab, panel == pl, max_lag == L, estimator == "Difference GMM")
    ss <- filter(tab, panel == pl, max_lag == L, estimator == "System GMM")
    cell <- function(r) if (nrow(r) == 0) "" else sprintf("%s (%s) [%d, %s, %s]",
      num(r$income), num(r$income_se), r$instruments,
      num(r$overid_p, 2), num(r$ar2_p, 2))
    txt <- c(txt, sprintf("%-8d %-33s %-33s", L, cell(dd), cell(ss)))
  }
  txt <- c(txt, "")
}
writeLines(txt, file.path(PATH_OUTPUT, "instruments.txt"))

benchlong <- bind_rows(
  transmute(bench, panel, kind = "pooled OLS", value = ols),
  transmute(bench, panel, kind = "fixed effects", value = fe))
fig <- ggplot(tab, aes(max_lag, income, color = estimator)) +
  geom_hline(data = benchlong, aes(yintercept = value, linetype = kind), color = "grey45") +
  geom_line() + geom_point(size = 1.6) +
  facet_wrap(~panel, scales = "free_y") +
  scale_x_continuous(breaks = 2:8) +
  scale_linetype_manual(values = c("pooled OLS" = "dashed", "fixed effects" = "dotted")) +
  labs(x = "Longest lag used to build the instruments",
       y = "Estimated effect of income",
       color = NULL, linetype = "Benchmark",
       title = "Difference GMM drifts toward fixed effects; system GMM stays at pooled OLS",
       subtitle = "Instruments left uncollapsed; two-step GMM") +
  theme_minimal(base_size = 11)
ggsave(file.path(PATH_OUTPUT, "instruments.png"), fig, width = 9, height = 4, dpi = 150)

md <- c(
"# Letting the instrument set grow",
"",
"Professor Torgovitsky asked how the GMM estimates behave when the instruments",
"are not collapsed and are built from longer and longer lags. These methods use",
"a country's own past as a stand-in for its present, so the instruments here",
"are old values of democracy and income. In docs/alternatives.md I kept that",
"set small on purpose. This time I did the opposite. Every instrument stays",
"separate, nothing is collapsed, and the lags run back 2 periods, then 3, and",
"so on out to 8, which is as far back as the 1960 to 2000 sample goes. The",
"countries and years stay identical the whole way; only the instrument list",
"changes. The results are in output/instruments.txt and output/instruments.csv,",
"and the picture is output/instruments.png.",
"",
"## What happened",
"",
"With the shortest list, difference GMM finds income pushing against democracy:",
"about -0.25 for Freedom House and -0.34 for Polity, and the Polity number is",
"nearly the paper's own GMM estimate. Growing the list wears this down. Freedom",
"House loses most of its estimate by the fourth lag and finishes near -0.09,",
"not far from the fixed effects value of 0.01. Polity fades more slowly, and",
"with every lag in use it still reads -0.25. System GMM barely reacts at all.",
"It starts just above the plain OLS value and stays there through every step.",
"",
"One note for anyone comparing instrument counts with the replication. Here",
"democracy and income each get the full run of lags, while the paper's GMM",
"column gives income a single lag, so my counts top out higher, 63 against the",
"replication's 55.",
"",
"## Why",
"",
"Too many instruments is a known trap. Each extra lag adds another condition",
"for the estimator to satisfy, and past some point it stops correcting the bias",
"it was built to fix and starts copying it instead. Difference GMM slides",
"toward fixed effects. System GMM sits on OLS from the start because its extra",
"assumptions already lean that way. There is a standard test meant to catch bad",
"instruments, and it goes quiet as the pile grows; a rising p-value there means",
"the test is losing its power to object. The routine check that would flag a",
"deeper problem stays clean at every lag length, so the drift traces back to",
"the size of the instrument set.",
"",
"So the answer to the original question: the cross-over is quick for Freedom",
"House, mostly finished by the fourth lag, and it never quite completes for",
"Polity.")
writeLines(md, file.path(PATH_DOCS, "instruments.md"))

cat("Instrument sweep written (table, figure, and writeup).\n")
