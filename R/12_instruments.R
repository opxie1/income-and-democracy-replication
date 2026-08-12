source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(plm))

# estimators
TRANSFORMATIONS <- c(d = "Difference GMM", ld = "System GMM")

# sweep
lag_sweep <- function(dep, inc, panel_label) {
  d <- prep_dynamic(FILE_P5, dep, inc)
  s <- filter(d, sample == 1)
  ols <- as.numeric(ce(fit_ols(s, dep, c("Ldep", "Linc"), FALSE), "Linc")["est"])
  fe  <- as.numeric(ce(fit_ols(s, dep, c("Ldep", "Linc"), TRUE),  "Linc")["est"])
  pd <- gmm_panel(s, dep, inc)
  rows <- list()
  for (wt in names(GMM_MODELS)) for (L in LAG_WINDOW)
    for (tr in names(TRANSFORMATIONS)) for (cl in c(FALSE, TRUE)) {
      m <- pgmm(gmm_formula(dep, inc, sprintf("2:%d", L)), data = pd, effect = "twoways",
                model = wt, transformation = tr, collapse = cl)
      ic <- unname(co_row(pgmm_co(m), inc))
      rows[[length(rows) + 1]] <- tibble(
        panel = panel_label,
        weighting = unname(GMM_MODELS[wt]),
        estimator = unname(TRANSFORMATIONS[tr]),
        instruments_kept = ifelse(cl, KEEP_LABEL[["lag"]], KEEP_LABEL[["none"]]),
        max_lag = L, income = ic[1], income_se = ic[2],
        obs = nobs(m), countries = pgmm_countries(m), n_inst = ncol(m$W[[1]]),
        ar1_p = tryCatch(mtest(m, 1)$p.value, error = function(e) NA_real_),
        ar2_p = tryCatch(mtest(m, 2)$p.value, error = function(e) NA_real_),
        overid_test = unname(OVERID_TEST[wt]),
        overid_df = tryCatch(unname(sargan(m, weights = wt)$parameter),
                             error = function(e) NA_real_),
        overid_p = tryCatch({
          sg <- sargan(m, weights = wt)
          if (is.na(sg$parameter) || sg$parameter < 1) NA_real_ else sg$p.value
        }, error = function(e) NA_real_))
    }
  tb <- bind_rows(rows)
  cells <- length(GMM_MODELS) * length(LAG_WINDOW) * length(TRANSFORMATIONS) * 2L
  key <- interaction(tb$estimator, tb$weighting, drop = TRUE)
  stopifnot(nrow(tb) == cells,
            all(tapply(tb$obs, key, function(x) length(unique(x))) == 1L),
            all(tapply(tb$countries, key, function(x) length(unique(x))) == 1L),
            all(tapply(tb$obs, tb$estimator, function(x) length(unique(x))) == 1L))
  list(tab = tb, bench = tibble(panel = panel_label, ols = ols, fe = fe))
}

# fits
res <- lapply(MEASURES, function(ms) lag_sweep(ms$dep, ms$inc, ms$label))
tab <- bind_rows(lapply(res, function(r) r$tab))
bench <- bind_rows(lapply(res, function(r) r$bench))
write_csv(tab, file.path(PATH_OUTPUT, "instruments.csv"))

# table
txt <- c()
for (pl in unique(tab$panel)) {
  b <- filter(bench, panel == pl)
  txt <- c(txt, sprintf("== %s (pooled OLS %.3f, fixed effects %.3f) ==", pl, b$ols, b$fe),
           "Each cell: estimate (SE) [instruments, overid p]",
           "'exact id' means there is nothing left over to test")
  for (wt in GMM_MODELS) for (es in TRANSFORMATIONS) {
    tst <- filter(tab, weighting == wt)$overid_test[1]
    txt <- c(txt, sprintf("-- %s, %s (overid p is %s's) --", wt, es, tst),
             sprintf("%-8s %-30s %-30s", "Max lag",
                     KEEP_LABEL[["none"]], KEEP_LABEL[["lag"]]))
    for (L in sort(unique(tab$max_lag))) {
      pick <- function(k) filter(tab, panel == pl, estimator == es, weighting == wt,
                                 max_lag == L, instruments_kept == k)
      cell <- function(r) if (nrow(r) == 0) "" else sprintf("%s (%s) [%d, %s]",
        num(r$income), num(r$income_se), r$n_inst,
        ifelse(is.na(r$overid_p), "exact id", num(r$overid_p, 2)))
      txt <- c(txt, sprintf("%-8d %-30s %-30s", L,
                            cell(pick(KEEP_LABEL[["none"]])), cell(pick(KEEP_LABEL[["lag"]]))))
    }
    txt <- c(txt, "")
  }
}
writeLines(txt, file.path(PATH_OUTPUT, "instruments.txt"))

# figure
benchlong <- bind_rows(
  transmute(bench, panel, kind = "pooled OLS", value = ols),
  transmute(bench, panel, kind = "fixed effects", value = fe))
fig <- ggplot(tab, aes(max_lag, income, color = instruments_kept)) +
  geom_hline(data = benchlong, aes(yintercept = value, linetype = kind), color = "grey45") +
  geom_line() + geom_point(size = 1.4) +
  facet_grid(estimator ~ panel + weighting, scales = "free_y") +
  scale_x_continuous(breaks = LAG_WINDOW) +
  scale_linetype_manual(values = c("pooled OLS" = "dashed", "fixed effects" = "dotted")) +
  labs(x = "Longest lag used to build the instruments",
       y = "Estimated effect of income",
       color = "Instruments", linetype = "Benchmark",
       title = "Collapsing holds the estimate still, but only under one-step weighting",
       subtitle = sprintf(paste("Within a row the vertical scale is shared, so the two",
                                "weighting columns are directly comparable.",
                                "\nAt the widest window collapsing cuts the instrument",
                                "count from %d to %d."),
         max(filter(tab, estimator == TRANSFORMATIONS[["d"]],
                    instruments_kept == KEEP_LABEL[["none"]])$n_inst),
         max(filter(tab, estimator == TRANSFORMATIONS[["d"]],
                    instruments_kept == KEEP_LABEL[["lag"]])$n_inst))) +
  theme_minimal(base_size = 11)
ggsave(file.path(PATH_OUTPUT, "instruments.png"), fig, width = 12, height = 6.5, dpi = 150)

# drift
drift <- tab |>
  group_by(panel, weighting, estimator, instruments_kept) |>
  arrange(max_lag, .by_group = TRUE) |>
  summarise(first = first(income), last = last(income),
            move = last(income) - first(income),
            travel = sum(abs(diff(income))),
            n_first = first(n_inst), n_last = last(n_inst), .groups = "drop")
write_csv(drift, file.path(PATH_OUTPUT, "instruments_drift.csv"))

# readers
pick <- function(pl, wt, es, ik) arrange(filter(tab, panel == pl, weighting == wt,
                                                estimator == es, instruments_kept == ik),
                                         max_lag)
ends <- function(d, col) c(d[[col]][1], d[[col]][nrow(d)])
DIFF <- TRANSFORMATIONS[["d"]]
SYS  <- TRANSFORMATIONS[["ld"]]
UNC  <- KEEP_LABEL[["none"]]
COL  <- KEEP_LABEL[["lag"]]
ONE  <- GMM_MODELS[["onestep"]]
TWO  <- GMM_MODELS[["twostep"]]

mv <- function(pl, wt, ik, col = "move")
  abs(filter(drift, panel == pl, weighting == wt, estimator == DIFF,
             instruments_kept == ik)[[col]])

fh <- MEASURES[[1]]$label
po <- MEASURES[[2]]$label
fhb <- filter(bench, panel == fh)
pob <- filter(bench, panel == po)

collapse_holds <- function(wt) {
  d <- filter(drift, weighting == wt, estimator == DIFF)
  all(abs(filter(d, instruments_kept == COL)$move) <
      abs(filter(d, instruments_kept == UNC)$move))
}
stopifnot(collapse_holds(ONE), !collapse_holds(TWO))

ratio <- function(pl, wt) mv(pl, wt, UNC) / mv(pl, wt, COL)
sys_worst <- max(filter(tab, estimator == SYS, panel == po,
                        weighting == TWO)$overid_p, na.rm = TRUE)
ar2_ok <- sum(is.finite(tab$ar2_p))
unc <- function(pl, wt)
  filter(tab, panel == pl, estimator == DIFF, weighting == wt,
         instruments_kept == UNC) |> arrange(max_lag)
one_max <- max(filter(tab, weighting == ONE)$overid_p, na.rm = TRUE)
two_rise <- function(pl) ends(unc(pl, TWO), "overid_p")
stopifnot(one_max < CI_LEVEL,
          vapply(c(fh, po), function(pl) {
            e <- two_rise(pl); e[2] > e[1] && e[2] > CI_LEVEL
          }, logical(1)))


# report
write_doc("instruments.md",
"# Collapsed instruments and a wider lag window",
sprintf(paste(
  "Professor Torgovitsky asked me to add collapsing as a second dial on the sweep",
  "from last round. The lag window and the collapsing choice now move together in one",
  "picture. These estimators use the past of a country as a stand-in for its present,",
  "so the instruments are old values of democracy and income. Without collapsing, the",
  "estimator gets one instrument for each pair of a lag and a time period. With",
  "collapsing, it gets one instrument per lag distance, so the count is much smaller.",
  "I ran both, for lag windows of %d through %d, for both GMM estimators and both",
  "democracy measures."), min(LAG_WINDOW), max(LAG_WINDOW)),
paste(
  "A third dial matters more than the lag window and the collapsing choice, so it is on",
  "the figure too. The sweep from last round used the two-step weighting matrix. Last",
  "round I did not see that weighting step as a choice. This round the weighting step",
  "is a dial on the figure. It is the choice that decides the answer. Results are in",
  "output/instruments.txt, output/instruments.csv and output/instruments_drift.csv, and",
  "the figure is output/instruments.png."),
sprintf(paste(
  "Two notes on how to read the output. Within each estimator the sample is fixed, so",
  "the rows of one block in output/instruments.txt differ only in the instrument set.",
  "Difference and system GMM do not use the same sample as each other. For %s",
  "they use %d and %d countries, because system GMM can also use the level equations.",
  "Their levels are therefore not strictly comparable, and only the trends carry",
  "information."),
  fh, filter(tab, panel == fh, estimator == DIFF)$countries[1],
  filter(tab, panel == fh, estimator == SYS)$countries[1]),
sprintf(paste(
  "The instrument count also increases quickly with the window, but not without",
  "limit. For difference GMM it runs %s across windows of %d to %d. Once the window",
  "is as long as the panel, the count is flat."),
  commas(pick(fh, ONE, DIFF, UNC)$n_inst), min(LAG_WINDOW), max(LAG_WINDOW)),

"## What happened",
paste(
  "Professor Torgovitsky expected the collapsed line to stay flat while the uncollapsed",
  "line drifted toward OLS or fixed effects. Under one-step weighting that is what",
  "happens, for both measures. Under two-step weighting it happens for one measure and",
  "not for the other. The sweep from last round used two-step weighting. So the puzzle",
  "I reported then is a fact about the weighting step and not about collapsing."),
sprintf(paste(
  "One-step difference GMM comes first, in the top row of the figure and its two",
  "one-step columns. Between the narrowest window and the widest, the collapsed",
  "estimate moves by %.3f for %s against %.3f uncollapsed. The collapsed estimate",
  "moves by %.3f for %s against %.3f. That is a factor of %.0f and %.0f. The",
  "uncollapsed lines drift toward the fixed-effects values of %.3f and %.3f. The",
  "collapsed lines do not."),
  mv(fh, ONE, COL), fh, mv(fh, ONE, UNC), mv(po, ONE, COL), po, mv(po, ONE, UNC),
  ratio(fh, ONE), ratio(po, ONE), fhb$fe, pob$fe),
sprintf(paste(
  "Two-step weighting breaks that pattern for %s. There the collapsed estimate moves",
  "by %.3f against %.3f uncollapsed, so collapsing does not slow the drift at all.",
  "The collapsed line runs from %.3f up to %.3f. The uncollapsed line runs from %.3f",
  "to %.3f. For %s, two-step weighting behaves like one-step, with a collapsed move",
  "of %.3f against %.3f uncollapsed."),
  fh, mv(fh, TWO, COL), mv(fh, TWO, UNC),
  ends(pick(fh, TWO, DIFF, COL), "income")[1], ends(pick(fh, TWO, DIFF, COL), "income")[2],
  ends(pick(fh, TWO, DIFF, UNC), "income")[1], ends(pick(fh, TWO, DIFF, UNC), "income")[2],
  po, mv(po, TWO, COL), mv(po, TWO, UNC)),
paste(
  "The reason is the two-step weighting matrix. The estimator builds that matrix from",
  "the residuals, and the matrix has one row and one column per instrument. With a",
  "wide window and a few dozen countries, the data are too thin to estimate the matrix",
  "well. The matrix then absorbs noise of its own. Collapsing makes that matrix",
  "smaller, but not small enough here."),
paste(
  "So the honest summary is that collapsing does its advertised job. The two-step",
  "weighting step is a separate form of overfitting, and collapsing does not fix it. A",
  "reader who looks only at the two-step columns will blame collapsing for a problem",
  "that belongs somewhere else."),
sprintf(paste(
  "System GMM barely moves under any setting. It starts a little above the pooled OLS",
  "value and stays there. The robust overidentification test also rejects it for %s",
  "everywhere. The largest p-value across the whole two-step sweep is %s. So those",
  "level conditions look invalid on this data, whatever the instrument count is."),
  po, num(sys_worst, 3)),

"## The standard errors",
paste(
  "The standard errors tell the same story from the other side, and they do it under",
  "both weighting choices."),
sprintf(paste(
  "In the uncollapsed one-step columns, the %s standard error decreases from %.3f to",
  "%.3f as the instrument count increases. The %s standard error decreases from %.3f to",
  "%.3f. In the collapsed one-step columns, both standard errors stay near their start.",
  "The %s standard error moves from %.3f to %.3f, and the %s standard error moves from",
  "%.3f to %.3f. The uncollapsed estimator therefore looks sharper as the window widens.",
  "The instrument count manufactures that sharpness, and the data do not earn it."),
  fh, ends(pick(fh, ONE, DIFF, UNC), "income_se")[1],
  ends(pick(fh, ONE, DIFF, UNC), "income_se")[2],
  po, ends(pick(po, ONE, DIFF, UNC), "income_se")[1],
  ends(pick(po, ONE, DIFF, UNC), "income_se")[2],
  fh, ends(pick(fh, ONE, DIFF, COL), "income_se")[1],
  ends(pick(fh, ONE, DIFF, COL), "income_se")[2],
  po, ends(pick(po, ONE, DIFF, COL), "income_se")[1],
  ends(pick(po, ONE, DIFF, COL), "income_se")[2]),
paste(
  "The weak-instrument work in docs/weak-instruments.md gives the same reading for the",
  "two-stage least squares columns. A conventional standard error assumes the",
  "instruments have real pull, and it decreases whether or not they do."),
sprintf(paste(
  "The instrument counts show the scale of this effect. For difference GMM, the widest",
  "window uses %d instruments in the collapsed set and %d instruments in the uncollapsed",
  "set."),
  ends(pick(fh, ONE, DIFF, COL), "n_inst")[2], ends(pick(fh, ONE, DIFF, UNC), "n_inst")[2]),

"## The overidentification test",
paste(
  "One caution comes before the numbers, and it took me a while to get this right. The",
  "plm package reports an overidentification test after either estimator. By default",
  "it builds the statistic from the robust weight matrix for both estimators. So the",
  "number it prints after a one-step fit is not the criterion of that estimator. I",
  "report the criterion of each estimator instead."),
sprintf(paste(
  "The %s statistic sits in the one-step rows, and it assumes homoskedastic errors.",
  "The %s statistic sits in the two-step rows. It is robust to heteroskedasticity, but",
  "it loses power as the instrument count increases. They are different tests, so a",
  "reader must compare p-values down a column and not across the weighting rows."),
  OVERID_TEST[["onestep"]], OVERID_TEST[["twostep"]]),
sprintf(paste(
  "The one-step rows say nothing useful, and the reason matters more than the numbers.",
  "The %s statistic rejects every specification in the sweep, and the largest p-value",
  "anywhere in that column is %s. The test rejects the collapsed sets and the",
  "uncollapsed sets, at every window, for both measures. That is what a",
  "homoskedasticity-assuming test does on data with heteroskedastic errors. The errors",
  "here are heteroskedastic enough that the paper clusters all of its own inference by",
  "country."),
  OVERID_TEST[["onestep"]], num(one_max, 3)),
paste(
  "The spread within that column is enormous and covers many orders of magnitude. The",
  "test breaks its own assumption everywhere, so I do not think that spread carries",
  "any weight."),
sprintf(paste(
  "The two-step rows carry the robust %s statistic, and that is the column to read.",
  "There the uncollapsed p-value increases from %s to %s for %s, and from %s to %s for",
  "%s. At the same time the instrument count increases from %d to %d. That is one",
  "statistic applied to a larger and larger instrument set, so the comparison is",
  "clean. The test loses power, and the instruments do not improve. So a high %s",
  "p-value on a large instrument set is not reassurance."),
  OVERID_TEST[["twostep"]],
  num(two_rise(fh)[1], 2), num(two_rise(fh)[2], 2), fh,
  num(two_rise(po)[1], 2), num(two_rise(po)[2], 2), po,
  ends(unc(fh, TWO), "n_inst")[1], ends(unc(fh, TWO), "n_inst")[2],
  OVERID_TEST[["twostep"]]),
sprintf(paste(
  "The collapsed sets behave differently. Their p-values stay high for most of the",
  "window, and then drop sharply at the widest window (%s for %s). I do not have a",
  "clean account of that last drop, and it is the weakest link in this section."),
  commas(pick(fh, TWO, DIFF, COL)$overid_p[-1], 2), fh),
sprintf(paste(
  "The second-order serial correlation test matters most for whether these instruments",
  "are valid at all. It passes everywhere. The smallest p-value across all %d fits is",
  "%.2f, so nothing here suggests that the lag-2 instruments are invalid."),
  ar2_ok, min(tab$ar2_p, na.rm = TRUE)),

"## How this fits with the other sweep",
sprintf(paste(
  "docs/aggregation.md runs the same idea through the one-step estimator written for",
  "the replication, and it reaches a tidier conclusion. Last round I did not know",
  "which of two differences explained the gap, the weighting step or the sample. The",
  "weighting dial answers that question, and the answer is the weighting step. The",
  "sample and the instrument design stay fixed, and only the weighting step changes.",
  "That one change turns the %s collapsed line from a %.3f drift into a %.3f one. The",
  "samples still differ between the two sweeps, but I no longer need that difference",
  "to explain anything."),
  fh, mv(fh, TWO, COL), mv(fh, ONE, COL)),

"## Other ways to collapse",
paste(
  "Professor Torgovitsky also asked whether Roodman describes other ways to collapse.",
  "Section V of the paper gives two levers rather than one, and one collapsing rule.",
  "Collapsing is a choice about which instrument columns the estimator adds together, so",
  "other rules exist and deserve a test. That comparison is in docs/aggregation.md."))

cat("Instrument sweep written (collapsed and uncollapsed, one-step and two-step).\n")
print(filter(drift, estimator == DIFF), n = 20)
