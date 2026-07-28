source(here::here("R", "00_setup.R"))

SCHEME_LABEL <- c(
  none = "One per lag and period (uncollapsed)",
  lag = "One per lag distance (Roodman)",
  period = "One per period, summed",
  period_mean = "One per period, averaged",
  full = "One instrument in total")

agg_sweep <- function(dep, inc, panel_label) {
  d <- add_lags(read_panel(FILE_P5), c(dep, inc, "year"), 1:2)
  d <- d |> mutate(
    Ldep = .data[[paste0(dep, "_l1")]], Linc = .data[[paste0(inc, "_l1")]],
    y = .data[[dep]] - .data[[paste0(dep, "_l1")]],
    dLdep = .data[[paste0(dep, "_l1")]] - .data[[paste0(dep, "_l2")]],
    dLinc = .data[[paste0(inc, "_l1")]] - .data[[paste0(inc, "_l2")]],
    L2inc = .data[[paste0(inc, "_l2")]])
  s <- filter(d, sample == 1)
  est <- complete_on(s, c("y", "dLdep", "dLinc", "L2inc"))
  fe <- as.numeric(ce(fit_ols(s, dep, c("Ldep", "Linc"), TRUE), "Linc")["est"])

  designs <- list(
    list(name = "Paper's design (lags of democracy only)", levels = dep, extra = "L2inc"),
    list(name = "Symmetric (lags of democracy and income)",
         levels = c(dep, inc), extra = character()))

  rows <- list()
  for (dg in designs) for (sc in AGG_SCHEMES) for (L in c(2:8, Inf)) {
    m <- tryCatch(fit_abgmm(d, est, dep_level = dg$levels, endog = c("dLdep", "dLinc"),
                            inst_extra = dg$extra, lag_max = L, scheme = sc),
                  error = function(e) NULL)
    if (is.null(m)) next
    rows[[length(rows) + 1]] <- tibble(
      panel = panel_label, design = dg$name,
      scheme = sc, scheme_label = unname(SCHEME_LABEL[sc]),
      lag_max = L, income = unname(m$coef["dLinc"]), income_se = unname(m$se["dLinc"]),
      dem = unname(m$coef["dLdep"]), n_inst = m$n_inst,
      obs = m$nobs, countries = m$n_country)
  }
  tb <- bind_rows(rows)
  stopifnot(length(unique(tb$obs)) == 1L, length(unique(tb$countries)) == 1L)
  list(tab = tb, bench = tibble(panel = panel_label, fe = fe))
}

res <- lapply(MEASURES, function(ms) agg_sweep(ms$dep, ms$inc, ms$label))
tab <- bind_rows(lapply(res, function(r) r$tab))
bench <- bind_rows(lapply(res, function(r) r$bench))

paper <- filter(tab, scheme == "none", is.infinite(lag_max), grepl("^Paper", design))
pub <- read_csv(file.path(PATH_DOCS, "published_values.csv"), show_col_types = FALSE) |>
  filter(table %in% c("2", "3"), column == 4, row == "inc")
stopifnot(nrow(paper) == 2L,
          all(abs(sort(round(paper$income, 3)) - sort(pub$value)) < 1e-9))
cat("Uncollapsed, all lags reproduces the published GMM column for both measures.\n")

write_csv(tab, file.path(PATH_OUTPUT, "aggregation.csv"))

txt <- c()
for (pl in unique(tab$panel)) {
  txt <- c(txt, sprintf("== %s (fixed effects %.3f) ==", pl, filter(bench, panel == pl)$fe),
           "Each cell: estimate (SE) [instruments]")
  for (dg in unique(tab$design)) {
    txt <- c(txt, sprintf("-- %s --", dg),
             sprintf("%-38s %-24s %-24s", "Aggregation of the lagged levels",
                     "Lags 2-4", "All lags"))
    for (sc in AGG_SCHEMES) {
      cell <- function(L) {
        r <- filter(tab, panel == pl, design == dg, scheme == sc,
                    if (is.infinite(L)) is.infinite(lag_max) else lag_max == L)
        if (nrow(r) == 0) "" else sprintf("%s (%s) [%d]",
                                          num(r$income), num(r$income_se), r$n_inst)
      }
      txt <- c(txt, sprintf("%-38s %-24s %-24s", SCHEME_LABEL[sc], cell(4), cell(Inf)))
    }
    txt <- c(txt, "")
  }
}
writeLines(txt, file.path(PATH_OUTPUT, "aggregation.txt"))

plotdat <- filter(tab, is.finite(lag_max))
fig <- ggplot(plotdat, aes(lag_max, income, color = scheme_label)) +
  geom_hline(data = bench, aes(yintercept = fe), color = "grey45", linetype = "dotted") +
  geom_line() + geom_point(size = 1.4) +
  facet_grid(design ~ panel, scales = "free_y") +
  scale_x_continuous(breaks = 2:8) +
  labs(x = "Longest lag used to build the instruments",
       y = "Estimated effect of income", color = "Instruments built as",
       title = "How the lagged levels are pooled into instruments changes the answer",
       subtitle = "Difference GMM, one-step; dotted line is the fixed-effects estimate") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")
ggsave(file.path(PATH_OUTPUT, "aggregation.png"), fig, width = 11, height = 6, dpi = 150)

at_full <- filter(tab, is.infinite(lag_max))
rng <- at_full |> group_by(panel, design) |>
  summarise(lo = min(income), hi = max(income), .groups = "drop")
papr <- filter(at_full, grepl("^Paper", design))
symm <- filter(at_full, grepl("^Symmetric", design))
worst_se <- at_full[which.max(at_full$income_se), ]

md <- c(
"# Other ways to pool the instruments",
"",
"Professor Torgovitsky asked whether Roodman describes other ways to collapse",
"the instruments, since there might be several ways to aggregate them. This is",
"what I found in the paper and what happens when the alternatives are tried on",
"this data. The numbers are in output/aggregation.txt and output/aggregation.csv,",
"and the figure is output/aggregation.png.",
"",
"## What Roodman actually describes",
"",
"Section V of the paper is the relevant part, and it lays out two levers rather",
"than one. The first is capping how far back the lags go. The second is",
"collapsing, which he describes as squeezing the instrument matrix sideways and",
"adding together the columns that end up on top of each other. He gives one",
"collapsing rule: make a single instrument for each lag distance. He also points",
"out that pulling both levers at once leaves a count that does not grow with the",
"length of the panel at all, and his Table 1 crosses the two levers the same way",
"the figure in docs/instruments.md does.",
"",
"Two further ideas appear in his footnotes rather than the main text. One is an",
"approach from Arellano that first models the instruments jointly with a vector",
"autoregression and uses the fitted coefficients as constraints, which Roodman",
"notes had not caught on. The other is drawing random subsets of the available",
"instruments repeatedly to see how much the answer moves around.",
"",
"## The alternatives I tried",
"",
"Collapsing is really just a decision about which instrument columns get added",
"together, so there is more than one way to do it. Roodman adds across time and",
"keeps the lag distances apart. The obvious alternative is to do it the other",
"way round, adding across lag distances and keeping the periods apart, which",
"gives one instrument per period. Adding and averaging are different here,",
"because the number of lags available grows over time, so averaging is a third",
"option. Pushing it to the end gives a fourth, where everything collapses into a",
"single instrument. I ran all four against the uncollapsed set.",
"",
"I also ran two instrument designs, because they turn out to matter. The paper",
"builds its lag blocks from democracy only and gives income a single lagged",
"level. The symmetric alternative gives both variables the full block of lags,",
"which is what the plm package does by default and what the sweep in",
"docs/instruments.md uses. The first design reproduces the published GMM column",
"exactly, and the code stops if it ever fails to.",
"",
"## What came out",
"",
"The pooling rule matters, and by more than I expected. Using every lag under",
sprintf("the paper's design, the estimates for Freedom House range from %.3f to %.3f",
        filter(rng, panel == "Freedom House", grepl("^Paper", design))$lo,
        filter(rng, panel == "Freedom House", grepl("^Paper", design))$hi),
sprintf("depending only on how the instruments are pooled, and for Polity from %.3f",
        filter(rng, panel == "Polity", grepl("^Paper", design))$lo),
sprintf("to %.3f. Same data, same lags, same estimator.",
        filter(rng, panel == "Polity", grepl("^Paper", design))$hi),
"",
"The uncollapsed set is the odd one out. In the figure it is the only line that",
"climbs steadily toward the fixed-effects value as the lag window widens; the",
"pooled versions stay far flatter. That is the same overfitting story as",
"before, and it shows up more sharply here than in the sweep in",
"docs/instruments.md because that sweep uses the symmetric design, where pooling",
"alone does not settle the Freedom House estimate down.",
"",
"Among the pooled rules, Roodman's is the best behaved. The two period-based",
"rules give noticeably more negative estimates, and collapsing everything into a",
"single instrument is the least stable of the lot: with the symmetric design it",
sprintf("produces a standard error of %.2f for %s, which is another way of saying",
        worst_se$income_se, worst_se$panel),
"the instrument has stopped carrying usable information. There is a middle",
"ground. Pool too little and the estimator overfits, pool too much and there is",
"not enough left to identify anything. That second failure is the weak",
"instrument problem, which is what docs/weak-instruments.md takes up.")
writeLines(md, file.path(PATH_DOCS, "aggregation.md"))

cat("Aggregation comparison written.\n")
print(filter(tab, is.infinite(lag_max)) |>
        transmute(panel, design = substr(design, 1, 9), scheme,
                  income = round(income, 3), se = round(income_se, 3), n_inst),
      n = 40)
