source(here::here("R", "00_setup.R"))

agg_sweep <- function(dep, inc, panel_label) {
  d <- prep_dynamic(FILE_P5, dep, inc)
  s <- filter(d, sample == 1)
  est <- complete_on(s, c("y", "dLdep", "dLinc", "L2inc"))
  fe <- as.numeric(ce(fit_ols(s, dep, c("Ldep", "Linc"), TRUE), "Linc")["est"])

  designs <- list(
    list(name = "Paper's design (lags of democracy only)", short = "paper",
         levels = dep, extra = "L2inc"),
    list(name = "Symmetric (lags of democracy and income)", short = "symmetric",
         levels = c(dep, inc), extra = character()))

  rows <- list()
  for (dg in designs) for (sc in AGG_SCHEMES) for (L in c(2:8, Inf)) {
    m <- tryCatch(fit_abgmm(d, est, dep_level = dg$levels, endog = c("dLdep", "dLinc"),
                            inst_extra = dg$extra, lag_max = L, scheme = sc),
                  error = function(e) NULL)
    if (is.null(m)) next
    rows[[length(rows) + 1]] <- tibble(
      panel = panel_label, design = dg$name, design_short = dg$short,
      scheme = sc, scheme_label = unname(SCHEME_LABEL[sc]),
      lag_max = L, income = unname(m$coef["dLinc"]), income_se = unname(m$se["dLinc"]),
      dem = unname(m$coef["dLdep"]), n_gmm = m$n_gmm, n_inst = m$n_inst,
      obs = m$nobs, countries = m$n_country)
  }
  tb <- bind_rows(rows)
  stopifnot(nrow(tb) == length(designs) * length(AGG_SCHEMES) * 8L,
            length(unique(tb$obs)) == 1L, length(unique(tb$countries)) == 1L)
  list(tab = tb, bench = tibble(panel = panel_label, fe = fe))
}

family_sweep <- function(dep, inc, panel_label) {
  d <- prep_dynamic(FILE_P5, dep, inc)
  s <- filter(d, sample == 1)
  est <- complete_on(s, c("y", "dLdep", "dLinc", "L2inc"))
  nper <- n_distinct(d$period)
  fit <- function(...) fit_abgmm(d, est, dep_level = dep, endog = c("dLdep", "dLinc"),
                                 inst_extra = "L2inc", ...)
  grab <- function(m) c(unname(m$coef["dLinc"]), unname(m$se["dLinc"]), m$n_gmm)

  ref_none <- grab(fit(scheme = "none"))
  ref_lag  <- grab(fit(scheme = "lag"))
  ref_per  <- grab(fit(scheme = "period"))
  stopifnot(
    max(abs(grab(fit(scheme = "block", block_size = 1L)) - ref_none)) < 1e-9,
    max(abs(grab(fit(scheme = "block", block_size = nper)) - ref_lag)) < 1e-9,
    max(abs(grab(fit(scheme = "period_geom", rho = 1)) - ref_per)) < 1e-9)

  rows <- list()
  for (b in 1:nper) {
    g <- grab(fit(scheme = "block", block_size = b))
    rows[[length(rows) + 1]] <- tibble(panel = panel_label, family = "Blocking the years",
      setting = b, label = sprintf("%d", b),
      income = g[1], income_se = g[2], n_gmm = as.integer(g[3]))
  }
  for (r in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) {
    g <- grab(fit(scheme = "period_geom", rho = r))
    rows[[length(rows) + 1]] <- tibble(panel = panel_label,
      family = "Fading out older lags", setting = r, label = sprintf("%.1f", r),
      income = g[1], income_se = g[2], n_gmm = as.integer(g[3]))
  }
  bind_rows(rows)
}

res <- lapply(MEASURES, function(ms) agg_sweep(ms$dep, ms$inc, ms$label))
tab <- bind_rows(lapply(res, function(r) r$tab))
bench <- bind_rows(lapply(res, function(r) r$bench))
fam <- bind_rows(lapply(MEASURES, function(ms) family_sweep(ms$dep, ms$inc, ms$label)))
cat("Both families reproduce the fixed rules at their endpoints.\n")
write_csv(fam, file.path(PATH_OUTPUT, "aggregation_families.csv"))

paper <- filter(tab, scheme == "none", is.infinite(lag_max), design_short == "paper")
pub <- read_csv(file.path(PATH_DOCS, "published_values.csv"), show_col_types = FALSE) |>
  filter(table %in% c("2", "3"), column == 4, row == "inc", type == "coef") |>
  mutate(panel = ifelse(table == "2", "Freedom House", "Polity"))
chk <- left_join(paper, select(pub, panel, published = value), by = "panel")
stopifnot(nrow(chk) == 2L, !anyNA(chk$published),
          all(abs(round(chk$income, 3) - chk$published) < 1e-9))
cat("Uncollapsed, all lags reproduces the published GMM column for both measures.\n")

write_csv(tab, file.path(PATH_OUTPUT, "aggregation.csv"))

txt <- c()
for (pl in unique(tab$panel)) {
  txt <- c(txt, sprintf("== %s (fixed effects %.3f) ==", pl, filter(bench, panel == pl)$fe),
           "Each cell: estimate (SE) [lagged-level instruments, before year dummies]")
  for (dg in unique(tab$design)) {
    txt <- c(txt, sprintf("-- %s --", dg),
             sprintf("%-38s %-24s %-24s", "Aggregation of the lagged levels",
                     "Lags 2-4", "All lags"))
    for (sc in AGG_SCHEMES) {
      cell <- function(L) {
        r <- filter(tab, panel == pl, design == dg, scheme == sc,
                    if (is.infinite(L)) is.infinite(lag_max) else lag_max == L)
        if (nrow(r) == 0) "" else sprintf("%s (%s) [%d]",
                                          num(r$income), num(r$income_se), r$n_gmm)
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

famfig <- ggplot(fam, aes(setting, income)) +
  geom_hline(data = bench, aes(yintercept = fe), color = "grey45", linetype = "dotted") +
  geom_line(color = "steelblue") + geom_point(color = "steelblue", size = 1.6) +
  facet_grid(panel ~ family, scales = "free") +
  labs(x = "Block size (left panel) or weight on each older lag (right panel)",
       y = "Estimated effect of income",
       title = "Turning the pooling up and down by degrees",
       subtitle = paste("Left: years grouped into blocks of this size; 1 is uncollapsed and",
                        "the largest is Roodman's rule.\nRight: each older lag multiplied",
                        "by this factor before adding; 1 is the plain sum.")) +
  theme_minimal(base_size = 11)
ggsave(file.path(PATH_OUTPUT, "aggregation_families.png"), famfig,
       width = 10, height = 5.5, dpi = 150)

famtxt <- c("Two families that turn pooling up and down by degrees.",
            "Blocking: years grouped into blocks of the given size (1 = uncollapsed,",
            "  largest = Roodman's rule). Fading: each older lag multiplied by this",
            "  factor before adding (1 = plain sum over lags within a year).",
            "Each cell: estimate (SE) [lagged-level instruments]")
for (pl in unique(fam$panel)) {
  famtxt <- c(famtxt, "", sprintf("== %s ==", pl))
  for (fm in unique(fam$family)) {
    r <- filter(fam, panel == pl, family == fm) |> arrange(setting)
    famtxt <- c(famtxt, sprintf("-- %s --", fm),
                paste(sprintf("%-8s", r$label), collapse = ""),
                paste(sprintf("%-8s", num(r$income)), collapse = ""),
                paste(sprintf("%-8s", paste0("[", r$n_gmm, "]")), collapse = ""))
  }
}
writeLines(famtxt, file.path(PATH_OUTPUT, "aggregation_families.txt"))

at_full <- filter(tab, is.infinite(lag_max))
rng <- at_full |> group_by(panel, design_short) |>
  summarise(lo = min(income), hi = max(income), .groups = "drop")
mv <- filter(tab, is.finite(lag_max)) |> group_by(panel, design_short, scheme) |>
  summarise(rng = max(income) - min(income), .groups = "drop")
worst_se <- tab[which.max(tab$income_se), ]

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
"The paper is Roodman (2009), \"A Note on the Theme of Too Many Instruments\",",
"Oxford Bulletin of Economics and Statistics 71(1), 135-158. Full details for",
"everything cited here are in docs/references.md.",
"",
"Section V is the relevant part, and it lays out two levers rather",
"than one. The first is capping how far back the lags go. The second is",
"collapsing, which he describes as squeezing the instrument matrix sideways and",
"adding together the columns that end up on top of each other. He gives one",
"collapsing rule: make a single instrument for each lag distance. He also points",
"out that pulling both levers at once leaves a count that does not grow with the",
"length of the panel at all, and his Table 1 crosses the two levers the same way",
"the figure in docs/instruments.md does.",
"",
"Two further ideas appear in his footnotes rather than the main text. One is an",
"approach from Arellano (2003) that first models the instruments jointly with a",
"vector autoregression and uses the fitted coefficients as constraints, which",
"Roodman notes had not caught on. The other is drawing random subsets of the",
"available instruments repeatedly to see how much the answer moves around.",
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
"single instrument. I ran all four against the uncollapsed set, and then two",
"further families with a knob on them, which are in the last section.",
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
        filter(rng, panel == "Freedom House", design_short == "paper")$lo,
        filter(rng, panel == "Freedom House", design_short == "paper")$hi),
sprintf("depending only on how the instruments are pooled, and for Polity from %.3f",
        filter(rng, panel == "Polity", design_short == "paper")$lo),
sprintf("to %.3f. Same data, same lags, same estimator.",
        filter(rng, panel == "Polity", design_short == "paper")$hi),
"",
"Under the paper's own instrument design, shown in the top row of the figure,",
"the uncollapsed set is the odd one out. It is the only line that climbs",
"steadily toward the fixed-effects value as the lag window widens (a range of",
sprintf("%.3f for Freedom House and %.3f for Polity), while the four pooled rules",
        filter(mv, design_short == "paper", panel == "Freedom House", scheme == "none")$rng,
        filter(mv, design_short == "paper", panel == "Polity", scheme == "none")$rng),
sprintf("all stay inside %.3f. That is the overfitting story again.",
        max(filter(mv, design_short == "paper", scheme != "none")$rng)),
"",
"The bottom row, where income also gets its own block of lags, is messier and",
"worth saying plainly rather than glossing. There the pooled rules are not",
"uniformly steadier. Roodman's rule is still the flattest line in every panel of",
sprintf("the figure, with a range of at most %.3f, but pooling by period swings more",
        max(filter(mv, scheme == "lag")$rng)),
sprintf("than the uncollapsed set does (%.3f against %.3f for Polity), and pooling",
        filter(mv, design_short == "symmetric", panel == "Polity", scheme == "period")$rng,
        filter(mv, design_short == "symmetric", panel == "Polity", scheme == "none")$rng),
"everything together is the least stable line on the whole figure. So the useful",
"statement is not that pooling always steadies the estimate. It is that pooling",
"by lag distance, which is the rule Roodman actually proposes, is the one that",
"holds up under both designs.",
"",
"Pooling too hard has its own failure mode. With the symmetric design, pooling",
sprintf("everything together produces a standard error of %.2f for %s, which is",
        worst_se$income_se, worst_se$panel),
"another way of saying the instrument has stopped carrying usable information.",
"There is a middle ground. Pool too little and the estimator overfits, pool too",
"much and there is not enough left to identify anything. That second failure is",
"the weak instrument problem, which is what docs/weak-instruments.md takes up.",
"",
"## Turning the dial by degrees",
"",
"The five rules above are all-or-nothing. Two more ways are not single rules but",
"families with a knob on them, which lets you ask what happens part of the way.",
"They are in output/aggregation_families.txt and output/aggregation_families.png.",
"",
"The first groups the years into blocks and pools within a block, so a block size",
"of one is the uncollapsed set and a block as wide as the panel is Roodman's",
"rule. The second keeps one instrument per year but multiplies each older lag by",
"a fading factor before adding, so a factor of one is the plain sum. Both",
"families have to land exactly on the fixed rules at their endpoints, and the",
"code stops if they do not.",
"",
"Neither did what I expected, in opposite directions.",
"",
"Blocking is not a smooth dial at all. The estimates bounce around and travel",
"well outside the two endpoints they connect. For Freedom House they run from",
sprintf("%.3f to %.3f, while the endpoints themselves are only %.3f and %.3f.",
        min(filter(fam, panel == "Freedom House", family == "Blocking the years")$income),
        max(filter(fam, panel == "Freedom House", family == "Blocking the years")$income),
        filter(fam, panel == "Freedom House", family == "Blocking the years",
               setting == 1)$income,
        filter(fam, panel == "Freedom House", family == "Blocking the years",
               setting == max(setting))$income),
"Even the instrument count refuses to fall as the blocks get wider. Block sizes",
"of four and five both give three blocks, but the wider blocks reach back to",
"years that have deeper lags on offer, so the count goes up rather than down.",
"Block size is therefore not a measure of how much pooling is going on, which is",
"why the line looks like noise.",
"",
"Fading the older lags does almost nothing. Across the whole range of factors",
sprintf("the estimate moves by %.3f for Freedom House and %.3f for Polity.",
        diff(range(filter(fam, panel == "Freedom House",
                          family == "Fading out older lags")$income)),
        diff(range(filter(fam, panel == "Polity",
                          family == "Fading out older lags")$income))),
"That is worth knowing rather than disappointing: it says the period-pooled",
"answer is not an artifact of counting every lag equally, since counting the",
"older ones for a tenth as much barely moves it.",
"",
"One note on reading the table against the figure. The table's second column",
"uses every available lag, which is past the right-hand edge of the figure; the",
"figure stops at a window of 8. The instrument counts in the table count only",
"the lagged levels, not the year dummies that also sit in the instrument set.")
writeLines(md, file.path(PATH_DOCS, "aggregation.md"))

cat("Aggregation comparison written.\n")
print(filter(tab, is.infinite(lag_max)) |>
        transmute(panel, design = substr(design, 1, 9), scheme,
                  income = round(income, 3), se = round(income_se, 3), n_inst),
      n = 40)
