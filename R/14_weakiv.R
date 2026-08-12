source(here::here("R", "00_setup.R"))

DEP <- DEP_WEAKIV

# columns
rows <- list(); curves <- list(); checks <- list()
for (sp in TSLS_SPECS) {
  dat <- spec_data(sp)
  k <- length(sp$inst)
  m <- fit_iv(dat, DEP, endog = "Linc", inst = sp$inst, exog = sp$exog, country_fe = TRUE)
  fs <- fit_first_stage(dat, endog = "Linc", inst = c(sp$inst, sp$exog), country_fe = TRUE)
  G <- mod_nc(m)
  est <- ce(m, "Linc")
  tcrit <- qt(1 - CI_LEVEL / 2, G - 1)

  want <- published(table == sp$tab, column == sp$col, row == "inc", type == "coef")$value
  stopifnot(length(want) == 1L, abs(round(est["est"], 3) - want) < 1e-9)

  fF <- wald_stat(fs, sp$inst, G)$F
  fEff <- effective_f(dat, "Linc", sp$inst, sp$exog, fs)
  if (k == 1L) stopifnot(abs(fEff - fF) < 1e-8)

  ar <- ar_pcurve(dat, DEP, "Linc", sp$inst, sp$exog, AR_GRID)
  clr <- clr_pcurve(dat, DEP, "Linc", sp$inst, sp$exog, CLR_GRID)
  clrr <- clr_robust_pcurve(dat, DEP, "Linc", sp$inst, sp$exog, CLR_GRID)

  if (k == 1L) {
    at_hat <- which.min(abs(AR_GRID - est["est"]))
    stopifnot(ar$p[at_hat] > 0.99,
              max(abs(clrr$p_at(CLR_GRID) - pchisq(ar$stat_at(CLR_GRID), 1,
                                                   lower.tail = FALSE))) < 1e-10)
  }

  ai <- set_info(AR_GRID, ar$p, ar$p_at)
  ci <- set_info(CLR_GRID, clr$p, clr$p_at)
  ri <- set_info(CLR_GRID, clrr$p, clrr$p_at)
  ar_bounded <- fF > qf(1 - CI_LEVEL, k, G - 1)
  stopifnot(identical(ai$grid_bounded, ar_bounded), ai$pieces == 1L, ri$pieces == 1L)

  lab <- spec_label(sp$tab, sp$col)
  rows[[length(rows) + 1]] <- tibble(
    spec = lab, table = sp$tab, column = sp$col,
    instruments = paste(sp$inst, collapse = " + "), k = k,
    countries = G, obs = mod_nobs(m),
    first_stage_F = fF, effective_F = fEff,
    stock_yogo = sy10(1L, k), mop_cv = MOP_CV_10,
    tsls = unname(est["est"]), tsls_se = unname(est["se"]),
    tsls_lo = unname(est["est"] - tcrit * est["se"]),
    tsls_hi = unname(est["est"] + tcrit * est["se"]),
    ar_set = ai$text, ar_lo = ai$lo, ar_hi = ai$hi, ar_bounded = ar_bounded,
    ar_p0 = ar$p_at(0),
    clrr_set = ri$text, clrr_lo = ri$lo, clrr_hi = ri$hi, clrr_p0 = clrr$p_at(0),
    clr_set = ci$text, clr_lo = ci$lo, clr_hi = ci$hi, clr_p0 = clr$p_at(0))
  curves[[length(curves) + 1]] <- bind_rows(
    tibble(spec = lab, beta = AR_GRID, p = ar$p, test = unname(TEST_LABEL["ar"])),
    tibble(spec = lab, beta = CLR_GRID, p = clrr$p, test = unname(TEST_LABEL["clr_robust"])),
    tibble(spec = lab, beta = CLR_GRID, p = clr$p, test = unname(TEST_LABEL["clr_hom"])))
}

# collect
tab <- bind_rows(rows)
curve <- bind_rows(curves) |> mutate(spec = factor(spec, levels = tab$spec))
write_csv(tab, file.path(PATH_OUTPUT, "weakiv.csv"))

# checks
CHECK_GRID <- c(-0.5, 0, 0.25)
for (sp in TSLS_SPECS) {
  dat <- spec_data(sp)
  a <- ar_pcurve(dat, DEP, "Linc", sp$inst, sp$exog, CHECK_GRID)
  px <- a$Zt %*% solve(crossprod(a$Zt), crossprod(a$Zt, a$xt))
  b_hand <- as.numeric(crossprod(px, a$yt) / crossprod(px, a$xt))
  b_fit <- unname(ce(fit_iv(dat, DEP, endog = "Linc", inst = sp$inst,
                            exog = sp$exog, country_fe = TRUE), "Linc")["est"])
  slow <- vapply(CHECK_GRID, function(b0) {
    d2 <- dat; d2$e0 <- d2[[DEP]] - b0 * d2$Linc
    mm <- fit_ols(d2, "e0", c(sp$inst, sp$exog), country_fe = TRUE)
    wald_stat(mm, sp$inst, mod_nc(mm))$p
  }, numeric(1))
  stopifnot(abs(b_hand - b_fit) < 1e-9, max(abs(slow - a$p)) < 1e-8)
}
cat("Partialling reproduces the 2SLS estimates, and the fast Anderson-Rubin curve\n",
    "matches the clustered regression it stands in for.\n", sep = "")

# anderson-hsiao
ah <- list()
for (ms in MEASURES) {
  s <- filter(prep_dynamic(FILE_P5, ms$dep, ms$inc), sample == 1)
  dat <- complete_on(s, c("y", AH_SPEC$endog, AH_SPEC$inst, "code"))
  tb <- if (ms$dep == MEASURES[[1]]$dep) AH_SPEC$tab[1] else AH_SPEC$tab[2]
  m <- fit_iv(dat, "y", endog = AH_SPEC$endog, inst = AH_SPEC$inst, country_fe = FALSE)
  want <- published(table == tb, column == AH_SPEC$col, row == "inc", type == "coef")$value
  stopifnot(abs(round(unname(ce(m, "dLinc")["est"]), 3) - want) < 1e-9)
  G <- mod_nc(m)
  per <- vapply(AH_SPEC$endog, function(v)
    wald_stat(fit_ols(dat, v, AH_SPEC$inst, country_fe = FALSE), AH_SPEC$inst, G)$F,
    numeric(1))
  ah[[length(ah) + 1]] <- tibble(
    spec = spec_label(tb, AH_SPEC$col), panel = ms$label, countries = G, obs = mod_nobs(m),
    f_dem = per[[1]], f_inc = per[[2]],
    cragg_donald = cragg_donald(dat, AH_SPEC$endog, AH_SPEC$inst, character(),
                                fe = FE_YEAR),
    stock_yogo = sy10(length(AH_SPEC$endog), length(AH_SPEC$inst)),
    est = unname(ce(m, "dLinc")["est"]), se = unname(ce(m, "dLinc")["se"]))
}
ahtab <- bind_rows(ah)
write_csv(ahtab, file.path(PATH_OUTPUT, "weakiv_anderson_hsiao.csv"))

# table
txt <- c(sprintf("%-16s %-14s %3s %9s %9s %-20s %-22s %-22s",
                 "Spec", "Instrument(s)", "k", "1st-st F", "eff. F",
                 "2SLS (SE)", "Conventional 95% CI", "Anderson-Rubin 95%"))
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  txt <- c(txt, sprintf("%-16s %-14s %3d %9.2f %9.2f %-20s %-22s %-22s",
    r$spec, r$instruments, r$k, r$first_stage_F, r$effective_F,
    sprintf("%s (%s)", num(r$tsls), num(r$tsls_se)),
    sprintf("[%.2f, %.2f]", r$tsls_lo, r$tsls_hi), r$ar_set))
}
txt <- c(txt, "",
         sprintf("Thresholds: Stock-Yogo 10%% maximal size = %s (by instrument count);",
                 commas(vapply(sort(unique(tab$k)), function(x) sy10(1L, x), numeric(1)), 2)),
         sprintf("Montiel Olea-Pflueger, worst case for 10%% Nagar bias = %.2f.", MOP_CV_10),
         "", "95% confidence sets:",
         sprintf("%-16s %-27s %-22s %-22s", "Spec", TEST_LABEL[["ar"]],
                 "CLR, cluster-robust", "CLR, homoskedastic"))
for (i in seq_len(nrow(tab)))
  txt <- c(txt, sprintf("%-16s %-27s %-22s %-22s", tab$spec[i], tab$ar_set[i],
                        tab$clrr_set[i], tab$clr_set[i]))
txt <- c(txt, "", "Anderson-Hsiao columns (two endogenous regressors, no CLR counterpart):",
         sprintf("%-16s %-16s %11s %11s %13s %9s", "Spec", "Panel", "F: democracy",
                 "F: income", "Cragg-Donald", "SY 10%"))
for (i in seq_len(nrow(ahtab))) {
  r <- ahtab[i, ]
  txt <- c(txt, sprintf("%-16s %-16s %11.2f %11.2f %13.2f %9.2f",
                        r$spec, r$panel, r$f_dem, r$f_inc, r$cragg_donald, r$stock_yogo))
}
writeLines(txt, file.path(PATH_OUTPUT, "weakiv.txt"))

# figure
fig <- ggplot(curve, aes(beta, p, color = test)) +
  geom_hline(yintercept = CI_LEVEL, color = "grey40", linetype = "dashed") +
  geom_line(linewidth = 0.5) +
  geom_vline(data = tab, aes(xintercept = tsls), color = "grey30",
             linetype = "dotted", linewidth = 0.4) +
  facet_wrap(~spec, ncol = 5) +
  coord_cartesian(xlim = c(-1, 0.6), ylim = c(0, 1)) +
  labs(x = "Value of the income coefficient being tested", y = "p-value",
       color = NULL,
       title = "Weak-instrument-robust p-values for the effect of income",
       subtitle = paste0("Dashed line is ", CI_LEVEL * 100,
                         "%; anything above it stays in the confidence set. ",
                         "Dotted line is the 2SLS estimate.")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top")
ggsave(file.path(PATH_OUTPUT, "weakiv.png"), fig, width = 11, height = 5.5, dpi = 150)

# summary
n_sy <- sum(tab$first_stage_F < tab$stock_yogo)
n_mop <- sum(tab$effective_F < tab$mop_cv)
n_zero <- sum(tab$ar_p0 >= CI_LEVEL)
n_single <- sum(tab$k == 1L)
n_multi <- nrow(tab) - n_single
widen <- tab |>
  mutate(conv = tsls_hi - tsls_lo, ar = ar_hi - ar_lo, ratio = ar / conv) |>
  arrange(desc(ratio))
stopifnot(all(widen$ratio >= 1), all(tab$ar_lo <= tab$tsls), all(tab$ar_hi >= tab$tsls))
ord <- arrange(widen, effective_F)
inv <- NULL
for (i in seq_len(nrow(ord) - 1L)) for (j in (i + 1L):nrow(ord))
  if (is.null(inv) && ord$ratio[i] < ord$ratio[j]) inv <- c(i, j)
strongest <- tab[which.max(tab$effective_F), ]
weakest <- tab[which.min(tab$effective_F), ]
multi <- filter(tab, k > 1L)
cal_gap <- max(abs(filter(tab, k == 1L)$ar_lo - filter(tab, k == 1L)$clrr_lo),
               abs(filter(tab, k == 1L)$ar_hi - filter(tab, k == 1L)$clrr_hi))
robust_wider <- with(multi, (clrr_hi - clrr_lo) / (clr_hi - clr_lo))
ah_cols <- paste(sprintf("%s (%s)", ahtab$spec, ahtab$panel), collapse = " and ")

# report
write_doc("weak-instruments.md",
"# Are the instruments strong enough?",
sprintf(paste(
  "Professor Torgovitsky asked for weak instrument diagnostics, then for",
  "Anderson-Rubin confidence sets and the Moreira conditional likelihood ratio",
  "approach. The natural place for all three is the two-stage least squares columns",
  "of Tables 5 and 6. Those columns instrument one variable, income, with a small",
  "number of instruments. The instrument is the savings rate in Table 5 and",
  "trade-weighted world income in Table 6. That gives %s columns."), spell(nrow(tab))),
paste(
  "The numbers are in output/weakiv.csv and output/weakiv.txt. The file",
  "output/weakiv.png draws the p-value curves. The tests come from Anderson and Rubin",
  "(1949), Moreira (2003) and Kleibergen (2005). The first-stage cutoffs come from",
  "Stock and Yogo (2005) and Montiel Olea and Pflueger (2013). Full details are in",
  "docs/references.md."),
paste(
  "The section \"What this does not reach\" covers the paper's other",
  "instrumental-variables columns. They are the Anderson-Hsiao and Arellano-Bond",
  "columns, which instrument two variables at once. The confidence sets here have no",
  "direct counterpart there. I report first-stage diagnostics for those columns."),

"## Why this matters",
paste(
  "A weak instrument barely moves the variable that it stands in for. With a weak",
  "instrument, the usual standard errors and confidence intervals are not",
  "trustworthy. They rest on the assumption that the instrument has real pull. The",
  "tests below stay honest whether the instrument is strong or weak. A comparison",
  "against the ordinary interval shows how much the ordinary interval relied on that",
  "assumption."),

"## How strong is the first stage",
sprintf(paste(
  "The first-stage F statistics run from %.1f to %.1f. I clustered them by country,",
  "the same way the paper clusters everything else. With one endogenous regressor this",
  "is the cluster-robust first-stage Wald F that standard software reports. Every",
  "column clears the old rule of thumb that an F of more than 10 is good enough.",
  "Against the stricter Stock and Yogo cutoff, %s of the %s columns fall short."),
  min(tab$first_stage_F), max(tab$first_stage_F), spell(n_sy), spell(nrow(tab))),
paste(
  "That last comparison is not strictly legitimate. Last round I said so and left it",
  "there. Stock and Yogo derive their cutoffs for well behaved, unclustered errors. I",
  "cluster these F statistics by country. Under clustering, the properly justified",
  "statistic is the effective F of Montiel Olea and Pflueger (2013). This round I",
  "implemented it, and it is in the table beside the ordinary one."),
sprintf(paste(
  "I have two observations about what the effective F changes. First, for the %s",
  "columns with a",
  "single instrument, the effective F is algebraically identical to the cluster-robust",
  "first-stage F. The script checks this identity, and does not assume it. So nothing",
  "moves there."), spell(n_single)),
sprintf(paste(
  "Second, for the %s columns with two instruments the effective F moves a lot, and in",
  "opposite directions. %s. The effective F divides by a trace that accounts for how",
  "unevenly the instruments contribute. The ordinary F does not."),
  spell(n_multi),
  paste(sprintf("%s goes from %.1f to %.1f", multi$spec, multi$first_stage_F,
                multi$effective_F), collapse = ". ")),
sprintf(paste(
  "The two cutoffs are not the same kind of object, so I report them separately rather",
  "than in one column. Stock and Yogo bound the size distortion of the conventional",
  "test. Montiel Olea and Pflueger bound the bias of the estimator relative to a worst",
  "case. Their conservative critical value of %.2f holds whatever the covariance",
  "structure turns out to be. Against that value, %s of the %s columns fall short."),
  MOP_CV_10, spell(n_mop), spell(nrow(tab))),
paste(
  "Either way the reading is the same. These instruments are not badly weak, and",
  "several are not comfortably strong. The two confidence sets below handle exactly",
  "that situation."),

"## Anderson-Rubin",
paste(
  "The Anderson-Rubin test asks a simple question. The test takes a candidate value",
  "for the effect of income. It subtracts what that value implies, then looks for",
  "anything left over that the instruments still explain. If something is left, the",
  "test rejects that candidate value. Every value that survives goes into the",
  "confidence set. That set stays valid no matter how weak the instruments are."),
paste(
  "I clustered the test by country to match the paper. That makes it exactly the",
  "cluster-robust Wald test on the auxiliary regression. The script checks it against",
  "an actual fit of that regression."),
sprintf(paste(
  "%s %s sets are bounded, and %s contain zero. The width of the grid I searched",
  "does not decide whether a set is bounded. With this form of the statistic, the",
  "Anderson-Rubin statistic tends to the first-stage F as the candidate value runs off",
  "to infinity. So the set is bounded exactly when that F exceeds the critical value.",
  "The script computes boundedness both ways. If the two answers disagree, the script",
  "stops."),
  ifelse(all(tab$ar_bounded), "All", paste(sum(tab$ar_bounded), "of the")),
  spell(nrow(tab)),
  ifelse(n_zero == nrow(tab), "all of them", paste(spell(n_zero), "of them"))),
sprintf(paste(
  "The first stage is strongest in %s, at an effective F of %.1f. There the set is %s,",
  "which is nearly the ordinary interval of [%.2f, %.2f]. The first stage is weakest",
  "in %s, at an effective F of %.1f. There the set is %s, against an ordinary interval",
  "of [%.2f, %.2f]. The Anderson-Rubin set is at least as wide as the ordinary interval",
  "in every column. So the ordinary intervals are optimistic throughout, not only at",
  "the extremes."),
  strongest$spec, strongest$effective_F, strongest$ar_set,
  strongest$tsls_lo, strongest$tsls_hi,
  weakest$spec, weakest$effective_F, weakest$ar_set,
  weakest$tsls_lo, weakest$tsls_hi),
sprintf(paste(
  "I must not push that comparison further than it goes. The columns do not form a",
  "tidy order by first-stage strength. The widening ratio runs from %.2f to %.2f. It",
  "does not decrease neatly as the first stage gets stronger."),
  min(widen$ratio), max(widen$ratio)),
sprintf(paste(
  "%s has the weaker first stage of the two, an effective F of %.1f against %.1f. It",
  "still widens less than %s does. The two widening ratios are %.2f against %.2f. These",
  "specifications differ in their controls and samples as well as in instrument",
  "strength. So the two ends of the range are suggestive rather than a clean",
  "relationship."),
  ord$spec[inv[1]], ord$effective_F[inv[1]], ord$effective_F[inv[2]],
  ord$spec[inv[2]], ord$ratio[inv[1]], ord$ratio[inv[2]]),

"## Moreira's conditional likelihood ratio",
paste(
  "The conditional likelihood ratio test can be sharper than Anderson-Rubin. It does",
  "not compare the statistic against a fixed cutoff. Instead it adjusts the cutoff with",
  "a quantity that carries the information about instrument strength. That adjustment",
  "buys back power that Anderson-Rubin gives away. I coded it directly, so one grid",
  "inversion serves all three tests here. R/15_ivcrosscheck.R checks it against the",
  "ivmodel package, which implements the same test."),
sprintf(paste(
  "The catch is that %s of the %s columns here have exactly one instrument. With one",
  "instrument there is nothing to buy back, and the two tests coincide. So the extra",
  "power only matters for the %s columns with a pair of instruments. For the",
  "single-instrument columns the cutoff is a chi-squared quantile with an exact form.",
  "The script uses that form rather than a simulation of a distribution it already",
  "knows."),
  spell(n_single), spell(nrow(tab)), spell(n_multi)),
paste(
  "Moreira's test assumes that the errors are well behaved and independent. The",
  "packages implement it that way. The paper clusters, and the Anderson-Rubin sets",
  "above cluster. So last round I reported the two side by side, with a caveat. My",
  "caveat was that the gap between them was mostly the clustering rather than the extra",
  "power. That was a conjecture."),
paste(
  "This round I also implemented the cluster-robust version, which is Kleibergen's",
  "(2005) generalization. It is the same conditional likelihood ratio statistic, built",
  "from moment conditions. I estimate the variance of those moment conditions by",
  "country cluster. All three sets are in output/weakiv.txt."),
sprintf(paste(
  "The cluster-robust version settles the question. For the %s single-instrument",
  "columns the cluster-robust statistic equals the clustered Anderson-Rubin statistic.",
  "It must do so, and the script checks that equality to ten decimal places at every",
  "candidate value on the grid. The two printed sets still differ by at most %.3f at an",
  "endpoint, because the two tests read that one statistic off different reference",
  "distributions. The likelihood ratio test uses chi-squared and the Anderson-Rubin test",
  "uses F. So for those columns the entire gap between the homoskedastic set and the",
  "Anderson-Rubin set was the clustering, as I guessed."),
  spell(n_single), cal_gap),
sprintf(paste(
  "For the %s columns with two instruments the answer is more interesting, and it is",
  "not what I guessed. The cluster-robust version widens the set by a factor of %.2f",
  "and %.2f. The result is still inside the Anderson-Rubin set rather than equal to it:",
  "%s against %s, and %s against %s. That remaining difference is the extra power the",
  "conditional likelihood ratio test buys, and it is real."),
  spell(n_multi), robust_wider[1], robust_wider[2],
  multi$clrr_set[1], multi$ar_set[1], multi$clrr_set[2], multi$ar_set[2]),
paste(
  "So the honest split is this. For eight columns the gap was all clustering. For the",
  "two overidentified ones it is part clustering and part genuine power."),

"## What it adds up to",
sprintf(paste(
  "The inference here does not lean on strong instruments. Under that inference, %s of",
  "the %s two-stage least squares columns leaves zero inside the confidence set. So",
  "none of them establishes an effect of income on democracy in either direction."),
  ifelse(n_zero == nrow(tab), "every one", spell(n_zero)), spell(nrow(tab))),
paste(
  "This does not overturn anything in the paper. The paper argues that the effect is",
  "not there once the fixed differences between countries come out. The result does",
  "mean that the instrumental-variables columns are uninformative about the sign of the",
  "effect. They are not evidence for a negative effect."),

"## What this does not reach",
sprintf(paste(
  "The paper's other instrumental-variables columns instrument two variables at once,",
  "lagged democracy and lagged income. The machinery above does not transfer here, for",
  "two reasons. A confidence set for the income coefficient alone needs a projection of",
  "a two-dimensional region rather than a grid inversion. Moreira also derives his test",
  "for a single endogenous regressor, so there is no conditional likelihood ratio",
  "counterpart at all.")),
paste(
  "The Arellano-Bond columns add a third reason. Their instrument counts run into the",
  "dozens against about a hundred countries. This sample cannot support a",
  "cluster-robust weight matrix of that size. That last point is not a side issue. It",
  "is the subject of docs/instruments.md, where the same instrument counts inflate the",
  "apparent precision of the estimates."),
sprintf(paste(
  "The first-stage diagnostic does transfer. So I report it for the Anderson-Hsiao",
  "columns, %s. With two endogenous regressors the right statistic is Cragg and",
  "Donald's minimum eigenvalue rather than a per-equation F. The Stock-Yogo cutoff",
  "changes with it, to %.2f. The Cragg-Donald statistics are %s and %s, so both clear",
  "it. The per-equation F statistics, clustered by country, are %.1f and %.1f for",
  "democracy and %.1f and %.1f for income."),
  ah_cols, ahtab$stock_yogo[1], num(ahtab$cragg_donald[1], 2),
  num(ahtab$cragg_donald[2], 2),
  ahtab$f_dem[1], ahtab$f_dem[2], ahtab$f_inc[1], ahtab$f_inc[2]),
paste(
  "The Cragg-Donald statistic is the homoskedastic one, so it carries the same caveat",
  "as the Stock-Yogo comparison above. There is no clustered version of the cutoff to",
  "hold it against. Taken together these are not the weakest instruments in the paper.",
  "The reason those columns still produce very wide intervals is the small number of",
  "moment conditions, not a first stage that fails."),

"## Checks",
paste(
  "This is fiddly enough that I did not want to rely on code I wrote myself. Several",
  "checks run every time the script does. The partialling step must reproduce the",
  "published two-stage least squares estimates. The fast Anderson-Rubin curve must",
  "match the slower clustered regression that it stands in for. Boundedness from the",
  "grid must agree with boundedness from the first-stage F."),
paste(
  "When there is one instrument, the cluster-robust conditional likelihood ratio test",
  "must equal the clustered Anderson-Rubin test. The effective F must equal the",
  "ordinary F in the same case. Every Anderson-Rubin set must contain its own point",
  "estimate. It must also be at least as wide as the conventional interval, which is",
  "the claim the section above leans on. If a check fails, the script stops."),
paste(
  "The script does not read confidence set endpoints from the grid. The grid locates",
  "the crossing, and then the script solves for the endpoint. So the printed set is the",
  "true set rather than the nearest grid point inside it."),
sprintf(paste(
  "R/15_ivcrosscheck.R also sends every column through the ivmodel package. If the",
  "ivmodel conditional likelihood ratio sets disagree with mine, the script stops. If",
  "the ivmodel two-stage least squares estimates do not match mine exactly, the script",
  "also stops. At present both comparisons pass for all %s columns. That script needs",
  "the ivmodel package, so it sits outside the main run. The command for it is",
  "`Rscript R/15_ivcrosscheck.R`."), spell(nrow(tab))))

cat("Weak-instrument diagnostics written.\n")
print(transmute(tab, spec, k, F = round(first_stage_F, 1), effF = round(effective_F, 1),
                tsls = round(tsls, 3), ar_set, clrr_set, clr_set), n = 20)
print(ahtab, n = 10)
