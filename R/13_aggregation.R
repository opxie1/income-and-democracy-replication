source(here::here("R", "00_setup.R"))
source(here::here("R", "_aggregation.R"))

# fits
res <- lapply(MEASURES, function(ms) agg_sweep(ms$dep, ms$inc, ms$label))
tab <- bind_rows(lapply(res, function(r) r$tab))
bench <- bind_rows(lapply(res, function(r) r$bench))
fam <- bind_rows(lapply(MEASURES, function(ms) family_sweep(ms$dep, ms$inc, ms$label)))
cat("Both families reproduce the fixed rules at their endpoints.\n")
write_csv(fam, file.path(PATH_OUTPUT, "aggregation_families.csv"))

# checks
stopifnot(all(filter(tab, scheme == "full")$n_inst == filter(tab, scheme == "full")$n_par),
          all((tab$n_inst - tab$n_par == 0L) == (tab$hansen_df == 0L)),
          all(is.na(tab$hansen_p) == (tab$hansen_df == 0L)))

# identities
chk_fr <- est_frame(MEASURES[[1]]$dep, MEASURES[[1]]$inc)
chk_fit <- function(model, scheme) fit_abgmm(chk_fr$full, chk_fr$est,
  dep_level = MEASURES[[1]]$dep, endog = c("dLdep", "dLinc"), inst_extra = "L2inc",
  scheme = scheme, model = model)
just1 <- chk_fit("onestep", "full"); just2 <- chk_fit("twostep", "full")
over2 <- chk_fit("twostep", "lag")
stopifnot(just1$n_inst == just1$n_par, just2$hansen_df == 0L,
          max(abs(just1$coef - just2$coef)) < 1e-8, just2$hansen < 1e-8,
          over2$hansen_df > 0L, over2$hansen > 0)
cat("Two-step reproduces one-step when the model is exactly identified,",
    "and its Hansen statistic is then zero.\n")

# published
paper <- filter(tab, scheme == "none", is.infinite(lag_max), design_short == "paper")
pub <- published(table %in% c("2", "3"), column == 4, row == "inc", type == "coef") |>
  mutate(panel = ifelse(table == "2", MEASURES[[1]]$label, MEASURES[[2]]$label))
chk <- left_join(paper, select(pub, panel, published = value), by = "panel")
stopifnot(nrow(chk) == length(MEASURES), !anyNA(chk$published),
          all(abs(round(chk$income, 3) - chk$published) < 1e-9))
cat("Uncollapsed, all lags reproduces the published GMM column for both measures.\n")

write_csv(tab, file.path(PATH_OUTPUT, "aggregation.csv"))

# subsets
sub <- bind_rows(lapply(MEASURES, function(ms) subset_sweep(ms$dep, ms$inc, ms$label)))
write_csv(sub, file.path(PATH_OUTPUT, "aggregation_subsets.csv"))


# outputs
out <- write_aggregation_outputs(tab, bench, fam, sub)
fixed_ref <- out$fixed_ref
subsum <- out$subsum

# summary
rng <- filter(tab, is.infinite(lag_max)) |> group_by(panel, design_short) |>
  summarise(lo = min(income), hi = max(income), .groups = "drop")
mv <- filter(tab, is.finite(lag_max)) |> group_by(panel, design, design_short, scheme) |>
  summarise(rng = max(income) - min(income), .groups = "drop")
flat <- mv |> group_by(panel, design_short) |> mutate(place = rank(rng)) |> ungroup()
n_panel <- n_distinct(paste(flat$panel, flat$design_short))
lag_first <- filter(flat, scheme == "lag", place == 1)
lag_lost <- filter(flat, scheme == "lag", place > 1) |>
  left_join(select(filter(flat, place == 1), panel, design_short,
                   win = scheme, win_rng = rng), by = c("panel", "design_short")) |>
  mutate(win_label = unname(SCHEME_LABEL[win]))
shape <- filter(tab, is.finite(lag_max)) |>
  left_join(bench, by = "panel") |> arrange(lag_max) |>
  group_by(panel, design_short, scheme) |>
  summarise(net = abs(first(income) - first(fe)) - abs(last(income) - last(fe)),
            turn = lag_max[which.max(abs(income - fe))],
            mono = all(diff(income) > 0), .groups = "drop")
un <- function(pl, col) filter(shape, design_short == "paper", panel == pl,
                               scheme == "none")[[col]]
worst_se <- tab[which.max(tab$income_se), ]

# overid
POOLED <- c("lag", "period", "period_mean")
pooled_p <- filter(tab, scheme %in% POOLED, hansen_df > 0L)
best_p <- pooled_p |> group_by(panel, design_short, lag_max) |>
  summarise(has_lag = any(scheme == "lag"),
            lag_p = if (any(scheme == "lag")) max(hansen_p[scheme == "lag"]) else NA_real_,
            rival = max(hansen_p[scheme != "lag"]), .groups = "drop") |>
  filter(has_lag) |> mutate(lag_wins = lag_p > rival)
stopifnot(nrow(best_p) > 0, mean(best_p$lag_wins) > 0.9)
lag_loss <- filter(best_p, !lag_wins)
rej <- pooled_p |> group_by(scheme) |>
  summarise(n_rej = sum(hansen_p < CI_LEVEL), n = n(), lo = min(hansen_p),
            hi = max(hansen_p), .groups = "drop")
rj <- function(sc, col) filter(rej, scheme == sc)[[col]]
unc_p <- filter(tab, scheme == "none", hansen_df > 0L)
fh <- MEASURES[[1]]$label
po <- MEASURES[[2]]$label
pr <- function(pl, sc, col = "rng")
  filter(mv, design_short == "paper", panel == pl, scheme == sc)[[col]]
sy <- function(pl, sc, col = "rng")
  filter(mv, design_short == "symmetric", panel == pl, scheme == sc)[[col]]
sub_wide <- subsum |> group_by(panel) |>
  summarise(small = sd[which.min(size)], big = sd[which.max(size)],
            med_small = med[which.min(size)], med_big = med[which.max(size)],
            lo_small = lo[which.min(size)], hi_small = hi[which.min(size)],
            se_small = med_se[which.min(size)], se_big = med_se[which.max(size)],
            s_small = min(size), s_big = max(size),
            rej_lo = min(rej), rej_hi = max(rej), .groups = "drop")
sw <- function(pl, col) filter(sub_wide, panel == pl)[[col]]
sub_total <- filter(tab, scheme == "none", is.infinite(lag_max),
                    design_short == "paper") |> select(panel, income, n_gmm)
st <- function(pl, col) filter(sub_total, panel == pl)[[col]]
stopifnot(all(abs(sub_wide$med_big - st(sub_wide$panel, "income")) < 0.02),
          all(sub_wide$big < sub_wide$small))
rej_series <- function(pl) commas(filter(subsum, panel == pl) |> arrange(size) |>
                                    pull(rej), 2)
near <- function(pl) {
  ref <- filter(fixed_ref, panel == pl, scheme == "lag")
  s <- filter(subsum, panel == pl)
  row <- s[which.min(abs(s$n_inst - ref$n_inst)), ]
  list(est = ref$income, n_ref = ref$n_inst, n_draw = row$n_inst,
       lo = row$lo, hi = row$hi, inside = ref$income >= row$lo && ref$income <= row$hi)
}

# report
write_doc("aggregation.md",
"# Other ways to collapse the instruments",
paste(
  "Professor Torgovitsky asked whether Roodman describes other ways to collapse the",
  "instruments. He asked because a collapse rule is a choice, and other choices exist.",
  "This file reports what the paper says. It also reports what the alternatives give on",
  "this data. The numbers are in output/aggregation.txt and output/aggregation.csv. The",
  "figure is output/aggregation.png."),

"## What Roodman describes",
paste(
  "The paper is Roodman (2009), \"A Note on the Theme of Too Many Instruments\",",
  "Oxford Bulletin of Economics and Statistics 71(1), 135-158. Full details for",
  "everything cited here are in docs/references.md."),
paste(
  "Section V is titled \"Techniques for reducing the instrument count\", and it opens",
  "with two techniques rather than one. The first technique uses only certain lags",
  "instead of all available ones. It still makes separate instruments for each period,",
  "but it caps the number per period. The count then grows in proportion to the length",
  "of the panel, and not with its square."),
paste(
  "Roodman describes this technique as a projection of the regressors onto the full",
  "instrument set, with the coefficients on certain lags held at zero. He attributes",
  sprintf("that reading to %s, his working paper on optimal instrumental variables.",
          CITE_ARELLANO_OPT)),
paste(
  "The second technique, which he calls the less common one, is \"to combine",
  "instruments through addition into smaller sets\". Its advantage is \"retaining more",
  "information, since no lags are actually dropped\". One sentence matters most here,",
  "because it says what collapsing does. It \"is equivalent to imposing the constraint",
  "in projecting regressors onto HENR instruments that certain subsets have the same",
  "coefficient\". The picture he gives is of \"squeezing the matrix ... horizontally and",
  "adding together formerly distinct columns\". He gives one rule for the choice of",
  "columns: \"collapsed instruments are straightforward conceptually: one is made for",
  "each lag distance, with 0 substituted for any missing values\"."),
paste(
  "Those two sentences license most of what follows. Where they do not, I say so. A cap",
  "on the lags sets some of the projection coefficients to zero. Collapsing sets some of",
  "them equal to each other."),
paste(
  "The choice of which coefficients to set equal is free, and the paper makes only one",
  "choice available. Any other partition of the instrument columns into groups is",
  "another collapsing rule in the same sense. So it is fair to ask how much the answer",
  "depends on his particular rule."),
paste(
  "The section closes with a statement of what these techniques are for. The statement",
  "deserves a quotation, because it sets the standard for this file. They \"provide the",
  "basis for some minimally arbitrary robustness and specification tests for Difference",
  "and System GMM: cut the instrument count in one of these ways and examine the",
  "behavior of the coefficient estimates and Hansen and Difference-in-Hansen tests\". So",
  "the coefficient on its own is only half of what he asks for. The table in",
  "output/aggregation.txt now carries an overidentification p-value in every cell,",
  "beside the estimate."),
paste(
  "I do not report the Difference-in-Hansen half here. That test compares one nested",
  "subset of instruments against the rest, and the rules compared here are not nested",
  "inside one another. It belongs in the difference-versus-system comparison, where the",
  "extra level conditions form a subset that I can add and remove. Both",
  "docs/alternatives.md and docs/instruments.md make that comparison."),
paste(
  "He also notes that the two techniques work together. The pair leaves a count that",
  "does not grow with the length of the panel at all. His Table 1 crosses them. It shows",
  "four variants of system GMM: the full instrument set, one-period lags only, the",
  "collapsed set, and both restrictions at once."),
paste(
  "Two further ideas appear in footnotes rather than the main text.",
  sprintf("Footnote 6 describes an approach from %s.", CITE_ARELLANO_OPT),
  "That approach first models the instrumenting variables as a group, as functions of",
  "their collective lags, with a vector autoregression. It then turns the fitted",
  "coefficients into constraints on the projection. Roodman notes that it \"has yet to",
  "enter common practice\"."),
paste(
  "Footnote 7 suggests \"repeatedly selecting random subsets from the collection of",
  "potential instruments and investigating how key results such as coefficients of",
  "interest and the p value on the J statistic vary with the number of instruments\".",
  "That idea is cheap to run, so I ran it. It is the last section here."),

"## The alternatives I tried",
paste(
  "Roodman adds across time and keeps the lag distances apart. The obvious alternative",
  "reverses this. It adds across lag distances and keeps the periods apart, which gives",
  "one instrument per period. Addition and the average differ here, because the number",
  "of available lags grows over time. So the average across lag distances is a third",
  "option. A fourth option pushes the rule to the end, where every column collapses into",
  "one instrument."),
paste(
  "The average is also the one rule here that steps outside the sentence quoted above.",
  "Its divisor is the count of lags that the country has in that period. So the divisor",
  "varies from country to country inside a single column. The result is not the original",
  "columns with some coefficients tied together, but a new instrument. The new instrument",
  "is still legitimate, because the divisor depends only on which lags exist and not on",
  "the outcome. The equal-coefficient reading does not reach it, and the same caveat",
  "applies to the fading family further down."),
sprintf(paste(
  "I ran all %s collapse rules against the uncollapsed set. I then ran %s further",
  "families with a knob on them. The section after next covers those families."),
  spell(length(AGG_SCHEMES) - 1L), spell(length(AGG_FAMILIES))),
paste(
  "I also ran two instrument designs, because the design matters. The paper builds its",
  "lag blocks from democracy only and gives income a single lagged level. The symmetric",
  "alternative gives both variables the full block of lags. The plm package uses that",
  "design by default, and so does the sweep in docs/instruments.md. The first design",
  "reproduces the published GMM column exactly. If it ever fails, the script stops."),
paste(
  "One reading note comes before the numbers. The collapse of every column into one",
  "leaves the model exactly identified here. That row therefore has no overidentifying",
  "restrictions left, and output/aggregation.txt prints `exact id` instead of a p-value.",
  "The estimate is still real, and it is the least stable one on the figure. It is not",
  "comparable to the others on anything that counts restrictions."),

"## What came out",
sprintf(paste(
  "The collapse rule matters, and by more than I expected. With every lag under the",
  "paper's design, the estimates for %s range from %.3f to %.3f. The estimates for %s",
  "range from %.3f to %.3f. Only the collapse rule changes across those numbers. The",
  "data, the lags and the estimator stay the same."),
  fh, filter(rng, panel == fh, design_short == "paper")$lo,
  filter(rng, panel == fh, design_short == "paper")$hi,
  po, filter(rng, panel == po, design_short == "paper")$lo,
  filter(rng, panel == po, design_short == "paper")$hi),
sprintf(paste(
  "The top row of the figure shows the paper's own instrument design. There the",
  "uncollapsed set is the exception. It travels much further than any of the collapse",
  "rules. Its range is %.3f for %s and %.3f for %s, against at most %.3f for the %s",
  "collapse rules. It also finishes closer to the fixed-effects value than it started,",
  "by %.3f and %.3f. That is the overfitting story again."),
  pr(fh, "none"), fh, pr(po, "none"), po,
  max(filter(mv, design_short == "paper", scheme != "none")$rng),
  spell(length(AGG_SCHEMES) - 1L), un(fh, "net"), un(po, "net")),
sprintf(paste(
  "The path is not a straight climb, and I do not want to say that it is. Both",
  "uncollapsed lines first move further from the fixed-effects value, out to a window of",
  "%d for %s and %d for %s. They reverse direction only after that point. The",
  "uncollapsed line is monotone in %s of the %s panels on the figure. %sSo the drift is",
  "a net direction over the whole window, and not a steady march."),
  un(fh, "turn"), fh, un(po, "turn"), po,
  {
    n_mono <- nrow(filter(shape, scheme == "none", mono))
    if (n_mono == 0) "none" else spell(n_mono)
  },
  spell(n_panel),
  {
    m <- filter(shape, scheme == "none", mono)
    dg <- ifelse(m$design_short == "paper", "the paper's", "the symmetric")
    if (nrow(m) == 0) ""
    else if (nrow(m) == 1) sprintf("That panel is %s under %s design. ", m$panel[1], dg[1])
    else sprintf("Those panels are %s. ",
                 commas(sprintf("%s under %s design", m$panel, dg)))
  }),
sprintf(paste(
  "The bottom row gives income its own block of lags, and the result is messier. That",
  "mess deserves a plain statement. There the collapse rules are not uniformly",
  "steadier. Roodman's rule is the flattest line in %s of the %s panels, with a range of",
  "at most %.3f. %s"),
  spell(nrow(lag_first)), spell(n_panel), max(filter(mv, scheme == "lag")$rng),
  if (nrow(lag_lost) == 0) "" else sprintf(
    paste("It loses only in the %s panel under the %s, where the exactly identified",
          "fully collapsed column happens to sit still (%.3f against %.3f). That column",
          "is not much of a rival, because it has nothing left to test."),
    lag_lost$panel[1],
    ifelse(lag_lost$design_short[1] == "paper", "paper's design", "symmetric design"),
    lag_lost$win_rng[1], lag_lost$rng[1])),
sprintf(paste(
  "Collapsing by period swings more than the uncollapsed set does (%.3f against %.3f",
  "for %s). Collapsing everything into one column is the least stable line on the whole",
  "figure. So the useful statement is not that a collapse always steadies the estimate.",
  "Collapsing by lag distance, the rule Roodman proposes, is the one that stays steady",
  "under both designs."),
  sy(po, "period"), sy(po, "none"), po),
sprintf(paste(
  "Too much collapsing has its own failure mode. With the symmetric design, collapsing",
  "everything into one column produces a standard error of %.2f for %s. That number says",
  "the instrument no longer carries usable information. A middle ground exists. Too",
  "little collapsing makes the estimator overfit, and too much collapsing leaves too",
  "little information for identification. That second failure is the weak instrument",
  "problem, and docs/weak-instruments.md examines it."),
  worst_se$income_se, worst_se$panel),

"## What the overidentification test says",
paste(
  "Last round I looked only at stability. Stability on its own does not settle much. A",
  "rule can sit still because it is well behaved, or because it no longer listens to the",
  "data. Roodman's own prescription is to examine the overidentification test as the",
  "instrument count falls. The table now carries one test in every cell. The test is the",
  "sharper tool of the two, and it points the same way."),
sprintf(paste(
  "Among the three collapse rules that leave anything to test, Roodman's rule has the",
  "highest p-value in %d of the %d overidentified cells. The one exception is %s under",
  "the %s at the widest window. There the three rules cluster close to the threshold on",
  "both sides. The gap is not close anywhere else. His p-values run from %s to %s, and",
  "the test rejects his rule at the %d%% level in %d of %d cells. The test rejects the",
  "collapse by period in %d of %d cells, and the collapse by period after averaging in",
  "%d of %d."),
  sum(best_p$lag_wins), nrow(best_p),
  if (nrow(lag_loss)) lag_loss$panel[1] else "none",
  if (nrow(lag_loss) && lag_loss$design_short[1] == "paper") "paper's design"
    else "symmetric design",
  num(rj("lag", "lo"), 3), num(rj("lag", "hi"), 3), as.integer(CI_LEVEL * 100),
  rj("lag", "n_rej"), rj("lag", "n"), rj("period", "n_rej"), rj("period", "n"),
  rj("period_mean", "n_rej"), rj("period_mean", "n")),
paste(
  "The other two rules are my own, and they turn Roodman's construction on its side.",
  "Those two rules are not merely less steady than his. The data reject them, and the",
  "data do not reject his rule."),
sprintf(paste(
  "The uncollapsed set needs the most careful comparison. Its p-values run from %s to",
  "%s, and the test rejects it in %d of %d cells. For %s the p-value climbs steadily",
  "with the count, from %s at the narrowest window to %s at the widest under the",
  "paper's design. For %s under the same design it stays low throughout, and it never",
  "rises above %s. Where the test does look comfortable, that comfort is not proof of a",
  "valid instrument set. Those cells have up to %d instruments against about %d",
  "countries, and that is exactly where the test loses its power."),
  num(min(unc_p$hansen_p), 2), num(max(unc_p$hansen_p), 2),
  sum(unc_p$hansen_p < CI_LEVEL), nrow(unc_p), fh,
  num(min(filter(unc_p, panel == fh, design_short == "paper")$hansen_p), 2),
  num(max(filter(unc_p, panel == fh, design_short == "paper")$hansen_p), 2), po,
  num(max(filter(unc_p, panel == po, design_short == "paper")$hansen_p), 2),
  max(unc_p$n_inst), filter(tab, panel == fh)$countries[1]),
sprintf(paste(
  "Roodman's rule reaches %s on as few as %d instruments. Against that mark, a p-value",
  "of %s on a set of %d means much less."),
  num(rj("lag", "hi"), 3), min(filter(tab, scheme == "lag", hansen_df > 0L)$n_inst),
  num(max(unc_p$hansen_p), 2), max(unc_p$n_inst)),

"## Turning the dial by degrees",
sprintf(paste(
  "The %s rules above are all-or-nothing. %s more ways are not single rules but",
  "families with a knob on them. With the knob I can ask what happens part of the way.",
  "The results are in output/aggregation_families.txt and",
  "output/aggregation_families.png."),
  spell(length(AGG_SCHEMES)), sentence_case(spell(length(AGG_FAMILIES)))),
paste(
  "The first family groups the years into blocks and collapses the columns inside each",
  "block. A block size of one is then the uncollapsed set, and a block as wide as the",
  "panel is Roodman's rule. The second family keeps one instrument per year, but it",
  "multiplies each older lag by a fading factor before the addition. A factor of one is",
  "then the plain sum. Both families must match the fixed rules exactly at their",
  "endpoints. If they do not, the script stops."),
"Neither family did what I expected, and they surprised me in opposite directions.",
sprintf(paste(
  "Blocking is not a smooth dial at all. The estimates bounce and travel well outside",
  "the two endpoints they connect. For %s they run from %.3f to %.3f, while the",
  "endpoints themselves are only %.3f and %.3f."),
  fh,
  min(filter(fam, panel == fh, family == "Blocking the years")$income),
  max(filter(fam, panel == fh, family == "Blocking the years")$income),
  filter(fam, panel == fh, family == "Blocking the years", setting == 1)$income,
  filter(fam, panel == fh, family == "Blocking the years",
         setting == max(setting))$income),
paste(
  "Even the instrument count refuses to fall for wider blocks. Block sizes of four and",
  "five both give three blocks. But the wider blocks reach back to years that have",
  "deeper lags on offer, so the count goes up rather than down. Block size is therefore",
  "not a measure of the amount of collapsing, and that is why the line looks like",
  "noise."),
sprintf(paste(
  "The fading family does almost nothing. Across the whole range of factors the estimate",
  "moves by %.3f for %s and %.3f for %s. That result is useful rather than",
  "disappointing. It says that the answer from the collapse by period is not an artifact",
  "of equal weight on every lag. A tenth of the weight on the older lags barely moves",
  "the estimate."),
  diff(range(filter(fam, panel == fh, family == "Fading out older lags")$income)), fh,
  diff(range(filter(fam, panel == po, family == "Fading out older lags")$income)), po),

"## Roodman's footnote 7: random subsets",
sprintf(paste(
  "The rules above all pick the instrument subsets deliberately. Roodman's footnote 7",
  "suggests a random pick instead. He asks how the coefficient and the overidentification",
  "p-value move as the count grows. I drew %d random subsets of the uncollapsed",
  "lagged-level columns at each of several sizes, under the paper's instrument design. I",
  "then refit the model."), SUBSET_DRAWS),
paste(
  "The coefficient comes from the one-step estimator, so it is comparable with the rest",
  "of this file. The Hansen test needs the two-step weight matrix, so that column is",
  "two-step. The results are in output/aggregation_subsets.txt,",
  "output/aggregation_subsets.csv and output/aggregation_subsets.png."),
sprintf(paste(
  "The first result shows how much of the answer comes from the analyst's choice rather",
  "than from the data. With only %d of the %d available lagged levels in play, the %s",
  "draws run from %.3f to %.3f. So on this data and this specification, the sign of the",
  "effect is not fixed at all. The sign depends on which instruments the draw contains.",
  "The spread then narrows as the count grows. The standard deviation falls from %.3f",
  "at %d instruments to %.3f at %d for %s, and from %.3f to %.3f for %s."),
  sw(fh, "s_small"), st(fh, "n_gmm"), fh, sw(fh, "lo_small"), sw(fh, "hi_small"),
  sw(fh, "small"), sw(fh, "s_small"), sw(fh, "big"), sw(fh, "s_big"), fh,
  sw(po, "small"), sw(po, "big"), po),
sprintf(paste(
  "The median moves with the spread, and it lands in the right place. At the widest",
  "draws it is %.3f for %s and %.3f for %s, against %.3f and %.3f for the full",
  "uncollapsed set. The script checks that. The median reported standard error also",
  "falls from %.3f to %.3f for %s. So the estimate looks more precise as the count",
  "grows. It converges on the uncollapsed answer that drifts toward fixed effects, which",
  "is the conclusion the deliberate rules reach by a different route."),
  sw(fh, "med_big"), fh, sw(po, "med_big"), po, st(fh, "income"), st(po, "income"),
  sw(fh, "se_small"), sw(fh, "se_big"), fh),
sprintf(paste(
  "The J statistic is the messier half, and I do not want to read more into it than it",
  "supports. The Hansen test rejects a share of the draws at the %d%% level. That share",
  "runs %s for %s across the %s sizes, and %s for %s. Neither series is monotone, and",
  "the two do not agree with each other. So this is not a clean demonstration of a loss",
  "of power in the test."),
  as.integer(CI_LEVEL * 100), rej_series(fh), fh,
  spell(n_distinct(sub$size)), rej_series(po), po),
sprintf(paste(
  "Part of the reason is itself the point of the paper. The two-step weight matrix",
  "behind the test has one row and column per instrument. About %d countries supply the",
  "estimate of that matrix. So by the widest draws the test leans on a matrix that the",
  "data cannot support. The safe reading is the one Roodman gives. A high Hansen p-value",
  "on a large instrument set is not evidence of anything."),
  filter(tab, panel == fh)$countries[1]),
sprintf(paste(
  "The diamonds on the figure put the deliberate rules on the same axes, which is the",
  "comparison worth having. Roodman's rule uses %d instruments and gives %.3f for %s.",
  "Random draws of about that size (%d instruments) run from %.3f to %.3f. His rule",
  "therefore sits %s that range. The range is wide enough that a position inside it is",
  "not much of a recommendation on its own. For %s the same comparison is %.3f against a",
  "range of %.3f to %.3f."),
  near(fh)$n_ref, near(fh)$est, fh, near(fh)$n_draw, near(fh)$lo, near(fh)$hi,
  ifelse(near(fh)$inside, "inside", "outside"),
  po, near(po)$est, near(po)$lo, near(po)$hi),
paste(
  "The deliberate rule does not give a different answer at a given count. It gives",
  "stability as the lag window widens. The figure earlier in this file shows that",
  "stability, and the random draws cannot."),
paste(
  sprintf(paste("I did not implement the other suggestion in footnote 6, the",
                "vector-autoregression restriction of %s."), CITE_ARELLANO_OPT),
  "It needs a first-stage model of the instruments themselves rather than a rule for",
  "the collapse of the columns. So it does not fit into the same comparison. Roodman",
  "himself notes that it did not enter common practice."),

"## A note on reading the table",
paste(
  "The second column of the table uses every available lag, which is past the right-hand",
  sprintf("edge of the figure. The figure stops at a window of %d.", max(LAG_WINDOW)),
  "The instrument counts in the table include only the lagged levels. They do not",
  "include the year dummies that also sit in the instrument set."))

cat("Collapse comparison written.\n")
print(filter(tab, is.infinite(lag_max)) |>
        transmute(panel, design = substr(design, 1, 9), scheme,
                  income = round(income, 3), se = round(income_se, 3), n_inst),
      n = 40)
print(subsum, n = 20)
