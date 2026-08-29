# report
write_mc_report <- function(sm, big, shrink, cals) {
  # table
  PARAM_LABEL <- c(alpha = "Persistence of democracy", beta = "Effect of income")
  txt <- c()
  for (mv in names(cals)) for (kd in MC_CALIB) for (dg in MC_DESIGNS) {
    b <- filter(sm, measure == mv, calibration == kd, design == dg)
    txt <- c(txt, sprintf("== %s | %s | %s ==", mv, kd, dg),
             sprintf("true persistence %.3f, true income effect %.3f",
                     b$alpha_true[1], b$beta_true[1]),
             sprintf("%-48s %8s %8s %8s %8s %8s", "Estimator", "bias", "sd", "mean se",
                     "coverage", "reject 0"))
    for (pm in names(PARAM_LABEL)) {
      txt <- c(txt, sprintf("-- %s --", PARAM_LABEL[[pm]]))
      r <- filter(b, param == pm) |> arrange(match(key, LADDER_ORDER))
      for (i in seq_len(nrow(r)))
        txt <- c(txt, sprintf("%-48s %8s %8s %8s %8s %8s", r$estimator[i],
          num(r$bias[i]), num(r$sd_est[i]), num(r$mean_se[i]),
          num(r$coverage[i], 2), num(r$reject0[i], 2)))
    }
    txt <- c(txt, "")
  }
  writeLines(txt, file.path(PATH_OUTPUT, "montecarlo.txt"))
  
  # figure
  plotdat <- sm |>
    mutate(param_label = unname(PARAM_LABEL[param]),
           short = ifelse(grepl("Mean", design), "stationary start", "start tied to effect"),
           estimator = fct_rev(estimator))
  mkfig <- function(v, xlab, ref, title, sub) {
    ggplot(plotdat, aes(.data[[v]], estimator, color = measure, shape = short)) +
      geom_vline(xintercept = ref, color = "grey40", linetype = "dashed") +
      geom_point(size = 2, alpha = 0.85,
                 position = position_dodge(width = 0.6)) +
      facet_grid(calibration ~ param_label, scales = "free_x",
                 labeller = label_wrap_gen(28)) +
      labs(x = xlab, y = NULL, color = NULL, shape = NULL,
           title = title, subtitle = sub) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "top", panel.grid.major.y = element_line(linewidth = 0.2))
  }
  ggsave(file.path(PATH_OUTPUT, "montecarlo_bias.png"),
         mkfig("bias", "Average error of the estimate", 0,
               "What each estimator gets wrong on data where the truth is known",
               sprintf("%d draws per design. The true income effect is zero.", MC_REPS)),
         width = 11, height = 7.5, dpi = 150)
  ggsave(file.path(PATH_OUTPUT, "montecarlo_coverage.png"),
         mkfig("coverage", "Share of draws whose 95% interval covers the truth", 1 - CI_LEVEL,
               "How often the reported confidence interval contains the true value",
               sprintf("%d draws per design. The dashed line is the nominal %d per cent.",
                       MC_REPS, as.integer((1 - CI_LEVEL) * 100))),
         width = 11, height = 7.5, dpi = 150)
  
  cat("Monte Carlo written.\n")
  pick <- function(mv, kd, dg, pm, k, col)
    filter(sm, measure == mv, calibration == unname(MC_CALIB[[kd]]),
           design == unname(MC_DESIGNS[[dg]]), param == pm, key == k)[[col]]
  fh <- names(cals)[1]; po <- names(cals)[2]
  cfh <- cals[[fh]]$cal$gmm; cpo <- cals[[po]]$cal$gmm
  
  base <- filter(sm, calibration == unname(MC_CALIB[["gmm"]]),
                 design == unname(MC_DESIGNS[["stationary"]]), param == "beta")
  worst_ols <- filter(base, key == "ols")
  best_rmse <- base |> filter(coverage >= 0.90) |>
    group_by(measure) |> slice_min(rmse, n = 1) |> ungroup()
  near_cov  <- base |> group_by(measure) |>
    slice_min(abs(coverage - (1 - CI_LEVEL)), n = 1) |> ungroup()
  sysfall <- filter(sm, param == "beta", key == "ld_two_coll",
                    calibration == unname(MC_CALIB[["gmm"]]))
  pull_chr <- function(d, col, mv) as.character(d[[col]][d$measure == mv])[1]
  pull_num <- function(d, col, mv) d[[col]][d$measure == mv][1]
  sysc <- function(mv, stat) sysfall$coverage[sysfall$measure == mv &
    (grepl("Mean", sysfall$design) == stat)][1]
  shrunk_n <- sum(shrink$shrunk, na.rm = TRUE)
  sh <- function(k, col) shrink[[col]][shrink$key == k][1]
  shr_pct <- function(k) 100 * (1 - abs(sh(k, "bias_big")) / abs(sh(k, "bias_small")))
  gmm_pct <- vapply(setdiff(shrink$key, c("ols", "fe", "abr", "ah")), shr_pct, numeric(1))
  bias_abr_fh <- pick(fh, "gmm", "stationary", "beta", "abr", "bias")
  bias_abr_po <- pick(po, "gmm", "stationary", "beta", "abr", "bias")
  
  write_doc("monte-carlo.md",
  "# A simulation with known answers",
  paste(
    "Professor Torgovitsky asked for a Monte Carlo study calibrated on this data. The",
    "idea is to build a dynamic panel that looks like the real one. It keeps the same",
    "number of countries and the same number of periods, and it takes persistence from",
    "the estimates here. In simulated data the true parameters are known, so the error of",
    "each estimator is measurable. The results give a guide to which estimates deserve",
    "trust in the real application."),
  sprintf(paste(
    "The numbers are in output/montecarlo.txt and output/montecarlo.csv. Every draw is in",
    "output/montecarlo_draws.csv.gz. The figures are output/montecarlo_bias.png and",
    "output/montecarlo_coverage.png. The code is R/16_montecarlo.R and R/_montecarlo.R.",
    "It runs %d draws for each of the %d designs."),
    MC_REPS, length(cals) * length(MC_CALIB) * length(MC_DESIGNS)),
  
  "## How the simulated data are built",
  sprintf(paste(
    "Each country gets a democracy equation and an income equation. Democracy depends on",
    "its own lag and on lagged income. Income depends on its own lag and on lagged",
    "democracy. Both equations carry a country effect, a period effect and a shock. The",
    "two shocks are drawn together, so a shock to one variable can move the other. The",
    "country effects are drawn together as well, with a correlation of %.2f for %s. That",
    "correlation is what makes rich countries democratic for fixed reasons. It is also",
    "what pooled OLS mistakes for an effect of income."),
    cfh$cor_fe, fh),
  sprintf(paste(
    "The simulated panel keeps the shape of the real one. It reproduces the observation",
    "pattern cell by cell, for each variable separately, so the estimators end up with",
    "the sample they have in the real data. For %s the estimation sample is %d",
    "observations on %d countries across %d periods. For %s it is %d observations on %d",
    "countries. The script checks these against the real counts and stops if they",
    "differ."),
    fh, cfh$gmm_obs, cfh$gmm_ctry, cfh$gmm_periods,
    po, cpo$gmm_obs, cpo$gmm_ctry),
  sprintf(paste(
    "The income equation is close to a random walk. Its persistence is %.3f for %s and",
    "%.3f for %s, both from difference GMM. That number matters more than it looks, and",
    "the results section returns to it."),
    cfh$rho, fh, cpo$rho, po),
  paste(
    "The true effect of income on democracy is zero in every design. That choice is",
    "deliberate. With a true zero, the share of draws in which an estimator reports a",
    "significant effect is its false-positive rate. It answers the question the paper",
    "asks, which is whether there is an effect at all."),
  
  "## The two persistence values",
  sprintf(paste(
    "Professor Torgovitsky asked for persistence set to what I estimate here. The trouble",
    "is that the estimates disagree. For %s, pooled OLS gives %.3f, fixed effects gives",
    "%.3f and difference GMM gives %.3f. Picking one of them presumes the answer that the",
    "simulation is meant to test."),
    fh, cals[[fh]]$target$ols, cals[[fh]]$target$fe, cals[[fh]]$alpha[["gmm"]]),
  sprintf(paste(
    "So the study runs two calibrations. The first sets persistence to the difference-GMM",
    "estimate, %.3f for %s and %.3f for %s. The second sets it to the value that makes",
    "the simulated fixed-effects estimate equal the real one, %.3f and %.3f. The second",
    "calibration assumes nothing about which estimator is right. It asks only what true",
    "persistence produces the fixed-effects number that the data show."),
    cals[[fh]]$alpha[["gmm"]], fh, cals[[po]]$alpha[["gmm"]], po,
    cals[[fh]]$alpha[["fe"]], cals[[po]]$alpha[["fe"]]),
  
  "## What the calibration matches, and what it misses",
  sprintf(paste(
    "The variance of the country effects is chosen so that the spread of countries within",
    "a period matches the real spread. For %s the target standard deviations are %.3f for",
    "democracy and %.3f for income."),
    fh, cfh$target_sd[1], cfh$target_sd[2]),
  paste(
    "One thing does not match, and it is worth stating. Pooled OLS is more biased in the",
    "simulation than in the data. A linear model with a country effect and a normal shock",
    "cannot match three things at once. Those three are the spread of countries, the",
    "fixed-effects estimate and the pooled OLS estimate. The simulation is therefore a",
    "slightly harsh world",
    "for the estimators that work in levels. Democracy in the real data also sits between",
    "0 and 1, and the simulated series has no such bound."),
  
  "## The two designs",
  sprintf(paste(
    "Each calibration runs under two starting conditions. In the first, countries start",
    "at their long-run mean, so the extra conditions that system GMM adds are valid. In",
    "the second, the starting point is tied to the country effect, so those extra",
    "conditions fail. Difference GMM stays valid under both. The second design is the one",
    "Roodman simulates, and docs/alternatives.md finds the same failure in the real data",
    "for %s."), po),
  
  "## What happens to the persistence estimate",
  sprintf(paste(
    "The classic results appear, which is the first sign that the simulation is built",
    "correctly. Under the difference-GMM calibration for %s the true persistence is %.3f.",
    "Pooled OLS averages %.3f and fixed effects averages %.3f. Their intervals cover the",
    "truth in %s and %s of draws."),
    fh, cals[[fh]]$alpha[["gmm"]],
    pick(fh, "gmm", "stationary", "alpha", "ols", "mean_est"),
    pick(fh, "gmm", "stationary", "alpha", "fe", "mean_est"),
    num(pick(fh, "gmm", "stationary", "alpha", "ols", "coverage"), 2),
    num(pick(fh, "gmm", "stationary", "alpha", "fe", "coverage"), 2)),
  sprintf(paste(
    "The collapsed difference GMM estimator is almost unbiased. Its average error is %.3f",
    "one-step and %.3f two-step. The price is noise, and the standard deviation across",
    "draws is %.3f. The uncollapsed estimator is the opposite. It is tighter, with a",
    "standard deviation of %.3f, but its average error is %.3f. Its intervals cover the",
    "truth only %s of the time. That is the cost of too many instruments, measured",
    "directly."),
    pick(fh, "gmm", "stationary", "alpha", "d_one_coll", "bias"),
    pick(fh, "gmm", "stationary", "alpha", "d_two_coll", "bias"),
    pick(fh, "gmm", "stationary", "alpha", "d_one_coll", "sd_est"),
    pick(fh, "gmm", "stationary", "alpha", "d_one_unc", "sd_est"),
    pick(fh, "gmm", "stationary", "alpha", "d_one_unc", "bias"),
    num(pick(fh, "gmm", "stationary", "alpha", "d_one_unc", "coverage"), 2)),
  
  "## What happens to the income effect",
  sprintf(paste(
    "The true effect is zero everywhere, so every number in this section is an error.",
    "Pooled OLS reports a positive effect of %.3f for %s and %.3f for %s. It calls that",
    "effect significant in %s and %s of draws, and its interval covers the true zero in",
    "%s and %s. A researcher who trusts pooled OLS on data of this shape finds an income",
    "effect that is not there, nearly every time."),
    pick(fh, "gmm", "stationary", "beta", "ols", "mean_est"), fh,
    pick(po, "gmm", "stationary", "beta", "ols", "mean_est"), po,
    num(pick(fh, "gmm", "stationary", "beta", "ols", "reject0"), 2),
    num(pick(po, "gmm", "stationary", "beta", "ols", "reject0"), 2),
    num(pick(fh, "gmm", "stationary", "beta", "ols", "coverage"), 2),
    num(pick(po, "gmm", "stationary", "beta", "ols", "coverage"), 2)),
  sprintf(paste(
    "The estimator in the GMM columns of the paper is difference GMM with uncollapsed",
    "instruments and every lag. Here it carries a negative bias of %.3f for %s and %.3f",
    "for %s, on data where the true effect is zero. The paper reports %.3f and %.3f in",
    "those columns. The bias therefore covers %.0f per cent of the published estimate for",
    "%s and %.0f per cent for %s. That does not show that the published numbers are only",
    "bias. It shows that a number of this kind arises on data of this shape when the true",
    "effect is zero."),
    bias_abr_fh, fh, bias_abr_po, po, cfh$beta_hat, cpo$beta_hat,
    100 * bias_abr_fh / cfh$beta_hat, fh, 100 * bias_abr_po / cpo$beta_hat, po),
  sprintf(paste(
    "I do not have a clean account of why the two measures differ this much. My first",
    "guess was the income process, which is close to a random walk and therefore gives",
    "weak instruments. That guess does not survive. Income is more persistent for %s",
    "(%.3f) than for %s (%.3f), so the weaker instruments belong to the measure with the",
    "smaller bias. The honest statement is that the bias is real and measure-specific,",
    "and that I cannot yet name its source."),
    po, cpo$rho, fh, cfh$rho),
  sprintf(paste(
    "Collapsing helps a little on this coefficient. The collapsed bias is %.3f against",
    "%.3f uncollapsed for %s. It costs a great deal of precision, and the standard",
    "deviation across draws rises from %.3f to %.3f. The gain from collapsing is far",
    "clearer on the persistence coefficient than on this one."),
    pick(fh, "gmm", "stationary", "beta", "d_one_coll", "bias"),
    pick(fh, "gmm", "stationary", "beta", "d_one_unc", "bias"), fh,
    pick(fh, "gmm", "stationary", "beta", "d_one_unc", "sd_est"),
    pick(fh, "gmm", "stationary", "beta", "d_one_coll", "sd_est")),
  
  "## Whether the confidence intervals are honest",
  sprintf(paste(
    "Coverage answers the question directly. A 95 per cent interval must contain the",
    "truth in 95 per cent of draws. Take the income effect, under the difference-GMM",
    "calibration and a stationary start. For %s the estimator closest to that target is",
    "%s at %s. Pooled OLS is the worst at %s."),
    fh, pull_chr(near_cov, "estimator", fh), num(pull_num(near_cov, "coverage", fh), 2),
    num(pick(fh, "gmm", "stationary", "beta", "ols", "coverage"), 2)),
  sprintf(paste(
    "The starting condition decides whether system GMM can be trusted. Take the two-step",
    "collapsed version. With a stationary start its intervals for the income effect cover",
    "%s for %s and %s for %s. With the start tied to the country effect the same",
    "intervals cover %s and %s. Difference GMM barely moves between the two designs. This",
    "is the warning that the overidentification test gives in the real data, and here the",
    "cost of ignoring it is measured."),
    num(sysc(fh, TRUE), 2), fh, num(sysc(po, TRUE), 2), po,
    num(sysc(fh, FALSE), 2), num(sysc(po, FALSE), 2)),
  
  "## What this says about the real estimates",
  sprintf(paste(
    "The answer depends on which coefficient is wanted, and that is the first thing to",
    "say. On persistence, pooled OLS and fixed effects are both hopeless. They miss by",
    "%.3f and %.3f for %s, and neither interval ever covers the truth. Collapsed",
    "difference GMM is the only estimator that is close to unbiased and honest about its",
    "own uncertainty at the same time."),
    pick(fh, "gmm", "stationary", "alpha", "ols", "bias"),
    pick(fh, "gmm", "stationary", "alpha", "fe", "bias"), fh),
  sprintf(paste(
    "On the income coefficient the two part company. Pooled OLS still fails, with",
    "coverage of %s. Fixed effects does not fail: its error is %.3f and its intervals",
    "cover the truth %s of the time. The bias of the fixed-effects estimator falls on the",
    "lagged dependent variable rather than on the other regressor. So a reader who cares",
    "only about the income coefficient loses less by using fixed effects than the",
    "persistence results suggest."),
    num(pick(fh, "gmm", "stationary", "beta", "ols", "coverage"), 2),
    pick(fh, "gmm", "stationary", "beta", "fe", "bias"),
    num(pick(fh, "gmm", "stationary", "beta", "fe", "coverage"), 2)),
  sprintf(paste(
    "The uncollapsed difference GMM column, which is the one the paper reports, is biased",
    "toward a negative income effect on data of this shape. Now take only the estimators",
    "whose intervals cover the truth at least 90 per cent of the time. Among those, the",
    "smallest root mean squared error on the income effect for %s belongs to %s at %.3f.",
    "Pooled OLS has a",
    "smaller root mean squared error still, at %.3f, and that is the trap in using that",
    "measure alone. A tight wrong answer beats a wide right one on root mean squared",
    "error, and coverage is what separates them."),
    fh, pull_chr(best_rmse, "estimator", fh), pull_num(best_rmse, "rmse", fh),
    pick(fh, "gmm", "stationary", "beta", "ols", "rmse")),
  paste(
    "The practical reading matches the weak-instrument work in docs/weak-instruments.md.",
    "The instrumental-variables columns of this paper are not informative about the sign",
    "of the income effect. The simulation adds a reason. On data with this many countries,",
    "this many periods and income this persistent, no estimator here separates a true zero",
    "from the effect the paper reports."),
  
  "## Checks",
  sprintf(paste(
    "Four checks run with the script. The estimator ladder is the same code that",
    "R/10_alternatives.R runs on the real data, so the simulated and the real columns",
    "cannot drift apart. The simulated estimation sample has to match the real one, and",
    "it does, at %d observations on %d countries for %s. Every estimator has to return",
    "an estimate in at least 90 per cent of draws. The bias has to shrink when the panel",
    "grows, and with %d times as many countries it falls for %d of the %d estimators. If",
    "a check fails, the script stops."),
    cfh$gmm_obs, cfh$gmm_ctry, fh, big$times[1], shrunk_n, nrow(shrink)),
  sprintf(paste(
    "That last check separates two kinds of error. The GMM errors fall by between %.0f",
    "and %.0f per cent as the panel grows, so they are finite-sample problems. Pooled",
    "OLS and fixed effects fall by about %.0f per cent, which is to say not at all,",
    "because neither is consistent here at any sample size. That is the difference",
    "between an estimator that needs a bigger panel and one that a bigger panel cannot",
    "save."),
    min(gmm_pct), max(gmm_pct), mean(c(shr_pct("ols"), shr_pct("fe")))),
  sprintf(paste(
    "One GMM row is slower than the rest, and it is worth naming. The estimator in the",
    "GMM columns of the paper falls only %.0f per cent, from %.3f to %.3f. Under the",
    "instrument design of the paper, income gets a single lagged level as its instrument,",
    "so that estimator has the least to work with. I read this as slow convergence rather",
    "than inconsistency, but the simulation here cannot separate the two."),
    shr_pct("abr"), sh("abr", "bias_small"), sh("abr", "bias_big")),
  paste(
    "One limitation stands out. The true income effect is zero in every design here, so",
    "these results measure bias and false positives. They do not measure power against a",
    "real effect. That is the natural next step."),
  paste(
    "A second limitation is one of scope. The ladder covers every estimator in",
    "docs/alternatives.md and the collapsed and uncollapsed variants from",
    "docs/instruments.md. It does not cover the two-stage least squares columns of",
    "Tables 5 and 6. Those columns need an instrument from outside the panel, the savings",
    "rate and trade-weighted world income, and the simulated data contain no such",
    "variable. Simulating one is a separate design and a separate question."))
  invisible(TRUE)
}
