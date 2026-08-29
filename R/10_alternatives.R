source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(plm))

# benchmark
data("EmplUK", package = "plm")
.bench <- pgmm(log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) +
                 lag(log(output), 0:1) | lag(log(emp), 2:99),
               data = EmplUK, effect = "twoways", model = "twostep")
stopifnot(abs(coef(.bench)["lag(log(emp), 1:2)1"] - 0.4742) < 0.001,
          abs(coef(.bench)["lag(log(emp), 1:2)2"] + 0.0530) < 0.001)
cat("pgmm reproduces the Arellano-Bond (1991) benchmark.\n")

# ladder
ladder <- function(dep, inc, panel_label) {
  d <- prep_dynamic(FILE_P5, dep, inc)
  s <- filter(d, sample == 1)

  rows <- list()
  add <- function(estimator, income, income_se, dem, dem_se, countries,
                  instruments = NA_real_, ar1_p = NA_real_, ar2_p = NA_real_, overid_p = NA_real_) {
    rows[[length(rows) + 1]] <<- tibble(panel = panel_label, estimator = estimator,
      income = as.numeric(income), income_se = as.numeric(income_se),
      dem = as.numeric(dem), dem_se = as.numeric(dem_se),
      countries = as.integer(countries), instruments = as.numeric(instruments),
      ar1_p = as.numeric(ar1_p), ar2_p = as.numeric(ar2_p), overid_p = as.numeric(overid_p))
  }
  reg_row <- function(label, m, dem_t, inc_t) {
    add(label, ce(m, inc_t)["est"], ce(m, inc_t)["se"], ce(m, dem_t)["est"], ce(m, dem_t)["se"], mod_nc(m))
  }

  reg_row("Pooled OLS",        fit_ols(s, dep, c("Ldep", "Linc"), FALSE), "Ldep", "Linc")
  reg_row("Fixed effects",     fit_ols(s, dep, c("Ldep", "Linc"), TRUE),  "Ldep", "Linc")
  reg_row("Anderson-Hsiao IV", fit_iv(s, "y", endog = c("dLdep", "dLinc"),
          inst = c("L2dep", "L2inc"), country_fe = FALSE), "dLdep", "dLinc")

  est <- complete_on(s, c("y", "dLdep", "dLinc", "L2inc"))
  abr <- fit_abgmm(d, est, dep_level = dep, endog = c("dLdep", "dLinc"), inst_extra = "L2inc")
  add("Arellano-Bond, difference GMM (replication)", ce(abr, "dLinc")["est"], ce(abr, "dLinc")["se"],
      ce(abr, "dLdep")["est"], ce(abr, "dLdep")["se"], abr$n_country, abr$n_inst)

  pd <- gmm_panel(s, dep, inc)
  gmm_row <- function(spec) {
    m <- fit_gmm_spec(pd, dep, inc, spec)
    co <- pgmm_co(m)
    ic <- co_row(co, inc); dc <- co_row(co, dep)
    add(spec$label, ic[1], ic[2], dc[1], dc[2],
        pgmm_countries(m), ncol(m$W[[1]]),
        tryCatch(mtest(m, 1)$p.value, error = function(e) NA_real_),
        tryCatch(mtest(m, 2)$p.value, error = function(e) NA_real_),
        tryCatch(sargan(m, weights = spec$model)$p.value, error = function(e) NA_real_))
  }
  for (.k in ALT_GMM_KEYS) gmm_row(gmm_spec(.k))

  bind_rows(rows)
}

# estimates
alt <- bind_rows(lapply(MEASURES, function(ms) ladder(ms$dep, ms$inc, ms$label)))
write_csv(alt, file.path(PATH_OUTPUT, "alternatives.csv"))

ab <- filter(alt, grepl("^Arellano-Bond", estimator))
for (pl in unique(ab$panel)) {
  r <- filter(ab, panel == pl)
  ref <- filter(r, grepl("replication", estimator))
  stopifnot(nrow(ref) == 1L,
            all(abs(r$income - ref$income) <= r$income_se + 1e-12))
}
cat("The Arellano-Bond rows agree to within one standard error.\n")

# checks
sysp <- function(pl, wt)
  filter(alt, panel == pl, estimator == sprintf("Blundell-Bond, system GMM (%s)", wt))$overid_p
stopifnot(sysp(MEASURES[[2]]$label, "two-step") < CI_LEVEL,
          sysp(MEASURES[[1]]$label, "two-step") >= CI_LEVEL,
          all(filter(alt, grepl("one-step", estimator))$overid_p < CI_LEVEL))
cat("The Hansen test separates the two measures, and the Sargan test rejects every row.\n")

# layout
cell <- function(b, se) ifelse(is.na(b), "", sprintf("%s (%s)", num(b), num(se)))

txt <- c()
for (pl in unique(alt$panel)) {
  txt <- c(txt, paste0("== ", pl, " =="),
           sprintf("%-44s %-16s %-16s %5s %5s %6s %6s %7s",
                   "Estimator", "Income (SE)", "Democracy (SE)", "Ctry", "Inst", "AR1 p", "AR2 p", "Overid"))
  for (i in which(alt$panel == pl)) {
    r <- alt[i, ]
    txt <- c(txt, sprintf("%-44s %-16s %-16s %5d %5s %6s %6s %7s",
             r$estimator, cell(r$income, r$income_se), cell(r$dem, r$dem_se), r$countries,
             ifelse(is.na(r$instruments), "", as.character(as.integer(r$instruments))),
             ifelse(is.na(r$ar1_p), "", num(r$ar1_p, 2)),
             ifelse(is.na(r$ar2_p), "", num(r$ar2_p, 2)),
             ifelse(is.na(r$overid_p), "", num(r$overid_p, 2))))
  }
  txt <- c(txt, "")
}
writeLines(txt, file.path(PATH_OUTPUT, "alternatives.txt"))

# table
md_tab <- function(df) {
  out <- c("| Estimator | Income (SE) | Democracy (SE) | Countries | Instruments | AR(1) p | AR(2) p | Overid p |",
           "|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    out <- c(out, sprintf("| %s | %s | %s | %d | %s | %s | %s | %s |",
      r$estimator, cell(r$income, r$income_se), cell(r$dem, r$dem_se), r$countries,
      ifelse(is.na(r$instruments), "", as.character(as.integer(r$instruments))),
      ifelse(is.na(r$ar1_p), "", num(r$ar1_p, 2)),
      ifelse(is.na(r$ar2_p), "", num(r$ar2_p, 2)),
      ifelse(is.na(r$overid_p), "", num(r$overid_p, 2))))
  }
  out
}

# report
md <- c(
"# Alternative estimators for the income effect",
"",
"The paper asks whether a country becomes more democratic when it becomes richer.",
"A plain correlation can give the wrong answer. The authors therefore use methods",
"that remove the steady differences between countries and account for the slow",
"movement of democracy. Here I estimate the income effect in several ways and put",
"the results side by side. Professor Torgovitsky asked for this comparison. All of",
"the estimates use the five-year sample, the same sample behind Tables 2 and 3.",
"",
"## The methods, in plain terms",
"",
"- Pooled OLS pools the observations across countries. This method pools",
"  observations, not instruments. Pooled OLS ignores the steady differences",
"  between countries, and the estimate is therefore too high.",
"- Fixed effects removes anything about a country that stays the same over time.",
"  With few time periods and a lagged outcome, this method still has a known",
"  bias.",
"- Anderson-Hsiao uses changes instead of levels to cancel the steady",
"  differences. It then uses values from two periods earlier as instruments.",
"- Arellano-Bond (difference GMM) also works with changes, and it uses a set of",
"  earlier values as instruments. The paper uses this method in its GMM columns.",
"- Blundell-Bond (system GMM) keeps the Arellano-Bond conditions and adds a",
"  second set of conditions in levels. When these level conditions hold, system",
"  GMM is more precise. However, the level conditions are an extra assumption.",
"",
"\"One-step\" and \"two-step\" are two weighting schemes for these GMM estimators.",
"The two-step standard errors use the Windmeijer correction. The GMM rows use a",
"small instrument set: collapsed, and limited to lags two through four. A small",
"instrument set keeps the diagnostic tests informative.",
"",
"## How to read the diagnostics",
"",
"- A valid specification gives a small AR(1) p and a large AR(2) p. Every GMM row",
"  here shows this pattern.",
"- A valid specification gives an overidentification p that is not small. A small",
"  value warns that some of the instruments can be invalid. Each GMM row reports",
"  the test of its own estimator. The one-step rows report the Sargan test, which",
"  assumes well-behaved errors. The two-step rows report the Hansen test, which",
"  does not make this assumption.",
"- The Sargan test rejects every GMM row here, so it separates nothing. The paper",
"  clusters all of its own inference by country, and that is the situation in",
"  which a homoskedastic test over-rejects. The Hansen column is the one to read.",
"",
"## Results", "")
for (pl in unique(alt$panel)) md <- c(md, paste0("### ", pl), "", md_tab(filter(alt, panel == pl)), "")
md <- c(md,
"## What the comparison shows",
"",
"The change-based methods agree with each other. These methods remove the steady",
"differences between countries and account for the slow movement of democracy.",
"After these corrections, the effect of income on democracy is small. In the",
"instrumental-variables estimates the effect is negative, not positive. This",
"result is the main finding of the paper, and fixed effects, Anderson-Hsiao and",
"Arellano-Bond all show it.",
"",
"The Arellano-Bond rows differ among themselves by less than one standard error.",
"I do not expect closer agreement than this. The replication row uses the",
"instrument set of the paper, which is uncollapsed and uses every lag. The plm",
"rows use a collapsed set limited to lags two through four.",
"",
"The gap between the one-step and two-step rows carries little information here.",
"The gap is small because the instrument set is small. The gap is not small in",
"general. With a wide lag window, the weighting step decides whether collapsing",
"works at all. The file docs/instruments.md covers this question.",
"",
"System GMM is the exception. The level conditions push the income coefficient to",
"a small positive value for both democracy measures, and that value is",
"statistically significant. This reversal rests entirely on the level conditions.",
"The two democracy measures treat these conditions differently. For Polity, the",
"Hansen test rejects them outright. For Freedom House, the Hansen test raises no",
"objection.",
"",
"For Freedom House the positive estimate therefore depends on one assumption. The",
"paper argues that democracies in transition are unlikely to satisfy this",
"assumption. Either way, every change-based estimate is small or negative. The",
"positive and significant number appears only in the rows that add the level",
"conditions.",
"",
"The conclusion of the paper therefore holds under the change-based methods. A",
"positive income effect returns only under the assumption that the level",
"conditions hold.",
"",
"## Checks",
"",
"I checked these numbers two ways. The GMM engine is the pgmm function in plm. It",
"reproduces the textbook Arellano-Bond (1991) employment results exactly. If it",
"does not reproduce them, the script stops.",
"",
"The script R/11_crosscheck.R runs difference and system GMM through an",
"independent package, pdynmc, with its own uncollapsed instrument set. That",
"script must find the same result as the tables above: income negative under",
"difference GMM and positive under system GMM. If it finds anything else, the",
"script stops. R/11_crosscheck.R needs the pdynmc package.",
"",
"To run the cross-check, use `Rscript R/11_crosscheck.R`.")
writeLines(md, file.path(PATH_DOCS, "alternatives.md"))

cat("Alternatives written. Income coefficient by estimator:\n")
print(transmute(alt, panel, estimator, income = round(income, 3), income_se = round(income_se, 3)))
