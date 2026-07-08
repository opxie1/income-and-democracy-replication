source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(plm))

data("EmplUK", package = "plm")
.bench <- pgmm(log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) +
                 lag(log(output), 0:1) | lag(log(emp), 2:99),
               data = EmplUK, effect = "twoways", model = "twostep")
stopifnot(abs(coef(.bench)["lag(log(emp), 1:2)1"] - 0.4742) < 0.001,
          abs(coef(.bench)["lag(log(emp), 1:2)2"] + 0.0530) < 0.001)
cat("pgmm reproduces the Arellano-Bond (1991) benchmark.\n")

ladder <- function(dep, inc, panel_label) {
  d <- add_lags(read_panel(FILE_P5), c(dep, inc, "year"), 1:2)
  d$Ldep  <- d[[paste0(dep, "_l1")]]
  d$Linc  <- d[[paste0(inc, "_l1")]]
  d$y     <- d[[dep]] - d[[paste0(dep, "_l1")]]
  d$dLdep <- d[[paste0(dep, "_l1")]] - d[[paste0(dep, "_l2")]]
  d$dLinc <- d[[paste0(inc, "_l1")]] - d[[paste0(inc, "_l2")]]
  d$L2dep <- d[[paste0(dep, "_l2")]]
  d$L2inc <- d[[paste0(inc, "_l2")]]
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
  gmm_row <- function(label, model, transformation) {
    m <- pgmm(gmm_formula(dep, inc, "2:4"), data = pd, effect = "twoways",
              model = model, transformation = transformation, collapse = TRUE)
    co <- pgmm_co(m)
    ic <- co_row(co, inc); dc <- co_row(co, dep)
    add(label, ic[1], ic[2], dc[1], dc[2],
        pgmm_countries(m), ncol(m$W[[1]]),
        tryCatch(mtest(m, 1)$p.value, error = function(e) NA_real_),
        tryCatch(mtest(m, 2)$p.value, error = function(e) NA_real_),
        tryCatch(sargan(m)$p.value, error = function(e) NA_real_))
  }
  gmm_row("Arellano-Bond, difference GMM (one-step)", "onestep", "d")
  gmm_row("Arellano-Bond, difference GMM (two-step)", "twostep", "d")
  gmm_row("Blundell-Bond, system GMM (one-step)",     "onestep", "ld")
  gmm_row("Blundell-Bond, system GMM (two-step)",     "twostep", "ld")

  bind_rows(rows)
}

alt <- bind_rows(lapply(MEASURES, function(ms) ladder(ms$dep, ms$inc, ms$label)))
write_csv(alt, file.path(PATH_OUTPUT, "alternatives.csv"))

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

md <- c(
"# Alternative estimators for the income effect",
"",
"The paper asks whether a country growing richer ends up more democratic. A plain",
"correlation can fool you, so the authors use methods that take out the steady",
"differences between countries and the fact that democracy changes slowly. Here I",
"estimate that income effect several different ways and put them side by side.",
"This is the follow-up Professor Torgovitsky asked for. All of these use the",
"five-year sample, the same one behind Tables 2 and 3.",
"",
"## The methods, in plain terms",
"",
"- Pooled regression: ignores the steady differences between countries, so it",
"  reads too high.",
"- Fixed effects: takes out anything about a country that stays the same over",
"  time. With only a few time periods and a lagged outcome, it still has a known",
"  bias.",
"- Anderson-Hsiao: looks at changes instead of levels to cancel the steady",
"  differences, then uses values from two periods earlier as instruments.",
"- Arellano-Bond (difference GMM): the same change-based idea, using a set of",
"  earlier values as instruments. This is the method in the paper's GMM columns.",
"- Blundell-Bond (system GMM): adds a second set of conditions, in levels, on top",
"  of Arellano-Bond. It is more precise when those conditions hold, but they are",
"  an extra assumption.",
"",
"\"One-step\" and \"two-step\" are two ways of weighting these GMM estimators. The",
"two-step standard errors use the Windmeijer correction. To keep the number of",
"instruments small, so the tests below mean something, the GMM rows use a small",
"instrument set: collapsed, and limited to lags two through four.",
"",
"## How to read the diagnostics",
"",
"- The AR(1) p should be small and the AR(2) p should be large. That is the",
"  pattern you want, and it holds throughout.",
"- The overidentification p (Hansen's J test) should not be small. A small",
"  value warns that the instruments may not all be valid.",
"",
"## Results", "")
for (pl in unique(alt$panel)) md <- c(md, paste0("### ", pl), "", md_tab(filter(alt, panel == pl)), "")
md <- c(md,
"## What the comparison shows",
"",
"The change-based methods agree. Once the steady differences and the slow movement",
"of democracy are accounted for, the effect of income on democracy is small, and",
"in the instrumental-variables estimates it is negative rather than positive. This",
"is the paper's main finding, and it shows up in fixed effects, Anderson-Hsiao,",
"and Arellano-Bond alike. The two Arellano-Bond versions here, the one I wrote by",
"hand for the replication and the one from the plm package, land in the same",
"place.",
"",
"System GMM is the exception. Adding the level conditions pushes the income",
"coefficient to a small positive, statistically significant value for both",
"democracy measures. That reversal leans entirely on the extra conditions, and",
"the two measures treat them differently. For Polity the overidentification",
"test rejects them outright. For Freedom House the test raises no objection, so",
"there the positive estimate stands or falls with an assumption the paper",
"argues democracies in transition are unlikely to satisfy. Either way, the",
"estimates built on changes alone are all small or negative; the positive,",
"significant number appears only once the added conditions come in.",
"",
"So the paper's conclusion holds up under the change-based methods. The only way",
"to get a positive income effect back is to assume the extra system-GMM conditions",
"hold.",
"",
"## Checks",
"",
"I checked these numbers two ways. The GMM engine (plm's pgmm) reproduces the",
"textbook Arellano-Bond (1991) employment results exactly, and the run stops if",
"it ever fails to. Separately, R/11_crosscheck.R runs difference and system GMM",
"through an independent package, pdynmc, with its own uncollapsed instrument",
"set, and it stops with an error unless it finds the same picture as the tables",
"above: income negative under difference GMM and positive under system GMM.",
"That script needs the pdynmc package; run it with Rscript R/11_crosscheck.R.")
writeLines(md, file.path(PATH_DOCS, "alternatives.md"))

cat("Alternatives written. Income coefficient by estimator:\n")
print(transmute(alt, panel, estimator, income = round(income, 3), income_se = round(income_se, 3)))
