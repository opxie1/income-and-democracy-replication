source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(plm))

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
  add <- function(estimator, income, income_se, dem, dem_se, n, countries,
                  instruments = NA_real_, hansen_p = NA_real_, ar2_p = NA_real_) {
    rows[[length(rows) + 1]] <<- tibble(panel = panel_label, estimator = estimator,
      income = as.numeric(income), income_se = as.numeric(income_se),
      dem = as.numeric(dem), dem_se = as.numeric(dem_se),
      n = as.integer(n), countries = as.integer(countries),
      instruments = as.numeric(instruments),
      hansen_p = as.numeric(hansen_p), ar2_p = as.numeric(ar2_p))
  }
  reg_row <- function(label, m, dem_t, inc_t) {
    add(label, ce(m, inc_t)["est"], ce(m, inc_t)["se"], ce(m, dem_t)["est"], ce(m, dem_t)["se"],
        mod_nobs(m), mod_nc(m))
  }

  reg_row("Pooled OLS",        fit_ols(s, dep, c("Ldep", "Linc"), FALSE), "Ldep", "Linc")
  reg_row("Fixed effects",     fit_ols(s, dep, c("Ldep", "Linc"), TRUE),  "Ldep", "Linc")
  reg_row("Anderson-Hsiao IV", fit_iv(s, "y", endog = c("dLdep", "dLinc"),
          inst = c("L2dep", "L2inc"), country_fe = FALSE), "dLdep", "dLinc")

  est <- complete_on(s, c("y", "dLdep", "dLinc", "L2inc"))
  abr <- fit_abgmm(d, est, dep_level = dep, endog = c("dLdep", "dLinc"), inst_extra = "L2inc")
  add("Arellano-Bond, difference GMM (replication)", ce(abr, "dLinc")["est"], ce(abr, "dLinc")["se"],
      ce(abr, "dLdep")["est"], ce(abr, "dLdep")["se"], abr$nobs, abr$n_country, abr$n_inst)

  pd <- pdata.frame(s[stats::complete.cases(s[, c(dep, inc)]), ], index = c("code", "year_numeric"))
  f <- as.formula(sprintf("%s ~ lag(%s, 1) + lag(%s, 1) | lag(%s, 2:99) + lag(%s, 2:99)",
                          dep, dep, inc, dep, inc))
  gmm_row <- function(label, model, transformation) {
    m <- pgmm(f, data = pd, effect = "twoways", model = model, transformation = transformation)
    co <- summary(m, robust = TRUE)$coefficients
    ir <- grep(paste0("lag\\(", inc), rownames(co))[1]
    dr <- grep(paste0("lag\\(", dep), rownames(co))[1]
    add(label, co[ir, 1], co[ir, 2], co[dr, 1], co[dr, 2],
        sum(lengths(m$residuals)), length(m$residuals), dim(m$W[[1]])[2],
        tryCatch(sargan(m)$p.value, error = function(e) NA_real_),
        tryCatch(mtest(m, order = 2)$p.value, error = function(e) NA_real_))
  }
  gmm_row("Arellano-Bond, difference GMM (one-step)", "onestep", "d")
  gmm_row("Arellano-Bond, difference GMM (two-step)", "twostep", "d")
  gmm_row("Blundell-Bond, system GMM (one-step)",     "onestep", "ld")
  gmm_row("Blundell-Bond, system GMM (two-step)",     "twostep", "ld")

  bind_rows(rows)
}

alt <- bind_rows(ladder("fhpolrigaug", "lrgdpch", "Freedom House"),
                 ladder("polity4", "lrgdpch", "Polity"))
write_csv(alt, file.path(PATH_OUTPUT, "alternatives.csv"))

num <- function(x, d = 3) ifelse(is.na(x), "", sprintf(paste0("%.", d, "f"), x))
cell <- function(b, se) ifelse(is.na(b), "", sprintf("%s (%s)", num(b), num(se)))
txt <- c()
for (pl in unique(alt$panel)) {
  txt <- c(txt, paste0("== ", pl, " =="),
           sprintf("%-44s %-16s %-16s %5s %5s %8s %6s",
                   "Estimator", "Income (SE)", "Democracy (SE)", "Ctry", "Inst", "Hansen p", "AR2 p"))
  sub <- filter(alt, panel == pl)
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    txt <- c(txt, sprintf("%-44s %-16s %-16s %5d %5s %8s %6s",
             r$estimator, cell(r$income, r$income_se), cell(r$dem, r$dem_se),
             r$countries,
             ifelse(is.na(r$instruments), "", as.character(as.integer(r$instruments))),
             ifelse(is.na(r$hansen_p), "", num(r$hansen_p, 2)),
             ifelse(is.na(r$ar2_p), "", num(r$ar2_p, 2))))
  }
  txt <- c(txt, "")
}
writeLines(txt, file.path(PATH_OUTPUT, "alternatives.txt"))

md_tab <- function(df) {
  out <- c("| Estimator | Income (SE) | Democracy (SE) | Countries | Instruments | Hansen p | AR(2) p |",
           "|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    out <- c(out, sprintf("| %s | %s | %s | %d | %s | %s | %s |",
      r$estimator, cell(r$income, r$income_se), cell(r$dem, r$dem_se), r$countries,
      ifelse(is.na(r$instruments), "", as.character(as.integer(r$instruments))),
      ifelse(is.na(r$hansen_p), "", num(r$hansen_p, 2)),
      ifelse(is.na(r$ar2_p), "", num(r$ar2_p, 2))))
  }
  out
}

md <- c(
"# Alternative estimators for the income effect",
"",
"The paper asks whether a country growing richer leads it to become more",
"democratic. A plain correlation is misleading, so the authors rely on panel",
"methods that strip out fixed differences between countries and the tendency of",
"democracy to persist. This note re-estimates that income effect with several",
"methods and lines them up, which is the follow-up Professor Torgovitsky asked",
"for. The estimates are for the baseline five-year sample, the same one behind",
"Tables 2 and 3.",
"",
"## The methods, in plain terms",
"",
"- Pooled regression: ignores fixed country differences, and reads high.",
"- Fixed effects: removes things about a country that do not change over time. In",
"  a short panel with a lagged outcome it still carries a known bias.",
"- Anderson-Hsiao: looks at changes rather than levels to cancel the fixed",
"  differences, then uses values from two periods earlier as instruments.",
"- Arellano-Bond (difference GMM): the same change-based idea, using the whole",
"  set of earlier values as instruments at once. This is the method in the",
"  paper's own GMM columns.",
"- Blundell-Bond (system GMM): adds a second set of conditions, in levels, on top",
"  of Arellano-Bond. It is more precise when those conditions hold, but they are",
"  an extra assumption.",
"",
"\"One-step\" and \"two-step\" are two ways of weighting the GMM estimators; the",
"two-step standard errors include the Windmeijer small-sample correction.",
"",
"## Results", "")
for (pl in unique(alt$panel)) md <- c(md, paste0("### ", pl), "", md_tab(filter(alt, panel == pl)), "")
md <- c(md,
"## What the comparison shows",
"",
"The change-based methods agree. Once fixed differences and persistence are",
"accounted for, the effect of income on democracy is small, and in the",
"instrumental-variables estimates it is negative rather than positive. That is",
"the paper's central finding, and it holds across fixed effects, Anderson-Hsiao,",
"and Arellano-Bond. The two Arellano-Bond implementations here, the hand-built",
"one used for the replication and the one from the plm package, land in the same",
"place, which is a useful check that the replication's estimator is doing what it",
"should.",
"",
"System GMM is the exception. Adding the level conditions pushes the income",
"coefficient to a small positive and statistically significant value for both",
"democracy measures. That looks like a reversal, but it rests entirely on the",
"extra assumption behind those level conditions, which is the kind of assumption",
"the paper is wary of. The overidentification (Hansen) test is borderline for",
"several of the one-step GMM runs, and the instrument counts are high, so the GMM",
"numbers deserve caution.",
"",
"So the paper's conclusion holds up under the change-based alternatives. The only",
"way to bring back a positive income effect is to assume the extra system-GMM",
"conditions hold.")
writeLines(md, file.path(PATH_DOCS, "alternatives.md"))

cat("Alternatives written. Income coefficient by estimator:\n")
print(transmute(alt, panel, estimator, income = round(income, 3), income_se = round(income_se, 3)))
