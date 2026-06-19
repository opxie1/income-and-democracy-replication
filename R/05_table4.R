# ---------------------------------------------------------------------------
# 05_table4.R  -- Table 4: Freedom House robustness checks (five-year data).
# Balanced panel, dropping socialist countries, and adding demographic and
# education covariates. Odd columns are fixed-effects OLS; even columns are
# Arellano-Bond GMM. Median age and the four age-structure shares are included
# but not displayed (as in the paper).
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

d5 <- read_parquet(FILE_P5) |> arrange(code, year_numeric) |> mutate(period = year_numeric)
d5 <- add_lags(d5, c("fhpolrigaug", "lrgdpch", "lpop", "medage", "education",
                     "age_veryyoung", "age_young", "age_midage", "age_old", "year"), 1:2)
d5 <- d5 |> mutate(
  Ldep = fhpolrigaug_l1, Linc = lrgdpch_l1,
  Llpop = lpop_l1, Lmed = medage_l1, Ledu = education_l1,
  Lavy = age_veryyoung_l1, Lay = age_young_l1, Lam = age_midage_l1, Lao = age_old_l1,
  y = fhpolrigaug - fhpolrigaug_l1,
  dLdep = fhpolrigaug_l1 - fhpolrigaug_l2, dLinc = lrgdpch_l1 - lrgdpch_l2, L2inc = lrgdpch_l2,
  dLpop = lpop_l1 - lpop_l2, dLmed = medage_l1 - medage_l2, dLedu = education_l1 - education_l2,
  dLavy = age_veryyoung_l1 - age_veryyoung_l2, dLay = age_young_l1 - age_young_l2,
  dLam = age_midage_l1 - age_midage_l2, dLao = age_old_l1 - age_old_l2)

AGE <- c("Lavy", "Lay", "Lam", "Lao"); dAGE <- c("dLavy", "dLay", "dLam", "dLao")

rows <- list()
push <- function(column, row, value, se = NA_real_, type = "coef")
  rows[[length(rows) + 1]] <<- tibble(column = column, row = row,
    value = as.numeric(value), se = as.numeric(se), type = type)
ols_col <- function(col, m, terms) {
  labs <- c(Ldep = "Democracy_t-1", Linc = "Log GDP per capita_t-1",
            Llpop = "Log population_t-1", Ledu = "Education_t-1")
  for (t in terms) push(col, labs[[t]], ce(m, t)["est"], ce(m, t)["se"])
  push(col, "Observations", mod_nobs(m), type = "count")
  push(col, "Countries", mod_nc(m), type = "count")
  push(col, "R-squared", mod_r2(m), type = "r2")
}
gmm_col <- function(col, m, terms) {
  labs <- c(dLdep = "Democracy_t-1", dLinc = "Log GDP per capita_t-1",
            dLpop = "Log population_t-1", dLedu = "Education_t-1")
  for (t in terms) push(col, labs[[t]], ce(m, t)["est"], ce(m, t)["se"])
  push(col, "Observations", mod_nobs(m), type = "count")
  push(col, "Countries", mod_nc(m), type = "count")
}

# C1-C2: balanced panel 1970-2000
ols_col(1, fit_ols(filter(d5, samplebalancefe == 1), "fhpolrigaug", c("Ldep", "Linc"), TRUE),
        c("Ldep", "Linc"))
est2 <- complete_on(filter(d5, samplebalancegmm == 1), c("y", "dLdep", "dLinc", "L2inc"))
gmm_col(2, fit_abgmm(filter(d5, year >= 1960), est2, "fhpolrigaug",
                     c("dLdep", "dLinc"), inst_extra = "L2inc"), c("dLdep", "dLinc"))

# C3-C4: exclude former socialist countries
nonsoc <- filter(d5, sample == 1, socialist != 1)
ols_col(3, fit_ols(nonsoc, "fhpolrigaug", c("Ldep", "Linc"), TRUE), c("Ldep", "Linc"))
est4 <- complete_on(nonsoc, c("y", "dLdep", "dLinc", "L2inc"))
gmm_col(4, fit_abgmm(d5, est4, "fhpolrigaug", c("dLdep", "dLinc"), inst_extra = "L2inc"),
        c("dLdep", "dLinc"))

# C5-C6: add log population, median age, age structure
base <- filter(d5, sample == 1)
ols_col(5, fit_ols(base, "fhpolrigaug", c("Ldep", "Linc", "Llpop", "Lmed", AGE), TRUE),
        c("Ldep", "Linc", "Llpop"))
est6 <- complete_on(base, c("y", "dLdep", "dLinc", "L2inc", "dLpop", "dLmed", dAGE))
gmm_col(6, fit_abgmm(d5, est6, "fhpolrigaug", c("dLdep", "dLinc"),
                     exog = c("dLpop", "dLmed", dAGE), inst_extra = "L2inc"),
        c("dLdep", "dLinc", "dLpop"))

# C7-C8: also add education
ols_col(7, fit_ols(base, "fhpolrigaug", c("Ldep", "Ledu", "Linc", "Llpop", "Lmed", AGE), TRUE),
        c("Ldep", "Linc", "Llpop", "Ledu"))
est8 <- complete_on(base, c("y", "dLdep", "dLinc", "L2inc", "dLedu", "dLpop", "dLmed", dAGE))
gmm_col(8, fit_abgmm(d5, est8, "fhpolrigaug", c("dLdep", "dLinc"),
                     exog = c("dLedu", "dLpop", "dLmed", dAGE), inst_extra = "L2inc"),
        c("dLdep", "dLinc", "dLpop", "dLedu"))

tab <- bind_rows(rows) |> mutate(table = "4", .before = 1)
ROW_ORDER <- c("Democracy_t-1", "Log GDP per capita_t-1", "Log population_t-1",
               "Education_t-1", "Observations", "Countries", "R-squared")
write_csv(tab, file.path(PATH_OUTPUT, "table_4.csv"))
writeLines(format_table_txt(tab, "Table 4 (FH robustness)", ROW_ORDER, ncol = 8),
           file.path(PATH_OUTPUT, "table_4.txt"))
cat("Table 4 written.\n")
print(filter(tab, row %in% c("Democracy_t-1", "Log GDP per capita_t-1")) |>
        transmute(column, row = substr(row, 1, 9), value = round(value, 3), se = round(se, 3)))
