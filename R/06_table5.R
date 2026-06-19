# ---------------------------------------------------------------------------
# 06_table5.R  -- Table 5: Two-stage least squares, Freedom House democracy,
# instrumenting log GDP per capita t-1 with the savings rate t-2.
# Panel A is the second stage; Panel B is the first stage. Columns 1-3 are OLS
# comparisons run on the 2SLS sample; column 6 is Arellano-Bond GMM.
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

d5 <- read_panel(FILE_P5)
d5 <- add_lags(d5, c("fhpolrigaug", "lrgdpch", "nsave", "laborshare", "year"), 1:3)
d5 <- d5 |> mutate(
  Linc = lrgdpch_l1, Ldep = fhpolrigaug_l1, L2dep = fhpolrigaug_l2, L3dep = fhpolrigaug_l3,
  z2 = nsave_l2, z3 = nsave_l3, Llabor = laborshare_l1,
  y = fhpolrigaug - fhpolrigaug_l1, dLdep = fhpolrigaug_l1 - fhpolrigaug_l2,
  dLinc = lrgdpch_l1 - lrgdpch_l2, dz = nsave_l2 - nsave_l3)
s5 <- filter(d5, sample == 1)

rows <- list()
push <- function(column, row, value, se = NA_real_, type = "coef")
  rows[[length(rows) + 1]] <<- tibble(column = column, row = row,
    value = as.numeric(value), se = as.numeric(se), type = type)
counts <- function(col, m, fs = NULL) {
  push(col, "Observations", mod_nobs(m), type = "count")
  push(col, "Countries", mod_nc(m), type = "count")
  if (!is.null(fs)) push(col, "First-stage R-squared", mod_r2(fs), type = "r2")
}

# ---- Columns 1-3: OLS comparisons, run on the 2SLS estimation sample ----
estIV  <- complete_on(s5, c("fhpolrigaug", "Linc", "z2"))
estIV3 <- complete_on(s5, c("fhpolrigaug", "Ldep", "Linc", "z2"))
m1 <- fit_ols(estIV, "fhpolrigaug", "Linc", FALSE)
push(1, "Log GDP per capita_t-1", ce(m1, "Linc")["est"], ce(m1, "Linc")["se"]); counts(1, m1)
m2 <- fit_ols(estIV, "fhpolrigaug", "Linc", TRUE)
push(2, "Log GDP per capita_t-1", ce(m2, "Linc")["est"], ce(m2, "Linc")["se"]); counts(2, m2)
m3 <- fit_ols(estIV3, "fhpolrigaug", c("Ldep", "Linc"), TRUE)
push(3, "Democracy_t-1", ce(m3, "Ldep")["est"], ce(m3, "Ldep")["se"])
push(3, "Log GDP per capita_t-1", ce(m3, "Linc")["est"], ce(m3, "Linc")["se"]); counts(3, m3)

# ---- helper for a standard 2SLS column (second stage + first stage) ----
sls <- function(col, dat, inst, exog = character(), second = c(), fs = c()) {
  m  <- fit_iv(dat, "fhpolrigaug", endog = "Linc", inst = inst, exog = exog, country_fe = TRUE)
  f  <- fit_first_stage(dat, endog = "Linc", inst = c(inst, exog), country_fe = TRUE)
  for (lab in names(second)) push(col, lab, ce(m, second[[lab]])["est"], ce(m, second[[lab]])["se"])
  for (lab in names(fs))     push(col, lab, ce(f, fs[[lab]])["est"], ce(f, fs[[lab]])["se"])
  counts(col, m, f); invisible(NULL)
}

# C4: base 2SLS
sls(4, estIV, inst = "z2",
    second = c("Log GDP per capita_t-1" = "Linc"),
    fs = c("First stage: Savings rate_t-2" = "z2"))
# C5: add lagged democracy
sls(5, estIV3, inst = "z2", exog = "Ldep",
    second = c("Democracy_t-1" = "Ldep", "Log GDP per capita_t-1" = "Linc"),
    fs = c("First stage: Democracy_t-1" = "Ldep", "First stage: Savings rate_t-2" = "z2"))

# C6: Arellano-Bond GMM (income instrumented by the differenced savings rate)
est6 <- complete_on(s5, c("y", "dLdep", "dLinc", "dz"))
m6 <- fit_abgmm(d5, est6, "fhpolrigaug", endog = c("dLdep", "dLinc"), inst_extra = "dz")
push(6, "Democracy_t-1", ce(m6, "dLdep")["est"], ce(m6, "dLdep")["se"])
push(6, "Log GDP per capita_t-1", ce(m6, "dLinc")["est"], ce(m6, "dLinc")["se"]); counts(6, m6)

# C7: add labor share
estL <- complete_on(s5, c("fhpolrigaug", "Llabor", "Linc", "z2"))
sls(7, estL, inst = "z2", exog = "Llabor",
    second = c("Log GDP per capita_t-1" = "Linc", "Labor share_t-1" = "Llabor"),
    fs = c("First stage: Labor share_t-1" = "Llabor", "First stage: Savings rate_t-2" = "z2"))

# C8: add three lags of democracy (report joint F-test p-values)
est8 <- complete_on(s5, c("fhpolrigaug", "Ldep", "L2dep", "L3dep", "Linc", "z2"))
m8 <- fit_iv(est8, "fhpolrigaug", endog = "Linc", inst = "z2",
             exog = c("Ldep", "L2dep", "L3dep"), country_fe = TRUE)
f8 <- fit_first_stage(est8, endog = "Linc", inst = c("z2", "Ldep", "L2dep", "L3dep"), country_fe = TRUE)
demlags <- c("Ldep", "L2dep", "L3dep")
push(8, "Democracy_t-1", wald_p(m8, demlags, mod_nc(m8)), type = "ftest_p")
push(8, "Log GDP per capita_t-1", ce(m8, "Linc")["est"], ce(m8, "Linc")["se"])
push(8, "First stage: Democracy_t-1", wald_p(f8, demlags, mod_nc(f8)), type = "ftest_p")
push(8, "First stage: Savings rate_t-2", ce(f8, "z2")["est"], ce(f8, "z2")["se"])
counts(8, m8, f8)

# C9: two savings instruments (t-2 and t-3)
est9 <- complete_on(s5, c("fhpolrigaug", "Linc", "z2", "z3"))
sls(9, est9, inst = c("z2", "z3"),
    second = c("Log GDP per capita_t-1" = "Linc"),
    fs = c("First stage: Savings rate_t-2" = "z2", "First stage: Savings rate_t-3" = "z3"))

tab <- bind_rows(rows) |> mutate(table = "5", .before = 1)
ROW_ORDER <- c("Democracy_t-1", "Log GDP per capita_t-1", "Labor share_t-1",
               "First stage: Democracy_t-1", "First stage: Labor share_t-1",
               "First stage: Savings rate_t-2", "First stage: Savings rate_t-3",
               "Observations", "Countries", "First-stage R-squared")
write_csv(tab, file.path(PATH_OUTPUT, "table_5.csv"))
writeLines(format_table_txt(tab, "Table 5 (2SLS, savings)", ROW_ORDER, ncol = 9),
           file.path(PATH_OUTPUT, "table_5.txt"))
cat("Table 5 written.\n")
