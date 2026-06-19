# ---------------------------------------------------------------------------
# 07_table6.R  -- Table 6: Two-stage least squares, Freedom House democracy,
# instrumenting log GDP per capita t-1 with trade-weighted world income.
# Same layout as Table 5. Columns 8-9 use deeper lags of the instrument; column
# 7 adds trade-weighted world democracy as a regressor.
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

d5 <- read_parquet(FILE_P5) |> arrange(code, year_numeric) |> mutate(period = year_numeric)
d5 <- add_lags(d5, c("fhpolrigaug", "lrgdpch", "worldincome", "year"), 1:2)
d5 <- d5 |> mutate(
  Linc = lrgdpch_l1, Ldep = fhpolrigaug_l1, wdem = worlddemocracy,
  z1 = worldincome_l1, z2wi = worldincome_l2,
  y = fhpolrigaug - fhpolrigaug_l1, dLdep = fhpolrigaug_l1 - fhpolrigaug_l2,
  dLinc = lrgdpch_l1 - lrgdpch_l2, dz_wi = worldincome_l1 - worldincome_l2)
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
sls <- function(col, dat, inst, exog = character(), second = c(), fs = c()) {
  m <- fit_iv(dat, "fhpolrigaug", endog = "Linc", inst = inst, exog = exog, country_fe = TRUE)
  f <- fit_first_stage(dat, endog = "Linc", inst = c(inst, exog), country_fe = TRUE)
  for (lab in names(second)) push(col, lab, ce(m, second[[lab]])["est"], ce(m, second[[lab]])["se"])
  for (lab in names(fs))     push(col, lab, ce(f, fs[[lab]])["est"], ce(f, fs[[lab]])["se"])
  counts(col, m, f); invisible(NULL)
}

# C1-C3: OLS comparisons on the 2SLS sample
estIV  <- complete_on(s5, c("fhpolrigaug", "Linc", "z1"))
estIV3 <- complete_on(s5, c("fhpolrigaug", "Ldep", "Linc", "z1"))
m1 <- fit_ols(estIV, "fhpolrigaug", "Linc", FALSE)
push(1, "Log GDP per capita_t-1", ce(m1, "Linc")["est"], ce(m1, "Linc")["se"]); counts(1, m1)
m2 <- fit_ols(estIV, "fhpolrigaug", "Linc", TRUE)
push(2, "Log GDP per capita_t-1", ce(m2, "Linc")["est"], ce(m2, "Linc")["se"]); counts(2, m2)
m3 <- fit_ols(estIV3, "fhpolrigaug", c("Ldep", "Linc"), TRUE)
push(3, "Democracy_t-1", ce(m3, "Ldep")["est"], ce(m3, "Ldep")["se"])
push(3, "Log GDP per capita_t-1", ce(m3, "Linc")["est"], ce(m3, "Linc")["se"]); counts(3, m3)

# C4-C5: base 2SLS and with lagged democracy
sls(4, estIV, inst = "z1",
    second = c("Log GDP per capita_t-1" = "Linc"),
    fs = c("First stage: Trade-weighted log GDP_t-1" = "z1"))
sls(5, estIV3, inst = "z1", exog = "Ldep",
    second = c("Democracy_t-1" = "Ldep", "Log GDP per capita_t-1" = "Linc"),
    fs = c("First stage: Democracy_t-1" = "Ldep",
           "First stage: Trade-weighted log GDP_t-1" = "z1"))

# C6: Arellano-Bond GMM (income instrumented by differenced world income)
est6 <- complete_on(s5, c("y", "dLdep", "dLinc", "dz_wi"))
m6 <- fit_abgmm(d5, est6, "fhpolrigaug", endog = c("dLdep", "dLinc"), inst_extra = "dz_wi")
push(6, "Democracy_t-1", ce(m6, "dLdep")["est"], ce(m6, "dLdep")["se"])
push(6, "Log GDP per capita_t-1", ce(m6, "dLinc")["est"], ce(m6, "dLinc")["se"]); counts(6, m6)

# C7: add trade-weighted world democracy
estW <- complete_on(s5, c("fhpolrigaug", "wdem", "Linc", "z1"))
sls(7, estW, inst = "z1", exog = "wdem",
    second = c("Log GDP per capita_t-1" = "Linc", "Trade-weighted democracy_t" = "wdem"),
    fs = c("First stage: Trade-weighted democracy" = "wdem",
           "First stage: Trade-weighted log GDP_t-1" = "z1"))

# C8: instrument with world income t-2
estIV8 <- complete_on(s5, c("fhpolrigaug", "Linc", "z2wi"))
sls(8, estIV8, inst = "z2wi",
    second = c("Log GDP per capita_t-1" = "Linc"),
    fs = c("First stage: Trade-weighted log GDP_t-2" = "z2wi"))

# C9: world income t-1 and t-2 as instruments
estIV9 <- complete_on(s5, c("fhpolrigaug", "Linc", "z1", "z2wi"))
sls(9, estIV9, inst = c("z1", "z2wi"),
    second = c("Log GDP per capita_t-1" = "Linc"),
    fs = c("First stage: Trade-weighted log GDP_t-1" = "z1",
           "First stage: Trade-weighted log GDP_t-2" = "z2wi"))

tab <- bind_rows(rows) |> mutate(table = "6", .before = 1)
ROW_ORDER <- c("Democracy_t-1", "Log GDP per capita_t-1", "Trade-weighted democracy_t",
               "First stage: Democracy_t-1", "First stage: Trade-weighted democracy",
               "First stage: Trade-weighted log GDP_t-1", "First stage: Trade-weighted log GDP_t-2",
               "Observations", "Countries", "First-stage R-squared")
write_csv(tab, file.path(PATH_OUTPUT, "table_6.csv"))
writeLines(format_table_txt(tab, "Table 6 (2SLS, world income)", ROW_ORDER, ncol = 9),
           file.path(PATH_OUTPUT, "table_6.txt"))
cat("Table 6 written.\n")
