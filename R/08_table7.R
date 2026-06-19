# ---------------------------------------------------------------------------
# 08_table7.R  -- Table 7: Polity democracy in the long run (Maddison data).
# Panel A is the 25-year sample (1875-2000); Panel B the 50-year sample
# (1900-2000). Columns: pooled OLS, fixed effects, Arellano-Bond GMM, fixed
# effects without lagged democracy, and fixed effects on non-extrapolated data.
# OLS/IV SEs cluster by Maddison aggregation id; the GMM column clusters by
# country, matching the authors' Stata.
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

build_t7_panel <- function(file, panel) {
  d <- read_panel(file)
  d <- add_lags(d, c("polity4", "lrgdpmad", "year"), 1:2)
  d <- d |> mutate(
    Ldep = polity4_l1, Linc = lrgdpmad_l1,
    y = polity4 - polity4_l1, dLdep = polity4_l1 - polity4_l2,
    dLinc = lrgdpmad_l1 - lrgdpmad_l2, L2inc = lrgdpmad_l2)
  s <- filter(d, sample == 1)

  rows <- list()
  push <- function(column, row, value, se = NA_real_, type = "coef")
    rows[[length(rows) + 1]] <<- tibble(panel = panel, column = column, row = row,
      value = as.numeric(value), se = as.numeric(se), type = type)
  ols_col <- function(col, m, dem = TRUE) {
    if (dem) push(col, "Democracy_t-1", ce(m, "Ldep")["est"], ce(m, "Ldep")["se"])
    push(col, "Log GDP per capita_t-1", ce(m, "Linc")["est"], ce(m, "Linc")["se"])
    push(col, "Observations", mod_nobs(m), type = "count")
    push(col, "Countries", mod_nc(m), type = "count")
    push(col, "R-squared", mod_r2(m), type = "r2")
  }

  ols_col(1, fit_ols(s, "polity4", c("Ldep", "Linc"), FALSE, cluster = "madid"))
  ols_col(2, fit_ols(s, "polity4", c("Ldep", "Linc"), TRUE,  cluster = "madid"))

  est3 <- complete_on(s, c("y", "dLdep", "dLinc", "L2inc"))
  m3 <- fit_abgmm(d, est3, "polity4", endog = c("dLdep", "dLinc"), inst_extra = "L2inc")
  push(3, "Democracy_t-1", ce(m3, "dLdep")["est"], ce(m3, "dLdep")["se"])
  push(3, "Log GDP per capita_t-1", ce(m3, "dLinc")["est"], ce(m3, "dLinc")["se"])
  push(3, "Observations", mod_nobs(m3), type = "count")
  push(3, "Countries", mod_nc(m3), type = "count")

  ols_col(4, fit_ols(s, "polity4", "Linc", TRUE, cluster = "madid"), dem = FALSE)
  ols_col(5, fit_ols(filter(s, noextrapolation == 1), "polity4", c("Ldep", "Linc"),
                     TRUE, cluster = "madid"))
  bind_rows(rows)
}

tabA <- build_t7_panel(FILE_P25, "A (25-year)")
tabB <- build_t7_panel(FILE_P50, "B (50-year)")
tab <- bind_rows(tabA, tabB) |> mutate(table = "7", .before = 1)

write_csv(tab, file.path(PATH_OUTPUT, "table_7.csv"))
ROW_ORDER <- c("Democracy_t-1", "Log GDP per capita_t-1", "Observations", "Countries", "R-squared")
txt <- c(format_table_txt(filter(tab, panel == "A (25-year)"), "Table 7 Panel A (25-year)", ROW_ORDER, 5),
         "", format_table_txt(filter(tab, panel == "B (50-year)"), "Table 7 Panel B (50-year)", ROW_ORDER, 5))
writeLines(txt, file.path(PATH_OUTPUT, "table_7.txt"))
cat("Table 7 written.\n")
