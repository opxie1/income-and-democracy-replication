# ---------------------------------------------------------------------------
# 03_table2.R  -- Table 2: Fixed-effects results, Freedom House democracy.
# Pooled OLS, fixed effects, Anderson-Hsiao IV, Arellano-Bond GMM, and the
# annual / ten-year / twenty-year columns. (Columns 3, 4, 8 are the IV/GMM
# estimates; the rest are OLS comparisons.)
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

ROW_ORDER <- c("Democracy_t-1", "Log GDP per capita_t-1",
               "Observations", "Countries", "R-squared")

tab <- build_dynamic_table("fhpolrigaug")
tab <- mutate(tab, table = "2", .before = 1)

write_csv(tab, file.path(PATH_OUTPUT, "table_2.csv"))
writeLines(format_table_txt(tab, "Table 2 (Freedom House)", ROW_ORDER),
           file.path(PATH_OUTPUT, "table_2.txt"))

cat("Table 2 written. Income coefficients by column:\n")
inc <- filter(tab, row == "Log GDP per capita_t-1")
print(transmute(inc, column, type, value = round(value, 3), se = round(se, 3)))
