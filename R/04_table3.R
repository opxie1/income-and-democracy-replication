# ---------------------------------------------------------------------------
# 04_table3.R  -- Table 3: Fixed-effects results, Polity democracy.
# Same nine-column structure as Table 2, with the Polity index as the
# dependent variable.
# ---------------------------------------------------------------------------

source(here::here("R", "00_setup.R"))

ROW_ORDER <- c("Democracy_t-1", "Log GDP per capita_t-1",
               "Observations", "Countries", "R-squared")

tab <- build_dynamic_table("polity4")
tab <- mutate(tab, table = "3", .before = 1)

write_csv(tab, file.path(PATH_OUTPUT, "table_3.csv"))
writeLines(format_table_txt(tab, "Table 3 (Polity)", ROW_ORDER),
           file.path(PATH_OUTPUT, "table_3.txt"))

cat("Table 3 written. Income coefficients by column:\n")
inc <- filter(tab, row == "Log GDP per capita_t-1")
print(transmute(inc, column, type, value = round(value, 3), se = round(se, 3)))
