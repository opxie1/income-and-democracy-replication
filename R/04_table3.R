source(here::here("R", "00_setup.R"))

# table3
tab <- build_dynamic_table("polity4")
tab <- mutate(tab, table = "3", .before = 1)

write_csv(tab, file.path(PATH_OUTPUT, "table_3.csv"))
writeLines(format_table_txt(tab, "Table 3 (Polity)", ROWS_DYNAMIC),
           file.path(PATH_OUTPUT, "table_3.txt"))

cat("Table 3 written. Income coefficients by column:\n")
inc <- filter(tab, row == LBL$inc)
print(transmute(inc, column, type, value = round(value, 3), se = round(se, 3)))
