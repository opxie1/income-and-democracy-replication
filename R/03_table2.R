source(here::here("R", "00_setup.R"))

# table2
tab <- build_dynamic_table("fhpolrigaug")
tab <- mutate(tab, table = "2", .before = 1)

write_csv(tab, file.path(PATH_OUTPUT, "table_2.csv"))
writeLines(format_table_txt(tab, "Table 2 (Freedom House)", ROWS_DYNAMIC),
           file.path(PATH_OUTPUT, "table_2.txt"))

cat("Table 2 written. Income coefficients by column:\n")
inc <- filter(tab, row == LBL$inc)
print(transmute(inc, column, type, value = round(value, 3), se = round(se, 3)))
