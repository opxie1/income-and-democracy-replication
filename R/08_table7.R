source(here::here("R", "00_setup.R"))

build_t7_panel <- function(file, panel) {
  DEP <- "polity4"
  d <- read_panel(file)
  d <- add_lags(d, c(DEP, "lrgdpmad", "year"), 1:2)
  d <- d |> mutate(
    Ldep = polity4_l1, Linc = lrgdpmad_l1,
    y = polity4 - polity4_l1, dLdep = polity4_l1 - polity4_l2,
    dLinc = lrgdpmad_l1 - lrgdpmad_l2, L2inc = lrgdpmad_l2)
  s <- filter(d, sample == 1)

  tb <- table_builder(DEP, panel = panel)
  push <- tb$push
  ols_col <- function(col, m, dem = TRUE) {
    if (dem) push(col, LBL$dem, ce(m, "Ldep")["est"], ce(m, "Ldep")["se"])
    push(col, LBL$inc, ce(m, "Linc")["est"], ce(m, "Linc")["se"])
    tb$counts(col, m, r2 = TRUE)
  }

  ols_col(1, fit_ols(s, DEP, c("Ldep", "Linc"), FALSE, cluster = "madid"))
  ols_col(2, fit_ols(s, DEP, c("Ldep", "Linc"), TRUE,  cluster = "madid"))

  est3 <- complete_on(s, c("y", "dLdep", "dLinc", "L2inc"))
  m3 <- fit_abgmm(d, est3, DEP, endog = c("dLdep", "dLinc"), inst_extra = "L2inc")
  push(3, LBL$dem, ce(m3, "dLdep")["est"], ce(m3, "dLdep")["se"])
  push(3, LBL$inc, ce(m3, "dLinc")["est"], ce(m3, "dLinc")["se"])
  tb$counts(3, m3)

  ols_col(4, fit_ols(s, DEP, "Linc", TRUE, cluster = "madid"), dem = FALSE)
  ols_col(5, fit_ols(filter(s, noextrapolation == 1), DEP, c("Ldep", "Linc"),
                     TRUE, cluster = "madid"))
  tb$collect()
}

tabA <- build_t7_panel(FILE_P25, "A (25-year)")
tabB <- build_t7_panel(FILE_P50, "B (50-year)")
tab <- bind_rows(tabA, tabB) |> mutate(table = "7", .before = 1)

write_csv(tab, file.path(PATH_OUTPUT, "table_7.csv"))
txt <- c(format_table_txt(filter(tab, panel == "A (25-year)"), "Table 7 Panel A (25-year)", ROWS_DYNAMIC, 5),
         "", format_table_txt(filter(tab, panel == "B (50-year)"), "Table 7 Panel B (50-year)", ROWS_DYNAMIC, 5))
writeLines(txt, file.path(PATH_OUTPUT, "table_7.txt"))
cat("Table 7 written.\n")
