source(here::here("R", "00_setup.R"))

# table4
DEP <- "fhpolrigaug"

d5 <- read_panel(FILE_P5)
d5 <- add_lags(d5, c(DEP, "lrgdpch", "lpop", "medage", "education",
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

tb <- table_builder(DEP)
push <- tb$push
labs <- c(Ldep = LBL$dem, dLdep = LBL$dem, Linc = LBL$inc, dLinc = LBL$inc,
          Llpop = LBL$logpop, dLpop = LBL$logpop, Ledu = LBL$educ, dLedu = LBL$educ)
col_out <- function(col, m, terms, r2 = TRUE) {
  for (t in terms) push(col, labs[[t]], ce(m, t)["est"], ce(m, t)["se"])
  tb$counts(col, m, r2 = r2)
}

col_out(1, fit_ols(filter(d5, samplebalancefe == 1), DEP, c("Ldep", "Linc"), TRUE),
        c("Ldep", "Linc"))
est2 <- complete_on(filter(d5, samplebalancegmm == 1), c("y", "dLdep", "dLinc", "L2inc"))
col_out(2, fit_abgmm(filter(d5, year >= 1960), est2, DEP,
                     c("dLdep", "dLinc"), inst_extra = "L2inc"), c("dLdep", "dLinc"),
        r2 = FALSE)

nonsoc <- filter(d5, sample == 1, socialist != 1 | is.na(socialist))
col_out(3, fit_ols(nonsoc, DEP, c("Ldep", "Linc"), TRUE), c("Ldep", "Linc"))
est4 <- complete_on(nonsoc, c("y", "dLdep", "dLinc", "L2inc"))
col_out(4, fit_abgmm(d5, est4, DEP, c("dLdep", "dLinc"), inst_extra = "L2inc"),
        c("dLdep", "dLinc"), r2 = FALSE)

base <- filter(d5, sample == 1)
col_out(5, fit_ols(base, DEP, c("Ldep", "Linc", "Llpop", "Lmed", AGE), TRUE),
        c("Ldep", "Linc", "Llpop"))
est6 <- complete_on(base, c("y", "dLdep", "dLinc", "L2inc", "dLpop", "dLmed", dAGE))
col_out(6, fit_abgmm(d5, est6, DEP, c("dLdep", "dLinc"),
                     exog = c("dLpop", "dLmed", dAGE), inst_extra = "L2inc"),
        c("dLdep", "dLinc", "dLpop"), r2 = FALSE)

col_out(7, fit_ols(base, DEP, c("Ldep", "Ledu", "Linc", "Llpop", "Lmed", AGE), TRUE),
        c("Ldep", "Linc", "Llpop", "Ledu"))
est8 <- complete_on(base, c("y", "dLdep", "dLinc", "L2inc", "dLedu", "dLpop", "dLmed", dAGE))
col_out(8, fit_abgmm(d5, est8, DEP, c("dLdep", "dLinc"),
                     exog = c("dLedu", "dLpop", "dLmed", dAGE), inst_extra = "L2inc"),
        c("dLdep", "dLinc", "dLpop", "dLedu"), r2 = FALSE)

tab <- tb$collect() |> mutate(table = "4", .before = 1)
ROW_ORDER <- c(LBL$dem, LBL$inc, LBL$logpop, LBL$educ, LBL$obs, LBL$ctry, LBL$r2)
write_csv(tab, file.path(PATH_OUTPUT, "table_4.csv"))
writeLines(format_table_txt(tab, "Table 4 (FH robustness)", ROW_ORDER, ncols = 8),
           file.path(PATH_OUTPUT, "table_4.txt"))
cat("Table 4 written.\n")
print(filter(tab, row %in% c(LBL$dem, LBL$inc)) |>
        transmute(column, row = substr(row, 1, 9), value = round(value, 3), se = round(se, 3)))
