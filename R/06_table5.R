source(here::here("R", "00_setup.R"))

DEP <- "fhpolrigaug"

d5 <- prep_savings_panel()
s5 <- filter(d5, sample == 1)

tb <- table_builder(DEP)
push <- tb$push; counts <- tb$counts; sls <- tb$sls

estIV  <- complete_on(s5, c(DEP, "Linc", "z2"))
estIV3 <- complete_on(s5, c(DEP, "Ldep", "Linc", "z2"))
m1 <- fit_ols(estIV, DEP, "Linc", FALSE)
push(1, LBL$inc, ce(m1, "Linc")["est"], ce(m1, "Linc")["se"]); counts(1, m1)
m2 <- fit_ols(estIV, DEP, "Linc", TRUE)
push(2, LBL$inc, ce(m2, "Linc")["est"], ce(m2, "Linc")["se"]); counts(2, m2)
m3 <- fit_ols(estIV3, DEP, c("Ldep", "Linc"), TRUE)
push(3, LBL$dem, ce(m3, "Ldep")["est"], ce(m3, "Ldep")["se"])
push(3, LBL$inc, ce(m3, "Linc")["est"], ce(m3, "Linc")["se"]); counts(3, m3)

sls(4, estIV, inst = "z2",
    second = setNames("Linc", LBL$inc),
    fs = setNames("z2", fstg(LBL$sav2)))
sls(5, estIV3, inst = "z2", exog = "Ldep",
    second = setNames(c("Ldep", "Linc"), c(LBL$dem, LBL$inc)),
    fs = setNames(c("Ldep", "z2"), c(fstg(LBL$dem), fstg(LBL$sav2))))

est6 <- complete_on(s5, c("y", "dLdep", "dLinc", "dz"))
m6 <- fit_abgmm(d5, est6, DEP, endog = c("dLdep", "dLinc"), inst_extra = "dz")
push(6, LBL$dem, ce(m6, "dLdep")["est"], ce(m6, "dLdep")["se"])
push(6, LBL$inc, ce(m6, "dLinc")["est"], ce(m6, "dLinc")["se"]); counts(6, m6)

estL <- complete_on(s5, c(DEP, "Llabor", "Linc", "z2"))
sls(7, estL, inst = "z2", exog = "Llabor",
    second = setNames(c("Linc", "Llabor"), c(LBL$inc, LBL$labor)),
    fs = setNames(c("Llabor", "z2"), c(fstg(LBL$labor), fstg(LBL$sav2))))

est8 <- complete_on(s5, c(DEP, "Ldep", "L2dep", "L3dep", "Linc", "z2"))
m8 <- fit_iv(est8, DEP, endog = "Linc", inst = "z2",
             exog = c("Ldep", "L2dep", "L3dep"), country_fe = TRUE)
f8 <- fit_first_stage(est8, endog = "Linc", inst = c("z2", "Ldep", "L2dep", "L3dep"), country_fe = TRUE)
demlags <- c("Ldep", "L2dep", "L3dep")
push(8, LBL$dem, wald_p(m8, demlags, mod_nc(m8)), type = "ftest_p")
push(8, LBL$inc, ce(m8, "Linc")["est"], ce(m8, "Linc")["se"])
push(8, fstg(LBL$dem), wald_p(f8, demlags, mod_nc(f8)), type = "ftest_p")
push(8, fstg(LBL$sav2), ce(f8, "z2")["est"], ce(f8, "z2")["se"])
counts(8, m8, f8)

est9 <- complete_on(s5, c(DEP, "Linc", "z2", "z3"))
sls(9, est9, inst = c("z2", "z3"),
    second = setNames("Linc", LBL$inc),
    fs = setNames(c("z2", "z3"), c(fstg(LBL$sav2), fstg(LBL$sav3))))

tab <- tb$collect() |> mutate(table = "5", .before = 1)
ROW_ORDER <- c(LBL$dem, LBL$inc, LBL$labor,
               fstg(LBL$dem), fstg(LBL$labor),
               fstg(LBL$sav2), fstg(LBL$sav3),
               LBL$obs, LBL$ctry, LBL$fsr2)
write_csv(tab, file.path(PATH_OUTPUT, "table_5.csv"))
writeLines(format_table_txt(tab, "Table 5 (2SLS, savings)", ROW_ORDER, ncols = 9),
           file.path(PATH_OUTPUT, "table_5.txt"))
cat("Table 5 written.\n")
