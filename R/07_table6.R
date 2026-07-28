source(here::here("R", "00_setup.R"))

DEP <- "fhpolrigaug"

d5 <- prep_worldincome_panel()
s5 <- filter(d5, sample == 1)

tb <- table_builder(DEP)
push <- tb$push; counts <- tb$counts; sls <- tb$sls

estIV  <- complete_on(s5, c(DEP, "Linc", "z1"))
estIV3 <- complete_on(s5, c(DEP, "Ldep", "Linc", "z1"))
m1 <- fit_ols(estIV, DEP, "Linc", FALSE)
push(1, LBL$inc, ce(m1, "Linc")["est"], ce(m1, "Linc")["se"]); counts(1, m1)
m2 <- fit_ols(estIV, DEP, "Linc", TRUE)
push(2, LBL$inc, ce(m2, "Linc")["est"], ce(m2, "Linc")["se"]); counts(2, m2)
m3 <- fit_ols(estIV3, DEP, c("Ldep", "Linc"), TRUE)
push(3, LBL$dem, ce(m3, "Ldep")["est"], ce(m3, "Ldep")["se"])
push(3, LBL$inc, ce(m3, "Linc")["est"], ce(m3, "Linc")["se"]); counts(3, m3)

sls(4, estIV, inst = "z1",
    second = setNames("Linc", LBL$inc),
    fs = setNames("z1", fstg(LBL$twgdp1)))
sls(5, estIV3, inst = "z1", exog = "Ldep",
    second = setNames(c("Ldep", "Linc"), c(LBL$dem, LBL$inc)),
    fs = setNames(c("Ldep", "z1"), c(fstg(LBL$dem), fstg(LBL$twgdp1))))

est6 <- complete_on(s5, c("y", "dLdep", "dLinc", "dz_wi"))
m6 <- fit_abgmm(d5, est6, DEP, endog = c("dLdep", "dLinc"), inst_extra = "dz_wi")
push(6, LBL$dem, ce(m6, "dLdep")["est"], ce(m6, "dLdep")["se"])
push(6, LBL$inc, ce(m6, "dLinc")["est"], ce(m6, "dLinc")["se"]); counts(6, m6)

estW <- complete_on(s5, c(DEP, "wdem", "Linc", "z1"))
sls(7, estW, inst = "z1", exog = "wdem",
    second = setNames(c("Linc", "wdem"), c(LBL$inc, LBL$twdem_t)),
    fs = setNames(c("wdem", "z1"), c(fstg(LBL$twdem), fstg(LBL$twgdp1))))

estIV8 <- complete_on(s5, c(DEP, "Linc", "z2wi"))
sls(8, estIV8, inst = "z2wi",
    second = setNames("Linc", LBL$inc),
    fs = setNames("z2wi", fstg(LBL$twgdp2)))

estIV9 <- complete_on(s5, c(DEP, "Linc", "z1", "z2wi"))
sls(9, estIV9, inst = c("z1", "z2wi"),
    second = setNames("Linc", LBL$inc),
    fs = setNames(c("z1", "z2wi"), c(fstg(LBL$twgdp1), fstg(LBL$twgdp2))))

tab <- tb$collect() |> mutate(table = "6", .before = 1)
ROW_ORDER <- c(LBL$dem, LBL$inc, LBL$twdem_t,
               fstg(LBL$dem), fstg(LBL$twdem),
               fstg(LBL$twgdp1), fstg(LBL$twgdp2),
               LBL$obs, LBL$ctry, LBL$fsr2)
write_csv(tab, file.path(PATH_OUTPUT, "table_6.csv"))
writeLines(format_table_txt(tab, "Table 6 (2SLS, world income)", ROW_ORDER, ncols = 9),
           file.path(PATH_OUTPUT, "table_6.txt"))
cat("Table 6 written.\n")
