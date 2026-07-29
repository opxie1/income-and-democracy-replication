source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(ivmodel))

DEP <- "fhpolrigaug"
TOL_BETA <- 1e-6
TOL_SET <- 0.02

s5 <- filter(prep_savings_panel(), sample == 1)
s6 <- filter(prep_worldincome_panel(), sample == 1)

SPECS <- list(
  list(tab = "5", col = 4, d = s5, inst = "z2",            exog = character()),
  list(tab = "5", col = 5, d = s5, inst = "z2",            exog = "Ldep"),
  list(tab = "5", col = 7, d = s5, inst = "z2",            exog = "Llabor"),
  list(tab = "5", col = 8, d = s5, inst = "z2",            exog = c("Ldep", "L2dep", "L3dep")),
  list(tab = "5", col = 9, d = s5, inst = c("z2", "z3"),   exog = character()),
  list(tab = "6", col = 4, d = s6, inst = "z1",            exog = character()),
  list(tab = "6", col = 5, d = s6, inst = "z1",            exog = "Ldep"),
  list(tab = "6", col = 7, d = s6, inst = "z1",            exog = "wdem"),
  list(tab = "6", col = 8, d = s6, inst = "z2wi",          exog = character()),
  list(tab = "6", col = 9, d = s6, inst = c("z1", "z2wi"), exog = character()))

mine <- read_csv(file.path(PATH_OUTPUT, "weakiv.csv"), show_col_types = FALSE)

cat("Cross-checking the hand-coded weak-instrument code against the ivmodel package.\n")
cat(sprintf("%-16s %-22s %-22s %s\n", "Spec", "my CLR set", "ivmodel CLR set", "2SLS agree"))
for (sp in SPECS) {
  lab <- sprintf("Table %s col %d", sp$tab, sp$col)
  row <- filter(mine, spec == lab)
  stopifnot(nrow(row) == 1L)

  dat <- complete_on(sp$d, c(DEP, "Linc", sp$inst, sp$exog, "code"))
  Xd <- model.matrix(reformulate(c("factor(year)", "factor(code)", sp$exog)), data = dat)[, -1]
  Zm <- as.matrix(dat[, sp$inst, drop = FALSE])
  dimnames(Zm) <- list(NULL, sp$inst); storage.mode(Zm) <- "double"
  iv <- ivmodel(Y = as.numeric(dat[[DEP]]), D = as.numeric(dat$Linc),
                Z = Zm, X = Xd, intercept = TRUE)
  ref_b <- unname(coef(iv)["TSLS", "Estimate"])
  ref_clr <- as.numeric(CLR(iv)$ci)

  stopifnot(abs(ref_b - row$tsls) < TOL_BETA,
            length(ref_clr) == 2L,
            abs(ref_clr[1] - row$clr_lo) < TOL_SET,
            abs(ref_clr[2] - row$clr_hi) < TOL_SET)

  if (length(sp$inst) == 1L) {
    ref_ar <- as.numeric(AR.test(iv)$ci)
    stopifnot(max(abs(ref_ar - ref_clr)) < TOL_SET)
  }

  cat(sprintf("%-16s %-22s %-22s %s\n", lab, row$clr_set,
              sprintf("[%.2f, %.2f]", ref_clr[1], ref_clr[2]), "yes"))
}
cat("\nAll ten columns agree: the hand-coded conditional likelihood ratio sets match\n",
    "the ivmodel package, and the two-stage least squares estimates are identical.\n", sep = "")
