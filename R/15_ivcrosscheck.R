source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(ivmodel))

DEP <- DEP_WEAKIV
TOL_BETA <- 1e-6
TOL_SET <- 0.005

# ivmodel
mine <- read_csv(file.path(PATH_OUTPUT, "weakiv.csv"), show_col_types = FALSE)

cat("Cross-checking the hand-coded weak-instrument code against the ivmodel package.\n")
cat(sprintf("%-16s %-22s %-22s %s\n", "Spec", "my CLR set", "ivmodel CLR set", "2SLS agree"))
for (sp in TSLS_SPECS) {
  lab <- spec_label(sp$tab, sp$col)
  row <- filter(mine, spec == lab)
  stopifnot(nrow(row) == 1L)

  dat <- spec_data(sp)
  Xd <- model.matrix(reformulate(c(FE_TERMS, sp$exog)), data = dat)[, -1]
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
cat(sprintf(paste0("\nAll %s columns agree: the hand-coded conditional likelihood ratio ",
                   "sets match\nthe ivmodel package to within %g, and the two-stage least ",
                   "squares estimates are identical.\n"),
            spell(length(TSLS_SPECS)), TOL_SET))
