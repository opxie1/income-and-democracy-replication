source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(pdynmc))

xcheck <- function(dep, inc, label) {
  lcol <- paste0(inc, "_l1")
  d0 <- add_lags(read_panel(FILE_P5), inc, 1) |> filter(sample == 1)
  d <- as.data.frame(d0[stats::complete.cases(d0[, c(dep, lcol)]), c("code", "year_numeric", dep, lcol)])
  names(d)[4] <- "Linc"
  mx <- n_distinct(d$year_numeric) - 2
  fit <- function(use.lev) {
    m <- pdynmc(dat = d, varname.i = "code", varname.t = "year_numeric",
      use.mc.diff = TRUE, use.mc.lev = use.lev, use.mc.nonlin = FALSE, inst.collapse = FALSE,
      include.y = TRUE, varname.y = dep, lagTerms.y = 1, maxLags.y = mx,
      include.x = TRUE, varname.reg.end = "Linc", lagTerms.reg.end = 0, maxLags.reg.end = mx,
      fur.con = FALSE, include.dum = TRUE, dum.diff = TRUE, dum.lev = use.lev,
      varname.dum = "year_numeric", w.mat = "iid.err", std.err = "corrected",
      estimation = "twostep", opt.meth = "none")
    co <- coef(m)
    c(income = unname(co[grep("Linc", names(co))[1]]), dem = unname(co[grep(dep, names(co))[1]]))
  }
  dd <- fit(FALSE); ss <- fit(TRUE)
  stopifnot(dd["income"] < 0, ss["income"] > 0, dd["dem"] > 0, ss["dem"] > 0)
  cat(sprintf("%-14s  difference GMM: income=%7.3f dem=%6.3f    system GMM: income=%7.3f dem=%6.3f\n",
              label, dd["income"], dd["dem"], ss["income"], ss["dem"]))
}

cat("pdynmc cross-check (compare signs and magnitudes with output/alternatives.txt):\n")
for (ms in MEASURES) xcheck(ms$dep, ms$inc, ms$label)
