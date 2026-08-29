# dynamic
derive_dynamic <- function(d, dep, inc = "lrgdpch", lags = 1:2) {
  d <- add_lags(d, c(dep, inc, "year"), lags)
  dl1 <- paste0(dep, "_l1"); dl2 <- paste0(dep, "_l2")
  il1 <- paste0(inc, "_l1"); il2 <- paste0(inc, "_l2")
  d$Ldep  <- d[[dl1]]
  d$Linc  <- d[[il1]]
  d$y     <- d[[dep]] - d[[dl1]]
  d$dLdep <- d[[dl1]] - d[[dl2]]
  d$dLinc <- d[[il1]] - d[[il2]]
  d$L2dep <- d[[dl2]]
  d$L2inc <- d[[il2]]
  d
}

prep_dynamic <- function(file, dep, inc = "lrgdpch", lags = 1:2)
  derive_dynamic(read_panel(file), dep, inc, lags)

# savings
prep_savings_panel <- function() {
  d <- add_lags(read_panel(FILE_P5),
                c("fhpolrigaug", "lrgdpch", "nsave", "laborshare", "year"), 1:3)
  mutate(d,
    Linc = lrgdpch_l1, Ldep = fhpolrigaug_l1, L2dep = fhpolrigaug_l2, L3dep = fhpolrigaug_l3,
    z2 = nsave_l2, z3 = nsave_l3, Llabor = laborshare_l1,
    y = fhpolrigaug - fhpolrigaug_l1, dLdep = fhpolrigaug_l1 - fhpolrigaug_l2,
    dLinc = lrgdpch_l1 - lrgdpch_l2, dz = nsave_l2 - nsave_l3)
}

# worldincome
prep_worldincome_panel <- function() {
  d <- add_lags(read_panel(FILE_P5),
                c("fhpolrigaug", "lrgdpch", "worldincome", "year"), 1:2)
  mutate(d,
    Linc = lrgdpch_l1, Ldep = fhpolrigaug_l1, wdem = worlddemocracy,
    z1 = worldincome_l1, z2wi = worldincome_l2,
    y = fhpolrigaug - fhpolrigaug_l1, dLdep = fhpolrigaug_l1 - fhpolrigaug_l2,
    dLinc = lrgdpch_l1 - lrgdpch_l2, dz_wi = worldincome_l1 - worldincome_l2)
}

# tables
build_dynamic_table <- function(dep) {
  tb <- table_builder(dep)
  push <- tb$push
  coefcol <- function(col, m, dem = NULL, inc = NULL) {
    if (!is.null(dem)) push(col, LBL$dem, ce(m, dem)["est"], ce(m, dem)["se"])
    if (!is.null(inc)) push(col, LBL$inc, ce(m, inc)["est"], ce(m, inc)["se"])
  }
  stats_col <- function(col, m, r2 = TRUE) tb$counts(col, m, r2 = r2)

  d5 <- prep_dynamic(FILE_P5, dep); s5 <- filter(d5, sample == 1)

  m1 <- fit_ols(s5, dep, c("Ldep", "Linc"), country_fe = FALSE); coefcol(1, m1, "Ldep", "Linc"); stats_col(1, m1)
  m2 <- fit_ols(s5, dep, c("Ldep", "Linc"), country_fe = TRUE);  coefcol(2, m2, "Ldep", "Linc"); stats_col(2, m2)
  m3 <- fit_iv(s5, "y", endog = c("dLdep", "dLinc"), inst = c("L2dep", "L2inc"),
               country_fe = FALSE);                              coefcol(3, m3, "dLdep", "dLinc"); stats_col(3, m3, r2 = FALSE)
  est4 <- complete_on(s5, c("y", "dLdep", "dLinc", "L2inc"))
  m4 <- fit_abgmm(d5, est4, dep_level = dep, endog = c("dLdep", "dLinc"), inst_extra = "L2inc")
  coefcol(4, m4, "dLdep", "dLinc"); stats_col(4, m4, r2 = FALSE)
  m5 <- fit_ols(s5, dep, c("Linc"), country_fe = TRUE);          coefcol(5, m5, inc = "Linc"); stats_col(5, m5)

  da <- add_lags(read_panel(FILE_PA), c(dep, "lrgdpch"), 1:5)
  deplags <- paste0(dep, "_l", 1:5); inclags <- paste0("lrgdpch_l", 1:5)
  sa <- filter(da, sample == 1)
  m6 <- fit_ols(sa, dep, c(deplags, inclags), country_fe = TRUE)
  push(6, LBL$dem, wald_p(m6, deplags, mod_nc(m6)), type = "ftest_p")
  push(6, LBL$inc, wald_p(m6, inclags, mod_nc(m6)), type = "ftest_p")
  stats_col(6, m6)

  d10 <- prep_dynamic(FILE_P10, dep); s10 <- filter(d10, sample == 1)
  m7 <- fit_ols(s10, dep, c("Ldep", "Linc"), country_fe = TRUE); coefcol(7, m7, "Ldep", "Linc"); stats_col(7, m7)
  est8 <- complete_on(s10, c("y", "dLdep", "dLinc", "L2inc"))
  m8 <- fit_abgmm(d10, est8, dep_level = dep, endog = c("dLdep", "dLinc"), inst_extra = "L2inc")
  coefcol(8, m8, "dLdep", "dLinc"); stats_col(8, m8, r2 = FALSE)

  d20 <- prep_dynamic(FILE_P20, dep); s20 <- filter(d20, sample == 1)
  m9 <- fit_ols(s20, dep, c("Ldep", "Linc"), country_fe = TRUE); coefcol(9, m9, "Ldep", "Linc"); stats_col(9, m9)

  tb$collect()
}

# layout
format_table_txt <- function(tab, title, row_order, ncols = max(tab$column)) {
  cell <- function(col, rw) {
    r <- tab[tab$column == col & tab$row == rw, ]
    if (nrow(r) == 0) return("")
    switch(r$type[1],
      coef    = sprintf("%.3f", r$value),
      ftest_p = sprintf("[%.2f]", r$value),
      count   = formatC(r$value, format = "d"),
      r2      = sprintf("%.2f", r$value),
      sprintf("%.3f", r$value))
  }
  se_cell <- function(col, rw) {
    r <- tab[tab$column == col & tab$row == rw, ]
    if (nrow(r) == 0 || r$type[1] != "coef" || is.na(r$se[1])) return("")
    sprintf("(%.3f)", r$se[1])
  }
  lw <- max(26L, nchar(c(title, row_order)) + 2L)
  hdr <- c(sprintf("%-*s", lw, title), sprintf("%9s", paste0("(", 1:ncols, ")")))
  out <- paste(hdr, collapse = "")
  for (rw in row_order) {
    line  <- c(sprintf("%-*s", lw, rw),
               vapply(1:ncols, function(cc) sprintf("%9s", cell(cc, rw)), character(1)))
    out <- c(out, paste(line, collapse = ""))
    ses <- vapply(1:ncols, function(cc) sprintf("%9s", se_cell(cc, rw)), character(1))
    if (any(nzchar(trimws(ses))))
      out <- c(out, paste(c(sprintf("%-*s", lw, ""), ses), collapse = ""))
  }
  out
}
