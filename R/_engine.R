# ---------------------------------------------------------------------------
# _engine.R  (sourced by 00_setup.R)
# Shared builders that assemble whole published tables from the fit_* helpers.
# Each builder returns a tidy tibble: column, row, value, se, type, where
#   type = "coef"  -> value is a coefficient, se its clustered SE
#          "ftest_p"-> value is the p-value of a joint F-test (annual columns)
#          "count" -> value is Observations or Countries
#          "r2"    -> value is an R-squared (or first-stage R-squared)
# ---------------------------------------------------------------------------

# Convenience: build the differenced/lagged columns a dynamic panel needs.
.prep_dynamic <- function(file, dep, lags = 1:2) {
  d <- read_parquet(file) |> arrange(code, year_numeric) |> mutate(period = year_numeric)
  d <- add_lags(d, c(dep, "lrgdpch", "year"), lags)
  dl1 <- paste0(dep, "_l1"); dl2 <- paste0(dep, "_l2")
  d$Ldep  <- d[[dl1]]
  d$Linc  <- d$lrgdpch_l1
  d$y     <- d[[dep]] - d[[dl1]]
  d$dLdep <- d[[dl1]] - d[[dl2]]
  d$dLinc <- d$lrgdpch_l1 - d$lrgdpch_l2
  d$L2dep <- d[[dl2]]
  d$L2inc <- d$lrgdpch_l2
  d
}

# Tables 2 (Freedom House) and 3 (Polity): nine columns spanning five-year,
# annual, ten-year, and twenty-year data.
build_dynamic_table <- function(dep) {
  rows <- list()
  push <- function(column, row, value, se = NA_real_, type = "coef") {
    rows[[length(rows) + 1]] <<- tibble(column = column, row = row,
      value = as.numeric(value), se = as.numeric(se), type = type)
  }
  coefcol <- function(col, m, dem = NULL, inc = NULL) {
    if (!is.null(dem)) push(col, "Democracy_t-1",         ce(m, dem)["est"], ce(m, dem)["se"])
    if (!is.null(inc)) push(col, "Log GDP per capita_t-1", ce(m, inc)["est"], ce(m, inc)["se"])
  }
  stats <- function(col, m, r2 = TRUE) {
    push(col, "Observations", mod_nobs(m), type = "count")
    push(col, "Countries",    mod_nc(m),   type = "count")
    if (r2) push(col, "R-squared", mod_r2(m), type = "r2")
  }

  # ---- five-year data: columns 1-5 ----
  d5 <- .prep_dynamic(FILE_P5, dep); s5 <- filter(d5, sample == 1)

  m1 <- fit_ols(s5, dep, c("Ldep", "Linc"), country_fe = FALSE); coefcol(1, m1, "Ldep", "Linc"); stats(1, m1)
  m2 <- fit_ols(s5, dep, c("Ldep", "Linc"), country_fe = TRUE);  coefcol(2, m2, "Ldep", "Linc"); stats(2, m2)
  m3 <- fit_iv(s5, "y", endog = c("dLdep", "dLinc"), inst = c("L2dep", "L2inc"),
               country_fe = FALSE);                              coefcol(3, m3, "dLdep", "dLinc"); stats(3, m3, r2 = FALSE)
  est4 <- complete_on(s5, c("y", "dLdep", "dLinc", "L2inc"))
  m4 <- fit_abgmm(d5, est4, dep_level = dep, endog = c("dLdep", "dLinc"), inst_extra = "L2inc")
  coefcol(4, m4, "dLdep", "dLinc"); stats(4, m4, r2 = FALSE)
  m5 <- fit_ols(s5, dep, c("Linc"), country_fe = TRUE);          coefcol(5, m5, inc = "Linc"); stats(5, m5)

  # ---- annual data: column 6 (five lags each; report F-test p-values) ----
  da <- read_parquet(FILE_PA) |> arrange(code, year_numeric) |> mutate(period = year_numeric)
  da <- add_lags(da, c(dep, "lrgdpch"), 1:5)
  deplags <- paste0(dep, "_l", 1:5); inclags <- paste0("lrgdpch_l", 1:5)
  sa <- filter(da, sample == 1)
  m6 <- fit_ols(sa, dep, c(deplags, inclags), country_fe = TRUE)
  push(6, "Democracy_t-1",          wald_p(m6, deplags, mod_nc(m6)), type = "ftest_p")
  push(6, "Log GDP per capita_t-1", wald_p(m6, inclags, mod_nc(m6)), type = "ftest_p")
  stats(6, m6)

  # ---- ten-year data: columns 7-8 ----
  d10 <- .prep_dynamic(FILE_P10, dep); s10 <- filter(d10, sample == 1)
  m7 <- fit_ols(s10, dep, c("Ldep", "Linc"), country_fe = TRUE); coefcol(7, m7, "Ldep", "Linc"); stats(7, m7)
  est8 <- complete_on(s10, c("y", "dLdep", "dLinc", "L2inc"))
  m8 <- fit_abgmm(d10, est8, dep_level = dep, endog = c("dLdep", "dLinc"), inst_extra = "L2inc")
  coefcol(8, m8, "dLdep", "dLinc"); stats(8, m8, r2 = FALSE)

  # ---- twenty-year data: column 9 ----
  d20 <- .prep_dynamic(FILE_P20, dep); s20 <- filter(d20, sample == 1)
  m9 <- fit_ols(s20, dep, c("Ldep", "Linc"), country_fe = TRUE); coefcol(9, m9, "Ldep", "Linc"); stats(9, m9)

  bind_rows(rows)
}

# Render a tidy table tibble as a fixed-width text block laid out like the paper.
format_table_txt <- function(tab, title, row_order, ncol = max(tab$column)) {
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
  hdr <- c(sprintf("%-26s", title), sprintf("%9s", paste0("(", 1:ncol, ")")))
  out <- paste(hdr, collapse = "")
  for (rw in row_order) {
    line  <- c(sprintf("%-26s", rw), sapply(1:ncol, function(cc) sprintf("%9s", cell(cc, rw))))
    out <- c(out, paste(line, collapse = ""))
    ses <- sapply(1:ncol, function(cc) se_cell(cc, rw))
    if (any(nzchar(ses)))
      out <- c(out, paste(c(sprintf("%-26s", ""), sprintf("%9s", ses)), collapse = ""))
  }
  out
}
