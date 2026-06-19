source(here::here("R", "00_setup.R"))

blank_na <- function(x) ifelse(is.na(x), "", as.character(x))

read_out <- function(file) {
  d <- read_csv(file, show_col_types = FALSE)
  if (!"panel" %in% names(d)) d$panel <- ""
  d |> transmute(table = as.character(table), panel = blank_na(panel),
                 column, row, type, mine = value, mine_se = se)
}

mine <- bind_rows(lapply(
  file.path(PATH_OUTPUT, paste0("table_", c(2:7), ".csv")), read_out))

pub <- read_csv(file.path(PATH_DOCS, "published_values.csv"), show_col_types = FALSE) |>
  mutate(table = as.character(table), panel = blank_na(panel),
         pub = value, pub_se = se) |>
  select(table, panel, column, row, type, pub, pub_se, note)

keys <- c("table", "panel", "column", "row", "type")
stopifnot(
  !anyDuplicated(mine[keys]),
  !anyDuplicated(pub[keys]),
  all(!is.na(pub$pub_se[pub$type == "coef"]))
)
cmp <- inner_join(mine, pub, by = keys)
stopifnot(nrow(cmp) == nrow(pub), nrow(cmp) == nrow(mine))

dec <- c(coef = 3, ftest_p = 2, count = 0, r2 = 2)
cmp <- cmp |>
  mutate(
    d_round  = dec[type],
    val_ok   = abs(round(mine * 10^d_round) / 10^d_round - pub) < 1e-6,
    se_ok    = is.na(pub_se) | (abs(round(mine_se, 3) - pub_se) < 1e-6),
    matched  = val_ok & se_ok,
    diff_val = abs(mine - pub),
    diff_se  = ifelse(is.na(pub_se), NA_real_, abs(mine_se - pub_se)),
    has_note = !is.na(note) & nzchar(note),
    failure  = !matched & !has_note
  )

summ <- cmp |>
  group_by(table) |>
  summarise(cells = n(), failures = sum(failure),
            max_diff = max(c(diff_val, diff_se[!is.na(diff_se)])), .groups = "drop")

for (tb in unique(cmp$table)) {
  out <- filter(cmp, table == tb) |>
    select(table, panel, column, row, type, mine, pub, diff_val, mine_se, pub_se, diff_se, matched, note)
  write_csv(out, file.path(PATH_DOCS, sprintf("diff_table_%s.csv", tb)))
}
write_csv(summ, file.path(PATH_DOCS, "diff_summary.csv"))
print(summ)

panel_files <- c("5-year" = FILE_P5, "annual" = FILE_PA, "10-year" = FILE_P10,
                 "20-year" = FILE_P20, "25-year" = FILE_P25, "50-year" = FILE_P50)
parq <- tibble()
for (nm in names(panel_files)) {
  d <- read_parquet(panel_files[[nm]])
  miss <- names(d)[vapply(d, function(x) {
    lb <- attr(x, "label"); is.null(lb) || !nzchar(lb) }, logical(1))]
  if (length(miss)) stop("Missing label(s) in ", nm, ": ", paste(miss, collapse = ", "))
  parq <- bind_rows(parq, tibble(panel = nm, rows = nrow(d), cols = ncol(d),
                                 kb = round(file.info(panel_files[[nm]])$size / 1024, 1)))
}

md <- c("# Replication check",
        "",
        "Every displayed cell of Tables 2-7 in Acemoglu, Johnson, Robinson, and",
        "Yared (2008) is compared against this code's output. Coefficients and",
        "standard errors are checked to three decimals, R-squared and F-test",
        "p-values to two, observation and country counts exactly.",
        "",
        "## Tables vs the published paper",
        "",
        "| Table | Cells | Mismatches | Max abs diff |",
        "|-------|-------|------------|--------------|")
for (i in seq_len(nrow(summ))) {
  md <- c(md, sprintf("| %s | %d | %d | %.3g |",
                      summ$table[i], summ$cells[i], summ$failures[i], summ$max_diff[i]))
}
total_fail <- sum(summ$failures)
md <- c(md, "",
        sprintf("Total cells checked: %d. Unexplained mismatches: %d.",
                nrow(cmp), total_fail),
        "")

noted <- filter(cmp, !matched & has_note)
if (nrow(noted)) {
  md <- c(md, "## Documented discrepancies", "")
  for (i in seq_len(nrow(noted))) {
    md <- c(md, sprintf("- Table %s, column %s, %s: this code gives %.3f (%.3f); the paper prints %.3f (%.3f). %s",
                        noted$table[i], noted$column[i], noted$row[i],
                        noted$mine[i], noted$mine_se[i], noted$pub[i], noted$pub_se[i], noted$note[i]))
  }
  md <- c(md, "")
}

md <- c(md, "## Method notes", "",
        "- OLS and fixed-effects columns: `lm_robust` with country dummies and",
        "  Stata-style clustered standard errors.",
        "- Anderson-Hsiao columns: `iv_robust` on the first-differenced equation,",
        "  instrumenting the lagged differences with twice-lagged levels.",
        "- Two-stage least squares columns (Tables 5-6): `iv_robust`, with the",
        "  first stage reported as a separate clustered regression.",
        "- Arellano-Bond columns: a one-step difference-GMM estimator written to",
        "  match Stata's `xtabond2 ... noleveleq robust` (uncollapsed GMM-style",
        "  instruments, level passthru instrument, H = first-difference weight,",
        "  one-step robust SEs). See fit_abgmm() in R/00_setup.R.",
        "",
        "## Parquet datasets", "",
        "| Panel | Rows | Cols | Size (KB) |",
        "|-------|------|------|-----------|")
for (i in seq_len(nrow(parq))) {
  md <- c(md, sprintf("| %s | %d | %d | %.1f |",
                      parq$panel[i], parq$rows[i], parq$cols[i], parq$kb[i]))
}
md <- c(md, "", "Every column in every panel carries a variable label.")
writeLines(md, file.path(PATH_DOCS, "replication_check.md"))

cat(sprintf("\nChecked %d cells across Tables 2-7. Unexplained mismatches: %d. Documented: %d.\n",
            nrow(cmp), total_fail, nrow(noted)))
if (total_fail > 0) {
  print(filter(cmp, failure) |> select(table, column, row, mine, pub, mine_se, pub_se))
  stop("Some cells did not match the paper; see docs/diff_table_*.csv.")
}
cat("All cells match the published paper (aside from the documented misprint).\n")
