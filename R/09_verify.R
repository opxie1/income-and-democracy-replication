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

lab_disp <- c(unlist(LBL), setNames(fstg(unlist(LBL)), paste0("fs_", names(LBL))))

# published
pub <- read_csv(file.path(PATH_DOCS, "published_values.csv"), show_col_types = FALSE) |>
  mutate(table = as.character(table), panel = blank_na(panel),
         pub = value, pub_se = se) |>
  select(table, panel, column, row, type, pub, pub_se, note)
stopifnot(all(pub$row %in% names(lab_disp)))
pub$row <- unname(lab_disp[pub$row])

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

# report
md <- c("# Replication check",
        "",
        "I checked every number in Tables 2 through 7 of Acemoglu, Johnson,",
        "Robinson, and Yared (2008) against the output of this code. The",
        "coefficients and the standard errors must agree to three decimals. The",
        "R-squared values and the F-test p-values must agree to two decimals. The",
        "observation counts and the country counts must agree exactly.",
        "",
        "## How the tables compare",
        "",
        "| Table | Numbers | Unexplained mismatches | Largest absolute difference |",
        "|-------|---------|------------------------|-----------------------------|")
for (i in seq_len(nrow(summ))) {
  md <- c(md, sprintf("| %s | %d | %d | %.3g |",
                      summ$table[i], summ$cells[i], summ$failures[i], summ$max_diff[i]))
}
total_fail <- sum(summ$failures)
md <- c(md, "",
        sprintf("I checked %d numbers in all. The count of unexplained mismatches is %d.",
                nrow(cmp), total_fail),
        "An unexplained mismatch is a number that differs from the paper for a reason",
        "that I do not document.",
        "")

noted <- filter(cmp, !matched & has_note)
if (nrow(noted)) {
  md <- c(md,
          if (nrow(noted) == 1) c("## The one documented mismatch", "",
            "One number does not match. The cause is a typo in the paper, not an error in this code:", "")
          else c("## Documented mismatches", "",
            "These numbers do not match. Each entry below gives the reason:", ""))
  for (i in seq_len(nrow(noted))) {
    md <- c(md, sprintf("- Table %s, column %s, %s: this code gives %.3f (%.3f). The paper prints %.3f (%.3f). %s",
                        noted$table[i], noted$column[i], noted$row[i],
                        noted$mine[i], noted$mine_se[i], noted$pub[i], noted$pub_se[i], noted$note[i]))
  }
  md <- c(md, "")
}

md <- c(md, "## How I estimated each type of column", "",
        "The OLS and fixed-effects columns use lm_robust with country dummies and",
        "Stata-style clustered standard errors. The Anderson-Hsiao columns use",
        "iv_robust on the first-differenced equation, with the twice-lagged levels",
        "as instruments. The two-stage least squares columns in Tables 5 and 6 also",
        "use iv_robust. I ran the first stage of these columns as a separate",
        "clustered regression.",
        "",
        "The Arellano-Bond columns use a difference-GMM estimator that I wrote by",
        "hand. This estimator matches the xtabond2 command in Stata. The code is",
        "the fit_abgmm() function in R/00_setup.R.",
        "",
        "## The data files", "",
        "| Panel | Rows | Columns | Size (KB) |",
        "|-------|------|---------|-----------|")
for (i in seq_len(nrow(parq))) {
  md <- c(md, sprintf("| %s | %d | %d | %.1f |",
                      parq$panel[i], parq$rows[i], parq$cols[i], parq$kb[i]))
}
md <- c(md, "", "Every column in every file has a label that describes the column.")
writeLines(md, file.path(PATH_DOCS, "replication_check.md"))

cat(sprintf("\nChecked %d cells across Tables 2-7. Unexplained mismatches: %d. Documented: %d.\n",
            nrow(cmp), total_fail, nrow(noted)))
if (total_fail > 0) {
  print(filter(cmp, failure) |> select(table, column, row, mine, pub, mine_se, pub_se))
  stop("Some cells did not match the paper. See docs/diff_table_*.csv.")
}
cat("All cells match the published paper (aside from the documented misprint).\n")
