# Replication check

Every displayed cell of Tables 2-7 in Acemoglu, Johnson, Robinson, and
Yared (2008) is compared against this code's output. Coefficients and
standard errors are checked to three decimals, R-squared and F-test
p-values to two, observation and country counts exactly.

## Tables vs the published paper

| Table | Cells | Mismatches | Max abs diff |
|-------|-------|------------|--------------|
| 2 | 41 | 0 | 0.00485 |
| 3 | 41 | 0 | 0.0359 |
| 4 | 42 | 0 | 0.00455 |
| 5 | 46 | 0 | 0.00367 |
| 6 | 44 | 0 | 0.00489 |
| 7 | 46 | 0 | 0.00475 |

Total cells checked: 260. Unexplained mismatches: 0.

## Documented discrepancies

- Table 3, column 3, Log GDP per capita_t-1: this code gives -0.413 (0.163); the paper prints -0.413 (0.127). paper misprints SE as 0.127 (duplicate of the democracy SE above); all standard methods give 0.163

## Method notes

- OLS and fixed-effects columns: `lm_robust` with country dummies and
  Stata-style clustered standard errors.
- Anderson-Hsiao columns: `iv_robust` on the first-differenced equation,
  instrumenting the lagged differences with twice-lagged levels.
- Two-stage least squares columns (Tables 5-6): `iv_robust`, with the
  first stage reported as a separate clustered regression.
- Arellano-Bond columns: a one-step difference-GMM estimator written to
  match Stata's `xtabond2 ... noleveleq robust` (uncollapsed GMM-style
  instruments, level passthru instrument, H = first-difference weight,
  one-step robust SEs). See fit_abgmm() in R/00_setup.R.

## Parquet datasets

| Panel | Rows | Cols | Size (KB) |
|-------|------|------|-----------|
| 5-year | 2321 | 23 | 161.4 |
| annual | 13293 | 9 | 63.2 |
| 10-year | 1477 | 9 | 17.9 |
| 20-year | 844 | 9 | 16.2 |
| 25-year | 175 | 10 | 10.3 |
| 50-year | 148 | 10 | 10.1 |

Every column in every panel carries a variable label.
