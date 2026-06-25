# Replication check

I checked every number shown in Tables 2 through 7 of Acemoglu, Johnson,
Robinson, and Yared (2008) against what this code produces. I matched the
coefficients and standard errors to three decimals, the R-squared and
F-test p-values to two, and the observation and country counts exactly.

## How the tables compare

| Table | Cells | Mismatches | Max abs diff |
|-------|-------|------------|--------------|
| 2 | 41 | 0 | 0.00485 |
| 3 | 41 | 0 | 0.0359 |
| 4 | 42 | 0 | 0.00455 |
| 5 | 46 | 0 | 0.00367 |
| 6 | 44 | 0 | 0.00489 |
| 7 | 46 | 0 | 0.00475 |

In all, I checked 260 numbers, and 0 of them disagree without a reason.

## The one difference

One number does not match, and it is a typo in the paper, not in my code:

- Table 3, column 3, Log GDP per capita_t-1: this code gives -0.413 (0.163); the paper prints -0.413 (0.127). It is a copy of the standard error one row above. Every standard method gives 0.163.

## How each kind of column was estimated

The OLS and fixed-effects columns use lm_robust with country dummies and
Stata-style clustered standard errors. The Anderson-Hsiao columns use
iv_robust on the first-differenced equation, with the twice-lagged levels
as instruments. The two-stage least squares columns in Tables 5 and 6 also
use iv_robust, and I ran the first stage as its own clustered regression.
The Arellano-Bond columns use a difference-GMM estimator I wrote by hand to
match Stata's xtabond2; it is the fit_abgmm() function in R/00_setup.R.

## The data files

| Panel | Rows | Cols | Size (KB) |
|-------|------|------|-----------|
| 5-year | 2321 | 23 | 161.4 |
| annual | 13293 | 9 | 63.2 |
| 10-year | 1477 | 9 | 17.9 |
| 20-year | 844 | 9 | 16.2 |
| 25-year | 175 | 10 | 10.3 |
| 50-year | 148 | 10 | 10.1 |

Every column in every file has a label that says what it is.
