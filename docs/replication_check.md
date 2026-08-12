# Replication check

I checked every number in Tables 2 through 7 of Acemoglu, Johnson,
Robinson, and Yared (2008) against the output of this code. The
coefficients and the standard errors must agree to three decimals. The
R-squared values and the F-test p-values must agree to two decimals. The
observation counts and the country counts must agree exactly.

## How the tables compare

| Table | Numbers | Unexplained mismatches | Largest absolute difference |
|-------|---------|------------------------|-----------------------------|
| 2 | 41 | 0 | 0.00485 |
| 3 | 41 | 0 | 0.0359 |
| 4 | 42 | 0 | 0.00455 |
| 5 | 46 | 0 | 0.00367 |
| 6 | 44 | 0 | 0.00489 |
| 7 | 46 | 0 | 0.00475 |

I checked 260 numbers in all. The count of unexplained mismatches is 0.
An unexplained mismatch is a number that differs from the paper for a reason
that I do not document.

## The one documented mismatch

One number does not match. The cause is a typo in the paper, not an error in this code:

- Table 3, column 3, Log GDP per capita_t-1: this code gives -0.413 (0.163). The paper prints -0.413 (0.127). It is a copy of the standard error one row above. Every standard method gives 0.163.

## How I estimated each type of column

The OLS and fixed-effects columns use lm_robust with country dummies and
Stata-style clustered standard errors. The Anderson-Hsiao columns use
iv_robust on the first-differenced equation, with the twice-lagged levels
as instruments. The two-stage least squares columns in Tables 5 and 6 also
use iv_robust. I ran the first stage of these columns as a separate
clustered regression.

The Arellano-Bond columns use a difference-GMM estimator that I wrote by
hand. This estimator matches the xtabond2 command in Stata. The code is
the fit_abgmm() function in R/00_setup.R.

## The data files

| Panel | Rows | Columns | Size (KB) |
|-------|------|---------|-----------|
| 5-year | 2321 | 23 | 161.4 |
| annual | 13293 | 9 | 63.2 |
| 10-year | 1477 | 9 | 17.9 |
| 20-year | 844 | 9 | 16.2 |
| 25-year | 175 | 10 | 10.3 |
| 50-year | 148 | 10 | 10.1 |

Every column in every file has a label that describes the column.
