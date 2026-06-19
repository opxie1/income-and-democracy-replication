# Income and Democracy (Acemoglu, Johnson, Robinson, Yared 2008): IV replication

This project reproduces the regression tables that report instrumental-variables
estimates in:

> Acemoglu, Daron, Simon Johnson, James A. Robinson, and Pierre Yared. 2008. "Income and Democracy." *American Economic Review* 98 (3): 808-842.

That means **Tables 2, 3, 4, 5, 6, and 7** in full. These are every table that
contains an IV or GMM estimate, with the OLS and fixed-effects columns kept
alongside them for comparison. Table 1 (cross-sectional correlations) and Table 8
(the 1500-2000 cross-section) have no IV estimates and are not included.

Each table is reproduced column by column, including:

- pooled OLS and fixed-effects OLS, with country-clustered standard errors;
- the **Anderson-Hsiao** instrumental-variables estimator (Tables 2 and 3, column 3);
- the **Arellano-Bond** one-step difference GMM estimator (the GMM columns of Tables 2-7);
- **two-stage least squares**, instrumenting income with the lagged savings rate
  (Table 5) and with trade-weighted world income (Table 6), with first stages.

## Setup

1. Install R (version 4 or newer). On Windows: `winget install RProject.R`.
2. Download the authors' data from [openICPSR project 113251](https://www.openicpsr.org/openicpsr/project/113251) (the file `Income-and-Democracy-Data-AER-adjustment.xls` and its readme) and unzip it into a folder called `replication-kit` in the project root. The workbook holds one sheet per sample; the Stata commands that generate each table live in a text box on each sheet.
3. Install the R packages:

```r
install.packages(c("tidyverse", "arrow", "readxl", "labelled", "estimatr", "plm", "here"))
```

`MASS` is also used (it ships with R). The pipeline does not need `fixest` or `pdftools`.

## How to run it

```sh
Rscript R/run_all.R
```

The scripts read the workbook, build one compressed parquet file per sample, make
Tables 2 through 7, and then check every number against the published paper. A
finished run prints a line reporting how many cells matched. On Windows, R can
still show a non-zero exit code after that line; it comes from a known shutdown
bug in the arrow package, not from these scripts, and every file is written by
then.

## How the code is organized

- `R/00_setup.R`: paths, variable labels, and the estimation helpers
  (`fit_ols`, `fit_iv`, `fit_abgmm`, lag construction, clustered Wald tests).
- `R/01_load_data.R`: read the six sample sheets from the workbook.
- `R/02_build_panels.R`: trim to the variables the tables use, label every
  column, and write one parquet per sample.
- `R/03_table2.R` ... `R/08_table7.R`: one script per table.
- `R/09_verify.R`: compare every cell to the paper and write the check report.
- `R/_engine.R`: shared builders for the structurally identical tables.

## What it produces

- `data/panel_5yr.parquet` and the annual, 10-, 20-, 25-, and 50-year panels:
  the compressed, labelled datasets, holding only the variables the tables use.
- `output/table_2.csv` ... `table_7.csv`: the tables in a tidy long format, and
  `table_*.txt`, laid out like the page in the paper.
- `docs/replication_check.md`: the cell-by-cell comparison against the paper.
- `docs/published_values.csv`: the published numbers, transcribed from the paper.
- `docs/diff_table_*.csv`: the comparison behind the report.
- `docs/sessionInfo.txt`: the R and package versions a run used.

## Checking the results

The check covers 260 displayed cells across Tables 2-7 (coefficients and standard
errors to three decimals, R-squared and F-test p-values to two, sample sizes
exactly). Every cell matches the published paper, with one documented exception.

**Table 3, column 3 (Anderson-Hsiao), log GDP per capita.** The paper prints the
standard error as 0.127, which is the same value shown one row above it for
lagged democracy. The coefficient (-0.413) matches exactly, and three independent
methods (`estimatr` with Stata clustering, `estimatr` with CR0, and `AER::ivreg`
with a clustered covariance) all give a standard error of about 0.163. The
printed 0.127 appears to be a typesetting duplication of the row above; this code
reports 0.163.

## Two details worth knowing

**Lags.** The authors set the time index with `tsset code_numeric year_numeric`,
where `year_numeric` advances by one each period (five years in the main panel,
ten/twenty/twenty-five/fifty in the others). Lags here are built on that index, on
the full panel before the sample restriction, so a 1960 observation can use its
1955 and 1950 values exactly as Stata does.

**Arellano-Bond GMM.** The authors use `xtabond2 ... noleveleq robust`. There is
no off-the-shelf R command that reproduces it here, because the spec mixes
uncollapsed GMM-style instruments with a level "passthru" instrument and a
one-step robust covariance. `fit_abgmm()` in `R/00_setup.R` builds that estimator
directly; it matches every published GMM coefficient and standard error to three
decimals.
