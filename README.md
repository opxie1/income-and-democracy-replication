# Income and Democracy: a replication

This project rebuilds, in R, the main results of a well-known economics paper and checks that every number comes out exactly as published:

> Acemoglu, Daron, Simon Johnson, James A. Robinson, and Pierre Yared. 2008. "Income and Democracy." American Economic Review 98 (3): 808-842.

The paper asks whether countries that grow richer go on to become more democratic. This code reproduces the six tables (Tables 2 through 7) where the authors work hardest to separate cause from coincidence, keeping the simpler comparison estimates that sit beside them.

## What you need

1. R, version 4 or newer. On Windows you can install it with `winget install RProject.R`.
2. The authors' data. Download it from [openICPSR project 113251](https://www.openicpsr.org/openicpsr/project/113251) and unzip it into a folder named `replication-kit` in this project.
3. A handful of R packages:

```r
install.packages(c("tidyverse", "arrow", "readxl", "estimatr", "plm", "here"))
```

## Running it

```sh
Rscript R/run_all.R
```

It reads the data, rebuilds each table, and compares every value against the paper. The whole run takes under a minute and ends by printing how many numbers matched. On Windows, R sometimes shows a harmless error code on the last line; that comes from the `arrow` package shutting down, and all the files are saved by then.

## What you get

- The rebuilt tables in `output/`, both as spreadsheets and laid out the way they appear in the paper.
- The cleaned data in `data/`.
- A short report, `docs/replication_check.md`, that compares every value to the published figure.

## How close it is

All 260 published numbers match, down to the last printed digit. The one exception is a single standard error in Table 3: the paper prints 0.127, which is a repeat of the number in the row directly above it, while the correct value (which three independent methods agree on) is 0.163. The code reports 0.163 and flags the difference.

## Trying other methods

Beyond reproducing the paper, the code re-estimates the income effect with a few alternative methods and lays them side by side. The table is in `output/alternatives.txt`, with a plain-language writeup in `docs/alternatives.md`. The short version: the methods that work by comparing changes within a country all agree with the paper that income has little or no positive effect on democracy, and the only way to bring a positive effect back is to use a method that leans on an extra assumption.
