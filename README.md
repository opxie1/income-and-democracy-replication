# Income and Democracy: a replication

This project redoes, in R, the main results of a famous economics paper, and checks that every number comes out the same as in the published version:

> Acemoglu, Daron, Simon Johnson, James A. Robinson, and Pierre Yared. 2008. "Income and Democracy." American Economic Review 98 (3): 808-842.

The paper asks whether countries that grow richer end up more democratic. I rebuilt the six tables (Tables 2 through 7) where the authors try hardest to tell whether income really causes democracy, and I kept the simpler estimates next to them so they are easy to compare.

## What you need

1. R, version 4 or newer. On Windows you can install it with `winget install RProject.R`.
2. The authors' data. Download it from [openICPSR project 113251](https://www.openicpsr.org/openicpsr/project/113251) and unzip it into a folder named `replication-kit` in this project.
3. A few R packages:

```r
install.packages(c("tidyverse", "arrow", "readxl", "estimatr", "plm", "here"))
```

## Running it

```sh
Rscript R/run_all.R
```

This reads the data, rebuilds each table, and compares every value to the paper. The whole run takes under a minute, and it prints how many numbers matched at the end. On Windows, R sometimes shows an error code on the last line. That is harmless. It comes from the `arrow` package closing, and the files are already saved by then.

## What you get

- The rebuilt tables in `output/`, both as spreadsheets and laid out the way they look in the paper.
- The cleaned data in `data/`.
- A short report, `docs/replication_check.md`, that compares every value to the published number.

## How close it is

All 260 published numbers match, down to the last printed digit. There is one exception. In Table 3, the paper prints a standard error of 0.127, but that is just a copy of the number in the row right above it. The correct value is 0.163, which three different methods all give, so the code reports 0.163 and points out the difference.

## Trying other methods

After the replication, I re-estimated the effect of income on democracy a few other ways and put them next to each other. The table is in `output/alternatives.txt`, and a plain writeup is in `docs/alternatives.md`. Here is what I found. The methods that work by comparing changes within a country all agree with the paper that income has little or no positive effect on democracy. The only method that gives a positive effect leans on an extra assumption, and once you test that assumption it does not hold up.

I checked these estimates two ways. The code reproduces a standard textbook result exactly, and it stops if it ever fails to. A second package gives the same answer, and that check is in `R/11_crosscheck.R`, which needs the `pdynmc` package.
