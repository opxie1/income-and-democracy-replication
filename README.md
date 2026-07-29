# Income and Democracy: a replication

This project redoes, in R, the main results of a famous economics paper, and checks that every number comes out the same as in the published version:

> Acemoglu, Daron, Simon Johnson, James A. Robinson, and Pierre Yared. 2008. "Income and Democracy." American Economic Review 98 (3): 808-842.

The paper asks whether countries that grow richer end up more democratic. I rebuilt the six tables (Tables 2 through 7) where the authors try hardest to tell whether income really causes democracy, and I kept the simpler estimates next to them so they are easy to compare.

## What you need

1. R, version 4 or newer. On Windows you can install it with `winget install RProject.R`.
2. The authors' data. Download it from [openICPSR project 113251](https://www.openicpsr.org/openicpsr/project/113251) and unzip it into a folder named `replication-kit` in this project.
3. A few R packages:

```r
install.packages(c("tidyverse", "arrow", "readxl", "estimatr", "plm", "here", "MASS"))
```

## Running it

```sh
Rscript R/run_all.R
```

This reads the data, rebuilds each table, and compares every value to the paper. The whole run takes under a minute, and it prints how many numbers matched. On Windows, R sometimes shows an error code on the last line. That is harmless. It comes from the `arrow` package closing, and the files are already saved by then.

## What you get

- The rebuilt tables in `output/`, both as spreadsheets and laid out the way they look in the paper.
- The cleaned data in `data/`.
- A short report, `docs/replication_check.md`, that compares every value to the published number.

## How close it is

Of the 260 published numbers, 259 match down to the last printed digit. The one that does not is a misprint in the paper: in Table 3 it prints a standard error of 0.127, a copy of the number in the row right above it. The correct value is 0.163, which three different methods all give, so the code reports 0.163 and points out the difference.

## Trying other methods

After the replication, I re-estimated the effect of income on democracy a few other ways and put them next to each other. The table is in `output/alternatives.txt`, and a plain writeup is in `docs/alternatives.md`. Here is what I found. The methods that work by comparing changes within a country all agree with the paper that income has little or no positive effect on democracy. The only method that gives a positive effect leans on an extra assumption. The data reject that assumption outright for one democracy measure, and for the other the positive result still hangs entirely on it.

I checked these estimates two ways. The code reproduces a standard textbook result exactly, and it stops if it ever fails to. A second package gives the same answer, and that check is in `R/11_crosscheck.R`, which needs the `pdynmc` package. Run it with `Rscript R/11_crosscheck.R`.

## Pushing on the instruments

These GMM methods let you use many past values as instruments, and using too many is a known trap. I re-ran both GMM estimators over a widening window of lags, once with the instruments left separate and once with them pooled together, to see how fast the estimates decay toward the naive ones. Pooling holds one of the two democracy measures steady but not the other, so it is not the whole answer. The clearer effect is on the error bars: left separate, they shrink by roughly two thirds as the instruments pile up, which is precision manufactured by the instrument count rather than earned from the data. The writeup is in `docs/instruments.md`.

Pooling can be done in more than one way, so I compared the obvious alternatives against each other in `docs/aggregation.md`. The rule matters. On the same data with the same lags, the answer for one measure moves by about a factor of two depending only on how the instruments are combined. I also ran two versions that can be turned up and down by degrees rather than switched on and off, which is where the more interesting result is: fading out the older lags barely changes anything, so that answer does not depend on counting every lag equally.

## Are the instruments strong enough

Instrumental-variables estimates go wrong quietly when the instrument barely moves the thing it stands in for. For the two-stage least squares columns I report first-stage F statistics and then two confidence sets built to stay valid even when the instrument is weak, one from Anderson and Rubin and one from Moreira. The instruments hold up reasonably, but the honest confidence sets are wider than the ordinary ones and all of them include zero. The writeup is in `docs/weak-instruments.md`.

I wrote those two tests by hand, so `R/15_ivcrosscheck.R` re-runs every column through the `ivmodel` package and stops if the answers disagree. It needs that package, which is why it sits outside the main run.

## Sources

Every paper and package this project leans on is listed in `docs/references.md`.

## The two optional scripts

`Rscript R/run_all.R` runs the replication and everything built on top of it. Two scripts sit outside it because each needs a package that is not otherwise required, and both are independent checks rather than part of the pipeline: `R/11_crosscheck.R` (needs `pdynmc`) and `R/15_ivcrosscheck.R` (needs `ivmodel`).
