# Income and Democracy: a replication

This project rebuilds the main results of a famous economics paper in R. It checks every number against the published version. The paper is:

> Acemoglu, Daron, Simon Johnson, James A. Robinson, and Pierre Yared. 2008. "Income and Democracy." American Economic Review 98 (3): 808-842.

The paper asks whether countries that grow richer become more democratic. I rebuilt the six tables (Tables 2 through 7) where the authors work hardest on whether income causes democracy. I kept the simpler estimates next to them for comparison.

## What you need

1. R, version 4 or newer. On Windows, install it with `winget install RProject.R`.
2. The data from the authors. Download it from [openICPSR project 113251](https://www.openicpsr.org/openicpsr/project/113251). Then unzip it into a folder named `replication-kit` in this project.
3. Install these R packages:

```r
install.packages(c("tidyverse", "arrow", "readxl", "estimatr", "plm", "here", "MASS"))
```

## How to run it

```sh
Rscript R/run_all.R
```

This command reads the data, rebuilds each table, and checks every value against the paper. It then runs the five follow-up studies described below. It prints how many numbers matched.

The replication itself takes a few seconds. The whole run takes about an hour. Almost all of that time goes to two simulation studies: the random-instrument draws and the Monte Carlo, both described further down.

On Windows, R sometimes shows an error code on the last line. That error is harmless. It comes from the `arrow` package as the package closes. The files are already saved at that point.

## What you get

- The rebuilt tables in `output/`, as spreadsheets and in the layout of the paper.
- The cleaned data in `data/`.
- A short report, `docs/replication_check.md`, that checks every value against the published number.

## How close it is

Of the 260 published numbers, 259 match to the last printed digit. The one number that does not match is a misprint in the paper. Table 3 prints a standard error of 0.127. That value is a copy of the number in the row above it. Three different methods all give 0.163 as the correct value. The code reports 0.163 and points out the difference.

## Other methods

After the replication, I re-estimated the effect of income on democracy in a few other ways. I put the results next to each other. The table is in `output/alternatives.txt`. A plain writeup is in `docs/alternatives.md`.

Here is what I found. Some of the methods compare changes within a country. These methods all agree with the paper that income has little or no positive effect on democracy. Only one method gives a positive effect, and that method depends on an extra assumption. The data reject that assumption outright for one democracy measure. For the other measure, the positive result still depends entirely on the assumption.

I checked these estimates two ways. The code reproduces a standard textbook result exactly. If the result ever fails to match, the script stops. A second package gives the same answer. This second check is in `R/11_crosscheck.R`, which needs the `pdynmc` package. The command `Rscript R/11_crosscheck.R` runs it.

## Too many instruments

These GMM methods allow many past values as instruments. Too many instruments is a known trap. I re-ran both GMM estimators over a series of lag windows, each one wider than the last. I ran each estimator twice: once with the instruments left uncollapsed, and once with them collapsed.

Collapsing holds the estimate still for both democracy measures. This result depends on one-step weighting of the conditions. Under two-step weighting, collapsing no longer holds the estimate still for one of the two measures. Last round I used two-step weighting without treating it as a choice. That result closes a loose end that I left open last time.

The error bars tell the same story from the other side. With the instruments uncollapsed, the error bars narrow by roughly half as the instrument count grows. That is precision manufactured by the instrument count rather than earned from the data. The writeup is in `docs/instruments.md`.

More than one rule can collapse the instruments. I compared the obvious alternatives against each other in `docs/aggregation.md`. The rule matters. On the same data with the same lags, the answer for one measure moves by about a factor of two. Only the collapsing rule differs between those two answers.

The comparison has a clear winner. It is Roodman's own rule. The overidentification test rejects that rule in 1 of 28 cells. I built two other rules by turning his construction on its side. The same test rejects them in 27 and 32 cells out of 32.

I also ran two versions that change by degrees rather than by an on-off switch. A fading weight on the older lags moves the estimate by at most 0.032, so that answer does not depend on counting every lag equally. Roodman suggests one more version in a footnote, and I missed it earlier. That version draws the instruments at random. On this data and this specification, with few instruments in hand, the particular draw decides even the sign of the estimate.

## Are the instruments strong enough?

When the instrument barely moves the variable that it stands for, an instrumental-variables estimate can be badly wrong with no warning. For the two-stage least squares columns, I report first-stage F statistics. One of these statistics is the version that is justified for the way the paper clusters the errors. I also report confidence sets that stay valid with a weak instrument, from Anderson and Rubin and from Moreira.

Moreira's test in its usual form assumes that the errors are well behaved. So I also built the version that allows for clustering. That version settles a question that I left as a guess last time.

These instruments are not badly weak, but several are not comfortably strong either. Four of the ten columns fall short of the Stock and Yogo threshold, and five fall short of the Montiel Olea and Pflueger value. The honest confidence sets are wider than the ordinary ones, and all ten include zero. The writeup is in `docs/weak-instruments.md`.

I wrote those tests by hand. For that reason, `R/15_ivcrosscheck.R` re-runs every column through the `ivmodel` package. If the answers disagree, the script stops. It needs that package, so it sits outside the main run.

## A simulation with known answers

The estimates above disagree with each other, and the data cannot say which one is right. So Professor Torgovitsky suggested building simulated data where the truth is known. I simulate a dynamic panel that reproduces the observation pattern of the real one cell by cell, so the estimators get the sample they have in the data: 838 observations on 127 countries for Freedom House. Persistence comes from the estimates here. The true effect of income on democracy is set to zero, so every reported effect is an error that can be measured. The study draws 500 panels for each of eight designs and runs the whole ladder of estimators on each. The writeup is in `docs/monte-carlo.md`.

Four results stand out. Pooled OLS reports a positive income effect in every draw, and its interval never covers the true zero. Fixed effects gets persistence wrong by 0.215, which is the known bias of that method with few periods. Collapsing the instruments repairs the persistence estimate: the error falls from 0.098 to 0.027, and the share of intervals that cover the truth rises from 0.73 to 0.96. System GMM is accurate when countries start at their long-run average and badly wrong when they do not, with coverage falling from 0.96 to 0.45.

The last of those matches what the overidentification test says about the real data. The simulation puts a number on the cost of ignoring it.

## Sources

The file `docs/references.md` lists every paper and package that this project uses.

## The two optional scripts

The command `Rscript R/run_all.R` runs the replication and everything built on it. Two scripts sit outside that run. Each one needs a package that the main run does not need. Both are independent checks rather than parts of the main run. They are `R/11_crosscheck.R`, which needs `pdynmc`, and `R/15_ivcrosscheck.R`, which needs `ivmodel`.
