# Income and Democracy: a replication

This project is about reproducing the main results of Acemoglu, Johnson, Robinson, and Yared (2008) in R, and then asking how much those results depend on the choices that the estimators require. The paper is:

> Acemoglu, Daron, Simon Johnson, James A. Robinson, and Pierre Yared. 2008. "Income and Democracy." *American Economic Review* 98 (3): 808-842.

The paper asks whether countries that grow richer become more democratic. I rebuild Tables 2 through 7, which are the tables where the authors work hardest on the causal question, and I keep the simpler estimates alongside them for comparison.

## What you need

1. R, version 4 or newer. On Windows, install it with `winget install RProject.R`.
2. The data from the authors. Download it from [openICPSR project 113251](https://www.openicpsr.org/openicpsr/project/113251) and unzip it into a folder named `replication-kit` in this project.
3. These R packages:

```r
install.packages(c("tidyverse", "arrow", "readxl", "estimatr", "plm", "here", "MASS"))
```

## How to run it

```sh
Rscript R/run_all.R
```

This reads the data, rebuilds each table, checks every value against the paper, and then runs the five follow up studies described below. It reports how many numbers matched.

The replication itself takes a few seconds. The whole run takes about an hour, and almost all of that hour goes to two simulation studies, the random instrument draws and the Monte Carlo, both of which are described further down.

On Windows, R sometimes prints an error code on the last line. That error is harmless. It comes from the `arrow` package as the package closes, and every file has already been written by that point.

## What you get

The rebuilt tables are in `output/`, as spreadsheets and in the layout of the paper. The cleaned data are in `data/`. A short report in `docs/replication_check.md` checks every value against the published number.

## How close it is

Of the 260 published numbers, 259 match to the last printed digit. The one number that does not match appears to be a misprint. Table 3 prints a standard error of 0.127, which is a copy of the number in the row above it. Three different methods all give 0.163 as the correct value. The code reports 0.163 and flags the difference.

## Other methods

After the replication, I estimate the effect of income on democracy in several other ways and put the results side by side. The table is in `output/alternatives.txt` and a plain description is in `docs/alternatives.md`.

The methods that compare changes within a country all agree with the paper that income has little or no positive effect on democracy. Only one method gives a positive effect, and that method depends on an extra assumption. For one democracy measure the data reject that assumption outright. For the other measure the positive result rests entirely on the assumption.

I check these estimates in two ways. The code reproduces a standard textbook result exactly, and the script stops if it ever fails to match. A second package reproduces the same pattern, with income negative under difference GMM and positive under system GMM. That check is in `R/11_crosscheck.R`, which needs the `pdynmc` package and runs with `Rscript R/11_crosscheck.R`.

## Too many instruments

These GMM methods allow many past values to serve as instruments, and too many instruments is a known trap. I run both GMM estimators again over a series of lag windows, each one wider than the last, and I run each estimator twice, once with the instruments left uncollapsed and once with them collapsed.

Collapsing holds the estimate still for both democracy measures. This result depends on one step weighting of the moment conditions. Under two step weighting, collapsing no longer holds the estimate still for one of the two measures. The sweep from the previous round used two step weighting without treating it as a choice, so this result closes a loose end that I had left open.

The error bars tell the same story from the other side. With the instruments uncollapsed, the error bars narrow by roughly half as the instrument count grows. That precision is manufactured by the instrument count rather than earned from the data. The description is in `docs/instruments.md`.

More than one rule can collapse the instruments, and I compare the obvious alternatives against each other in `docs/aggregation.md`. The rule matters. On the same data with the same lags, the answer for one measure moves by about a factor of two, and only the collapsing rule differs between those two answers.

The comparison has a clear winner, which is Roodman's own rule. The overidentification test rejects that rule in 1 of 28 cells. I build two other rules by turning his construction on its side, and the same test rejects them in 27 and 32 cells out of 32.

I also run two versions that change by degrees rather than by an on and off switch. A fading weight on the older lags moves the estimate by at most 0.032, so the answer does not depend on counting every lag equally. Roodman suggests one more version in a footnote, which I had missed earlier. That version draws the instruments at random. On this data and this specification, with few instruments in hand, the particular draw decides even the sign of the estimate.

## Are the instruments strong enough?

When an instrument barely moves the variable that it stands for, an instrumental variables estimate can be badly wrong with no warning. For the two stage least squares columns I report first stage F statistics, one of which is the version that is justified for the way the paper clusters its errors. I also report confidence sets that stay valid with a weak instrument, from Anderson and Rubin and from Moreira.

Moreira's test in its usual form assumes that the errors are well behaved, so I also implement the version that allows for clustering. That version settles a question that I had left as a guess.

These instruments are not badly weak, although several of them are not comfortably strong either. Four of the ten columns fall short of the Stock and Yogo threshold, and five fall short of the Montiel Olea and Pflueger value. The honest confidence sets are wider than the ordinary ones, and all ten of them include zero. The description is in `docs/weak-instruments.md`.

I wrote those tests by hand, so `R/15_ivcrosscheck.R` sends every column through the `ivmodel` package and stops if the answers disagree. It needs that package, so it sits outside the main run.

## A simulation with known answers

The estimates above disagree with each other, and the data cannot say which one is right. Professor Torgovitsky suggested building simulated data in which the truth is known. I simulate a dynamic panel that reproduces the observation pattern of the real one cell by cell, so that the estimators face the sample they face in the data, which is 838 observations on 127 countries for Freedom House. Persistence comes from the estimates here. The true effect of income on democracy is set to zero, so every reported effect is an error whose size can be measured. The study draws 500 panels for each of eight designs and runs the whole ladder of estimators on each one. The description is in `docs/monte-carlo.md`.

Four results stand out. Pooled OLS reports a positive income effect in every draw, and its interval never covers the true zero. Fixed effects misses persistence by 0.215, which is the known bias of that method with few periods. Collapsing the instruments repairs the persistence estimate, since the error falls from 0.098 to 0.027 and the share of intervals that cover the truth rises from 0.73 to 0.96. System GMM is accurate when countries start at their long run average and quite inaccurate when they do not, with coverage falling from 0.96 to 0.45.

That last result matches what the overidentification test says about the real data. The simulation puts a number on the cost of ignoring it.

## Sources

The file `docs/references.md` lists every paper and package that this project uses.

## The two optional scripts

`Rscript R/run_all.R` runs the replication and everything built on it. Two scripts sit outside that run, and each one needs a package that the main run does not need. Both are independent checks. They are `R/11_crosscheck.R`, which needs `pdynmc`, and `R/15_ivcrosscheck.R`, which needs `ivmodel`.
