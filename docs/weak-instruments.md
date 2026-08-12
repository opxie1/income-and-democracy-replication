# Are the instruments strong enough?

Professor Torgovitsky asked for weak instrument diagnostics, then for
Anderson-Rubin confidence sets and the Moreira conditional likelihood ratio
approach. The natural place for all three is the two-stage least squares
columns of Tables 5 and 6. Those columns instrument one variable, income,
with a small number of instruments. The instrument is the savings rate in
Table 5 and trade-weighted world income in Table 6. That gives ten columns.

The numbers are in output/weakiv.csv and output/weakiv.txt. The file
output/weakiv.png draws the p-value curves. The tests come from Anderson and
Rubin (1949), Moreira (2003) and Kleibergen (2005). The first-stage cutoffs
come from Stock and Yogo (2005) and Montiel Olea and Pflueger (2013). Full
details are in docs/references.md.

The section "What this does not reach" covers the paper's other
instrumental-variables columns. They are the Anderson-Hsiao and Arellano-Bond
columns, which instrument two variables at once. The confidence sets here
have no direct counterpart there. I report first-stage diagnostics for those
columns.

## Why this matters

A weak instrument barely moves the variable that it stands in for. With a
weak instrument, the usual standard errors and confidence intervals are not
trustworthy. They rest on the assumption that the instrument has real pull.
The tests below stay honest whether the instrument is strong or weak. A
comparison against the ordinary interval shows how much the ordinary interval
relied on that assumption.

## How strong is the first stage

The first-stage F statistics run from 12.1 to 39.5. I clustered them by
country, the same way the paper clusters everything else. With one endogenous
regressor this is the cluster-robust first-stage Wald F that standard
software reports. Every column clears the old rule of thumb that an F of more
than 10 is good enough. Against the stricter Stock and Yogo cutoff, four of
the ten columns fall short.

That last comparison is not strictly legitimate. Last round I said so and
left it there. Stock and Yogo derive their cutoffs for well behaved,
unclustered errors. I cluster these F statistics by country. Under
clustering, the properly justified statistic is the effective F of Montiel
Olea and Pflueger (2013). This round I implemented it, and it is in the table
beside the ordinary one.

I have two observations about what the effective F changes. First, for the
eight columns with a single instrument, the effective F is algebraically
identical to the cluster-robust first-stage F. The script checks this
identity, and does not assume it. So nothing moves there.

Second, for the two columns with two instruments the effective F moves a lot,
and in opposite directions. Table 5 col 9 goes from 12.1 to 21.6. Table 6 col
9 goes from 17.8 to 13.2. The effective F divides by a trace that accounts
for how unevenly the instruments contribute. The ordinary F does not.

The two cutoffs are not the same kind of object, so I report them separately
rather than in one column. Stock and Yogo bound the size distortion of the
conventional test. Montiel Olea and Pflueger bound the bias of the estimator
relative to a worst case. Their conservative critical value of 23.11 holds
whatever the covariance structure turns out to be. Against that value, five
of the ten columns fall short.

Either way the reading is the same. These instruments are not badly weak, and
several are not comfortably strong. The two confidence sets below handle
exactly that situation.

## Anderson-Rubin

The Anderson-Rubin test asks a simple question. The test takes a candidate
value for the effect of income. It subtracts what that value implies, then
looks for anything left over that the instruments still explain. If something
is left, the test rejects that candidate value. Every value that survives
goes into the confidence set. That set stays valid no matter how weak the
instruments are.

I clustered the test by country to match the paper. That makes it exactly the
cluster-robust Wald test on the auxiliary regression. The script checks it
against an actual fit of that regression.

All ten sets are bounded, and all of them contain zero. The width of the grid
I searched does not decide whether a set is bounded. With this form of the
statistic, the Anderson-Rubin statistic tends to the first-stage F as the
candidate value runs off to infinity. So the set is bounded exactly when that
F exceeds the critical value. The script computes boundedness both ways. If
the two answers disagree, the script stops.

The first stage is strongest in Table 6 col 7, at an effective F of 39.5.
There the set is [-0.48, 0.06], which is nearly the ordinary interval of
[-0.46, 0.06]. The first stage is weakest in Table 6 col 9, at an effective F
of 13.2. There the set is [-0.71, 0.07], against an ordinary interval of
[-0.51, 0.08]. The Anderson-Rubin set is at least as wide as the ordinary
interval in every column. So the ordinary intervals are optimistic
throughout, not only at the extremes.

I must not push that comparison further than it goes. The columns do not form
a tidy order by first-stage strength. The widening ratio runs from 1.05 to
1.45. It does not decrease neatly as the first stage gets stronger.

Table 6 col 9 has the weaker first stage of the two, an effective F of 13.2
against 21.6. It still widens less than Table 5 col 9 does. The two widening
ratios are 1.32 against 1.45. These specifications differ in their controls
and samples as well as in instrument strength. So the two ends of the range
are suggestive rather than a clean relationship.

## Moreira's conditional likelihood ratio

The conditional likelihood ratio test can be sharper than Anderson-Rubin. It
does not compare the statistic against a fixed cutoff. Instead it adjusts the
cutoff with a quantity that carries the information about instrument
strength. That adjustment buys back power that Anderson-Rubin gives away. I
coded it directly, so one grid inversion serves all three tests here.
R/15_ivcrosscheck.R checks it against the ivmodel package, which implements
the same test.

The catch is that eight of the ten columns here have exactly one instrument.
With one instrument there is nothing to buy back, and the two tests coincide.
So the extra power only matters for the two columns with a pair of
instruments. For the single-instrument columns the cutoff is a chi-squared
quantile with an exact form. The script uses that form rather than a
simulation of a distribution it already knows.

Moreira's test assumes that the errors are well behaved and independent. The
packages implement it that way. The paper clusters, and the Anderson-Rubin
sets above cluster. So last round I reported the two side by side, with a
caveat. My caveat was that the gap between them was mostly the clustering
rather than the extra power. That was a conjecture.

This round I also implemented the cluster-robust version, which is
Kleibergen's (2005) generalization. It is the same conditional likelihood
ratio statistic, built from moment conditions. I estimate the variance of
those moment conditions by country cluster. All three sets are in
output/weakiv.txt.

The cluster-robust version settles the question. For the eight
single-instrument columns the cluster-robust statistic equals the clustered
Anderson-Rubin statistic. It must do so, and the script checks that equality
to ten decimal places at every candidate value on the grid. The two printed
sets still differ by at most 0.011 at an endpoint, because the two tests read
that one statistic off different reference distributions. The likelihood
ratio test uses chi-squared and the Anderson-Rubin test uses F. So for those
columns the entire gap between the homoskedastic set and the Anderson-Rubin
set was the clustering, as I guessed.

For the two columns with two instruments the answer is more interesting, and
it is not what I guessed. The cluster-robust version widens the set by a
factor of 1.38 and 1.25. The result is still inside the Anderson-Rubin set
rather than equal to it: [-0.21, 0.21] against [-0.28, 0.26], and [-0.58,
0.01] against [-0.71, 0.07]. That remaining difference is the extra power the
conditional likelihood ratio test buys, and it is real.

So the honest split is this. For eight columns the gap was all clustering.
For the two overidentified ones it is part clustering and part genuine power.

## What it adds up to

The inference here does not lean on strong instruments. Under that inference,
every one of the ten two-stage least squares columns leaves zero inside the
confidence set. So none of them establishes an effect of income on democracy
in either direction.

This does not overturn anything in the paper. The paper argues that the
effect is not there once the fixed differences between countries come out.
The result does mean that the instrumental-variables columns are
uninformative about the sign of the effect. They are not evidence for a
negative effect.

## What this does not reach

The paper's other instrumental-variables columns instrument two variables at
once, lagged democracy and lagged income. The machinery above does not
transfer here, for two reasons. A confidence set for the income coefficient
alone needs a projection of a two-dimensional region rather than a grid
inversion. Moreira also derives his test for a single endogenous regressor,
so there is no conditional likelihood ratio counterpart at all.

The Arellano-Bond columns add a third reason. Their instrument counts run
into the dozens against about a hundred countries. This sample cannot support
a cluster-robust weight matrix of that size. That last point is not a side
issue. It is the subject of docs/instruments.md, where the same instrument
counts inflate the apparent precision of the estimates.

The first-stage diagnostic does transfer. So I report it for the
Anderson-Hsiao columns, Table 2 col 3 (Freedom House) and Table 3 col 3
(Polity). With two endogenous regressors the right statistic is Cragg and
Donald's minimum eigenvalue rather than a per-equation F. The Stock-Yogo
cutoff changes with it, to 7.03. The Cragg-Donald statistics are 12.68 and
10.99, so both clear it. The per-equation F statistics, clustered by country,
are 34.9 and 26.3 for democracy and 10.9 and 9.4 for income.

The Cragg-Donald statistic is the homoskedastic one, so it carries the same
caveat as the Stock-Yogo comparison above. There is no clustered version of
the cutoff to hold it against. Taken together these are not the weakest
instruments in the paper. The reason those columns still produce very wide
intervals is the small number of moment conditions, not a first stage that
fails.

## Checks

This is fiddly enough that I did not want to rely on code I wrote myself.
Several checks run every time the script does. The partialling step must
reproduce the published two-stage least squares estimates. The fast
Anderson-Rubin curve must match the slower clustered regression that it
stands in for. Boundedness from the grid must agree with boundedness from the
first-stage F.

When there is one instrument, the cluster-robust conditional likelihood ratio
test must equal the clustered Anderson-Rubin test. The effective F must equal
the ordinary F in the same case. Every Anderson-Rubin set must contain its
own point estimate. It must also be at least as wide as the conventional
interval, which is the claim the section above leans on. If a check fails,
the script stops.

The script does not read confidence set endpoints from the grid. The grid
locates the crossing, and then the script solves for the endpoint. So the
printed set is the true set rather than the nearest grid point inside it.

R/15_ivcrosscheck.R also sends every column through the ivmodel package. If
the ivmodel conditional likelihood ratio sets disagree with mine, the script
stops. If the ivmodel two-stage least squares estimates do not match mine
exactly, the script also stops. At present both comparisons pass for all ten
columns. That script needs the ivmodel package, so it sits outside the main
run. The command for it is `Rscript R/15_ivcrosscheck.R`.
