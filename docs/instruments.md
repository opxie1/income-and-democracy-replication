# Collapsed instruments and a wider lag window

Professor Torgovitsky asked me to add collapsing as a second dial on the
sweep from last round. The lag window and the collapsing choice now move
together in one picture. These estimators use the past of a country as a
stand-in for its present, so the instruments are old values of democracy and
income. Without collapsing, the estimator gets one instrument for each pair
of a lag and a time period. With collapsing, it gets one instrument per lag
distance, so the count is much smaller. I ran both, for lag windows of 2
through 8, for both GMM estimators and both democracy measures.

A third dial matters more than the lag window and the collapsing choice, so
it is on the figure too. The sweep from last round used the two-step
weighting matrix. Last round I did not see that weighting step as a choice.
This round the weighting step is a dial on the figure. It is the choice that
decides the answer. Results are in output/instruments.txt,
output/instruments.csv and output/instruments_drift.csv, and the figure is
output/instruments.png.

Two notes on how to read the output. Within each estimator the sample is
fixed, so the rows of one block in output/instruments.txt differ only in the
instrument set. Difference and system GMM do not use the same sample as each
other. For Freedom House they use 124 and 133 countries, because system GMM
can also use the level equations. Their levels are therefore not strictly
comparable, and only the trends carry information.

The instrument count also increases quickly with the window, but not without
limit. For difference GMM it runs 21, 33, 43, 51, 57, 61, 63 across windows
of 2 to 8. Once the window is as long as the panel, the count is flat.

## What happened

Professor Torgovitsky expected the collapsed line to stay flat while the
uncollapsed line drifted toward OLS or fixed effects. Under one-step
weighting that is what happens, for both measures. Under two-step weighting
it happens for one measure and not for the other. The sweep from last round
used two-step weighting. So the puzzle I reported then is a fact about the
weighting step and not about collapsing.

One-step difference GMM comes first, in the top row of the figure and its two
one-step columns. Between the narrowest window and the widest, the collapsed
estimate moves by 0.033 for Freedom House against 0.201 uncollapsed. The
collapsed estimate moves by 0.091 for Polity against 0.426. That is a factor
of 6 and 5. The uncollapsed lines drift toward the fixed-effects values of
0.010 and -0.006. The collapsed lines do not.

Two-step weighting breaks that pattern for Freedom House. There the collapsed
estimate moves by 0.179 against 0.169 uncollapsed, so collapsing does not
slow the drift at all. The collapsed line runs from -0.234 up to -0.055. The
uncollapsed line runs from -0.253 to -0.085. For Polity, two-step weighting
behaves like one-step, with a collapsed move of 0.027 against 0.094
uncollapsed.

The reason is the two-step weighting matrix. The estimator builds that matrix
from the residuals, and the matrix has one row and one column per instrument.
With a wide window and a few dozen countries, the data are too thin to
estimate the matrix well. The matrix then absorbs noise of its own.
Collapsing makes that matrix smaller, but not small enough here.

So the honest summary is that collapsing does its advertised job. The
two-step weighting step is a separate form of overfitting, and collapsing
does not fix it. A reader who looks only at the two-step columns will blame
collapsing for a problem that belongs somewhere else.

System GMM barely moves under any setting. It starts a little above the
pooled OLS value and stays there. The robust overidentification test also
rejects it for Polity everywhere. The largest p-value across the whole
two-step sweep is 0.047. So those level conditions look invalid on this data,
whatever the instrument count is.

## The standard errors

The standard errors tell the same story from the other side, and they do it
under both weighting choices.

In the uncollapsed one-step columns, the Freedom House standard error
decreases from 0.173 to 0.101 as the instrument count increases. The Polity
standard error decreases from 0.249 to 0.092. In the collapsed one-step
columns, both standard errors stay near their start. The Freedom House
standard error moves from 0.147 to 0.127, and the Polity standard error moves
from 0.202 to 0.210. The uncollapsed estimator therefore looks sharper as the
window widens. The instrument count manufactures that sharpness, and the data
do not earn it.

The weak-instrument work in docs/weak-instruments.md gives the same reading
for the two-stage least squares columns. A conventional standard error
assumes the instruments have real pull, and it decreases whether or not they
do.

The instrument counts show the scale of this effect. For difference GMM, the
widest window uses 21 instruments in the collapsed set and 63 instruments in
the uncollapsed set.

## The overidentification test

One caution comes before the numbers, and it took me a while to get this
right. The plm package reports an overidentification test after either
estimator. By default it builds the statistic from the robust weight matrix
for both estimators. So the number it prints after a one-step fit is not the
criterion of that estimator. I report the criterion of each estimator
instead.

The Sargan statistic sits in the one-step rows, and it assumes homoskedastic
errors. The Hansen statistic sits in the two-step rows. It is robust to
heteroskedasticity, but it loses power as the instrument count increases.
They are different tests, so a reader must compare p-values down a column and
not across the weighting rows.

The one-step rows say nothing useful, and the reason matters more than the
numbers. The Sargan statistic rejects every specification in the sweep, and
the largest p-value anywhere in that column is 0.006. The test rejects the
collapsed sets and the uncollapsed sets, at every window, for both measures.
That is what a homoskedasticity-assuming test does on data with
heteroskedastic errors. The errors here are heteroskedastic enough that the
paper clusters all of its own inference by country.

The spread within that column is enormous and covers many orders of
magnitude. The test breaks its own assumption everywhere, so I do not think
that spread carries any weight.

The two-step rows carry the robust Hansen statistic, and that is the column
to read. There the uncollapsed p-value increases from 0.01 to 0.17 for
Freedom House, and from 0.01 to 0.13 for Polity. At the same time the
instrument count increases from 21 to 63. That is one statistic applied to a
larger and larger instrument set, so the comparison is clean. The test loses
power, and the instruments do not improve. So a high Hansen p-value on a
large instrument set is not reassurance.

The collapsed sets behave differently. Their p-values stay high for most of
the window, and then drop sharply at the widest window (0.25, 0.37, 0.59,
0.62, 0.71, 0.04 for Freedom House). I do not have a clean account of that
last drop, and it is the weakest link in this section.

The second-order serial correlation test matters most for whether these
instruments are valid at all. It passes everywhere. The smallest p-value
across all 112 fits is 0.17, so nothing here suggests that the lag-2
instruments are invalid.

## How this fits with the other sweep

docs/aggregation.md runs the same idea through the one-step estimator written
for the replication, and it reaches a tidier conclusion. Last round I did not
know which of two differences explained the gap, the weighting step or the
sample. The weighting dial answers that question, and the answer is the
weighting step. The sample and the instrument design stay fixed, and only the
weighting step changes. That one change turns the Freedom House collapsed
line from a 0.179 drift into a 0.033 one. The samples still differ between
the two sweeps, but I no longer need that difference to explain anything.

## Other ways to collapse

Professor Torgovitsky also asked whether Roodman describes other ways to
collapse. Section V of the paper gives two levers rather than one, and one
collapsing rule. Collapsing is a choice about which instrument columns the
estimator adds together, so other rules exist and deserve a test. That
comparison is in docs/aggregation.md.
