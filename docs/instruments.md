# Collapsing the instruments, and letting the lag window grow

Professor Torgovitsky asked me to take the sweep from last round and make
collapsing a second dial, so the lag window and the collapsing choice move
together in one picture. These estimators use a country's own past as a
stand-in for its present, so the instruments are old values of democracy and
income. Left separate, the estimator gets its own instrument for every pairing
of a lag with a time period. Collapsed, it gets one instrument per lag
distance, so the count stays much smaller. I ran both, for lag windows of 2
through 8, for both GMM estimators and both democracy measures. Results are in
output/instruments.txt and output/instruments.csv, and the figure is
output/instruments.png.

Two notes on reading it. Within each estimator the sample is fixed, so only the
instrument list changes down a column, but difference and system GMM do not use
the same sample as each other: for Freedom House they use 124 and 133 countries
respectively, because system GMM can also use the level equations. Their levels
are therefore not strictly comparable, only their trends. Also, the instrument
count grows quickly with the window but not without limit. For difference GMM
it runs 21, 33, 43, 51, 57, 61, 63 across windows of 2 to 8,
flattening out once the window reaches back as far as the panel goes.

## What happened

Half of the expectation held and half of it did not, so it is worth separating
the two measures.

For Polity, collapsing does what it is supposed to. The collapsed estimate
starts at -0.429 and ends at -0.457, drifting slightly away from zero rather
than toward it, while the uncollapsed one climbs from -0.343 to -0.249, moving
toward the fixed-effects value of -0.006.

For Freedom House it does not. Both lines drift, and by almost the same amount:
collapsed runs from -0.234 to -0.055 and uncollapsed from -0.253 to -0.085, both
heading toward the fixed-effects value of 0.010. Collapsing slows the drift
here but does not stop it. The reason is that collapsing does not freeze the
instrument count, it only slows its growth, so a wide enough window still piles
up enough instruments to bite.

System GMM barely moves under either setting. It starts a little above the
pooled OLS value and stays there, collapsed or not.

The standard errors are the clearer signal. Uncollapsed, the Freedom House
standard error falls from 0.186 to 0.059 as the instruments pile up, and the
Polity one from 0.238 to 0.098. Collapsed, they hold up much better: 0.147 to
0.102 and 0.202 to 0.192. The uncollapsed estimator looks like it is getting
sharper as the window widens, and that apparent sharpness is manufactured by
the instrument count rather than earned from the data.

The counts show the scale of it. For difference GMM the widest window uses 21
instruments collapsed against 63 uncollapsed.

The overidentification test broadly loses its bite as the count grows, which is
the test running out of power rather than the instruments improving, but it is
not a clean trend and I do not want to oversell it. Uncollapsed for Freedom
House it goes 0.01, 0.11, 0.13, 0.20, 0.16, 0.19, 0.17.
Collapsed, the narrowest window is exactly identified and leaves nothing to
test; from there the p-value rises most of the way, 0.25, 0.37, 0.59, 0.62, 0.71,
before dropping back sharply to 0.04 at the widest window, which is also where
that estimate makes its largest jump. I do not have a clean account of that
last point, and it is the weakest link in the Freedom House story above.

The test that matters most for whether these instruments are allowed at all is
the second-order serial correlation check, and it passes everywhere: the
smallest p-value across all 56 fits is 0.17, so nothing here suggests the lag-2
instruments are invalid.

One caveat on reading this next to docs/aggregation.md, which runs the same
idea through the hand-built estimator and reaches a tidier conclusion. The two
differ in three ways at once: this sweep uses two-step GMM through the plm
package and that one uses the one-step estimator written for the replication,
the samples differ (this one keeps fewer countries), and the levels differ
noticeably as a result. I checked whether the instrument design was the
explanation and it is not, since the symmetric design stays flat there too. So
the honest summary is that the Freedom House drift under collapsing shows up
with the two-step plm estimator on its sample and not with the one-step
estimator on the paper's sample, and I have not isolated which of those two
differences is doing the work.

## Other ways to aggregate

He also asked whether Roodman describes other ways to collapse. Section V of
the paper lays out two levers, not one, and gives a single collapsing rule.
Since collapsing is really just a choice about which instrument columns get
added together, there are other rules worth trying. That comparison is in
docs/aggregation.md.
