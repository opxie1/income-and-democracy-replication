# Collapsing the instruments, and letting the lag window grow

Professor Torgovitsky asked me to take the sweep from last round and make
collapsing a second dial, so the lag window and the collapsing choice move
together in one picture. These estimators use a country's own past as a
stand-in for its present, so the instruments are old values of democracy and
income. Left separate, the estimator gets its own instrument for every pairing
of a lag with a time period, and the count grows with the square of the lag
window. Collapsed, it gets one instrument per lag distance, and the count grows
in step with the window instead. I ran both, for lag windows of 2 through 8,
for both GMM estimators and both democracy measures. The countries and years
stay the same throughout; only the instrument list changes. Results are in
output/instruments.txt and output/instruments.csv, and the figure is
output/instruments.png.

## What happened

Half of the expectation held and half of it did not, so it is worth separating
the two measures.

For Polity, collapsing does what it is supposed to. The collapsed estimate
starts at -0.429 and ends at -0.457, drifting slightly away from zero rather
than toward it, while the uncollapsed one
climbs from -0.343 to -0.249, moving toward the fixed-effects value of -0.006.

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
instruments collapsed against 63 uncollapsed. The overidentification p-value
climbs as that count grows, which is the test losing its power to object rather
than the instruments getting better.

One caveat on reading this next to the other writeups. This sweep gives both
democracy and income their own block of lags, which is what the plm package
does by default. The paper builds its own GMM column differently, from lags of
democracy only. Under the paper's design, collapsing does hold the Freedom
House estimate steady, which is in docs/aggregation.md.

## Other ways to aggregate

He also asked whether Roodman describes other ways to collapse. Section V of
the paper lays out two levers, not one, and gives a single collapsing rule.
Since collapsing is really just a choice about which instrument columns get
added together, there are other rules worth trying. That comparison is in
docs/aggregation.md.
