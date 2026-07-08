# Letting the instrument set grow

Professor Torgovitsky asked how the GMM estimates behave when the instruments
are not collapsed and are built from longer and longer lags. These methods use
a country's own past as a stand-in for its present, so the instruments here
are old values of democracy and income. In docs/alternatives.md I kept that
set small on purpose. This time I did the opposite. Every instrument stays
separate, nothing is collapsed, and the lags run back 2 periods, then 3, and
so on out to 8, which is as far back as the 1960 to 2000 sample goes. The
countries and years stay identical the whole way; only the instrument list
changes. The results are in output/instruments.txt and output/instruments.csv,
and the picture is output/instruments.png.

## What happened

With the shortest list, difference GMM finds income pushing against democracy:
about -0.25 for Freedom House and -0.34 for Polity, and the Polity number is
nearly the paper's own GMM estimate. Growing the list wears this down. Freedom
House loses most of its estimate by the fourth lag and finishes near -0.09,
not far from the fixed effects value of 0.01. Polity fades more slowly, and
with every lag in use it still reads -0.25. System GMM barely reacts at all.
It starts just above the plain OLS value and stays there through every step.

One note for anyone comparing instrument counts with the replication. Here
democracy and income each get the full run of lags, while the paper's GMM
column gives income a single lag, so my counts top out higher, 63 against the
replication's 55.

## Why

Too many instruments is a known trap. Each extra lag adds another condition
for the estimator to satisfy, and past some point it stops correcting the bias
it was built to fix and starts copying it instead. Difference GMM slides
toward fixed effects. System GMM sits on OLS from the start because its extra
assumptions already lean that way. There is a standard test meant to catch bad
instruments, and it goes quiet as the pile grows; a rising p-value there means
the test is losing its power to object. The routine check that would flag a
deeper problem stays clean at every lag length, so the drift traces back to
the size of the instrument set.

So the answer to the original question: the cross-over is quick for Freedom
House, mostly finished by the fourth lag, and it never quite completes for
Polity.
