# Other ways to pool the instruments

Professor Torgovitsky asked whether Roodman describes other ways to collapse
the instruments, since there might be several ways to aggregate them. This is
what I found in the paper and what happens when the alternatives are tried on
this data. The numbers are in output/aggregation.txt and output/aggregation.csv,
and the figure is output/aggregation.png.

## What Roodman actually describes

Section V of the paper is the relevant part, and it lays out two levers rather
than one. The first is capping how far back the lags go. The second is
collapsing, which he describes as squeezing the instrument matrix sideways and
adding together the columns that end up on top of each other. He gives one
collapsing rule: make a single instrument for each lag distance. He also points
out that pulling both levers at once leaves a count that does not grow with the
length of the panel at all, and his Table 1 crosses the two levers the same way
the figure in docs/instruments.md does.

Two further ideas appear in his footnotes rather than the main text. One is an
approach from Arellano that first models the instruments jointly with a vector
autoregression and uses the fitted coefficients as constraints, which Roodman
notes had not caught on. The other is drawing random subsets of the available
instruments repeatedly to see how much the answer moves around.

## The alternatives I tried

Collapsing is really just a decision about which instrument columns get added
together, so there is more than one way to do it. Roodman adds across time and
keeps the lag distances apart. The obvious alternative is to do it the other
way round, adding across lag distances and keeping the periods apart, which
gives one instrument per period. Adding and averaging are different here,
because the number of lags available grows over time, so averaging is a third
option. Pushing it to the end gives a fourth, where everything collapses into a
single instrument. I ran all four against the uncollapsed set.

I also ran two instrument designs, because they turn out to matter. The paper
builds its lag blocks from democracy only and gives income a single lagged
level. The symmetric alternative gives both variables the full block of lags,
which is what the plm package does by default and what the sweep in
docs/instruments.md uses. The first design reproduces the published GMM column
exactly, and the code stops if it ever fails to.

## What came out

The pooling rule matters, and by more than I expected. Using every lag under
the paper's design, the estimates for Freedom House range from -0.249 to -0.118
depending only on how the instruments are pooled, and for Polity from -0.554
to -0.351. Same data, same lags, same estimator.

The uncollapsed set is the odd one out. In the figure it is the only line that
climbs steadily toward the fixed-effects value as the lag window widens; the
pooled versions stay far flatter. That is the same overfitting story as
before, and it shows up more sharply here than in the sweep in
docs/instruments.md because that sweep uses the symmetric design, where pooling
alone does not settle the Freedom House estimate down.

Among the pooled rules, Roodman's is the best behaved. The two period-based
rules give noticeably more negative estimates, and collapsing everything into a
single instrument is the least stable of the lot: with the symmetric design it
produces a standard error of 1.61 for Freedom House, which is another way of saying
the instrument has stopped carrying usable information. There is a middle
ground. Pool too little and the estimator overfits, pool too much and there is
not enough left to identify anything. That second failure is the weak
instrument problem, which is what docs/weak-instruments.md takes up.
