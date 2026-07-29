# Other ways to pool the instruments

Professor Torgovitsky asked whether Roodman describes other ways to collapse
the instruments, since there might be several ways to aggregate them. This is
what I found in the paper and what happens when the alternatives are tried on
this data. The numbers are in output/aggregation.txt and output/aggregation.csv,
and the figure is output/aggregation.png.

## What Roodman actually describes

The paper is Roodman (2009), "A Note on the Theme of Too Many Instruments",
Oxford Bulletin of Economics and Statistics 71(1), 135-158. Full details for
everything cited here are in docs/references.md.

Section V is the relevant part, and it lays out two levers rather
than one. The first is capping how far back the lags go. The second is
collapsing, which he describes as squeezing the instrument matrix sideways and
adding together the columns that end up on top of each other. He gives one
collapsing rule: make a single instrument for each lag distance. He also points
out that pulling both levers at once leaves a count that does not grow with the
length of the panel at all, and his Table 1 crosses the two levers the same way
the figure in docs/instruments.md does.

Two further ideas appear in his footnotes rather than the main text. One is an
approach from Arellano (2003) that first models the instruments jointly with a
vector autoregression and uses the fitted coefficients as constraints, which
Roodman notes had not caught on. The other is drawing random subsets of the
available instruments repeatedly to see how much the answer moves around.

## The alternatives I tried

Collapsing is really just a decision about which instrument columns get added
together, so there is more than one way to do it. Roodman adds across time and
keeps the lag distances apart. The obvious alternative is to do it the other
way round, adding across lag distances and keeping the periods apart, which
gives one instrument per period. Adding and averaging are different here,
because the number of lags available grows over time, so averaging is a third
option. Pushing it to the end gives a fourth, where everything collapses into a
single instrument. I ran all four against the uncollapsed set, and then two
further families with a knob on them, which are in the last section.

I also ran two instrument designs, because they turn out to matter. The paper
builds its lag blocks from democracy only and gives income a single lagged
level. The symmetric alternative gives both variables the full block of lags,
which is what the plm package does by default and what the sweep in
docs/instruments.md uses. The first design reproduces the published GMM column
exactly, and the code stops if it ever fails to.

## What came out

The pooling rule matters, and by more than I expected. Using every lag under
the paper's design, the estimates for Freedom House range from -0.243 to -0.118
depending only on how the instruments are pooled, and for Polity from -0.553
to -0.351. Same data, same lags, same estimator.

Under the paper's own instrument design, shown in the top row of the figure,
the uncollapsed set is the odd one out. It is the only line that climbs
steadily toward the fixed-effects value as the lag window widens (a range of
0.084 for Freedom House and 0.245 for Polity), while the four pooled rules
all stay inside 0.044. That is the overfitting story again.

The bottom row, where income also gets its own block of lags, is messier and
worth saying plainly rather than glossing. There the pooled rules are not
uniformly steadier. Roodman's rule is still the flattest line in every panel of
the figure, with a range of at most 0.054, but pooling by period swings more
than the uncollapsed set does (0.477 against 0.332 for Polity), and pooling
everything together is the least stable line on the whole figure. So the useful
statement is not that pooling always steadies the estimate. It is that pooling
by lag distance, which is the rule Roodman actually proposes, is the one that
holds up under both designs.

Pooling too hard has its own failure mode. With the symmetric design, pooling
everything together produces a standard error of 1.61 for Freedom House, which is
another way of saying the instrument has stopped carrying usable information.
There is a middle ground. Pool too little and the estimator overfits, pool too
much and there is not enough left to identify anything. That second failure is
the weak instrument problem, which is what docs/weak-instruments.md takes up.

## Turning the dial by degrees

The five rules above are all-or-nothing. Two more ways are not single rules but
families with a knob on them, which lets you ask what happens part of the way.
They are in output/aggregation_families.txt and output/aggregation_families.png.

The first groups the years into blocks and pools within a block, so a block size
of one is the uncollapsed set and a block as wide as the panel is Roodman's
rule. The second keeps one instrument per year but multiplies each older lag by
a fading factor before adding, so a factor of one is the plain sum. Both
families have to land exactly on the fixed rules at their endpoints, and the
code stops if they do not.

Neither did what I expected, in opposite directions.

Blocking is not a smooth dial at all. The estimates bounce around and travel
well outside the two endpoints they connect. For Freedom House they run from
-0.274 to -0.105, while the endpoints themselves are only -0.129 and -0.154.
Even the instrument count refuses to fall as the blocks get wider. Block sizes
of four and five both give three blocks, but the wider blocks reach back to
years that have deeper lags on offer, so the count goes up rather than down.
Block size is therefore not a measure of how much pooling is going on, which is
why the line looks like noise.

Fading the older lags does almost nothing. Across the whole range of factors
the estimate moves by 0.032 for Freedom House and 0.015 for Polity.
That is worth knowing rather than disappointing: it says the period-pooled
answer is not an artifact of counting every lag equally, since counting the
older ones for a tenth as much barely moves it.

One note on reading the table against the figure. The table's second column
uses every available lag, which is past the right-hand edge of the figure; the
figure stops at a window of 8. The instrument counts in the table count only
the lagged levels, not the year dummies that also sit in the instrument set.
