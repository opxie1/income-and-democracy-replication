# Other ways to collapse the instruments

Professor Torgovitsky asked whether Roodman describes other ways to collapse
the instruments. He asked because a collapse rule is a choice, and other
choices exist. This file reports what the paper says. It also reports what
the alternatives give on this data. The numbers are in output/aggregation.txt
and output/aggregation.csv. The figure is output/aggregation.png.

## What Roodman describes

The paper is Roodman (2009), "A Note on the Theme of Too Many Instruments",
Oxford Bulletin of Economics and Statistics 71(1), 135-158. Full details for
everything cited here are in docs/references.md.

Section V is titled "Techniques for reducing the instrument count", and it
opens with two techniques rather than one. The first technique uses only
certain lags instead of all available ones. It still makes separate
instruments for each period, but it caps the number per period. The count
then grows in proportion to the length of the panel, and not with its square.

Roodman describes this technique as a projection of the regressors onto the
full instrument set, with the coefficients on certain lags held at zero. He
attributes that reading to Arellano (2003b), his working paper on optimal
instrumental variables.

The second technique, which he calls the less common one, is "to combine
instruments through addition into smaller sets". Its advantage is "retaining
more information, since no lags are actually dropped". One sentence matters
most here, because it says what collapsing does. It "is equivalent to
imposing the constraint in projecting regressors onto HENR instruments that
certain subsets have the same coefficient". The picture he gives is of
"squeezing the matrix ... horizontally and adding together formerly distinct
columns". He gives one rule for the choice of columns: "collapsed instruments
are straightforward conceptually: one is made for each lag distance, with 0
substituted for any missing values".

Those two sentences license most of what follows. Where they do not, I say
so. A cap on the lags sets some of the projection coefficients to zero.
Collapsing sets some of them equal to each other.

The choice of which coefficients to set equal is free, and the paper makes
only one choice available. Any other partition of the instrument columns into
groups is another collapsing rule in the same sense. So it is fair to ask how
much the answer depends on his particular rule.

The section closes with a statement of what these techniques are for. The
statement deserves a quotation, because it sets the standard for this file.
They "provide the basis for some minimally arbitrary robustness and
specification tests for Difference and System GMM: cut the instrument count
in one of these ways and examine the behavior of the coefficient estimates
and Hansen and Difference-in-Hansen tests". So the coefficient on its own is
only half of what he asks for. The table in output/aggregation.txt now
carries an overidentification p-value in every cell, beside the estimate.

I do not report the Difference-in-Hansen half here. That test compares one
nested subset of instruments against the rest, and the rules compared here
are not nested inside one another. It belongs in the difference-versus-system
comparison, where the extra level conditions form a subset that I can add and
remove. Both docs/alternatives.md and docs/instruments.md make that
comparison.

He also notes that the two techniques work together. The pair leaves a count
that does not grow with the length of the panel at all. His Table 1 crosses
them. It shows four variants of system GMM: the full instrument set,
one-period lags only, the collapsed set, and both restrictions at once.

Two further ideas appear in footnotes rather than the main text. Footnote 6
describes an approach from Arellano (2003b). That approach first models the
instrumenting variables as a group, as functions of their collective lags,
with a vector autoregression. It then turns the fitted coefficients into
constraints on the projection. Roodman notes that it "has yet to enter common
practice".

Footnote 7 suggests "repeatedly selecting random subsets from the collection
of potential instruments and investigating how key results such as
coefficients of interest and the p value on the J statistic vary with the
number of instruments". That idea is cheap to run, so I ran it. It is the
last section here.

## The alternatives I tried

Roodman adds across time and keeps the lag distances apart. The obvious
alternative reverses this. It adds across lag distances and keeps the periods
apart, which gives one instrument per period. Addition and the average differ
here, because the number of available lags grows over time. So the average
across lag distances is a third option. A fourth option pushes the rule to
the end, where every column collapses into one instrument.

The average is also the one rule here that steps outside the sentence quoted
above. Its divisor is the count of lags that the country has in that period.
So the divisor varies from country to country inside a single column. The
result is not the original columns with some coefficients tied together, but
a new instrument. The new instrument is still legitimate, because the divisor
depends only on which lags exist and not on the outcome. The
equal-coefficient reading does not reach it, and the same caveat applies to
the fading family further down.

I ran all four collapse rules against the uncollapsed set. I then ran two
further families with a knob on them. The section after next covers those
families.

I also ran two instrument designs, because the design matters. The paper
builds its lag blocks from democracy only and gives income a single lagged
level. The symmetric alternative gives both variables the full block of lags.
The plm package uses that design by default, and so does the sweep in
docs/instruments.md. The first design reproduces the published GMM column
exactly. If it ever fails, the script stops.

One reading note comes before the numbers. The collapse of every column into
one leaves the model exactly identified here. That row therefore has no
overidentifying restrictions left, and output/aggregation.txt prints `exact
id` instead of a p-value. The estimate is still real, and it is the least
stable one on the figure. It is not comparable to the others on anything that
counts restrictions.

## What came out

The collapse rule matters, and by more than I expected. With every lag under
the paper's design, the estimates for Freedom House range from -0.243 to
-0.118. The estimates for Polity range from -0.553 to -0.351. Only the
collapse rule changes across those numbers. The data, the lags and the
estimator stay the same.

The top row of the figure shows the paper's own instrument design. There the
uncollapsed set is the exception. It travels much further than any of the
collapse rules. Its range is 0.084 for Freedom House and 0.245 for Polity,
against at most 0.044 for the four collapse rules. It also finishes closer to
the fixed-effects value than it started, by 0.039 and 0.190. That is the
overfitting story again.

The path is not a straight climb, and I do not want to say that it is. Both
uncollapsed lines first move further from the fixed-effects value, out to a
window of 4 for Freedom House and 5 for Polity. They reverse direction only
after that point. The uncollapsed line is monotone in one of the four panels
on the figure. That panel is Polity under the symmetric design. So the drift
is a net direction over the whole window, and not a steady march.

The bottom row gives income its own block of lags, and the result is messier.
That mess deserves a plain statement. There the collapse rules are not
uniformly steadier. Roodman's rule is the flattest line in three of the four
panels, with a range of at most 0.054. It loses only in the Polity panel
under the paper's design, where the exactly identified fully collapsed column
happens to sit still (0.010 against 0.021). That column is not much of a
rival, because it has nothing left to test.

Collapsing by period swings more than the uncollapsed set does (0.477 against
0.332 for Polity). Collapsing everything into one column is the least stable
line on the whole figure. So the useful statement is not that a collapse
always steadies the estimate. Collapsing by lag distance, the rule Roodman
proposes, is the one that stays steady under both designs.

Too much collapsing has its own failure mode. With the symmetric design,
collapsing everything into one column produces a standard error of 1.61 for
Freedom House. That number says the instrument no longer carries usable
information. A middle ground exists. Too little collapsing makes the
estimator overfit, and too much collapsing leaves too little information for
identification. That second failure is the weak instrument problem, and
docs/weak-instruments.md examines it.

## What the overidentification test says

Last round I looked only at stability. Stability on its own does not settle
much. A rule can sit still because it is well behaved, or because it no
longer listens to the data. Roodman's own prescription is to examine the
overidentification test as the instrument count falls. The table now carries
one test in every cell. The test is the sharper tool of the two, and it
points the same way.

Among the three collapse rules that leave anything to test, Roodman's rule
has the highest p-value in 27 of the 28 overidentified cells. The one
exception is Polity under the symmetric design at the widest window. There
the three rules cluster close to the threshold on both sides. The gap is not
close anywhere else. His p-values run from 0.031 to 0.998, and the test
rejects his rule at the 5% level in 1 of 28 cells. The test rejects the
collapse by period in 27 of 32 cells, and the collapse by period after
averaging in 32 of 32.

The other two rules are my own, and they turn Roodman's construction on its
side. Those two rules are not merely less steady than his. The data reject
them, and the data do not reject his rule.

The uncollapsed set needs the most careful comparison. Its p-values run from
0.00 to 0.39, and the test rejects it in 12 of 32 cells. For Freedom House
the p-value climbs steadily with the count, from 0.00 at the narrowest window
to 0.26 at the widest under the paper's design. For Polity under the same
design it stays low throughout, and it never rises above 0.13. Where the test
does look comfortable, that comfort is not proof of a valid instrument set.
Those cells have up to 99 instruments against about 127 countries, and that
is exactly where the test loses its power.

Roodman's rule reaches 0.998 on as few as 12 instruments. Against that mark,
a p-value of 0.39 on a set of 99 means much less.

## Turning the dial by degrees

The five rules above are all-or-nothing. Two more ways are not single rules
but families with a knob on them. With the knob I can ask what happens part
of the way. The results are in output/aggregation_families.txt and
output/aggregation_families.png.

The first family groups the years into blocks and collapses the columns
inside each block. A block size of one is then the uncollapsed set, and a
block as wide as the panel is Roodman's rule. The second family keeps one
instrument per year, but it multiplies each older lag by a fading factor
before the addition. A factor of one is then the plain sum. Both families
must match the fixed rules exactly at their endpoints. If they do not, the
script stops.

Neither family did what I expected, and they surprised me in opposite
directions.

Blocking is not a smooth dial at all. The estimates bounce and travel well
outside the two endpoints they connect. For Freedom House they run from
-0.274 to -0.105, while the endpoints themselves are only -0.129 and -0.154.

Even the instrument count refuses to fall for wider blocks. Block sizes of
four and five both give three blocks. But the wider blocks reach back to
years that have deeper lags on offer, so the count goes up rather than down.
Block size is therefore not a measure of the amount of collapsing, and that
is why the line looks like noise.

The fading family does almost nothing. Across the whole range of factors the
estimate moves by 0.032 for Freedom House and 0.015 for Polity. That result
is useful rather than disappointing. It says that the answer from the
collapse by period is not an artifact of equal weight on every lag. A tenth
of the weight on the older lags barely moves the estimate.

## Roodman's footnote 7: random subsets

The rules above all pick the instrument subsets deliberately. Roodman's
footnote 7 suggests a random pick instead. He asks how the coefficient and
the overidentification p-value move as the count grows. I drew 200 random
subsets of the uncollapsed lagged-level columns at each of several sizes,
under the paper's instrument design. I then refit the model.

The coefficient comes from the one-step estimator, so it is comparable with
the rest of this file. The Hansen test needs the two-step weight matrix, so
that column is two-step. The results are in output/aggregation_subsets.txt,
output/aggregation_subsets.csv and output/aggregation_subsets.png.

The first result shows how much of the answer comes from the analyst's choice
rather than from the data. With only 4 of the 45 available lagged levels in
play, the Freedom House draws run from -0.366 to 0.879. So on this data and
this specification, the sign of the effect is not fixed at all. The sign
depends on which instruments the draw contains. The spread then narrows as
the count grows. The standard deviation falls from 0.152 at 4 instruments to
0.017 at 44 for Freedom House, and from 0.076 to 0.029 for Polity.

The median moves with the spread, and it lands in the right place. At the
widest draws it is -0.133 for Freedom House and -0.354 for Polity, against
-0.129 and -0.351 for the full uncollapsed set. The script checks that. The
median reported standard error also falls from 0.144 to 0.077 for Freedom
House. So the estimate looks more precise as the count grows. It converges on
the uncollapsed answer that drifts toward fixed effects, which is the
conclusion the deliberate rules reach by a different route.

The J statistic is the messier half, and I do not want to read more into it
than it supports. The Hansen test rejects a share of the draws at the 5%
level. That share runs 0.41, 0.76, 0.69, 0.32, 0.00, 0.00 for Freedom House
across the six sizes, and 0.32, 0.36, 0.33, 0.28, 0.48, 0.94 for Polity.
Neither series is monotone, and the two do not agree with each other. So this
is not a clean demonstration of a loss of power in the test.

Part of the reason is itself the point of the paper. The two-step weight
matrix behind the test has one row and column per instrument. About 127
countries supply the estimate of that matrix. So by the widest draws the test
leans on a matrix that the data cannot support. The safe reading is the one
Roodman gives. A high Hansen p-value on a large instrument set is not
evidence of anything.

The diamonds on the figure put the deliberate rules on the same axes, which
is the comparison worth having. Roodman's rule uses 19 instruments and gives
-0.154 for Freedom House. Random draws of about that size (22 instruments)
run from -0.293 to 0.057. His rule therefore sits inside that range. The
range is wide enough that a position inside it is not much of a
recommendation on its own. For Polity the same comparison is -0.512 against a
range of -0.649 to -0.362.

The deliberate rule does not give a different answer at a given count. It
gives stability as the lag window widens. The figure earlier in this file
shows that stability, and the random draws cannot.

I did not implement the other suggestion in footnote 6, the
vector-autoregression restriction of Arellano (2003b). It needs a first-stage
model of the instruments themselves rather than a rule for the collapse of
the columns. So it does not fit into the same comparison. Roodman himself
notes that it did not enter common practice.

## A note on reading the table

The second column of the table uses every available lag, which is past the
right-hand edge of the figure. The figure stops at a window of 8. The
instrument counts in the table include only the lagged levels. They do not
include the year dummies that also sit in the instrument set.
