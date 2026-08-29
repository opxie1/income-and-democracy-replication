# A simulation with known answers

Professor Torgovitsky asked for a Monte Carlo study calibrated on this data.
The idea is to build a dynamic panel that looks like the real one. It keeps
the same number of countries and the same number of periods, and it takes
persistence from the estimates here. In simulated data the true parameters
are known, so the error of each estimator is measurable. The results give a
guide to which estimates deserve trust in the real application.

The numbers are in output/montecarlo.txt and output/montecarlo.csv. Every
draw is in output/montecarlo_draws.csv.gz. The figures are
output/montecarlo_bias.png and output/montecarlo_coverage.png. The code is
R/16_montecarlo.R and R/_montecarlo.R. It runs 500 draws for each of the 8
designs.

## How the simulated data are built

Each country gets a democracy equation and an income equation. Democracy
depends on its own lag and on lagged income. Income depends on its own lag
and on lagged democracy. Both equations carry a country effect, a period
effect and a shock. The two shocks are drawn together, so a shock to one
variable can move the other. The country effects are drawn together as well,
with a correlation of 0.59 for Freedom House. That correlation is what makes
rich countries democratic for fixed reasons. It is also what pooled OLS
mistakes for an effect of income.

The simulated panel keeps the shape of the real one. It reproduces the
observation pattern cell by cell, for each variable separately, so the
estimators end up with the sample they have in the real data. For Freedom
House the estimation sample is 838 observations on 127 countries across 9
periods. For Polity it is 747 observations on 114 countries. The script
checks these against the real counts and stops if they differ.

The income equation is close to a random walk. Its persistence is 0.915 for
Freedom House and 0.958 for Polity, both from difference GMM. That number
matters more than it looks, and the results section returns to it.

The true effect of income on democracy is zero in every design. That choice
is deliberate. With a true zero, the share of draws in which an estimator
reports a significant effect is its false-positive rate. It answers the
question the paper asks, which is whether there is an effect at all.

## The two persistence values

Professor Torgovitsky asked for persistence set to what I estimate here. The
trouble is that the estimates disagree. For Freedom House, pooled OLS gives
0.697, fixed effects gives 0.364 and difference GMM gives 0.489. Picking one
of them presumes the answer that the simulation is meant to test.

So the study runs two calibrations. The first sets persistence to the
difference-GMM estimate, 0.489 for Freedom House and 0.590 for Polity. The
second sets it to the value that makes the simulated fixed-effects estimate
equal the real one, 0.600 and 0.702. The second calibration assumes nothing
about which estimator is right. It asks only what true persistence produces
the fixed-effects number that the data show.

## What the calibration matches, and what it misses

The variance of the country effects is chosen so that the spread of countries
within a period matches the real spread. For Freedom House the target
standard deviations are 0.348 for democracy and 1.025 for income.

One thing does not match, and it is worth stating. Pooled OLS is more biased
in the simulation than in the data. A linear model with a country effect and
a normal shock cannot match three things at once. Those three are the spread
of countries, the fixed-effects estimate and the pooled OLS estimate. The
simulation is therefore a slightly harsh world for the estimators that work
in levels. Democracy in the real data also sits between 0 and 1, and the
simulated series has no such bound.

## The two designs

Each calibration runs under two starting conditions. In the first, countries
start at their long-run mean, so the extra conditions that system GMM adds
are valid. In the second, the starting point is tied to the country effect,
so those extra conditions fail. Difference GMM stays valid under both. The
second design is the one Roodman simulates, and docs/alternatives.md finds
the same failure in the real data for Polity.

## What happens to the persistence estimate

The classic results appear, which is the first sign that the simulation is
built correctly. Under the difference-GMM calibration for Freedom House the
true persistence is 0.489. Pooled OLS averages 0.778 and fixed effects
averages 0.274. Their intervals cover the truth in 0.00 and 0.00 of draws.

The collapsed difference GMM estimator is almost unbiased. Its average error
is -0.027 one-step and -0.029 two-step. The price is noise, and the standard
deviation across draws is 0.117. The uncollapsed estimator is the opposite.
It is tighter, with a standard deviation of 0.071, but its average error is
-0.098. Its intervals cover the truth only 0.73 of the time. That is the cost
of too many instruments, measured directly.

## What happens to the income effect

The true effect is zero everywhere, so every number in this section is an
error. Pooled OLS reports a positive effect of 0.040 for Freedom House and
0.025 for Polity. It calls that effect significant in 1.00 and 0.99 of draws,
and its interval covers the true zero in 0.00 and 0.01. A researcher who
trusts pooled OLS on data of this shape finds an income effect that is not
there, nearly every time.

The estimator in the GMM columns of the paper is difference GMM with
uncollapsed instruments and every lag. Here it carries a negative bias of
-0.058 for Freedom House and -0.002 for Polity, on data where the true effect
is zero. The paper reports -0.129 and -0.351 in those columns. The bias
therefore covers 45 per cent of the published estimate for Freedom House and
1 per cent for Polity. That does not show that the published numbers are only
bias. It shows that a number of this kind arises on data of this shape when
the true effect is zero.

I do not have a clean account of why the two measures differ this much. My
first guess was the income process, which is close to a random walk and
therefore gives weak instruments. That guess does not survive. Income is more
persistent for Polity (0.958) than for Freedom House (0.915), so the weaker
instruments belong to the measure with the smaller bias. The honest statement
is that the bias is real and measure-specific, and that I cannot yet name its
source.

Collapsing helps a little on this coefficient. The collapsed bias is -0.063
against -0.091 uncollapsed for Freedom House. It costs a great deal of
precision, and the standard deviation across draws rises from 0.172 to 0.523.
The gain from collapsing is far clearer on the persistence coefficient than
on this one.

## Whether the confidence intervals are honest

Coverage answers the question directly. A 95 per cent interval must contain
the truth in 95 per cent of draws. Take the income effect, under the
difference-GMM calibration and a stationary start. For Freedom House the
estimator closest to that target is Blundell-Bond, system GMM (two-step) at
0.96. Pooled OLS is the worst at 0.00.

The starting condition decides whether system GMM can be trusted. Take the
two-step collapsed version. With a stationary start its intervals for the
income effect cover 0.96 for Freedom House and 0.96 for Polity. With the
start tied to the country effect the same intervals cover 0.45 and 0.79.
Difference GMM barely moves between the two designs. This is the warning that
the overidentification test gives in the real data, and here the cost of
ignoring it is measured.

## What this says about the real estimates

The answer depends on which coefficient is wanted, and that is the first
thing to say. On persistence, pooled OLS and fixed effects are both hopeless.
They miss by 0.289 and -0.215 for Freedom House, and neither interval ever
covers the truth. Collapsed difference GMM is the only estimator that is
close to unbiased and honest about its own uncertainty at the same time.

On the income coefficient the two part company. Pooled OLS still fails, with
coverage of 0.00. Fixed effects does not fail: its error is -0.036 and its
intervals cover the truth 0.90 of the time. The bias of the fixed-effects
estimator falls on the lagged dependent variable rather than on the other
regressor. So a reader who cares only about the income coefficient loses less
by using fixed effects than the persistence results suggest.

The uncollapsed difference GMM column, which is the one the paper reports, is
biased toward a negative income effect on data of this shape. Now take only
the estimators whose intervals cover the truth at least 90 per cent of the
time. Among those, the smallest root mean squared error on the income effect
for Freedom House belongs to Fixed effects at 0.058. Pooled OLS has a smaller
root mean squared error still, at 0.041, and that is the trap in using that
measure alone. A tight wrong answer beats a wide right one on root mean
squared error, and coverage is what separates them.

The practical reading matches the weak-instrument work in
docs/weak-instruments.md. The instrumental-variables columns of this paper
are not informative about the sign of the income effect. The simulation adds
a reason. On data with this many countries, this many periods and income this
persistent, no estimator here separates a true zero from the effect the paper
reports.

## Checks

Four checks run with the script. The estimator ladder is the same code that
R/10_alternatives.R runs on the real data, so the simulated and the real
columns cannot drift apart. The simulated estimation sample has to match the
real one, and it does, at 838 observations on 127 countries for Freedom
House. Every estimator has to return an estimate in at least 90 per cent of
draws. The bias has to shrink when the panel grows, and with 6 times as many
countries it falls for 11 of the 11 estimators. If a check fails, the script
stops.

That last check separates two kinds of error. The GMM errors fall by between
45 and 98 per cent as the panel grows, so they are finite-sample problems.
Pooled OLS and fixed effects fall by about 1 per cent, which is to say not at
all, because neither is consistent here at any sample size. That is the
difference between an estimator that needs a bigger panel and one that a
bigger panel cannot save.

One GMM row is slower than the rest, and it is worth naming. The estimator in
the GMM columns of the paper falls only 8 per cent, from -0.058 to -0.054.
Under the instrument design of the paper, income gets a single lagged level
as its instrument, so that estimator has the least to work with. I read this
as slow convergence rather than inconsistency, but the simulation here cannot
separate the two.

One limitation stands out. The true income effect is zero in every design
here, so these results measure bias and false positives. They do not measure
power against a real effect. That is the natural next step.

A second limitation is one of scope. The ladder covers every estimator in
docs/alternatives.md and the collapsed and uncollapsed variants from
docs/instruments.md. It does not cover the two-stage least squares columns of
Tables 5 and 6. Those columns need an instrument from outside the panel, the
savings rate and trade-weighted world income, and the simulated data contain
no such variable. Simulating one is a separate design and a separate
question.
