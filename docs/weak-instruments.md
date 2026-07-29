# Are the instruments strong enough?

Professor Torgovitsky asked for weak instrument diagnostics, then for
Anderson-Rubin confidence sets and the Moreira conditional likelihood ratio
approach. All three apply to the two-stage least squares columns of Tables 5
and 6, which are the parts of the paper with one variable being instrumented,
income, and a small number of instruments for it: the savings rate in Table 5
and trade-weighted world income in Table 6. That gives ten columns to check.
The numbers are in output/weakiv.csv and output/weakiv.txt, and the p-value
curves are drawn in output/weakiv.png. The two tests are from Anderson and
Rubin (1949) and Moreira (2003), the first-stage cutoffs from Stock and Yogo
(2005), and full details are in docs/references.md.

## Why this matters

An instrument is weak when it barely moves the variable it is supposed to be
standing in for. When that happens the usual standard errors and confidence
intervals are not trustworthy, because they are built on the assumption that
the instrument has real pull. The tests below are built to stay honest whether
the instrument is strong or weak, so comparing them against the ordinary
interval shows how much the ordinary one was relying on that assumption.

## How strong is the first stage

The first-stage F statistics run from 12.1 to 39.5, clustered by country the
same way the paper clusters everything else. Every column clears the old rule
of thumb that an F above 10 is good enough. Against the stricter Stock and
Yogo threshold, 4 of the 10 columns fall short. So the instruments are not
badly weak, but a few are not comfortably strong either, which is exactly the
situation the next two tests are designed for.

One caveat on that comparison. The Stock and Yogo cutoffs were worked out for a
setting with well behaved, unclustered errors, and the F statistics here are
clustered by country. Lining the two up is what most applied work does, and it
is what the standard software prints, but it is not strictly justified. I use
the threshold as a rough marker rather than a real test, which is another
reason to lean on the confidence sets below instead. The properly justified
cutoff for clustered errors would be the effective F statistic of Montiel Olea
and Pflueger (2013), which I have not implemented.

## Anderson-Rubin

The Anderson-Rubin test asks a simple question. Pick a candidate value for the
effect of income, subtract off what that value would imply, and see whether
anything is left over that the instruments can still explain. If something is,
that candidate value is rejected. Collecting every value that survives gives a
confidence set, and it stays valid no matter how weak the instruments are. I
clustered it by country to match the paper.

All 10 sets come out bounded, and all of them contain zero. Where the first stage
is strongest (Table 6 col 7, F of 39.5) the set is [-0.48, 0.06], which is
nearly the same as the ordinary interval of [-0.46, 0.06]. Where it is
weakest (Table 5 col 9, F of 12.1) the set widens to [-0.28, 0.26] against an
ordinary interval of [-0.17, 0.20]. The Anderson-Rubin set is at least as wide
as the ordinary interval in every column, so the ordinary intervals are
optimistic throughout, not just at the extremes.

I should not push that comparison further than it goes. Lining the ten columns
up by first-stage F does not give a tidy ordering. The widening ratio runs
from 1.05 to 1.45, and it does not fall neatly as the first stage gets
stronger. Table 6 col 8 has the weaker first stage of the two,
F of 14.2 against 17.8, yet widens less, 1.21 against 1.31, than Table 6 col 9 does.
These specifications
differ in their controls and samples as well as in instrument strength, so the
two ends of the range are suggestive rather than a clean relationship.

## Moreira's conditional likelihood ratio

The conditional likelihood ratio test can be sharper than Anderson-Rubin.
Rather than comparing against a fixed cutoff, it adjusts the cutoff using a
quantity that carries the information about how strong the instruments are,
which buys back power that Anderson-Rubin gives away. I coded it directly,
since no package for it was installed.

The catch is that 8 of the 10 columns here have exactly one instrument, and
with one instrument there is nothing to buy back: the two tests are the same
test. So the extra power is only in play for the two columns with a pair of
instruments. For the single-instrument columns the cutoff is a chi-squared
quantile that can be written down exactly, and the code uses that rather than
simulating a distribution it already knows.

There is a caveat worth stating plainly. The version of this test I implemented
assumes the errors are well behaved and independent, and it does not cluster by
country. The paper clusters, and the Anderson-Rubin sets above cluster. So the
conditional likelihood ratio sets come out narrower than the Anderson-Rubin
ones, but that gap is mostly the clustering, not the extra power of the test.
The code checks this: for the columns with a single instrument the conditional
likelihood ratio test and the unclustered Anderson-Rubin test agree, and the
run stops if they ever do not. The right comparison is therefore
Anderson-Rubin, and I treat the other as a benchmark.

## What it adds up to

Under inference that does not lean on the instruments being strong, every one
of the 10 two-stage least squares columns leaves zero inside the confidence
set. So none of them establishes an effect of income on democracy in either
direction.

This does not overturn anything in the paper, whose argument is that the effect
is not there once the fixed differences between countries are taken out. It
does mean the instrumental-variables columns should be read as uninformative
about the sign of the effect rather than as evidence for a negative one.

## Checks

This is fiddly enough that I did not want to rely on code I wrote myself. Three
checks run every time the script does. The partialling step has to reproduce
the published two-stage least squares estimates, the fast Anderson-Rubin curve
has to match the slower clustered regression it stands in for, and the
conditional likelihood ratio test has to collapse onto Anderson-Rubin for the
columns with a single instrument. Any failure stops the run.

On top of that, R/15_ivcrosscheck.R sends all ten columns through the ivmodel
package and stops unless its conditional likelihood ratio sets agree with mine
and its two-stage least squares estimates match exactly. They do, for all ten.
That script needs the ivmodel package, so it sits outside the main run; use
Rscript R/15_ivcrosscheck.R.
