# Are the instruments strong enough?

Professor Torgovitsky asked for weak instrument diagnostics, then for
Anderson-Rubin confidence sets and the Moreira conditional likelihood ratio
approach. All three apply to the two-stage least squares columns of Tables 5
and 6, which are the parts of the paper with one variable being instrumented,
income, and a small number of instruments for it: the savings rate in Table 5
and trade-weighted world income in Table 6. That gives ten columns to check.
The numbers are in output/weakiv.csv and output/weakiv.txt, and the p-value
curves are drawn in output/weakiv.png.

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

## Anderson-Rubin

The Anderson-Rubin test asks a simple question. Pick a candidate value for the
effect of income, subtract off what that value would imply, and see whether
anything is left over that the instruments can still explain. If something is,
that candidate value is rejected. Collecting every value that survives gives a
confidence set, and it stays valid no matter how weak the instruments are. I
clustered it by country to match the paper.

All of the 10 sets come out bounded, and all of them contain zero. Where the first
stage is strongest (Table 6 col 7, F of 39.5) the set is [-0.48, 0.06], which is
nearly the same as the ordinary interval of [-0.46, 0.06]. Where it is
weakest (Table 5 col 9, F of 12.1) the set widens to [-0.28, 0.26] against an
ordinary interval of [-0.17, 0.20]. That is the pattern you would hope to
see: the two agree when the instrument is strong and separate when it is not,
so the ordinary intervals understate the uncertainty in precisely the columns
where the first stage is thinnest.

## Moreira's conditional likelihood ratio

The conditional likelihood ratio test is the sharper of the two. Rather than
comparing against a fixed cutoff, it adjusts the cutoff using a quantity that
carries the information about how strong the instruments are, which buys back
power that Anderson-Rubin gives away. I coded it directly, including the
conditional cutoff by simulation, since no package for it was installed.

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
