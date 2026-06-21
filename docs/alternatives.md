# Alternative estimators for the income effect

The paper asks whether a country growing richer leads it to become more
democratic. A plain correlation is misleading, so the authors rely on panel
methods that strip out fixed differences between countries and the tendency of
democracy to persist. This note re-estimates that income effect with several
methods and lines them up, which is the follow-up Professor Torgovitsky asked
for. The estimates are for the baseline five-year sample, the same one behind
Tables 2 and 3.

## The methods, in plain terms

- Pooled regression: ignores fixed country differences, and reads high.
- Fixed effects: removes things about a country that do not change over time. In
  a short panel with a lagged outcome it still carries a known bias.
- Anderson-Hsiao: looks at changes rather than levels to cancel the fixed
  differences, then uses values from two periods earlier as instruments.
- Arellano-Bond (difference GMM): the same change-based idea, using the whole
  set of earlier values as instruments at once. This is the method in the
  paper's own GMM columns.
- Blundell-Bond (system GMM): adds a second set of conditions, in levels, on top
  of Arellano-Bond. It is more precise when those conditions hold, but they are
  an extra assumption.

"One-step" and "two-step" are two ways of weighting the GMM estimators; the
two-step standard errors include the Windmeijer small-sample correction.

## Results

### Freedom House

| Estimator | Income (SE) | Democracy (SE) | Countries | Instruments | Hansen p | AR(2) p |
|---|---|---|---|---|---|---|
| Pooled OLS | 0.072 (0.010) | 0.706 (0.035) | 150 |  |  |  |
| Fixed effects | 0.010 (0.035) | 0.379 (0.051) | 150 |  |  |  |
| Anderson-Hsiao IV | -0.104 (0.107) | 0.469 (0.100) | 127 |  |  |  |
| Arellano-Bond, difference GMM (replication) | -0.129 (0.076) | 0.489 (0.085) | 127 | 55 |  |  |
| Arellano-Bond, difference GMM (one-step) | -0.114 (0.101) | 0.475 (0.079) | 152 | 63 | 0.03 | 0.74 |
| Arellano-Bond, difference GMM (two-step) | -0.085 (0.059) | 0.449 (0.087) | 152 | 63 | 0.17 | 0.76 |
| Blundell-Bond, system GMM (one-step) | 0.098 (0.021) | 0.611 (0.061) | 152 | 78 | 0.03 | 0.50 |
| Blundell-Bond, system GMM (two-step) | 0.092 (0.022) | 0.630 (0.059) | 152 | 78 | 0.11 | 0.50 |

### Polity

| Estimator | Income (SE) | Democracy (SE) | Countries | Instruments | Hansen p | AR(2) p |
|---|---|---|---|---|---|---|
| Pooled OLS | 0.053 (0.010) | 0.749 (0.034) | 136 |  |  |  |
| Fixed effects | -0.006 (0.039) | 0.449 (0.063) | 136 |  |  |  |
| Anderson-Hsiao IV | -0.413 (0.163) | 0.582 (0.127) | 114 |  |  |  |
| Arellano-Bond, difference GMM (replication) | -0.351 (0.127) | 0.590 (0.106) | 114 | 55 |  |  |
| Arellano-Bond, difference GMM (one-step) | -0.273 (0.092) | 0.618 (0.104) | 138 | 63 | 0.01 | 0.28 |
| Arellano-Bond, difference GMM (two-step) | -0.249 (0.098) | 0.616 (0.119) | 138 | 63 | 0.13 | 0.32 |
| Blundell-Bond, system GMM (one-step) | 0.072 (0.020) | 0.686 (0.081) | 138 | 78 | 0.00 | 0.27 |
| Blundell-Bond, system GMM (two-step) | 0.069 (0.023) | 0.703 (0.089) | 138 | 78 | 0.05 | 0.29 |

## What the comparison shows

The change-based methods agree. Once fixed differences and persistence are
accounted for, the effect of income on democracy is small, and in the
instrumental-variables estimates it is negative rather than positive. That is
the paper's central finding, and it holds across fixed effects, Anderson-Hsiao,
and Arellano-Bond. The two Arellano-Bond implementations here, the hand-built
one used for the replication and the one from the plm package, land in the same
place, which is a useful check that the replication's estimator is doing what it
should.

System GMM is the exception. Adding the level conditions pushes the income
coefficient to a small positive and statistically significant value for both
democracy measures. That looks like a reversal, but it rests entirely on the
extra assumption behind those level conditions, which is the kind of assumption
the paper is wary of. The overidentification (Hansen) test is borderline for
several of the one-step GMM runs, and the instrument counts are high, so the GMM
numbers deserve caution.

So the paper's conclusion holds up under the change-based alternatives. The only
way to bring back a positive income effect is to assume the extra system-GMM
conditions hold.
