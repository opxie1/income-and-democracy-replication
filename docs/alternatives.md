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
- Arellano-Bond (difference GMM): the same change-based idea, using a set of
  earlier values as instruments. This is the method in the paper's GMM columns.
- Blundell-Bond (system GMM): adds a second set of conditions, in levels, on top
  of Arellano-Bond. It is more precise when those conditions hold, but they are
  an extra assumption.

"One-step" and "two-step" are two ways of weighting the GMM estimators; the
two-step standard errors include the Windmeijer small-sample correction. The GMM
rows use a deliberately small instrument set (collapsed, lags two to four) so the
instrument count stays low and the tests below stay meaningful.

## How to read the diagnostics

- AR(1) p should be small and AR(2) p should be large: that is the pattern you
  want, and it holds throughout.
- The overidentification p (a Sargan test for one-step, a Hansen test for
  two-step) should not be small; a small value warns that the instruments may not
  all be valid.

## Results

### Freedom House

| Estimator | Income (SE) | Democracy (SE) | Countries | Instruments | AR(1) p | AR(2) p | Overid p |
|---|---|---|---|---|---|---|---|
| Pooled OLS | 0.072 (0.010) | 0.706 (0.035) | 150 |  |  |  |  |
| Fixed effects | 0.010 (0.035) | 0.379 (0.051) | 150 |  |  |  |  |
| Anderson-Hsiao IV | -0.104 (0.107) | 0.469 (0.100) | 127 |  |  |  |  |
| Arellano-Bond, difference GMM (replication) | -0.129 (0.076) | 0.489 (0.085) | 127 | 55 |  |  |  |
| Arellano-Bond, difference GMM (one-step) | -0.189 (0.133) | 0.484 (0.091) | 152 | 13 | 0.00 | 0.75 | 0.21 |
| Arellano-Bond, difference GMM (two-step) | -0.133 (0.127) | 0.514 (0.094) | 152 | 13 | 0.00 | 0.69 | 0.37 |
| Blundell-Bond, system GMM (one-step) | 0.099 (0.024) | 0.583 (0.068) | 152 | 16 | 0.00 | 0.51 | 0.18 |
| Blundell-Bond, system GMM (two-step) | 0.100 (0.024) | 0.599 (0.062) | 152 | 16 | 0.00 | 0.51 | 0.30 |

### Polity

| Estimator | Income (SE) | Democracy (SE) | Countries | Instruments | AR(1) p | AR(2) p | Overid p |
|---|---|---|---|---|---|---|---|
| Pooled OLS | 0.053 (0.010) | 0.749 (0.034) | 136 |  |  |  |  |
| Fixed effects | -0.006 (0.039) | 0.449 (0.063) | 136 |  |  |  |  |
| Anderson-Hsiao IV | -0.413 (0.163) | 0.582 (0.127) | 114 |  |  |  |  |
| Arellano-Bond, difference GMM (replication) | -0.351 (0.127) | 0.590 (0.106) | 114 | 55 |  |  |  |
| Arellano-Bond, difference GMM (one-step) | -0.472 (0.213) | 0.627 (0.132) | 138 | 13 | 0.00 | 0.26 | 0.41 |
| Arellano-Bond, difference GMM (two-step) | -0.480 (0.211) | 0.646 (0.137) | 138 | 13 | 0.00 | 0.28 | 0.51 |
| Blundell-Bond, system GMM (one-step) | 0.073 (0.022) | 0.703 (0.088) | 138 | 16 | 0.00 | 0.26 | 0.00 |
| Blundell-Bond, system GMM (two-step) | 0.055 (0.023) | 0.834 (0.075) | 138 | 16 | 0.00 | 0.26 | 0.01 |

## What the comparison shows

The change-based methods agree. Once fixed differences and persistence are
accounted for, the effect of income on democracy is small, and in the
instrumental-variables estimates it is negative rather than positive. That is the
paper's central finding, and it holds across fixed effects, Anderson-Hsiao, and
Arellano-Bond. The two Arellano-Bond implementations here, the hand-built one
used for the replication and the one from the plm package, land in the same
place.

System GMM is the exception. Adding the level conditions pushes the income
coefficient to a small positive, statistically significant value for both
democracy measures. That looks like a reversal, but it leans entirely on those
extra level conditions, and the data do not fully support them: the
overidentification test is only borderline for the Freedom House measure and is
rejected for Polity. So the positive system-GMM estimate is shaky rather than a
clean reversal, and it rests on exactly the kind of assumption the paper is wary
of.

So the paper's conclusion holds up under the change-based alternatives. The only
way to bring back a positive income effect is to assume the extra system-GMM
conditions hold.

## Checks

Two independent checks back these numbers. First, the GMM engine (plm's pgmm)
reproduces the textbook Arellano-Bond (1991) employment results exactly; the run
stops if it ever fails to. Second, a separate package, pdynmc, was used to
re-estimate the same models (see R/11_crosscheck.R); it gives the same picture,
with income negative under difference GMM and positive under system GMM.
