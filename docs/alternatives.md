# Alternative estimators for the income effect

The paper asks whether a country growing richer ends up more democratic. A plain
correlation can fool you, so the authors use methods that take out the steady
differences between countries and the fact that democracy changes slowly. Here I
estimate that income effect several different ways and put them side by side.
This is the follow-up Professor Torgovitsky asked for. All of these use the
five-year sample, the same one behind Tables 2 and 3.

## The methods, in plain terms

- Pooled regression: ignores the steady differences between countries, so it
  reads too high.
- Fixed effects: takes out anything about a country that stays the same over
  time. With only a few time periods and a lagged outcome, it still has a known
  bias.
- Anderson-Hsiao: looks at changes instead of levels to cancel the steady
  differences, then uses values from two periods earlier as instruments.
- Arellano-Bond (difference GMM): the same change-based idea, using a set of
  earlier values as instruments. This is the method in the paper's GMM columns.
- Blundell-Bond (system GMM): adds a second set of conditions, in levels, on top
  of Arellano-Bond. It is more precise when those conditions hold, but they are
  an extra assumption.

"One-step" and "two-step" are two ways of weighting these GMM estimators. The
two-step standard errors use the Windmeijer correction. To keep the number of
instruments small, so the tests below mean something, the GMM rows use a small
instrument set: collapsed, and limited to lags two through four.

## How to read the diagnostics

- The AR(1) p should be small and the AR(2) p should be large. That is the
  pattern you want, and it holds throughout.
- The overidentification p (Hansen's J test) should not be small. A small
  value warns that the instruments may not all be valid.

## Results

### Freedom House

| Estimator | Income (SE) | Democracy (SE) | Countries | Instruments | AR(1) p | AR(2) p | Overid p |
|---|---|---|---|---|---|---|---|
| Pooled OLS | 0.072 (0.010) | 0.706 (0.035) | 150 |  |  |  |  |
| Fixed effects | 0.010 (0.035) | 0.379 (0.051) | 150 |  |  |  |  |
| Anderson-Hsiao IV | -0.104 (0.107) | 0.469 (0.100) | 127 |  |  |  |  |
| Arellano-Bond, difference GMM (replication) | -0.129 (0.076) | 0.489 (0.085) | 127 | 55 |  |  |  |
| Arellano-Bond, difference GMM (one-step) | -0.189 (0.133) | 0.484 (0.091) | 124 | 13 | 0.00 | 0.75 | 0.21 |
| Arellano-Bond, difference GMM (two-step) | -0.133 (0.127) | 0.514 (0.094) | 124 | 13 | 0.00 | 0.69 | 0.37 |
| Blundell-Bond, system GMM (one-step) | 0.099 (0.024) | 0.583 (0.068) | 133 | 16 | 0.00 | 0.51 | 0.18 |
| Blundell-Bond, system GMM (two-step) | 0.100 (0.024) | 0.599 (0.062) | 133 | 16 | 0.00 | 0.51 | 0.30 |

### Polity

| Estimator | Income (SE) | Democracy (SE) | Countries | Instruments | AR(1) p | AR(2) p | Overid p |
|---|---|---|---|---|---|---|---|
| Pooled OLS | 0.053 (0.010) | 0.749 (0.034) | 136 |  |  |  |  |
| Fixed effects | -0.006 (0.039) | 0.449 (0.063) | 136 |  |  |  |  |
| Anderson-Hsiao IV | -0.413 (0.163) | 0.582 (0.127) | 114 |  |  |  |  |
| Arellano-Bond, difference GMM (replication) | -0.351 (0.127) | 0.590 (0.106) | 114 | 55 |  |  |  |
| Arellano-Bond, difference GMM (one-step) | -0.472 (0.213) | 0.627 (0.132) | 111 | 13 | 0.00 | 0.26 | 0.41 |
| Arellano-Bond, difference GMM (two-step) | -0.480 (0.211) | 0.646 (0.137) | 111 | 13 | 0.00 | 0.28 | 0.51 |
| Blundell-Bond, system GMM (one-step) | 0.073 (0.022) | 0.703 (0.088) | 120 | 16 | 0.00 | 0.26 | 0.00 |
| Blundell-Bond, system GMM (two-step) | 0.055 (0.023) | 0.834 (0.075) | 120 | 16 | 0.00 | 0.26 | 0.01 |

## What the comparison shows

The change-based methods agree. Once the steady differences and the slow movement
of democracy are accounted for, the effect of income on democracy is small, and
in the instrumental-variables estimates it is negative rather than positive. This
is the paper's main finding, and it shows up in fixed effects, Anderson-Hsiao,
and Arellano-Bond alike. The two Arellano-Bond versions here, the one I wrote by
hand for the replication and the one from the plm package, land in the same
place.

System GMM is the exception. Adding the level conditions pushes the income
coefficient to a small positive, statistically significant value for both
democracy measures. That reversal leans entirely on the extra conditions, and
the two measures treat them differently. For Polity the overidentification
test rejects them outright. For Freedom House the test raises no objection, so
there the positive estimate stands or falls with an assumption the paper
argues democracies in transition are unlikely to satisfy. Either way, the
estimates built on changes alone are all small or negative; the positive,
significant number appears only once the added conditions come in.

So the paper's conclusion holds up under the change-based methods. The only way
to get a positive income effect back is to assume the extra system-GMM conditions
hold.

## Checks

I checked these numbers two ways. The GMM engine (plm's pgmm) reproduces the
textbook Arellano-Bond (1991) employment results exactly, and the run stops if
it ever fails to. Separately, R/11_crosscheck.R runs difference and system GMM
through an independent package, pdynmc, with its own uncollapsed instrument
set, and it stops with an error unless it finds the same picture as the tables
above: income negative under difference GMM and positive under system GMM.
That script needs the pdynmc package; run it with Rscript R/11_crosscheck.R.
