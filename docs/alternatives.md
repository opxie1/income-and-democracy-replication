# Alternative estimators for the income effect

The paper asks whether a country becomes more democratic when it becomes richer.
A plain correlation can give the wrong answer. The authors therefore use methods
that remove the steady differences between countries and account for the slow
movement of democracy. Here I estimate the income effect in several ways and put
the results side by side. Professor Torgovitsky asked for this comparison. All of
the estimates use the five-year sample, the same sample behind Tables 2 and 3.

## The methods, in plain terms

- Pooled OLS pools the observations across countries. This method pools
  observations, not instruments. Pooled OLS ignores the steady differences
  between countries, and the estimate is therefore too high.
- Fixed effects removes anything about a country that stays the same over time.
  With few time periods and a lagged outcome, this method still has a known
  bias.
- Anderson-Hsiao uses changes instead of levels to cancel the steady
  differences. It then uses values from two periods earlier as instruments.
- Arellano-Bond (difference GMM) also works with changes, and it uses a set of
  earlier values as instruments. The paper uses this method in its GMM columns.
- Blundell-Bond (system GMM) keeps the Arellano-Bond conditions and adds a
  second set of conditions in levels. When these level conditions hold, system
  GMM is more precise. However, the level conditions are an extra assumption.

"One-step" and "two-step" are two weighting schemes for these GMM estimators.
The two-step standard errors use the Windmeijer correction. The GMM rows use a
small instrument set: collapsed, and limited to lags two through four. A small
instrument set keeps the diagnostic tests informative.

## How to read the diagnostics

- A valid specification gives a small AR(1) p and a large AR(2) p. Every GMM row
  here shows this pattern.
- A valid specification gives an overidentification p that is not small. A small
  value warns that some of the instruments can be invalid. Each GMM row reports
  the test of its own estimator. The one-step rows report the Sargan test, which
  assumes well-behaved errors. The two-step rows report the Hansen test, which
  does not make this assumption.
- The Sargan test rejects every GMM row here, so it separates nothing. The paper
  clusters all of its own inference by country, and that is the situation in
  which a homoskedastic test over-rejects. The Hansen column is the one to read.

## Results

### Freedom House

| Estimator | Income (SE) | Democracy (SE) | Countries | Instruments | AR(1) p | AR(2) p | Overid p |
|---|---|---|---|---|---|---|---|
| Pooled OLS | 0.072 (0.010) | 0.706 (0.035) | 150 |  |  |  |  |
| Fixed effects | 0.010 (0.035) | 0.379 (0.051) | 150 |  |  |  |  |
| Anderson-Hsiao IV | -0.104 (0.107) | 0.469 (0.100) | 127 |  |  |  |  |
| Arellano-Bond, difference GMM (replication) | -0.129 (0.076) | 0.489 (0.085) | 127 | 55 |  |  |  |
| Arellano-Bond, difference GMM (one-step) | -0.189 (0.133) | 0.484 (0.091) | 124 | 13 | 0.00 | 0.75 | 0.00 |
| Arellano-Bond, difference GMM (two-step) | -0.133 (0.127) | 0.514 (0.094) | 124 | 13 | 0.00 | 0.69 | 0.37 |
| Blundell-Bond, system GMM (one-step) | 0.099 (0.024) | 0.583 (0.068) | 133 | 16 | 0.00 | 0.51 | 0.00 |
| Blundell-Bond, system GMM (two-step) | 0.100 (0.024) | 0.599 (0.062) | 133 | 16 | 0.00 | 0.51 | 0.30 |

### Polity

| Estimator | Income (SE) | Democracy (SE) | Countries | Instruments | AR(1) p | AR(2) p | Overid p |
|---|---|---|---|---|---|---|---|
| Pooled OLS | 0.053 (0.010) | 0.749 (0.034) | 136 |  |  |  |  |
| Fixed effects | -0.006 (0.039) | 0.449 (0.063) | 136 |  |  |  |  |
| Anderson-Hsiao IV | -0.413 (0.163) | 0.582 (0.127) | 114 |  |  |  |  |
| Arellano-Bond, difference GMM (replication) | -0.351 (0.127) | 0.590 (0.106) | 114 | 55 |  |  |  |
| Arellano-Bond, difference GMM (one-step) | -0.472 (0.213) | 0.627 (0.132) | 111 | 13 | 0.00 | 0.26 | 0.00 |
| Arellano-Bond, difference GMM (two-step) | -0.480 (0.211) | 0.646 (0.137) | 111 | 13 | 0.00 | 0.28 | 0.51 |
| Blundell-Bond, system GMM (one-step) | 0.073 (0.022) | 0.703 (0.088) | 120 | 16 | 0.00 | 0.26 | 0.00 |
| Blundell-Bond, system GMM (two-step) | 0.055 (0.023) | 0.834 (0.075) | 120 | 16 | 0.00 | 0.26 | 0.01 |

## What the comparison shows

The change-based methods agree with each other. These methods remove the steady
differences between countries and account for the slow movement of democracy.
After these corrections, the effect of income on democracy is small. In the
instrumental-variables estimates the effect is negative, not positive. This
result is the main finding of the paper, and fixed effects, Anderson-Hsiao and
Arellano-Bond all show it.

The Arellano-Bond rows differ among themselves by less than one standard error.
I do not expect closer agreement than this. The replication row uses the
instrument set of the paper, which is uncollapsed and uses every lag. The plm
rows use a collapsed set limited to lags two through four.

The gap between the one-step and two-step rows carries little information here.
The gap is small because the instrument set is small. The gap is not small in
general. With a wide lag window, the weighting step decides whether collapsing
works at all. The file docs/instruments.md covers this question.

System GMM is the exception. The level conditions push the income coefficient to
a small positive value for both democracy measures, and that value is
statistically significant. This reversal rests entirely on the level conditions.
The two democracy measures treat these conditions differently. For Polity, the
Hansen test rejects them outright. For Freedom House, the Hansen test raises no
objection.

For Freedom House the positive estimate therefore depends on one assumption. The
paper argues that democracies in transition are unlikely to satisfy this
assumption. Either way, every change-based estimate is small or negative. The
positive and significant number appears only in the rows that add the level
conditions.

The conclusion of the paper therefore holds under the change-based methods. A
positive income effect returns only under the assumption that the level
conditions hold.

## Checks

I checked these numbers two ways. The GMM engine is the pgmm function in plm. It
reproduces the textbook Arellano-Bond (1991) employment results exactly. If it
does not reproduce them, the script stops.

The script R/11_crosscheck.R runs difference and system GMM through an
independent package, pdynmc, with its own uncollapsed instrument set. That
script must find the same result as the tables above: income negative under
difference GMM and positive under system GMM. If it finds anything else, the
script stops. R/11_crosscheck.R needs the pdynmc package.

To run the cross-check, use `Rscript R/11_crosscheck.R`.
