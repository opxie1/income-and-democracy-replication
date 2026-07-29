# References

The paper being replicated.

- Acemoglu, Daron, Simon Johnson, James A. Robinson, and Pierre Yared. 2008.
  "Income and Democracy." *American Economic Review* 98 (3): 808-842.
  Data from openICPSR project 113251.

On too many instruments and on collapsing them. This is the paper Professor
Torgovitsky sent, and Section V is the part on reducing the instrument count.

- Roodman, David. 2009. "A Note on the Theme of Too Many Instruments."
  *Oxford Bulletin of Economics and Statistics* 71 (1): 135-158.
- Roodman, David. 2009. "How to Do xtabond2: An Introduction to Difference and
  System GMM in Stata." *Stata Journal* 9 (1): 86-136. This is the companion
  piece describing the software whose behaviour R/00_setup.R reproduces.
- Arellano, Manuel. 2003. *Panel Data Econometrics*. Oxford University Press.
  Cited in a footnote of the first Roodman paper as the source of a
  vector-autoregression approach to constraining the instrument set.

The estimators being compared.

- Anderson, T. W., and Cheng Hsiao. 1981. "Estimation of Dynamic Models with
  Error Components." *Journal of the American Statistical Association* 76 (375):
  598-606.
- Arellano, Manuel, and Stephen Bond. 1991. "Some Tests of Specification for
  Panel Data." *Review of Economic Studies* 58 (2): 277-297. The employment
  results in this paper are the benchmark R/10_alternatives.R reproduces on
  every run.
- Blundell, Richard, and Stephen Bond. 1998. "Initial Conditions and Moment
  Restrictions in Dynamic Panel Data Models." *Journal of Econometrics* 87 (1):
  115-143.
- Holtz-Eakin, Douglas, Whitney Newey, and Harvey S. Rosen. 1988. "Estimating
  Vector Autoregressions with Panel Data." *Econometrica* 56 (6): 1371-1395.
  The source of the instrument matrix that collapsing operates on.
- Windmeijer, Frank. 2005. "A Finite Sample Correction for the Variance of
  Linear Efficient Two-Step GMM Estimators." *Journal of Econometrics* 126 (1):
  25-51.

Weak instruments.

- Anderson, T. W., and Herman Rubin. 1949. "Estimation of the Parameters of a
  Single Equation in a Complete System of Stochastic Equations." *Annals of
  Mathematical Statistics* 20 (1): 46-63.
- Moreira, Marcelo J. 2003. "A Conditional Likelihood Ratio Test for Structural
  Models." *Econometrica* 71 (4): 1027-1048.
- Andrews, Donald W. K., Marcelo J. Moreira, and James H. Stock. 2006.
  "Optimal Two-Sided Invariant Similar Tests for Instrumental Variables
  Regression." *Econometrica* 74 (3): 715-752.
- Stock, James H., and Motohiro Yogo. 2005. "Testing for Weak Instruments in
  Linear IV Regression." In *Identification and Inference for Econometric
  Models*, edited by Donald W. K. Andrews and James H. Stock, 80-108. Cambridge
  University Press. The cutoffs used here are the 10 percent maximal size
  values for one endogenous regressor: 16.38, 19.93, 22.30 and 24.58 for one
  through four instruments.
- Montiel Olea, José Luis, and Carolin Pflueger. 2013. "A Robust Test for Weak
  Instruments." *Journal of Business and Economic Statistics* 31 (3): 358-369.
  The effective F statistic here would be the properly justified cutoff for
  clustered standard errors, which is the caveat noted in
  docs/weak-instruments.md.

R packages doing the heavy lifting.

- Croissant, Yves, and Giovanni Millo. 2008. "Panel Data Econometrics in R: The
  plm Package." *Journal of Statistical Software* 27 (2).
- Fritsch, Markus, Joachim Schnurbus, and Andrew Adrian Yu Pua. 2021. "pdynmc:
  A Package for Estimating Linear Dynamic Panel Data Models Based on Nonlinear
  Moment Conditions." *R Journal* 13 (1): 218-231.
- Kang, Hyunseung, Yang Jiang, Qingyuan Zhao, and Dylan S. Small. 2021.
  "ivmodel: An R Package for Inference and Sensitivity Analysis of Instrumental
  Variables Models with One Endogenous Variable." *Observational Studies* 7 (2):
  1-24.
