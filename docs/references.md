# References

The paper that this project replicates.

- Acemoglu, Daron, Simon Johnson, James A. Robinson, and Pierre Yared. 2008.
  "Income and Democracy." *American Economic Review* 98 (3): 808-842.
  Data from openICPSR project 113251.

On too many instruments and on collapsing them. Professor Torgovitsky sent the
first paper in this group.

- Roodman, David. 2009. "A Note on the Theme of Too Many Instruments."
  *Oxford Bulletin of Economics and Statistics* 71 (1): 135-158. Section V,
  "Techniques for reducing the instrument count", is the relevant part. It
  gives lag limiting and collapsing as the two techniques in use. It crosses
  them in its Table 1. It adds two more ideas in footnotes 6 and 7.
- Roodman, David. 2009. "How to Do xtabond2: An Introduction to Difference and
  System GMM in Stata." *Stata Journal* 9 (1): 86-136. This is the companion
  piece. It describes the software, and R/00_setup.R reproduces the behavior of
  that software.
- Arellano, Manuel. 2003. "Modeling Optimal Instrumental Variables for Dynamic
  Panel Data Models." Working Paper 0310, Centro de Estudios Monetarios y
  Financieros, Madrid. This is Roodman's "Arellano (2003b)". He cites it twice.
  In Section V he cites it for the reading of lag limiting as a zero
  restriction on the projection coefficients. In footnote 6 he cites it as the
  source of the vector-autoregression approach that constrains the instrument
  set.
- Arellano, Manuel. 2003. *Panel Data Econometrics*. Oxford University Press.
  This is Roodman's "Arellano (2003a)". He cites it once, at page 171. He cites
  it only for the textbook remark that IV estimators perform poorly with many
  instruments. This entry separates the two, because both are from 2003.

The estimators that this project compares.

- Anderson, T. W., and Cheng Hsiao. 1981. "Estimation of Dynamic Models with
  Error Components." *Journal of the American Statistical Association* 76 (375):
  598-606.
- Arellano, Manuel, and Stephen Bond. 1991. "Some Tests of Specification for
  Panel Data." *Review of Economic Studies* 58 (2): 277-297. The employment
  results in this paper are the benchmark. R/10_alternatives.R reproduces them
  on every run.
- Blundell, Richard, and Stephen Bond. 1998. "Initial Conditions and Moment
  Restrictions in Dynamic Panel Data Models." *Journal of Econometrics* 87 (1):
  115-143.
- Holtz-Eakin, Douglas, Whitney Newey, and Harvey S. Rosen. 1988. "Estimating
  Vector Autoregressions with Panel Data." *Econometrica* 56 (6): 1371-1395.
  This paper is the source of the instrument matrix that collapsing operates
  on.
- Nickell, Stephen. 1981. "Biases in Dynamic Models with Fixed Effects."
  *Econometrica* 49 (6): 1417-1426. The source of the downward bias in the
  fixed-effects estimate of persistence, which the simulation in
  docs/monte-carlo.md measures directly.
- Windmeijer, Frank. 2005. "A Finite Sample Correction for the Variance of
  Linear Efficient Two-Step GMM Estimators." *Journal of Econometrics* 126 (1):
  25-51.

Weak instruments.

- Anderson, T. W., and Herman Rubin. 1949. "Estimation of the Parameters of a
  Single Equation in a Complete System of Stochastic Equations." *Annals of
  Mathematical Statistics* 20 (1): 46-63.
- Moreira, Marcelo J. 2003. "A Conditional Likelihood Ratio Test for Structural
  Models." *Econometrica* 71 (4): 1027-1048.
- Kleibergen, Frank. 2005. "Testing Parameters in GMM Without Assuming That
  They Are Identified." *Econometrica* 73 (4): 1103-1123. This is the source of
  the cluster-robust version of the conditional likelihood ratio test in
  R/_weakiv.R. I built that version from his K and J statistics. It uses a rank
  statistic in place of Moreira's conditioning statistic.
- Cragg, John G., and Stephen G. Donald. 1993. "Testing Identifiability and
  Specification in Instrumental Variable Models." *Econometric Theory* 9 (2):
  222-240. This is the source of the minimum-eigenvalue statistic. I report
  that statistic for the two columns with two endogenous regressors.
- Andrews, Donald W. K., Marcelo J. Moreira, and James H. Stock. 2006.
  "Optimal Two-Sided Invariant Similar Tests for Instrumental Variables
  Regression." *Econometrica* 74 (3): 715-752.
- Stock, James H., and Motohiro Yogo. 2005. "Testing for Weak Instruments in
  Linear IV Regression." In *Identification and Inference for Econometric
  Models*, edited by Donald W. K. Andrews and James H. Stock, 80-108. Cambridge
  University Press. The cutoffs used here are the 10 percent maximal size
  values for one endogenous regressor. They are 16.38, 19.93, 22.30 and 24.58
  for one through four instruments.
- Montiel Olea, José Luis, and Carolin Pflueger. 2013. "A Robust Test for Weak
  Instruments." *Journal of Business and Economic Statistics* 31 (3): 358-369.
  The effective F statistic is the properly justified first-stage measure for
  clustered errors. I report it next to the ordinary one. The conservative
  critical value used here is 23.109. It is their worst case over covariance
  structures for two-stage least squares, at a 10 percent Nagar-bias threshold
  and a 5 percent test.

The R packages that do the main work.

- Croissant, Yves, and Giovanni Millo. 2008. "Panel Data Econometrics in R: The
  plm Package." *Journal of Statistical Software* 27 (2).
- Fritsch, Markus, Joachim Schnurbus, and Andrew Adrian Yu Pua. 2021. "pdynmc:
  A Package for Estimating Linear Dynamic Panel Data Models Based on Nonlinear
  Moment Conditions." *R Journal* 13 (1): 218-231.
- Kang, Hyunseung, Yang Jiang, Qingyuan Zhao, and Dylan S. Small. 2021.
  "ivmodel: An R Package for Inference and Sensitivity Analysis of Instrumental
  Variables Models with One Endogenous Variable." *Observational Studies* 7 (2):
  1-24.
