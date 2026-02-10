# fEGarch 1.0.6 (2026-02-10)
- GJR-GARCH, TGARCH, APARCH, FIGJR-GARCH, FITGARCH and FIAPARCH 
  pre-sample values have been stabilized slightly.
- In dual models, where for some models in the volatility part pre-sample
  quantities are set through empirical values, a preliminary
  ARMA or FARIMA estimation step has been implemented to obtain their
  residuals and subsequently obtain the pre-sample values based on those
  residuals; this stabilizes estimated conditional standard deviations
  at the beginning of the training period; dual models with
  long-memory EGF part for the variance model stay fully unaffected.
- fixed a bug in `garch_sim`, where the simulation failed, if no
  setting of `phi` in the list for argument `pars` was provided
  manually.

# fEGarch 1.0.5 (2026-02-02)
- a bug in the backtesting functions for output of 
  `measure_risk,fEGarch_distr_est-method` was fixed that was introduced
  with package version 1.0.4; this did not affect the `measure_risk` function.
- step size different from 1 in rolling volatility forecasts under model
  refitting now available for GARCH-type models different from the exponential
  GARCH family (was already and is available for exponential GARCH family
  models).

# fEGarch 1.0.4 (2026-01-07)
- the method `predict_roll()` now allows for refitting of the underlying
  model via the new argument `refit_after` during the test period.
- `measure_risk()` and backtest functions take into account the potential
  refitting to compute risk measures and backtest quantities.
- all fitting functions now include an argument `skip_vcov` to allow for
  skipping the variance-covariance matrix computation

# fEGarch 1.0.3 (2025-11-07)
- numerical stability of estimators improved regarding parameter mu in cases
  with strong outliers
- some parameter default starting values were adjusted slightly to further
  increase numerical stability
- typos in the documentation for `trafflight_test,fEGarch_risk-method` were
  fixed
- some unnecessary elements in the manual for internal functions have been
  removed from the manual for improved clarity

# fEGarch 1.0.2 (2025-09-11)
- required package update for RcppArmadillo update to Armadillo 15.0.*
- a bug was fixed, where for the nonparametric scale estimation always the
  automated bandwidth selection under long-memory errors was employed,
  even for short-memory models; now the functions select the appropriate
  bandwidth selection algorithm correctly for short- and long-memory
  models
- a safety measure was included to avoid that the log-likelihood can
  become -Inf
- math formulas were fully removed from the README, because they were not
  displayed properly on CRAN
- some minor typos were corrected in the README

# fEGarch 1.0.1 (2025-06-20)
- a bug was fixed, where for most models except for EGARCH-type models the
  mean estimate was not considered properly in volatility forecasts
- a bug was fixed, where for EGARCH-family models rolling point forecasts
  under the (skewed) ALD were erroneous
- the mathematical formula in the README was not displayed correctly on the
  CRAN servers and was therefore removed entirely
- the package description was adjusted to contain also information on the
  the more advanced models and backtesting capabilities of the package,
  making it easier for users to see what to expect from the package
- the package title was adjusted to contain more information on the
  package's contents, making it easier for users to see what to expect from
  the package
- a bug was fixed in the computation of the hessian matrix, where now
  it is caught if a hessian is not invertible through solve() and
  the fitting functions therefore now don't stop due to an error anymore
- an additional dataset `UKinflation` was added to the package to show
  applications of the package beyond return data
- some examples in the README were adjusted / added
