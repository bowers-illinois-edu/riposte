# Covariance-adjustment specification

Bundles the covariates, the learner, and whether to refit per
permutation. Pass the result as `adjust =` to
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md),
or pass a bare one-sided covariate formula and let
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
wrap it with the default learner.

## Usage

``` r
riposte_adjust(covariates, learner = riposte_ridge_learner(), refit = TRUE)
```

## Arguments

- covariates:

  a one-sided formula of covariates, e.g. `~ x1 + x2`.

- learner:

  a learner `(X, Y, control) -> Yhat`; ridge by default. See
  [`riposte_ridge_learner()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_ridge_learner.md).

- refit:

  refit the learner inside every re-randomization (`TRUE`, exact for any
  learner) or fit once on the observed controls and permute the fixed
  residuals (`FALSE`, the shortcut that breaks exactness when the
  learner overfits — kept for comparison).

## Value

a `riposte_adjust` specification.
