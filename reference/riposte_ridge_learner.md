# Controls-only ridge learner

Returns a learner that fits ridge regression on the control units and
predicts for all units. Ridge has a closed form (no iterative fitting)
and its shrinkage keeps predictions stable when there are many
covariates — the regime where an unregularized fit overfits and the
fit-once shortcut loses exactness.

## Usage

``` r
riposte_ridge_learner(lambda = NULL)
```

## Arguments

- lambda:

  ridge penalty. If `NULL`, chosen by leave-one-out CV on the control
  units each time the learner is called. A fixed `lambda` (for instance
  the CV value from the observed controls) is faster and equally exact
  under refit-per-permutation.

## Value

a function `(X, Y, control) -> Yhat`.
