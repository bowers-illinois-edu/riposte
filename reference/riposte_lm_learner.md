# Controls-only ordinary-least-squares learner

Returns a learner that fits OLS on the control units. With many
covariates relative to the number of controls, OLS overfits — which is
exactly the case that makes the fit-once / fixed-residual shortcut
anti-conservative and that refit-per-permutation handles correctly.
Useful for demonstrating that contrast; ridge is the better default in
practice.

## Usage

``` r
riposte_lm_learner()
```

## Value

a function `(X, Y, control) -> Yhat`.
