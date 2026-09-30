# Control parameters for the quadratic-vs-Cauchy screen

Control parameters for the quadratic-vs-Cauchy screen

## Usage

``` r
riposte_screen_control(
  method = c("shrink", "threshold"),
  threshold = 100,
  lambda = NULL
)
```

## Arguments

- method:

  `"shrink"` (smooth, the default) or `"threshold"` (hard switch).

- threshold:

  condition-number cutoff for `method = "threshold"`.

- lambda:

  optional fixed shrinkage; if `NULL`, chosen from the conditioning of
  Sigma.

## Value

a list of screen settings.
