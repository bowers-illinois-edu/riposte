# Stephenson rank scores

For an outcome vector, the Stephenson score of a unit at ascending rank
`j` is `choose(j - 1, r - 1)`. Ties take average ranks. Larger `r`
concentrates weight on the top ranks; `r=2` reduces to the (rank-1)
Wilcoxon score.

## Usage

``` r
riposte_stephenson_scores(y, r)
```

## Arguments

- y:

  numeric vector (a single block's outcomes).

- r:

  the Stephenson tuning value, a whole number; `r=2` is the Wilcoxon
  (rank - 1) score and `r=1` gives a constant (uninformative) score.
  Must be a whole number: a fractional `r` would silently round inside
  [`choose()`](https://rdrr.io/r/base/Special.html).

## Value

numeric vector of scores the same length as `y`.
