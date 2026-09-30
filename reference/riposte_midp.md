# Two-sided mid-p permutation p-values

Given a statistic evaluated across the observed assignment and its
re-randomizations (a single vector `v`, length 1 + number of draws),
returns the two-sided mid-p permutation p-value of each entry relative
to the whole set. The most extreme value (largest `|v|`) gets `0.5 / n`;
the least extreme gets `(n - 0.5) / n`. Every value is strictly inside
(0, 1), so the Cauchy transform never reaches a pole.

## Usage

``` r
riposte_midp(v)
```

## Arguments

- v:

  numeric vector of a statistic across assignments (observed plus
  draws). Two-sidedness is by `abs(v)`.

## Value

numeric vector of mid-p values, the same length as `v`, all in (0, 1).

## Details

The mid-p value is the ordinary upper-tail count minus half the mass at
the observed value:
`(# strictly more extreme) + 0.5 * (# tied, including self)`, divided by
`n`. With `ties.method = "average"` the rank-based form below gives
exactly that. The linear statistics riposte combines are continuous, so
ties are negligible in practice.
