# Mid-p permutation p-values

Given a statistic evaluated across the observed assignment and its
re-randomizations (a single vector `v`, length 1 + number of draws),
returns the mid-p permutation p-value of each entry relative to the
whole set: two-sided by default, or one-sided with
`alternative = "greater"` (large values are extreme) or `"less"` (small
values are extreme). The most extreme value (largest `|v|`) gets
`0.5 / n`; the least extreme gets `(n - 0.5) / n`. Every value is
strictly inside (0, 1), so the Cauchy transform never reaches a pole.

## Usage

``` r
riposte_midp(v, alternative = c("two.sided", "greater", "less"))
```

## Arguments

- v:

  numeric vector of a statistic across assignments (observed plus
  draws). Two-sidedness is by `abs(v)`.

- alternative:

  `"two.sided"` (default), `"greater"`, or `"less"`.

## Value

numeric vector of mid-p values, the same length as `v`, all in (0, 1).
Without ties, the `"greater"` and `"less"` values add to one.

## Details

The mid-p value is the ordinary upper-tail count minus half the mass at
the observed value:
`(# strictly more extreme) + 0.5 * (# tied, including self)`, divided by
`n`. Values within about `1.5e-8` of each other, relative to their size
when it exceeds one, count as tied, because rank scores make many
re-randomizations tie exactly, and floating-point summation in a
different order can separate such ties in their last bits.
