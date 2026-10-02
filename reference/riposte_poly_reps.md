# Polynomial rank-score representations at several tuning values

Returns a named list of representation functions, one per value of
`zeta`, for
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
or
[`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md).
Combining them gives a test that responds to gains concentrated at the
top of each block across several weightings instead of one.

## Usage

``` r
riposte_poly_reps(zeta = c(2, 7, 12, 17, 22), tail = c("upper", "lower"))
```

## Arguments

- zeta:

  vector of tuning values, each `>= 1`.

- tail:

  `"upper"` (default) or `"lower"`, passed to
  [`riposte_poly_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_scores.md).

## Value

a named list of representation functions, named `poly2`, `poly7`, ...
for the upper tail and `polylow2`, `polylow7`, ... for the lower tail.

## Details

The default, `c(2, 7, 12, 17, 22)`, is the set in Bowers and Burton's
rank-score tables, chosen for blocks of 50 units with 25 treated. For
other designs, supply your own values; the vignette shows how to choose
them. In brief: the largest `zeta` near `n_b / m - 1` targets gains
confined to about `m` units per block; and two scores at `zeta` and
`zeta'` correlate at about `2 sqrt(u v) / (u + v)` with `u = 2 zeta - 1`
and `v = 2 zeta' - 1`, so values whose `2 zeta - 1` grow by a constant
factor (for example `zeta = 2, 5, 14, 41`, where `2 zeta - 1` triples)
are evenly spaced, while equal steps in `zeta` give neighbours whose
correlation approaches 1.

`zeta = 2` is the rank rescaled. When it is added to
[`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md),
which already holds the rank, leave it out:
`c(riposte_reps_default(), riposte_poly_reps(c(7, 12, 17, 22)))`.

With `tail = "lower"` the scores weight the bottom of each block, for a
treatment that lowers the outcomes of a few units; see
[`riposte_poly_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_scores.md)
and
[`vignette("harm")`](https://bowers-illinois-edu.github.io/riposte/articles/harm.md).

## See also

[`riposte_poly_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_scores.md),
[`riposte_stephenson_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_stephenson_reps.md).
