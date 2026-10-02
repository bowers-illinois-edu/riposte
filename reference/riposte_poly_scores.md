# Polynomial rank scores

For an outcome vector of length `n`, the polynomial rank score of a unit
at ascending rank `k` is `(k / (n + 1))^(zeta - 1)`, with ties at their
average rank. This is the score of Kim, Su, Bowers, and Li (arXiv
2605.08027). At `zeta = 2` it is the rank divided by `n + 1`, which
gives the Wilcoxon rank-sum test. Larger `zeta` puts more of the block's
total score on its top ranks: treating `k / (n + 1)` as uniform on (0,
1), the score-weighted average distance from the top of the block is
about `n / (zeta + 1)` units.

## Usage

``` r
riposte_poly_scores(y, zeta)
```

## Arguments

- y:

  numeric vector (a single block's outcomes).

- zeta:

  the tuning value, a single number `>= 1`. It need not be a whole
  number. `zeta = 1` gives the constant 1, which carries no information.

## Value

numeric vector of scores the same length as `y`.

## Details

Unlike
[`riposte_stephenson_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_stephenson_scores.md),
the score lies in `[0, 1)` whatever the block size, so blocks of
different sizes enter on the same scale.

## See also

[`riposte_poly_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_reps.md),
[`riposte_stephenson_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_stephenson_scores.md).
