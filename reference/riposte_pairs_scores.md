# Signed-rank-type scores for matched-pair differences

Scores `h_i` that depend on the rank of `|y_i|`, for use with
[`riposte_pairs_bound()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_pairs_bound.md)
and
[`riposte_pairs_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_pairs_test.md).

## Usage

``` r
riposte_pairs_scores(
  y,
  type = c("wilcoxon", "brown", "noether", "ustat"),
  m = 8,
  m1 = 7,
  m2 = 8,
  q1 = NULL,
  q2 = NULL,
  exact = TRUE
)
```

## Arguments

- y:

  numeric vector of treated-minus-control matched-pair differences.

- type:

  `"wilcoxon"` (ranks of `|y|`), `"brown"` (Brown / Markowski-
  Hettmansperger combined quantile average, quantiles 1/3 and 2/3),
  `"noether"` (0/1 for ranks in the top third), or `"ustat"`
  (Rosenbaum's 2011 U-statistic, needs `m`, `m1`, `m2`).

- m, m1, m2:

  U-statistic indices with `1 <= m1 <= m2 <= m`; `(8, 7, 8)` and
  `(8, 6, 7)` are the two highlighted in Rosenbaum (2012).

- q1, q2:

  quantiles for `"brown"` (defaults 1/3, 2/3) and `"noether"` (`q1`
  only, default 2/3).

- exact:

  for `"ustat"`, use the integer form
  `sum choose(a-1, l-1) * choose(I-a, m-l)` of the paper (`TRUE`) or the
  smooth approximation (`FALSE`, what
  [`testtwice::multrnk()`](https://rdrr.io/pkg/testtwice/man/multrnk.html)
  does by default).

## Value

a numeric vector of scores, scaled to have maximum 1 (the standardized
statistics do not depend on the scale), with 0 for zero differences.
