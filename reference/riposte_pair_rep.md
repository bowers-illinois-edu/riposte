# Matched-pair representation (Rosenbaum 2012)

A representation for
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
and
[`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md)
that scores each pair by the rank of its absolute difference among all
pairs, and gives that score to the unit with the larger outcome. Every
block must be a pair (two units, one treated). See
[`riposte_pairs_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_pairs_scores.md)
for the scores.

## Usage

``` r
riposte_pair_rep(
  type = c("wilcoxon", "brown", "noether", "ustat"),
  m = 8,
  m1 = 7,
  m2 = 8,
  exact = TRUE
)
```

## Arguments

- type:

  `"wilcoxon"` (ranks of `|y|`), `"brown"` (Brown / Markowski-
  Hettmansperger combined quantile average, quantiles 1/3 and 2/3),
  `"noether"` (0/1 for ranks in the top third), or `"ustat"`
  (Rosenbaum's 2011 U-statistic, needs `m`, `m1`, `m2`).

- m, m1, m2:

  U-statistic indices with `1 <= m1 <= m2 <= m`; `(8, 7, 8)` and
  `(8, 6, 7)` are the two highlighted in Rosenbaum (2012).

- exact:

  for `"ustat"`, use the integer form
  `sum choose(a-1, l-1) * choose(I-a, m-l)` of the paper (`TRUE`) or the
  smooth approximation (`FALSE`, what
  [`testtwice::multrnk()`](https://rdrr.io/pkg/testtwice/man/multrnk.html)
  does by default).

## Value

a function of `(y, block)` of class `"riposte_pooled_rep"`, carrying the
pair-level score function as attribute `pair_scores`.
