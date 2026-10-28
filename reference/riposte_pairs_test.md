# Test one hypothesis with several matched-pair statistics

Convenience wrapper: builds the score columns with
[`riposte_pairs_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_pairs_scores.md)
and calls
[`riposte_pairs_bound()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_pairs_bound.md).

## Usage

``` r
riposte_pairs_test(
  y,
  statistics = c("u878", "u867"),
  gamma = 1,
  alternative = c("greater", "less"),
  dose = NULL,
  algorithm = NULL
)
```

## Arguments

- y:

  matched-pair differences (treated minus control).

- statistics:

  character vector from `"wilcoxon"`, `"brown"`, `"noether"`, `"u858"`,
  `"u888"`, `"u878"`, `"u868"`, `"u867"` (U-statistics with `m = 8`).
  Default: the two U-statistics `u878` and `u867` highlighted in
  Rosenbaum (2012).

- gamma:

  sensitivity parameter.

- alternative:

  `"greater"` or `"less"`.

- dose:

  optional nonnegative weights (e.g. doses) multiplying every score.

- algorithm:

  passed to
  [`riposte_pairs_bound()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_pairs_bound.md).

## Value

a
[`riposte_pairs_bound()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_pairs_bound.md)
list.
