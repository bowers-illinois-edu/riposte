# Rosenbaum's matched-pair representations

Wilcoxon, Brown, Noether and the two U-statistics `(8, 7, 8)` and
`(8, 6, 7)` highlighted in Rosenbaum (2012), as a named list for the
`representations` argument of
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md).

## Usage

``` r
riposte_pair_reps(which = c("wilcoxon", "brown", "noether", "u878", "u867"))
```

## Arguments

- which:

  names to keep, from `"wilcoxon"`, `"brown"`, `"noether"`, `"u878"`,
  `"u867"`, `"u858"`, `"u888"`, `"u868"`.

## Value

a named list of
[`riposte_pair_rep()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_pair_rep.md)
functions.
