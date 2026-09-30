# Stephenson rank-score representations at several tuning values

Returns a named list of representation functions, one per tuning value
`r`, suitable for
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
or
[`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md).
Combining them (the default quadratic or the screen) gives the combined
Stephenson rank test, which targets the upper tail across a range of
weightings instead of committing to one. Add the default representations
with `c(riposte_reps_default(), riposte_stephenson_reps())` to combine
distribution-shape and upper-tail views together.

## Usage

``` r
riposte_stephenson_reps(r = c(2, 6, 10))
```

## Arguments

- r:

  vector of whole-number tuning values; `c(2, 6, 10)` by default
  (Wilcoxon plus two heavier upper-tail weightings), following the
  proposal.

## Value

a named list of representation functions, named `S2`, `S6`, ...

## See also

[`riposte_stephenson_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_stephenson_scores.md),
[`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md).
