# Sign-pattern test for the maximum of L statistics, with a sensitivity bound

Under the sharp null in a paired randomized experiment, each
difference's sign is equally likely to be either, independently across
pairs. Under Rosenbaum's sensitivity model with hidden bias at most
`gamma`, each sign is instead positive with some chance no larger than
`theta = gamma / (1 + gamma)`. Because the scores are nonnegative, the
chance that the maximum standardized statistic reaches its observed
value is largest when every sign is positive with chance `theta`, so
that chance is an upper bound on the p-value (and the exact p-value at
`gamma = 1`). This function computes it from the distribution of the
signs: by listing every sign pattern with its chance when `nresample` is
`NULL` (at most 16 pairs with nonzero difference), otherwise from
`nresample` random patterns with the observed signs counted as one of
them.
[`riposte_pairs_bound()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_pairs_bound.md)
gives the normal approximation to the same bound.

## Usage

``` r
riposte_pairs_perm_test(y, H, nresample = 9999L, gamma = 1)
```

## Arguments

- y:

  matched-pair differences.

- H:

  numeric matrix, one row per pair, one column per statistic (`2 <= L`).

- nresample:

  number of random sign vectors, or `NULL` for full enumeration
  (requires at most 16 nonzero pairs).

- gamma:

  sensitivity parameter, `>= 1`; `1` (the default) gives the
  randomization test of the sharp null.

## Value

a list with `p.value` (an upper bound when `gamma > 1`) and the observed
`statistic`, the maximum of the statistics standardized at `theta`.

## References

Rosenbaum, P. R. (2012) Testing one hypothesis twice in observational
studies. Biometrika 99(4), 763-774.
