# Maximum-of-L sensitivity bound (Rosenbaum 2012, procedure 6)

Given L score columns `H` for the same matched pairs, computes each
statistic `T_l = sum B_i H[i, l]`, its standardized worst-case deviate
at sensitivity `gamma`, and the upper bound on the one-sided p-value of
their maximum from the multivariate normal orthant probability implied
by `crossprod(H)`.

## Usage

``` r
riposte_pairs_bound(y, H, gamma = 1, algorithm = NULL)
```

## Arguments

- y:

  matched-pair differences.

- H:

  numeric matrix, one row per pair, one column per statistic (`2 <= L`).

- gamma:

  sensitivity parameter, `>= 1`.

- algorithm:

  an `mvtnorm` algorithm object. Default:
  [`mvtnorm::Miwa()`](https://rdrr.io/pkg/mvtnorm/man/algorithms.html)
  for `L <= 5` (what `testtwice` uses), else
  [`mvtnorm::GenzBretz()`](https://rdrr.io/pkg/mvtnorm/man/algorithms.html)
  (Miwa is exponentially slow in L; GenzBretz is Monte Carlo, so set a
  seed).

## Value

a list with `p.value`, `deviates`, and `cor`.

## Details

This is a normal approximation (as in the paper and in
[`testtwice::tt()`](https://rdrr.io/pkg/testtwice/man/tt.html)).
[`riposte_pairs_perm_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_pairs_perm_test.md)
computes the same bound from the distribution of the signs instead,
which is what
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
does unless `engine = "asymptotic"`.

## References

Rosenbaum, P. R. (2012) Testing one hypothesis twice in observational
studies. Biometrika 99(4), 763-774. See also the CRAN package
`testtwice`.
