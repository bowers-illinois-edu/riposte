# Rosenbaum's p-value for the exact adaptive test

The smallest level `alpha` at which the adaptive test, with its cutoffs
chosen by conditions (4)-(6) of Rosenbaum (2012), rejects the observed
counts. Rosenbaum defines the test's p-value this way: "The P-value ...
for the adaptive test is the smallest alpha in (4)-(6) that leads to
rejection" (p. 100). Because the cutoffs are chosen separately at each
level, this number can fall at or below a level more often than that
level allows when there is no effect: with 11 untied pairs at
`gamma = 1`, it is at most 16/256 with chance 20/256.
[`riposte_adaptive_pvalue()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_adaptive_pvalue.md)
gives a p-value that keeps its level, and
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
prints both.

## Usage

``` r
riposte_adaptive_smallest_level(y, gamma = 1, lambda = c(1/3, 2/3))
```

## Arguments

- y:

  matched-pair differences (treated minus control). Zero differences are
  dropped. Ties among `|y|` are broken by order
  (`ties.method = "first"`) with a warning, because the paper assumes
  untied ranks.

- gamma:

  sensitivity parameter.

- lambda:

  the two band fractions `0 < lambda[1] < lambda[2] <= 1` (default
  `c(1/3, 2/3)`).

## Value

the smallest rejecting level, in (0, 1\]; 1 when no level below 1
rejects.

## References

Rosenbaum, P. R. (2012) An exact adaptive test with superior design
sensitivity in an observational study of treatments for ovarian cancer.
Annals of Applied Statistics 6(1), 83-105.
