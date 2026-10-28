# Exact adaptive test for matched pairs (Brown + Noether)

Exact adaptive test for matched pairs (Brown + Noether)

## Usage

``` r
riposte_adaptive_test(y, gamma = 1, alpha = 0.05, lambda = c(1/3, 2/3))
```

## Arguments

- y:

  matched-pair differences (treated minus control). Zero differences are
  dropped. Ties among `|y|` are broken by order
  (`ties.method = "first"`) with a warning, because the paper assumes
  untied ranks.

- gamma:

  sensitivity parameter.

- alpha:

  level.

- lambda:

  the two band fractions `0 < lambda[1] < lambda[2] <= 1` (default
  `c(1/3, 2/3)`).

## Value

a list with `reject`, `B1`, `T`, `kB`, `kT`, `I1`, `I2`, `size`.

## References

Rosenbaum, P. R. (2012) An exact adaptive test with superior design
sensitivity in an observational study of treatments for ovarian cancer.
Annals of Applied Statistics 6(1), 83-105.
