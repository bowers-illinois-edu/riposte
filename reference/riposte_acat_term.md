# Pole-aware Cauchy (ACAT) transform of p-values

`tan((0.5 - p) * pi)`, with the accurate small-/large-p limits
substituted near the poles (this is what Liu et al.'s ACAT
implementation does): `1 / (p * pi)` as `p -> 0` and
`-1 / ((1 - p) * pi)` as `p -> 1`. With mid-p inputs the substitution
almost never fires; it makes the transform safe if riposte is ever asked
to combine p-values produced elsewhere that reach 0 or 1.

## Usage

``` r
riposte_acat_term(p)
```

## Arguments

- p:

  numeric vector of p-values.

## Value

numeric vector of Cauchy-transformed values, the same length as `p`.

## Examples

``` r
# The paper's Hybrid test gives equal weight to seven p-values:
# six individual representations and their quadratic combination.
dat <- data.frame(
  outcome = c(-3, -1, 0, 0.5, 2, 4, 7, 9,
              -2, -1, 0, 1, 2, 3, 4, 5, 8, 13),
  treated = c(1, 0, 1, 0, 0, 1, 0, 0,
              1, 0, 1, 1, 0, 1, 0, 1, 1, 0),
  block = factor(rep(1:2, c(8, 10)))
)
six <- riposte_test(outcome ~ treated | block, dat,
                    statistic = "cauchy", engine = "asymptotic")
quad <- riposte_test(outcome ~ treated | block, dat,
                     statistic = "quadratic", engine = "asymptotic")
hybrid_inputs <- c(six$component_p, quadratic = quad$p.value)
hybrid_statistic <- mean(riposte_acat_term(hybrid_inputs))
pcauchy(hybrid_statistic, lower.tail = FALSE)
#> [1] 0.7383777

# These are seven inputs, not the two combined p-values with equal weight.
# The final Cauchy tail is an approximation, not a permutation p-value.
```
