# Cosine of the within-block rank

Each unit gets \\\cos(\pi r / (n_b + 1))\\, where \\r\\ is the rank of
its outcome among the \\n_b\\ outcomes of its block (ties at their
average rank). The score runs from near 1 at the lowest outcome to near
-1 at the highest. It is the score of the paper's cosine statistic
(Bowers and Burton, Section 4.7), which riposte computes as a block
statistic through
`riposte_minp_control(block = list(cosine = riposte_cosine))`. Within a
block of fifty units it correlates with the rank at -0.99, so as one of
the representations added over blocks it nearly repeats the rank.

## Usage

``` r
riposte_cosine(y)
```

## Arguments

- y:

  numeric vector (a single block's outcomes).

## Value

numeric vector the same length as `y`.

## Examples

``` r
riposte_cosine(c(3.2, -1, 0.5, 7))
#> [1] -0.309017  0.809017  0.309017 -0.809017
```
