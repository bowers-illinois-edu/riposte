# Mean pairwise L1 distance computed on the within-block ranks

Like
[`riposte_mean_dist()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_mean_dist.md)
but applied to the mid-ranks of `y`, so it measures how far each unit
sits from the others in rank space. Robust to the outcome's scale and to
outliers.

## Usage

``` r
riposte_mean_rank_dist(y)
```

## Arguments

- y:

  numeric vector (a single block's outcomes).

## Value

numeric vector the same length as `y`.
