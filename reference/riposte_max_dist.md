# Maximum distance of each unit to the block extremes

For each unit, `max(y_i - min(y), max(y) - y_i)`: how far the unit is
from the nearer end of its block's range. Translation invariant.

## Usage

``` r
riposte_max_dist(y)
```

## Arguments

- y:

  numeric vector (a single block's outcomes).

## Value

numeric vector the same length as `y`.
