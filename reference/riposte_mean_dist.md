# Mean pairwise L1 distance of each unit to the others (within a block)

For an outcome vector `y` of length n, returns, for each unit i, the
mean absolute distance `mean_{j != i} |y_i - y_j|` with denominator
`n - 1`. This is the energy/distance representation: units far from the
rest of their block (in either tail) score high.

## Usage

``` r
riposte_mean_dist(y)
```

## Arguments

- y:

  numeric vector (a single block's outcomes).

## Value

numeric vector the same length as `y`. A block of length \< 2 has no
distances and returns zeros.

## Details

The computation is the sorted closed form (no n-by-n distance matrix),
made numerically stable by first subtracting the median. L1 distances
are translation invariant, so the subtraction changes nothing about the
answer; it keeps the prefix sums small so the fast form matches the
definitional sum-of-distances on the same input to machine precision
even when `y` carries a large location (dollars, populations). It does
not recover precision already lost when `y` is stored as a large offset
plus a small spread — no algorithm can — but unlike a plain
cumsum-then-subtract it adds no error of its own.
