# The default set of outcome representations

Returns the six representations of Bowers and Burton, "A More Powerful
Test for Randomized Experiments" (working paper), as a named list of
functions, each mapping a block's outcome vector to a per-unit score
vector of the same length:

- raw:

  the outcome itself (identity)

- rank:

  within-block mid-ranks,
  [`riposte_rank()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_rank.md)

- mean_dist:

  mean L1 distance,
  [`riposte_mean_dist()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_mean_dist.md)

- mean_rank_dist:

  mean L1 distance on ranks,
  [`riposte_mean_rank_dist()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_mean_rank_dist.md)

- max_dist:

  max distance to the block extremes,
  [`riposte_max_dist()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_max_dist.md)

- huber:

  Huber's psi of the MAD-standardised outcome,
  [`riposte_huber()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_huber.md)

## Usage

``` r
riposte_reps_default()
```

## Value

a named list of representation functions.

## Details

Pass your own list (or this one extended) to
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
or
[`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md)
to combine a different set. Each element must be a function of one
numeric vector returning a numeric vector of the same length; it is
applied within each block.

## See also

[`riposte_reps_proposal()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_proposal.md)
for the earlier set with raw `tanh(y)`.
