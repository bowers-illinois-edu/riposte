# Critical values of the exact adaptive test

With `I1` pairs in the top `lambda[1]` fraction of ranks of `|y|` and
`I2` in the next band, `B1`, `B2` count positive differences, and
`T = 2*B1 + B2`. Under the sensitivity bound they are independent
binomials with success probability `kappa = gamma / (1 + gamma)`. The
test rejects when `B1 >= kB` or `T >= kT`, where `(kB, kT)` satisfy
conditions (4)-(6) of the paper: size at most `alpha`, neither constant
can be lowered by one without exceeding `alpha`, and
`|P(B1 >= kB) - P(T >= kT)|` is minimized. Ties in the minimization go
to the smaller `kB`, as in the paper.

## Usage

``` r
riposte_adaptive_critical(I1, I2, gamma = 1, alpha = 0.05)
```

## Arguments

- I1, I2:

  band sizes.

- gamma:

  sensitivity parameter, `>= 1`.

- alpha:

  level.

## Value

a list with `kB`, `kT`, `size` (P(reject)), `pB` (P(B1 \>= kB)), `pT`.
