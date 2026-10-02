# Polynomial rank scores

For an outcome vector of length `n`, the polynomial rank score of a unit
at ascending rank `k` is `(k / (n + 1))^(zeta - 1)`, with ties at their
average rank. This is the score of Kim, Su, Bowers, and Li (arXiv
2605.08027). At `zeta = 2` it is the rank divided by `n + 1`, which
gives the Wilcoxon rank-sum test. Larger `zeta` puts more of the block's
total score on its top ranks: treating `k / (n + 1)` as uniform on (0,
1), the score-weighted average distance from the top of the block is
about `n / (zeta + 1)` units.

## Usage

``` r
riposte_poly_scores(y, zeta, tail = c("upper", "lower"))
```

## Arguments

- y:

  numeric vector (a single block's outcomes).

- zeta:

  the tuning value, a single number `>= 1`. It need not be a whole
  number. `zeta = 1` gives a constant, which carries no information.

- tail:

  `"upper"` (default) weights the top of the block; `"lower"` weights
  the bottom.

## Value

numeric vector of scores the same length as `y`.

## Details

Unlike
[`riposte_stephenson_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_stephenson_scores.md),
the score lies in `[0, 1)` whatever the block size, so blocks of
different sizes enter on the same scale.

With `tail = "lower"` the score weights the bottom of the block instead:
the unit at ascending rank `k` gets
`-((n + 1 - k) / (n + 1))^(zeta - 1)`, which lies in `(-1, 0]`. The
minus sign makes the score rise with the outcome, as the upper-tail
score does, so that with `alternative = "less"` in \[riposte_test()\]
both tails ask whether treated outcomes are lower. Without ties, the
lower-tail score is `-riposte_poly_scores(-y, zeta)`.

\[0, 1)\` whatever the block size, so blocks of different sizes enter on
the same scale.

With `tail = "lower"` the score weights the bottom of the block instead:
the unit at ascending rank `k` gets
`-((n + 1 - k) / (n + 1))^(zeta - 1)`, which lies in \`(-1, 0\]:
R:0,%201)%60%20whatever%20the%0Ablock%20size,%20so%20blocks%20of%20different%20sizes%20enter%20on%20the%20same%20scale.%0A%0AWith%20%60tail%20=%20%22lower%22%60%20the%20score%20weights%20the%20bottom%20of%20the%20block%20instead:%20the%0Aunit%20at%20ascending%20rank%20%60k%60%20gets%20%60-((n%20+%201%20-%20k)%20/%20(n%20+%201))%5E(zeta%20-%201)%60,%20which%0Alies%20in%20%60(-1,%200
\[riposte_test()\]: R:riposte_test()

## See also

[`riposte_poly_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_reps.md),
[`riposte_stephenson_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_stephenson_scores.md).
