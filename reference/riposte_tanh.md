# Bounded (tanh) transform of a block's outcomes, robustly standardised

A bounded representation: `tanh((y - median(y)) / s)`, where `s` is a
robust scale (the MAD, falling back to the SD and then to 1 when those
are zero).

## Usage

``` r
riposte_tanh(y)
```

## Arguments

- y:

  numeric vector (a single block's outcomes).

## Value

numeric vector the same length as `y`.

## Details

The standardisation matters. Raw `tanh(y)` saturates to +/-1 once `|y|`
exceeds about 3, so on an outcome measured in dollars or raw test points
every value maps to +/-1 and the representation carries no information.
Centring and scaling first keeps the transform sensitive in the body of
each block's distribution while still bounding the influence of extreme
values. To reproduce the proposal's raw `tanh(y)` exactly, use
[`riposte_reps_proposal()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_proposal.md).

This transform is not in either built-in set:
[`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md)
uses
[`riposte_huber()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_huber.md)
in the sixth place. Add it to a list of representations to use it.
