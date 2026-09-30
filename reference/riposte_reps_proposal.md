# The proposal's representation set (raw tanh)

[`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md)
with Huber's psi replaced by the raw `tanh(y)` that the "Power for
Policy" proposal used (as do `manytestsr`'s `pIndepDist()` and
`pCombCauchyDist()`), for reproducing the proposal's tables exactly.
Prefer
[`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md)
for real outcomes, whose scale `tanh(y)` does not respect.

## Usage

``` r
riposte_reps_proposal()
```

## Value

a named list of representation functions.
