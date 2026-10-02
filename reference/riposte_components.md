# Per-representation p-values, covariance, and condition number

Per-representation p-values, covariance, and condition number

## Usage

``` r
riposte_components(
  formula,
  data,
  blocks = NULL,
  clusters = NULL,
  representations = riposte_reps_default(),
  nresample = 1999L,
  cluster_agg = mean,
  seed = NULL,
  alternative = c("two.sided", "greater", "less"),
  ...
)
```

## Arguments

- formula:

  `Y ~ treatment` for complete randomization, or `Y ~ treatment | block`
  for randomization within blocks (or use `blocks`).

- data:

  a data.frame/data.table.

- blocks:

  optional block specification (see
  [`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)).

- clusters:

  optional cluster specification (see
  [`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md));
  the design is collapsed to the cluster level when given.

- representations:

  a named list of representation functions;
  [`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md)
  by default.

- nresample:

  number of re-randomizations for the marginal p-values.

- cluster_agg:

  how to aggregate the outcome within a cluster; the cluster mean by
  default.

- seed:

  optional integer seed (L'Ecuyer-CMRG) for reproducibility.

- alternative:

  `"two.sided"` (default), `"greater"`, or `"less"`: the direction of
  each representation's mid-p value; see
  [`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md).

- ...:

  reserved.

## Value

an object of class `riposte_components` with the representations' mid-p
permutation p-values in the requested direction, the closed-form
covariance `Sigma`, its `condition` number, the kept/dropped
representations, and the score matrix.
