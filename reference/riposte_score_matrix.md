# Build the within-block-centred score matrix from representations

Applies each representation function to the outcome within each block,
then centres each resulting score within block (subtract the block
mean). Returns the matrix of centred scores, one column per surviving
representation. A representation that is constant within every block
contributes nothing to a permutation statistic and is dropped, with a
recorded reason.

## Usage

``` r
riposte_score_matrix(y, block, representations = riposte_reps_default())
```

## Arguments

- y:

  numeric outcome vector.

- block:

  factor or vector of block labels, the same length as `y`.

- representations:

  a named list of representation functions; see
  [`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md).

## Value

a list with

- scores:

  numeric matrix, `length(y)` rows, one column per kept representation,
  centred within block.

- kept:

  character vector of representation names retained.

- dropped:

  character vector of representation names dropped for being
  within-block constant.
