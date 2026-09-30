# Within-block re-randomizations preserving each block's treated count

Builds a 0/1 assignment matrix whose first column is the observed
treatment and whose remaining `nresample` columns are independent
re-randomizations, each drawn by
[`sample()`](https://rdrr.io/r/base/sample.html) within block so that
every block keeps its observed number of treated units. This is the
randomization reference for the Monte Carlo combinations.

## Usage

``` r
riposte_block_draws(z, block, nresample = 1999L)
```

## Arguments

- z:

  integer/numeric 0/1 vector of observed treatment assignment.

- block:

  factor or vector of block labels, the same length as `z`.

- nresample:

  number of re-randomizations (columns 2..nresample+1).

## Value

a numeric matrix with `length(z)` rows and `nresample + 1` columns;
column 1 is `z`. Row order matches the input.

## Details

For a cluster-randomized design, pass cluster-level vectors (`z` and
`block` indexed by cluster, one entry per cluster); the permutation is
then over cluster assignments, which is the correct reference for
clusters.
