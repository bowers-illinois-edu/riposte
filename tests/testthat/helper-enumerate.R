## helper-enumerate.R --- exhaustive within-block randomization enumerator.
## Sourced automatically by testthat before the tests. Used by the exactness
## tests (block- and cluster-level) to build the FULL randomization distribution
## on tiny designs and check riposte's p-values against the exact value.

## All within-block assignments for a balanced design of blocks of size nb with
## nb/2 treated each: the Cartesian product of each block's choose(nb, nb/2) ways.
## Returns a matrix with one column per assignment.
enumerate_assignments <- function(block) {
  block <- as.factor(block)
  per_block <- lapply(levels(block), function(b) {
    ix <- which(block == b); nb <- length(ix); k <- nb %/% 2L
    combn(nb, k, function(t) { z <- integer(nb); z[t] <- 1L; z }, simplify = FALSE)
  })
  grid <- do.call(expand.grid, lapply(per_block, seq_along))
  vapply(seq_len(nrow(grid)), function(r) {
    z <- integer(length(block))
    for (j in seq_along(per_block)) z[block == levels(block)[j]] <- per_block[[j]][[grid[r, j]]]
    z
  }, numeric(length(block)))
}
