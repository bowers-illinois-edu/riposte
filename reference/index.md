# Package index

## Tests

The user-facing test and per-representation diagnostic.

- [`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
  : Combined randomization-based test of the sharp null
- [`riposte_components()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_components.md)
  : Per-representation p-values, covariance, and condition number

## Outcome representations

The views of the outcome that get combined; an open, user-supplied set.

- [`riposte_reps_default()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_default.md)
  : The default set of outcome representations
- [`riposte_reps_proposal()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_reps_proposal.md)
  : The proposal's representation set (raw tanh)
- [`riposte_mean_dist()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_mean_dist.md)
  : Mean pairwise L1 distance of each unit to the others (within a
  block)
- [`riposte_mean_rank_dist()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_mean_rank_dist.md)
  : Mean pairwise L1 distance computed on the within-block ranks
- [`riposte_max_dist()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_max_dist.md)
  : Maximum distance of each unit to the block extremes
- [`riposte_rank()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_rank.md)
  : Mid-ranks of a block's outcomes
- [`riposte_huber()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_huber.md)
  : Huber's psi of a block's outcomes, standardised by the block's MAD
- [`riposte_tanh()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_tanh.md)
  : Bounded (tanh) transform of a block's outcomes, robustly
  standardised
- [`riposte_stephenson_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_stephenson_scores.md)
  : Stephenson rank scores
- [`riposte_stephenson_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_stephenson_reps.md)
  : Stephenson rank-score representations at several tuning values
- [`riposte_poly_scores()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_scores.md)
  : Polynomial rank scores
- [`riposte_poly_reps()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_poly_reps.md)
  : Polynomial rank-score representations at several tuning values

## Combination screen

Choosing or blending the quadratic and the Cauchy from the conditioning
of the covariance.

- [`riposte_screen()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_screen.md)
  : Choose (or blend) the combination from the conditioning of Sigma
- [`riposte_screen_control()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_screen_control.md)
  : Control parameters for the quadratic-vs-Cauchy screen

## Covariance adjustment

Controls-only learners, refit inside every re-randomization.

- [`riposte_adjust()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_adjust.md)
  : Covariance-adjustment specification
- [`riposte_ridge_learner()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_ridge_learner.md)
  : Controls-only ridge learner
- [`riposte_lm_learner()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_lm_learner.md)
  : Controls-only ordinary-least-squares learner

## Building blocks

The permutation engine and the closed-form moments behind the tests.

- [`riposte_sw_moments()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_sw_moments.md)
  : Closed-form Strasser-Weber permutation mean and covariance
- [`riposte_block_draws()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_block_draws.md)
  : Within-block re-randomizations preserving each block's treated count
- [`riposte_score_matrix()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_score_matrix.md)
  : Build the within-block-centred score matrix from representations
- [`riposte_midp()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_midp.md)
  : Two-sided mid-p permutation p-values
- [`riposte_acat_term()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_acat_term.md)
  : Pole-aware Cauchy (ACAT) transform of p-values
