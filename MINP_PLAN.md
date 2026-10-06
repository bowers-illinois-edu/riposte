# Plan: add single-step min-p as a fifth combination

riposte combines its six representations in four ways today: the
quadratic, the Cauchy, the max, and the hybrid. We will add a fifth,
single-step min-p. On every re-randomization it computes each
representation’s own p-value and keeps the smallest. Each
representation’s p-value is its mid-p value. The mid-p value counts
re-randomizations tied with the one in hand as half, and it is the value
riposte’s Cauchy and hybrid combinations already use (decision of 5
October 2026). The ordinary p-value counts ties in full. On the Wallsten
and Nteta outcomes it gave the same false positive rate and power as the
mid-p value, within simulation error. The mid-p value keeps every
per-representation p-value riposte reports the same quantity across
combinations. The combined p-value is the share of re-randomizations
whose smallest p-value is at most the observed smallest p-value.
Westfall and Young call this procedure single-step min-p. It is the
procedure in the EGAP guide “10 Things to Know About Multiple
Comparisons.” Rosenbaum’s 2012 adaptive test in the *Annals of Applied
Statistics* computes the same thing exactly for two statistics.

We will write the tests first. Then we will add the new combination to
the code that computes the combinations and let
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
call it, with and without covariance adjustment. Last, we will update
the documentation. We will include it in the simulations and studies. We
do not yet recommend it over the max.

## Why the max is not enough

The max takes the largest absolute standardized statistic. Suppose the
six absolute standardized statistics have the same null distribution, so
that a given value is reached in the same share of re-randomizations
whichever statistic reaches it. Then the smallest p-value is a
decreasing function of the largest absolute standardized statistic, so
min-p and the max give the same p-value. They give different p-values
when the null distributions differ. That happens when one
representation’s statistic takes only a few values, for example the
maximum distance on an outcome with five answer categories. In the
example in block_test_power’s multiple-comparisons memo, a statistic
that takes only three absolute values reaches a standardized value of
2.03 in 125 of 1,000 re-randomizations. A normal statistic reaches 2.03
or more in 37 of the same 1,000. The max treats the two values of 2.03
as equally extreme. Min-p compares their p-values, 0.125 and 0.037,
instead.

## What will be built

1.  One new function in `R/combine.R`, beside the function that computes
    the max from the statistics on every re-randomization
    (`riposte_max_from_T()`). It takes the same inputs and an
    `alternative`. For each representation it computes the mid-p value
    on every re-randomization in the requested direction, with the
    function the Cauchy combination already uses
    ([`riposte_midp()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_midp.md)).
    On each re-randomization it keeps the smallest of those p-values.
    The combined p-value is the share of all re-randomizations, observed
    assignment included, whose smallest p-value is at most the observed
    one. The existing permutation p-value function
    (`riposte_perm_pvalue()`) computes that share when it is given the
    negated smallest p-values, because it counts values at least as
    large as the observed one. The function returns four things: 1.1 the
    observed smallest p-value, as the statistic; 1.2 the combined
    p-value; 1.3 the name of the representation that gave the observed
    smallest p-value; 1.4 the six observed p-values, as `component_p`,
    the name the hybrid already uses.
2.  A matching wrapper, `riposte_minp()`, beside `riposte_max()`, for
    the tests to call directly.
3.  In
    [`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md):
    3.1 `"minp"` joins the choices for `statistic`. 3.2 The combination
    is called from both places that choose among the combinations: the
    path without covariance adjustment and
    `riposte_adjusted_combination()` in `R/adjust.R`. Both already
    compute the statistics on every re-randomization. Min-p needs
    nothing more. 3.3 A one-sided `alternative` is allowed, as for the
    max and the Cauchy, because each representation’s p-value can be
    computed in one direction. 3.4 The large-sample and saddlepoint
    engines stop with an error, as they do for the max. Under the normal
    approximation every standardized statistic is standard normal, so
    min-p and the max are the same test there and a large-sample min-p
    would add nothing.
4.  The p-value is exact in finite samples for any number of
    re-randomizations. Each representation’s p-value on each
    re-randomization is computed against all the re-randomizations,
    observed assignment included. So the collection of smallest p-values
    is the same whichever of those assignments was the observed one.
    When the treatment had no effect, the observed assignment is equally
    likely to be any of them, so its smallest p-value is equally likely
    to hold any rank in that collection. The permutation Cauchy
    combination is exact by the same argument.

## Tests, written before the code

Each test states the principle it checks in a comment.

1.  Exactness. On a small design whose assignments are all listed by
    `helper-enumerate.R`, the min-p p-value equals a direct count over
    every assignment, two-sided and one-sided. The one-sided max test in
    `test-one-sided.R` is the model.
2.  The identity. Take the raw outcome as one representation. As the
    second, take the same outcome values shuffled among the units within
    each block. The two have the same null distribution, because the sum
    over the treated units draws from the same values in each block, but
    they take different values on each assignment. On the full listing
    of assignments, min-p and the max give the same p-value for every
    assignment. With one representation, both equal that
    representation’s own p-value, as `test-exactness.R` already checks
    for the max.
3.  Where they differ. On the same listing, with the raw outcome and a
    representation that is 1 for the largest outcome in each block and 0
    otherwise, at least one assignment gives min-p and the max different
    p-values.
4.  The Bonferroni bracket. With outcomes that have no ties, the min-p
    p-value is at least the observed smallest p-value and at most six
    times it. The lower bound says min-p’s p-value is never smaller than
    the p-value of the best single representation tested alone. The
    upper bound says it is never larger than the Bonferroni-adjusted
    p-value.
5.  Level under the null, including a representation that takes two
    values, following `test-size.R`.
6.  Level under the null with covariance adjustment, in its own test in
    `test-adjust-validity.R` with 199 re-randomizations. The test beside
    it uses 49, and at 49 min-p cannot reject at 0.05: each
    representation’s most extreme re-randomization ties for the smallest
    mid-p value, so the min-p p-value cannot fall below about six
    divided by 50, or 0.12.
7.  The large-sample and saddlepoint engines give their errors. A
    one-sided call runs and reports `"minp"` as its combination.

## Documentation

We will add min-p wherever the combinations are listed: the description
of the `statistic` argument of
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md),
the comments at the top of `R/test.R` and `R/combine.R`, `NEWS.md`,
`README.md`, and `dev/riposte_spec.md`. The version goes from 0.0.0.9008
to 0.0.0.9009.

## Outside riposte

The power simulation in block_test_power,
`Analysis/truncation_review/multiplicity_sim.R`, computes its own
min-p. Once riposte has min-p, that script will call riposte instead and
be re-run, so the numbers in the multiple-comparisons memo come from the
package.

## Checkpoints

1.  After the tests are written, before any code.
2.  After the code is written, before `devtools::check()`.
3.  Whenever a design question arises that this plan does not settle.
