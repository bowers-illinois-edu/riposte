# Exact level for the covariance-adjusted combined test

Status: settled, with proof and simulation evidence (2026-06-24). Resolves the
two open theory items previously listed in the development notes ("adjusted-quadratic moments"
and "screen under adjustment"). Written in the design-based, finite-population
tradition (Rosenbaum; Hansen).

## The question

riposte's unadjusted quadratic uses the closed-form Strasser-Weber permutation
moments, which are exact. Under covariance adjustment the outcome is residualized
against a learner that is refit on each re-randomization, so the representation
scores change from draw to draw and the closed-form moments no longer apply. The
package estimates the metric (mean and covariance) from the draws instead. The
worry was that this Monte-Carlo metric makes the adjusted quadratic
anti-conservative. It does --- but only because of HOW the metric was estimated
(from the draws alone, excluding the observed assignment). Estimating it from the
pooled set (observed and draws together) restores exact level. This note states
and proves that, and reports the simulations that confirm it.

## Setup (object and measure ledger)

- Finite population of units `i = 1, ..., N`, partitioned into blocks. In block
  `b`, `m_b` of `n_b` units are treated.
- Potential outcomes `y_i(0), y_i(1)` and covariates `x_i` are FIXED constants.
- Sharp null `H0`: `y_i(1) = y_i(0)` for every `i`. Under `H0` the observed
  outcome `Y_i = y_i(0)` does not depend on the assignment; `Y = (Y_1,...,Y_N)`
  is a fixed vector. (For a constant-shift null `H0(t0)`: `y_i(1) = y_i(0) + t0`,
  replace `Y` by the imputed `Y - t0 z`, also fixed; take `t0 = 0` below.)
- Randomness comes ONLY from the assignment. `P` is the block-randomization law:
  independently across blocks, the treated set is uniform among the `C(n_b, m_b)`
  choices. Let `Z0` be the realized (observed) assignment and `Z1, ..., ZB` be `B`
  independent draws from `P`, independent of `Z0`.
- A learner `L` maps `(X, Y, C)` to fitted values for all units, fit on the
  control set `C`. `L` is any MEASURABLE map; no smoothness, consistency, or
  correct specification is assumed. The residual at assignment `z` is
  `e(z) = Y - L(X, Y, {i : z_i = 0})`.
- Representations `f_1, ..., f_C`. The within-block-centred score for
  representation `c` at assignment `z` is `s_{c,i}(z) = f_c(e(z) restricted to
  i's block) - (its block mean)`. The linear statistic is
  `T_c(z) = sum_i z_i s_{c,i}(z)`, and `T(z) = (T_1(z), ..., T_C(z))` in `R^C`.
- Pooled metric over the `B + 1` assignments `Z0, ..., ZB`:
  `mu_hat = (B+1)^{-1} sum_{j=0}^B T(Zj)` and
  `Sigma_hat = B^{-1} sum_{j=0}^B (T(Zj) - mu_hat)(T(Zj) - mu_hat)'`.
- Combined statistic. For the quadratic,
  `Q_j = (T(Zj) - mu_hat)' Sigma_hat^+ (T(Zj) - mu_hat)` (`^+` is the
  Moore-Penrose inverse). The max and Cauchy combinations are defined analogously
  from the same pooled set (see remark 4).
- Monte-Carlo p-value: `p = (B+1)^{-1} sum_{j=0}^B 1{Q_j >= Q_0}`.

## Claim

**Proposition.** Under `H0` and the block-randomization design above, for any
`B >= 1`, any measurable learner `L`, and any representations `f_c`, the
Monte-Carlo p-value `p` is exactly level-valid: `P(p <= a) <= a` for all
`a in (0,1)`, where `P` is the randomization law.

The same holds for the max and Cauchy combinations, and for the screen
(remark 4).

## Assumptions, by role

- (A1) Sharp null. LOAD-BEARING. Makes `Y` fixed so that `T(z)` is a function of
  `z` alone. Under a non-sharp null `Y` is random and the argument fails at
  step 1.
- (A2) Refit per permutation. LOAD-BEARING. `e(z)` (hence `T(z)`) must be
  recomputed at each `z`. If instead the residuals are fixed at `Z0`
  (`e = Y - L(X, Y, {i: Z0_i = 0})`), then `T(z) = sum_i z_i s_{c,i}` with `s`
  depending on `Z0`; `T` is then a function of `(z, Z0)`, not of `z` alone, and
  step 1 fails. Failure mode: anti-conservative when `L` overfits (verified).
- (A3) Symmetric (pooled) metric AND symmetric column selection. LOAD-BEARING.
  Everything computed from the statistics --- the metric `mu_hat, Sigma_hat` AND
  the choice of which representations to keep (drop a representation only when its
  linear statistic is constant across the pooled orbit) --- must be a symmetric
  function of all `B + 1` statistics INCLUDING the observed. The draws-only metric
  (excluding `Z0`) breaks symmetry at step 4 (verified: size 0.116 at `B = 49`,
  `C = 6`). Selecting the kept representations from the OBSERVED assignment alone
  breaks it the same way: the surviving column set then depends on `Z0`, so the map
  `(Z0,...,ZB) -> (Q_0,...,QB)` is not symmetric. riposte selects symmetrically
  (it drops a representation only on zero pooled variance across the orbit).
- (A4) Measurability of `L` and `f_c`. TECHNICAL. Needed only so that `T(.)` is a
  well-defined measurable map; no other regularity is used. This is the sense in
  which exactness holds for "any learner."

## Proof

We show `(Q_0, ..., Q_B)` is exchangeable; level validity then follows from the
uniform rank of `Q_0`.

1. Determinism under the null. By (A1), `Y` is fixed; with `X` fixed, (A2) makes
   `e(.)`, `s_c(.)`, and hence `T(.)` deterministic measurable functions of the
   assignment: `T(z) = T(z; Y, X)`. Write the fixed map `phi(z) = T(z)`.

2. Exchangeability of the assignments. `Z0` is, under the design, a draw from `P`;
   `Z1, ..., ZB` are independent draws from the same `P`. Hence
   `(Z0, Z1, ..., ZB)` is an i.i.d. (so exchangeable) sequence.

3. Exchangeability of the statistics. Applying the fixed map `phi` coordinatewise,
   `(phi(Z0), ..., phi(ZB)) = (T(Z0), ..., T(ZB))` is exchangeable (a fixed
   function of an exchangeable sequence is exchangeable).

4. Permutation-equivariance of the combination. Take `T(z)` to be the full vector
   of all `C` representations. Two things are computed from the statistics, and
   both depend on `(T(Z0), ..., T(ZB))` only through the unordered multiset `M`:
   (i) the kept-column set `K(M)` (drop coordinate `c` iff its values are constant
   across `M`, i.e. zero pooled variance), and (ii) the pooled mean and covariance
   `mu(M), Sigma(M)` of the kept coordinates --- means, variances, and a
   constant-across-the-multiset test are all symmetric in their arguments. Define
   `g(t; M) = (t_{K(M)} - mu(M))' Sigma(M)^+ (t_{K(M)} - mu(M))`. Then
   `Q_j = g(T(Zj); M)` with the SAME `M` for every `j`. The map
   `(t_0, ..., t_B) -> (g(t_0; M), ..., g(t_B; M))`, `M = {t_0, ..., t_B}`, is
   symmetric: permuting the inputs leaves `M` --- hence `K`, `mu`, `Sigma` ---
   unchanged and merely relabels the `t_j`. (Selecting `K` from the observed `Z0`
   alone would make `K` depend on the labelling and break this symmetry; riposte
   selects from the pooled orbit, which does not.)

5. Exchangeability of the combined statistics. A symmetric function of an
   exchangeable vector is exchangeable, so `(Q_0, ..., Q_B)` is exchangeable.

6. Uniform rank and level. By exchangeability the rank of `Q_0` among
   `(Q_0, ..., Q_B)` (ties broken by the `>=` convention, which is conservative)
   is stochastically no larger than uniform on `{1, ..., B+1}`. Since
   `p = (B+1)^{-1} #{j : Q_j >= Q_0}`, for any `a`,
   `P(p <= a) = P(#{j : Q_j >= Q_0} <= a(B+1)) <= a`. This is the standard
   Monte-Carlo permutation-test bound with the observed assignment included
   (Phipson and Smyth 2010). QED.

## Remarks (scope and the two open items)

1. Finite-sample, design-based. The guarantee is exact in finite samples under the
   assignment mechanism; no asymptotics, no consistency of `L`, no smoothness.
   This is the Rosenbaum guarantee that "the test statistic is the whole
   algorithm" (the refit included): a deterministic function of `(Z, Y, X)` has an
   exact randomization distribution. The novelty is small and worth stating
   plainly: exactness for any measurable learner is Rosenbaum's; here we only add
   that the COMBINED statistic with a POOLED metric inherits it, and that the
   screen does too.

2. "Exact" means exact LEVEL, not the enumeration value. `p` is a Monte-Carlo
   permutation p-value: valid for any `B`, and equal to the full-enumeration
   p-value only as `B -> inf`. The pooled `Sigma_hat` is a (consistent-as-`B`)
   estimate of the assignment-induced covariance of `T`; its sampling noise
   affects POWER, not level. So open item #1 is resolved not by deriving the
   closed-form adjusted moments (which would only help power) but by showing the
   estimated, pooled moments already give exact level.

3. The unadjusted case is a special case with a better metric. There `T` is linear
   in `z` and the closed-form Strasser-Weber `mu, Sigma` are the exact permutation
   moments --- constants, trivially symmetric --- so the same step 4--6 argument
   applies and, additionally, the metric carries no Monte-Carlo noise. That is why
   riposte keeps the closed form for the unadjusted quadratic and the screen.

4. Other combinations and the screen (resolves open item #2). The max combination
   uses `g(t; M) = max_c |t_c - mu_c(M)| / sigma_c(M)` with `sigma_c` from the
   pooled `Sigma`; symmetric in `M`, so step 4 holds. The Cauchy uses
   `g(t; M) = mean_c ACAT(midp_c(t; M))` where `midp_c(t; M)` is the mid-rank of
   `t_c` in the pooled column `{t_{0,c}, ..., t_{B,c}}` --- again symmetric (the
   implementation already ranks the observed within the full column, so the Cauchy
   was already exact). The SCREEN selects the combination, or the shrinkage
   intensity `lambda`, from `kappa(Sigma_hat)` (and produces a shrunk
   `Sigma_hat(lambda)`); all are symmetric functions of `M`. A symmetric selection
   composed with symmetric combinations is still symmetric, so step 5 holds and
   the adjusted screen is exact. Under adjustment the condition number is no longer
   ancillary in the closed-form sense (it depends on the realized draws), but it is
   a symmetric function of the exchangeable pooled set, which is what the proof
   needs.

5. Still open. The shrinkage rule for the screen (how `lambda` should depend on the
   conditioning) is a power/efficiency question, not a level question, and remains
   the Project 1 research deliverable. Nothing above depends on the choice of
   `lambda`: exact level holds for any symmetric `lambda` rule.

## Simulation evidence

Design-based size under `H0` (no effect; covariates unrelated to `Y`),
block-randomized, ridge learner refit per draw, six default representations.

- Pooled vs draws-only metric, `B = 49`, `C = 6`, 1500 simulations (size SE about
  0.006): quadratic size 0.116 (draws-only) vs 0.045 (pooled); max 0.063
  (draws-only) vs 0.043 (pooled). The pooled metric removes the anti-conservatism;
  draws-only is badly anti-conservative at small `B`.
- Full enumeration, two blocks of four with two treated each (36 assignments),
  ridge refit per assignment, pooled quadratic, 3000 simulations: the p-value is
  uniform under `H0` --- `P(p <= k/36)` matched `k/36` at `k = 1, 3, 6, 9, 18`
  (e.g. 0.166 vs 0.167 at `k = 6`; mean p 0.518 vs uniform 0.514). This is the
  exact-level claim with no Monte-Carlo metric noise in the reference.

Both checks are reproduced as `skip_on_cran` tests in
`tests/testthat/test-adjust-validity.R`.

## What changed in the code

`riposte_adjusted_combination()` (`R/adjust.R`) now computes `mu, Sigma` from the
full `Tmat` (observed + draws) instead of `Tmat[-1, ]` (draws only). One line; the
proof above is the justification.
