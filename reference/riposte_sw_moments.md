# Closed-form Strasser-Weber permutation mean and covariance

Computes the exact permutation mean vector and covariance matrix of the
representation linear statistics `T_c = sum_i z_i * s_{c,i}`, from the
scores and the design, without re-randomizing. The permutation fixes
each block's treated count, so within block `z` is a simple random
sample of size `m_b` from `n_b` units. For that sampling, with
`p_b = m_b / n_b`,

- `E[T_c] = sum_b p_b * sum_{i in b} s_{c,i}`, and

- `Cov(T_a, T_b) = sum_b w_b * sum_{i in b} (s_{a,i} - sbar_a)(s_{b,i} - sbar_b)`,
  where `w_b = m_b (n_b - m_b) / (n_b (n_b - 1))` and `sbar` is the
  within-block mean.

The covariance uses within-block-centred cross products, so it is the
same whether `scores` are already centred (as from
[`riposte_score_matrix()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_score_matrix.md))
or not; the function centres internally. When `scores` are centred, the
mean is zero, as riposte's combinations assume.

## Usage

``` r
riposte_sw_moments(scores, z, block)
```

## Arguments

- scores:

  numeric matrix of scores (rows = units, columns = representations),
  e.g. `riposte_score_matrix()$scores`.

- z:

  0/1 treatment vector (its per-block sums give the treated counts).

- block:

  factor of block labels.

## Value

a list with `mu` (length = number of representations) and `Sigma`
(representation-by-representation covariance matrix).

## Details

Verified against `coin`'s `expectation()` and `covariance()` (the same
linear-statistic framework) to machine precision; see the tests.
