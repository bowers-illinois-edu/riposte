# Truncated Cauchy combination of p-values

Combines p-values with the truncated Cauchy conversion of Gui, Jiang,
and Wang (2025). Each p-value \\p\\ is converted to \\\tan((0.5 - t
p)\pi)\\, where \\t\\ is `truncation`; with \\n\\ p-values whose
converted values sum to \\S\\, the combined p-value is \\\min(1, n P(C
\> S) / t)\\, with \\C\\ a standard Cauchy variable.

## Usage

``` r
riposte_truncated_cauchy(p, truncation = 0.9)
```

## Arguments

- p:

  numeric vector of p-values in \[0, 1\].

- truncation:

  the share \\t\\ of the Cauchy distribution kept, in (0, 1\]; 0.9 by
  default.

## Value

the combined p-value.

## Details

Liu and Xie's Cauchy combination converts \\p\\ to \\\tan((0.5 -
p)\pi)\\, which is minus infinity at \\p = 1\\: one p-value of 1 then
makes the combined p-value 1 whatever the others are, and a p-value near
1 does almost as much (0.999 converts to about -318). On outcomes with
few distinct values a representation's treated sum can equal its mean
exactly, and its large-sample p-value is then 1. With \\t = 0.9\\, the
value Gui, Jiang, and Wang recommend, \\p = 1\\ converts to
\\\tan(-0.4\pi) = -3.08\\. `truncation = 1` gives Liu and Xie's
untruncated combination, referred to the standard Cauchy as
[`riposte_test()`](https://bowers-illinois-edu.github.io/riposte/reference/riposte_test.md)
with `engine = "asymptotic"` did before riposte 0.0.0.9008.

## References

Gui, L., Jiang, Y., and Wang, J. (2025). Aggregating dependent signals
with heavy-tailed combination tests. *Biometrika*, 112(4), asaf038.
[doi:10.1093/biomet/asaf038](https://doi.org/10.1093/biomet/asaf038)

Fang, Y., Chang, C., Park, Y., and Tseng, G. C. (2023). Heavy-tailed
distribution for combining dependent p-values with asymptotic
robustness. *Statistica Sinica*, 33, 1115-1142.
[doi:10.5705/ss.202022.0046](https://doi.org/10.5705/ss.202022.0046)

## Examples

``` r
# four small p-values and one representation that finds no difference
p <- c(0.01, 0.01, 0.01, 0.01, 1)
riposte_truncated_cauchy(p)                  # 0.013
#> [1] 0.01278124
riposte_truncated_cauchy(p, truncation = 1)  # Liu and Xie's version: 1
#> [1] 1
```
