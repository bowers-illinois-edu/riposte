# The point of riposte

riposte is an R package for researchers who analyze randomized
experiments — completely randomized, block-randomized, or
cluster-randomized — and want to test the sharp null hypothesis that
treatment changed no one’s outcome. A difference of means can miss an
effect that changes the spread of the outcomes or changes the outcomes
of only a few units. riposte combines several representations of the
outcome (raw values, ranks, distances, a bounded transform) into one
test with power against such effects. It compares the observed test
statistic with its distribution over re-randomizations of treatment, so
the chance of rejecting a true sharp null hypothesis is at most the
level the researcher chose, in samples of any size. When a user picks a
large-sample approximation instead, the printed result says that it is
an approximation.

For matched pairs, riposte also offers Rosenbaum’s sensitivity analysis
for observational studies. The researcher chooses gamma: an unmeasured
difference between the two units of a pair could have made one unit’s
odds of treatment as much as gamma times the other’s. riposte reports an
upper bound on the p-value, the largest p-value that hidden bias no
larger than gamma could produce. At gamma = 1 that bound is the p-value
of a paired randomized experiment. Under the sharp null hypothesis only
the signs of the treated-minus-control differences depend on which unit
was treated, so riposte computes the bound from the chance of each
pattern of signs, by listing every pattern or by drawing patterns at
random. It uses a normal approximation only when the user asks for one.

The package is the software for Bowers and Burton, “A More Powerful Test
for Randomized Experiments”: a reader of the paper can run its tests on
their own experiment.

riposte tests hypotheses. It does not estimate the size of effects
(propertee), locate effects among blocks (manytestsr), or count the
units harmed (CMRSS).
