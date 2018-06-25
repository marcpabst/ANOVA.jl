# Analysis of Variance in Julia
Calculate ANOVA tables for linear models. Currently supports type I and type III.

Usage:

```julia
anova(lm(@formula(Y ~ A + B), data))
```

