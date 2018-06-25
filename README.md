# Analysis of Variance in Julia
[![Build Status](https://travis-ci.org/marcpabst/ANOVA.jl.svg?branch=master)](https://travis-ci.org/marcpabst/ANOVA.jl)
Calculate ANOVA tables for linear models. Currently supports type I and type III.

Usage:

```julia
anova(lm(@formula(Y ~ A + B), data))
```

