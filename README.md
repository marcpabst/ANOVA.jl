# Analysis of Variance in Julia
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Build Status](https://travis-ci.org/marcpabst/ANOVA.jl.svg?branch=master)](https://travis-ci.org/marcpabst/ANOVA.jl)
[![Coverage Status](https://coveralls.io/repos/github/marcpabst/ANOVA.jl/badge.svg?x=b&branch=master)](https://coveralls.io/github/marcpabst/ANOVA.jl?branch=master)


Calculate ANOVA tables for linear models. If no `anovatype` argument is provided, a type III ANOVA will be calculated. Type I and II are also supported for compatibility. Support for mixed models and a more convenient way to create ANOVAs (similar to the `ez` package in R) is planned.

Important: Make sure to use `EffectsCoding` on all your predictors, or results won't be meaningful.

Minimal Example:

```julia
using RDatasets, GLM, ANOVA

data = dataset("datasets", "ToothGrowth")
categorical!(data, :Dose)

model = fit(LinearModel,
            @formula(Len ~  Supp + Dose),
            data,
            contrasts = Dict(:Supp => EffectsCoding(), :Dose => EffectsCoding()))
Anova(model)
ANOVA table for linear model
             DF      SS     MSS       F      p
Supp        1.0  205.35  205.35 14.0166 0.0004
Dose        2.0 2426.43 1213.22 82.8109 <1e-16
Residuals  56.0 820.425 14.6504     0.0 <1e-99
```
