# Analysis of Variance in Julia
[![Build Status](https://travis-ci.org/marcpabst/ANOVA.jl.svg?branch=master)](https://travis-ci.org/marcpabst/ANOVA.jl)
[![Coverage Status](https://coveralls.io/repos/github/marcpabst/ANOVA.jl/badge.svg?x=b&branch=master)](https://coveralls.io/github/marcpabst/ANOVA.jl?branch=master)


Calculate ANOVA tables for linear models. Currently supports type I and type III.

Minimal Example:

```julia
using ANOVA
using GLM
using RDatasets

data = dataset("datasets", "ToothGrowth")
data[:Dose] = categorical(data[:Dose])

model = fit(LinearModel,
            @formula(Len ~  Supp + Dose), 
            data, 
            contrasts = Dict(:Supp => EffectsCoding(),:Dose => EffectsCoding()))
anova(model)
 ```
 Output
 ```
ANOVA table for linear model
3×6 DataFrames.DataFrame
│ Row │ Source    │ DF   │ SS      │ MSS     │ F       │ p           │
├─────┼───────────┼──────┼─────────┼─────────┼─────────┼─────────────┤
│ 1   │ Supp      │ 1.0  │ 205.35  │ 205.35  │ 14.0166 │ 0.000429279 │
│ 2   │ Dose      │ 2.0  │ 2426.43 │ 1213.22 │ 82.8109 │ 1.87116e-17 │
│ 3   │ Residuals │ 56.0 │ 820.425 │ 14.6504 │ -1.0    │ -1.0        │
```

