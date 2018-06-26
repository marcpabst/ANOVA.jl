using ANOVA
using RDatasets
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# # write your own tests here

data = dataset("datasets", "ToothGrowth")
data[:Dose] = categorical(data[:Dose])

model = fit(LinearModel,
            @formula(Len ~  Supp + Dose), 
            data, 
            contrasts = Dict(:Supp => EffectsCoding(),:Dose => EffectsCoding()))

@test ANOVA.anova(model).p[1] ≈ 0.00042927927678531035
@test ANOVA.anova(model).p[2] ≈ 1.8711626362930449e-17
