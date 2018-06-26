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

            (Intercept) 2326.91  1 158.828 < 2.2e-16 ***
supp         205.35  1  14.017 0.0004293 ***
dose        2426.43  2  82.811 < 2.2e-16 ***
Residuals    820.43 56  

@test ANOVA.anova(model).p[1] ≈ 0.00042927927678531035
@test ANOVA.anova(model).p[2] ≈ 1.8711626362930449e-17
