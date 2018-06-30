import ANOVA
using DataFrames
using GLM

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# # write your own tests here
data = convert(DataFrame,[
"1" "GT" "10" "enabled" 19.3557287740925;
"2" "GT" "1.5" "disabled" 23.3806943448596;
"3" "MX" "1.5" "enabled" 24.0077846814578;
"4" "MX" "10" "enabled" 25.0714195337438;
"5" "LS" "10" "disabled" 26.3983294771175;
"6" "GT" "10" "enabled" 18.6088838586582;
"7" "LS" "1.5" "enabled" 23.2992470891206;
"8" "LS" "1.5" "disabled" 26.1332521553406;
"9" "LS" "10" "disabled" 31.024110406057;
"10" "MX" "10" "disabled" 26.175152381339;
"11" "MX" "10" "disabled" 27.0498523179508;
"12" "GT" "10" "disabled" 24.6421568043518;
"13" "LS" "1.5" "enabled" 23.8080587240068;
"14" "GT" "1.5" "disabled" 25.9001814923499;
"15" "MX" "10" "enabled" 24.0457318442931;
"16" "MX" "1.5" "enabled" 24.8136216797396;
"17" "GT" "10" "disabled" 26.0228457322805;
"18" "MX" "1.5" "disabled" 27.2711342214877;
"19" "LS" "1.5" "disabled" 29.4439325048592;
"20" "LS" "10" "enabled" 24.0498649478329;
"21" "GT" "1.5" "enabled" 21.9492468934002;
"22" "MX" "1.5" "disabled" 26.2485867993077;
"23" "GT" "1.5" "enabled" 22.4742018940028;
"24" "LS" "10" "enabled" 21.0579535297715;
"25" "LS" "10" "enabled" 20.0579535297715])

names!(data, [:ID, :Tire, :Tread, :ABS, :Distance])

data[:Tread] = categorical(data[:Tread])
data[:Tire] = categorical(data[:Tire])
data[:ABS] = categorical(data[:ABS])
data[:Distance] = convert(Array{Float64,1}, data[:Distance])

# prepare test for ANOVA type III

anova3_expected = ANOVA.AnovaObject( 
          vcat("ABS", "Tread", "Tire","ABS & Tread", "ABS & Tire", "Tread & Tire", "&(ABS, Tread, Tire)", "Residuals"),
          vcat(1 , 1, 2, 1, 2, 2, 2, 13),
          vcat(102.1138211805702980, 2.0356309732669118, 40.1164208238708397, 6.9792255733881632, 12.3635146287293978, 1.6557631208168218, 4.7578414794463484, 31.237427095361188),
          vcat(102.1138211805702980, 2.0356309732669118, 40.1164208238708397, 6.9792255733881632, 12.3635146287293978, 1.6557631208168218, 4.7578414794463484, 31.237427095361188)./vcat(1 , 1, 2, 1, 2, 2, 2, 13),
          vcat(42.49644733207060199, 0.84716332659802784, 8.34757403543913412, 2.90452642520995008, 2.57264610306765462, 0.34453734785690932, 0.99002934915192542, 0),
          vcat(1.9453276992336526e-05, 3.7412049964586458e-01, 4.6577488940770014e-03, 1.1209905804371272e-01, 1.1446323312399301e-01, 7.1482574202441140e-01, 3.9791832675247707e-01, 0)
          )
model3 = fit(LinearModel,
            @formula(Distance ~ ABS * Tread * Tire), 
            data, 
            contrasts = Dict(:ABS => EffectsCoding(),:Tread => EffectsCoding(), :Tire => EffectsCoding()))
anova3_actual  = ANOVA.anova(model3, anovatype = 3)


# prepare test for ANOVA type I
anova1_expected = ANOVA.AnovaObject( 
          vcat("ABS", "Tread", "Tire","ABS & Tread", "ABS & Tire", "Tread & Tire", "&(ABS, Tread, Tire)", "Residuals"),
          vcat(1 , 1, 2, 1, 2, 2, 2, 13),
          vcat(106.5991986929681445,   2.3651163588824207,  38.2807197340170120,   7.8051389666327555,  12.5399359441607761,   1.6429152697879308, 4.7578414794463368,31.2374270953612125),
          vcat(106.5991986929681445,   2.3651163588824207,  38.2807197340170120,   7.8051389666327555,  12.5399359441607761,   1.6429152697879308, 4.7578414794463368,31.2374270953612125)./vcat(1 , 1, 2, 1, 2, 2, 2, 13),
          vcat(44.36311539929538128,  0.98428441534601785,  7.96559452580719363,  3.24824468598099525,  2.60935650648222861,  0.34186391923448084,  0.99002934915192564, 0),
          vcat(1.5624161055791751e-05, 3.3925131671067088e-01, 5.5176093522618852e-03, 9.4722296495528976e-02, 1.1149792216959341e-01, 7.1664324192550988e-01, 3.9791832675247707e-01, 0)
          )
model1 = fit(LinearModel,
            @formula(Distance ~ ABS * Tread * Tire), 
            data, 
            contrasts = Dict(:ABS => EffectsCoding(),:Tread => EffectsCoding(), :Tire => EffectsCoding()))
anova1_actual  = ANOVA.anova(model1, anovatype = 1)

@test anova3_actual ≈ anova3_expected
@test anova1_actual ≈ anova1_expected

function Base.isapprox(a::ANOVA.AnovaObject, b::ANOVA.AnovaObject)
    all( [a.Source == b.Source,
    all(a.DF .≈ b.DF),
    all(a.SS .≈ b.SS),
    all(a.MSS .≈ b.MSS),
    all(a.F .≈ b.F),
    all( a.p .≈ b.p)])
end
