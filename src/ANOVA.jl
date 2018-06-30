module ANOVA

using GLM, DataFrames, Distributions, CategoricalArrays, ArgCheck, StatsModels
import StatsBase

struct AnovaObject
  Source::Array{String,1}
  DF::Array{AbstractFloat,1}
  SS::Array{AbstractFloat,1}
  MSS::Array{AbstractFloat,1}
  F::Array{AbstractFloat,1}
  p::Array{AbstractFloat,1}
end


struct AnovaDataFrameRegressionModel{M} <: StatsModels.RegressionModel
  model::M
  mf::ModelFrame
  mm::ModelMatrix
  anova::AnovaObject
end

struct AnovaOptions
  anovatype::Int
end


# function StatsBase.fit(::Type{LinearModel}, f::Formula, df::AbstractDataFrame, options::AnovaOptions,args...; contrasts::Dict = Dict(), kwargs...)

#   @argcheck options.anovatype == 1 || options.anovatype == 3 "'anovatype' in 'options' must be either '1' or '3'"

#   trms = StatsModels.Terms(f)
#   StatsModels.drop_intercept(LinearModel) && (trms.intercept = true)
#   mf = ModelFrame(trms, df, contrasts=contrasts)
#   StatsModels.drop_intercept(LinearModel) && (mf.terms.intercept = false)
#   mm = ModelMatrix(mf)
#   y = StatsModels.model_response(mf)

#   model = fit(LinearModel, mm.m, y, args...; kwargs...)
#   AnovaDataFrameRegressionModel{LinearModel}(model, mf, mm, anova(model, mf, mm, anovatype = options.anovatype))
# end

## helper functions ##

# calculate effects for type I ANOVA (I stole this somewhere, but don't remember where :( )
function effects(mod)
    return (mod.pp.X / cholfact!(mod.pp)[:U])' * mod.rr.y
  end

# this will drop variables from the matrix using the mask argument and fit a new linear model
function droptermbymask(mod, mask)
    matrix = mod.pp.X[:,.!mask]
    y = mod.rr.y
    return lm(matrix,y)
end

# calculate ANOVA
function anova(mod::LinearModel, mf::ModelFrame, mm::ModelMatrix; anovatype = 3)

  eff = effects(mod) # get effects 

  T = eltype(mod.pp.X) # get types
  response = mf.terms.eterms[1] 
  terms = mf.terms.terms

  ## calculate variables for residuals (this is the full modell) ##
  DFres = dof_residual(mod)           # degrees of freedom
  #SSres = sum(abs2.(residuals(mod))) # used deviance() instead
  SSres = deviance(mod)               # calculate sum of squares of weighted residuals
  MSSres = SSres / DFres              # calculate mean sum of squares for residuals

  # check if there is an intercept
  v = all(mod.pp.X[:,1] .== 1.0) ? 1 : 0

  # get assigned numbers for factors
  factors_assig = unique(mm.assign)

  ## calculate some variables for computation ##
  n = length(predict(mod))              # sample size
  Nfactors = length(factors_assig) - v  # number of factors

  ## create arrays for ANOVA factors ##
  DF = zeros(T,Nfactors)   # will contain degrees of freedoms for factor
  SS = zeros(T,Nfactors)   # will contain sum of squares for factor
  MSS = zeros(T,Nfactors)  # will contain mean sum of squares for factor
  RSS = zeros(T,Nfactors)  # 
  F = zeros(T,Nfactors)    # will contain F values for factors
  p = zeros(T,Nfactors)    # will contain pvalue for factors


for i in 1:Nfactors
  mask = mm.assign .== factors_assig[i+v]
  if anovatype == 3 # if this is a type III ANOVA, we use a similar procedure to R's drop1() function
    print(mask)
    newmod = droptermbymask(mod, mask)
    RSS[i] = deviance(newmod)
    SS[i] = RSS[i] - deviance(mod)
    DF[i] = sum(mask)
    MSS[i] = SS[i]/DF[i]
    F[i] = ftest(mod,newmod).fstat[1]
  else # type I ANOVA
    SS[i] = sum(abs2.(eff[mask]))
    RSS[i] = 0
    DF[i] = sum(mask)
    MSS[i] = SS[i]/DF[i]
    F[i] = MSS[i]/MSSres
  end

  p[i] = ccdf(FDist(DF[i], DFres), F[i])
end

AnovaObject( 
          vcat(string.(terms),"Residuals"),
          vcat(DF,DFres),
          vcat(SS, SSres),
          vcat(MSS, MSSres),
          vcat(F, 0),
          vcat(p, 0)
          )
end

anova(mod::StatsModels.DataFrameRegressionModel; anovatype = 3) = anova(mod.model, mod.mf, mod.mm, anovatype = anovatype)



StatsBase.predict(m::AnovaDataFrameRegressionModel) = predict(m.model)
StatsBase.residuals(m::AnovaDataFrameRegressionModel) = residuals(m.model)

function Base.show(io::IO, m::AnovaDataFrameRegressionModel)
  println("ANOVA table for linear model")
  print(m.anova)
end

function Base.show(io::IO, x::AnovaObject)
  print(DataFrame( 
          Source = x.Source,
          DF = x.DF,
          SS = x.SS,
          MSS = x.MSS,
          F = x.F,
          p = x.p
          ))
end

export AnovaDataFrameRegressionModel, AnovaOptions, anova, AnovaObject

end # module
