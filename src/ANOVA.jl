module ANOVA

using GLM, DataFrames, Distributions, CategoricalArrays, ArgCheck, StatsModels, LinearAlgebra
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

## helper functions ##

function merge_bool_array(a, b)
  return ifelse.(b,b,a)
end

# calculate effects for type I ANOVA
function effects(mod)
    return (mod.pp.X / cholesky!(mod.pp).U)' * mod.rr.y
end

# this will drop variables from the matrix using the mask argument then fit a new linear model
function droptermbymask(mod, mask)
  matrix = mod.pp.X[:,.!mask]
  y = mod.rr.y
  return lm(matrix,y)
end


# calculate ANOVA

"""
    anova(mod, mf, mm[, anovatype = 3])

Compute an ANOVA table for the given linear model.
"""
function anova(mod::LinearModel, mf::ModelFrame, mm::ModelMatrix; anovatype = 2)

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
    fullmod = mod
    newmod = droptermbymask(fullmod, mask)
    RSS[i] = deviance(newmod)
    SS[i] = RSS[i] - deviance(fullmod)
    DF[i] = sum(mask)
    MSS[i] = SS[i]/DF[i]
    F[i] = ftest(fullmod,newmod).fstat[1]
  elseif anovatype == 2 # if this is a type II ANOVA
    full_model_mask = falses(length(mm.assign))
    for j in 1:Nfactors
      term = terms[j]
      if(typeof(term) == Expr && typeof(terms[i])==Symbol) # this is not an interaction term
        if(terms[i] in term.args)
          full_model_mask = merge_bool_array(full_model_mask, mm.assign .== factors_assig[j+v])
        end      
       elseif(typeof(term) == Expr && typeof(terms[i])==Expr) # this is an interaction term
        if(all(in.(terms[i].args, (term.args,) )) && term != terms[i])
          full_model_mask = merge_bool_array(full_model_mask, mm.assign .== factors_assig[j+v])
        end
      end
    end
    fullmod = droptermbymask(mod, full_model_mask) # we have to drop all interaction terms containing the current term in the full modell
    newmod = droptermbymask(mod, merge_bool_array(full_model_mask, mask))
    RSS[i] = deviance(newmod)
    SS[i] = RSS[i] - deviance(fullmod)
    DF[i] = sum(mask)
    MSS[i] = SS[i]/DF[i]
    F[i] = MSS[i]/MSSres
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

function Base.isapprox(a::ANOVA.AnovaObject, b::ANOVA.AnovaObject)
  all( [a.Source == b.Source,
  all(a.DF .≈ b.DF),
  all(a.SS .≈ b.SS),
  all(a.MSS .≈ b.MSS),
  all(a.F .≈ b.F),
  all( a.p .≈ b.p)])
end

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
