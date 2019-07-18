module ANOVA
import GLM: ccdf, cholesky!, CoefTable, deviance, dof_residual, FDist, fit,
            ftest, LinearModel, nobs, drop_term, StatsModels.hasintercept, term,
            width

import Base: show, summary

"""
    Anova(mod, anovatype = 3)

Compute an ANOVA table for the given linear model.
"""
struct Anova
  Source::Vector{String}
  DF::Vector{Float64}
  SS::Vector{Float64}
  MSS::Vector{Float64}
  F::Vector{Float64}
  p::Vector{Float64}
  function Anova(mod; anovatype = 3)
    mf = mod.mf
    mm = mod.mm
    mod = mod.model

    eff = effects(mod) # get effects
    response = mf.f.lhs
    terms = drop_term(mf.f.rhs, term(1))

    ## calculate variables for residuals (this is the full modell) ##
    DFres = dof_residual(mod)           # degrees of freedom
    SSres = deviance(mod)               # calculate sum of squares of weighted residuals
    MSSres = SSres / DFres              # calculate mean sum of squares for residuals

    # check if there is an intercept
    v = hasintercept(mf.f)

    # Temporary hack to fix https://github.com/JuliaStats/StatsModels.jl/issues/133
    mm.assign .= mapreduce((it)->(first(it)-1)*ones(width(last(it))), append!, enumerate(mf.f.rhs.terms), init=Int[])
    # get assigned numbers for factors
    factors_assig = unique(mm.assign)

    ## calculate some variables for computation ##
    n = nobs(mod)                                # sample size
    Nfactors = count(x -> x > 0, factors_assig)  # number of factors

    ## create arrays for ANOVA factors ##
    Source = fill("", Nfactors + 1)
    DF = zeros(Nfactors + 1)
    SS = zeros(Nfactors + 1)
    MSS = zeros(Nfactors + 1)
    RSS = zeros(Nfactors + 1)
    F = zeros(Nfactors + 1)
    p = zeros(Nfactors + 1)

    for i in 1:Nfactors
      Source[i] = string(terms.terms[i])
      mask = mm.assign .== factors_assig[i+v]
      if anovatype == 3 # if this is a type III ANOVA, we use a similar procedure to R's drop1() function
        fullmod = mod
        newmod = droptermbymask(fullmod, mask)
        RSS[i] = deviance(newmod)
        SS[i] = RSS[i] - deviance(fullmod)
        DF[i] = count(mask)
        MSS[i] = SS[i]/DF[i]
        F[i] = ftest(fullmod,newmod).fstat[1]
      elseif anovatype == 2 # if this is a type II ANOVA
        full_model_mask = falses(length(mm.assign))
        for j in 1:Nfactors
          term = terms[j]
          if (isa(term, Expr) && isa(terms[i], Symbol)) # this is not an interaction term
            if terms[i] ∈ term.args
              full_model_mask = merge_bool_array(full_model_mask, mm.assign .== factors_assig[j+v])
            end
          elseif (isa(term, Expr) && isa(terms[i], Expr)) # this is an interaction term
            if (all(x -> x ∈ term.args, terms[i].args) && term ≠ terms[i])
              full_model_mask = merge_bool_array(full_model_mask, mm.assign .== factors_assig[j+v])
            end
          end
        end
        fullmod = droptermbymask(mod, full_model_mask) # we have to drop all interaction terms containing the current term in the full modell
        newmod = droptermbymask(mod, merge_bool_array(full_model_mask, mask))
        RSS[i] = deviance(newmod)
        SS[i] = RSS[i] - deviance(fullmod)
        DF[i] = count(mask)
        MSS[i] = SS[i] / DF[i]
        F[i] = MSS[i] / MSSres
      else # type I ANOVA
        SS[i] = sum(abs2, eff[mask])
        RSS[i] = 0
        DF[i] = count(mask)
        MSS[i] = SS[i] / DF[i]
        F[i] = MSS[i] / MSSres
      end
      p[i] = ccdf(FDist(DF[i], DFres), F[i])
    end
    Source[end] = "Residuals"
    DF[end] = DFres
    SS[end] = SSres
    MSS[end] = MSSres
    F[end] = NaN
    p[end] = NaN
    new(Source, DF, SS, MSS, F, p)
  end
  Anova(Source::AbstractVector{<:AbstractString},
        DF::AbstractVector{<:Real},
        SS::AbstractVector{<:Real},
        F::AbstractVector{<:Real},
        p::AbstractVector{<:Real}) = new(Source, DF, SS, SS ./ DF, F, p)
end

## helper functions ##

merge_bool_array(a, b) = ifelse.(b, b, a)

# calculate effects for type I ANOVA
effects(mod) = (mod.pp.X / cholesky!(mod.pp).U)' * mod.rr.y

# this will drop variables from the matrix using the mask argument then fit a new linear model
droptermbymask(mod, mask) = fit(LinearModel, mod.pp.X[:,.!mask], mod.rr.y)

summary(io::IO, obj::Anova) =
  println(io, "ANOVA table for linear model")

show(io::IO, obj::Anova) =
  println(io,
          summary(obj),
          CoefTable(hcat(obj.DF, obj.SS, obj.MSS, obj.F, obj.p),
                    ["DF", "SS", "MSS", "F", "p"],
                    obj.Source,
                    5))

export Anova

end # module
