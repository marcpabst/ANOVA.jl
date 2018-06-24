module ANOVA

using GLM, DataFrames, Distributions, CategoricalArrays


## helper functions
function effects(mod::LinearModel)
    return (mod.pp.X / cholfact!(mod.pp)[:U])' * mod.rr.y
  end

function droptermbymask(mod, mask)
    matrix = mod.model.pp.X[:,.!mask]
    y = mod.model.rr.y
    return lm(matrix,y)
end


function anova(mod; anovatype = 3)
    

  eff = effects(mod.model) # get effect sizes 

  T = eltype(mod.model.pp.X)
  response = mod.mf.terms.eterms[1]
  terms = mod.mf.terms.terms
  assign = mod.mm.assign

  ## calculate variables for residuals (this is the full / saturated modell) ##
  DFres = dof_residual(mod)           # degrees of freedom
  #SSres = sum(abs2.(residuals(mod))) # used deviance() instead
  SSres = deviance(mod)               # calculate sum of squares of weighted residuals
  MSSres = SSres / DFres              # calculate mean sum of squares for residuals

  # check if there is an intercept
  k = all(mod.model.pp.X[:,1] .== 1.0) ? 1 : 0

  # get assigned numbers for factors
  factors_assig = unique(assign)

  ## calculate some variables for computation ##
  n = length(predict(mod))              # sample size
  Ncomb = length(assign)                # number of combinations
  Nfactors = length(factors_assig) - k  # number of factors

  ## create arrays for ANOVA factors ##
  DF = zeros(T,Nfactors)   # will contain degrees of freedoms for factor
  SS = zeros(T,Nfactors)   # will contain sum of squares for factor
  MSS = zeros(T,Nfactors)  # will contain mean sum of squares for factor
  RSS = zeros(T,Nfactors)  # 
  F = zeros(T,Nfactors)    # will contain F values for factors
  pval = zeros(T,Nfactors) # will contain pvalue for factors


for i in 1:Nfactors
  v = factors_assig[i+k]
  mask = assign .== v
  if anovatype == 3 # if this is a type III ANOVA, we use a similar procedure to R's drop1() function
    newmod = droptermbymask(mod, mask)
    RSS[i] = deviance(newmod)
    SS[i] = RSS[i] - deviance(mod)

    DF[i] = sum(mask)
    MSS[i] = SS[i]/DF[i]
    F[i] = ftest(mod.model,newmod).fstat[1]
  else
    SS[i] = sum(abs2.(eff[mask]))
    RSS[i] = 0
    DF[i] = sum(mask)
    MSS[i] = SS[i]/DF[i]
    F[i] = MSS[i]/MSSres
  end

  pval[i] = ccdf(FDist(DF[i], DFres), F[i])
  print(terms[i])
end


output2 = DataFrame( 
          Source = vcat(terms,"Residuals", "Total"),
          DF = vcat(DF,DFres, "-"),
          SS = vcat(SS, SSres, "-"),
          MSS = vcat(MSS, MSSres, "-"),
          F = vcat(F, "", ""),
          p = vcat(pval, "", "")
          )
end

export anova

end # module
