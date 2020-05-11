## ----echo=FALSE---------------------------------------------------------------
library(selfisher)

## ----scale--------------------------------------------------------------------
m0=selfisher(nwide/(nfine+nwide)~scale(Lengths), total=nfine+nwide, haddock)

## ----NA function, eval=FALSE--------------------------------------------------
#  Warning in nlminb(start = par, objective = fn, gradient = gr) :
#    NA/NaN function evaluation

## ----Cholmod, eval=FALSE------------------------------------------------------
#  Cholmod warning 'matrix not positive definite'

## ----lgamma, eval=FALSE-------------------------------------------------------
#  Warning in f(par, order = order, ...) : value out of range in 'lgamma'

