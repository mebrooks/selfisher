pkgname <- "selfisher"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('selfisher')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("VarCorr.selfisher")
### * VarCorr.selfisher

flush(stderr()); flush(stdout())

### Name: VarCorr.selfisher
### Title: Extract variance and correlation components
### Aliases: VarCorr.selfisher VarCorr
### Keywords: internal

### ** Examples

## Comparing variance-covariance matrix with manual computation



cleanEx()
nameEx("fixef")
### * fixef

flush(stderr()); flush(stdout())

### Name: fixef
### Title: Extract fixed-effects estimates
### Aliases: fixef fixef.selfisher
### Keywords: models

### ** Examples

data(haddock)
dat=transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))
fixef(selfisher(prop~Lengths, p=~1, psplit=TRUE, total=tot, dat))



cleanEx()
nameEx("formFuns")
### * formFuns

flush(stderr()); flush(stdout())

### Name: inForm
### Title: test formula: does it contain a particular element?
### Aliases: inForm extractForm dropHead drop.special2
### Keywords: internal

### ** Examples

inForm(z~.,quote(.))
inForm(z~y,quote(.))
inForm(z~a+b+c,quote(c))
inForm(z~a+b+(d+e),quote(c))
f <- ~ a + offset(x)
f2 <- z ~ a
inForm(f,quote(offset))
inForm(f2,quote(offset))
extractForm(~a+offset(b),quote(offset))
extractForm(~c,quote(offset))
extractForm(~a+offset(b)+offset(c),quote(offset))
dropHead(~a+offset(b),quote(offset))
dropHead(~a+poly(x+z,3)+offset(b),quote(offset))



cleanEx()
nameEx("getGrpVar")
### * getGrpVar

flush(stderr()); flush(stdout())

### Name: getGrpVar
### Title: Get Grouping Variable
### Aliases: getGrpVar
### Keywords: internal

### ** Examples

data(cbpp,package="lme4")
cbpp$obs <- factor(seq(nrow(cbpp)))
rt <- lme4::glFormula(cbind(size,incidence-size)~(1|herd)+(1|obs),
  data=cbpp,family=binomial)$reTrms
getGrpVar(rt$flist)



cleanEx()
nameEx("getReStruc")
### * getReStruc

flush(stderr()); flush(stdout())

### Name: getReStruc
### Title: Calculate random effect structure Calculates number of random
###   effects, number of parameters, blocksize and number of blocks.
###   Mostly for internal use.
### Aliases: getReStruc

### ** Examples

data(sleepstudy, package="lme4")
rt <- lme4::lFormula(Reaction~Days+(1|Subject)+(0+Days|Subject),
                    sleepstudy)$reTrms
rt2 <- lme4::lFormula(Reaction~Days+(Days|Subject),
                    sleepstudy)$reTrms
getReStruc(rt)



cleanEx()
nameEx("numFactor")
### * numFactor

flush(stderr()); flush(stdout())

### Name: numFactor
### Title: Factor with numeric interpretable levels.
### Aliases: numFactor parseNumLevels

### ** Examples

## 1D example
numFactor(sample(1:5,20,TRUE))
## 2D example
coords <- cbind( sample(1:5,20,TRUE), sample(1:5,20,TRUE) )
(f <- numFactor(coords))
parseNumLevels(levels(f)) ## Sorted
## Used as part of a model.matrix
model.matrix( ~f )
## parseNumLevels( colnames(model.matrix( ~f )) )
## Error: 'Failed to parse numeric levels: (Intercept)'
parseNumLevels( colnames(model.matrix( ~ f-1 )) )



cleanEx()
nameEx("predict.selfisher")
### * predict.selfisher

flush(stderr()); flush(stdout())

### Name: predict.selfisher
### Title: prediction
### Aliases: predict.selfisher

### ** Examples

data(haddock)
dat <- transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))
m1 <- selfisher(prop~Lengths, p=~1, psplit=TRUE, total=tot, dat)
nd <- data.frame(Lengths=20:50, tot=100)
predict(m1, newdata=nd, se.fit=TRUE)



cleanEx()
nameEx("ranef.selfisher")
### * ranef.selfisher

flush(stderr()); flush(stdout())

### Name: ranef.selfisher
### Title: Extract Random Effects
### Aliases: ranef.selfisher ranef

### ** Examples

data(ccmhsdat)
ranef(selfisher(prop~length*type+(1|haul), total=total, ccmhsdat))
print(ranef(selfisher(prop~length*type+(1|haul), total=total, ccmhsdat)), simplify=FALSE)



cleanEx()
nameEx("selfisher")
### * selfisher

flush(stderr()); flush(stdout())

### Name: selfisher
### Title: Fit gear selectivity models with TMB
### Aliases: selfisher

### ** Examples

dat <- transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))
m0 <- selfisher(prop~Lengths, pformula=~0, psplit=TRUE, total=tot, dat)
m1 <- selfisher(prop~Lengths, pformula=~1, psplit=TRUE, total=tot, dat)



cleanEx()
nameEx("splitForm")
### * splitForm

flush(stderr()); flush(stdout())

### Name: addForm
### Title: Combine right-hand sides of an arbitrary number of formulas
### Aliases: addForm splitForm noSpecials
### Keywords: internal

### ** Examples

splitForm(~x+y)                     ## no specials or RE
splitForm(~x+y+(f|g))               ## no specials
splitForm(~x+y+diag(f|g))           ## one special
splitForm(~x+y+(diag(f|g)))         ## 'hidden' special
splitForm(~x+y+(f|g)+cs(1|g))       ## combination
splitForm(~x+y+(1|f/g))             ## 'slash'; term
splitForm(~x+y+(1|f/g/h))             ## 'slash'; term
splitForm(~x+y+(1|(f/g)/h))             ## 'slash'; term
splitForm(~x+y+(f|g)+cs(1|g)+cs(a|b,stuff))  ## complex special
splitForm(~(((x+y))))               ## lots of parentheses

noSpecials(y~1+us(1|f))
noSpecials(y~1+us(1|f),delete=FALSE)
noSpecials(y~us(1|f))
noSpecials(y~us+1)  ## should *not* delete unless head of a function



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
