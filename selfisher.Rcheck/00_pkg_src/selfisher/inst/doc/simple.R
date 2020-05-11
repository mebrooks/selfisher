## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----libs---------------------------------------------------------------------
library(selfisher)
library(bbmle) #for AICtab

## ----dat----------------------------------------------------------------------
data(haddock)
head(haddock)
dat=transform(haddock, 
              tot=nfine+nwide, 
              prop=nwide/(nfine+nwide))

## ----mod0---------------------------------------------------------------------
m0=selfisher(prop~Lengths, pformula=~0, psplit=TRUE, total=tot, dat)

## ----mod1---------------------------------------------------------------------
m1=selfisher(prop~Lengths, pformula=~1, psplit=TRUE, total=tot, dat)

## ----aic----------------------------------------------------------------------
AICtab(m0, m1)

## ----sum----------------------------------------------------------------------
summary(m1)

## ----pred---------------------------------------------------------------------
plot(dat$Length, residuals(m1, type="deviance"), type='h', ylim=c(-1,2))
abline(0,0)

