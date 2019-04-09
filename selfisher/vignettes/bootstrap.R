## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----preliminaries-------------------------------------------------------
library(selfisher)
library(plyr)
library(ggplot2); theme_set(theme_bw())
data("ccmhsdat")
head(ccmhsdat)

## ----aggregate-----------------------------------------------------------
sumdat=ddply(ccmhsdat, ~length+type, summarize, prop=sum(test)/sum(total), total=sum(total))

## ----fits----------------------------------------------------------------
mod_both=selfisher(prop~length*type, total=total, ccmhsdat, haul=haul)

## ----pred----------------------------------------------------------------
newdata=expand.grid(length=unique(ccmhsdat$length),
                total=1,
                haul=1,
                type=c("baseline", "stimulation"))

newdata$prop=predict(mod_both, newdata=newdata, type="response")

## ----ci------------------------------------------------------------------
bs=bootSel(mod_both, nsim=100, parallel = "multicore", ncpus = 4, FUN=function(mod){predict(mod, newdata=newdata, type="response")})

quants=apply(bs$t, 2, quantile, c(0.025, 0.5, 0.975))
newdata[,c("lo", "mid", "hi")]=t(quants)

## ----plot----------------------------------------------------------------
ggplot(sumdat, aes(length, prop, colour=type))+geom_point()+
   geom_line(data=newdata)+
   geom_ribbon(data=newdata, aes(ymin=lo, ymax=hi, fill=type), alpha=0.2)

