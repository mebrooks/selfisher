## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache=TRUE)

## ----libs----------------------------------------------------------------
library(selfisher)
library(plyr)
library(ggplot2); theme_set(theme_bw())
library(bbmle)
library(knitr)
library(splines) 

## ----data, echo=TRUE-----------------------------------------------------
data(ccmhsampdat)
dat=ccmhsampdat
head(dat)

## ----scale---------------------------------------------------------------

mu=mean(rep(dat$length, dat$total))
var=var(rep(dat$length, dat$total))
  
dat=transform(dat, 
				q=sampling_test1/sampling_test2,
				sl=(length-mu)/sqrt(var))

## ----aggregate-----------------------------------------------------------
sumdat=ddply(dat, ~length+sl, summarize, 
				prop=sum(test1)/sum(total), 
				total=sum(total),
				raised_prop=sum(test1/sampling_test1)/
						(sum(test1/sampling_test1)+sum(test2/sampling_test2)),
				raised_total=sum(test1/sampling_test1)+sum(test2/sampling_test2),
				raised_ratio = sum(test1/sampling_test1)/sum(test2/sampling_test2))

## ----fits, warning=FALSE-------------------------------------------------
m1=selfisher(prop~offset(log(q))+sl+I(sl^2)+I(sl^3)+I(sl^4), total=total, dat, haul=haul)

## ----pred----------------------------------------------------------------
newdata=data.frame(length=unique(dat$length))
newdata=transform(newdata,
				sl=(length-mu)/sqrt(var),
				total=1,
				haul=1,
				q=1
)
newdata$raised_ratio=predict(m1, newdata=newdata, type="ratio")

## ----plot----------------------------------------------------------------
p1=ggplot(sumdat, aes(length, raised_ratio))+
	geom_point(aes(size=total))+
   geom_line(data=newdata)+
	geom_hline(yintercept = 1, lty=2)
p1

## ----ci------------------------------------------------------------------
bs=bootSel(m1, nsim=1000, parallel = "multicore", ncpus = 4, FUN=function(mod){predict(mod, newdata=newdata, type="ratio")})

quants=apply(bs$t, 2, quantile, c(0.025, 0.5, 0.975))
newdata[,c("lo", "mid", "hi")]=t(quants)

ggplot(sumdat, aes(length, raised_ratio))+
	geom_point(aes(size=total))+
   geom_line(data=newdata)+
   geom_ribbon(data=newdata, aes(ymin=lo, ymax=hi), alpha=0.2)+
	geom_hline(yintercept = 1, lty=2)

## ----spline_fits---------------------------------------------------------
s3=selfisher(prop~offset(log(q))+ns(sl, 3), psplit=FALSE, total=total, dat, haul=haul)
s4=selfisher(prop~offset(log(q))+ns(sl, 4), psplit=FALSE, total=total, dat, haul=haul)
s5=selfisher(prop~offset(log(q))+ns(sl, 5), psplit=FALSE, total=total, dat, haul=haul)
AICtab(s3, s4, s5, m1)

## ----spline_pred---------------------------------------------------------
newdata=data.frame(length=unique(dat$length))
newdata=transform(newdata,
          sl=(length-mean(dat$length))/sqrt(var(dat$length)),
          total=1,
          haul=1,
					q=1
)
newdata$raised_ratio=predict(s3, newdata=newdata, type="ratio")

## ----spline_ci-----------------------------------------------------------
bs=bootSel(s3, nsim=100, parallel = "multicore", ncpus = 4, FUN=function(mod){predict(mod, newdata=newdata, type="ratio")})

quants=apply(bs$t, 2, quantile, c(0.025, 0.5, 0.975))
newdata[,c("lo", "mid", "hi")]=t(quants)

## ----spline_plot---------------------------------------------------------
p2=ggplot(sumdat, aes(length, raised_ratio))+
	geom_point(aes(size=total))+
   geom_line(data=newdata)+
   geom_ribbon(data=newdata, aes(ymin=lo, ymax=hi), alpha=0.2)+
	geom_hline(yintercept = 1, lty=2)
p2

