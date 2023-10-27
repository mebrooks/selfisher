## ----setup, include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache=FALSE)


## ----pkgs----------------------------------------------------------------------------
library(boot) #for inv.logit
library(ggplot2); theme_set(theme_bw())
library(plyr)
library(selfisher)
set.seed(11)


## ----pars----------------------------------------------------------------------------
length=1:100 #measured length classes
pars1=c(-6, .1) #selection curve parameters for gear 1
pars2=c(-11, .3) #selection curve parameters for gear 2
nhaul=10 #number of hauls with each type of gear
q1=runif(nhaul, 0.6, 0.8) #sampling fractions for gear 1
q2=runif(nhaul, 0.2, 0.4) #sampling fractions for gear 2
psplit=0.5 # psplit will be obsorbed by the intercept
nfish=c(rep(0, 20), 
				rep(5, 20), 
				rep(10, 20), 
				rep(20, 20), 
				rep(10, 10),
				rep(1, 10)) # avg number of fish encountered in each length class in each haul


## ----dat-----------------------------------------------------------------------------
bdat=expand.grid(length = length, haul=1:nhaul)
bdat$nfish=nfish

bdat=transform(bdat, 
			sampling_test1 =q1[haul], 
			sampling_test2 =q2[haul], 
			r1=inv.logit(pars1[1]+ pars1[2]* length),
			r2=inv.logit(pars2[1]+ pars2[2]* length),
			simtotal=rpois(nrow(bdat), nfish)
			)


## ----plot1---------------------------------------------------------------------------
ggplot(bdat)+
	geom_line(aes(length, r1))+
	geom_line(aes(length, r2), lty=2)+
	ylab("true retention probabilities")


## ----sim_retention-------------------------------------------------------------------
dat=transform(bdat, length= length, haul=haul,
	test1=rbinom(nrow(bdat), size= simtotal, prob=r1*q1*psplit),
	test2=rbinom(nrow(bdat), size= simtotal, prob=r2*q2*(1-psplit)))


## ----rescale-------------------------------------------------------------------------
mu=mean(rep(dat$length, dat$test1+dat$test2))
var=var(rep(dat$length, dat$test1+dat$test2))
  
dat=transform(dat, 
      sl=(length-mu)/sqrt(var),
      q_ratio=sampling_test1/sampling_test2,
      prop=test1/(test1+test2),
      total=test1+test2,
      true_raised_ratio=r1/r2)

dat=subset(dat, !is.na(prop))


## ----sumdat--------------------------------------------------------------------------
sumdat=ddply(dat, ~length + sl + true_raised_ratio, summarize, 
      prop=sum(test1)/sum(total), 
      total=sum(total),
      raised_prop=sum(test1/sampling_test1)/
                        (sum(test1/sampling_test1)+sum(test2/sampling_test2)),
      raised_total=sum(test1/sampling_test1)+sum(test2/sampling_test2),
      raised_ratio = sum(test1/sampling_test1)/sum(test2/sampling_test2))


## ----plotmod-------------------------------------------------------------------------
p1=ggplot(sumdat, aes(x=length,))+
  geom_point(aes(y= raised_ratio, size=total), alpha=.5)+
  geom_hline(yintercept = 1, lty=2)+
  geom_line(data=sumdat, aes(length, true_raised_ratio), colour="red", alpha=.5, lwd=2)+
  ylab("ratio")+
  coord_cartesian(ylim=c(0, 3))
p1


## ----model---------------------------------------------------------------------------
m1=selfisher(prop~offset(log(q_ratio))+sl+I(sl^2)+I(sl^3)+I(sl^4), total=total, dat, haul=haul)


## ----preds---------------------------------------------------------------------------
newdata=data.frame(length=unique(dat$length))
newdata=transform(newdata,
    sl=(length-mu)/sqrt(var),
    total=1,
    haul=NA,
    q_ratio=1
)

newdata$est_raised_ratio=predict(m1, newdata=newdata, type="ratio")


## ----boot, message=FALSE, warning=FALSE----------------------------------------------
bs=bootSel(m1, nsim=1000, parallel = "multicore", ncpus = 4,
        FUN=function(mod){predict(mod, newdata=newdata, type="ratio")})



## ----bootwin, eval=FALSE-------------------------------------------------------------
## library(snow)
## ncpus = 4
## cl = makeCluster(rep("localhost", ncpus)) clusterExport(cl, "newdata")
## bs = bootSel(m1, nsim=1000, parallel = "snow", cl=cl,
##     FUN=function(mod){predict(mod, newdata=newdata, type="ratio")})
## stopCluster(cl)


## ----plotboot, message=FALSE, warning=FALSE------------------------------------------
#apply quantile function to bootstraps and match them with the newdata used for predictions 
quants=apply(bs$t, 2, quantile, c(0.025, 0.5, 0.975)) 

newdata[,c("lo", "mid", "hi")]=t(quants)

p1+
	geom_ribbon(data=newdata, aes(ymin=lo, ymax=hi), alpha=0.2)+
	geom_line(data=newdata, aes(y=est_raised_ratio))
	


## ----model_sp------------------------------------------------------------------------
library(splines)

m1_sp=selfisher(prop~offset(log(q_ratio))+bs(sl, df=4), total=total, dat, haul=haul)


## ----preds_sp------------------------------------------------------------------------
newdata$raised_ratio_sp=predict(m1_sp, newdata=newdata, type="ratio")


## ----boot_sp, message=FALSE, warning=FALSE-------------------------------------------

bs=bootSel(m1_sp, nsim=1000, parallel = "multicore", ncpus = 4,
        FUN=function(mod){predict(mod, newdata=newdata, type="ratio")})

quants=apply(bs$t, 2, quantile, c(0.025, 0.5, 0.975))

newdata[,c("lo", "mid", "hi")]=t(quants)

p1+
  geom_line(data=newdata, aes(length, raised_ratio_sp), lty=3)+
  geom_line(data=newdata, aes(y=est_raised_ratio))+
  geom_ribbon(data=newdata, aes(ymin=lo, ymax=hi), alpha=0.2)



## ----aic-----------------------------------------------------------------------------
library(bbmle)
AICtab(m1, m1_sp)

