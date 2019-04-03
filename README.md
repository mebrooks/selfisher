# selfisher

`selfisher` is an R package for estimating the selectivity of fisheries gear, built using code from [glmmTMB](https://github.com/glmmTMB/glmmTMB), which is built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on CppAD and Eigen. It is currently pre-alpha or alpha software. Fixed and random effects models can be specified for the selectivity and  relative fishing power components of the model. Possible links for the selectivity model include logit (i.e. logistic), probit (i.e. normal probability ogive), complimentary log-log (i.e. negative extreme value), log-log (i.e. extreme value/Gompertz), and Richards.

## Installation 

You can install `selfisher` via
```
devtools::install_github("mebrooks/selfisher/selfisher")
```
(this string denotes "Github user `mebrooks`, repository `selfisher`, subdirectory `selfisher`"). If the install fails at the vignette-building step, try specifying `build_vignettes=FALSE` within the `install_github` call. Alternatively you can use `install_github()` from the `ghit` package, which has fewer dependencies. You'll need to have development tools (compilers etc.) installed: `devtools::dr_devtools()` and the [RStudio devtools docs](https://www.rstudio.com/products/rpackages/devtools/) should help. Installing the source version will ensure that you get the very latest version of the package: since the package is in rapid development, that's a good idea. 
