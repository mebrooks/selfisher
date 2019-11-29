# selfisher

`selfisher` is an R package for estimating the selectivity of fisheries gear, built using code from [glmmTMB](https://github.com/glmmTMB/glmmTMB), which is built on [Template Model Builder](https://github.com/kaskr/adcomp), which is in turn built on CppAD and Eigen. It has been compared to other software used for the same analyses. Fixed and random effects models can be specified for the selectivity and  relative fishing power components of the model. Possible links for the selectivity model include logit (i.e. logistic), probit (i.e. normal probability ogive), complimentary log-log (i.e. negative extreme value), log-log (i.e. extreme value/Gompertz), and Richards.

## Installation 

You can install `selfisher` via
```
remotes::install_github("mebrooks/selfisher/selfisher", build_vignette = TRUE)
```
If you haven't done so before, you may need to install Rtools on Windows from [here](https://cran.r-project.org/bin/windows/Rtools/). You'll need to have development tools (compilers etc.) installed: `devtools::dr_devtools()`. You may also need to install other packages that are used in the vignettes and examples
```
install.packages(c("knitr", "rmarkdown", "ggplot2", "bbmle", "plyr"))
install.packages('TMB', type = 'source')
```

If you get a warning saying `In checkMatrixPackageVersion() : Package version inconsistency detected.`, then you may need to first install TMB from source using the command `install.packages('TMB', type = 'source')`, then reinstall selfisher. Installing the source version will ensure that you get the very latest version of the package; since the package is in rapid development, that's a good idea. 
