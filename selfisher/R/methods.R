##' Extract the fixed-effects estimates
##'
##' Extract the estimates of the fixed-effects parameters from a fitted model.
##' @name fixef
##' @title Extract fixed-effects estimates
##' @aliases fixef fixef.selfisher
##' @docType methods
##' @param object any fitted model object from which fixed effects estimates can
##' be extracted.
##' @param \dots optional additional arguments. Currently none are used in any
##' methods.
##' @return a named, numeric vector of fixed-effects estimates.
##' @keywords models
##' @examples
##' data(haddock)
##' dat=transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))
##' fixef(selfisher(prop~Lengths, p=~1, psplit=TRUE, total=tot, dat))
##' @importFrom nlme fixef
##' @export fixef
##' @export
fixef.selfisher <- function(object, ...) {
  pl <- object$obj$env$parList(object$fit$par, object$fit$parfull)
  structure(list(r = setNames(pl$betar,   colnames(getME(object, "Xr"))),
                 p    = setNames(pl$betap, colnames(getME(object, "Xp"))),
                 d = setNames(pl$betad, colnames(getME(object, "Xd")))),
            class =  "fixef.selfisher")
}

## general purpose matching between component names and printable names
cNames <- list(r = "Selectivity model",
               p = "Relative fishing power model",
               d = "Richards exponent model")

trivialDisp <- function(object) {
    ## This version works on summary object or fitted model object
    object$modelInfo$family$link!="Richards" ||(
      identical(deparse(object$call$dformula),"~1") &
      object$modelInfo$family$link=="Richards")
}
trivialFixef <- function(xnm,nm) {
    length(xnm)==0 ||
        (nm %in% c('d') && identical(xnm,'(Intercept)'))
}


##' @method print fixef.selfisher
##' @export
print.fixef.selfisher <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  for(nm in names(x)) {
      if (!trivialFixef(names(x[[nm]]),nm)) {
          cat(sprintf("\n%s:\n", cNames[[nm]]))
          print.default(format(x[[nm]], digits=digits), print.gap = 2L, quote = FALSE)
      }
  }
  invisible(x)
}

##' Extract Random Effects
##'
##' Generic function to extract random effects from \code{selfisher} models, both
##' for the selectivity (i.e. retention) model and relative fising power model.
##'
##' @param object a \code{selfisher} model.
##' @param ... some methods for this generic function require additional
##'   arguments.
##'
##' @return Object of class \code{ranef.selfisher} with two components:
##'   \item{r}{a list of data frames, containing random effects
##'     for the selectivity (i.e. retention) model.}
##'   \item{p}{a list of data frames, containing random effects for
##'     the relative fising power model.}
##'
##' @note When a model has no model of relative fishing power, the default behavior of
##'   \code{ranef} is to simplify the printed format of the random effects. To
##'   show the full list structure, run \code{print(ranef(model),
##'   simplify=FALSE)}. In all cases, the full list structure is used to access
##'   the data frames (see example).
##'
##' @seealso \code{\link{fixef.selfisher}}.
##'##'
##' @examples
##' data(ccmhsdat)
##' ranef(selfisher(prop~length*type+(1|haul), total=total, ccmhsdat))
##' print(ranef(selfisher(prop~length*type+(1|haul), total=total, ccmhsdat)), simplify=FALSE)
##' @aliases ranef ranef.selfisher
##' @importFrom nlme ranef
##' @export ranef
##' @export
ranef.selfisher <- function(object, ...) {
  ## The arrange() function converts a vector of random effects to a list of
  ## data frames, in the same way as lme4 does.
  arrange <- function(x, listname)
  {
    cnms <- object$modelInfo$reTrms[[listname]]$cnms
    flist <- object$modelInfo$reTrms[[listname]]$flist
    if (!is.null(cnms)) {
      levs <- lapply(fl <- flist, levels)
      asgn <- attr(fl, "assign")
      nc <- vapply(cnms, length, 1L)
      nb <- nc * vapply(levs, length, 1L)[asgn]
      nbseq <- rep.int(seq_along(nb), nb)
      ml <- split(x, nbseq)
      for (i in seq_along(ml))
        ml[[i]] <- matrix(ml[[i]], ncol=nc[i], byrow=TRUE,
                          dimnames=list(NULL, cnms[[i]]))
      x <- lapply(seq_along(fl), function(i)
        data.frame(do.call(cbind, ml[asgn==i]), row.names=levs[[i]],
                   check.names=FALSE))
      names(x) <- names(fl)
      x
    }
    else {
      list()
    }
  }

  pl <- getParList(object)
  structure(list(r = arrange(pl$br, "r"),
                 p    = arrange(pl$bp, "p")),
            class = "ranef.selfisher")
}

##' @method print ranef.selfisher
##' @export
print.ranef.selfisher <- function(x, simplify=TRUE, ...) {
    print(if (simplify && length(x$p) == 0L)
              unclass(x$r) else unclass(x),
          ...)
    invisible(x)
}


##' Extract or Get Generalize Components from a Fitted Mixed Effects Model
##'
##' @aliases getME
##' @param object a fitted \code{selfisher} object
##' @param name of the component to be retrieved
##' @param \dots ignored, for method compatibility
##'
##' @seealso \code{\link[lme4]{getME}}
##' Get generic and re-export:
##' @importFrom lme4 getME
##' @export getME
##'
##' @method getME selfisher
##' @export
getME.selfisher <- function(object,
                          name = c("Xr", "Xp","Zr", "Zp", "Xd",
                                   "betar", "betap", "betad",
                                   "br", "bp", "thetar", "thetap"),
                          ...)
{
  if(missing(name)) stop("'name' must not be missing")
  ## Deal with multiple names -- "FIXME" is inefficiently redoing things
  if (length(name <- as.character(name)) > 1) {
    names(name) <- name
    return(lapply(name, getME, object = object))
  }
  if(name == "ALL") ## recursively get all provided components
      return(sapply(eval(formals()$name),
                    getME.selfisher, object=object, simplify=FALSE))

  stopifnot(inherits(object, "selfisher"))
  name <- match.arg(name)

  oo.env <- object$obj$env
  ### Start of the switch
  switch(name,
         "Xr"    = oo.env$data$Xr,
         "Xp"    = oo.env$data$Xp,
         "Zr"    = oo.env$data$Zr,
         "Zp"    = oo.env$data$Zp,
         "Xd"    = oo.env$data$Xd,
         "betar" = oo.env$parList(object$fit$par, object$fit$parfull)$betar,
         "betap" = oo.env$parList(object$fit$par, object$fit$parfull)$betap,
         "betad" = oo.env$parList(object$fit$par, object$fit$parfull)$betad,
         "br" = oo.env$parList(object$fit$par, object$fit$parfull)$br,
         "bp" = oo.env$parList(object$fit$par, object$fit$parfull)$bp,
         "thetar" = oo.env$parList(object$fit$par, object$fit$parfull)$thetar,
         "thetap" = oo.env$parList(object$fit$par, object$fit$parfull)$thetap,

         "..foo.." = # placeholder!
           stop(gettextf("'%s' is not implemented yet",
                         sprintf("getME(*, \"%s\")", name))),
         ## otherwise
         stop(sprintf("Mixed-Effects extraction of '%s' is not available for class \"%s\"",
                      name, class(object))))
}## {getME}

## FIXME: (1) why is this non-standard (containing nobs, nall?)
##        (2) do we really need to document it??
## Extract the log likelihood of a selfisher model
##
## @return object of class \code{logLik} with attributes
## \item{val}{log likelihood}
## \item{nobs,nall}{number of non NA observations initially supplied to TMB}
## \item{df}{number of parameters}
##' @importFrom stats logLik
##' @export
logLik.selfisher <- function(object, ...) {
  if(!is.null(object$sdr)){
    val <- if(object$sdr$pdHess){-object$fit$objective}else{NA}
  }else val <- -object$fit$objective

  nobs <- nobs.selfisher(object)
  structure(val, nobs = nobs, nall = nobs, df = length(object$fit$par),
            class = "logLik")
}

##' @importFrom stats nobs
##' @export
nobs.selfisher <- function(object, ...) sum(!is.na(object$obj$env$data$yobs))

##' @importFrom stats df.residual
##' @method df.residual selfisher
##' @export
##  TODO: not clear whether the residual df should be based
##  on p=length(beta) or p=length(c(theta,beta)) ... but
##  this is just to allow things like aods3::gof to work ...
##  Taken from LME4, including the todo
##
df.residual.selfisher <- function(object, ...) {
  nobs(object)-length(object$fit$par)
}


##' Calculate Variance-Covariance Matrix for a Fitted selfisher model
##'
##' @param object a \dQuote{selfisher} fit
##' @param full return a full variance-covariance matrix?
##' @param \dots ignored, for method compatibility
##' @return By default (\code{full==FALSE}), a list of separate variance-covariance matrices for each model component (conditional, zero-inflation, dispersion).  If \code{full==TRUE}, a single square variance-covariance matrix for \emph{all} top-level model parameters (conditional, dispersion, and variance-covariance parameters)
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats vcov
##' @export
vcov.selfisher <- function(object, full=FALSE, ...) {
  if(is.null(sdr <- object$sdr)) {
    warning("Calculating sdreport. Use se=TRUE in selfisher to avoid repetitive calculation of sdreport")
    sdr <- sdreport(object$obj)
  }
  keepTag <- if (full) { "."
             } else if (!trivialDisp(object)) { "beta*"
             } else "beta($|[^d])"
  to_keep <- grep(keepTag,colnames(sdr$cov.fixed)) # only keep betas
  covF <- sdr$cov.fixed[to_keep,to_keep,drop=FALSE]

  mkNames <- function(tag) {
      X <- getME(object,paste0("X",tag))
      if (trivialFixef(nn <- colnames(X),tag) &&
          ## if 'full', keep d even if trivial
          !(full && tag =="d")) character(0)
      else paste(tag,nn,sep="~")
  }

  nameList <- setNames(list(colnames(getME(object,"Xr")),
                       mkNames("p"),
                       mkNames("d")),
                names(cNames))

  if(full) {
      ## FIXME: haven't really decided if we should drop the
      ##   trivial variance-covariance dispersion parameter ??
      ## if (trivialDisp(object))
      ##    res <- covF[-nrow(covF),-nrow(covF)]

      reNames <- function(tag) {
          re <- object$modelInfo$reStruc[[paste0(tag,"ReStruc")]]
          nn <- mapply(function(n,L) paste(n,seq(L),sep="."),
                 names(re),
                 sapply(re,"[[","blockNumTheta"))
          if (length(nn)==0) return(nn)
          return(paste("theta",gsub(" ","",nn),sep="_"))
      }
      nameList <- c(nameList,list(reNames("r"),reNames("p")))

      colnames(covF) <- rownames(covF) <- unlist(nameList)
      res <- covF        ## return just a matrix in this case
  } else {
      splitMat <- function(x) {
          ss <- split(seq_along(colnames(x)),
                      colnames(x))
          lapply(ss,function(z) x[z,z,drop=FALSE])
      }
      covList <- splitMat(covF)
      names(covList) <-
          names(cNames)[match(names(covList),c("betar","betap","betad"))]
      for (nm in names(covList)) {
          if (length(xnms <- nameList[[nm]])==0) {
              covList[[nm]] <- NULL
          }
          else dimnames(covList[[nm]]) <- list(xnms,xnms)
      }
      res <- covList
      ##  FIXME: should vcov always return a three-element list
      ## (with NULL values for trivial models)?
      class(res) <- c("vcov.selfisher","matrix")
  }
  return(res)
}

##' @method print vcov.selfisher
##' @export
print.vcov.selfisher <- function(x,...) {
    for (nm in names(x)) {
        cat(cNames[[nm]],":\n",sep="")
        print(x[[nm]])
        cat("\n")
    }
    invisible(x)
}

cat.f <- function(...) cat(..., fill = TRUE)

.prt.call.selfisher <- function(call, long = TRUE) {
  pass <- 0
  if (!is.null(cc <- call$rformula)){
    cat.f("Selectivity formula:         ", deparse(cc))
    rhs <- cc[[2]]
    if (!is.null(rhs)) {
        pass<-nchar(deparse(rhs))
    }
  }
  if(!identical(cc <- deparse(call$pformula),"~0") & call$psplit)
    cat.f("Relative fishing power formula:  ",rep(' ',pass+2), cc, sep='')
  if(!identical(cc <- deparse(call$dformula),"~1"))
    cat.f("Richards exponent:      ",rep(' ',pass+2), cc, sep='')
  if (!is.null(cc <- call$data))
    cat.f("Data:", deparse(cc))
  if (!is.null(cc <- call$total))
    cat.f("Total:", deparse(cc))
  if (!is.null(cc <- call$offset))
    cat.f(" Offset:", deparse(cc))
#  if (long && length(cc <- call$control) &&
#      !identical((dc <- deparse(cc)), "lmerControl()"))
    ## && !identical(eval(cc), lmerControl()))
#    cat.f("Control:", dc)
#  if (!is.null(cc <- call$subset))
#    cat.f(" Subset:", deparse(cc))
}

.prt.retention <- function(ret, SR) {
	if (!is.null(ret)){
		cat("\nSize at retention probability:\n")
		print(ret, rownames=FALSE)
		cat("\nSelectivity range (SR):\n")
		print(SR)
	}
}

### FIXME: attempted refactoring ...
cat.f2 <- function(call,component,label,lwid,fwid=NULL,cind=NULL) {
    if (!is.null(cc <- call[[component]])) {
        if (!is.null(cind)) {
            ## try to extract component (of formula)
            if (!is.null(ccc <- cc[[cind]]))
                cc <- ccc
        }
        f1 <- format(paste0(label,":"),width=lwid,justify="right")
        f2 <- deparse(cc)
        if (!is.null(fwid)) {
            f2 <- format(f2,width=fwid,justify="right")
        }
        cat(f1,f2,fill=TRUE)
    }
}

## don't use ##' until we're ready to generate a man page
## @param s delta (results of delta(x) for original object
printDispersion <- function(s) {

        dname <- "Richards exponent parameter"
        sname <- "delta"
        sval <- s
        cat(sprintf("\n%s (%s): %s",
                    dname,sname,
                    formatC(sval,digits=3)),"\n")
}

##' @importFrom lme4 .prt.aictab
##' @method print selfisher
##' @export
print.selfisher <-
    function(x, digits = max(3, getOption("digits") - 3),
             correlation = NULL, symbolic.cor = FALSE,
             signif.stars = getOption("show.signif.stars"),
             longCall = TRUE, ranef.comp = "Std.Dev.", ...)
{
  ## Type Of Model fit --- REML? ---['class']  & Family & Call
  .prt.call.selfisher(x$call, long=longCall)
  ## the 'digits' argument should have an action here
  aictab <- c(AIC = AIC(x), BIC = BIC(x), logLik = logLik(x),
              df.resid = df.residual(x))
  .prt.aictab(aictab, digits=digits+1)
  ## varcorr
  if (!all(sapply(vc <- VarCorr(x),is.null))) {
      cat("Random-effects (co)variances:\n")
      print(VarCorr(x), digits=digits, comp = ranef.comp)
  }
  ## ngroups
  gvec <- list(obs=sprintf("\nNumber of obs: %d",nobs(x)))
  ng <- ngrps.selfisher(x)
  for (i in seq_along(ng)) {
      if (length(ng[[i]])>0) {
          nm <- names(ng)[i]
          gvec[[nm]] <- paste0(cNames[nm],": ",
                      paste(paste(names(ng[[i]]), ng[[i]], sep=", "), collapse="; "))
      }
  }
  cat(do.call(paste,c(gvec,list(sep=" / "))),fill=TRUE)

  if(trivialDisp(x) & x$modelInfo$link=="Richards") {# if trivial print here, else below(~x) or none(~0)
    printDispersion(Richardsdelta(x))
  }
  ## Fixed effects:
  if(length(cf <- fixef(x)) > 0) {
    cat("\nFixed Effects:\n")
    print(cf, ...)
  } else
    cat("No fixed effect coefficients\n")
  invisible(x)
}

##' @export
model.frame.selfisher <- function(formula, ...) {
    formula$frame
}


##' Compute residuals for a selfisher object
##'
##' @param object a \dQuote{selfisher} object
##' @param type (character) residual type
##' @param \dots ignored, for method compatibility
##' @importFrom stats fitted model.response residuals
##' @export
residuals.selfisher <- function(object, type=c("response", "pearson", "deviance"), ...) {
#     stop("residuals are not implemented yet")
    type <- match.arg(type)
    r <- model.response(object$frame)-fitted(object)
    switch(type,
           response=r,
           pearson={
#               if (is.null(v <- family(object)$variance))
#                   stop("variance function undefined for family ",
#                        sQuote(family(object)$family),"; cannot compute",
#                        " Pearson residuals")
#               vv <- switch(length(formals(v)),
#                            v(fitted(object)),
#                            v(fitted(object),delta(object)),
#                            stop("variance function should take 1 or 2 arguments"))
#               r/sqrt(vv)
               y <- model.response(object$frame)
               yhat <- fitted(object)
               n <- model.total(object$frame)
               sqrt(n)*(y-yhat)/sqrt(yhat*(1-yhat))
           },
           deviance={
               y <- model.response(object$frame)
               yhat <- fitted(object)
               n <- model.total(object$frame)
               sign(y-yhat)*(binomial()$dev.resids(y, yhat, n))^.5
           }
       )
}

## Helper to get CI of simple *univariate monotone* parameter
## function, i.e. a function of 'fit$par' and/or 'fit$parfull'.
## Examples: 'sigma.glmmTMB' and some parts of 'VarCorr.glmmTMB'.

##' @importFrom stats qchisq
.CI_univariate_monotone <- function(object, f, reduce=NULL,
                                    level=0.95,
                                    name.prepend=NULL,
                                    estimate = TRUE) {
  x <- object
  par <- x$fit$par
  i <- seq_along(x$fit$parfull) ## Pointers into long par vector
  r <- x$obj$env$random
  if(!is.null(r)) i <- i[-r]    ## Pointers into short par subset
  sdr <- x$sdr
  sdpar <- summary(sdr, "fixed")[,2]
  q <- sqrt(qchisq(level, df=1))
  ans <- list()
  x$fit$parfull[i] <- x$fit$par <- par - q * sdpar
  ans$lower <- f(x)
  x$fit$parfull[i] <- x$fit$par <- par + q * sdpar
  ans$upper <- f(x)
  if (estimate) {
    ans$Estimate <- f(object)
  }
  if(is.null(reduce)) reduce <- function(x) x
  ans <- lapply(ans, reduce)
  nm <- names(ans)
  tmp <- cbind(ans$lower, ans$upper)
  if (is.null(tmp) || nrow(tmp) == 0L) return (NULL)
  sort2 <- function(x) if(any(is.nan(x))) x * NaN else sort(x)
  ans <- cbind( t( apply(tmp, 1, sort2) ) , ans$Estimate )
  colnames(ans) <- nm
  if (!is.null(name.prepend))
    name.prepend <- rep(name.prepend, length.out = nrow(ans))
  rownames(ans) <- paste(name.prepend,
                         rownames(ans), sep="")
  ans
}

## copied from 'stats'

format.perc <- function (probs, digits) {
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
    "%")
}

##' Calculate confidence intervals
##'
##' @details
##' Available methods are
##' \describe{
##' \item{wald}{These intervals are based on the standard errors
##' calculated for parameters on the scale
##' of their internal parameterization depending on the family. Derived
##' quantities such as standard deviation parameters and dispersion
##' parameters are backtransformed. It follows that confidence
##' intervals for these derived quantities are asymmetric.}
##' \item{profile}{This method computes a likelihood profile
##' for the specified parameter(s) using \code{profile.glmmTMB};
##' fits a spline function to each half of the profile; and
##' inverts the function to find the specified confidence interval.}
##' \item{uniroot}{This method uses the \code{\link{uniroot}}
##' function to find critical values of one-dimensional profile
##' functions for each specified parameter.}
##' }
##' @param object \code{selfisher} fitted object.
##' @param parm Specification of a parameter subset \emph{after}
##'     \code{component} subset has been applied.
##' @param level Confidence level.
##' @param method 'wald', 'profile', or 'uniroot': see Details
##' function)
##' @param component Which of the three components 'r', 'p' or
##'     'other' to select. Default is to select 'all'.
##' @param estimate (logical) add a third column with estimate ?
##' @param parallel method (if any) for parallel computation
##' @param ncpus number of CPUs/cores to use for parallel computation
##' @param cl cluster to use for parallel computation
##' @param ... arguments may be passed to \code{\link{profile.selMod}} or
##' \code{\link{tmbroot}}
##' @importFrom stats qnorm
##' @importFrom stats confint
##' @export
confint.selfisher <- function (object, parm, level = 0.95,
                             method=c("wald",
                                      "Wald",
                                      "profile",
                                      "uniroot"),
                             component = c("all", "r", "p", "other"),
                             estimate = TRUE,
                             parallel = c("no", "multicore", "snow"),
                             ncpus = getOption("profile.ncpus", 1L),
                             cl = NULL,
                             ...)
{
    method <- tolower(match.arg(method))
    if (method=="wald") {
        dots <- list(...)
        if (length(dots)>0) {
            if (is.null(names(dots))) {
                warning("extra (unnamed) arguments ignored")
            } else {
                warning(paste("extra arguments ignored: ",
                              paste(names(dots),collapse=", ")))
            }
        }
    }
    components <- match.arg(component, several.ok = TRUE)
    components.has <- function(x)
        any(match(c(x, "all"), components, nomatch=0L)) > 0L
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- format.perc(a, 3)
    fac <- qnorm(a)
    estimate <- as.logical(estimate)
    ci <- matrix(NA, 0, 2 + estimate,
                 dimnames=list(NULL,
                               c(pct, "Estimate")
                               [c(TRUE, TRUE, estimate)] ))
    if (tolower(method)=="wald") {
        for (component in c("r", "p") ) {
            if (components.has(component)) {
                cf <- unlist(fixef(object)[component])
                vv <- vcov(object)[component]
                ss <- unlist(lapply(vv,diag))
                ses <- sqrt(ss)
                ci.tmp <- cf + ses %o% fac
                if (estimate) ci.tmp <- cbind(ci.tmp, cf)
                ci <- rbind(ci, ci.tmp)
                ## VarCorr -> stddev
                reduce <- function(VC) sapply(VC[[component]],
                                              function(x)attr(x, "stddev"))
                ci.sd <- .CI_univariate_monotone(object,
                                                 VarCorr,
                                                 reduce = reduce,
                                                 level = level,
                                                 name.prepend=paste(component,
                                                                    "Std.Dev.",
                                                                    sep="."),
                                                 estimate = estimate)
                ci <- rbind(ci, ci.sd)
            }
        }
        if (components.has("other")) {
            ## sigma
            ff <- object$modelInfo$family$family
            if (object$modelInfo$link=="Richards") {
                ci.sigma <- .CI_univariate_monotone(object,
                                                    sigma,
                                                    reduce = NULL,
                                                    level=level,
                                                    name.prepend="sigma",
                                                    estimate = estimate)
                ci <- rbind(ci, ci.sigma)
            }
        }
        ## Take subset
        if (!missing(parm)) {
            ## FIXME: DRY/refactor with confint.profile
            ## FIXME: beta_ not well defined; sigma parameters not
            ## distinguishable (all called "sigma")
            ## for non-trivial dispersion model
            theta_parms <- grep("\\.(Std\\.Dev|Cor)\\.",rownames(ci))
            ## if non-trivial disp, keep disp parms for "beta_"
            disp_parms <- if (!trivialDisp(object)) numeric(0) else grep("^sigma",rownames(ci))
            if (identical(parm,"theta_")) {
                parm <- theta_parms
            } else if (identical(parm,"beta_")) {
                parm <- seq(nrow(ci))[-c(theta_parms,disp_parms)]
            }
            ci <- ci[parm, , drop=FALSE]
        }
        ## end Wald method
    } else if (tolower(method=="uniroot")) {
        ## FIXME: allow greater flexibility in specifying different
        ##  ranges, etc. for different parameters
        if (missing(parm)) {
            parm <- seq_along(names(object$obj$par))
        }
        plist <- parallel_default(parallel,ncpus)
        parallel <- plist$parallel
        do_parallel <- plist$do_parallel
        FUN <- function(n) {
            tmbroot(obj=object$obj, name=n, target=0.5*qchisq(level,df=1),
                    ...)
        }
        if (do_parallel) {
            if (parallel == "multicore") {
                L <- parallel::mclapply(parm, FUN, mc.cores = ncpus)
            } else if (parallel=="snow") {
                if (is.null(cl)) {
                    ## start cluster
                    new_cl <- TRUE
                    cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                }
                ## run
                L <- parallel::clusterApply(cl, parm, FUN)
                if (new_cl) {
                    ## stop cluster
                    parallel::stopCluster(cl)
                }
            }
        } else { ## non-parallel
            L <- lapply(as.list(parm), FUN)
        }
        L <- do.call(rbind,L)
        rownames(L) <- rownames(vcov(object,full=TRUE))[parm]
        if (estimate) {
            ee <- object$obj$env
            par <- ee$last.par.best
            if (!is.null(ee$random))
                par <- par[-ee$random]
            par <- par[parm]
            L <- cbind(L,par)
        }
        ci <- rbind(ci,L) ## really just adding column names!
    }
    else {  ## profile CIs
        pp <- profile(object, parm=parm, level_max=level,
                      parallel=parallel,ncpus=ncpus,
                      ...)
        ci <- confint(pp)
    }
    return(ci)
}


##' @export
## FIXME: establish separate 'terms' components for
##   each model component (selectivity, random, fishing power, dispersion ...)
terms.selfisher <- function(x, component="r", part="fixed", ...) {
    if (part != "fixed") stop("only fixed terms currently available")
    return(x$modelInfo$reTrms[[component]]$terms[[part]])
    ## terms(x$frame)
}

##' @export
extractAIC.selfisher <- function(fit, scale, k = 2, ...) {
    L <- logLik(fit)
    edf <- attr(L,"df")
    return(c(edf,c(-2*L + k*edf)))
}

## deparse(.) returning \bold{one} string
## copied from lme4/R/utilities.R
## Protects against the possibility that results from deparse() will be
##       split after 'width.cutoff' (by default 60, maximally 500)
safeDeparse <- function(x, collapse=" ") paste(deparse(x, 500L), collapse=collapse)

abbrDeparse <- function(x, width=60) {
    r <- deparse(x, width)
    if(length(r) > 1) paste(r[1], "...") else r
}

##' @importFrom methods is
##' @importFrom stats var getCall pchisq anova
##' @export
anova.selfisher <- function (object, ..., model.names = NULL)
{
    mCall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
    ## detect multiple models, i.e. models in ...
    modp <- as.logical(vapply(dots, is, NA, "selfisher"))
    if (any(modp)) {
        mods <- c(list(object), dots[modp])
        nobs.vec <- vapply(mods, nobs, 1L)
        if (var(nobs.vec) > 0)
            stop("models were not all fitted to the same size of dataset")
        if (is.null(mNms <- model.names))
            mNms <- vapply(as.list(mCall)[c(FALSE, TRUE, modp)],
                           safeDeparse, "")
        if (any(duplicated(mNms))) {
            warning("failed to find unique model names, assigning generic names")
            mNms <- paste0("MODEL", seq_along(mNms))
        }
        if (length(mNms) != length(mods))
            stop("model names vector and model list have different lengths")
        names(mods) <- sub("@env$", "", mNms)
        llks <- lapply(mods, logLik)
        ii <- order(Df <- vapply(llks, attr, FUN.VALUE = numeric(1),
            "df"))
        mods <- mods[ii]
        llks <- llks[ii]
        Df <- Df[ii]
        calls <- lapply(mods, getCall)
        data <- lapply(calls, `[[`, "data")
        if (!all(vapply(data, identical, NA, data[[1]])))
            stop("all models must be fit to the same data object")
        header <- paste("Data:", abbrDeparse(data[[1]]))
        subset <- lapply(calls, `[[`, "subset")
        if (!all(vapply(subset, identical, NA, subset[[1]])))
            stop("all models must use the same subset")
        if (!is.null(subset[[1]]))
            header <- c(header, paste("Subset:", abbrDeparse(subset[[1]])))
        llk <- unlist(llks)
        chisq <- 2 * pmax(0, c(NA, diff(llk)))
        dfChisq <- c(NA, diff(Df))
        val <- data.frame(Df = Df, AIC = .sapply(llks, AIC),
            BIC = .sapply(llks, BIC), logLik = llk, deviance = -2 *
                llk, Chisq = chisq, `Chi Df` = dfChisq, `Pr(>Chisq)` = pchisq(chisq,
                dfChisq, lower.tail = FALSE), row.names = names(mods),
            check.names = FALSE)
        class(val) <- c("anova", class(val))
        forms <- lapply(lapply(calls, `[[`, "formula"), deparse)
        structure(val, heading = c(header, "Models:", paste(rep(names(mods),
            times = lengths(forms)), unlist(forms), sep = ": ")))
    } else stop("no single-model anova() method for selfisher")
}

#' @importFrom stats predict
#' @export
fitted.selfisher <- function(object, ...) {
    predict(object)
}


##' Simulate from a selfisher fitted model
##' @method simulate selfisher
##' @param object selfisher fitted model
##' @param nsim number of response lists to simulate. Defaults to 1.
##' @param seed random number seed
##' @param ... extra arguments
##' @details Random effects are also simulated from their estimated distribution.
##' Currently, it is not possible to condition on estimated random effects.
##' @return returns a list of vectors. The list has length \code{nsim}.
##' Each simulated vector of observations is the same size as the vector of response variables in the original data set.
##' @importFrom stats simulate
##' @export
simulate.selfisher<-function(object, nsim=1, seed=NULL, ...){
    if(!is.null(seed)) set.seed(seed)
    ret <- replicate(nsim, object$obj$simulate()$yobs, simplify=FALSE)
    ret
}

##' refit the same model to a new response
##' @method refit selfisher
##' @param object a fitted \code{selfisher} object
##' @param newdata a data set with the same predictors used in the model
##' @importFrom lme4 refit
##' @export
refit.selfisher <- function(object, newdata, ...) {
  cc <- getCall(object)
  cc$data <- newdata
  return(eval(cc))
}

##' read in data from a single haul
##' @param name part of the file name that stays the same
##' @param x possibly a number or other indicator of the unique haul
##' @param extension what type of file is it
##' @param raising name of raising factor if there is one e.g. "RAISING_FACTOR"
##' @param sampling name of sampled fraction if there is one e.g. "SAMPLING"
##' @details the name of the file where the data is stored is paste0(name, x, extension)
##' @export
read_in_haul=function(x, name="Haul", extension=".txt", raising=NULL, sampling=NULL){
  y=read.table(paste0(name, x, extension), header=TRUE, fill=TRUE)
  lengths=subset(y, !is.na(as.numeric(as.character(y$LENGTH)))) #length data
  lengths$LENGTH=as.numeric(as.character(lengths$LENGTH))
  extra=subset(y, is.na(as.numeric(as.character(y$LENGTH))))#extra lines after length data

  #Add columns of raising or sampling
  if(!is.null(raising)) {
    RFindex= which(extra$LENGTH==raising)
    RF=extra[RFindex,]
    colnames(RF)=paste0("raising_",colnames(RF))
    for(i in 2:ncol(RF)) {lengths[,colnames(RF)[i]]=RF[i]}
    extra=extra[-RFindex,] #remove that row from extras
  } else
    if(!is.null(sampling)) {
      SFindex=which(extra$LENGTH==sampling)
      SF=extra[SFindex,]
      colnames(SF)=paste0("sampling_",colnames(SF))
      for(i in 2:ncol(SF)) {lengths[,colnames(SF)[i]]=SF[i]}
      extra=extra[-SFindex,] #remove that row from extras
    }

  #Add columns of covariates that are left in extra lines
  covs=subset(extra, !is.na(extra[,2]) )[,1:2] #covariates
  if(nrow(covs)>0)
  {
    for (i in 1:nrow(covs))
    {
      lengths[,as.character(covs[i,1])]=covs[i,2]
    }
  }

  lengths$haul=x

  return(lengths)
}
