##' Extract info from formulas, reTrms, etc., format for TMB
##' @param rformula selectivity formula
##' @param pformula relative fishing power formula
##' @param dformula Richards delta parameter formula
##' @param mf call to model frame
##' @param fr frame
##' @param yobs observed y
##' @param offset offset
##' @param total total
##' @param link character
##' @param cc (logical) covered codend model (i.e. big fish go in experimental net and small fish go in covered codend)
##' @param pPredictCode relative fishing power code
##' @param doPredict flag to enable sds of predictions
##' @param whichPredict which observations in model frame represent predictions
##' @param x0 vector of initial values for the size selectivity model
##' @keywords internal
##' @importFrom stats model.offset
mkTMBStruc <- function(rformula, pformula, dformula,
                       mf, fr,
                       yobs, offset, total,
                       family, link_char, cc,
                       pPredictCode="selection",
                       doPredict=0,
                       whichPredict=integer(0),
                       x0=NULL) {

  mapArg <- NULL

  ## p=0.5 for equal fishing power in test and control codend
  if(pformula == ~0) {
    pformula[] <- ~1
    betap_init <- 0 #logit(.5) #p always has logit link
    mapArg <- c(mapArg, list(betap = factor(NA))) ## Fix betap
  }
  if(cc) {
  	pformula[] <- ~0 #no p in cc models
  }

  ## n.b. eval.parent() chain needs to be preserved because
  ## we are going to try to eval(mf) at the next level down,
  ## need to be able to find data etc.
  rList  <- eval.parent(getXReTrms(rformula, mf, fr))
  pList  <- eval.parent(getXReTrms(pformula, mf, fr))
  dList  <- eval.parent(getXReTrms(dformula, mf, fr,
                                        ranOK=FALSE, "dispersion"))

  rReStruc <- with(rList, getReStruc(reTrms, ss))
  pReStruc <- with(pList, getReStruc(reTrms, ss))

  grpVar <- with(rList, getGrpVar(reTrms$flist))

  nobs <- nrow(fr)
  ## FIXME: deal with offset in formula
  ##if (grepl("offset", safeDeparse(formula)))
  ##  stop("Offsets within formulas not implemented. Use argument.")

  if (is.null(offset <- model.offset(fr)))
      offset <- rep(0,nobs)

  if (is.null(total <- fr[["(total)"]]))
    total <- rep(1,nobs) #needed for predict function

  Lindex = grep("length", colnames(rList$X), ignore.case=TRUE)-1
  if(length(Lindex)!=1) Lindex = -1 #flag for complex function => no L50 or SR
  data.tmb <- namedList(
    Xr = rList$X,
    Zr = rList$Z,
    Xp = pList$X,
    Zp = pList$Z,
    Xd = dList$X,
    yobs,
    offset,
    total,
    ## information about random effects structure
    termsr = rReStruc,
    termsp = pReStruc,
    link = .valid_link[link_char],
    pPredictCode = .valid_ppredictcode[pPredictCode],
    doPredict = doPredict,
    Lindex = Lindex,
    cc = as.numeric(cc),
    whichPredict = whichPredict
  )
  getVal <- function(obj, component)
    vapply(obj, function(x) x[[component]], numeric(1))

  if(is.null(x0)) {
    betar    = with(data.tmb, c(interceptinit(link_char), rep(0, ncol(Xr)-1)))
  } else {
    betar = x0
  }
  parameters <- with(data.tmb,
                     list(
                       betar    = betar,
                       br       = rep(0, ncol(Zr)),
                       betap    = rep(0, ncol(Xp)),
                       bp       = rep(0, ncol(Zp)),
                       thetar   = rep(0, sum(getVal(rReStruc,"blockNumTheta"))),
                       thetap   = rep(0, sum(getVal(pReStruc,  "blockNumTheta"))),
                       betad    = rep(0, ncol(Xd))# d=1 makes richards become logisitc
                    ))
  randomArg <- c(if(ncol(data.tmb$Zr) > 0) "br",
                 if(ncol(data.tmb$Zp) > 0) "bp")
  namedList(data.tmb, parameters, mapArg, randomArg, grpVar,
            rList, pList, dList, rReStruc, pReStruc)
}

##' Initialize the intercept based on the link funciton
##' Assuming catchability of length 0 indivs is near 0
##' @param link character
interceptinit <- function(link) {
  r0 = 1e-12
  switch(link,
         "logit"    = log(r0/(1-r0)),
         "probit"   = qnorm(r0),
         "cloglog"  = log(-log(1-r0)),
         "loglog"   = -log(-log(r0)),
         "richards" = log(r0/(1-r0)))
}
##' Create X and random effect terms from formula
##' @param formula current formula, containing both fixed & random effects
##' @param mf matched call
##' @param fr full model frame
##' @param ranOK random effects allowed here?
##' @param type label for model type
##' @return a list composed of
##' \item{X}{design matrix for fixed effects}
##' \item{Z}{design matrix for random effects}
##' \item{reTrms}{output from \code{\link{mkReTrms}} from \pkg{lme4}}
##'
##' @importFrom stats model.matrix contrasts
##' @importFrom methods new
##' @importFrom lme4 findbars nobars
getXReTrms <- function(formula, mf, fr, ranOK=TRUE, type="") {
    ## fixed-effects model matrix X -
    ## remove random effect parts from formula:
    fixedform <- formula
    RHSForm(fixedform) <- nobars(RHSForm(fixedform))

    nobs <- nrow(fr)
    ## check for empty fixed form

    if (identical(RHSForm(fixedform), ~  0) ||
        identical(RHSForm(fixedform), ~ -1)) {
        X <- NULL
    } else {
        mf$formula <- fixedform

        terms_fixed <- terms(eval.parent(mf))
        
        ## FIXME: make model matrix sparse?? i.e. Matrix:::sparse.model.matrix()
        X <- model.matrix(fixedform, fr, contrasts)
        ## will be 0-column matrix if fixed formula is empty

        terms <- list(fixed=terms(terms_fixed))

    }
    ## ran-effects model frame (for predvars)
    ## important to COPY formula (and its environment)?
    ranform <- formula
    if (is.null(findbars(ranform))) {
        reTrms <- NULL
        Z <- new("dgCMatrix",Dim=c(as.integer(nobs),0L)) ## matrix(0, ncol=0, nrow=nobs)
        ss <- integer(0)
    } else {

        ## FIXME: check whether predvars are carried along correctly in terms
        if (!ranOK) stop("no random effects allowed in ", type, " term")
        RHSForm(ranform) <- subbars(RHSForm(reOnly(formula)))

        mf$formula <- ranform
        reTrms <- mkReTrms(findbars(RHSForm(formula)), fr)

        ss <- splitForm(formula)
        ss <- unlist(ss$reTrmClasses)

        Z <- t(reTrms$Zt)   ## still sparse ...
    }

    ## if(is.null(rankX.chk <- control[["check.rankX"]]))
    ## rankX.chk <- eval(formals(lmerControl)[["check.rankX"]])[[1]]
    ## X <- chkRank.drop.cols(X, kind=rankX.chk, tol = 1e-7)
    ## if(is.null(scaleX.chk <- control[["check.scaleX"]]))
    ##     scaleX.chk <- eval(formals(lmerControl)[["check.scaleX"]])[[1]]
    ## X <- checkScaleX(X, kind=scaleX.chk)

    ## list(fr = fr, X = X, reTrms = reTrms, family = family, formula = formula,
    ##      wmsgs = c(Nlev = wmsgNlev, Zdims = wmsgZdims, Zrank = wmsgZrank))

    namedList(X, Z, reTrms, ss, terms)
}

##' Extract grouping variables for random effect terms from a factor list
##' @title Get Grouping Variable
##' @param x "flist" object; a data frame of factors including an \code{assign} attribute
##' matching columns to random effect terms
##' @return character vector of grouping variables
##' @keywords internal
##' @examples
##' data(cbpp,package="lme4")
##' cbpp$obs <- factor(seq(nrow(cbpp)))
##' rt <- lme4::glFormula(cbind(size,incidence-size)~(1|herd)+(1|obs),
##'   data=cbpp,family=binomial)$reTrms
##' getGrpVar(rt$flist)
##' @export
getGrpVar <- function(x)
{
  assign <- attr(x,"assign")
  names(x)[assign]
}

##' Calculate random effect structure
##' Calculates number of random effects, number of parameters,
##' blocksize and number of blocks.  Mostly for internal use.
##' @param reTrms random-effects terms list
##' @param ss a character string indicating a valid covariance structure. 
##' Must be one of \code{names(selfisher:::.valid_covstruct)};
##' default is to use an unstructured  variance-covariance
##' matrix (\code{"us"}) for all blocks).
##' @return a list
##' \item{blockNumTheta}{number of variance covariance parameters per term}
##' \item{blockSize}{size (dimension) of one block}
##' \item{blockReps}{number of times the blocks are repeated (levels)}
##' \item{covCode}{structure code}
##' @examples
##' data(sleepstudy, package="lme4")
##' rt <- lme4::lFormula(Reaction~Days+(1|Subject)+(0+Days|Subject),
##'                     sleepstudy)$reTrms
##' rt2 <- lme4::lFormula(Reaction~Days+(Days|Subject),
##'                     sleepstudy)$reTrms
##' getReStruc(rt)
##' @importFrom stats setNames dist
##' @export
getReStruc <- function(reTrms, ss=NULL) {

  ## information from ReTrms is contained in cnms, flist
  ## cnms: list of column-name vectors per term
  ## flist: data frame of grouping variables (factors)
  ##   'assign' attribute gives match between RE terms and factors
    if (is.null(reTrms)) {
        list()
    } else {
        ## Get info on sizes of RE components

        assign <- attr(reTrms$flist,"assign")
        nreps <- vapply(assign,
                          function(i) length(levels(reTrms$flist[[i]])),
                          0)
        blksize <- diff(reTrms$Gp) / nreps
        ## figure out number of parameters from block size + structure type

        if (is.null(ss)) {
            ss <- rep("us",length(blksize))
        }

        covCode <- .valid_covstruct[ss]

        parFun <- function(struc, blksize) {
            switch(as.character(struc),
                   "0" = blksize, # diag
                   "1" = blksize * (blksize+1) / 2, # us
                   "2" = blksize + 1, # cs
                   "3" = 2,  # ar1
                   "4" = 2,  # ou
                   "5" = 2,  # exp
                   "6" = 2,  # gau
                   "7" = 3,  # mat
                   "8" = 2 * blksize - 1) # toep
        }
        blockNumTheta <- mapply(parFun, covCode, blksize, SIMPLIFY=FALSE)

        ans <-
            lapply( seq_along(ss), function(i) {
                tmp <-
                    list(blockReps = nreps[i],
                         blockSize = blksize[i],
                         blockNumTheta = blockNumTheta[[i]],
                         blockCode = covCode[i]
                         )
                if(ss[i] == "ar1"){
                    ## FIXME: Keep this warning ?
                    if (any(reTrms$cnms[[i]][1] == "(Intercept)") )
                        warning("AR1 not meaningful with intercept")
                }
                if(ss[i] == "ou"){
                    times <- parseNumLevels(reTrms$cnms[[i]])
                    if (ncol(times) != 1)
                        stop("'ou' structure is for 1D coordinates only.")
                    if (is.unsorted(times, strictly=TRUE))
                        stop("'ou' is for strictly sorted times only.")
                    tmp$times <- drop(times)
                }
                if(ss[i] %in% c("exp", "gau", "mat")){
                    coords <- parseNumLevels(reTrms$cnms[[i]])
                    tmp$dist <- as.matrix( dist(coords) )
                }
                tmp
            })
        setNames(ans, names(reTrms$Ztlist))
    }
}

## select only desired pieces from results of getXReTrms
stripReTrms <- function(xrt, whichReTrms = c("cnms","flist"), which="terms") {
  c(xrt$reTrms[whichReTrms],setNames(xrt[which],which))
}

##' Fit models with TMB
##' @param rformula combined fixed and random effects formula for the selectivity model, following lme4
##'     syntax. The left-hand side of the formula should be the proportion of fish entering the test gear.
##' @param pformula a \emph{one-sided} (i.e., no response variable) formula for
##'     the ralaive fishing power of the test versus the control gear combining fixed and random effects:
##' \code{~0} can be used to specify equal fishing power (p=0.5).
##' The relative fishing power model uses a logit link.
##' @param link A character indicating the link function for the selectivity model. 
##' \code{"logit"} is the default, but other options can be used (use \code{getCapabilities()} to see options).
##' @param dformula a formula for the delta parameter in Richards selection curve. Ignored unless \code{link="richards"}.
##' @param cc (logical) covered codend model (i.e. big fish go in experimental net and small fish go in cover)
##' @param x0 vector of initial values for the size selectivity model
##' @param data data frame
##' @param total The number of total fish caught in the test and control gear.
##' @param offset offset
##' @param se whether to return standard errors
##' @param verbose logical indicating if some progress indication should be printed to the console.
##' @param debug whether to return the preprocessed data and parameter objects,
##'     without fitting the model
##' @importFrom stats binomial nlminb as.formula terms model.weights
##' @importFrom lme4 subbars findbars mkReTrms nobars
##' @importFrom Matrix t
##' @importFrom TMB MakeADFun sdreport
##' @details
##' \itemize{
##' \item in all cases \code{selfisher} returns maximum likelihood estimates - random effects variance-covariance matrices are not REML (so use \code{REML=FALSE} when comparing with \code{lme4::lmer}), and residual standard deviations (\code{\link{sigma}}) are not bias-corrected. Because the \code{\link{df.residual}} method for \code{selfisher} currently counts the dispersion parameter, one would need to multiply by \code{sqrt(nobs(fit)/(1+df.residual(fit)))} when comparing with \code{lm} ...
##' }
##' @useDynLib selfisher
##' @importFrom stats update
##' @export
##' @examples
selfisher <- function (
    rformula,
    pformula = ~1,
    dformula = ~1,
    cc = FALSE,
    x0 = NULL,
    data = NULL,
    link = "logit",
    total=NULL,
    offset=NULL,
    se=TRUE,
    verbose=FALSE,
    debug=FALSE
    )
{
    ## FIXME: check for offsets in pformula/dispformula, throw an error

    call <- mf <- mc <- match.call()

    family <- binomial()
    familyStr <- "binomial"

    ## lme4 function for warning about unused arguments in ...
    ## ignoreArgs <- c("start","verbose","devFunOnly",
    ##   "optimizer", "control", "nAGQ")
    ## l... <- list(...)
    ## l... <- l...[!names(l...) %in% ignoreArgs]
    ## do.call(checkArgs, c(list("glmer"), l...))

    # substitute evaluated versions
    ## FIXME: denv leftover from lme4, not defined yet

    environment(rformula) <- parent.frame()
    call$rformula <- mc$rformula <- rformula

    environment(pformula) <- environment(rformula)
    call$pformula <- pformula

    environment(dformula) <- environment(rformula)
    call$dformula <- dformula

    if(link!="richards") dformula = ~0

    ## now work on evaluating model frame
    m <- match(c("data", "subset", "total", "na.action", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")

    ## want the model frame to contain the union of all variables
    ## used in any of the terms
    ## combine all formulas
    formList <- list(rformula, pformula)
    formList <- lapply(formList,
                   function(x) noSpecials(subbars(x), delete=FALSE))
                       ## substitute "|" by "+"; drop special
    combForm <- do.call(addForm,formList)
    environment(combForm) <- environment(rformula)
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in ...),
    ## so they have to be put there ...
 #   for (i in c("total", "offset")) {
	  for (i in c("offset")) {
        if (!eval(bquote(missing(x=.(i)))))
            assign(i, get(i, parent.frame()), environment(combForm))
    }

    mf$formula <- combForm
    fr <- eval(mf,envir=environment(formula),enclos=parent.frame())
    
    ## FIXME: throw an error *or* convert character to factor
    ## convert character vectors to factor (defensive)
    ## fr <- factorize(fr.form, fr, char.only = TRUE)
    ## store full, original formula & offset
    ## attr(fr,"formula") <- combForm  ## unnecessary?
    nobs <- nrow(fr)
    total <- as.vector(model.total(fr))

    if(is.null(total)) {
      stop("The total number of fish caught in the test and control gear must be specified using 'total' argument.")
    }
    
    ## sanity checks (skipped!)
    ## wmsgNlev <- checkNlevels(reTrms$ flist, n=n, control, allow.n=TRUE)
    ## wmsgZdims <- checkZdims(reTrms$Ztlist, n=n, control, allow.n=TRUE)
    ## wmsgZrank <- checkZrank(reTrms$Zt, n=n, control, nonSmall=1e6, allow.n=TRUE)

    ## store info on location of response variable
    respCol <- attr(terms(fr), "response")
    names(respCol) <- names(fr)[respCol]
     
    ## extract response variable
    ## (name *must* be 'y' to match guts of family()$initialize
    y <- fr[,respCol]

    ## eval.parent() necessary because we will try to eval(mf) down below
    TMBStruc <- eval.parent(
        mkTMBStruc(rformula, pformula, dformula,
                   mf, fr,
                   yobs=y, offset, total,
                   family=familyStr, link_char=link, cc=cc, x0=x0))

    ## short-circuit
    if(debug) return(TMBStruc)

    obj <- with(TMBStruc,
                MakeADFun(data.tmb,
                     parameters,
                     map = mapArg,
                     random = randomArg,
                     profile = NULL,
                     silent = !verbose,
                     DLL = "selfisher"))

    optTime <- system.time(fit <- with(obj, nlminb(start=par, objective=fn,
                                                   gradient=gr)))

    fit$parfull <- obj$env$last.par.best ## This is in sync with fit$par

    fitted <- NULL

    if (se) {
        sdr <- sdreport(obj)
        ## FIXME: assign original rownames to fitted?
    } else {
        sdr <- NULL
    }
    if(!is.null(sdr$pdHess)) {
      if(!sdr$pdHess) {
        warning(paste0("Model convergence problem; ",
                       "non-positive-definite Hessian matrix. ", 
                       "See vignette('troubleshooting')"))
      } else {
        eigval <- try(1/eigen(sdr$cov.fixed)$values, silent=TRUE)
        if( is(eigval, "try-error") || ( min(eigval) < .Machine$double.eps*10 ) ) {
          warning(paste0("Model convergence problem; ",
                       "extreme or very small eigen values detected. ", 
                       "See vignette('troubleshooting')"))
        }
      }
    }

    modelInfo <- with(TMBStruc,
                      namedList(nobs, respCol, grpVar, familyStr, family, link, cc,
                                ## FIXME:apply condList -> cond earlier?
                                reTrms = lapply(list(r=rList, p=pList),
                                                stripReTrms),
                                reStruc = namedList(rReStruc, pReStruc),
                                allForm = namedList(combForm, rformula,
                                                    pformula, dformula)))
    ## FIXME: are we including obj and frame or not?
    ##  may want model= argument as in lm() to exclude big stuff from the fit
    ## If we don't include obj we need to get the basic info out
    ##    and provide a way to regenerate it as necessary
    ## If we don't include frame, then we may have difficulty
    ##    with predict() in its current form
    structure(namedList(obj, fit, sdr, call, frame=fr, modelInfo,
                        fitted),
              class = "selfisher")
}

##' @importFrom stats AIC BIC
llikAIC <- function(object) {
    llik <- logLik(object)
    AICstats <- 
        c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
          deviance = -2*llik, ## FIXME:
          df.resid = df.residual(object))
    list(logLik = llik, AICtab = AICstats)
}

## FIXME: export/import from lme4?
ngrps <- function(object, ...) UseMethod("ngrps")

ngrps.default <- function(object, ...) stop("Cannot extract the number of groups from this object")

ngrps.selfisher <- function(object, ...) {
    res <- lapply(object$modelInfo$reTrms,
           function(x) vapply(x$flist, nlevels, 1))
    ## FIXME: adjust reTrms names for consistency rather than hacking here
    names(res) <- gsub("List$","",names(res))
    return(res)
    
}

ngrps.factor <- function(object, ...) nlevels(object)


##' @importFrom stats pnorm
##' @method summary selfisher
##' @export
summary.selfisher <- function(object,...)
{
    if (length(list(...)) > 0) {
        ## FIXME: need testing code
        warning("additional arguments ignored")
    }
    ## figure out useSc
    sig <- richardsdelta(object)

    famL <- family(object)$family
    link <- object$modelInfo$link

    mkCoeftab <- function(coefs,vcov) {
        p <- length(coefs)
        coefs <- cbind("Estimate" = coefs,
                       "Std. Error" = sqrt(diag(vcov)))
        if (p > 0) {
            coefs <- cbind(coefs, (cf3 <- coefs[,1]/coefs[,2]),
                           deparse.level = 0)
            ## statType <- if (useSc) "t" else "z"
            statType <- "z"
            ## ??? should we provide Wald p-values???
            coefs <- cbind(coefs, 2*pnorm(abs(cf3), lower.tail = FALSE))
            colnames(coefs)[3:4] <- c(paste(statType, "value"),
                                      paste0("Pr(>|",statType,"|)"))
        }
        coefs
    }

    ff <- fixef(object)
    vv <- vcov(object)
    coefs <- setNames(lapply(names(ff),
            function(nm) if (trivialFixef(names(ff[[nm]]),nm)) NULL else
                             mkCoeftab(ff[[nm]],vv[[nm]])),
                      names(ff))

    llAIC <- llikAIC(object)
                   
    ## FIXME: You can't count on object@re@flist,
    ##	      nor compute VarCorr() unless is(re, "reTrms"):
    varcor <- VarCorr(object)
					# use S3 class for now
    structure(list(logLik = llAIC[["logLik"]],
                   family = famL, link = link,
		   ngrps = ngrps(object),
                   nobs = nobs(object),
		   coefficients = coefs, delta = sig,
		   vcov = vcov(object),
		   varcor = varcor, # and use formatVC(.) for printing.
		   AICtab = llAIC[["AICtab"]], call = object$call
                   ## residuals = residuals(object,"pearson",scaled = TRUE),
		   ## fitMsgs = .merMod.msgs(object),
                   ## optinfo = object@optinfo
		   ), class = "summary.selfisher")
               
}

## copied from lme4:::print.summary.merMod (makes use of
##' @importFrom lme4 .prt.family .prt.resids .prt.VC .prt.grps
##' @importFrom stats printCoefmat
##' @method print summary.selfisher
##' @export
print.summary.selfisher <- function(x, digits = max(3, getOption("digits") - 3),
                                 signif.stars = getOption("show.signif.stars"),
                                 ranef.comp = c("Variance", "Std.Dev."),
                                 show.resids = FALSE, ...)
{
    .prt.family(x)
    .prt.call.selfisher(x$call); cat("\n")
    .prt.aictab(x$AICtab); cat("\n")
    if (show.resids)
        .prt.resids(x$residuals, digits = digits)

    if (any(whichRE <- !sapply(x$varcor,is.null))) {
        cat("Random effects:\n")
        for (nn in names(x$varcor[whichRE])) {
            cat("\n",cNames[[nn]],":\n",sep="")
            ## lme4:::.prt.VC is not quite what we want here
            print(formatVC(x$varcor[[nn]],
                           digits = digits,
                           comp = ranef.comp),
                  quote=FALSE, digits=digits)
            ## FIXME: redundant nobs output
            .prt.grps(x$ngrps[[nn]],nobs=x$nobs)
        }
    }
   if((x$call$dformula== ~1)&(x$link=="richards")) {# if trivial print here, else below(~x) or none(~0)
    printDispersion(x$delta)  
   }
    for (nn in names(x$coefficients)) {
        cc <- x$coefficients[[nn]]
        p <- length(cc)
        if (p > 0) {
            cat("\n",cNames[[nn]],":\n",sep="")
            printCoefmat(cc, zap.ind = 3, #, tst.ind = 4
                         digits = digits, signif.stars = signif.stars)
        } ## if (p>0)
    }

    invisible(x)
}## print.summary.selfisher


