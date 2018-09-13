##' Extract info from formulas, reTrms, etc., format for TMB
##' @param rformula selectivity formula
##' @param pformula relative fishing power formula
##' @param dformula Richards delta parameter formula
##' @param combForm combined formula
##' @param mf call to model frame
##' @param fr frame
##' @param yobs observed y
##' @param respCol response column
##' @param total total
##' @param link character
##' @param cover (logical) covered codend model (i.e. big fish go in experimental net and small fish go in covered codend)
##' @param Lp controls calculation of length (L) at retention prob (p)
##' @param pPredictCode relative fishing power code
##' @param doPredict flag to enable sds of predictions
##' @param whichPredict which observations in model frame represent predictions
##' @param x0 vector of initial values for the size selectivity model
##' @keywords internal
##' @importFrom stats model.offset
mkTMBStruc <- function(rformula, pformula, dformula,
                       combForm,
                       mf, fr,
                       yobs, total,
                       family, link_char, cover, Lp,
                       pPredictCode="selection",
                       doPredict=0,
                       whichPredict=integer(0),
                       x0=NULL,
                       call=NULL,
                       respCol) {

  mapArg <- NULL
  pformula.orig <- pformula ## Restore when done

  ## p=0.5 for equal fishing power in test and control codend
  if(pformula == ~0) {
    pformula[] <- ~1
    betap_init <- 0 #logit(.5) #p always has logit link
    mapArg <- c(mapArg, list(betap = factor(NA))) ## Fix betap
  }
  if(cover) {
    pformula[] <- ~0 #no p in cover models
  }

  rList  <- getXReTrms(rformula, mf, fr)
  pList  <- getXReTrms(pformula, mf, fr)
  dList  <- getXReTrms(dformula, mf, fr,
                                  ranOK=FALSE, "dispersion")

  rReStruc <- with(rList, getReStruc(reTrms, ss))
  pReStruc <- with(pList, getReStruc(reTrms, ss))

  grpVar <- with(rList, getGrpVar(reTrms$flist))

  nobs <- nrow(fr)

  if (is.null(total)) total <- rep(1,nobs)

  ## store info on location of response variable
  respCol <- attr(terms(fr), "response")
  names(respCol) <- names(fr)[respCol]

  #FLAGS
  #FIXME: might be 'width' instead of 'length'
        #does it ever need to output L50 for 2 or more types?
#  Lindex = grep("length", colnames(rList$X), ignore.case=TRUE)-1
  Lp <- match.arg(Lp, c("basic", "none", "full", "100"))
  Lpflag = switch(Lp, "basic"=1, "none"=0, "full"=2, "100"=3)


  data.tmb <- namedList(
    Xr = rList$X,
    Zr = rList$Z,
    Xp = pList$X,
    Zp = pList$Z,
    Xd = dList$X,
    yobs,
    respCol,
    offset = rList$offset,
    total,
    ## information about random effects structure
    termsr = rReStruc,
    termsp = pReStruc,
    link = .valid_link[link_char],
    pPredictCode = .valid_ppredictcode[pPredictCode],
    doPredict = doPredict,
    Lpflag = Lpflag,
    cover = as.numeric(cover),
    whichPredict = whichPredict
  )
  getVal <- function(obj, component)
    vapply(obj, function(x) x[[component]], numeric(1))

  if(is.null(x0)) {
    if(with(data.tmb,ncol(Xr))==2) {
      betar = with(data.tmb, c(interceptinit(link_char), .3))
    } else {
      betar = with(data.tmb, c(interceptinit(link_char), rep(0, ncol(Xr)-1)))
    }
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
  pformula <- pformula.orig ## May have changed - restore
  namedList(data.tmb, parameters, mapArg, randomArg, grpVar,
            rList, pList, dList, rReStruc, pReStruc,
            respCol,
            allForm=namedList(combForm,rformula,pformula,dformula),
            fr, call)
}

##' Initialize the intercept based on the link funciton
##' Assuming catchability of length 0 indivs is near 0
##' @param link character
interceptinit <- function(link) {
  r0 = 1e-10
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
##' \item{offset}{offset vector, or vector of zeros if offset not specified}
##'
##' @importFrom stats model.matrix
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
    terms_fixed <- terms(eval(mf,envir=environment(fixedform)))
    ## FIXME: make model matrix sparse?? i.e. Matrix:::sparse.model.matrix()
#    X <- model.matrix(drop.special2(fixedform), fr)
    X <- model.matrix(fixedform, fr)
    ## will be 0-column matrix if fixed formula is empty

    offset <- rep(0,nobs)
    terms <- list(fixed=terms(terms_fixed))
    if (inForm(fixedform,quote(offset))) {
      ## hate to match offset terms with model frame names
      ##  via deparse, but since that what was presumably done
      ##  internally to get the model frame names in the first place ...
      for (o in extractForm(fixedform,quote(offset))) {
        offset_nm <- deparse(o)
        ## don't think this will happen, but ...
        if (length(offset_nm)>1) {
          stop("trouble reconstructing offset name")
        }
        offset <- offset + fr[[offset_nm]]
      }
    }
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

  namedList(X, Z, reTrms, ss, terms, offset)
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

##' Fit gear selectivity models with TMB
##' @param rformula combined fixed and random effects formula for the selectivity model, following lme4
##'     syntax. The left-hand side of the formula should be the proportion of fish entering the test gear.
##' @param pformula a \emph{one-sided} (i.e., no response variable) formula for
##'     the ralaive fishing power of the test versus the control gear combining fixed and random effects:
##' \code{~0} can be used to specify equal fishing power (p=0.5).
##' The relative fishing power model uses a logit link.
##' @param link A character indicating the link function for the selectivity model.
##' \code{"logit"}(logistic) is the default, but other options are "probit" (i.e. normal probability ogiv), "cloglog" (i.e. negative extreme value), "loglog" (i.e. extreme value/Gompert), or "richards"
##' @param dformula a formula for the delta parameter in Richards selection curve. Ignored unless \code{link="richards"}.
##' @param cover (logical) covered codend model (i.e. big fish go in experimental net and small fish go in cover)
##' @param x0 vector of initial values for the size selectivity model
##' @param data data frame
##' @param total The number of total fish caught in the test and control gear.
##' @param haul Name of column representing different hauls. Needed for double bootstrap.
##' @param Lp controls calculation of length (L) at retention prob (p), see details
##' @param se whether to return standard errors
##' @param verbose logical indicating if some progress indication should be printed to the console.
##' @param debug whether to return the preprocessed data and parameter objects,
##'     without fitting the model
##' @param optControl control parameters passed to \code{nlminb}
##' @importFrom stats binomial nlminb as.formula terms model.weights
##' @importFrom lme4 subbars findbars mkReTrms nobars
##' @importFrom Matrix t
##' @importFrom TMB MakeADFun sdreport
##' @details
##' \itemize{
##' \item in all cases \code{selfisher} returns maximum likelihood estimates.
##' \item You only need to specify \code{haul} in models that are going to be bootstraped with type="double".
##' \item Lp="basic" will return values for L50 and SR.
##' \item Lp="none" supresses calculation of L50 and SR to save time.
##' \item Lp="full" will return values of Lp for p=5 to 95 as well as SR
##' \item Lp="100" will return values of Lp for p=1 to 100 as well as SR
##' \item Use \code{getCapabilities()} to see options for links and RE
##' }
##' @useDynLib selfisher
##' @importFrom stats update
##' @export
##' @examples
##' dat <- transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))
##' m0 <- selfisher(prop~Lengths, p=~0, total=tot, dat)
##' m1 <- selfisher(prop~Lengths, p=~1, total=tot, dat)
selfisher <- function (
    rformula,
    data = NULL,
    pformula = ~1,
    dformula = ~1,
    cover = FALSE,
    x0 = NULL,
    link = "logit",
    total=NULL,
    haul=NULL,
    offset=NULL,
    Lp="basic",
    se=TRUE,
    verbose=FALSE,
    debug=FALSE,
    optControl=list(iter.max=300, eval.max=400)
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
    ## add offset-specified-as-argument to formula as + offset(...)
    ## need evaluate offset within envi
    if (!is.null(eval(substitute(offset),data,
                      enclos=environment(rformula)))) {
      rformula <- addForm0(rformula,makeOp(substitute(offset),op=quote(offset)))
    }

    environment(pformula) <- environment(rformula)
    call$pformula <- pformula

    environment(dformula) <- environment(rformula)
    call$dformula <- dformula

    call$cover <- cover

    if(link!="richards") dformula[] <- ~0

    ## now work on evaluating model frame
    m <- match(c("data", "subset", "total", "haul", "offset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")

    ## want the model frame to contain the union of all variables
    ## used in any of the terms
    ## combine all formulas
    formList <- list(rformula, pformula, dformula)
    for (i in seq_along(formList)) {
      f <- formList[[i]] ## abbreviate
      ## substitute "|" by "+"; drop specials
      f <- noSpecials(subbars(f),delete=FALSE)
      formList[[i]] <- f
    }
    combForm <- do.call(addForm,formList)
    environment(combForm) <- environment(rformula)
    ## model.frame.default looks for these objects in the environment
    ## of the *formula* (see 'extras', which is anything passed in ...),
    ## so they have to be put there ...
#    for (i in c("weights", "offset")) {
#        if (!eval(bquote(missing(x=.(i)))))
#            assign(i, get(i, parent.frame()), environment(combForm))
#    }

    mf$formula <- combForm
    fr <- eval(mf,envir=environment(rformula),enclos=parent.frame())

    ## FIXME: throw an error *or* convert character to factor
    ## convert character vectors to factor (defensive)
    ## fr <- factorize(fr.form, fr, char.only = TRUE)
    ## store full, original formula
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

    TMBStruc <-
        mkTMBStruc(rformula=rformula, pformula=pformula, dformula=dformula,
                   combForm = combForm,
                   mf=mf, fr=fr,
                   yobs=y, total=total,
                   family=familyStr, link_char=link, cover=cover, x0=x0, Lp=Lp,
                   call=call, respCol=respCol)

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
                                                   gradient=gr,
                                                   control=optControl)))

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
                      namedList(nobs, respCol, grpVar, familyStr, family, link, cover,
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
    structure(namedList(obj, fit, sdr, call=TMBStruc$call,
                        frame=TMBStruc$fr, modelInfo,
                        fitted),
              class = "selfisher")
}

##' @importFrom stats AIC BIC
llikAIC <- function(object) {
    llik <- logLik(object)
    AICstats <-
        c(AIC = AIC(llik), BIC = BIC(llik), logLik = c(llik),
          deviance = sum(na.omit((residuals(object,type="deviance"))^2)),
          Pearson.ChiSq=sum((residuals(object,type="pearson"))^2),
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
            coefs <- cbind(coefs, 2*pnorm(abs(cf3), lower.tail = FALSE))
            colnames(coefs)[3:4] <- c(paste(statType, "value"),
                                      paste0("Pr(>|",statType,"|)"))
        }
        coefs
    }

    ff <- fixef(object)
    vv <- vcov(object)
    coefs <- setNames(lapply(names(ff),
            function(nm) if (trivialFixef(names(ff[[nm]]),nm) | is.null(vv[[nm]])) NULL else
                             mkCoeftab(ff[[nm]],vv[[nm]])),
                      names(ff))

    llAIC <- llikAIC(object)

    ## FIXME: You can't count on object@re@flist,
    ##	      nor compute VarCorr() unless is(re, "reTrms"):
    varcor <- VarCorr(object)
    #If the model is simple, extract Lp and SR
#    if(all(names(object$fit$par)=="betar")) {

    if("SR" %in% rownames(summary(object$sdr))) {
      SR <- summary(object$sdr, "report")["SR",]
      retention <- data.frame(p=object$obj$report()$retp,
            Lp.Est=summary(object$sdr, "report")[1:length(object$obj$report()$retp),1],
            Lp.Std.Err=summary(object$sdr, "report")[1:length(object$obj$report()$retp),2])
    } else {
      SR <- NULL
      retention <- NULL
    }
    structure(list(logLik = llAIC[["logLik"]],
                   family = famL, link = link,
                   ngrps = ngrps(object),
                   nobs = nobs(object),
                   coefficients = coefs, delta = sig,
                   vcov = vcov(object),
                   SR = SR,
                   retention = retention,
                   varcor = varcor, # and use formatVC(.) for printing.
                   AICtab = llAIC[["AICtab"]],
                   call = object$call
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
    .prt.retention(x$retention, x$SR)

    invisible(x)
}## print.summary.selfisher


