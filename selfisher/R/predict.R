##' prediction
##' @param object a \code{selfisher} object
##' @param newdata new data for prediction
##' @param se.fit return the standard errors of the predicted values?
##' @param type
##' return expected value ("response": pr/(pr+1-p)),
##' the  predicted selection curve ("selection": r),
##' or the relative fishing power of the test gear ("prob")?
##' @param debug (logical) return the \code{TMBStruc} object that will be
##' used internally for debugging?
##' @param re.form (not yet implemented) specify which random effects to condition on when predicting
##' @param allow.new.levels (not yet implemented) allow previously unobserved levels in random-effects grouping variables?
##' @param \dots unused - for method compatibility

##' @examples data(haddock)
##' dat <- transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))
##' m1 <- selfisher(prop~Lengths, p=~1, total=tot, dat)
##' nd <- data.frame(Lengths=20:50, tot=100)
##' predict(m1, newdata=nd, se.fit=TRUE)
##' @importFrom TMB sdreport
##' @importFrom stats optimHess
##' @export
predict.selfisher <- function(object,newdata=NULL,
                            se.fit=FALSE,
                            re.form, allow.new.levels=FALSE,
                            type = c("response","selection","prob"),
                            debug=FALSE,
                            ...)
{
  ## FIXME: add re.form
  ## FIXME: deal with napredict stuff ...

  type <- match.arg(type)
  if (!missing(re.form)) stop("re.form not yet implemented")
  ##if (allow.new.levels) stop("allow.new.levels not yet implemented")
  mc <- mf <- object$call
  ## FIXME: DRY so much
  ## now work on evaluating model frame
  ## do we want to re-do this part???

  ## need to 'fix' call to proper model.frame call whether or not
  ## we have new data, because ... (??)
  m <- match(c("subset", "total","haul", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]

  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf$formula <- RHSForm(object$modelInfo$allForm$combForm, as.form=TRUE)

  if (is.null(newdata)) {
    mf$data <- mc$data ## restore original data
    newFr <- object$fr
  } else {
    mf$data <- newdata
    newFr <- eval.parent(mf)
  }

  omi <- object$modelInfo  ## shorthand ("**o**bject$**m**odel**I**nfo")

  respCol <- match(respNm <- names(omi$respCol),names(newFr))
  ## create *or* overwrite response column for prediction data with NA
  newFr[[respNm]] <- NA

  ## FIXME: not yet handling population-level predictions (re.form
  ##  or new levels/allow.new.levels)

  ## append to existing model frame
  augFr <- rbind(object$fr,newFr)

  ## Pointers into 'new rows' of augmented data frame.
  w <- nrow(object$fr) + seq_len(nrow(newFr))

  ## Variety of possible binomial inputs are taken care of by
  ## 'mkTMBStruc' further down.
  yobs <- augFr[[names(omi$respCol)]]


  ## match zitype arg with internal name
  PredNm <- switch(type,
                       response="response",
                       selection="selection",
                       prob="prob",
                       stop("unknown type ",type))
  PredCode <- .valid_ppredictcode[PredNm]

  ## need eval.parent() because we will do eval(mf) down below ...
  TMBStruc <-
        ## with() interfering with eval.parent() ?
        eval.parent(mkTMBStruc(rformula=RHSForm(omi$allForm$rformula,as.form=TRUE),
                               pformula=omi$allForm$pformula,
                               dformula=omi$allForm$dformula,
                               combForm=omi$allForm$combForm,
                               mf=mf,
                               fr=augFr,
                               yobs=yobs,
                               total=model.total(augFr),
                               family=omi$familyStr,
                               link=omi$link,
                               cover=omi$cover,
                               Lp=omi$Lp,
                               pPredictCode=PredNm,
                               doPredict=as.integer(se.fit),
                               whichPredict=w))

  ## short-circuit
  if(debug) return(TMBStruc)

  newObj <- with(TMBStruc,
                 MakeADFun(data.tmb,
                           parameters,
                           map = mapArg,
                           random = randomArg,
                           profile = NULL, # TODO: Optionally "beta"
                           silent = TRUE,
                           DLL = "selfisher"))

  oldPar <- object$fit$par
  newObj$fn(oldPar)  ## call once to update internal structures
  lp <- newObj$env$last.par

  if (!se.fit) {
    pred <- newObj$report(lp)$mu_predict
  } else {
    H <- with(object,optimHess(oldPar,obj$fn,obj$gr))
    ## FIXME: Eventually add 'getReportCovariance=FALSE' to this sdreport
    ##        call to fix memory issue (requires recent TMB version)
    ## Fixed! (but do we want a flag to get it ? ...)
    sdr <- sdreport(newObj,oldPar,hessian.fixed=H,getReportCovariance=FALSE)
    sdrsum <- summary(sdr, "report") ## TMB:::summary.sdreport(sdr, "report")
    pred <- sdrsum[,"Estimate"]
    se <- sdrsum[,"Std. Error"]
  }
  if (!se.fit) return(pred) else return(list(fit=pred,se.fit=se))
}
