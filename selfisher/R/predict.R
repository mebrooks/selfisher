## Helper function for predict.
## Assert that we can use old model (data.tmb0) as basis for
## predictions using the new data (data.tmb1):
assertIdenticalModels <- function(data.tmb1, data.tmb0, allow.new.levels=FALSE)
{
  ## Check terms. Only 'blockReps' and 'blockSize' are allowed to
  ## change.  Note that we allow e.g. spatial covariance matrices to
  ## change, while e.g. an unstrucured covariance must remain the
  ## same.
  checkTerms <- function(t1, t0) {
    ## Defensive check:
    stopifnot(identical(names(t1), names(t0)))
    ## *Never* allowed to differ:
    testIdentical <- function(checkNm) {
      unlist( Map( function(x,y)
        identical(x[checkNm], y[checkNm]), t0, t1) )
    }
    ok <- testIdentical( c("blockNumTheta", "blockCode") )
    if ( ! all(ok) ) {
      msg <- c("Prediction is not possible for terms: ",
               paste(names(t1)[!ok], collapse=", "), "\n",
               "Probably some factor levels in 'newdata' require fitting a new model.")
      stop(msg)
    }
    ## Sometimes allowed to differ:
    if ( ! allow.new.levels ) {
      ok <- testIdentical( c( "blockReps", "blockSize") )
      if ( ! all(ok) ) {
        msg <- c("Predicting new random effect levels for terms: ",
                 paste(names(t1)[!ok], collapse=", "), "\n",
                 "Disable this warning with 'allow.new.levels=TRUE'")
        ## FIXME: warning or error ?
        warning(msg)
      }
    }
  }
  checkTerms( data.tmb1$termsr,   data.tmb0$termsr )
  checkTerms( data.tmb1$termsp, data.tmb0$termsp )
  ## Fixed effect parameters must be identical
  checkModelMatrix <- function(X1, X0) {
    if( !identical(colnames(X1), colnames(X0)) ) {
      msg <- c("Prediction is not possible for unknown fixed effects: ",
               paste( setdiff(colnames(X1), colnames(X0)), collapse=", "), "\n",
               "Probably some factor levels in 'newdata' require fitting a new model.")
      stop(msg)
    }
  }
  checkModelMatrix(data.tmb1$Xr, data.tmb0$Xr)
  checkModelMatrix(data.tmb1$Xp, data.tmb0$Xp)
  NULL
}

##' prediction
##' @param object a \code{selfisher} object
##' @param newdata new data for prediction
##' @param se.fit return the standard errors of the predicted values?
##' @param type
##' \itemize{
##' \item return expected response value ("response": see details below),
##' \item predicted selection curve ("selection": r),
##' \item relative fishing power of the test gear ("prob"),
##' \item catch ratio from catch comparison models ("ratio": r/(1-r)).
##' }
##' Some types (ratio, link) might not make sense to use with psplit models.
##' @param debug (logical) return the \code{TMBStruc} object that will be
##' used internally for debugging?
##' @param re.form (not yet implemented) specify which random effects to condition on when predicting. For now, all random effects are included.
##' @param allow.new.levels (not yet implemented) allow previously unobserved levels in random-effects grouping variables?
##' @param \dots unused - for method compatibility
##' @details Predicting with type="response" returns values comparable to the
##' response variable (the left-hand side of the model's rformula);
##' that is  pr/(pr+1-p) in a model with psplit=TRUE
##' or r in a model with psplit=FALSE.
##' This function could work with random effects, but is untested.
##' @examples data(haddock)
##' dat <- transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))
##' m1 <- selfisher(prop~Lengths, p=~1, psplit=TRUE, total=tot, dat)
##' nd <- data.frame(Lengths=20:50, tot=100)
##' predict(m1, newdata=nd, se.fit=TRUE)
##' @importFrom TMB sdreport
##' @importFrom stats optimHess
##' @export
predict.selfisher <- function(object,newdata=NULL,
                            se.fit=FALSE,
                            re.form, allow.new.levels=FALSE,
                            type = c("response","selection","prob", "ratio", "link"),
                            na.action = na.pass,
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
    mf$na.action <- na.action
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
                       ratio="ratio",
                       link="link",
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
                               psplit=omi$psplit,
                               Lp=omi$Lp,
                               pPredictCode=PredNm,
                               doPredict=as.integer(se.fit),
                               whichPredict=w))

  ## short-circuit
  if(debug) return(TMBStruc)

  ## Check that the model specification is unchanged:
  assertIdenticalModels(TMBStruc$data.tmb,
                        object$obj$env$data, allow.new.levels)

  ## Check that the neccessary predictor variables are finite (not NA nor NaN)
  if(se.fit) {
    with(TMBStruc$data.tmb, if(any(!is.finite(Xr)) |
                               any(!is.finite(Zr@x)) |
                               any(!is.finite(Xp)) |
                               any(!is.finite(Zp@x)) |
                               any(!is.finite(Xd))
    ) stop("Some variables in newdata needed for predictions contain NAs or NaNs.
           This is currently incompatible with se.fit=TRUE."))
  }

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
