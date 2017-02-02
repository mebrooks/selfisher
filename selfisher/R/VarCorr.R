## returns binomial()
##' @importFrom stats family
##' @export
family.selfisher <- function(object, ...) {
    object$modelInfo$family
}
## returns the link used for the selectivity model
##' @export
link.selfisher <- function(object, ...) {
    object$modelInfo$link
}

## don't quite know why this (rather than just ...$parList()) is
## necessary -- used in ranef.selfisher and sigma.selfisher
getParList <- function(object) {
    object$obj$env$parList(object$fit$par, object$fit$parfull)
}


##' Extract Richards exponent parameter
##'
##' @param object a \dQuote{selfisher} fitted object
##' @param \dots (ignored; for method compatibility)
##' @export delta
##' @method delta selfisher
##' @export
delta.selfisher <- function(object, ...) {
    pl <- getParList(object)
    ff <- object$modelInfo$link
    if (ff=="richards") return(exp(pl$betad))
}


mkVC <- function(cor, sd, cnms, sc, useSc) {
    stopifnot(length(cnms) == (nc <- length(cor)),  nc == length(sd),
              is.list(cnms), is.list(cor), is.list(sd),
              is.character(nnms <- names(cnms)), nzchar(nnms))
    ##
    ## FIXME: do we want this?  Maybe not.
    ## Potential problem: the names of the elements of the VarCorr() list
    ##  are not necessarily unique (e.g. fm2 from example("lmer") has *two*
    ##  Subject terms, so the names are "Subject", "Subject".  The print method
    ##  for VarCorrs handles this just fine, but it's a little awkward if we
    ##  want to dig out elements of the VarCorr list ... ???
    if (anyDuplicated(nnms))
        nnms <- make.names(nnms, unique = TRUE)
    ##
    ## cov :=  F(sd, cor) :
    do1cov <- function(sd, cor, n = length(sd)) {
        sd * cor * rep(sd, each = n)
    }
    docov <- function(sd,cor,nm) {
        ## diagonal model:
        diagmodel <- identical(dim(cor),c(0L,0L))
        if (diagmodel) cor <- diag(length(sd))
        cov <- do1cov(sd, cor)
        names(sd) <- nm
        dimnames(cov) <- dimnames(cor) <- list(nm,nm)
        structure(cov,stddev=sd,correlation=cor)
    }
    ss <- setNames(mapply(docov,sd,cor,cnms,SIMPLIFY=FALSE),nnms)
    ## ONLY first element -- otherwise breaks formatVC
    ## FIXME: do we want a message/warning here, or elsewhere,
    ##   when the 'Residual' var parameters are truncated?
    attr(ss,"sc") <- sc[1]
    attr(ss,"useSc") <- useSc
    ss
}


##' Extract variance and correlation components
##'
##' @aliases VarCorr
##' @param x a fitted \code{selfisher} model
##' @param sigma residual standard deviation (usually set automatically from internal information)
##' @param extra arguments (for consistency with generic method)
##' @importFrom nlme VarCorr
## and re-export the generic:
##' @export VarCorr
##' @export
##' @examples
##' ## Comparing variance-covariance matrix with manual computation
##' @keywords internal
VarCorr.selfisher <- function(x, sigma = 1, ... )
{
    ## FIXME:: add type=c("varcov","sdcorr","logs" ?)
    ## FIXME:: do we need 'sigma' any more (now that nlme generic
    ##         doesn't have it?)
    stopifnot(is.numeric(sigma), length(sigma) == 1)
    xrep <- x$obj$env$report(x$fit$parfull)
    reT <- x$modelInfo$reTrms
    familyStr <- family(x)$family
    useSc <- if (missing(sigma)) {
        sigma <- sigma(x)
        familyStr=="gaussian"
        ## *only* report residual variance for Gaussian family ...
        ## usesDispersion(familyStr)
    } else TRUE

    vc.cond <- if(length(cn <- reT$cond$cnms))
        mkVC(cor = xrep$corr,  sd = xrep$sd,   cnms = cn,
             sc = sigma, useSc = useSc)
    vc.zi   <- if(length(cn <- reT$zi$cnms))
        mkVC(cor = xrep$corrzi, sd = xrep$sdzi, cnms = cn,
             sc = sigma, useSc = useSc)
    structure(list(cond = vc.cond, zi = vc.zi),
	      sc = usesDispersion(familyStr), ## 'useScale'
	      class = "VarCorr.selfisher")
}

##' Printing The Variance and Correlation Parameters of a \code{selfisher}
##' @method print VarCorr.selfisher
##' @export
##' @importFrom lme4 formatVC
##  document as it is a method with "surprising arguments":
##' @param x a result of \code{\link{VarCorr}(<selfisher>)}.
##' @param digits number of significant digits to use.
##' @param comp a string specifying the component to format and print.
##' @param formatter a \code{\link{function}}.
##' @param ... optional further arguments, passed the next \code{\link{print}} method.
print.VarCorr.selfisher <- function(x, digits = max(3, getOption("digits") - 2),
				  comp = "Std.Dev.", formatter = format, ...)
{
    for (cc in names(x))  if(!is.null(x[[cc]])) {
	cat(sprintf("\n%s:\n", cNames[[cc]]))
	print(formatVC(x[[cc]],
		       digits = digits, comp = comp, formatter = formatter),
	      quote = FALSE, ...)
    }
    invisible(x)
}

