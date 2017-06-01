### Copied from bootMer() 
.simpleCap <- function(x) {
  paste0(toupper(substr(x, 1,1)), substr(x, 2, 1000000L), collapse=" ")
}

##' a function taking a fitted \code{selfisher} object as input and returning the
##' L50 and SR estimates as a named numeric vector.
##' @param x a fitted \code{selfisher} object
##' @export
L50SR <- function(x) {
	L50 <- summary(x$sdr, "report")[which(x$obj$report()$retp==0.5),1]
	SR <- summary(x$sdr, "report")["SR",1]
	return(c("L50"=L50,"SR"=SR)) 
}

##' a function that refits the same model to a new data set
##' @param x a fitted \code{selfisher} object
##' @param newdata
##' @export
refit.selfisher <- funciton(x, newdata) {
	
}

##' Perform bootstrap.
##'
##' @param x a fitted \code{selfisher} object
##' @param FUN a function taking a fitted 
##' \code{selfisher} object as input and returning the
##' \emph{statistic} of interest, which must be a (possibly named) numeric vector.
##' @param nsim number of simulations, positive integer
##' @param seed optional argument to \code{\link{set.seed}}
##' @param type character string specifying the type of
##' bootstrap, \code{double}(the defualt) as defined in gear selectivity literature (e.g. Millar 1993),
##' \code{"parametric"} or \code{"nonparametric"}; partial matching is allowed.
##' @details the default bootstrap type "double" is specific to fisheries literature.
##' @export
##' @importFrom lme4 refit
bootSel <- function(x, FUN = L50SR, nsim = 1, seed = NULL,
                   type=c("double", "parametric", "nonparameteric"),
                   verbose = FALSE,
                   .progress = "none", PBargs=list(),
                   parallel = c("no", "multicore", "snow"),
                   ncpus = getOption("boot.ncpus", 1L), cl = NULL)
{
    stopifnot((nsim <- as.integer(nsim[1])) > 0)
    if (.progress!="none") { ## progress bar
        pbfun <- get(paste0(.progress,"ProgressBar"))
        setpbfun <- get(paste0("set",.simpleCap(.progress),"ProgressBar"))
        pb <- do.call(pbfun,PBargs)
    }
    if (missing(parallel)) parallel <- getOption("boot.parallel", "no")
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    do_parallel <- (parallel != "no" && ncpus > 1L)
    if (do_parallel) {
        if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow") have_snow <- TRUE
        if (!(have_mc || have_snow))
        do_parallel <- FALSE # (only for "windows")
    }
    if (do_parallel && .progress != "none")
        message("progress bar disabled for parallel operations")

    FUN <- match.fun(FUN)
    if(!is.null(seed)) set.seed(seed)
    else if(!exists(".Random.seed", envir = .GlobalEnv))
        runif(1) # initialize the RNG if necessary

    mc <- match.call()
    t0 <- FUN(x)
    if (!is.numeric(t0))
        stop("bootSel currently only handles functions that return numeric vectors")

    mle <- x$obj$env$parList(x$fit$par, x$fit$parfull)

    type <- match.arg(type)
    if (type=="parametric") {
       argList <- list(x, nsim=nsim)
       ss <- do.call(simulate,argList)
    } else {
        yobs <- x$frame[,x$modelInfo$respCol]
        if (type=="double") {
           #resample hauls
           hauls <- unique(x$frame[,"(haul)"])
           if(length(hauls)<=1) stop("Double bootstrap is only useful for multiple hauls. Maybe you want 'nonparameteric'.")
           splith <- split(x$frame, x$frame[,"(haul)"])
           #create nsim new frames in ss
           newhauls <- replicate(nsim, sample(hauls, length(hauls), replace=TRUE))

           #within hauls, resample obs for each length class
           newframe <- apply(newhauls, 2, function(i){ do.call(rbind, splith[i])})
           ss <- lapply(newframe, function(z) {
                   #overwrite the response variable
                   z[,x$modelInfo$respCol] <- rbinom(nrow(z), size=z[,"(total)"], prob=z[,x$modelInfo$respCol])/z[,"(total)"]
                   return(z)
                 })
        } else {
            if (type=="nonparameteric") {
              
              ss <- replicate(nsim, function() {
                      z  <- x$frame
                      #overwrite the response variable
                      z[,z$modelInfo$respCol] <- rbinom(length(z), z[,"(total)"], z[,z$modelInfo$respCol])/z[,"(total)"]
                      return(z)
                    })

            } else {
            stop("unknown 'type' specified in call to bootSel")
              }
          }
      }

    # define ffun as a closure containing the referenced variables
    # in its scope to avoid explicit clusterExport statement
    # in the PSOCKcluster case
    ffun <- local({
      FUN
      refit
      x
      ss
      verbose
      do_parallel
      length.t0 <- length(t0)
      function(i) {
        ret <- tryCatch(FUN(refit(x,ss[[i]])), error=function(e)e)
        if(verbose) { cat(sprintf("%5d :",i)); str(ret) }
        if (!do_parallel && .progress!="none") { setpbfun(pb,i/nsim) }
        if (inherits(ret, "error"))
            structure(rep(NA, length.t0), "fail.msgs" = ret$message)
        else
            ret
    }})

    simvec <- seq_len(nsim)
     res <- if (do_parallel) {
        if (have_mc) {
            parallel::mclapply(simvec, ffun, mc.cores = ncpus)
        } else if (have_snow) {
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
                ## explicit export of the lme4 namespace since most FUNs will probably
                ## use some of them
                parallel::clusterExport(cl, varlist=getNamespaceExports("selfisher"))
                if(RNGkind()[1L] == "L'Ecuyer-CMRG")
                    parallel::clusterSetRNGStream(cl)
                res <- parallel::parLapply(cl, simvec, ffun)
                parallel::stopCluster(cl)
                res
            } else parallel::parLapply(cl, simvec, ffun)
        }
    } else lapply(simvec, ffun)

    t.star <- do.call(cbind,res)
    rownames(t.star) <- names(t0)
    if ((numFail <- sum(bad.runs <- apply(is.na(t.star),2,all)))>0) {
        warning("some bootstrap runs failed (",numFail,"/",nsim,")")
        fail.msgs <- vapply(res[bad.runs],FUN=attr,FUN.VALUE=character(1),
                            "fail.msgs")
    } else fail.msgs <- NULL
    ## boot() ends with the equivalent of
    ## structure(list(t0 = t0, t = t.star, R = R, data = data, seed = seed,
    ##		      statistic = statistic, sim = sim, call = call,
    ##		      ran.gen = ran.gen, mle = mle),
    ##		 class = "boot")
    s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = x$frame,
		   seed = .Random.seed,
		   statistic = FUN, call = mc,
		   ## these two are dummies
		   ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle),
	      class = "boot")
    attr(s,"bootFail") <- numFail
    attr(s,"boot.fail.msgs") <- fail.msgs
    s
} ## {bootSel}

##' @S3method as.data.frame boot
as.data.frame.boot <- function(x,...) {
  as.data.frame(x$t)
}

