% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.selfisher}
\alias{predict.selfisher}
\title{prediction}
\usage{
\method{predict}{selfisher}(
  object,
  newdata = NULL,
  se.fit = FALSE,
  re.form,
  allow.new.levels = FALSE,
  type = c("response", "selection", "prob", "ratio", "link"),
  na.action = na.pass,
  debug = FALSE,
  ...
)
}
\arguments{
\item{object}{a \code{selfisher} object}

\item{newdata}{new data for prediction}

\item{se.fit}{return the standard errors of the predicted values?}

\item{re.form}{(not yet implemented) specify which random effects to condition on when predicting. For now, all random effects are included.}

\item{allow.new.levels}{(not yet implemented) allow previously unobserved levels in random-effects grouping variables?}

\item{type}{\itemize{
\item return expected response value ("response": see details below),
\item predicted selection curve ("selection": r),
\item relative fishing power of the test gear ("prob"),
\item catch ratio from catch comparison models ("ratio": r/(1-r)).
}
Some types (ratio, link) might not make sense to use with psplit models.}

\item{debug}{(logical) return the \code{TMBStruc} object that will be
used internally for debugging?}

\item{\dots}{unused - for method compatibility}
}
\description{
prediction
}
\details{
Predicting with type="response" returns values comparable to the
response variable (the left-hand side of the model's rformula);
that is  pr/(pr+1-p) in a model with psplit=TRUE
or r in a model with psplit=FALSE.
\itemize{
\item Prediction of new random effect levels is possible as long as the model specification (fixed effects and parameters) is kept constant.
However, to ensure intentional usage, a warning is triggered if \code{allow.new.levels=FALSE} (the default).
\item Prediction using "data-dependent bases" (variables whose scaling or transformation depends on the original data, e.g. \code{\link{poly}}, \code{\link[splines]{ns}}, or \code{\link{poly}}) should work properly; however, users are advised to check results extra-carefully when using such variables. Models with different versions of the same data-dependent basis type in different components (e.g. \code{formula= y ~ poly(x,3), dispformula= ~poly(x,2)}) will probably \emph{not} produce correct predictions.
}
##' @examples data(haddock)
dat <- transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))
m1 <- selfisher(prop~Lengths, p=~1, psplit=TRUE, total=tot, dat)
nd <- data.frame(Lengths=20:50, tot=100)
predict(m1, newdata=nd, se.fit=TRUE)
}
