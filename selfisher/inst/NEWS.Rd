\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}

\name{NEWS}
\title{selfisher News}
\encoding{UTF-8}
\section{CHANGES IN VERSION 1.1.0}{
%  \subsection{NEW FEATURES}{
%    \itemize{
%      \item
%      }
%  } % new features
  \subsection{BUG FIXES}{
   \itemize{
       \item \code{predict} was not handling data-dependent predictors (e.g. \code{poly}, \code{spline}, \code{scale}) correctly
    }
  } % bug fixes
%  \subsection{USER-VISIBLE CHANGES}{
%    \itemize{
%      \item
%    }
%  }
} % 1.1.0


\section{CHANGES IN VERSION 1.0.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added data sets plus examples to go with upcoming manuscript
      }
  } % new features
%  \subsection{BUG FIXES}{
%   \itemize{
%      \item
%    }
%  } % bug fixes
%  \subsection{USER-VISIBLE CHANGES}{
%    \itemize{
%      \item
%    }
%  }
} % 1.0.0

\section{CHANGES IN VERSION 0.3.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item New function `inits` to give initial values based on guesses for l50 and SR.
      \item Now possible to include sampling fractions in trouser trawl models, paired gear models, and models with non-logit links.
    }
  } % new features
  \subsection{BUG FIXES}{
    \itemize{
      \item Nonparametric bootstrapping should be working now.
    }
  } % bug fixes
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item
    }
  }
} % 0.3.0


\section{CHANGES IN VERSION 0.2.0}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item It is now possible to do the double bootstrap while sampling hauls from separate pools and maintaining the same number of samples as in the original pools. This may be useful for gillnet analyses (see Herrmann et al. 2017).
    }
  } % new features
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed double bootstrapping algorithm
    }
  } % bug fixes
  \subsection{USER-VISIBLE CHANGES}{
    \itemize{
      \item NAs in the response variable are no longer allowed in bootstrapping
    }
  }
} % 0.2.0

