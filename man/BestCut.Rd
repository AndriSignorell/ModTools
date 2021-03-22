\name{BestCut}
\alias{BestCut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Best Cutpoint for a ROC Curve
%%  ~~function to do ... ~~
}
\description{Returns the best cutpoint for a given classification model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BestCut(x, method = c("youden", "closest.topleft"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a roc object from the roc function
%%     ~~Describe \code{x} here~~
}
  \item{method}{one of \code{"youden"} or \code{"closest.topleft"}, controls how the
  optimal threshold is determined. See details.
%%     ~~Describe \code{method} here~~
}
}
\details{
  The \code{method} argument controls how the
  optimal threshold is determined.
  \describe{
    \item{'\code{youden}'}{
      Youden's J statistic (Youden, 1950) is employed. The optimal cut-off is the threshold
      that maximizes the distance to the identity (diagonal) line. Can
      be shortened to \dQuote{y}.

      The optimality criterion is:
      \deqn{max(sensitivities + specificities)}{max(sensitivities + specificities)}
    }
    \item{'\code{closest.topleft}'}{
      The optimal threshold is the point closest to the top-left part of
      the plot with perfect sensitivity or specificity. Can be shortened
      to \dQuote{c} or \dQuote{t}.

      The optimality criterion is:
      \deqn{min((1 - sensitivities)^2 + (1- specificities)^2)}{min((1 - sensitivities)^2 + (1- specificities)^2)}
    }
}
%%  ~~ If necessary, more details than the description above ~~
}
\value{the threshold value
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \href{http://dx.doi.org/10.1186/1471-2105-12-77}{10.1186/1471-2105-12-77}.
%% ~put references to the literature/web site here ~
}
\author{Robin Xavier <pROC-cran@xavier.robin.name>,
Andri Signorell <andri@signorell.net> (interface)
%%  ~~who you are~~
}

\seealso{\code{\link{ROC}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
r.glm <- FitMod(diabetes ~ ., data = d.pima, fitfn="logit")

ROC(r.glm)
BestCut(ROC(r.glm))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
