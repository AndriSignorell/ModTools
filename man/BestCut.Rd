\name{BestCut}
\alias{BestCut}
\title{Best Cutpoint for a ROC Curve
}
\description{Returns the best cutpoint for a given classification model.
}
\usage{
BestCut(x, method = c("youden", "closest.topleft"))
}
\arguments{
  \item{x}{a roc object from the roc function
}
  \item{method}{one of \code{"youden"} or \code{"closest.topleft"}, controls how the
  optimal threshold is determined. See details.
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
}
\value{the threshold value
}
\references{
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  \doi{10.1186/1471-2105-12-77}.
}
\author{Robin Xavier <pROC-cran@xavier.robin.name>,
Andri Signorell <andri@signorell.net> (interface)
}

\seealso{\code{\link{ROC}}
}
\examples{
r.glm <- FitMod(diabetes ~ ., data = d.pima, fitfn="logit")

ROC(r.glm)
BestCut(ROC(r.glm))
}
