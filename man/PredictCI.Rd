\name{PredictCI}
\alias{PredictCI}
\title{Confidence Intervals for Predictions of a GLM
}
\description{Provides confidence intervals for predictions of a GLM.
}
\usage{
PredictCI(mod, newdata, conf.level = 0.95)
}
\arguments{
  \item{mod}{the binomial model
}
  \item{newdata}{the data to be predicted
}
  \item{conf.level}{confidence level of the interval. Default is 0.95.
}
}
\details{The confidence intervals for predictions are calculated with the se of the model and the normal quantile.
}
\value{a matrix with 3 columns for the fit, the lower confidence interval and the upper confidence interval
}
\references{\url{https://stackoverflow.com/questions/14423325/confidence-intervals-for-predictions-from-logistic-regression}
}
\author{Andri Signorell <andri@signorell.net>
}
\seealso{\code{\link{FitMod}}
}
\examples{
r.logit <- FitMod(diabetes ~ age, d.pima, fitfn = "logit")
head(PredictCI(r.logit, newdata=d.pima))
}
