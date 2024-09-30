\name{CoeffDiffCI}
\alias{CoeffDiffCI}
\title{Confidence Interval for the Difference of Two Coefficients in a Linear Model
}
\description{Calculate the confidence interval for the difference of two coefficients in a linear model.
}
\usage{
CoeffDiffCI(x, coeff, conf.level = 0.95, sides = c("two.sided", "left", "right"))
}
\arguments{
  \item{x}{the linear model object
}
  \item{coeff}{a vector of length two, containing either the names or the index of the two coefficients whose difference should be used
}
  \item{conf.level}{confidence level of the interval.
}
\item{sides}{a character string specifying the side of the confidence interval, must be one of \code{"two.sided"} (default), \code{"left"} or \code{"right"}. You can specify just the initial letter. \code{"left"} would be analogue to a hypothesis of \code{"greater"} in a \code{t.test}.}

}
\details{This is quite useful in the course of the modelling process.
}
\value{a numeric vector with 3 elements:
  \item{mean}{mean}
  \item{lwr.ci}{lower bound of the confidence interval}
  \item{upr.ci}{upper bound of the confidence interval}
}
\author{Andri Signorell <andri@signorell.net>}

\seealso{\code{\link[car]{linearHypothesis}()}
}
\examples{
# get some model first...
r.lm <- FitMod(Fertility ~ ., data=swiss, fitfn="lm")

# calculate the confidence interval for the difference of the
# coefficients Examination and Education
CoeffDiffCI(r.lm, c("Examination", "Education"))

# the test could be calculated as
car::linearHypothesis(r.lm, "Education = Examination")
}
