\name{CoeffDiffCI}
\alias{CoeffDiffCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence Interval for the Difference of Two Coefficients in a Linear Model
%%  ~~function to do ... ~~
}
\description{Calculate the confidence interval for the difference of two coefficients in a linear model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CoeffDiffCI(x, coeff, conf.level = 0.95, sides = c("two.sided", "left", "right"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the linear model object
%%     ~~Describe \code{x} here~~
}
  \item{coeff}{a vector of length two, containing either the names or the index of the two coefficients whose difference should be used
%%     ~~Describe \code{coeff} here~~
}
  \item{conf.level}{confidence level of the interval.
%%     ~~Describe \code{conf.level} here~~
}
\item{sides}{a character string specifying the side of the confidence interval, must be one of \code{"two.sided"} (default), \code{"left"} or \code{"right"}. You can specify just the initial letter. \code{"left"} would be analogue to a hypothesis of \code{"greater"} in a \code{t.test}.}

}
\details{This is quite useful in the course of the modelling process.
%%  ~~ If necessary, more details than the description above ~~
}
\value{a numeric vector with 3 elements:
  \item{mean}{mean}
  \item{lwr.ci}{lower bound of the confidence interval}
  \item{upr.ci}{upper bound of the confidence interval}
}
\author{Andri Signorell <andri@signorell.net>}

\seealso{\code{\link[car]{linearHypothesis}()}
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
