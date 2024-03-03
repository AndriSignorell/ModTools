\name{PredictCI}
\alias{PredictCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Confidence Intervals for Predictions of a GLM
%%  ~~function to do ... ~~
}
\description{Provides confidence intervals for predictions of a GLM.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PredictCI(mod, newdata, conf.level = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{the binomial model
%%     ~~Describe \code{mod} here~~
}
  \item{newdata}{the data to be predicted
%%     ~~Describe \code{newdata} here~~
}
  \item{conf.level}{confidence level of the interval. Default is 0.95.
%%     ~~Describe \code{conf.level} here~~
}
}
\details{The confidence intervals for predictions are calculated with the se of the model and the normal quantile.
%%  ~~ If necessary, more details than the description above ~~
}
\value{a matrix with 3 columns for the fit, the lower confidence interval and the upper confidence interval
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{\url{https://stackoverflow.com/questions/14423325/confidence-intervals-for-predictions-from-logistic-regression}
%% ~put references to the literature/web site here ~
}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}
\seealso{\code{\link{FitMod}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
r.logit <- FitMod(diabetes ~ age, d.pima, fitfn = "logit")
head(PredictCI(r.logit, newdata=d.pima))
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
