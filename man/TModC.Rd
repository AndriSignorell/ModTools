\name{TModC}
\alias{TModC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Puts together results for classification models
%%  ~~function to do ... ~~
}
\description{Puts together results for classification models
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TModC(..., newdata = NULL, reference = NULL,
        type = C("ROC", "measures", "varimp"), col = NULL,
        args.legend = NULL, ord = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{a list with models
%%     ~~Describe \code{\dots} here~~
}
  \item{newdata}{the data to use for predicting
%%     ~~Describe \code{newdata} here~~
}
  \item{reference}{the reference values
%%     ~~Describe \code{reference} here~~
}
  \item{type}{type of result one out of ROC, measures, varimp
%%     ~~Describe \code{type} here~~
}
  \item{col}{the color for the plotlines in ROC
%%     ~~Describe \code{col} here~~
}
  \item{args.legend}{the legend in the ROC plot
%%     ~~Describe \code{args.legend} here~~
}
  \item{ord}{the order of the results, one out of
%%     ~~Describe \code{ord} here~~
}
}
\details{Put together the results to compare classificators.
%%  ~~ If necessary, more details than the description above ~~
}
\value{plot and table
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}
\seealso{\code{\link[DescTools]{TMod}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
d.pim <- SplitTrainTest(d.pima, p = 0.2)
mdiab <- formula(diabetes ~ pregnant + glucose + pressure + triceps +
                            insulin + mass + pedigree + age)

r.glm <- FitMod(mdiab, data=d.pim$train, fitfn="logit")
r.rp <- FitMod(mdiab, data=d.pim$train, fitfn="rpart")
mods <- list(glm=r.glm, rp=r.rp)

TModC(mods, ord="a")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}% use one of  RShowDoc("KEYWORDS")
