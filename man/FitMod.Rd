\name{FitMod}
\alias{FitMod}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Wrapper for Several Model Functions
%%  ~~function to do ... ~~
}
\description{Wrap some classifiers to give them a homogen interface.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FitMod(formula, data, ..., subset, na.action = na.pass, fitfn = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula expression as for FitModession models, of the form \code{response ~ predictors}. The response should be a factor or a matrix with K columns, which will be interpreted as counts for each of K classes. See the documentation of \code{\link{formula}()} for other details.
%%     ~~Describe \code{formula} here~~
}
  \item{data}{an optional data frame in which to interpret the variables occurring in formula.
%%     ~~Describe \code{data} here~~
}
  \item{\dots}{additional arguments for nnet
%%     ~~Describe \code{\dots} here~~
}
  \item{subset}{expression saying which subset of the rows of the data should be used in the fit. All observations are included by default.
%%     ~~Describe \code{subset} here~~
}
  \item{na.action}{a function to filter missing data.
%%     ~~Describe \code{na.action} here~~
}
  \item{fitfn}{fitting function to use for regression or classifying, one out of \code{lm}, \code{glm}, \code{polr}, \code{rpart}, \code{randomForest}, \code{logit}
%%     ~~Describe \code{fitfn} here~~
}
}
\details{Needed to make it easy to use
%%  ~~ If necessary, more details than the description above ~~
}
\value{the model
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}

\seealso{
\code{\link{lm}}, \code{\link{rpart}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
r.lm <- FitMod(Fertility ~ ., data=swiss, fitfn="lm")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
