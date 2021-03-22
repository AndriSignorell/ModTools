\name{FitMod}
\alias{FitMod}
\alias{print.FitMod}
\alias{plot.FitMod}
\alias{predict.FitMod}
\alias{summary.FitMod}
\alias{drop1.FitMod}

\title{Wrapper for Several Model Functions
}
\description{Popular implementations of algorithms are characterized by partly unconventional implementations of the operating standards in R. For example, the function \code{e1057::SVM()} returns the predicted values as attributes! \cr
\code{FitMod()} is designed as a wrapping function to offer a consistent interface for a selection of most often used classification and regression models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FitMod(formula, data, ..., subset, na.action = na.pass, fitfn = NULL)

\method{predict}{FitMod}(object, ...)
\method{plot}{FitMod}(x, ...)
\method{summary}{FitMod}(object, ...)
\method{drop1}{FitMod}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{a fitted object of class \code{"FitMod"}.}
  \item{formula}{a formula expression as for classification and regression models, of the form \code{response ~ predictors}. The response should be a factor or a matrix with K columns, which will be interpreted as counts for each of K classes. See the documentation of \code{\link{formula}()} for other details.
%%     ~~Describe \code{formula} here~~
}
  \item{data}{an optional data frame in which to interpret the variables occurring in formula.
%%     ~~Describe \code{data} here~~
}
  \item{subset}{expression saying which subset of the rows of the data should be used in the fit. All observations are included by default.
%%     ~~Describe \code{subset} here~~
}
  \item{na.action}{a function to filter missing data.
%%     ~~Describe \code{na.action} here~~
}
  \item{fitfn}{code for the fitting function to be used for regression or classifying. So far implemented are: \code{lm}, \code{lmrob}, \code{poisson}, \code{quasipoisson}, \code{gamma}, \code{negbin}, \code{poisson}, \code{polr}, \code{tobit}, \code{zeroinfl}, \code{multinom}, \code{poisson}, \code{rpart}, \code{randomForest}, \code{logit}, \code{nnet}, \code{C5.0}, \code{lda}, \code{qda}, \code{svm}, \code{naive_bayes}, \code{lb}.
%%     ~~Describe \code{fitfn} here~~
}
\item{object}{the model object.}

\item{\dots}{further arguments passed to the underlying functions.}
}
\details{The function will in general return the original object, extended by a further class \code{FitMod}, which allows to capture the output and plot routines.

The classifying algorithms will at the minimum offer the predicting options \code{type = c("class", "prob")} additionally to those implemented by the underlying function.
}
\value{model object as returned by the calculating function extended with the \code{FitMod} class.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Andri Signorell <andri@signorell.net>
}

\seealso{
\code{\link{lm}}, \code{\link{rpart}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
r.lm <- FitMod(Fertility ~ ., data=swiss, fitfn="lm")

r.logit <- FitMod(diabetes ~ glucose + pressure + mass + age,
                  data=d.pima, fitfn="logit")
r.svm <- FitMod(diabetes ~ glucose + pressure + mass + age,
                  data=d.pima, fitfn="svm")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
