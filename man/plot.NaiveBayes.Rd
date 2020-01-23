\name{plot.NaiveBayes}
\alias{plot.NaiveBayes}
\title{Plot Method for NaiveBayes Object}
\usage{
\method{plot}{NaiveBayes}(x, which = NULL, ask = FALSE, legend = TRUE,
  legend.box = FALSE, arg.num = list(), arg.cat = list(), ...)
}
\arguments{
\item{x}{object of class inheriting from \code{"NaiveBayes"}.}

\item{which}{variables to be plotted (all by default). This can be any valid indexing vector or vector containing names of variables.}

\item{ask}{logical; if \code{TRUE}, the user is asked before each plot, see \code{\link[graphics]{par}(ask=.)}.}

\item{legend}{logical; if \code{TRUE} a \code{\link{legend}} will be be plotted.}

\item{legend.box}{logical; if \code{TRUE} a box will be drawn around the legend.}

\item{arg.num}{other parameters to be passed as a named list to \code{\link[graphics]{matplot}}.}

\item{arg.cat}{other parameters to be passed as a named list to \code{\link[graphics]{mosaicplot}}.}

\item{...}{not used.}
}
\description{
Plot method for objects of class \code{"NaiveBayes"} designed for a quick look at the marginal probabilities of predictor variables given the class.
}
\details{
Conditional probabilities are visualised by \code{\link[graphics]{matplot}} (for metric predictors) and  \code{\link[graphics]{mosaicplot}} (for categorical predictors). In case of non parametric estimation of class-conditional densities, the bandwidths are reported. Nothing is returned.
}
\examples{
data(iris)
iris2 <- cbind(iris, New = sample(letters[1:3], 150, TRUE))
nb <- NaiveBayes(Species ~ ., data = iris2)
plot(nb)
plot(nb, which = c(1, 2), arg.num = list(col = 1:3, lty = 1,
     main = "Naive Bayes Plot"))
plot(nb, which = "Newb", arg.cat = list(color = 4:7))
}
\author{
Michal Majka, \email{michalmajka@hotmail.com}
}
\seealso{
\code{\link{NaiveBayes}}, \code{\link{predict.NaiveBayes}}, \code{\link[naivebayes]{tables}}
}

