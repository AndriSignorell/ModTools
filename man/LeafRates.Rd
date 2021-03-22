\name{LeafRates}
\alias{Purity}
\alias{LeafRates}
\alias{plot.LeafRates}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Leafrates for the Nodes of an 'rpart' Tree
%%  ~~function to do ... ~~
}
\description{Return the frequencies of correct and wrong classifications in given node(s) in tabular form. The 'purity', denoting the relative frequency of correctly classified elements, is a useful information for the interpretation of regression and classification trees and a measure for its quality.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LeafRates(x)

\method{plot}{LeafRates}(x, col = NULL, which = c("rel", "abs"),
                         layout = NULL, ylim = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{fitted model object of class \code{rpart}.}
  \item{col}{color for the bars in the plot}
  \item{which}{one out of \code{"rel"} or \code{"abs"}, denoting whether relative or absolute frequencies should be used for the plot.}
  \item{layout}{vector defining the layout}
  \item{ylim}{the y limits of the plot.}
  \item{\dots}{further arguments (not used).}
}

\details{The result comprises absolute and relative frequencies per leaf.
}
\value{A list with 5 elements consisting of:
\item{node }{the node id (of the leaf)}
\item{freq }{the absolute frequency of correct and wrong classifications}
\item{p.row }{the relative frequency of correct and wrong classifications}
\item{mfreq}{the total number of cases}
\item{mperc}{the percentage of the sample in the leaf}
}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}

\seealso{\code{\link{Node}}, \code{\link{Rules}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
r.rp <- FitMod(Species ~ ., data=iris, fitfn="rpart")
LeafRates(r.rp)

plot(LeafRates(r.rp))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }% __ONLY ONE__ keyword per line
