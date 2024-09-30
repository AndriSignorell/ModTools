\name{LeafRates}
\alias{Purity}
\alias{LeafRates}
\alias{plot.LeafRates}

\title{Leafrates for the Nodes of an 'rpart' Tree
}
\description{Return the frequencies of correct and wrong classifications in given node(s) in tabular form. The 'purity', denoting the relative frequency of correctly classified elements, is a useful information for the interpretation of regression and classification trees and a measure for its quality.
}
\usage{
LeafRates(x)

\method{plot}{LeafRates}(x, col = NULL, which = c("rel", "abs"),
                         layout = NULL, ylim = NULL, \dots)
}
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
}

\seealso{\code{\link{Node}}, \code{\link{Rules}}
}
\examples{
r.rp <- FitMod(Species ~ ., data=iris, fitfn="rpart")
LeafRates(r.rp)

plot(LeafRates(r.rp))
}
\keyword{ tree }
