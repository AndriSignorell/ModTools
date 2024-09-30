\name{Rules}
\alias{Rules}
\title{Extract Rules from 'rpart' Object
}
\description{Extract rules from an rpart object. This can be useful, if the rules must be implemented in another system. The rules contain all the criteria for the binary splits of an rpart tree from the root node down to the specified leaf.
}
\usage{
Rules(x, node = NULL, leafonly = FALSE)
}
\arguments{
  \item{x}{the rpart object to extract the rules from
}
\item{node}{integer vector, defining the nodes whose details are required.}
\item{leafonly}{boolean, defining if only the rules leading to end nodes ("leafs") should be returned.}
}
\details{The function builds upon the original function \code{\link[rpart]{path.rpart}}, which is bulky in some situations.
}
\value{a list with the rules
\item{frame }{the frame of the rpart}
\item{ylevels}{the y values of the node}
\item{ds.size}{the size of the dataset}
\item{path}{a list of character vecotrs containing the rules}
}
\author{Andri Signorell <andri@signorell.net>
}

\seealso{
 \code{\link[rpart]{rpart}}, \code{\link[rpart]{path.rpart}}
 }
\examples{
r.rp <- FitMod(diabetes ~ ., data=d.pima, fitfn="rpart")
Rules(r.rp)
}
\keyword{ models}
