\name{Node}
\alias{Node}
\alias{Splits}

\title{Nodes and Splits in an rpart Tree
}
\description{The \code{rpart} result object has a complex and compact design. This can make practical use tedious for occasional users as it is difficult to figure out how to access some specific information. The function \code{Node()} is designed as accessor to the most important properties of a node, being a 'split' or a 'leaf' (aka. 'endnode'). It also serves as base for further convenience functions as e.g. \code{\link{LeafRates}()}.
}
\usage{
Node(x, node = NULL, type = c("all", "split", "leaf"), digits = 3)
}
\arguments{
  \item{x}{
fitted model object of class \code{rpart}.
}
  \item{node}{integer vector, defining the nodes whose details are required.
}
  \item{type}{one out of \code{"all"} (default), \code{"split"}, \code{"leaf"}, where the latter two restrict the result set to splits or end nodes only. Can be abbreviated.
}
\item{digits}{the number of digits for numeric values}
}
\details{\code{Node()} returns detailed information for a single node in the tree. It reports all the data in the summary of a node, but with the option to provide a nodelist. The structure of the result is organised as a list.
}
\value{A list containing:
\item{id}{int, id of the node}
\item{vname }{character, one out of \code{'leaf'} or \code{'split'}}
\item{isleaf }{logical, \code{TRUE} for leaves \code{FALSE} else}
\item{nobs}{integer, number of observation in the node}
\item{group  }{character, the predicted class for the node}
\item{ycount}{numeric, the number of observation per class in the node}
\item{yprob}{numeric, the relative frequencies for the each class}
\item{nodeprob }{the global probability for an observation to fall in the node}
\item{complexity}{numeric, the complexity parameter for the node}
\item{tprint}{character, the text to be printed}
}

\author{Andri Signorell <andri@signorell.net>}

\seealso{\code{\link{LeafRates}}, \code{\link{Rules}}}

\examples{
r.rpart <- FitMod(Species ~ ., data=iris, fitfn="rpart")
# return Node nr. 3
Node(r.rpart, node=3)

r.rp <- FitMod(Type ~ ., data = d.glass, fitfn="rpart")
# return all the splits
Node(r.rpart, type="split")
}
\keyword{ tree }

