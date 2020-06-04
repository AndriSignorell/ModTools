\name{Node}
\alias{Node}
\alias{Splits}
\alias{LeafRates}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Nodes, Splits and Leafrates in an rpart Tree
%%  ~~function to do ... ~~
}
\description{The \code{rpart} object has an extremely compact and abstract design. This sometimes makes practical use  tedious as it can be difficult to figure out how to extract some specific information. The function \code{Node()} is designed as accessor to the most important properties of a node, being a split or a leaf. It serves as base for further convenience functions as e.g. \code{LeafRates()}, which returns the misclassification rates in all end nodes.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Node(x, node=NULL, type=c("all", "split", "leaf"), digits=3)

LeafRates(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
fitted model object of class \code{rpart}.
%%     ~~Describe \code{x} here~~
}
  \item{node}{integer vector, defining the nodes whose details are required.
%%     ~~Describe \code{node} here~~
}
  \item{type}{one out of \code{"all"} (default), \code{"split"}, \code{"leaf"}, where the latter two restrict the result set to splits or end nodes only. Can be abbreviated.
%%     ~~Describe \code{node} here~~
}
}
\details{The function returns detailed information for a single node in the tree. It reports all the data in the summary of a node, but with the option to provide a nodelist. The structure of the result is organised as a list and therefore should be easier to understand than the original.

The row.names of frame contain the (unique) node numbers that follow a binary ordering indexed by node depth. Columns of frame include var, a factor giving the names of the variables used in the split at each node (leaf nodes are denoted by the level "<leaf>"), n, the number of observations reaching the node, wt, the sum of case weights for observations reaching the node, dev, the deviance of the node, yval, the fitted value of the response at the node, and splits, a two column matrix of left and right split labels for each node. Also included in the frame are complexity, the complexity parameter at which this split will collapse, ncompete, the number of competitor splits recorded, and nsurrogate, the number of surrogate splits recorded.
%%  ~~ If necessary, more details than the description above ~~
}
\value{An object of class \code{rpart}. See \code{\link{rpart.object}}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}

\seealso{\code{\link{Purity}}, \code{\link{Rules}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
r.rpart <- FitMod(Species ~ ., data=iris, fitfn="rpart")
Node(r.rpart, node=3)

r.rp <- FitMod(Type ~ ., data = d.glass, fitfn="rpart")
Node(r.rpart, type="splits")

LeafRates(r.rp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }% __ONLY ONE__ keyword per line
