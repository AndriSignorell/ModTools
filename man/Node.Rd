\name{Node}
\alias{Node}
\alias{Splits}
\alias{LeafRates}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Nodes, Splits and Leafrates in an rpart Tree
%%  ~~function to do ... ~~
}
\description{The \code{rpart} object is designed extremely densely, which sometimes makes practical use a little tedious. The functions \code{Node()} and \code{Splits()} extract the specific information from a \code{rpart} object for easier access. \code{LeafRates()} returns the misclassification rates in all end nodes.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Node(x, node, digits = getOption("digits"))
Splits(x)

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
  \item{digits}{the number of digits of numbers to print.
%%     ~~Describe \code{digits} here~~
}
}
\details{The function returns detailed information for single nodes in the tree. It reports all the data in the summary of a node, but with the option to provide a node.
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
Splits(r.rp)
LeafRates(r.rp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }% __ONLY ONE__ keyword per line
