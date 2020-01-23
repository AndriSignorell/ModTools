\name{Purity}
\alias{Purity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Purity in Nodes
%%  ~~function to do ... ~~
}
\description{Displays the frequencies of classifications in a given node. The purity denotes the relative frequency of correctly classified elements.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Purity(x, leaves = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
fitted model object of class \code{rpart}. This is assumed to be the result of some function that produces an object with the same named components as that returned by the \code{rpart} function.
%%     ~~Describe \code{x} here~~
}
  \item{leaves}{boolean, if set to \code{TRUE}, only purity of the end nodes will be displayed, else all the nodes will be included.
%%     ~~Describe \code{node} here~~
}
}
\details{The purity is an interesting measure for the quality of a tree.
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
Purity(r.rpart)

# all the nodes
Purity(r.rpart, leaves=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }% __ONLY ONE__ keyword per line
