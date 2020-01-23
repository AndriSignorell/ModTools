\name{Rules}
\alias{Rules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract Rules from a rpart Object
%%  ~~function to do ... ~~
}
\description{Extract Rules from a rpart Object
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Rules(x, node = NULL, leafonly = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the rpart object to extract the rules from
%%     ~~Describe \code{x} here~~
}
\item{node}{integer vector, defining the nodes whose details are required.}
\item{leafonly}{boolean, defining if only end nodes ("leafs") should be returned.}
}
\value{a list with the rules
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
 \code{\link{rpart}}
 }
\examples{
r.rp <- FitMod(diabetes ~ ., data=d.pima, fitfn="rpart")
Rules(r.rp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models}% use one of  RShowDoc("KEYWORDS")
