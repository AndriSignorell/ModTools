\name{Response}
\alias{Response}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{extracts the response from several models
%%  ~~function to do ... ~~
}
\description{extracts the response from several models
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Response(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the model to use
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{more arguments
%%     ~~Describe \code{\dots} here~~
}
}
\details{extracts the response from several models
%%  ~~ If necessary, more details than the description above ~~
}
\value{the response
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{Andri Signorell <andri@signorell.net>
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{model.frame}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
r.rpart <- FitMod(diabetes ~ ., d.pima, fitfn="rpart")
Response(r.rpart)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
